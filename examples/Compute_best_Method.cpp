#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"


int main( int argc, char *argv[] )
{
    std::cout << std::endl;
    std::cout << "-----------RBM-Clyde Adaptive Model starts-------------" << std::endl << std::endl;
    //std::cout << "-----------Computing best method-------" << std::endl;

    std::cout << " Initializing common variables " << std::endl << std::endl;

    std::string filecfg = argv[1];
    // std::string mode = argv[2];
    prob_settings settings;
    int ndim = 2;   //dimension of the problem (1D/2D/3D)

    //Reading configuration file
    Read_cfg( filecfg, settings );
    // Config_stream ( settings );
    int Nm_POD, Nm_DMD, Nm_mrDMD, Nm_RDMD;
    std::vector<int> Nm_SPOD = {};

    // Calculate number of grid points
    int Nr = N_gridpoints ( settings.in_file );
    std::cout << "Number of grid points : " << Nr << std::endl;
    
    // Reading coordinates
    std::cout << "Reading Coordinates ... \t ";
    Eigen::MatrixXd Coords = read_col( settings.in_file, Nr, settings.Cols_coords );
    std::cout << "Done " << std::endl;

    // Create matrix of snapshots
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings.Ns, settings.Ds, settings.nstart,
                                        settings.Cols,
                                        settings.in_file,
                                        settings.flag_prob);

    double tol_rec = 0.05;
    double tol = 1e-1;
    double t_0 = (double)settings.nstart*settings.Dt_cfd;
    Eigen::VectorXd t_vec( settings.Ns );
    t_vec(0) = t_0;

    for ( int i = 1; i < settings.Ns; i++ )
        t_vec(i) = t_vec(i-1) + settings.Dt_cfd*settings.Ds;

    std::cout << std::endl;
    std::cout << "Initialized vector of times " << std::endl;

    Eigen::VectorXd mean = sn_set.rowwise().mean();
    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(settings.Ns);

    for ( int i = 0; i < settings.Ns; i ++ )
    {
        norm_sn_set(i) = sn_set.col(i).norm();
    }


    
    Eigen::VectorXd Err_POD_time = Eigen::VectorXd::Zero(settings.Ns);

    std::vector<int> Nf(4);
    Nf[0] = std::ceil(settings.Ns/10.0);
    Nf[1] = std::ceil(settings.Ns/2.0);
    Nf[2] = std::ceil(2.0*settings.Ns/3.0);
    Nf[3] = settings.Ns;
    Eigen::MatrixXd Err_SPOD_time = Eigen::MatrixXd::Zero(settings.Ns,Nf.size());


//Defining common scope for POD and SPOD
    {
        Eigen::VectorXd lambda(settings.Ns);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
//------------------>POD
        {
            for ( int i = 0; i < settings.Ns; i++ )
                sn_set.col(i) -= mean;

            std::cout << "Extracting  POD basis ... " << "\t";        

            Eigen::MatrixXd Phi = SPOD_basis( sn_set,
                                    lambda, K_pc, eig_vec,
                                    0,
                                    settings.flag_bc, 
                                    settings.flag_filter,  
                                    settings.sigma);

            std::cout << " Done! " << std::endl << std::endl;




            int Nrec = Nmod( settings.En, K_pc);
            Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nrec, Nrec);
            std::cout << "Number of modes for the desired energy content : " << Nrec << std::endl;

            // for ( int i = 0; i < settings.Ns; i++ )
            //     sn_set.col(i) += mean;

            for ( int i = 0; i < Nrec; i++ )
                Sig(i,i) = std::sqrt(lambda(i));

            Eigen::MatrixXd Err_POD_map = sn_set - Phi.leftCols(Nrec)*Sig*eig_vec.transpose().topRows(Nrec);

            for ( int i = 0; i < settings.Ns; i++ )
            {
                for ( int j = 0; j < ndim*Nr; j++ )
                    if ( sn_set(j,i) < tol )
                        Err_POD_time(i) += Err_POD_map(j,i)*Err_POD_map(j,i);
                    else    
                        Err_POD_time(i) += Err_POD_map(j,i)*Err_POD_map(j,i);// /(sn_set(j,i)*sn_set(j,i));


                Err_POD_time(i) = std::sqrt(Err_POD_time(i))/norm_sn_set(i);
            }


            Nm_POD = Nrec;

        }

//------------------>SPOD
        {
            int Nrec = 0;

            for ( int nfj = 0; nfj < Nf.size(); nfj++ )
            {
                // std::cout << "Subtracting mean from snapshots ... " << std::endl << std::endl;

                // for ( int nt = 0; nt < settings.Ns; nt++ )
                //     sn_set.col(nt) -= mean;

                std::cout << "Extracting  SPOD " << Nf[nfj] << " basis ... " << "\t";        

                Eigen::MatrixXd Phi = SPOD_basis( sn_set,
                                        lambda, K_pc, eig_vec,
                                        Nf[nfj],
                                        settings.flag_bc, 
                                        settings.flag_filter,  
                                        settings.sigma);

                std::cout << " Done! " << std::endl << std::endl;
                
                Nrec = Nmod( settings.En, K_pc);
                Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nrec, Nrec);
                std::cout << "Number of modes for the desired energy content : " << Nrec << std::endl;
    
                // for ( int i = 0; i < settings.Ns; i++ )
                //     sn_set.col(i) += mean;
    
                for ( int i = 0; i < Nrec; i++ )
                    Sig(i,i) = std::sqrt(lambda(i));
    
                std::cout << "Computing error in time ..." << std::endl;
    
                Eigen::MatrixXd Err_SPOD_map = sn_set - Phi.leftCols(Nrec)*Sig*eig_vec.transpose().topRows(Nrec);
    
                for ( int i = 0; i < settings.Ns; i++ )
                {
                    for ( int j = 0; j < ndim*Nr; j++ )
                        if ( sn_set(j,i) < tol )    
                            Err_SPOD_time(i,nfj) += Err_SPOD_map(j,i)*Err_SPOD_map(j,i);
                        else
                            Err_SPOD_time(i,nfj) += Err_SPOD_map(j,i)*Err_SPOD_map(j,i); // /(sn_set(j,i)*sn_set(j,i));  
    
                    Err_SPOD_time(i,nfj) = std::sqrt(Err_SPOD_time(i,nfj))/norm_sn_set(i);
                }

                Nm_SPOD.push_back(Nrec);
            }

        }

    }

//Adding again mean over snapshots to perform DMD and its variants
    for ( int i = 0; i < settings.Ns; i++ )
    sn_set.col(i) += mean;

    Eigen::VectorXd Err_DMD_time = Eigen::VectorXd::Zero(settings.Ns);
//Defining scope for DMD ( Rank=SVHT, Coeffs = OPT )
    {

        Eigen::VectorXd lambda_POD;
        Eigen::MatrixXd eig_vec_POD;
        Eigen::VectorXcd lambda_DMD;
        Eigen::MatrixXcd eig_vec_DMD;

        std::cout << "Extracting basis DMD using SVHT ... " << "\t";        
        Eigen::MatrixXcd Phi;
        Eigen::VectorXcd alfa;

        Phi = DMD_basis( sn_set,
                        lambda_DMD,
                        eig_vec_DMD,
                        lambda_POD,
                        eig_vec_POD,
                        0 );

        std::cout << " Done! " << std::endl << std::endl;
        
        int Nm = Phi.cols();
        std::cout << "Number of modes extracted : " << Nm << std::endl;

        // Eigen::VectorXcd omega(Nm);
        // for ( int i = 0; i < Nm; i++ )
        //         omega(i) = std::log(lambda_DMD(i))/(settings.Dt_cfd*settings.Ds);

        std::cout << "Calculating coefficients DMD ... " << "\t";
            
        alfa = Calculate_Coefs_DMD_exact ( sn_set.leftCols(settings.Ns-1),  
                                                            lambda_DMD, 
                                                            Phi );
            
        std::cout << " Done! " << std::endl << std::endl;

        Eigen::MatrixXcd V_and(lambda_DMD.size(), settings.Ns);
        
        for ( int i = 0; i < lambda_DMD.size(); i++ )
        {
            for ( int j = 0; j < settings.Ns; j++ )
                V_and(i,j) = std::pow(lambda_DMD(i), (double)j);                                                                                         
        }        

        Eigen::MatrixXcd Psi = Eigen::MatrixXd::Zero(alfa.size(), settings.Ns);
        for ( int i = 0; i < settings.Ns; i++ )
            Psi.col(i) = alfa.cwiseProduct(V_and.col(i));
        
        Eigen::MatrixXcd D_dmd = Phi*Psi;
        Eigen::MatrixXd Err_DMD_map = sn_set - D_dmd.real();

        for ( int i = 0; i < settings.Ns; i++ )
        {
            for ( int j = 0; j < ndim*Nr; j++ )
                if ( sn_set(j,i) < tol )
                    Err_DMD_time(i) += Err_DMD_map(j,i)*Err_DMD_map(j,i);
                else    
                    Err_DMD_time(i) += Err_DMD_map(j,i)*Err_DMD_map(j,i); // /(sn_set(j,i)*sn_set(j,i));  

            Err_DMD_time(i) = std::sqrt(Err_DMD_time(i))/norm_sn_set(i);
        }
        
        Nm_DMD = Nm;
    }

    Eigen::VectorXd Err_mrDMD_time = Eigen::VectorXd::Zero(settings.Ns);
//Defining scope for mrDMD
    {

        std::cout << "Extracting basis mrDMD ..." << std::endl;
        std::cout << "Computing nodes for Multi-Resolution Analysis " << std::endl <<std::endl;        

        int level = 0, bin_num = 0, offset = 0, max_levels = settings.max_levels, max_cycles = settings.max_cycles;
        std::vector<node_mrDMD> nodes = {}; 
        nodes = mrDMD_basis( sn_set,      
                            nodes,  
                            settings.r,                    
                            settings.Dt_cfd*settings.Ds,                      
                            t_0,                
                            level,                          
                            bin_num,
                            offset,
                            max_levels,
                            max_cycles,
                            "OPT");

        std::cout << " Done! " << std::endl << std::endl;
        
        int sum = 0;
        for ( int n_nodes = 0; n_nodes < nodes.size(); n_nodes++ )
            sum += nodes[n_nodes].Modes.cols();

        std::cout << "Number of nodes stored : " << nodes.size() << std::endl;
        std::cout << "Number of total modes stored : " << sum << std::endl << std::endl;


//try to compute multi resolution reconstruction in a more clever and straightforward way that DOESN'T WORK!

        // int n_levels = std::floor(std::log((double)settings.Ns/((double)settings.max_cycles*4.0))/std::log(2.0));

        // std::vector<Eigen::MatrixXd> LRec;

        // int count;
        // for ( int i = 0; i < n_levels; i++ )
        // {
        //     count = 0;
        //     Eigen::MatrixXd Rec(ndim*Nr, settings.Ns);
        //     for ( int j = 0; j < nodes.size(); j++ )
        //     {
        //         if ( nodes[j].l == i )
        //         {
        //             while ( count != nodes[j].bin_num )
        //             {
        //                 Rec.middleCols(count*nodes[j].bin_size, nodes[j].bin_size-1) = Eigen::MatrixXd::Zero(ndim*Nr, nodes[j].bin_size-1);         
        //                 count++;
        //             }

        //             if ( nodes[j].bin_num == (std::pow(2,nodes[j].l) - 1) )
        //             {
        //                 Eigen::MatrixXcd temp = nodes[j].Modes*nodes[j].Psi;
        //                 Rec.middleCols(count*nodes[j].bin_size, nodes[j].bin_size) = temp.real();    
        //             }
        //             else
        //             {
        //             Eigen::MatrixXcd temp = nodes[j].Modes*nodes[j].Psi.leftCols(nodes[j].Psi.cols()-1);
        //             Rec.middleCols(count*nodes[j].bin_size, nodes[j].bin_size-1) = temp.real();
        //             count++;
        //             }
        //         }
            
        //     }
        //     LRec.push_back(Rec);
        // }

        // Eigen::MatrixXcd mrDMD_Rec;
        // for ( int i = 0; i < LRec.size(); i++ )
        //     mrDMD_Rec += LRec[i];


        Eigen::MatrixXcd mrDMD_Rec = Eigen::MatrixXcd::Zero(ndim*Nr, settings.Ns);
        
        for ( int i = 0; i < settings.Ns; i++ )
        {
            mrDMD_Rec.col(i) = Reconstruction_mrDMD ( t_vec(i), settings.Dt_cfd*settings.Ds,
                                                    nodes,
                                                    "SCALAR" );
        }

        Eigen::MatrixXd Err_mrDMD_map = sn_set - mrDMD_Rec.real();

        for ( int i = 0; i < settings.Ns; i++ )
        {
            for ( int j = 0; j < ndim*Nr; j++ )
                if ( sn_set(j,i) < tol )
                    Err_mrDMD_time(i) += Err_mrDMD_map(j,i)*Err_mrDMD_map(j,i);
                else    
                    Err_mrDMD_time(i) += Err_mrDMD_map(j,i)*Err_mrDMD_map(j,i); // /(sn_set(j,i)*sn_set(j,i));  

            Err_mrDMD_time(i) = std::sqrt(Err_mrDMD_time(i))/norm_sn_set(i);
        }

        Nm_mrDMD = sum;

    }

    Eigen::VectorXd Err_RDMD_time = Eigen::VectorXd::Zero(settings.Ns);
//Defining scope for RDMD
    {

        Eigen::VectorXd lambda = Eigen::VectorXd::Zero(3*settings.Ns);
        Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(3*settings.Ns, settings.Ns); 
        Eigen::MatrixXd Phi;

        if ( argc == 1 )
        {
            std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        

        
            Phi = RDMD_modes_coefs ( sn_set,
                                    Coefs,
                                    lambda,     
                                    settings.r,
                                    settings.r_RDMD,
                                    settings.En );
        }
        else
        {
            std::cout << "Reading basis and extracting Coeffs RDMD ... " << "\t";
            std::string file_modes = argv[2];
            Phi = read_modes( file_modes, ndim*Nr, settings.r_RDMD );
            Coefs = Phi.transpose()*sn_set;

        }

        std::cout << " Done! " << std::endl << std::endl;

        Eigen::MatrixXd Err_RDMD_map = sn_set - Phi*Coefs;

        for ( int i = 0; i < settings.Ns; i++ )
        {
            for ( int j = 0; j < ndim*Nr; j++ )
                if ( sn_set(j,i) < tol )
                    Err_RDMD_time(i) += Err_RDMD_map(j,i)*Err_RDMD_map(j,i);
                else    
                    Err_RDMD_time(i) += Err_RDMD_map(j,i)*Err_RDMD_map(j,i); // /(sn_set(j,i)*sn_set(j,i));  

            Err_RDMD_time(i) = std::sqrt(Err_RDMD_time(i))/norm_sn_set(i);
        }

        Nm_RDMD = settings.r_RDMD;

    }

    std::cout << "Writing Error file " << "\t";

    std::ofstream error_methods;
    error_methods.open("Error_RBM.dat");

    error_methods << "Time(s)" << "\t";
    error_methods << "Err_POD" << "\t";
    error_methods << "Err_SPOD_" << Nf[0] << "\t";
    error_methods << "Err_SPOD_" << Nf[1] << "\t";
    error_methods << "Err_SPOD_" << Nf[2] << "\t";
    error_methods << "Err_SPOD_" << Nf[3] << "\t";
    error_methods << "Err_DMD" << "\t";
    error_methods << "Err_mrDMD" << "\t";
    error_methods << "Err_RDMD" << "\n";
    
   
    for( int j = 0; j < settings.Ns; j++ ) 
    {
        error_methods << std::setprecision(8) << t_vec(j) << "\t";
        error_methods <<  std::setprecision(8) << Err_POD_time(j) << "\t";
        
        for ( int k = 0; k < Nf.size(); k++ )
            error_methods << std::setprecision(8) << Err_SPOD_time(j,k) << "\t";

        error_methods << std::setprecision(8) << Err_DMD_time(j) << "\t";
        error_methods << std::setprecision(8) << Err_mrDMD_time(j) << "\t";
        error_methods << std::setprecision(8) << Err_RDMD_time(j);

        error_methods << std::endl;

    }
   
    error_methods.close();

    std::cout << "Done" << std::endl;

    int n_err = Nf.size() + 4;
    Eigen::MatrixXd Err_RBM( settings.Ns, n_err );
    Err_RBM.col(0) = Err_POD_time;
    Err_RBM.middleCols(1, Nf.size()) = Err_SPOD_time;
    Err_RBM.col(Nf.size()+1) = Err_DMD_time;
    Err_RBM.col(Nf.size()+2) = Err_mrDMD_time;
    Err_RBM.col(Nf.size()+3) = Err_RDMD_time;
    
    Eigen::VectorXi N_modes_RBM(n_err);

    N_modes_RBM(0) = Nm_POD;
    for ( int i = 0; i < Nf.size(); i ++)
        N_modes_RBM(i+1) = Nm_SPOD[i];
    
    N_modes_RBM(Nf.size()+1) = Nm_DMD;
    N_modes_RBM(Nf.size()+2) = Nm_mrDMD;
    N_modes_RBM(Nf.size()+3) = Nm_RDMD;

    Eigen::VectorXd Err_interp(n_err);
    int index1, index2;

    std::vector<int> pos = {};
    int best_method_idx;

//Adaptive reconstruction on each selected time step
    int Nf_SPOD = 0;
    
    for ( int i = 0; i < settings.t_rec.size(); i++ )
    {
        std::cout << " Adaptive reconstruction at time : " << settings.t_rec[i] << std::endl;

        index1 = 0; 
        index2 = 0;
        for ( int nt = 0; nt < t_vec.size()-1; nt ++ )
        {
            if ( (settings.t_rec[i] > t_vec(nt)) && (settings.t_rec[i] <= t_vec(nt+1)) )
            {
                index1 = nt;
                index2 = nt+1;
                break;
            }
        }

        if ( index1 == index2 )
        {
            std::cout << "Time for reconstruction out of interval!" << std::endl;
            continue;
        }

        int count = 0;
        for ( int k = 0; k < n_err; k ++ )
        {
            Err_interp(k) = Err_RBM(index1,k) + (Err_RBM(index2,k) - Err_RBM(index1,k))/
                            (settings.Dt_cfd*settings.Ds)*(t_vec[index1] - settings.t_rec[i]);
        
            if ( Err_interp(k) < tol_rec )
            {
                pos.push_back(k);
                count ++;
            }

        }

        if ( count == 0 )
        {
            double min_val = Err_interp.minCoeff( &best_method_idx );
            std::cout << " No methods under the prescribed tolerance\n Smallest Error selected " << std::endl;
        }

        if ( count == 1 )
        {
            best_method_idx = pos[0];
            std::cout << " One methods under the prescribed tolerance " << std::endl;
        }

        if ( count > 1 )
        {
            std::cout << count << " methods under the prescribed tolerance\n Selecting the one with less modes " << std::endl;
            Eigen::VectorXi Nm_new(pos.size());
            for ( int j = 0; j < pos.size(); j++ )
                Nm_new(j) = N_modes_RBM(pos[j]);
            
            int temp;
            int min_val = Nm_new.minCoeff( &temp );
            best_method_idx = pos[temp];
            

        }

        std::string method = method_selected ( best_method_idx, Nf_SPOD, Nf );
        std::cout << " Error : " << Err_interp(best_method_idx) << " using method " << method 
                    << " with Nf ( value meaningful only for SPOD ) : " << Nf_SPOD << std::endl;
        
        std::cout << "Computing Reconstruction using selected method " << std::endl;
        
        if ( method == "SPOD" )
        {
            Eigen::VectorXd lambda(settings.Ns);
            Eigen::VectorXd K_pc(settings.Ns);
            Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);

            for ( int kt = 0; kt < settings.Ns; kt++ )
                sn_set.col(kt) -= mean;        

            Eigen::MatrixXd Phi = SPOD_basis( sn_set,
                                    lambda, K_pc, eig_vec,
                                    Nf_SPOD,
                                    settings.flag_bc, 
                                    settings.flag_filter,  
                                    settings.sigma);

            int Nrec = Nmod( settings.En, K_pc);

            std::vector<double> t_v( settings.Ns );
            t_v[0] = settings.nstart*settings.Dt_cfd;

            for ( int kt = 1; kt < settings.Ns; kt++ )
                t_v[kt] = t_v[kt-1] + settings.Dt_cfd*settings.Ds;

            Eigen::MatrixXd Rec = Reconstruction_S_POD ( t_v,
                                K_pc, lambda, eig_vec.transpose(),
                                Phi, settings.t_rec[i],
                                settings.En,
                                settings.flag_prob,
                                settings.flag_interp ) ;

            for ( int kt = 0; kt < Rec.cols(); kt++)
                Rec.col(kt) = Rec.col(kt) + mean.segment(kt*Nr, Nr);

            std::cout << "Writing reconstructed field ..." << "\t";

            write_Reconstructed_fields ( Rec, Coords,
                                    settings.out_file,
                                    settings.flag_prob, i );

            std::cout << "Done" << std::endl << std::endl;
            
        }

        if ( method == "DMD" )
        {

        Eigen::VectorXd lambda_POD;
        Eigen::MatrixXd eig_vec_POD;
        Eigen::VectorXcd lambda_DMD;
        Eigen::MatrixXcd eig_vec_DMD;      
        Eigen::MatrixXcd Phi;
        Eigen::VectorXcd alfa;    

        Phi = DMD_basis( sn_set,
                        lambda_DMD,
                        eig_vec_DMD,
                        lambda_POD,
                        eig_vec_POD,
                        0 );

        alfa = Calculate_Coefs_DMD_exact ( sn_set.leftCols(settings.Ns-1),  
                                                            lambda_DMD,  
                                                            Phi );                    
                
        Eigen::MatrixXcd Rec = Reconstruction_DMD ( settings.t_rec[i],
                                                settings.Dt_cfd*settings.Ds,
                                                alfa,
                                                Phi,
                                                lambda_DMD,
                                                settings.flag_prob );

        std::cout << "Writing reconstructed field ..." << "\t";

        write_Reconstructed_fields ( Rec.real(), Coords,
                                settings.out_file,
                                settings.flag_prob, i );

        std::cout << "Done" << std::endl << std::endl;
        
        }


        if ( method == "mrDMD" )
        {

            double dts = settings.Dt_cfd*settings.Ds;

            int level = 0, bin_num = 0, offset = 0, max_levels = settings.max_levels, max_cycles = settings.max_cycles;
            std::vector<node_mrDMD> nodes = {}; 

            nodes = mrDMD_basis( sn_set,      
                                nodes,  
                                -1,                    
                                dts,                      
                                t_0,                
                                level,                          
                                bin_num,
                                offset,
                                max_levels,
                                max_cycles,
                                "OPT");

            Eigen::MatrixXcd Rec = Reconstruction_mrDMD ( settings.t_rec[i],                                                                                                                                                                                                                                                                                                                              
                                                        dts,       
                                                        nodes,     
                                                        settings.flag_prob );  

            std::cout << "Writing reconstructed field ..." << "\t";

            write_Reconstructed_fields ( Rec.real(), Coords,
                                    settings.out_file,
                                    settings.flag_prob, i );

            std::cout << "Done" << std::endl << std::endl;

        }

        if ( method == "RDMD" )
        {

            // for ( int i = 0; i < settings.Ns; i++ )
                // sn_set.col(i) -= mean;        
            Eigen::VectorXd lambda = Eigen::VectorXd::Zero(3*settings.Ns);
            Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(3*settings.Ns, settings.Ns);
            Eigen::MatrixXd Phi;

            if ( argc == 1 )
            {
                std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        
                //You can define rank DMD at each time step from the config file ( use -1 for the adaptive study adviced)
                Phi = RDMD_modes_coefs ( sn_set,
                                        Coefs,
                                        lambda,     
                                        settings.r,
                                        settings.r_RDMD,
                                        settings.En );
            }
            else
            {

                std::string file_modes = argv[2];
                Phi = read_modes( file_modes, ndim*Nr, settings.r_RDMD );
                Coefs = Phi.transpose()*sn_set;

            }

            std::vector<double> t_st_vec(settings.Ns);
            t_st_vec[0] = t_0;

            for ( int i = 1; i < settings.Ns; i++ )
                t_st_vec[i] = t_st_vec[i-1] + settings.Dt_cfd*settings.Ds;

            Eigen::MatrixXd Rec = Reconstruction_RDMD ( settings.t_rec[i],
                                    t_st_vec,
                                    Coefs,
                                    Phi,
                                    settings.flag_prob,
                                    settings.flag_interp );

                        // for ( int i = 0; i < Rec.cols(); i++)
                            // Rec.col(i) = Rec.col(i) + mean.segment(i*Nr, Nr);

            std::cout << "Writing reconstructed field ..." << "\t";

            write_Reconstructed_fields ( Rec, Coords,
                                    settings.out_file,
                                    settings.flag_prob, i );

            std::cout << "Done" << std::endl << std::endl;            

        }


    }

    return 0;

}