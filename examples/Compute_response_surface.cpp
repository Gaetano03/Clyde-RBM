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

    std::cout << "Initializing common variables " << std::endl << std::endl;

    std::string filecfg = argv[1];
    // std::string mode = argv[2];
    double gamma = 1.4;     //ratio of specific heats;
    double R = 287.058;     //gas constant;
    double Mach = 0.2;      //Mach number;
    double T = 300;         //Temperature;

    double Vinf = Mach*std::sqrt(gamma*R*T);

    prob_settings settings;

    //Reading configuration file
    Read_cfg( filecfg, settings );
    int ndim = settings.ndim;   //dimension of the problem (1D/2D/3D)
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

    double tol_rec = settings.tol;
    double tol = 1e-1;
    double t_0 = (double)settings.nstart*settings.Dt_cfd;
    Eigen::VectorXd t_vec( settings.Ns );
    t_vec(0) = t_0;

    for ( int i = 1; i < settings.Ns; i++ )
        t_vec(i) = t_vec(i-1) + settings.Dt_cfd*settings.Ds;

    std::cout << std::endl;
    std::cout << "Initialized vector of times " << std::endl << std::endl;

    Eigen::VectorXd mean = sn_set.rowwise().mean();
    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(settings.Ns);

    for ( int i = 0; i < settings.Ns; i ++ )
    {
        norm_sn_set(i) = sn_set.col(i).norm();
    }

    std::vector<int> Nf(5);
    Nf[0] = 0;
    Nf[1] = std::ceil(settings.Ns/10.0);
    Nf[2] = std::ceil(settings.Ns/2.0);
    Nf[3] = std::ceil(2.0*settings.Ns/3.0);
    Nf[4] = settings.Ns;
    std::vector<Eigen::MatrixXd> Err_SPOD_Nf_Nm_time;
    std::vector<Eigen::MatrixXd> J_SPOD_Nf_Nm_time;

//Defining common scope for POD-SPOD
    {
        Eigen::VectorXd lambda(settings.Ns);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);

        int Nrec = 0;

        for ( int nt = 0; nt < settings.Ns; nt++ )
            sn_set.col(nt) -= mean;

        for ( int nfj = 0; nfj < Nf.size(); nfj++ )
        {

            std::cout << "Extracting  SPOD " << Nf[nfj] << " basis ... " << "\t";        

            Eigen::MatrixXd Phi = SPOD_basis( sn_set,
                                    lambda, K_pc, eig_vec,
                                    Nf[nfj],
                                    settings.flag_bc, 
                                    settings.flag_filter,  
                                    settings.sigma);

            std::cout << " Done! " << std::endl;
            
            Nrec = Phi.cols();
            Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nrec, Nrec);

            // for ( int i = 0; i < settings.Ns; i++ )
            //     sn_set.col(i) += mean;

            for ( int i = 0; i < Nrec; i++ )
                Sig(i,i) = std::sqrt(lambda(i));

            std::cout << "Computing error and Jaccard index surface..." << "\t";

            Eigen::MatrixXd Err_SPOD_map( sn_set.rows(), sn_set.cols() );
            Eigen::MatrixXd Err_SPOD_Nm_time( settings.Ns, Phi.cols() );
            Eigen::MatrixXd J_SPOD_Nm_time( settings.Ns, Phi.cols() );

            for ( int nm = 1; nm <= Nrec; nm ++ )
            {
                Err_SPOD_map = sn_set - Phi.leftCols(nm)*Sig.block(0,0,nm,nm)*eig_vec.transpose().topRows(nm);

                for ( int i = 0; i < settings.Ns; i++ )
                {

                    int count = 0;
                    for ( int j = 0; j < ndim*Nr; j++ )
                    {
                        Err_SPOD_Nm_time(i,nm-1) += Err_SPOD_map(j,i)*Err_SPOD_map(j,i);
                    
                        if ( Err_SPOD_map(j,i) < 0.1)
                            count++;
                    }
                    
                    Err_SPOD_Nm_time(i,nm-1) = std::sqrt(Err_SPOD_Nm_time(i,nm-1))/norm_sn_set(i);
                    J_SPOD_Nm_time(i,nm-1) = (double)count/((double)ndim*(double)Nr);
                
                }
            
            }

            Err_SPOD_Nf_Nm_time.push_back(Err_SPOD_Nm_time);
            J_SPOD_Nf_Nm_time.push_back(J_SPOD_Nm_time);

            std::cout << "Done" << std::endl;
        }

    }

//Adding again mean over snapshots to perform DMD and its variants
    for ( int i = 0; i < settings.Ns; i++ )
    sn_set.col(i) += mean;

    Eigen::MatrixXd Err_DMD_Nm_time;
    Eigen::MatrixXd J_DMD_Nm_time;

//Defining scope for DMD ( Rank=-1 preferable, Coeffs = OPT )
    {

        Eigen::VectorXd lambda_POD;
        Eigen::MatrixXd eig_vec_POD;
        Eigen::VectorXcd lambda_DMD;
        Eigen::MatrixXcd eig_vec_DMD;

        std::cout << "Extracting basis DMD using rank " << settings.r << "\t";        
        Eigen::MatrixXcd Phi;
        Eigen::VectorXcd alfa;

        Phi = DMD_basis( sn_set,
                        lambda_DMD,
                        eig_vec_DMD,
                        lambda_POD,
                        eig_vec_POD,
                        settings.r );

        std::cout << " Done! " << std::endl;
        int Nm = Phi.cols();

        Eigen::VectorXcd omega(Nm);
        for ( int i = 0; i < Nm; i++ )
            omega(i) = std::log(lambda_DMD(i))/(settings.Dt_cfd*settings.Ds);

        std::cout << "Calculating coefficients DMD ... " << "\t";
            
        alfa = Calculate_Coefs_DMD_exact ( sn_set.leftCols(settings.Ns-1),  
                                                            lambda_DMD, 
                                                            Phi );
            
        std::cout << " Done! " << std::endl;

        std::cout << "Reordering modes DMD ... " << "\t";
        Eigen::VectorXd En = Eigen::VectorXd::Zero(Phi.cols());
        double T = t_vec(t_vec.size()-1);

        for ( int i = 0 ; i < Phi.cols(); i ++ )
        {

            double alfa_i = alfa(i).imag();
            double alfa_r = alfa(i).real();
            double sigma = omega(i).real();

            En(i) = (alfa_r*alfa_r + alfa_i*alfa_i)*(std::exp(2.0*sigma*T) - 1.0)/(2.0*sigma);

            std::cout << "En at mode " << i << " " <<  En(i) << std::endl;

        }

        dmd_sort( En, Phi, lambda_DMD, alfa);

        std::cout << "Done" << std::endl;

        std::cout << "Computing error and Jaccard index surface ... " << std::endl << std::endl;

        Eigen::MatrixXcd V_and(lambda_DMD.size(), settings.Ns);

      
        for ( int i = 0; i < lambda_DMD.size(); i++ )
        {
            for ( int j = 0; j < settings.Ns; j++ )
                V_and(i,j) = std::pow(lambda_DMD(i), (double)j);                                                                                         
        }        

        Eigen::MatrixXcd Psi = Eigen::MatrixXcd::Zero(alfa.size(), settings.Ns);
        for ( int i = 0; i < settings.Ns; i++ )
            Psi.col(i) = alfa.cwiseProduct(V_and.col(i));
  
        Eigen::MatrixXd Err_DMD_map(sn_set.rows(), sn_set.cols());
        Err_DMD_Nm_time = Eigen::MatrixXd::Zero(settings.Ns, Phi.cols());
        J_DMD_Nm_time = Eigen::MatrixXd::Zero(settings.Ns, Phi.cols());

        for ( int nm = 1; nm <= Nm; nm ++ )
        {

            Eigen::MatrixXcd D_dmd = Phi.leftCols(nm)*Psi.topRows(nm);
            Err_DMD_map = sn_set - D_dmd.real();

            for ( int i = 0; i < settings.Ns; i++ )
            {
                int count = 0;
                for ( int j = 0; j < ndim*Nr; j++ )
                { 
                    Err_DMD_Nm_time(i,nm-1) += Err_DMD_map(j,i)*Err_DMD_map(j,i);
                    if ( Err_DMD_map(j,i) < 0.1 )
                        count++;
                }

                Err_DMD_Nm_time(i,nm-1) = std::sqrt(Err_DMD_Nm_time(i,nm-1))/norm_sn_set(i);
                J_DMD_Nm_time(i,nm-1) = (double)count/((double)ndim*(double)Nr);
            }

        }
    }


    Eigen::MatrixXd Err_RDMD_Nm_time;
    Eigen::MatrixXd J_RDMD_Nm_time;
//Defining scope for RDMD
    {

        Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns); 
        Eigen::MatrixXd Phi;

        if ( argc == 2 )
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

        std::cout << " Done! " << std::endl;

        std::cout << "Computing error and Jaccard index surface ... " << std::endl << std::endl;
        Eigen::MatrixXd Err_RDMD_map;
        Err_RDMD_Nm_time = Eigen::MatrixXd::Zero(settings.Ns, Phi.cols());
        J_RDMD_Nm_time = Eigen::MatrixXd::Zero(settings.Ns, Phi.cols());

        for ( int nm = 1; nm <= Phi.cols(); nm ++ )
        {

            Err_RDMD_map = sn_set - Phi.leftCols(nm)*Coefs.topRows(nm);
            for ( int i = 0; i < settings.Ns; i++ )
            {
                int count = 0;

                for ( int j = 0; j < ndim*Nr; j++ )
                {
                    Err_RDMD_Nm_time(i,nm-1) += Err_RDMD_map(j,i)*Err_RDMD_map(j,i);
                    if ( Err_RDMD_map(j,i) < 0.1 )
                        count++;
                }

                Err_RDMD_Nm_time(i,nm-1) = std::sqrt(Err_RDMD_Nm_time(i,nm-1))/norm_sn_set(i);
                J_RDMD_Nm_time(i,nm-1) = (double)count/((double)ndim*(double)Nr);
            }

        }

    }

    std::cout << "Writing surface error files ... " << "\t";
    // -------- Nm ---------
    // |                    |
    // |                    |
    // t         E          |
    // |                    |
    // |                    |
    //----------------------

    for ( int i = 0; i < Nf.size(); i++ )
    {
    
        std::string filename = "Error_SPOD_Nf_" + std::to_string(i) + ".dat";
        std::ofstream error_SPOD;
        error_SPOD.open(filename);
    
        for( int j = 0; j < settings.Ns; j++ ) 
        {
            for ( int nm = 0; nm < Err_SPOD_Nf_Nm_time[i].cols(); nm ++ )
                error_SPOD <<  std::setprecision(8) << Err_SPOD_Nf_Nm_time[i](j,nm) << "\t";

            error_SPOD << std::endl;

        }
    
        error_SPOD.close();

    }

    std::ofstream error_DMD;
    error_DMD.open("Error_DMD.dat");

    for( int j = 0; j < settings.Ns; j++ ) 
    {
        for ( int nm = 0; nm < Err_DMD_Nm_time.cols(); nm ++ )
            error_DMD <<  std::setprecision(8) << Err_DMD_Nm_time(j,nm) << "\t";

        error_DMD << std::endl;

    }

    error_DMD.close();

    std::ofstream error_RDMD;
    error_RDMD.open("Error_RDMD.dat");

    for( int j = 0; j < settings.Ns; j++ ) 
    {
        for ( int nm = 0; nm < Err_RDMD_Nm_time.cols(); nm ++ )
            error_RDMD <<  std::setprecision(8) << Err_RDMD_Nm_time(j,nm) << "\t";

        error_RDMD << std::endl;

    }

    error_RDMD.close();
    std::cout << "Done" << std::endl << std::endl;

    std::cout << "Writing Jaccard index files ... " << "\t";
    //--------- Nm ----------
    // |                    |
    // |                    |
    // t          J         |
    // |                    |
    // |                    |
    //-----------------------

    for ( int i = 0; i < Nf.size(); i++ )
    {
    
        std::string filename = "Jaccard_SPOD_Nf_" + std::to_string(i) + ".dat";
        std::ofstream Jaccard_SPOD;
        Jaccard_SPOD.open(filename);
    
        for( int j = 0; j < settings.Ns; j++ ) 
        {
            for ( int nm = 0; nm < J_SPOD_Nf_Nm_time[i].cols(); nm ++ )
                Jaccard_SPOD <<  std::setprecision(8) << J_SPOD_Nf_Nm_time[i](j,nm) << "\t";

            Jaccard_SPOD << std::endl;

        }
    
        Jaccard_SPOD.close();

    }

    std::ofstream Jaccard_DMD;
    Jaccard_DMD.open("Jaccard_DMD.dat");

    for( int j = 0; j < settings.Ns; j++ ) 
    {
        for ( int nm = 0; nm < J_DMD_Nm_time.cols(); nm ++ )
            Jaccard_DMD <<  std::setprecision(8) << J_DMD_Nm_time(j,nm) << "\t";

        Jaccard_DMD << std::endl;

    }

    Jaccard_DMD.close();

    std::ofstream Jaccard_RDMD;
    Jaccard_RDMD.open("Jaccard_RDMD.dat");

    for( int j = 0; j < settings.Ns; j++ ) 
    {
        for ( int nm = 0; nm < J_RDMD_Nm_time.cols(); nm ++ )
            Jaccard_RDMD <<  std::setprecision(8) << J_RDMD_Nm_time(j,nm) << "\t";

        Jaccard_RDMD << std::endl;

    }

    Jaccard_RDMD.close();
    std::cout << "Done" << std::endl << std::endl;



    // if ( settings.flag_rec == "YES")
    // {

    //     int n_err = Nf.size() + 4;
    //     Eigen::MatrixXd Err_RBM( settings.Ns, n_err );
    //     Err_RBM.col(0) = Err_POD_time;
    //     Err_RBM.middleCols(1, Nf.size()) = Err_SPOD_time;
    //     Err_RBM.col(Nf.size()+1) = Err_DMD_time;
    //     Err_RBM.col(Nf.size()+2) = Err_mrDMD_time;
    //     Err_RBM.col(Nf.size()+3) = Err_RDMD_time;
        
    //     Eigen::VectorXi N_modes_RBM(n_err);

    //     N_modes_RBM(0) = Nm_POD;
    //     for ( int i = 0; i < Nf.size(); i ++)
    //         N_modes_RBM(i+1) = Nm_SPOD[i];
        
    //     N_modes_RBM(Nf.size()+1) = Nm_DMD;
    //     N_modes_RBM(Nf.size()+2) = Nm_mrDMD;
    //     N_modes_RBM(Nf.size()+3) = Nm_RDMD;

    //     Eigen::VectorXd Err_interp(n_err);
    //     int index1, index2;


    //     int best_method_idx;

    // //Adaptive reconstruction on each selected time step
    //     int Nf_SPOD = 0;
    
    //     for ( int i = 0; i < settings.t_rec.size(); i++ )
    //     {
    //         std::vector<int> pos = {};
    //         std::cout << " Adaptive reconstruction at time : " << settings.t_rec[i] << std::endl;

    //         index1 = 0;
    //         index2 = 0;
    //         for ( int nt = 0; nt < t_vec.size()-1; nt ++ )
    //         {
    //             if ( (settings.t_rec[i] > t_vec(nt)) && (settings.t_rec[i] <= t_vec(nt+1)) )
    //             {
    //                 index1 = nt;
    //                 index2 = nt+1;
    //                 break;
    //             }
    //         }

    //         if ( index1 == index2 )
    //         {
    //             std::cout << "Time for reconstruction out of interval!" << std::endl;
    //             continue;
    //         }

    //         int count = 0;
    //         for ( int k = 0; k < n_err; k ++ )
    //         {
    //             Err_interp(k) = Err_RBM(index1,k) + (Err_RBM(index2,k) - Err_RBM(index1,k))/
    //                             (settings.Dt_cfd*settings.Ds)*(settings.t_rec[i] - t_vec[index1]);
            
    //             if ( Err_interp(k) < tol_rec )
    //             {
    //                 pos.push_back(k);
    //                 count ++;
    //             }

    //         }

    // // std::cout << "Err_interp for current time step : " << Err_interp << std::endl;
    // // std::cout << " positions for best methods (under tolerance " <<  settings.tol << ")" << std::endl;
    // // for ( int kk = 0; kk < pos.size(); kk++)
    // //     std::cout << pos[kk] << "\t";

    // //     std::cout << std::endl;

    //         if ( count == 0 )
    //         {
    //             double min_val = Err_interp.minCoeff( &best_method_idx );
    //             std::cout << " No methods under the prescribed tolerance\n Smallest Error selected " << std::endl;
    //         }

    //         if ( count == 1 )
    //         {
    //             best_method_idx = pos[0];
    //             std::cout << " One methods under the prescribed tolerance " << std::endl;
    //         }

    //         if ( count > 1 )
    //         {
    //             std::cout << count << " methods under the prescribed tolerance\n Selecting the one with less modes " << std::endl;
    //             Eigen::VectorXi Nm_new(pos.size());
    //             for ( int j = 0; j < pos.size(); j++ )
    //                 Nm_new(j) = N_modes_RBM(pos[j]);
                
    //             int temp;
    //             int min_val = Nm_new.minCoeff( &temp );
    //             best_method_idx = pos[temp];
                

    //         }

    //         std::string method = method_selected ( best_method_idx, Nf_SPOD, Nf );
    //         std::cout << " Error : " << Err_interp(best_method_idx) << " using method " << method 
    //                     << " with Nf ( value meaningful only for SPOD ) : " << Nf_SPOD << std::endl;
            
    //         std::cout << "Computing Reconstruction using selected method " << std::endl;
            
    //         if ( method == "SPOD" )
    //         {
    //             Eigen::VectorXd lambda(settings.Ns);
    //             Eigen::VectorXd K_pc(settings.Ns);
    //             Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);

    //             for ( int kt = 0; kt < settings.Ns; kt++ )
    //                 sn_set.col(kt) -= mean;        

    //             Eigen::MatrixXd Phi = SPOD_basis( sn_set,
    //                                     lambda, K_pc, eig_vec,
    //                                     Nf_SPOD,
    //                                     settings.flag_bc, 
    //                                     settings.flag_filter,  
    //                                     settings.sigma);

    //             int Nrec = Nmod( settings.En, K_pc);

    //             std::vector<double> t_v( settings.Ns );
    //             t_v[0] = settings.nstart*settings.Dt_cfd;

    //             for ( int kt = 1; kt < settings.Ns; kt++ )
    //                 t_v[kt] = t_v[kt-1] + settings.Dt_cfd*settings.Ds;

    //             Eigen::MatrixXd Rec = Reconstruction_S_POD ( t_v,
    //                                 K_pc, lambda, eig_vec.transpose(),
    //                                 Phi, settings.t_rec[i],
    //                                 settings.En,
    //                                 settings.flag_prob,
    //                                 settings.flag_interp ) ;

    //             for ( int kt = 0; kt < Rec.cols(); kt++)
    //                 Rec.col(kt) = Rec.col(kt) + mean.segment(kt*Nr, Nr);

    //             std::cout << "Writing reconstructed field ..." << "\t";

    //             write_Reconstructed_fields ( Rec, Coords,
    //                                     settings.out_file,
    //                                     settings.flag_prob, i );

    //             std::cout << "Done" << std::endl << std::endl;
                
    //         }

    //         if ( method == "DMD" )
    //         {

    //         Eigen::VectorXd lambda_POD;
    //         Eigen::MatrixXd eig_vec_POD;
    //         Eigen::VectorXcd lambda_DMD;
    //         Eigen::MatrixXcd eig_vec_DMD;      
    //         Eigen::MatrixXcd Phi;
    //         Eigen::VectorXcd alfa;    

    //         Phi = DMD_basis( sn_set,
    //                         lambda_DMD,
    //                         eig_vec_DMD,
    //                         lambda_POD,
    //                         eig_vec_POD,
    //                         0 );

    //         alfa = Calculate_Coefs_DMD_exact ( sn_set.leftCols(settings.Ns-1),  
    //                                                             lambda_DMD,  
    //                                                             Phi );                    
                    
    //         Eigen::MatrixXcd Rec = Reconstruction_DMD ( settings.t_rec[i],
    //                                                 settings.Dt_cfd*settings.Ds,
    //                                                 alfa,
    //                                                 Phi,
    //                                                 lambda_DMD,
    //                                                 settings.flag_prob );

    //         std::cout << "Writing reconstructed field ..." << "\t";

    //         write_Reconstructed_fields ( Rec.real(), Coords,
    //                                 settings.out_file,
    //                                 settings.flag_prob, i );

    //         std::cout << "Done" << std::endl << std::endl;
            
    //         }


    //         if ( method == "mrDMD" )
    //         {

    //             double dts = settings.Dt_cfd*settings.Ds;

    //             int level = 0, bin_num = 0, offset = 0, max_levels = settings.max_levels, max_cycles = settings.max_cycles;
    //             std::vector<node_mrDMD> nodes = {}; 

    //             nodes = mrDMD_basis( sn_set,      
    //                                 nodes,  
    //                                 -1,                    
    //                                 dts,                      
    //                                 t_0,                
    //                                 level,                          
    //                                 bin_num,
    //                                 offset,
    //                                 max_levels,
    //                                 max_cycles,
    //                                 "OPT");

    //             Eigen::MatrixXcd Rec = Reconstruction_mrDMD ( settings.t_rec[i],                                                                                                                                                                                                                                                                                                                              
    //                                                         dts,       
    //                                                         nodes,     
    //                                                         settings.flag_prob );  

    //             std::cout << "Writing reconstructed field ..." << "\t";

    //             write_Reconstructed_fields ( Rec.real(), Coords,
    //                                     settings.out_file,
    //                                     settings.flag_prob, i );

    //             std::cout << "Done" << std::endl << std::endl;

    //         }

    //         if ( method == "RDMD" )
    //         {

    //             // for ( int i = 0; i < settings.Ns; i++ )
    //                 // sn_set.col(i) -= mean;        
    //             Eigen::VectorXd lambda = Eigen::VectorXd::Zero(3*settings.Ns);
    //             Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(3*settings.Ns, settings.Ns);
    //             Eigen::MatrixXd Phi;

    //             if ( argc == 2 )
    //             {
    //                 std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        
    //                 //You can define rank DMD at each time step from the config file ( use -1 for the adaptive study adviced)
    //                 Phi = RDMD_modes_coefs ( sn_set,
    //                                         Coefs,
    //                                         lambda,     
    //                                         settings.r,
    //                                         settings.r_RDMD,
    //                                         settings.En );
    //             }
    //             else
    //             {
    //                 std::cout << "Reading basis and extracting Coeffs RDMD ... " << "\t"; 
    //                 std::string file_modes = argv[2];
    //                 Phi = read_modes( file_modes, ndim*Nr, settings.r_RDMD );
    //                 Coefs = Phi.transpose()*sn_set;

    //             }

    //             std::vector<double> t_st_vec(settings.Ns);
    //             t_st_vec[0] = t_0;

    //             for ( int i = 1; i < settings.Ns; i++ )
    //                 t_st_vec[i] = t_st_vec[i-1] + settings.Dt_cfd*settings.Ds;

    //             Eigen::MatrixXd Rec = Reconstruction_RDMD ( settings.t_rec[i],
    //                                     t_st_vec,
    //                                     Coefs,
    //                                     Phi,
    //                                     settings.flag_prob,
    //                                     settings.flag_interp );

    //                         // for ( int i = 0; i < Rec.cols(); i++)
    //                             // Rec.col(i) = Rec.col(i) + mean.segment(i*Nr, Nr);

    //             std::cout << "Writing reconstructed field ..." << "\t";

    //             write_Reconstructed_fields ( Rec, Coords,
    //                                     settings.out_file,
    //                                     settings.flag_prob, i );

    //             std::cout << "Done" << std::endl << std::endl;            

    //         }


    //     }
    // }


    return 0;

}