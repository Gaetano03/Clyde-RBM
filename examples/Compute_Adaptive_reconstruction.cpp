#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"


int main( int argc, char *argv[] )
{
    std::cout << "Adaptive Reconstruction RBM-Clyde starts " << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];
    
    //Reading configuration file
    Read_cfg( filecfg, settings );
    double t_0 = settings.nstart*settings.Dt_cfd;

    int s_Nf = 5;
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
    Nf[1] = std::ceil(settings.Ns/10.0);
    Nf[2] = std::ceil(settings.Ns/2.0);
    Nf[3] = std::ceil(2.0*settings.Ns/3.0);
    Nf[4] = settings.Ns;
    int Nf_SPOD = 0;

    int Nr = N_gridpoints ( settings.in_file );
    std::cout << "Number of grid points : " << Nr << std::endl;

    std::cout << "Reading Coordinates ... \t ";
    Eigen::MatrixXd Coords = read_col( settings.in_file, Nr, settings.Cols_coords );
    std::cout << "Done " << std::endl;

    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings.Ns, settings.Ds, settings.nstart,
                                        settings.Cols,
                                        settings.in_file,
                                        settings.flag_prob);

    Eigen::VectorXd mean = sn_set.rowwise().mean();
    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(settings.Ns);

    for ( int nt = 0; nt < settings.Ns; nt++ )
        sn_set.col(nt) -= mean;

    std::cout << "Initializing error surface ... " << std::endl;
    std::vector<Eigen::MatrixXd> Err_RBM_Nm_time;
    std::vector<Eigen::MatrixXd> J_RBM_Nm_time;
    std::vector<Eigen::MatrixXd> EJ_RBM_Nm_time;

    std::string filename1, filename2;
    for ( int i = 0; i < s_Nf; i++ )
    {
        filename1 = "Error_SPOD_Nf_" + std::to_string(i) + ".dat";
        Eigen::MatrixXd Err_map = read_err_j ( filename1, settings.Ns );

        filename2 = "Jaccard_SPOD_Nf_" + std::to_string(i) + ".dat";
        Eigen::MatrixXd J_map = read_err_j ( filename2, settings.Ns );

        Err_RBM_Nm_time.push_back(Err_map);
        J_RBM_Nm_time.push_back(J_map);
        EJ_RBM_Nm_time.push_back(Err_map.cwiseQuotient(J_map));

    }

    filename1 = "Error_DMD.dat";
    Eigen::MatrixXd Err_map = read_err_j ( filename1, settings.Ns );

    filename2 = "Jaccard_DMD.dat";
    Eigen::MatrixXd J_map = read_err_j ( filename2, settings.Ns );

    Err_RBM_Nm_time.push_back(Err_map);
    J_RBM_Nm_time.push_back(J_map);
    EJ_RBM_Nm_time.push_back(Err_map.cwiseQuotient(J_map));

    filename1 = "Error_RDMD.dat";
    Err_map = read_err_j ( filename1, settings.Ns );

    filename2 = "Jaccard_RDMD.dat";
    J_map = read_err_j ( filename2, settings.Ns );

    Err_RBM_Nm_time.push_back(Err_map);
    J_RBM_Nm_time.push_back(J_map);
    EJ_RBM_Nm_time.push_back(Err_map.cwiseQuotient(J_map));      

//Adaptive reconstruction on each selected time step
    int best_method_idx;

    std::cout << "Initializing Vector of time ... " << std::endl; 
    Eigen::VectorXd t_vec( settings.Ns );
    t_vec(0) = settings.nstart*settings.Dt_cfd;
    for ( int i = 1; i < settings.Ns; i++ )
        t_vec(i) = t_vec(i-1) + settings.Dt_cfd*settings.Ds;
    
    double tol = settings.tol;
    int index1, index2;
    Eigen::MatrixXd Err_interp(EJ_RBM_Nm_time.size(), settings.Ns);

    for ( int i = 0; i < settings.t_rec.size(); i++ )
    {
        std::vector<int> pos = {};
        std::cout << " Adaptive reconstruction at time : " << settings.t_rec[i] << std::endl;

        index1 = 0;
        index2 = 0;
        for ( int nt = 0; nt < t_vec.size()-1; nt ++ )
        {
            if ( (settings.t_rec[i] >= t_vec(nt)) && (settings.t_rec[i] <= t_vec(nt+1)) )
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
        for ( int k = 0; k < EJ_RBM_Nm_time.size(); k ++ )
        {
            for ( int nm = 0; nm < settings.Ns; nm ++)
                Err_interp(k,nm) = EJ_RBM_Nm_time[k](index1,nm) + (EJ_RBM_Nm_time[k](index2,nm) - EJ_RBM_Nm_time[k](index1,nm))/
                                    (settings.Dt_cfd*settings.Ds)*(settings.t_rec[i] - t_vec[index1]);
        }

        double eps = 1.0;
        int ncount = 0;

        while ( eps > tol )
        {
            eps = Err_interp.col(ncount).minCoeff( &best_method_idx );
            ncount++;
        }

        std::string method = method_selected ( best_method_idx, Nf_SPOD, Nf );
        std::cout << "Best method is " << method << " with number of modes " << ncount - 1
            << " and Nf ( value meaningful only for SPOD ) : " << Nf_SPOD << std::endl;
        
        std::cout << " Error : " << Err_interp(best_method_idx, ncount-1) << std::endl;
                        
            
        std::cout << "Computing Reconstruction using selected method " << std::endl;
        
        if ( method == "SPOD" )
        {
            Eigen::VectorXd lambda(settings.Ns);
            Eigen::VectorXd K_pc(settings.Ns);
            Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);        

            Eigen::MatrixXd Phi = SPOD_basis( sn_set,
                                    lambda, K_pc, eig_vec,
                                    Nf_SPOD,
                                    settings.flag_bc, 
                                    settings.flag_filter,  
                                    settings.sigma);

            int Nrec = ncount - 1;

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
                            settings.r );

            alfa = Calculate_Coefs_DMD_exact ( sn_set.leftCols(settings.Ns-1),  
                                                                lambda_DMD,  
                                                                Phi );                    
                    
            Eigen::MatrixXcd Rec = Reconstruction_DMD ( settings.t_rec[i],
                                                    settings.Dt_cfd*settings.Ds,
                                                    alfa,
                                                    Phi,
                                                    lambda_DMD,
                                                    settings.flag_prob );

            for ( int kt = 0; kt < Rec.cols(); kt++)
                Rec.col(kt) = Rec.col(kt) + mean.segment(kt*Nr, Nr);

            std::cout << "Writing reconstructed field ..." << "\t";

            write_Reconstructed_fields ( Rec.real(), Coords,
                                    settings.out_file,
                                    settings.flag_prob, i );

            std::cout << "Done" << std::endl << std::endl;
        
        }


        if ( method == "RDMD" )
        {
        
            Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
            Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
            Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
            Eigen::MatrixXd Phi;

            if ( argc == 2 )
            {
                std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        
                //You can define rank DMD at each time step from the config file ( use -1 for the adaptive study adviced)
                Phi = RDMD_modes_coefs ( sn_set,
                                        Coefs,
                                        lambda,
                                        K_pc,     
                                        settings.r,
                                        settings.r_RDMD,
                                        settings.En );
            }
            else
            {
                std::cout << "Reading basis and extracting Coeffs RDMD ... " << "\t"; 
                std::string file_modes = argv[2];
                std::string file_coefs = argv[3];
                Phi = read_modes( file_modes, settings.ndim*Nr, settings.r_RDMD );
                Coefs = read_coefs( file_coefs, settings.Ns, settings.r_RDMD );

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

            for ( int kt = 0; kt < Rec.cols(); kt++)
                Rec.col(kt) = Rec.col(kt) + mean.segment(kt*Nr, Nr);

            std::cout << "Writing reconstructed field ..." << "\t";

            write_Reconstructed_fields ( Rec, Coords,
                                    settings.out_file,
                                    settings.flag_prob, i );

            std::cout << "Done" << std::endl << std::endl << std::endl;            

        }

    }

    std::cout << "Adaptive Reconstruction RBM-Clyde ends " << std::endl;

    return 0;
}






