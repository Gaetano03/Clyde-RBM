#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"


int main( int argc, char *argv[] )
{

    std::cout << "-----------RBM-Clyde error surface calculation starts-------------" << std::endl << std::endl;
    
    std::cout << "Initializing common variables ... " << std::endl << std::endl;
    prob_settings settings;
    std::string filecfg = argv[1];
    Read_cfg( filecfg, settings );
    int s_Nf = 5;   //Number of values for the SPOD filter
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
    Nf[1] = std::ceil(settings.Ns/10.0);
    Nf[2] = std::ceil(settings.Ns/2.0);
    Nf[3] = std::ceil(2.0*settings.Ns/3.0);
    Nf[4] = settings.Ns;
    std::vector<Eigen::MatrixXd> Err_RBM_Nm_time;
    std::vector<Eigen::MatrixXd> J_RBM_Nm_time;

    // Calculate number of grid points
    int Nr = N_gridpoints ( settings.in_file );
    std::cout << "Number of grid points : " << Nr << std::endl;

    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings.Ns, settings.Ds, settings.nstart,
                                        settings.Cols,
                                        settings.in_file,
                                        settings.flag_prob);

    std::cout << "Initializing Vector of time ... " << std::endl; 
    Eigen::VectorXd t_vec( settings.Ns );
    t_vec(0) = settings.nstart*settings.Dt_cfd;
    for ( int i = 1; i < settings.Ns; i++ )
        t_vec(i) = t_vec(i-1) + settings.Dt_cfd*settings.Ds;

    std::cout << "Computing mean and norm of CFD solution ... " << std::endl;
    Eigen::VectorXd mean = sn_set.rowwise().mean();
    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(settings.Ns);
    for ( int i = 0; i < settings.Ns; i ++ )
    {
        norm_sn_set(i) = sn_set.col(i).norm();
    }

    for ( int nt = 0; nt < settings.Ns; nt++ )
        sn_set.col(nt) -= mean;

//Defining common scope for POD-SPOD
    {
        Eigen::VectorXd lambda(settings.Ns);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        int Nrec;
    
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
            std::cout << "Number of modes extracted: " << Nrec << std::endl;
            Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nrec, Nrec);

            for ( int i = 0; i < Nrec; i++ )
                Sig(i,i) = std::sqrt(lambda(i));

            std::cout << "Computing error and Jaccard index surface..." << "\t";

            Eigen::MatrixXd Err_SPOD_map = Eigen::MatrixXd::Zero( sn_set.rows(), sn_set.cols() );
            Eigen::MatrixXd Err_SPOD_Nm_time = Eigen::MatrixXd::Zero( settings.Ns, settings.Ns );
            Eigen::MatrixXd J_SPOD_Nm_time = Eigen::MatrixXd::Zero( settings.Ns, settings.Ns );

            for ( int nm = 1; nm <= Nrec; nm ++ )
            {
                Err_SPOD_map = sn_set - Phi.leftCols(nm)*Sig.block(0,0,nm,nm)*eig_vec.transpose().topRows(nm);

                for ( int i = 0; i < settings.Ns; i++ )
                {

                    int count = 0;
                    for ( int j = 0; j < settings.ndim*Nr; j++ )
                    {
                        Err_SPOD_Nm_time(i,nm-1) += Err_SPOD_map(j,i)*Err_SPOD_map(j,i);
                    
                        if ( Err_SPOD_map(j,i) < 0.1)
                            count++;
                    }
                    
                    Err_SPOD_Nm_time(i,nm-1) = std::sqrt(Err_SPOD_Nm_time(i,nm-1))/norm_sn_set(i);
                    J_SPOD_Nm_time(i,nm-1) = (double)count/((double)settings.ndim*(double)Nr);
                
                }
            
            }

            if ( Phi.cols() < settings.Ns )
            {
                for ( int plus = 0; plus < (settings.Ns - Phi.cols()); plus++ )
                {
                    Err_SPOD_Nm_time.col(Phi.cols() + plus) = Err_SPOD_Nm_time.col(Phi.cols() - 1);
                    J_SPOD_Nm_time.col(Phi.cols() + plus) = J_SPOD_Nm_time.col(Phi.cols() - 1);
                }
            }

            Err_RBM_Nm_time.push_back(Err_SPOD_Nm_time);
            J_RBM_Nm_time.push_back(J_SPOD_Nm_time);

            std::cout << "Done" << std::endl;
        }
    
    }


    
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
        std::cout << "Number of modes extracted : " << Nm << std::endl;

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
        Eigen::MatrixXd Err_DMD_Nm_time = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
        Eigen::MatrixXd J_DMD_Nm_time = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);

        for ( int nm = 1; nm <= Nm; nm ++ )
        {

            Eigen::MatrixXcd D_dmd = Phi.leftCols(nm)*Psi.topRows(nm);
            Err_DMD_map = sn_set - D_dmd.real();

            for ( int i = 0; i < settings.Ns; i++ )
            {
                int count = 0;
                for ( int j = 0; j < settings.ndim*Nr; j++ )
                { 
                    Err_DMD_Nm_time(i,nm-1) += Err_DMD_map(j,i)*Err_DMD_map(j,i);
                    if ( Err_DMD_map(j,i) < 0.1 )
                        count++;
                }

                Err_DMD_Nm_time(i,nm-1) = std::sqrt(Err_DMD_Nm_time(i,nm-1))/norm_sn_set(i);
                J_DMD_Nm_time(i,nm-1) = (double)count/((double)settings.ndim*(double)Nr);
            }

        }
        
        if ( Phi.cols() < settings.Ns )
        {
            for ( int plus = 0; plus < (settings.Ns - Phi.cols()); plus++ )
            {
                Err_DMD_Nm_time.col(Phi.cols() + plus) = Err_DMD_Nm_time.col(Phi.cols() - 1);
                J_DMD_Nm_time.col(Phi.cols() + plus) = J_DMD_Nm_time.col(Phi.cols() - 1);
            }
        }

        Err_RBM_Nm_time.push_back(Err_DMD_Nm_time);
        J_RBM_Nm_time.push_back(J_DMD_Nm_time);

    }
    

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
            Phi = read_modes( file_modes, settings.ndim*Nr, settings.r_RDMD );
            Coefs = Phi.transpose()*sn_set;
        }

        std::cout << " Done! " << std::endl;

        std::cout << "Computing error and Jaccard index surface ... " << std::endl << std::endl;
        Eigen::MatrixXd Err_RDMD_map;
        Eigen::MatrixXd Err_RDMD_Nm_time = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);
        Eigen::MatrixXd J_RDMD_Nm_time = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns);

        for ( int nm = 1; nm <= settings.r_RDMD; nm ++ )
        {

            Err_RDMD_map = sn_set - Phi.leftCols(nm)*Coefs.topRows(nm);
            for ( int i = 0; i < settings.Ns; i++ )
            {
                int count = 0;

                for ( int j = 0; j < settings.ndim*Nr; j++ )
                {
                    Err_RDMD_Nm_time(i,nm-1) += Err_RDMD_map(j,i)*Err_RDMD_map(j,i);
                    if ( Err_RDMD_map(j,i) < 0.1 )
                        count++;
                }

                Err_RDMD_Nm_time(i,nm-1) = std::sqrt(Err_RDMD_Nm_time(i,nm-1))/norm_sn_set(i);
                J_RDMD_Nm_time(i,nm-1) = (double)count/((double)settings.ndim*(double)Nr);
            }

        }

        // if ( Phi.cols() < settings.Ns )
        // {
        //     for ( int plus = 0; plus < (settings.Ns - Phi.cols()); plus++ )
        //     {    
        //         Err_RDMD_Nm_time.col(Phi.cols() + plus) = Err_RDMD_Nm_time.col(Phi.cols() - 1);
        //         J_RDMD_Nm_time.col(Phi.cols() + plus) = J_RDMD_Nm_time.col(Phi.cols() - 1);
        //     }
        // }

        Err_RBM_Nm_time.push_back(Err_RDMD_Nm_time);
        J_RBM_Nm_time.push_back(J_RDMD_Nm_time);

    }


    std::cout << "Writing error and Jaccard index error files ... " << "\t";
    // -------- Nm ---------
    // |                    |
    // |                    |
    // t        E,J         |
    // |                    |
    // |                    |
    // ----------------------
    for ( int i = 0; i < Nf.size(); i++ )
    {
        
        std::string filename1 = "Error_SPOD_Nf_" + std::to_string(i) + ".dat";
        write_err_j( Err_RBM_Nm_time[i], filename1 );
        std::string filename2 = "Jaccard_SPOD_Nf_" + std::to_string(i) + ".dat";
        write_err_j( J_RBM_Nm_time[i], filename2 );

    }

    std::string filename1 = "Error_DMD.dat";
    write_err_j( Err_RBM_Nm_time[s_Nf], filename1 );
    std::string filename2 = "Jaccard_DMD.dat";
    write_err_j( J_RBM_Nm_time[s_Nf], filename2 );


    filename1 = "Error_RDMD.dat";
    write_err_j( Err_RBM_Nm_time[s_Nf+1], filename1 );
    filename2 = "Jaccard_RDMD.dat";
    write_err_j( J_RBM_Nm_time[s_Nf+1], filename2 );

    return 0;

}