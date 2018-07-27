#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"

int main(int argc, char *argv[]) {

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "-----------RBM-Clyde start-------------" << std::endl << std::endl;

    std::string filecfg = argv[1];
    // std::string mode = argv[2];
    prob_settings settings;

    //Reading configuration file
    Read_cfg( filecfg, settings );
    std::cout << std::endl;
    
    // Calculate number of grid points
    int Nr = N_gridpoints ( settings.in_file );
    std::cout << "Number of grid points : " << Nr << std::endl;
    
    // Reading coordinates
    std::cout << "Reading Coordinates ... \t ";
    Eigen::MatrixXd Coords = read_col( settings.in_file, Nr, settings.Cols_coords );
    std::cout << "Done " << std::endl;


    // Create matrix of snapshots
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings.Ns, settings.Ds,
                                        settings.Cols,
                                        settings.in_file,
                                        settings.flag_prob);


    Eigen::VectorXd mean = sn_set.rowwise().mean();

    if ( settings.flag_method == "SPOD")
    {

        std::vector<double> t_vec( settings.Ns );
        t_vec[0] = 0.0;

        for ( int i = 1; i < settings.Ns; i++ )
            t_vec[i] = t_vec[i-1] + settings.Dt_cfd*settings.Ds;

        std::cout << std::endl;
        std::cout << "Initialized vector of times " << std::endl;

        Eigen::VectorXd lambda(settings.Ns);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);

        if ( settings.flag_mean == "YES" )
        {
            std::cout << "Subtracting mean from snapshots ... " << std::endl << std::endl;
            for ( int i = 0; i < settings.Ns; i++ )
                sn_set.col(i) -= mean;
        }


        std::cout << "Extracting basis ... " << "\t";        

        Eigen::MatrixXd Phi = SPOD_basis( sn_set,
                                lambda, K_pc, eig_vec,
                                settings.Nf,
                                settings.flag_bc, 
                                settings.flag_filter,  
                                settings.sigma);

        std::cout << " Done! " << std::endl << std::endl;

        std::cout << "Number of non-zero modes : " << Phi.cols() << std::endl;

        int Nrec = Nmod( settings.En, K_pc);
        std::cout << "Number of modes for the desired energy content : " << Nrec << std::endl;

        if ( settings.flag_wdb_be == "YES" )
        {
            
            std::cout << "Writing modes ..." << "\t";
            write_modes_sPOD ( Phi.leftCols(Nrec), Coords, settings.flag_prob );
            std::cout << "Complete!" << std::endl;

            std::cout << "Writing Coefficients ..." << "\t";
            write_coeffs_sPOD ( eig_vec.transpose(), t_vec, lambda );
            std::cout << "Complete!" << std::endl;
            std::cout << std::endl;

        }

        if ( settings.flag_rec == "YES" )
        {

            std::cout << "Reconstructing field at time : " << settings.t_rec << "\t";

            Eigen::MatrixXd Rec = Reconstruction_S_POD ( t_vec,
                                K_pc, lambda, eig_vec.transpose(),
                                Phi, settings.t_rec,
                                settings.En,
                                settings.flag_prob,
                                settings.flag_interp ) ;

            std::cout << "Done" << std::endl;

            if ( settings.flag_mean == "YES" )
            {

                for ( int i = 0; i < Rec.cols(); i++)
                    Rec.col(i) = Rec.col(i) + mean.segment(i*Nr, Nr);

            }

            std::cout << "Writing reconstructed field ..." << "\t";

            write_Reconstructed_fields ( Rec, Coords,
                                    settings.out_file,
                                    settings.flag_prob );

            std::cout << "Done" << std::endl << std::endl;

        }


    }

    if ( settings.flag_method == "DMD" )
    {


        Eigen::VectorXd t_vec( settings.Ns );
        t_vec(0) = 0.0;

        for ( int i = 1; i < settings.Ns; i++ )
            t_vec(i) = t_vec(i-1) + settings.Dt_cfd*settings.Ds;

        std::cout << std::endl;
        std::cout << "Initialized vector of times " << std::endl;

        Eigen::VectorXd lambda_POD(settings.Ns - 1);
        Eigen::VectorXd K_pc(settings.Ns - 1);
        Eigen::MatrixXd eig_vec_POD(settings.Ns - 1, settings.Ns - 1);
        Eigen::VectorXcd lambda_DMD(settings.Ns - 1);
        Eigen::MatrixXcd eig_vec_DMD(settings.Ns - 1, settings.Ns - 1);
        

        if ( settings.flag_mean == "YES" )
        {
            std::cout << "Subtracting mean from snapshots ... " << std::endl << std::endl;
            for ( int i = 0; i < settings.Ns; i++ )
                sn_set.col(i) -= mean;
        }

        std::cout << "Extracting basis ... " << "\t";        

        Eigen::MatrixXcd Phi = DMD_basis( sn_set,
                                        lambda_DMD,
                                        eig_vec_DMD,
                                        lambda_POD,
                                        eig_vec_POD,
                                        K_pc,
                                        settings.En );

        int Nm = Nmod( settings.En, K_pc);

        std::cout << " Done! " << std::endl << std::endl;

        std::cout << "Calculating coefficients DMD ... " << "\t";

        Eigen::MatrixXcd alfa = Calculate_Coefs_DMD ( eig_vec_DMD,
                                            eig_vec_POD,
                                            lambda_DMD,
                                            lambda_POD,
                                            Nm );

        std::cout << " Done! " << std::endl << std::endl;

        if ( settings.flag_wdb_be == "YES" )
        {
            
            std::cout << "Writing modes ..." << "\t";
            write_modes_DMD ( Phi, Coords, settings.flag_prob );
            std::cout << "Complete!" << std::endl;

            Eigen::VectorXcd omega(Nm);
            for ( int i = 0; i < Nm; i++ )
                 omega(i) = std::log(lambda_DMD(i))/(settings.Dt_cfd*settings.Ds);

            std::cout << "Writing time dynamics ..." << "\t";
            write_TimeDynamics_DMD ( omega, alfa, t_vec );
            std::cout << "Complete!" << std::endl;
            std::cout << std::endl;

        }


        if ( settings.flag_rec == "YES" )
        {

            std::cout << "Reconstructing field at time : " << settings.t_rec << "\t";

            Eigen::MatrixXcd Rec = Reconstruction_DMD ( settings.t_rec,
                                                    settings.Dt_cfd*settings.Ds,
                                                    alfa,
                                                    Phi,
                                                    lambda_DMD,
                                                    settings.flag_prob );

            std::cout << "Done" << std::endl;

            if ( settings.flag_mean == "YES" )
            {

                for ( int i = 0; i < Rec.cols(); i++)
                    Rec.col(i) = Rec.col(i) + mean.segment(i*Nr, Nr);

            }

            std::cout << "Writing reconstructed field ..." << "\t";

            write_Reconstructed_fields ( Rec.real(), Coords,
                                    settings.out_file,
                                    settings.flag_prob );

            std::cout << "Done" << std::endl << std::endl;

        }

    }

    std::cout << std::endl;    
    std::cout << "-----------RBM-Clyde end-------------" << std::endl << std::endl;

    return 0;

}
