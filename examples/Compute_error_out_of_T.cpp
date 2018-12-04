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
    std::vector<Eigen::VectorXd> Err_RBM_Nm_time;
    std::vector<Eigen::VectorXd> J_RBM_Nm_time;
    std::vector<Eigen::VectorXd> EN;

    // Calculate number of grid points
    int Nr = N_gridpoints ( settings.in_file );
    std::cout << "Number of grid points : " << Nr << std::endl;

    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings.Ns, settings.Ds, settings.nstart,
                                        settings.Cols,
                                        settings.in_file,
                                        settings.flag_prob);

    std::cout << "Initializing Vector of time ... " << std::endl; 
    std::vector<double> t_vec( settings.Ns );
    t_vec[0] = settings.nstart*settings.Dt_cfd;
    for ( int i = 1; i < settings.Ns; i++ )
        t_vec[i] = t_vec[i-1] + settings.Dt_cfd*settings.Ds;
    
    std::vector<double> t(settings.Ns-1);
    t[0] = settings.nstart*settings.Dt_cfd + settings.Dt_cfd*(double)settings.Ds/2.0;
    for ( int i = 1; i < t.size(); i++)
        t[i] = t[i-1] + settings.Dt_cfd*(double)settings.Ds;

    std::cout << "Computing mean of CFD solution ... " << std::endl;
    Eigen::VectorXd mean = sn_set.rowwise().mean();
    
    Eigen::VectorXd norm_sn_set = Eigen::VectorXd::Zero(settings.Ns-1);

    std::cout << "Reading matrix for check ... " << std::endl;
    Eigen::MatrixXd sn_set_check = generate_snap_matrix( Nr, settings.Ns-1, settings.Ds, settings.nstart + settings.Ds/2,
                                                        settings.Cols,
                                                        settings.in_file,
                                                        settings.flag_prob);
    
    Eigen::VectorXd mean_check = sn_set.rowwise().mean();

    for ( int i = 0; i < settings.Ns-1; i ++ )
    {
        norm_sn_set(i) = sn_set_check.col(i).norm();
    }

// //Defining common scope for POD-SPOD
    {
        for ( int nt = 0; nt < settings.Ns; nt++ )
            sn_set.col(nt) -= mean;

        Eigen::VectorXd lambda(settings.Ns);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        int Nm;

        for ( int nfj = 0; nfj < Nf.size(); nfj++ )
        {

            std::cout << "Extracting SPOD " << Nf[nfj] << " basis ... " << "\t";        
 
            Eigen::MatrixXd Phi = SPOD_basis( sn_set,
                                    lambda, K_pc, eig_vec,
                                    Nf[nfj],
                                    settings.flag_bc, 
                                    settings.flag_filter,  
                                    settings.sigma);

            std::cout << " Done! " << std::endl;            

            Nm = Nmod(settings.En, K_pc);
            std::cout << "Number of modes for desired energetic content: " << Nm << std::endl;

            std::vector<rbf> surr_coefs =  getSurrCoefs (t_vec,
                                                        eig_vec,
                                                        settings.flag_interp);
            
            Eigen::MatrixXd coef_t(settings.Ns-1, Nm);


            std::vector<double> tr(1);
            for ( int j = 0; j < settings.Ns - 1; j++ )
            {    
                tr[0] = t[j];
                for ( int i = 0; i < Nm; i++ )
                    surr_coefs[i].evaluate(tr, coef_t(j,i));
            }
            

            Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);

            for ( int i = 0; i < Nm; i++ )
                Sig(i,i) = std::sqrt(lambda(i));

            std::cout << "Computing error and Jaccard index vector..." << "\t";

            Eigen::MatrixXd Err_SPOD_map = Eigen::MatrixXd::Zero( sn_set.rows(), sn_set_check.cols() );
            Eigen::VectorXd Err_SPOD_Nm_time = Eigen::VectorXd::Zero(settings.Ns-1);
            Eigen::VectorXd J_SPOD_Nm_time = Eigen::VectorXd::Zero(settings.Ns-1);
            
            Err_SPOD_map = sn_set_check - Phi.leftCols(Nm)*Sig*coef_t.transpose();
            
            for ( int i = 0; i < Err_SPOD_map.cols(); i++ )
            {
                Err_SPOD_map.col(i) -= mean; 
            }

            for ( int i = 0; i < settings.Ns-1; i++ )
            {

                int count = 0;
                for ( int j = 0; j < settings.ndim*Nr; j++ )
                {
                    Err_SPOD_Nm_time(i) += Err_SPOD_map(j,i)*Err_SPOD_map(j,i);
                
                    if ( Err_SPOD_map(j,i) < 0.1)
                        count++;
                }
                
                Err_SPOD_Nm_time(i) = std::sqrt(Err_SPOD_Nm_time(i))/norm_sn_set(i);
                J_SPOD_Nm_time(i) = (double)count/((double)settings.ndim*(double)Nr);
            
            }     

            Err_RBM_Nm_time.push_back(Err_SPOD_Nm_time);
            J_RBM_Nm_time.push_back(J_SPOD_Nm_time);
            EN.push_back(K_pc);

            std::cout << "Done" << std::endl;
        }
    
    }


    for ( int nt = 0; nt < settings.Ns; nt++ )
        sn_set.col(nt) += mean;
    
// //Defining scope for DMD ( Rank=-1 preferable, Coeffs = OPT )
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
//         int Nm = Phi.cols();
//         std::cout << "Number of modes extracted : " << Nm << std::endl;

        Eigen::VectorXcd omega(Phi.cols());
        for ( int i = 0; i < Phi.cols(); i++ )
            omega(i) = std::log(lambda_DMD(i))/(settings.Dt_cfd*settings.Ds);

        std::cout << "Calculating coefficients DMD ... " << "\t";            
        alfa = Calculate_Coefs_DMD_exact ( sn_set.leftCols(settings.Ns-1),  
                                                            lambda_DMD, 
                                                            Phi );
        std::cout << " Done! " << std::endl;

        std::cout << "Reordering modes DMD ... " << "\t";
        Eigen::VectorXd En = Eigen::VectorXd::Zero(Phi.cols());
        double T = t_vec[t_vec.size()-1];

        for ( int i = 0 ; i < Phi.cols(); i ++ )
        {

            double alfa_i = alfa(i).imag();
            double alfa_r = alfa(i).real();
            double sigma = omega(i).real();
            En(i) = (alfa_r*alfa_r + alfa_i*alfa_i)*(std::exp(2.0*sigma*T) - 1.0)/(2.0*sigma);

        }

        dmd_sort( En, Phi, lambda_DMD, alfa);
        std::cout << "Done" << std::endl;

        double sum = 0;
        Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
        for (int i = 0; i < Phi.cols(); i++)
        {
            sum += En(i)/En.sum();
            K_pc(i) = sum;
        }

        int Nm = Nmod(settings.En, K_pc);
        std::cout << "Number of modes for the desired energetic content : " << Nm << std::endl;

        std::cout << "Computing error and Jaccard index surface ... " << std::endl << std::endl;

        Eigen::MatrixXcd V_and(lambda_DMD.size(), settings.Ns-1);      
        for ( int i = 0; i < lambda_DMD.size(); i++ )
        {
            for ( int j = 0; j < settings.Ns-1; j++ )
                V_and(i,j) = std::pow(lambda_DMD(i), (double)j + 0.5);                                                                                         
        }        
        Eigen::MatrixXcd Psi = Eigen::MatrixXcd::Zero(alfa.size(), settings.Ns-1);
        for ( int i = 0; i < settings.Ns-1; i++ )
            Psi.col(i) = alfa.cwiseProduct(V_and.col(i));
  
        Eigen::MatrixXd Err_DMD_map(sn_set.rows(), sn_set_check.cols());
        Eigen::VectorXd Err_DMD_Nm_time = Eigen::VectorXd::Zero(settings.Ns-1);
        Eigen::VectorXd J_DMD_Nm_time = Eigen::VectorXd::Zero(settings.Ns-1);


        Eigen::MatrixXcd D_dmd = Phi.leftCols(Nm)*Psi.topRows(Nm);
        Err_DMD_map = sn_set_check - D_dmd.real();

        for ( int i = 0; i < settings.Ns-1; i++ )
        {
            int count = 0;
            for ( int j = 0; j < settings.ndim*Nr; j++ )
            { 
                Err_DMD_Nm_time(i) += Err_DMD_map(j,i)*Err_DMD_map(j,i);
                if ( Err_DMD_map(j,i) < 0.1 )
                    count++;
            }

            Err_DMD_Nm_time(i) = std::sqrt(Err_DMD_Nm_time(i))/norm_sn_set(i);
            J_DMD_Nm_time(i) = (double)count/((double)settings.ndim*(double)Nr);
        }

        
        

        Err_RBM_Nm_time.push_back(Err_DMD_Nm_time);
        J_RBM_Nm_time.push_back(J_DMD_Nm_time);
        EN.push_back(K_pc);
    }
    
    for ( int nt = 0; nt < settings.Ns; nt++ )
        sn_set.col(nt) -= mean;

// //Defining scope for RDMD
    {

        Eigen::VectorXd lambda = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::VectorXd K_pc = Eigen::VectorXd::Zero(settings.Ns);
        Eigen::MatrixXd Coefs = Eigen::MatrixXd::Zero(settings.Ns, settings.Ns); 
        Eigen::MatrixXd Phi;

        if ( argc == 2 )
        {
            std::cout << "Extracting basis and Coeffs RDMD ... " << "\t";        
        
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
            std::string file_En = argv[4];
            Phi = read_modes( file_modes, settings.ndim*Nr, settings.r_RDMD );
            Coefs = read_coefs( file_coefs, settings.Ns, settings.r_RDMD );

            std::ifstream En_data;
            En_data.open( file_En );
            if ( !En_data.is_open() )
            {
                std::cout << "File : " << file_En << " not found" << std::endl;    
                exit (EXIT_FAILURE);
            }
            std::string line_flow_data ;
            getline( En_data, line_flow_data );
            std::istringstream iss(line_flow_data);
            std::string token;

            int count = 0;
            while( getline( iss, token, ' ') )
            {
                K_pc(count) = std::stod(token);
                count ++;
            } 
            En_data.close();
        }

        std::cout << " Done! " << std::endl;

        int Nm = Nmod(settings.En, K_pc);
        std::cout << "number of modes for the desired energetic content " << Nm << std::endl;

        std::vector<rbf> surr_coefs =  getSurrCoefs (t_vec,
                                                    Coefs.transpose(),
                                                    settings.flag_interp);
        
        Eigen::MatrixXd coef_t(settings.Ns-1, Nm);

        std::vector<double> tr(1);
        for ( int j = 0; j < settings.Ns - 1; j++ )
        {    
            tr[0] = t[j];
            for ( int i = 0; i < Nm; i++ )
                surr_coefs[i].evaluate(tr, coef_t(j,i));
        }


        std::cout << "Computing error and Jaccard index surface ... " << std::endl << std::endl;
        Eigen::MatrixXd Err_RDMD_map;
        Eigen::VectorXd Err_RDMD_Nm_time = Eigen::VectorXd::Zero(settings.Ns-1);
        Eigen::VectorXd J_RDMD_Nm_time = Eigen::VectorXd::Zero(settings.Ns-1);

 
        Err_RDMD_map = sn_set - Phi.leftCols(Nm)*coef_t.transpose();
        for ( int i = 0; i < settings.Ns-1; i++ )
        {
            int count = 0;

            for ( int j = 0; j < settings.ndim*Nr; j++ )
            {
                Err_RDMD_Nm_time(i) += Err_RDMD_map(j,i)*Err_RDMD_map(j,i);
                if ( Err_RDMD_map(j,i) < 0.1 )
                    count++;
            }

            Err_RDMD_Nm_time(i) = std::sqrt(Err_RDMD_Nm_time(i))/norm_sn_set(i);
            J_RDMD_Nm_time(i) = (double)count/((double)settings.ndim*(double)Nr);
        }
 
        Err_RBM_Nm_time.push_back(Err_RDMD_Nm_time);
        J_RBM_Nm_time.push_back(J_RDMD_Nm_time);
        EN.push_back(K_pc);
    }


    std::cout << "Writing error and Jaccard index error files ... " << std::endl;

    std::ofstream errfile;
    errfile.open("Err_RBM.dat");

    for ( int nm = 0; nm < settings.Ns-1; nm ++ )    
    {
        for( int j = 0; j < Err_RBM_Nm_time.size(); j++ ) 
            errfile <<  std::setprecision(8) << Err_RBM_Nm_time[j](nm) << "\t";

        errfile << std::endl;

    }

    errfile.close();


    std::cout << "Writing energetic content ... " << std::endl;

    std::ofstream datafile;
    datafile.open("Encontent_RBM.dat");

    for( int j = 0; j < EN.size(); j++ ) 
    {
        for ( int nm = 0; nm < settings.Ns; nm ++ )
            datafile <<  std::setprecision(8) << EN[j](nm) << "\t";

        datafile << std::endl;

    }

    datafile.close();

    std::cout << "RBM-Clyde error surface calculation ends" << std::endl << std::endl;

    return 0;

}

