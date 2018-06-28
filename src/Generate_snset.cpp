#include "read_Inputs.hpp"
#include "Generate_snset.hpp"

Eigen::MatrixXd generate_snap_matrix( const int Nr, const int Ns, const int ds,
                                        std::vector<int> Cols,
                                        std::string root_inputfile,
                                        std::string input_format,
                                        std::string flag_prob, 
                                        std::string flag_dim )
{

    Eigen::MatrixXd field(Nr, Cols.size());
    std::string file_temp;
    int k = 0;

    if ( flag_prob == "VECTOR")
    {

        if ( flag_dim == "2D"){
            
            Eigen::MatrixXd snap(2*Nr, Ns);

            for( int i = 0; i < Ns*ds; i += ds )
            {

                std::stringstream buffer;
                buffer << std::setfill('0') << std::setw(5) << std::to_string(i);
                file_temp = root_inputfile + buffer.str() + input_format;
                std::cout << "Reading fields from : " << file_temp << "\t";
                field = read_col(file_temp, Nr, Cols, field);
                std::cout << "Complete!" << std::endl;

                Eigen::VectorXd gx = field.col(1);
                Eigen::VectorXd gy = field.col(2);

                snap.col(k) << gx,
                                gy;
                
                k++;
            }

            return snap;

        } else if ( flag_dim == "3D")
        {

            Eigen::MatrixXd snap(3*Nr, Ns);

            for( int i = 0; i < Ns*ds; i += ds )
            {

                std::stringstream buffer;
                buffer << std::setfill('0') << std::setw(5) << std::to_string(i);
                file_temp = root_inputfile + buffer.str() + input_format;
                std::cout << "Reading fields from : " << file_temp << "\t";
                field = read_col(file_temp, Nr, Cols, field);
                std::cout << "Complete!" << std::endl;

                Eigen::VectorXd gx = field.col(1);
                Eigen::VectorXd gy = field.col(2);
                Eigen::VectorXd gz = field.col(3);

                snap.col(k) << gx,
                                gy,
                                gz;

                k++;

            }   

        return snap;

        } else 
        {

            std::cout << " Set well problem dimension! Exiting ... " << std::endl;
            exit (EXIT_FAILURE);
            
        } 


    } else if ( flag_prob == "VELOCITY" ) 
    {


        if ( flag_dim == "2D" )
        {

            Eigen::MatrixXd snap(2*Nr, Ns);


            for( int i = 0; i < Ns*ds; i += ds )
            {

                std::stringstream buffer;
                buffer << std::setfill('0') << std::setw(5) << std::to_string(i);
                file_temp = root_inputfile + buffer.str() + input_format;
                std::cout << "Reading fields from : " << file_temp << "\t";
                field = read_col(file_temp, Nr, Cols, field);
                std::cout << "Complete!" << std::endl;

                Eigen::VectorXd rho = field.col(0);
                Eigen::VectorXd rho_u = field.col(1);
                Eigen::VectorXd rho_v = field.col(2);
                Eigen::VectorXd u = rho_u.cwiseQuotient(rho);
                Eigen::VectorXd v = rho_v.cwiseQuotient(rho);

                snap.col(k) << u,
                                v;
                
                k++;
            }

            return snap;

        } else if ( flag_dim == "3D" )
        {

            Eigen::MatrixXd snap(3*Nr, Ns);

            for( int i = 0; i < Ns*ds; i += ds )
            {

                std::stringstream buffer;
                buffer << std::setfill('0') << std::setw(5) << std::to_string(i);
                file_temp = root_inputfile + buffer.str() + input_format;
                std::cout << "Reading fields from : " << file_temp << "\t";
                field = read_col(file_temp, Nr, Cols, field);
                std::cout << "Complete!" << std::endl;

                Eigen::VectorXd rho = field.col(0);
                Eigen::VectorXd rho_u = field.col(1);
                Eigen::VectorXd rho_v = field.col(2);
                Eigen::VectorXd rho_w = field.col(3);
                Eigen::VectorXd u = rho_u.cwiseQuotient(rho);
                Eigen::VectorXd v = rho_v.cwiseQuotient(rho);
                Eigen::VectorXd w = rho_w.cwiseQuotient(rho);

                snap.col(k) << u,
                                v,
                                w;

                k++;

            }

        return snap;

        } else 
        {
            
            std::cout << " Set well problem dimension! Exiting ... " << std::endl;
            exit (EXIT_FAILURE);            

        }
    
    } else if ( flag_prob == "SCALAR" )
    {

            Eigen::MatrixXd snap(Nr, Ns);

            for( int i = 0; i < Ns*ds; i += ds )
            {

                std::stringstream buffer;
                buffer << std::setfill('0') << std::setw(5) << std::to_string(i);
                file_temp = root_inputfile + buffer.str() + input_format;
                std::cout << "Reading fields from : " << file_temp << "\t";
                field = read_col(file_temp, Nr, Cols, field);
                std::cout << "Complete!" << std::endl;

                snap.col(k) = field.col(0);

                k++;

            }

        return snap;


    }

}




