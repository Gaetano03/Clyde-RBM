#include "Extract_Basis.hpp"

void eig_sort( Eigen::VectorXd &lam, Eigen::MatrixXd &eig_vec ) {

    unsigned int swap_count = 1;
    double temp;
    Eigen::VectorXd temp_vec(eig_vec.rows());

    while (swap_count > 0)
    {

        swap_count = 0;

        for(unsigned int index = 1; index < lam.size(); index++)
        {

            if ( lam(index) > lam(index-1) )
            {

                temp = lam(index-1);
                lam(index-1) = lam(index);
                lam(index) = temp;

                temp_vec = eig_vec.col(index-1);
                eig_vec.col(index-1) = eig_vec.col(index);
                eig_vec.col(index) = temp_vec;

                swap_count++;

            }
        }
    }
}





Eigen::MatrixXd SPOD_basis( const Eigen::MatrixXd &snap_set,
                                const int Ns, const int Nr,     
                                std::string bc_flag , 
                                std::string filter_flag , 
                                std::string meanflag , 
                                double sigma )
                                {

                                    int count;

                                    Eigen::MatrixXd R_f(Ns, Ns);
                                    Eigen::MatrixXd phi_c(Nr, Ns);
                                    Eigen::VectorXd mean(Nr);
                                    Eigen::VectorXd g(2*Nf+1);

                                    R_f.setZero(Ns, Ns);

                                    mean = snap_set.rowwise().mean();

                                    if ( meanflag == "YES" )
                                    {
                                        for ( int i = 0; i < Ns; i++ )
                                            snap_set.col(i) -= mean;
                                    }

                                    if ( Nf == 0) 
                                    {                           //Calculating R POD
                                        R_f = snap_set.transpose()*snap_set;

                                    } else
                                    {                                      //Calculating R SPOD
                                    
                                        if ( filter_flag == "BOX" )
                                        {             //with BOX filter  

                                            for (int k = 0; k <= 2*Nf; k++ )
                                                g(k) = 1.0/(2.0*(double)Nf+1.0); 

                                        } else if ( filter_flag == "GAUSSIAN" )
                                        {             //with GAUSSIAN filter
                                            
                                            if ( sigma == 0.0)
                                            {
                                                std::cout << "sigma = 0 then only POD could be performed" << std::endl;
                                                exit (EXIT_FAILURE);
                                            }
 
                                            for ( int k = 0; k <= 2*Nf; k++ )
                                                g(k) = exp(-k*k/(2.0*sigma*sigma));
                                            

                                        } else
                                        {
                                            std::cout << "Filter selected not available" << std::endl;
                                            exit (EXIT_FAILURE);
                                        }

                                        Eigen::MatrixXd R = snap_set.transpose()*snap_set;

                                        //Zero-padded boundary conditions
                                        if (bc_flag == "ZERO"){

                                            for ( int i = 0; i < Ns; i++ )
                                            {
                                                for ( int j = 0; j < Ns; j++ )
                                                {

                                                    count = 0;

                                                    for ( int k = -Nf; k <= Nf; k++ ){
                                                        if ( i + k < 0 || j + k < 0 || i + k >= Ns || j + k >= Ns )
                                                            R_f(i,j) += 0.0;
                                                        else
                                                            R_f(i,j) += g(count)*R(i+k, j+k);
                                                        
                                                        count++; 
                                        
                                                    }
                                                }
                                            }
                                        } else 
                                        {
                                            std::cout << "Booundary condition not implemented " << std::endl;
                                            exit (EXIT_FAILURE);
                                        }
                                    }


                                    Eigen::EigenSolver<Eigen::MatrixXd> es(R_f); 
                                    Eigen::VectorXd lam = es.eigenvalues().real();
                                    Eigen::MatrixXd eig_vec = es.eigenvectors().real();
                                    eig_sort( lam, eig_vec);

                                    double sum = 0;

                                    for (int i = 0; i < Ns; i++){
                                        sum += lam(i)/lam.sum();
                                        K_pc(i) = sum;
                                    }

                                    double tol = lam(0)*1e-12;
                                    phi_c = snap_set*eig_vec;

                                    count = 0;
                                    while ( count < lam.size() && abs(lam(count)) > tol)
                                            count++;

                                    Eigen::MatrixXd phi(Nr,count);
                                    for ( int i = 0 ; i < count ; i++ )
                                        phi.col(i) = phi_c.col(i)/sqrt(lam(i));


                                    return phi;   


                                }