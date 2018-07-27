#include "Extract_Basis.hpp"


int Nmod ( double En, Eigen::VectorXd K_pc )
{

    double En_content = 0.0;
    int count = 0;

    while ( En_content < En && count < K_pc.size() )
    {
        En_content = K_pc(count);
        count++;
    }

    return count;

}



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

// Should work only with singular values ( Not eigenvalues problems )
int SVHT ( Eigen::VectorXd lam, int m, int n )
{

    double median;
    int n_sv = lam.size();

    for ( int i = 0; i < n_sv; i++ ) 
        lam(i) = std::sqrt(lam(i));

    std::cout << "lambda :\n" << lam << std::endl;

    if ( n_sv%2 != 0 )
        median = lam((n_sv-1)/2);
    else
        median = 0.5*(lam(n_sv/2) + lam(n_sv/2-1));
    
    std::cout << "median :\n" << median << std::endl;

    double beta = std::min(m,n)/std::max(n,m);
    beta = 1.0;
    double omega = 0.56*std::pow(beta,3.0) - 0.95*std::pow(beta,2.0) 
                    + 1.82*beta + 1.43;
    
    double tau = omega*median;

    std::cout << "tau :\n" << tau << std::endl;

    double eps = 1000;
    int Nm = 0;

    do
    {
        eps = lam(Nm);
        Nm++;
    } while ( eps > tau );

    return Nm-1;

}


Eigen::MatrixXd SPOD_basis( const Eigen::MatrixXd &snap_set,
                                Eigen::VectorXd &lam,
                                Eigen::VectorXd &K_pc,
                                Eigen::MatrixXd &eig_vec,
                                const int Nf,    
                                std::string bc_flag , 
                                std::string filter_flag ,  
                                double sigma )
                                {

                                    int count;
                                    int Nr = snap_set.rows();
                                    int Ns = snap_set.cols();

                                    Eigen::MatrixXd R_f(Ns, Ns);
                                    Eigen::MatrixXd phi_c(Nr, Ns);
                                    Eigen::VectorXd mean(Nr);
                                    Eigen::VectorXd g(2*Nf+1);

                                    R_f.setZero(Ns, Ns);

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
                                            std::cout << "Boundary condition not implemented " << std::endl;
                                            exit (EXIT_FAILURE);
                                        }
                                    }


                                    Eigen::EigenSolver<Eigen::MatrixXd> es(R_f); 
                                    lam = es.eigenvalues().real();
                                    eig_vec = es.eigenvectors().real();
                                    eig_sort( lam, eig_vec);

                                    double sum = 0;

                                    for (int i = 0; i < Ns; i++){
                                        sum += lam(i)/lam.sum();
                                        K_pc(i) = sum;
                                    }

                                    // double tol = lam(0)*1e-12;
                                    double tol = 1e-16;
                                    phi_c = snap_set*eig_vec;

                                    count = 0;
                                    while ( count < lam.size() && std::abs(lam(count)) > tol)
                                            count++;

                                    Eigen::MatrixXd phi(Nr,count);
                                    for ( int i = 0 ; i < count ; i++ )
                                        phi.col(i) = phi_c.col(i)/sqrt(lam(i));


                                    return phi;   


                                }


Eigen::MatrixXcd DMD_basis ( const Eigen::MatrixXd &snap_set,
                            Eigen::VectorXcd &lam,
                            Eigen::MatrixXcd &eig_vec,
                            Eigen::VectorXd &lam_POD,
                            Eigen::MatrixXd &eig_vec_POD,
                            Eigen::VectorXd &K_pc,
                            const double En )
                            {   

                                int Ns = snap_set.cols() - 1;
 
                                Eigen::MatrixXd U = SPOD_basis(snap_set.leftCols(Ns), lam_POD, K_pc, eig_vec_POD );
                                Eigen::MatrixXd Sig_inv = Eigen::MatrixXd::Zero(U.cols(), U.cols());

                                std::cout << "Number of non-zero modes : " << U.cols() << std::endl;

                                for ( int i = 0; i < U.cols(); i++ )
                                    Sig_inv(i, i) = 1.0/std::sqrt(lam_POD(i)); 

                                int Nm = Nmod( En, K_pc );

                                if ( En == 1.0 )
                                    Nm = U.cols();
                                
                                //Nm = SVHT ( lam_POD, Ns, snap_set.rows() );
                                // Nm = 49;


                                Eigen::MatrixXd Atilde = U.leftCols(Nm).transpose()*snap_set.rightCols(Ns)*
                                                            eig_vec_POD.leftCols(Nm)*Sig_inv.block(0,0,Nm,Nm);

                                Eigen::EigenSolver<Eigen::MatrixXd> es(Atilde); 
                                lam = es.eigenvalues();
                                eig_vec = es.eigenvectors();

                                Eigen::MatrixXcd appo = snap_set.rightCols(Ns)*
                                                        eig_vec_POD.leftCols(Nm)*Sig_inv.block(0,0,Nm,Nm);

                                

                                // Non divido per lambda
                                // return appo*eig_vec;

                                // Divido per lambda

                                Eigen::MatrixXcd phi(snap_set.rows(), Nm);
                                for (int i = 0; i < Nm; i++)
                                    phi.col(i) = 1.0/lam(i)*appo*eig_vec.col(i);
// 
                                return phi;
                                
                                //Standard DMD
                                // return U*eig_vec;
                            }


