#include "Basis_Extraction.hpp"


/**
 * Default constructor - Empty 
 */
basis_extraction::basis_extraction() :
    m_Ns(0), m_Nr(0), m_Nf(0)
{
}


basis_extraction::basis_extraction( const int &Ns, 
                                    const int &Nr, 
                                    const int &Nf,
                                    const Eigen::MatrixXd sn_set ) :
    m_Ns(Ns), m_Nr(Nr), m_Nf(Nf)
{

    m_snset = sn_set;
    m_lam = Eigen::VectorXd(Ns);
    m_K_pc = Eigen::VectorXd(Ns);
    m_eig_vec = Eigen::MatrixXd(Ns, Ns);

}

// Destructor
basis_extraction::~basis_extraction() {

}



Eigen::MatrixXd basis_extraction::SPOD_basis( std::string bc_flag, 
                                        std::string filter_flag, 
                                        std::string meanflag, 
                                        double sigma){

    int count;

    Eigen::MatrixXd R_f(m_Ns, m_Ns);
    Eigen::MatrixXd phi_c(m_Nr, m_Ns);
    Eigen::VectorXd mean(m_Nr);
    Eigen::VectorXd g(2*m_Nf+1);

    R_f.setZero(m_Ns, m_Ns);

    mean = m_snset.rowwise().mean();

    if ( meanflag == "YES" ){
        for (int i = 0; i < m_Ns; i++)
            m_snset.col(i) -= mean;
    }

    if ( m_Nf == 0) {                           //Performing pure POD
        R_f = m_snset.transpose()*m_snset;

    }else{                                      //Performing SPOD
    
        if ( filter_flag == "BOX"){             //with BOX filter  

            for (int k = 0; k <= 2*m_Nf; k++ )
                g(k) = 1.0/(2.0*(double)m_Nf+1.0); 

        }else if ( filter_flag == "GAUSSIAN"){  //with GAUSSIAN filter
            
            if ( sigma == 0.0){
                std::cout << "sigma = 0 then only POD could be performed" << std::endl;
                exit (EXIT_FAILURE);
                }

            count = 0; 
            for (int k = -m_Nf; k <= m_Nf; k++ ){
                    g(count) = exp(-k*k/(2.0*sigma*sigma));
                    count++;
            }

        }else{
            std::cout << "Filter selected not available" << std::endl;
            exit (EXIT_FAILURE);
        }

        Eigen::MatrixXd R = m_snset.transpose()*m_snset;

        //Zero-padded boundary conditions
        if (bc_flag == "ZERO"){

            for (int i = 0; i < m_Ns; i++){
                for (int j = 0; j < m_Ns; j++){

                    count = 0;

                    for (int k = -m_Nf; k <= m_Nf; k++){
                        if (i + k < 0 || j + k < 0 || i + k >= m_Ns || j + k >= m_Ns)
                            R_f(i,j) += 0.0;
                        else
                            R_f(i,j) += g(count)*R(i+k, j+k);
                        
                        count++; 
        
                    }
                }
            }
        }else{
            std::cout << "Booundary condition not implemented " << std::endl;
            exit (EXIT_FAILURE);
        }
    }


    Eigen::EigenSolver<Eigen::MatrixXd> es(R_f); 
    Eigen::VectorXd m_lam = es.eigenvalues().real();
    Eigen::MatrixXd m_eig_vec = es.eigenvectors().real();
    eig_sort();

    double sum = 0;

    for (int i = 0; i < m_Ns; i++){
        sum += m_lam(i)/m_lam.sum();
        m_K_pc(i) = sum;
    }

    double tol = m_lam(0)*1e-12;
    phi_c = m_snset*m_eig_vec;

    count = 0;
    while ( count < m_lam.size() && abs(m_lam(count)) > tol)
            count++;

    Eigen::MatrixXd phi(m_Nr,count);
    for ( int i = 0 ; i < count ; i++ )
        phi.col(i) = phi_c.col(i)/sqrt(m_lam(i));



    return phi;

}


void basis_extraction::eig_sort() {

    unsigned int swap_count = 1;
    double temp;
    Eigen::VectorXd temp_vec(m_eig_vec.rows());

    while (swap_count > 0){

        swap_count = 0;

        for(unsigned int index = 1; index < m_lam.size(); index++){

            if (m_lam(index) > m_lam(index-1)){

                temp = m_lam(index-1);
                m_lam(index-1) = m_lam(index);
                m_lam(index) = temp;

                temp_vec = m_eig_vec.col(index-1);
                m_eig_vec.col(index-1) = m_eig_vec.col(index);
                m_eig_vec.col(index) = temp_vec;

                swap_count++;
            }
        }
    }
}


Eigen::VectorXd basis_extraction::get_eigValues() {
    return m_lam;
}

Eigen::MatrixXd basis_extraction::get_eigVectors() {
    return m_eig_vec;
}

Eigen::VectorXd basis_extraction::get_Kpc() {
    return m_K_pc;
}







