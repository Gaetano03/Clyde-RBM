#include "Reconstruction.hpp"


smartuq::surrogate::RBF_FUNCTION get_key_rbf ( const std::string &key_string ) 
{

    if(key_string == "LINEAR")
        return smartuq::surrogate::LINEAR;
    else if(key_string == "CUBIC")
        return smartuq::surrogate::CUBIC;
    else if(key_string == "GAUSSIAN")
        return smartuq::surrogate::GAUSSIAN;
    else if(key_string == "THIN_PLATE")
        return smartuq::surrogate::THIN_PLATE;
    else if(key_string == "MULTIQUADRATICS")
        return smartuq::surrogate::MULTIQUADRATICS;


}



Eigen::MatrixXd Reconstruction_S_POD ( const std::vector<double> &t_vec,
                                const Eigen::VectorXd &K_pc,
                                const Eigen::VectorXd &lam,
                                const Eigen::MatrixXd &Coeffs,
                                const Eigen::MatrixXd &phi,
                                const double time,
                                const double En,
                                std::string flag_prob,
                                std::string flag_interp ) 
                                {

                                    std::vector< std::vector<double> > T( t_vec.size(), std::vector<double>(1));
                                    std::vector<double> t(1, time);

                                    double avgDt = 0.0;

                                    for ( int i = 0; i < t_vec.size(); i++ ) {
                                        T[i][0] = t_vec[i];
                                    }

                                    for ( int i = 1; i < t_vec.size(); i++ ) {
                                        avgDt += t_vec[i] - t_vec[i-1];
                                    }

                                    avgDt = avgDt/(double)(t_vec.size()-1);

                                    //Vector of surrogate coefficients
                                    std::vector<double> coefs_intrp( t_vec.size() );

                                    // Create surrogates for coefficients
                                    std::vector<rbf> surr_coefs;
                                    RBF_CONSTANTS rbf_const {avgDt, 0.0};

                                    // Define the number of modes Nrec to use in the reconstruction
                                    int Nrec = Nmod(En, K_pc);

                                    for ( int i = 0; i < Nrec; i++ ){
                                        
                                        std::vector<double> coefs ;
                                        for (int j = 0 ; j < t_vec.size() ; j++)
                                            coefs.push_back(Coeffs(i,j));

                                        surr_coefs.push_back( rbf(T, coefs, get_key_rbf( flag_interp ), rbf_const) );               
                                        surr_coefs[i].build();
                                        surr_coefs[i].evaluate(t, coefs_intrp[i]);
                                
                                    }

                                    Eigen::VectorXd coefs_t(t_vec.size());

                                    for (int i = 0; i < t_vec.size() ; i++)
                                        coefs_t(i) = coefs_intrp[i];

                                    Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nrec, Nrec);

                                    for ( int i = 0; i < Nrec; i++ )
                                        Sig(i, i) = sqrt(lam(i));


                                    if ( flag_prob == "SCALAR" ) 
                                    {

                                        Eigen::MatrixXd Rec_field = phi.leftCols(Nrec)*Sig*coefs_t.head(Nrec);  
                                        return Rec_field;

                                    } else if ( flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D" ) 
                                    {

                                        Eigen::MatrixXd Rec_field (phi.rows()/2, 2);
                                        Rec_field.col(0) = phi.topLeftCorner(phi.rows()/2, Nrec)*Sig*coefs_t.head(Nrec);
                                        Rec_field.col(1) = phi.bottomLeftCorner(phi.rows()/2, Nrec)*Sig*coefs_t.head(Nrec);
                                        return Rec_field;


                                    } else if ( flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D" ) 
                                    {

                                        Eigen::MatrixXd Rec_field (phi.rows()/3, 3);
                                        Rec_field.col(0) = phi.topLeftCorner(phi.rows()/3, Nrec)*Sig*coefs_t.head(Nrec); 
                                        Rec_field.col(1) = phi.block(phi.rows()/3, 0, phi.rows()/3, Nrec)*Sig*coefs_t.head(Nrec);
                                        Rec_field.col(2) = phi.bottomLeftCorner(phi.rows()/3, Nrec)*Sig*coefs_t.head(Nrec);                                        
                                        return Rec_field;

                                    } else 
                                    {

                                            std::cout << " Set well flag_prob! Now Exiting ..." << std::endl;
                                            exit (EXIT_FAILURE);

                                    }

                                }


Eigen::VectorXcd Calculate_Coefs_DMD ( const Eigen::MatrixXcd &eig_vec,
                                    const Eigen::MatrixXcd &eig_vec_POD,
                                    const Eigen::VectorXcd &lam,
                                    const Eigen::VectorXcd &lam_POD,
                                    const int Nm )
                                    {


                                        Eigen::MatrixXcd Sigma = Eigen::MatrixXcd::Zero(Nm, Nm);
                                        Eigen::MatrixXcd V_and(Nm, lam_POD.size()); 

                                        for ( int i = 0; i < Nm; i++)
                                            Sigma(i, i) = std::sqrt(lam_POD(i));
                                      
                                        for ( int i = 0; i < Nm; i++ )
                                        {
                                            for ( int j = 0; j < lam_POD.size(); j++ )
                                                V_and(i,j) = std::pow(lam(i), j);                                                                                         
                                        }

                                        Eigen::MatrixXcd Y_sq = eig_vec.conjugate().transpose()*eig_vec;
                                        Eigen::MatrixXcd V_and_sq = V_and*V_and.conjugate().transpose();
                                        V_and_sq = V_and_sq.conjugate();

                                        Eigen::MatrixXcd M1(Nm, Nm);

                                        for ( int i = 0; i < Nm; i++ )
                                        {

                                            for ( int j = 0; j < Nm; j++)
                                                M1(i,j) = Y_sq(i,j)*V_and_sq(i,j);

                                        }

                                        Eigen::MatrixXcd M2 = V_and*eig_vec_POD.leftCols(Nm)*Sigma.conjugate()*eig_vec;
                                        Eigen::VectorXcd dum = M2.diagonal();
                                        dum = dum.conjugate();

                                        return M1.inverse()*dum;
                                       

                                    } 


Eigen::MatrixXcd Reconstruction_DMD ( const double time, const double dt, 
                                    const Eigen::VectorXcd &alfa,
                                    const Eigen::MatrixXcd &Phi,
                                    const Eigen::VectorXcd &lam,
                                    const std::string flag_prob )
                                    {

                                        int Nm = lam.size();
                                        Eigen::VectorXcd omega(Nm);

                                        for ( int i = 0; i < Nm; i ++)
                                        {
/*Remove after check*/                      std::cout << " Lambda_" << i << " = " << lam(i) << "\t"; 
                                            omega(i) = std::log(lam(i))/dt;
/*Remove after check*/                      std::cout << " omega_" << i << " = " << omega(i) << std::endl;
                                        }


                                        Eigen::MatrixXcd v_rec (Phi.rows(),1);

                                        for ( int i = 0; i < Nm; i++ )
                                        {
                                            
                                            v_rec += std::exp(omega(i)*time)*alfa(i)*Phi.col(i);
        
                                        }

                                        if ( flag_prob == "SCALAR" )
                                        {

                                            return v_rec;

                                        } else if ( flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D" )
                                        {

                                            Eigen::MatrixXcd Rec(Phi.rows()/2,2);
                                            Rec.col(0) = v_rec.topRows(Phi.rows()/2);
                                            Rec.col(1) = v_rec.bottomRows(Phi.rows()/2);
                                            return Rec;

                                        } else if ( flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D" )
                                        {

                                            Eigen::MatrixXcd Rec(Phi.rows()/3,3);
                                            Rec.col(0) = v_rec.topRows(Phi.rows()/3);
                                            Rec.col(1) = v_rec.middleRows(Phi.rows()/3,Phi.rows()/3);
                                            Rec.col(2) = v_rec.bottomRows(Phi.rows()/3);
                                            return Rec;

                                        } else 
                                        {

                                            std::cout << "Set well problem flag! Exiting ... " << std::endl;
                                            exit (EXIT_FAILURE);

                                        }

                                    } 



Eigen::MatrixXcd TimeEvo_DMD ( Eigen::VectorXd &time,
                            double dt,
                            const Eigen::VectorXcd &alfa,
                            const Eigen::MatrixXcd &Phi,
                            const Eigen::VectorXcd &lam )
                            {

                                int Nm = lam.size();
                                Eigen::VectorXcd omega(Nm);

                                for ( int i = 0; i < Nm; i ++)
                                {
/*Remove after check*/              std::cout << " Lambda_" << i << " = " << lam(i) << "\t"; 
                                    omega(i) = std::log(lam(i))/dt;
/*Remove after check*/              std::cout << " omega_" << i << " = " << omega(i) << std::endl;
                                }

                                Eigen::MatrixXcd v_rec = Eigen::MatrixXcd::Zero(Phi.rows(),time.size());

                                for ( int k = 0; k < time.size(); k++ )
                                {
                                    for ( int i = 0; i < Nm; i++ )
                                    {
                                        v_rec.col(k) += std::exp(omega(i)*time(k))*alfa(i)*Phi.col(i);
                                        // std::cout << " exp(omega_" << i << "*t) : " << std::exp(omega(i)*time(k)) << std::endl;
                                    }
                                }

                                return v_rec;

                            }