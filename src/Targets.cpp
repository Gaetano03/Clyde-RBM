#include "Targets.hpp"

Targets::Targets()
{
}




Targets::Targets( const Eigen::MatrixXd Phi,
            const Eigen::MatrixXd Coeffs,
            const int t_vec,
            const int t,
            const std::string flag_interp ) {

                m_Phi = Phi;
                m_Coeffs = Coeffs;
                m_t_vec = t_vec;
                m_time = t;
                m_flag_interp = flag_interp;

            }


Targets::~Targets()
{    
}


void Targets::set_target_time( double t ) {

    m_time = t;
}


key Targets::get_key ( const std::string &key_string ){

       if(key_string == "LINEAR")
        return LINEAR;
    else if(key_string == "CUBIC")
        return CUBIC;
    else if(key_string == "GAUSSIAN")
        return GAUSSIAN;
    else if(key_string == "THIN_PLATE")
        return THIN_PLATE;
    else if(key_string == "MULTIQUADRATICS")
        return MULTIQUADRATICS;


}


Eigen::VectorXd Targets::Reconstruction () {

    std::vector< std::vector<double> > T( m_t_vec.size(), std::vector<double>(1));
    std::vector<double> t(1, m_time);

    double avgDt = 0.0;

    for ( int i = 0; i < m_t_vec.size(); i++ ) {
        T[i][0] = m_t_vec[i];
    }

    for ( int i = 1; i < m_t_vec.size(); i++ ) {
        Dt += m_t_vec[i] - m_t_vec[i-1];
    }

    Dt = Dt/(double)(m_t_vec.size()-1);

    //Vector of surrogate coefficients
    vector<double> coefs_intrp( m_t_vec.size() );

    // Create surrogates for coefficients
    vector<rbf> surr_coefs;

    
    RBF_CONSTANTS rbf_const {Dt, 0.0};

    for ( int i = 0; i < m_t_vec.size(); i++ ){
        
        vector<double> coefs ;
        for (int j = 0 ; j < Ns ; j++)
            coefs.push_back(Coeffs(i,j));

        // cout << "Coefficients to interpolate : \n" << Coeffs.row(i) << endl;
        // CUBIC, GAUSSIAN, LINEAR, MULTIQUADRATICS




        surr_coefs.push_back( rbf(T, coefs, LINEAR, rbf_const) );
        surr_coefs[i].build();
        surr_coefs[i].evaluate(t,coefs_intrp[i]);

        // cout << "Value of the interpolated coefficient : " << coefs_intrp[i] << endl;
   
    }

    VectorXd coefs_t(Ns);

    for (int i = 0; i < Ns; i++)
        coefs_t(i) = coefs_intrp[i]; 

    //Reconstruction of the field for which RBM has been applied

    return coefs_t;


}