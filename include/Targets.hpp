#ifndef TARGETS_HPP
#define TARGETS_HPP

#include <iostream>
#include <stdio.h> 
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include "LinearAlgebra/Eigen/Dense"
#include "LinearAlgebra/Eigen/Eigenvalues"
#include "smartmath.h"
#include "smartuq.h"


//List of keywords for interpolation technique of RBF function



class Targets {

protected:

    //Matrix of modes
    Eigen::MatrixXd m_Phi;

    //Matrix of coefficients
    Eigen::MatrixXd m_Coeffs;

    //Vector of input time
    std::vector<double> m_t_vec;

    //Target time instant
    double m_time; 

    //Flag for interpolation methods in RBF
    std::string m_flag_interp;

    //Keywords for RBM interp
    enum key { 
        LINEAR,
        CUBIC,
        GAUSSIAN,
        THIN_PLATE,
        MULTIQUADRATICS }

public:


    /*
    **Default constructor 
     */
    Target ();

    //Constructor with modes set, Coeffs set
    Target ( const Eigen::MatrixXd Phi,
            const Eigen::MatrixXd Coeffs,
            const int t_vec,
            const int t,
            const std::string flag_interp = "LINEAR" );

    //Destructor
    ~Target ();

    //set target time
    void set_target_time( double t );


    key get_key ( const std::string &key_string );

    //Reconstruction
    Eigen::VectorXd Reconstruction ();





};









#endif //TARGETS_HPP