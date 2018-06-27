#ifndef BASIS_EXTRACTION_HPP
#define BASIS_EXTRACTION_HPP

#include "ctime"
#include <iostream>
#include <stdio.h> 
#include <vector>
#include <fstream>
#include <sstream>
#include <complex> 
#include <cmath>
#include <iomanip>
#include <string>
//#include "smartmath.h"
//#include "smartuq.h"
#include "LinearAlgebra/Eigen/Dense"
#include "LinearAlgebra/Eigen/Eigenvalues"
//#include "LinearAlgebra/Eigen/Core"
//#include "LinearAlgebra/Eigen/Householder"


class basis_extraction {
    
protected:

    // Number of snapshots
    int m_Ns;

    //Number of grid points
    int m_Nr;

    // Filter size for SPOD
    int m_Nf;

    // Vector of eigenvalues
    Eigen::VectorXd m_lam;

    // Vector of energetic content
    Eigen::VectorXd m_K_pc;

    // Eigenvector matrix
    Eigen::MatrixXd m_eig_vec;  

    // Matrix of snaps
    Eigen::MatrixXd m_snset;



public: 

    /**
     * Default constructor - Empty 
     */
    basis_extraction();

    // Input number of snapshot, Number of grid points, SPOD filter size
    basis_extraction( const int &Ns, 
                    const int &Nr,
                    const int &Nf,
                    const int &snset );

    // Destructor
    virtual ~basis_extraction();


    // Order eigenvalues and eigenvectors columns in descending order
    void eig_sort();

    // Get Eingenvalues of the correlation matrix
    Eigen::VectorXd get_eigValues();

    // Get EigenVectors of the correlation matrix
    Eigen::MatrixXd get_eigVectors();

    // Get Vector of the energetic content
    Eigen::VectorXd get_Kpc();

    //Calculate POD-SPOD modes (POD when Nf = 0)
    /**Performing by default a SPOD with zero-padded boundary condition,
    box filter and subtracting the mean to the data**/
    Eigen::MatrixXd SPOD_basis( std::string bc_flag = "ZERO", 
                                            std::string filter_flag = "BOX", 
                                            std::string meanflag = "YES", 
                                            double sigma = 1.0);



};

#endif // BASIS_EXTRACTION_HPP