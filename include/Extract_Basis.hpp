#ifndef EXTRACT_BASIS_HPP
#define EXTRACT_BASIS_HPP


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

#include "LinearAlgebra/Eigen/Dense"
#include "LinearAlgebra/Eigen/Eigenvalues"


int Nmod ( double En, Eigen::VectorXd K_pc );

void eig_sort( Eigen::VectorXd &lam, Eigen::MatrixXd &eig_vec );

int SVHT ( Eigen::VectorXd lam, int m, int n );



Eigen::MatrixXd SPOD_basis( const Eigen::MatrixXd &snap_set,
                                Eigen::VectorXd &lam,
                                Eigen::VectorXd &K_pc,
                                Eigen::MatrixXd &eig_vec,
                                const int Nf = 0,
                                std::string bc_flag = "ZERO", 
                                std::string filter_flag = "BOX",  
                                double sigma = 1.0);



Eigen::MatrixXcd DMD_basis( const Eigen::MatrixXd &snap_set,
                            Eigen::VectorXcd &lam,
                            Eigen::MatrixXcd &eig_vec,
                            Eigen::VectorXd &lam_POD,
                            Eigen::MatrixXd &eig_vec_POD,
                            Eigen::VectorXd &K_pc,
                            const double En );


class Mr_DMD {

protected:

    //Number of levels for multi-resolution analysis
    int m_L;

    //Number of snapshots
    int m_Ns;

    //mr_DMD modes
    std::vector<Eigen::MatrixXcd> Modes;

    //mr_DMD Coefficients
    std::vector<Eigen::MatrixXcd> Coefs;
    

public:


    //Default constructor
    Mr_DMD();

    //Constructor
    Mr_DMD( const int L,
            const int Ns,
            const double En,
            const Eigen::MatrixXd &sn_set );

    
    

};


#endif //EXTRACT_BASIS_HPP