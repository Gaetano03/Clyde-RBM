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
#include "LinearAlgebra/Eigen/MatrixFunctions"

int Nmod ( double En, Eigen::VectorXd K_pc );


void eig_sort( Eigen::VectorXd &lam, Eigen::MatrixXd &eig_vec );


int SVHT ( Eigen::VectorXd lam, int m, int n );


struct node_mrDMD
{
    
    int l;                              //level number
    int bin_num;                        //time bin number
    int bin_size;                       //time bin size
    int step;                           //Delta between snaps based on nyquist frequency
    int start, stop;                    //starting and stopping index for snapshots
    double rho;                         //cut-off frequency
    double t_begin, t_end;              //initial and final instant of each window with respect to the whole interval
    double dt;
    int r;                              //rank reduction
    int n;                              //number of slow modes
    Eigen::MatrixXcd Modes;             //Matrix of modes
    Eigen::VectorXcd Coefs;             //Vector of optimized coefs
    Eigen::VectorXcd lam;               //Vector of eigenvalues
    Eigen::MatrixXcd Psi;               //Time evolution matrix

};


Eigen::MatrixXd nullspace( Eigen::VectorXd s, Eigen::MatrixXd vh, const double atol = 1e-13, const double rtol = 0 );


bool check_linear_consistency( Eigen::MatrixXd X, Eigen::MatrixXd Y, Eigen::MatrixXd Nullspace, const double atol = 1e-13, const double rtol = 0);


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
                            const int r = 0 );


Eigen::VectorXcd Calculate_Coefs_DMD ( const Eigen::MatrixXcd &eig_vec,
                                    const Eigen::MatrixXcd &eig_vec_POD,
                                    const Eigen::VectorXcd &lam,
                                    const Eigen::VectorXcd &lam_POD,
                                    const int Ns );


Eigen::VectorXcd Calculate_Coefs_DMD_exact ( const Eigen::MatrixXd &sn_set,  //matrix of first Ns-1 snaps 
                                            const Eigen::VectorXcd &lam,  //slow eigenvalues
                                            const Eigen::MatrixXcd &Phi ); //slow exact DMD modes


Eigen::MatrixXcd Calculate_Coefs_Matrix_DMD ( const Eigen::MatrixXd &sn_set,
                                            const Eigen::MatrixXcd &Phi,
                                            const Eigen::VectorXcd &omega,
                                            const double t_0,
                                            const double dt_dmd);


std::vector<node_mrDMD> mrDMD_basis( Eigen::MatrixXd &snap_set,       //Initial set of snapshots
                                    std::vector<node_mrDMD> &nodes,          //Initialize as empty vector
                                    const int r,                                  //DMD rank
                                    double dts,                             //dt between initial snapshots
                                    double t_0 = 0.0,                       //time instant of the first snapshot of the series
                                    int level = 0,                          
                                    int bin_num = 0,
                                    int offset = 0,
                                    int max_levels = 7,
                                    int max_cycles = 2,
                                    std::string flag_coefs = "OPT" );



Eigen::MatrixXd RDMD_modes_coefs ( const Eigen::MatrixXd &sn_set,
                                    Eigen::MatrixXd &Coefs,     //Define N_mod RDMD through the dimension of matrix Coefs
                                    Eigen::VectorXd &lambda,
                                    const int r );              //rank of pure DMD at each step



Eigen::MatrixXcd fbDMD_basis ( const Eigen::MatrixXd &snap_set,
                            Eigen::VectorXcd &lam,
                            Eigen::MatrixXcd &eig_vec,
                            const int r );

#endif //EXTRACT_BASIS_HPP