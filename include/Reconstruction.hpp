#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include <iostream>
#include <sstream>
#include <stdio.h> 
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <complex>
#include "Extract_Basis.hpp"
#include "LinearAlgebra/Eigen/Dense"
#include "LinearAlgebra/Eigen/Eigenvalues"
#include "smartmath.h"
#include "smart-uq/Surrogates/rbf.h"

using namespace smartuq::surrogate;


smartuq::surrogate::RBF_FUNCTION get_key_rbf ( const std::string &key_string ); 


Eigen::MatrixXd Reconstruction_S_POD ( const std::vector<double> &t_vec,
                                const Eigen::VectorXd &K_pc,
                                const Eigen::VectorXd &lam,
                                const Eigen::MatrixXd &Coeffs,
                                const Eigen::MatrixXd &phi,
                                const double time,
                                const double En,
                                std::string flag_prob = "VELOCITY-2D",
                                std::string flag_interp = "LINEAR");


Eigen::MatrixXcd Reconstruction_DMD ( const double time, const double dt, 
                                    const Eigen::VectorXcd &alfa,
                                    const Eigen::MatrixXcd &Phi,
                                    const Eigen::VectorXcd &lam,
                                    const std::string flag_prob );

Eigen::MatrixXcd TimeEvo_DMD ( Eigen::VectorXd &time,
                            double dt,
                            const Eigen::VectorXcd &alfa,
                            const Eigen::MatrixXcd &Phi,
                            const Eigen::VectorXcd &lam );
 

void nodes_mrDMD_sort( std::vector<node_mrDMD> &nodes ); //Sorting nodes in ordered levels


Eigen::MatrixXcd Reconstruction_mrDMD ( const double time,                      //time desired for reconstruction                                                                                                                                                                                                                                                                                                                              
                                    const double dts,                           //time between initial set of snapshots
                                    const std::vector<node_mrDMD> &nodes,       //nodes
                                    const std::string flag_prob );              

#endif // RECONSTRUCTION_HPP

