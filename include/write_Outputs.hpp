#ifndef WRITE_OUTPUTS_HPP
#define WRITE_OUTPUTS_HPP

#include <iostream>
#include <stdio.h> 
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "LinearAlgebra/Eigen/Dense"
#include "Reconstruction.hpp"
#include "read_Inputs.hpp"

void Config_stream ( prob_settings settings );

void write_modes_sPOD ( const Eigen::MatrixXd &Phi_cut, 
                        const Eigen::MatrixXd &Coords,
                        std::string flag_prob );

void write_modes_DMD ( const Eigen::MatrixXcd &Phi_cut,
                    const Eigen::MatrixXd &Coords, 
                    std::string flag_prob );

void write_coeffs_sPOD ( const Eigen::MatrixXd &Coeffs,
                        const std::vector<double> &t_vec,
                        const Eigen::VectorXd &lam );

void write_TimeDynamics_DMD ( const Eigen::VectorXcd omega,
                            const Eigen::VectorXcd alfa,
                            const Eigen::VectorXd t);

void write_CoefsDynamics_mrDMD( std::vector<node_mrDMD> &nodes, 
                                const int level, 
                                const int ns, 
                                const int max_levels );

void write_Reconstructed_fields ( const Eigen::MatrixXd Rec,
                                    const Eigen::MatrixXd &Coords,
                                    std::string filename,
                                    std::string flag_prob,
                                    const int nt );

void write_modes ( const Eigen::MatrixXd &Phi_cut ); //write modes RDMD

void write_coefs ( const Eigen::MatrixXd &Coefs ); //write coefs RDMD

void write_err_j ( const Eigen::MatrixXd data, std::string filename); //write error/jaccard surface for RBM method


#endif // WRITE_OUTPUTS_HPP