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


Eigen::MatrixXd SPOD_basis( const Eigen::MatrixXd snap_set,
                                const int Ns, const int Nr,
                                std::string bc_flag = "ZERO", 
                                std::string filter_flag = "BOX", 
                                std::string meanflag = "YES", 
                                double sigma = 1.0);




#endif //EXTRACT_BASIS_HPP