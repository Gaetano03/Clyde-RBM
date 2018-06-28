#ifndef GENERATE_SNSET_HPP
#define GENERATE_SNSET_HPP


#include "read_Inputs.hpp"

Eigen::MatrixXd generate_snap_matrix( const int Nr, const int Ns, const int ds,
                                        std::vector<int> Cols,
                                        std::string root_inputfile,
                                        std::string input_format,
                                        std::string flag_prob = "VELOCITY", 
                                        std::string flag_dim = "2D");



#endif // GENERATE_SNSET_HPP