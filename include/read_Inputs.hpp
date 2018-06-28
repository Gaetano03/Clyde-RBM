#ifndef READ_INPUTS_HPP
#define READ_INPUTS_HPP

#include <iostream>
#include <sstream>
#include <stdio.h> 
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include "LinearAlgebra/Eigen/Dense"
#include "LinearAlgebra/Eigen/Eigenvalues"



// Structure to be filled with information from cfg file
struct prob_settings 
{

    int Ns;
    int Nf;
    double En;
    double Ds;
    std::string in_file_root;
    std::string in_file_format;
    std::string out_file;
    std::string dim_prob;
    std::string flag_filter;
    std::string flag_mean;
    std::vector<int> Cols;

};


// List of Keywords in config file
enum keywords 
            { 
                NS, DS, EN, PROB_DIM, NF, 
                INPUT_FILE_ROOT, 
                INPUT_FILE_FORMAT, 
                COLS_FIELDS 
            };


// Compare keywords with input string
keywords read_keyword_type( const std::string &key_string );


// Read config file and store info in prob_settings
void Read_cfg ( prob_settings &settings );


//Get Number of grid points
int N_gridpoints ( const std::string file_in );


// Get fields on the specified columns 
Eigen::MatrixXd read_col( std::string filename );


#endif //READ_INPUTS_HPP
