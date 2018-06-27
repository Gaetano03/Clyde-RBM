#ifndef CREATE_SNAPSET_HPP
#define CREATE_SNAPSET_HPP

#include <iostream>
#include <stdio.h> 
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip> //for setting the output stream format (setprecision, setfill, setw, etc...)
//#include "smartmath.h"
//#include "smartuq.h"
#include "LinearAlgebra/Eigen/Dense"
#include "LinearAlgebra/Eigen/Eigenvalues"
//#include "LinearAlgebra/Eigen/Core"
//#include "LinearAlgebra/Eigen/Householder"

class Create_snapset {

protected:

    // Number of grid points
    int m_Nr;

    // Number of snapshots
    int m_Ns;

    // Delta between snapshots
    int m_ds;

    // Root of filename with fields
    std::string m_root_inputfile;

    // Columns to be read 
    std::vector<int> m_Cols;

    // Flag for problem dimension
    std::string m_flag1;

    // Flag for derived quantities
    std::string m_flag2;

    // Flag for input format
    std::string m_input_format;

public:

    /**
     * Default Constructor
    **/
   Create_snapset();

   // Constructor with 
   Create_snapset( const int &Ns,                               //Number of snaps 
                const int &ds,                                  //Delta between snaps 
                const std::vector<int> Cols,                    //Columns to be read
                const std::string filename,                     //Root filename of fields to be read
                const std::string flag1,                        //flag for problem dimension
                const std::string flag2 = "VELOCITY",           //flag for derived fields
                const std::string input_format = ".dat" );      //file format to read


    // Destructor
    ~Create_snapset();

    // Set the number of grid points
    void set_gridpoints();    

    // Get the number of grid points
    int get_gridpoints();

    // Get fields on the specified columns 
    Eigen::MatrixXd read_col( std::string filename );

    // Generate matrix with desired fields
    Eigen::MatrixXd generate_snap_matrix();


};

#endif //CREATE_SNAPSET_HPP