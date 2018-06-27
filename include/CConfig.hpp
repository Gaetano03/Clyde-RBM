#ifndef CCONFIG_HPP
#define CCONFIG_HPP


#include <iostream>
#include <sstream>
#include <stdio.h> 
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

// Structure to be filled with information from cfg file
struct prob_settings {

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
enum keywords { NS, DS, EN, PROB_DIM, NF, INPUT_FILE_ROOT, INPUT_FILE_FORMAT, COLS_FIELDS};


// Class which manage the reading of config file
class CConfig {

protected:

    // Structure with problem settings
    prob_settings m_settings;

    // Config file name
    std::string m_filename;

public:

    /**
     * Default constructor
     **/
    CConfig();

    // Input config filename
    CConfig( std::string filename );

    // Destructor
    ~CConfig();

    // Get all problem settings (structure prob_settings)
    prob_settings get_settings();

    // Compare keywords with input string
    keywords read_keyword_type(const std::string &key_string);

    // Read config file and store info in prob_settings
    void Read_cfg ();


};

#endif // CCONFIG_HPP