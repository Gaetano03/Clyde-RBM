#include "CConfig.hpp"


CConfig::CConfig( std::string filename ) {
    m_filename = filename;
}

CConfig::~CConfig() {
}


prob_settings CConfig::get_settings() {
    return m_settings;
}


keywords CConfig::read_keyword_type(const std::string &key_string) {
    
    if(key_string == "NS")
        return NS;
    else if(key_string == "DS")
        return DS;
    else if(key_string == "EN")
        return EN;
    else if(key_string == "PROB_DIM")
        return PROB_DIM;
    else if(key_string == "INPUT_FILE_ROOT")
        return INPUT_FILE_ROOT;
    else if(key_string == "INPUT_FILE_FORMAT")
        return INPUT_FILE_FORMAT;
    else if(key_string == "NF")
        return NF;
    else if(key_string == "COLS_FIELDS")
        return COLS_FIELDS;
    else{
        std::cout << "Something wrong in cfg file" << std::endl;
        exit (EXIT_FAILURE);       
    }

}


void CConfig::Read_cfg (){

    
    std::ifstream cFile (m_filename);
    if (cFile.is_open()){

        size_t delimiterPos; 
        std::string name, value;
        std::string line;

        while(getline(cFile, line)){
            line.erase(remove_if(line.begin(), line.end(), isspace),
                                 line.end());  //include ::  in front of isspace if using namespace std
            if(line[0] == '#' || line.empty())
                continue;
            delimiterPos = line.find("=");
            name = line.substr(0, delimiterPos);
            value = line.substr(delimiterPos + 1);

            switch (read_keyword_type(name)){

                case PROB_DIM:
                {
                    m_settings.dim_prob = value;
                    std::cout << "Problem dimension : " << value << std::endl;
                    break;
                }

                case NS:
                {
                    m_settings.Ns = std::stoi(value);
                    std::cout << "Number of snapshots : " << value << std::endl;
                    break;
                }

                case DS:
                {
                    m_settings.Ds = std::stod(value);
                    std::cout << "Delta between selected snapshots (Equispaced) : " << value << std::endl; 
                    break;
                }

                case EN:
                {
                    m_settings.En = std::stod(value);
                    std::cout << "Energy content used in the reconstruction : " << value << std::endl;
                    break;
                }

                case INPUT_FILE_ROOT:
                {
                    m_settings.in_file_root = value;
                    std::cout << "Input file root name : " << value << std::endl;
                    break;
                }

                case INPUT_FILE_FORMAT:
                {
                    m_settings.in_file_format = value;
                    std::cout << "Input file format : " << value << std::endl;
                    break;
                }

                case NF:
                {
                    m_settings.Nf = std::stoi(value);
                    std::cout << "Filter size for feature extraction : " << value << std::endl;
                    break;
                }

                case COLS_FIELDS:
                {
                    std::string str = value;
                    std::stringstream ss(str);

                    int i;

                    while ( ss >> i ) {

                        m_settings.Cols.push_back(i);

                        if (ss.peek() != ',' || ss.peek() != ' ')
                            ss.ignore();

                    }

                    std::cout << "Number of columns to process: \t";
                    
                    for (   i = 0; i < m_settings.Cols.size(); i++ )
                        std::cout << m_settings.Cols[i] << "\t";

                    std::cout << std::endl;

                    break;
                }

                default:
                {
                    break;
                }

            }

        }
        cFile.close();
    }
    else
        std::cout << "Unable to open config file." << '\n';


}