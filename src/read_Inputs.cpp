#include "read_Inputs.hpp"


keywords read_keyword_type( const std::string &key_string )
{

    if( key_string == "NS" )
        return NS;
    else if( key_string == "DS" )
        return DS;
    else if( key_string == "EN" )
        return EN;
    else if( key_string == "PROB_DIM" )
        return PROB_DIM;
    else if( key_string == "INPUT_FILE_ROOT" )
        return INPUT_FILE_ROOT;
    else if( key_string == "INPUT_FILE_FORMAT" )
        return INPUT_FILE_FORMAT;
    else if( key_string == "NF" )
        return NF;
    else if( key_string == "COLS_FIELDS" )
        return COLS_FIELDS;
    else
    {
        std::cout << "Something wrong in cfg file" << std::endl;
        exit (EXIT_FAILURE);       
    }

}


void Read_cfg ( const std::string filename, prob_settings &settings )
{

    
    std::ifstream cFile ( filename );
    if ( cFile.is_open() ){

        size_t delimiterPos; 
        std::string name, value;
        std::string line;

        while( getline(cFile, line) )
        {
            line.erase( remove_if(line.begin(), line.end(), isspace),
                                 line.end() );  //include ::  in front of isspace if using namespace std
            if( line[0] == '#' || line.empty() )
                continue;
            delimiterPos = line.find("=");
            name = line.substr(0, delimiterPos);
            value = line.substr(delimiterPos + 1);

            switch (read_keyword_type(name))
            {

                case PROB_DIM:
                {
                    settings.dim_prob = value;
                    std::cout << "Problem dimension : " << value << std::endl;
                    break;
                }

                case NS:
                {
                    settings.Ns = std::stoi(value);
                    std::cout << "Number of snapshots : " << value << std::endl;
                    break;
                }

                case DS:
                {
                    settings.Ds = std::stod(value);
                    std::cout << "Delta between selected snapshots (Equispaced) : " << value << std::endl; 
                    break;
                }

                case EN:
                {
                    settings.En = std::stod(value);
                    std::cout << "Energy content used in the reconstruction : " << value << std::endl;
                    break;
                }

                case INPUT_FILE_ROOT:
                {
                    settings.in_file_root = value;
                    std::cout << "Input file root name : " << value << std::endl;
                    break;
                }

                case INPUT_FILE_FORMAT:
                {
                    settings.in_file_format = value;
                    std::cout << "Input file format : " << value << std::endl;
                    break;
                }

                case NF:
                {
                    settings.Nf = std::stoi(value);
                    std::cout << "Filter size for feature extraction : " << value << std::endl;
                    break;
                }

                case COLS_FIELDS:
                {
                    std::string str = value;
                    std::stringstream ss(str);

                    int i;

                    while ( ss >> i ) {

                        settings.Cols.push_back(i);

                        if (ss.peek() != ',' || ss.peek() != ' ')
                            ss.ignore();

                    }

                    std::cout << "Number of columns to process: \t";
                    
                    for ( i = 0; i < settings.Cols.size(); i++ )
                        std::cout << settings.Cols[i] << "\t";

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
    {
        std::cout << "Unable to open config file! Terminating ... " << '\n';
        exit (EXIT_FAILURE);
    }

}


int N_gridpoints( const std::string file_in) {

    std::string line_flow_data;
    std::ifstream flow_data;
    flow_data.open(file_in.c_str());

    if( !flow_data.is_open() )
    {

        std::cout << " While getting number of grid points, \nFile: " 
            << file_in << " not found" << std::endl;
        exit (EXIT_FAILURE);

    }

    int n_row = 0;


    while(getline( flow_data, line_flow_data )) 
    {
           
        if ( line_flow_data.compare(0,1,"A") == 0 ) //if reading SU2 native restart file
            break;
            
        n_row++;

   } 

    flow_data.close();

    return (n_row-1);

}





Eigen::MatrixXd read_col( std::string filename, int Nr, std::vector<int> Cols, Eigen::MatrixXd &field )
{


    std::ifstream flow_data;
    flow_data.open( filename );    

    if ( !flow_data.is_open() )
    {
        std::cout << "File : " << filename << " not found" << std::endl;    
        exit (EXIT_FAILURE);
    }

    std::string line_flow_data ;

    // Read row of headers
    getline( flow_data, line_flow_data );

    int n_row = 0;

    while ( getline( flow_data, line_flow_data ) )
    {

        Eigen::RowVectorXd point(Cols.size());
        std::istringstream iss(line_flow_data);
        std::string token;
        long double rubbish;
        int count = 0, c = 0; 

        if ( filename.compare( filename.size()-3, 3, "csv") == 0 ) 
        {

            while( getline( iss, token, ',') ) 
            {

                rubbish = std::stold(token);

                //This trick is due to some bad values I found in the restart SU2 file
                if ( (rubbish!=0.0) && ((std::abs(rubbish) < 1.0e-300) || (std::abs(rubbish) > 1.0e300)) )
                {
                    std::cout << " Rubbish value : " << std::setprecision(17) << std::scientific << 
                        rubbish <<  " on the row : "<< n_row << std::endl;
                    rubbish = 0.0;
                }
                
                if ( count == Cols[c])
                {

                    point(c) = rubbish;
                    c++;
                }
                count ++;
            }

        } else if ( filename.compare( filename.size()-3, 3, "dat") == 0 ) 
        {

            while( getline( iss, token, '\t') )
            {

                if ( token.compare(0,1,"A") == 0 ) //if reading SU2 restart file
                    break;

                rubbish = std::stold(token);
                

                //This trick is due to some bad values I found in the restart SU2 file
                if ( (rubbish!=0.0) && ((std::abs(rubbish) < 1.0e-300) || (std::abs(rubbish) > 1.0e300)))
                {
                    std::cout << " Rubbish value : " << std::setprecision(17) << std::scientific << 
                        rubbish <<  " on the row : "<< n_row << std::endl;
                    rubbish = 0.0;
                }

                if ( count == Cols[c])
                {

                    point(c) = rubbish;
                    c++;
                }
                count ++;

            } 

        } else
        {
            std::cout << "Bad file format" << std::endl;
            exit (EXIT_FAILURE);
        }

        if ( token.compare(0,1,"A") == 0 ) //if reading SU2 restart file
            break;

        field.row(n_row) = point; 
        n_row++;

    }

    flow_data.close();

    return field;
  
}