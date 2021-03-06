
#include "read_Inputs.hpp"

keywords read_keyword_type( const std::string &key_string )
{

    if( key_string == "NS" )
        return NS;
    else if( key_string == "DS" )
        return DS;
    else if( key_string == "EN" )
        return EN;
    else if( key_string == "SIGMA" )
        return SIGMA;
    else if( key_string == "NSTART" )
        return NSTART;
    else if( key_string == "NDIM" )
        return NDIM;
    else if( key_string == "DT_CFD" )
        return DT_CFD;
    else if( key_string == "FLAG_DIM" )
        return FLAG_DIM;
    else if( key_string == "FLAG_PROB" )
        return FLAG_PROB;
    else if( key_string == "INPUT_FILE" )
        return INPUT_FILE;
    else if( key_string == "OUTPUT_FILE" )
        return OUTPUT_FILE;
    else if( key_string == "NF" )
        return NF;
    else if( key_string == "COLS_COORDS" )
        return COLS_COORDS;
    else if( key_string == "COLS_FIELDS" )
        return COLS_FIELDS;
    else if( key_string == "FLAG_METHOD" )
        return FLAG_METHOD;
    else if( key_string == "FLAG_MEAN" )
        return FLAG_MEAN;
    else if( key_string == "FLAG_BC" )
        return FLAG_BC;
    else if( key_string == "FLAG_FILTER" )
        return FLAG_FILTER;
    else if( key_string == "FLAG_WDB_BE" )
        return FLAG_WDB_BE;
    else if( key_string == "FLAG_REC" )
        return FLAG_REC;
    else if( key_string == "FLAG_INTERP" )
        return FLAG_INTERP;
    else if( key_string == "T_REC" )
        return T_REC;
    else if( key_string == "RANK_RDMD" )
        return RANK_RDMD;
    else if( key_string == "RANK" )
        return RANK;
    else if( key_string == "HO_D" )
        return HO_D;
    else if( key_string == "DMD_COEF_FLAG" )
        return DMD_COEF_FLAG;
    else if( key_string == "MAX_CYCLES" )
        return MAX_CYCLES;
    else if( key_string == "MAX_LEVELS" )
        return MAX_LEVELS;
    else if( key_string == "TOL" )
        return TOL;
    else
    {
        std::cout << "Something wrong in cfg file" << std::endl;
        exit (EXIT_FAILURE);       
    }

}


void Read_cfg ( const std::string filename, prob_settings &settings )
{

    
    std::ifstream cFile ( filename );
    if ( cFile.is_open() )
    {

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

                case FLAG_PROB:
                {
                    settings.flag_prob = value;
                    //std::cout << "Problem flag : " << value << std::endl;
                    break;
                }

                case FLAG_DIM:
                {
                    settings.flag_dim = value;
                    //std::cout << "Dimension flag : " << value << std::endl;
                    break;
                }

                case NS:
                {
                    settings.Ns = std::stoi(value);
                    //std::cout << "Number of snapshots : " << value << std::endl;
                    break;
                }

                case DS:
                {
                    settings.Ds = std::stod(value);
                    //std::cout << "Delta between selected snapshots (Equispaced) : " << value << std::endl; 
                    break;
                }

                case EN:
                {
                    settings.En = std::stod(value);
                    //std::cout << "Energy content used in the reconstruction : " << value << std::endl;
                    break;
                }

                case SIGMA:
                {
                    settings.sigma = std::stod(value);
                    //std::cout << "Sigma for SPOD gaussian filter : " << value << std::endl;
                    break;
                }

                case TOL:
                {
                    settings.tol = std::stod(value);
                    //std::cout << "Sigma for SPOD gaussian filter : " << value << std::endl;
                    break;
                }                

                case NSTART:
                {
                    settings.nstart = std::stoi(value);
                    //std::cout << "Initial snapshot number : " << value << std::endl;
                    break;
                }

                case NDIM:
                {
                    settings.ndim = std::stoi(value);
                    //std::cout << "Initial snapshot number : " << value << std::endl;
                    break;
                }
                
                case DT_CFD:
                {
                    settings.Dt_cfd = std::stod(value);
                    //std::cout << "Dt used in CFD simulation : " << value << std::endl;
                    break;
                }

                case INPUT_FILE:
                {
                    settings.in_file = value;
                    //std::cout << "Input file root name and format : " << value << std::endl;
                    break;
                }

                case OUTPUT_FILE:
                {
                    settings.out_file = value;
                    //std::cout << "Output file root name and format : " << value << std::endl;
                    break;
                }

                case NF:
                {
                    settings.Nf = std::stoi(value);
                    //std::cout << "Filter size for feature extraction : " << value << std::endl;
                    break;
                }

                case FLAG_METHOD:
                {
                    settings.flag_method = value;
                    //std::cout << "Method for feature extraction : " << value << std::endl;
                    break;
                }

                case FLAG_MEAN:
                {
                    settings.flag_mean = value;
                    //std::cout << "Mean subtraction : " << value << std::endl;
                    break;
                }

                case FLAG_BC:
                {
                    settings.flag_bc = value;
                    //std::cout << "Boundary consition for correlation matrix : " << value << std::endl;
                    break;
                }

                case FLAG_FILTER:
                {
                    settings.flag_filter = value;
                    //std::cout << "Filter type for SPOD : " << value << std::endl;
                    break;
                }

                case FLAG_WDB_BE:
                {
                    settings.flag_wdb_be = value;
                    //std::cout << "Write database basis extraction (modes and coefficients) : " << value << std::endl;
                    break;
                }

                case FLAG_REC:
                {
                    settings.flag_rec = value;
                    //std::cout << "Compute reconstructed field : " << value << std::endl;
                    break;
                }

                case FLAG_INTERP:
                {
                    settings.flag_interp = value;
                    //std::cout << "Interpolation technique for rbf : " << value << std::endl;
                    break;
                }

                case RANK:
                {
                    settings.r = std::stoi(value);
                    //std::cout << "DMD rank : " << value << std::endl;
                    break;
                }

                case RANK_RDMD:
                {
                    settings.r_RDMD = std::stoi(value);
                    //std::cout << "Recursive-DMD rank : " << value << std::endl;
                    break;
                }

                case DMD_COEF_FLAG:
                {
                    settings.dmd_coef_flag = value;
                    //std::cout << "DMD coefs method : " << value << std::endl;
                    break;
                }

                case HO_D:
                {
                    settings.d = std::stoi(value);
                    //std::cout << "DMD coefs method : " << value << std::endl;
                    break;
                }

                case MAX_CYCLES:
                {
                    settings.max_cycles = std::stoi(value);
                    //std::cout << "Max_cycles for mrDMD : " << value << std::endl;
                    break;
                }

                case MAX_LEVELS:
                {
                    settings.max_levels = std::stoi(value);
                    //std::cout << "Max_levels for mrDMD : " << value << std::endl;
                    break;
                }

                case T_REC:
                {
                    std::string str = value;
                    std::stringstream ss(str);

                    double i;

                    while ( ss >> i ) {

                        settings.t_rec.push_back(i);

                        if (ss.peek() != ',' || ss.peek() != ' ')
                            ss.ignore();

                    }

                    //std::cout << "Times desired for reconstruction: \t";
                    
                    //for ( i = 0; i < settings.t_rec.size(); i++ )
                        //std::cout << settings.t_rec[i] << "\t";

                    //std::cout << std::endl;

                    break;
                }

                case COLS_COORDS:
                {
                    std::string str = value;
                    std::stringstream ss(str);

                    int i;

                    while ( ss >> i ) {

                        settings.Cols_coords.push_back(i);

                        if (ss.peek() != ',' || ss.peek() != ' ')
                            ss.ignore();

                    }

                    //std::cout << "Number of columns with coordinates: \t";
                    
                    //for ( i = 0; i < settings.Cols_coords.size(); i++ )
                        //std::cout << settings.Cols_coords[i] << "\t";

                    ///std::cout << std::endl;

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

                    //std::cout << "Number of columns to process: \t";
                    
                    //for ( i = 0; i < settings.Cols.size(); i++ )
                        //std::cout << settings.Cols[i] << "\t";

                    //std::cout << std::endl;

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
            << file_in << " not found " << std::endl;
        exit (EXIT_FAILURE);

    }

    int n_row = 0;


    while(getline( flow_data, line_flow_data )) 
    {
           
        if ( line_flow_data.compare(0,1,"E") == 0 || line_flow_data.compare(0,1,"A") == 0 ) //if reading SU2 native restart file
            break;
            
        n_row++;

   } 

    flow_data.close();

    return (n_row-1);

}




Eigen::MatrixXd read_col( std::string filename, int Nr, std::vector<int> Cols )
{


    Eigen::MatrixXd field (Nr, Cols.size());
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

                if ( token.compare(0,1,"E") == 0 || token.compare(0,1,"A") == 0 ) //if reading SU2 restart file
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

        if ( token.compare(0,1,"E") == 0 || token.compare(0,1,"A") == 0 ) //if reading SU2 restart file
            break;

        field.row(n_row) = point; 
        n_row++;

    }

    flow_data.close();

    return field;
  
}


Eigen::MatrixXd read_modes( std::string filename, int Nr, int r_RDMD )
{
    Eigen::MatrixXd f(Nr, r_RDMD);

    std::ifstream Modes_data;
    Modes_data.open( filename );

        if ( !Modes_data.is_open() )
    {
        std::cout << "File : " << filename << " not found" << std::endl;    
        exit (EXIT_FAILURE);
    }

    std::string line_flow_data ;

    // Read row of headers
    getline( Modes_data, line_flow_data );

    int n_row = 0;

    while ( getline( Modes_data, line_flow_data ) )
    {

        Eigen::RowVectorXd point(r_RDMD);
        std::istringstream iss(line_flow_data);
        std::string token;
        double rubbish;
        int count = 0, c = 0; 

        while( getline( iss, token, ' ') )
        {
            rubbish = std::stod(token);
            if ( count < r_RDMD )
                point(count) = rubbish;
            
            count ++;
        } 

        f.row(n_row) = point; 
        n_row++;

    }

    Modes_data.close();

    return f;

}


Eigen::MatrixXd read_coefs( std::string filename, int Ns, int r_RDMD )
{

    Eigen::MatrixXd f(r_RDMD, Ns);

    std::ifstream Coefs_data;
    Coefs_data.open( filename );

    if ( !Coefs_data.is_open() )
    {
        std::cout << "File : " << filename << " not found" << std::endl;    
        exit (EXIT_FAILURE);
    }

    std::string line_flow_data ;

    // Read row of headers
    getline( Coefs_data, line_flow_data );

    int n_row = 0;

    int count = 0;
    while ( getline( Coefs_data, line_flow_data ) && n_row < r_RDMD )
    {

        Eigen::RowVectorXd point(Ns);
        std::istringstream iss(line_flow_data);
        std::string token;
        double rubbish;
        int count = 0, c = 0; 

        while( getline( iss, token, ' ') && count < Ns )
        {
            rubbish = std::stod(token);
            point(count) = rubbish;
            
            count ++;
        } 

        f.row(n_row) = point; 
        n_row++;

    }

    Coefs_data.close();

    return f;
}



Eigen::MatrixXd read_err_j ( std::string filename, int Ns )
{
        std::ifstream file_data;
        file_data.open( filename );

            if ( !file_data.is_open() )
        {
            std::cout << "File : " << filename << " not found" << std::endl;    
            exit (EXIT_FAILURE);
        }

        std::string line_flow_data ;

        int n_row = 0, count = 0;
        Eigen::MatrixXd Err_map = Eigen::MatrixXd::Zero(Ns, Ns); 

        while ( getline( file_data, line_flow_data ) )
        {

            
            std::istringstream iss(line_flow_data);
            std::string token;
            double err;
            count = 0; 

            while( getline( iss, token, '\t') )
            {
                err = std::stod(token);
                Err_map(n_row, count) = err;

                count ++;
            } 
 
            n_row++;

        }

        file_data.close();

        return Err_map;

}