#include "Create_snapset.hpp"


Create_snapset::Create_snapset(const int &Ns, 
                            const int &ds,
                            const std::vector<int> Cols,
                            const std::string filename,
                            const std::string flag1,
                            const std::string flag2,
                            const std::string input_format ) {

                                m_Ns = Ns;
                                m_ds = ds;
                                m_root_inputfile = filename;
                                m_input_format = input_format;
                                m_Cols = Cols;
                                m_flag1 = flag1;
                                m_flag2 = flag2;
                                set_gridpoints();

                            }

Create_snapset::~Create_snapset() {
}


void Create_snapset::set_gridpoints() {

    std::string line_flow_data;
    std::string file_in = m_root_inputfile + "00000" + m_input_format;
    std::ifstream flow_data;
    flow_data.open(file_in.c_str());

    if(!flow_data.is_open()){

        std::cout << " While getting number of grid points, \nFile: " 
            << file_in << " not found" << std::endl;
        exit (EXIT_FAILURE);

    }

    int n_row = 0;


    while(getline( flow_data, line_flow_data )) {
           
        if ( line_flow_data.compare(0,1,"A") == 0 ) //if reading SU2 native restart file
            break;
            
        n_row++;

   } 

    flow_data.close();

    m_Nr = n_row-1;

}


int Create_snapset::get_gridpoints() {
    return m_Nr;
}


Eigen::MatrixXd Create_snapset::read_col( std::string filename ) {

    std::ifstream flow_data;
    flow_data.open( filename );    
    Eigen::MatrixXd field (m_Nr, m_Cols.size());

    if (!flow_data.is_open()){
        std::cout << "File : " << filename << " not found" << std::endl;    
        exit (EXIT_FAILURE);
    }

    std::string line_flow_data ;

    // Read row of headers
    getline( flow_data, line_flow_data );

    int n_row = 0;

    while ( getline( flow_data, line_flow_data ) ){

        Eigen::RowVectorXd point(m_Cols.size());
        std::istringstream iss(line_flow_data);
        std::string token;
        long double rubbish;
        int count = 0, c = 0; 

        if ( m_input_format == ".csv" ) {

            while( getline( iss, token, ',') ) {

                rubbish = std::stold(token);

                //This trick is due to some bad values I found in the restart SU2 file
                if ( (rubbish!=0.0) && ((std::abs(rubbish) < 1.0e-300) || (std::abs(rubbish) > 1.0e300))){
                    std::cout << " Rubbish value : " << std::setprecision(17) << std::scientific << 
                        rubbish <<  " on the row : "<< n_row << std::endl;
                    rubbish = 0.0;
                }
                
                if ( count == m_Cols[c]){

                    point(c) = rubbish;
                    c++;
                }
                count ++;
            }

        }else if ( m_input_format == ".dat" ) {

            while( getline( iss, token, '\t') ){

                if ( token.compare(0,1,"A") == 0 ) //if reading SU2 restart file
                    break;

                rubbish = std::stold(token);
                

                //This trick is due to some bad values I found in the restart SU2 file
                if ( (rubbish!=0.0) && ((std::abs(rubbish) < 1.0e-300) || (std::abs(rubbish) > 1.0e300))){
                    std::cout << " Rubbish value : " << std::setprecision(17) << std::scientific << 
                        rubbish <<  " on the row : "<< n_row << std::endl;
                    rubbish = 0.0;
                }

                if ( count == m_Cols[c]){

                    point(c) = rubbish;
                    c++;
                }
                count ++;

            } 

        }else{
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



Eigen::MatrixXd Create_snapset::generate_snap_matrix() {

    Eigen::MatrixXd field(m_Nr, m_Cols.size());
    std::string file_temp;
    int k = 0;

    if ( m_flag2 == "VELOCITY"){

        if ( m_flag1 == "2D"){
            
            Eigen::MatrixXd snap(2*m_Nr, m_Ns);

            for( int i = 0; i < m_Ns*m_ds; i += m_ds ){

                std::stringstream buffer;
                buffer << std::setfill('0') << std::setw(5) << std::to_string(i);
                file_temp = m_root_inputfile + buffer.str() + m_input_format;
                std::cout << "Reading fields from : " << file_temp << "\t";
                field = read_col(file_temp);
                std::cout << "Complete!" << std::endl;

                Eigen::VectorXd rho = field.col(0);
                Eigen::VectorXd rho_u = field.col(1);
                Eigen::VectorXd rho_v = field.col(2);
                Eigen::VectorXd u = rho_u.cwiseQuotient(rho);
                Eigen::VectorXd v = rho_v.cwiseQuotient(rho);

                snap.col(k) << u,
                                v;
                
                k++;
            }

            return snap;

        }else if ( m_flag1 == "3D"){

            Eigen::MatrixXd snap(3*m_Nr, m_Ns);

            for( int i = 0; i < m_Ns*m_ds; i += m_ds ){

                std::stringstream buffer;
                buffer << std::setfill('0') << std::setw(5) << std::to_string(i);
                file_temp = m_root_inputfile + buffer.str() + m_input_format;
                std::cout << "Reading fields from : " << file_temp << "\t";
                field = read_col(file_temp);
                std::cout << "Complete!" << std::endl;

                Eigen::VectorXd rho = field.col(0);
                Eigen::VectorXd rho_u = field.col(1);
                Eigen::VectorXd rho_v = field.col(2);
                Eigen::VectorXd rho_w = field.col(3);
                Eigen::VectorXd u = rho_u.cwiseQuotient(rho);
                Eigen::VectorXd v = rho_v.cwiseQuotient(rho);
                Eigen::VectorXd w = rho_w.cwiseQuotient(rho);

                snap.col(k) << u,
                                v,
                                w;

                k++;

        }

        return snap;

        }else {

            std::cout << " Set well problem dimension! Exiting ... " << std::endl;
            exit (EXIT_FAILURE);
            
        }

    }

}
    
