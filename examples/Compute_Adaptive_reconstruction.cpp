#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"


int main( int argc, char *argv[] )
{
    std::cout << "Adaptive Reconstruction RBM-Clyde starts " << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];
    
    //Reading configuration file
    Read_cfg( filecfg, settings );

    int s_Nf = 5;
    std::vector<Eigen::MatrixXd> Err_RBM_Nm_time;
    std::vector<Eigen::MatrixXd> J_RBM_Nm_time;
    std::vector<Eigen::MatrixXd> EJ_RBM_Nm_time;

    std::string filename1, filename2;
    for ( int i = 0; i < s_Nf; i++ )
    {
        filename1 = "Error_SPOD_Nf_" + std::to_string(i) + ".dat";
        Eigen::MatrixXd Err_map = read_err_j ( filename1, settings.Ns );

        filename2 = "Jaccard_SPOD_Nf_" + std::to_string(i) + ".dat";
        Eigen::MatrixXd J_map = read_err_j ( filename2, settings.Ns );

        Err_RBM_Nm_time.push_back(Err_map);
        J_RBM_Nm_time.push_back(J_map);
        J_RBM_Nm_time.push_back(Err_map.cwiseQuotient(J_map));

    }

    std::cout << "Test stampa error POD:\n" << Err_RBM_Nm_time[0] << std::endl;
    std::cout << "Test stampa J POD:\n" << J_RBM_Nm_time[0] << std::endl;

    return 0;
}






