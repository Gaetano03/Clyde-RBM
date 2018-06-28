#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"


int main(int argc, char *argv[]) {

    std::cout << "-----------RBM-Clyde start-------------" << std::endl << std::endl;

    std::string filecfg = argv[1];
    std::string prob_dim;
    std::string file_temp; 
    std::string input_format, input_root;
    std::vector<int> Col_fields;
    int Ns, delta_s, Nf;
    double Ec;

    prob_settings settings;
    Read_cfg( filecfg, settings );


    // prob_dim = settings.dim_prob;
    // Ns = settings.Ns;
    // Nf = settings.Nf;
    // delta_s = settings.Ds;
    // Ec = settings.En;
    // Col_fields = settings.Cols;

    // input_format = settings.in_file_format;
    // input_root = settings.in_file_root;

    // Create_snapset sn_set ( Ns, delta_s, Col_fields, input_root,
    //                          prob_dim );

    // int Nr = sn_set.get_gridpoints();

    // std::cout << "Number of grid points : " << Nr << std::endl;


    // Eigen::MatrixXd snapM = sn_set.generate_snap_matrix();

    // basis_extraction Base_set ( Ns, Nr, Nf, sn_set);


    return 0;
}
