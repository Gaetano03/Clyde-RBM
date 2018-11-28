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
        EJ_RBM_Nm_time.push_back(Err_map.cwiseQuotient(J_map));

    }

    filename1 = "Error_DMD.dat";
    Eigen::MatrixXd Err_map = read_err_j ( filename1, settings.Ns );

    filename2 = "Jaccard_DMD.dat";
    Eigen::MatrixXd J_map = read_err_j ( filename2, settings.Ns );

    Err_RBM_Nm_time.push_back(Err_map);
    J_RBM_Nm_time.push_back(J_map);
    EJ_RBM_Nm_time.push_back(Err_map.cwiseQuotient(J_map));

    filename1 = "Error_RDMD.dat";
    Err_map = read_err_j ( filename1, settings.Ns );

    filename2 = "Jaccard_RDMD.dat";
    J_map = read_err_j ( filename2, settings.Ns );

    Err_RBM_Nm_time.push_back(Err_map);
    J_RBM_Nm_time.push_back(J_map);
    EJ_RBM_Nm_time.push_back(Err_map.cwiseQuotient(J_map));      

//Adaptive reconstruction on each selected time step
    int best_method_idx;

    std::cout << "Initializing Vector of time ... " << std::endl; 
    Eigen::VectorXd t_vec( settings.Ns );
    t_vec(0) = settings.nstart*settings.Dt_cfd;
    for ( int i = 1; i < settings.Ns; i++ )
        t_vec(i) = t_vec(i-1) + settings.Dt_cfd*settings.Ds;
    
    double tol = settings.tol;
    int index1, index2;
    Eigen::MatrixXd Err_interp(EJ_RBM_Nm_time.size(), settings.Ns);

    for ( int i = 0; i < settings.t_rec.size(); i++ )
    {
        std::vector<int> pos = {};
        std::cout << " Adaptive reconstruction at time : " << settings.t_rec[i] << std::endl;

        index1 = 0;
        index2 = 0;
        for ( int nt = 0; nt < t_vec.size()-1; nt ++ )
        {
            if ( (settings.t_rec[i] >= t_vec(nt)) && (settings.t_rec[i] <= t_vec(nt+1)) )
            {
                index1 = nt;
                index2 = nt+1;
                break;
            }
        }

        if ( index1 == index2 )
        {
            std::cout << "Time for reconstruction out of interval!" << std::endl;
            continue;
        }

        int count = 0;
        for ( int k = 0; k < EJ_RBM_Nm_time.size(); k ++ )
        {
            for ( int nm = 0; nm < settings.Ns; nm ++)
                Err_interp(k,nm) = EJ_RBM_Nm_time[k](index1,nm) + (EJ_RBM_Nm_time[k](index2,nm) - EJ_RBM_Nm_time[k](index1,nm))/
                                    (settings.Dt_cfd*settings.Ds)*(settings.t_rec[i] - t_vec[index1]);
        }

        double eps = 1.0;
        int ncount = 0;

        while ( eps > tol )
        {
            eps = Err_interp.col(ncount).minCoeff( &best_method_idx );
            ncount++;
        }

        std::cout << "Best method is " << best_method_idx << " with number of modes " << ncount << std::endl;
        
    }


    return 0;
}






