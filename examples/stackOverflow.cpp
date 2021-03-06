#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"



int main ( int argc, char *argv[] )
{

//Input string Method - TestNumber - Number of snapshots - rank - reconstruction

std::cout << "Main start" << std::endl;
std::string method = argv[1];
int test = std::atoi(argv[2]);
int Nsnap = std::atoi(argv[3]);
int r = std::atoi(argv[4]);
int fac_rec = std::atoi(argv[5]);

double pi = 3.14159;
Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(50,-5.0,5.0);
Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(50,-5.0,5.0);
Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(fac_rec*Nsnap,0.0,5.0);

//Generate gaussian noise
// const double mean = 0.0;
// const double stddev = 0.3;
// std::default_random_engine generator;
// std::normal_distribution<double> dist(mean, stddev);


int Ns = t.size();
int Npx = x.size();
int Npy = y.size();

double dt = t(2)-t(1);
double t_0 = t(0);

// double kl = 1.0, omega = 1.0;

//Define parameter of the first gaussian
double x01 = -4.0, y01 = -4.0;
double sigma_x1 = 0.5, sigma_y1 = 0.5;
double vel1 = 3.5; //Velocity or frequency
double A1 = 1.0;

//Define parameter of the second gaussian
double x02 = -4.0, y02 = 4.0;
double sigma_x2 = 0.5, sigma_y2 = 0.5;
double vel2 = 1.0;
double A2 = 1.0;

Eigen::MatrixXd f (Npx, Npy);
Eigen::MatrixXd data(Npx*Npy,Ns);


for ( int k = 0; k < Ns; k++ )
{

    for ( int i = 0; i < Npx; i++ )
    {
        for ( int j = 0; j < Npy; j++)
        {   
            switch ( test )
            {
                case 1: 
                {
                //double sin wave with noise
                //f(i,j) = std::sin(0.2*pi*x(i))*std::sin(0.5*t(j)) + dist(generator);
                
                //gaussian 1D traveling wave
                // f(i,j) = std::exp(-std::pow((x(i)-t(j)+5.0)/2.0, 2.0));
                
                //Sinusoidal 2D traveling wave
                //f(i,j) = std::sin(kl*(x(i)+y(j))-omega*t(k));

                //Two gaussian 2D traveling wave at different velocities
                    f(i,j) = A1*std::exp( -((x(i) - x01 - vel1*t(k))*(x(i) - x01 - vel1*t(k))/2.0/sigma_x1/sigma_x1 +
                        (y(j) - y01)*(y(j) - y01)/2.0/sigma_y1/sigma_y1)) + A2*std::exp( -((x(i) - x02)*(x(i) - x02)/2.0/sigma_x2/sigma_x2 +
                        (y(j) - y02 - vel2*t(k))*(y(j) - y02 - vel2*t(k))/2.0/sigma_y2/sigma_y2));
                break;
                }
                case 2:
                {
                //Two gaussian 2D traveling wave at different velocities ( non linear )
                    f(i,j) = A1*std::exp( -((x(i) - x01 - vel1*t(k)*t(k))*(x(i) - x01 - vel1*t(k)*t(k))/2.0/sigma_x1/sigma_x1 +
                        (y(j) - y01)*(y(j) - y01)/2.0/sigma_y1/sigma_y1)) + A2*std::exp( -((x(i) - x02)*(x(i) - x02)/2.0/sigma_x2/sigma_x2 +
                        (y(j) - y02 + vel2*std::sqrt(t(k)))*(y(j) - y02 + vel2*std::sqrt(t(k)))/2.0/sigma_y2/sigma_y2));

                break;
                }
                case 3:
                {
                //Two gaussian 2D oscillating wave with different frequencies
                    f(i,j) = A1*std::exp( -((x(i) - x01 - std::sin(vel1*t(k)))*(x(i) - x01 - std::sin(vel1*t(k)))/2.0/sigma_x1/sigma_x1 +
                        (y(j) - y01)*(y(j) - y01)/2.0/sigma_y1/sigma_y1)) + A2*std::exp( -((x(i) - x02)*(x(i) - x02)/2.0/sigma_x2/sigma_x2 +
                        (y(j) - y02 - std::sin(vel2*t(k)))*(y(j) - y02 - std::sin(vel2*t(k)))/2.0/sigma_y2/sigma_y2));
                break;
                }
                //One gaussian 2D oscillating wave with assigned frequency
                // f(i,j) = A1*std::exp( -((x(i) - x01 - std::sin(vel1*t(k)))*(x(i) - x01 - std::sin(vel1*t(k)))/2.0/sigma_x1/sigma_x1 +
                //         (y(j) - y01)*(y(j) - y01)/2.0/sigma_y1/sigma_y1));
                default:
                    break;
            }
        }
    }

    Eigen::Map<Eigen::RowVectorXd> ft(f.data(),f.size());
    data.col(k) = ft.transpose();

}

std::cout << "Writing original field ..." << std::endl; 

std::ofstream f_data;
f_data.open("Original_field.dat");  

for ( int k = 0; k < Npx*Npy; k++ )
{

    for( int j = 0; j < Ns; j++ ) 
        f_data << data(k,j) << " ";   

f_data << std::endl; 
}


f_data.close();

Eigen::VectorXcd lam;
Eigen::MatrixXcd eig_vec;
Eigen::VectorXd lam_POD;
Eigen::MatrixXd eig_vec_POD;


Eigen::VectorXd tnew(fac_rec*Nsnap); //times for reconstruction
tnew(0) = 0.0;
for ( int i = 1; i < fac_rec*Nsnap ; i ++ )
{
    tnew(i) = tnew(i-1) + dt/fac_rec;
}

Eigen::MatrixXcd Rec(Npx*Npy, tnew.size());

Eigen::MatrixXd sn_set(Npx*Npy, Nsnap);

for ( int i = 0; i < Nsnap; i ++ )
{
    sn_set.col(i) = data.col(i*fac_rec);

}


std::string flag_coef = "OPT";


if( method == "HODMD")
{

    Eigen::VectorXcd Coefs;
    double tol = 1e-16;
    int d = 5;
    Eigen::MatrixXcd Phi =  HODMD_basis( sn_set,
                                        lam,
                                        eig_vec,
                                        Coefs,
                                        tol,
                                        d );


    std::cout << "Computing Reconstruction... " << std::endl;

    Rec = TimeEvo_DMD ( tnew,
                        dt,
                        Coefs,
                        Phi,
                        lam );

    std::cout << "Done" << std::endl;
    


}


if( method == "mrDMD")
{
    int max_levels = 10;
    int max_cycles = 8;
    std::vector<node_mrDMD> nodes = {};
    nodes = mrDMD_basis( sn_set,       
                        nodes,        
                        r,                                  
                        t(1) - t(0), 0.0, 0, 0, 0, max_levels, max_cycles, flag_coef);

    std::cout << "Number of nodes : " << nodes.size() << std::endl << std::endl;

    //Reconstruction when only one node is present (standard DMD)
    // Eigen::MatrixXcd Rec = TimeEvo_DMD ( tnew,
                                        // dt,
                                        // nodes[0].Coefs,
                                        // nodes[0].Modes,
                                        // nodes[0].lam );
    
    std::cout << "Write coefs dynamics : " << std::endl;


    for ( int i = 0; i < max_levels; i ++ )
    {
        std::cout << "Level " << i << "\t";
        write_CoefsDynamics_mrDMD( nodes, i, tnew.size()/std::pow(2,i), 7);
        std::cout << "Done" << std::endl;
    }

    for ( int i = 0; i < tnew.size(); i++ )
    {
        Rec.col(i) = Reconstruction_mrDMD ( tnew(i), t(1) - t(0),
                                        nodes,
                                        "SCALAR" );

    }


    if ( fac_rec == 1 )
    {
        Eigen::MatrixXd Error_mrDMD_map = sn_set - Rec.real();
        Eigen::VectorXd Error_mrDMD = Eigen::VectorXd::Zero(sn_set.cols());

        for ( int i = 0; i < sn_set.cols(); i++ )
        {
            for ( int j = 0; j < sn_set.rows(); j++ )
            {
                Error_mrDMD(i) += Error_mrDMD_map(i,j)*Error_mrDMD_map(i,j);
            }
            Error_mrDMD(i) = std::sqrt(Error_mrDMD(i))/sn_set.col(i).norm();
            std::cout << "Error at time step " << t(i) << " : " << Error_mrDMD(i) << std::endl;
        }
    
    }

std::cout << " Computed Reconstruction" << std::endl;

}


if( method == "DMD" )
{
    Eigen::MatrixXcd Phi = DMD_basis( sn_set,
                                    lam,
                                    eig_vec,
                                    lam_POD,
                                    eig_vec_POD, 
                                    r);

    for ( int i = 0; i < Phi.cols(); i ++ )
    {
        Phi.col(i) = Phi.col(i)/Phi.col(i).norm();
    } 


    // std::cout << "Singular values:\n " << lam_POD << std::endl;
    // std::cout << "Eigen-values DMD : " << std::endl;

    // for ( int i = 0; i < lam.size(); i++)
    //     std::cout << std::log(lam(i))/dt << std::endl;

    //Calculate coefficients lsq
    // Eigen::VectorXcd b = Eigen::VectorXcd::Zero(data.rows()); 
    // for ( int k = 0; k < data.rows(); k++ )
    //     b(k).real(data(k,0));  
     
    // Eigen::VectorXcd alfa = Phi.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

    // std::cout << "Coefficients : " << alfa << std::endl;

    //Calculate optimized coefficients
    Eigen::VectorXcd alfa = Calculate_Coefs_DMD_exact (sn_set.leftCols(Nsnap-1), lam, Phi);

    Eigen::VectorXcd omega(lam.size());
    for ( int i = 0; i < lam.size(); i++ )
        omega(i) = std::log(lam(i))/dt;

    //Computing energy for each mode
    Eigen::VectorXd En = Eigen::VectorXd::Zero(Phi.cols());
    double T = t(t.size()-1);

    for ( int i = 0 ; i < Phi.cols(); i ++ )
    {

        double alfa_i = alfa(i).imag();
        double alfa_r = alfa(i).real();
        double sigma = omega(i).real();

        En(i) = (alfa_r*alfa_r + alfa_i*alfa_i)*(std::exp(2.0*sigma*T) - 1.0)/(2.0*sigma);

        std::cout << "En at mode " << i << " " <<  En(i) << std::endl;

    }

    dmd_sort( En, Phi, lam, alfa);

    std::cout << "En reordered :\n " << En << std::endl;
    std::cout << "Lam reordered : \n" << lam << std::endl;


// std::cout << "Done line 230" << std::endl;

    // std::cout << "DMD eigenvalues :\n" << lam << std::endl;
    // std::cout << "DMD omega :\n" << omega << std::endl;


    //Calculate Hybridized coefficients
    std::cout << "Calculating cofficients ... " << std::endl;
    Eigen::MatrixXcd Alfas = Calculate_Coefs_Matrix_DMD ( sn_set,
                                                        Phi,
                                                        omega,
                                                        t_0,
                                                        dt );
    
    std::cout << "Done" << std::endl; 

    std::cout << "Writing training points ..." << std::endl;

    std::ofstream train_real;
    train_real.open("train_real.dat");


    for ( int k = 0; k < Nsnap; k++ )
    {

        for( int j = 0; j < Alfas.cols(); j++ ) 
            train_real << Alfas(k,j).real() << " ";   

    train_real << std::endl;

    }

    train_real.close();


    std::ofstream train_imag;
    train_imag.open("train_imag.dat");


    for ( int k = 0; k < Nsnap; k++ )
    {

        for( int j = 0; j < Alfas.cols(); j++ ) 
            train_imag << Alfas(k,j).imag() << " ";   

    train_imag << std::endl;

    }

    train_imag.close();

    std::cout << "Done" << std::endl;

    std::vector<double> t_vec;
    for ( int i = 0; i < Nsnap; i++ )
        t_vec.push_back(t(i));

    // std::vector<rbf> surr_coefs_real = getSurrCoefs ( t_vec,
    //                                                 Alfas.real(),
    //                                                 "LINEAR" );


    // std::vector<rbf> surr_coefs_imag = getSurrCoefs ( t_vec,
    //                                                 Alfas.imag(),
    //                                                 "LINEAR" );

    // std::cout << "Writing file with interpolating coefficients..." << std::endl;

    // std::ofstream coefs_real;
    // coefs_real.open("Coefs__interp_real.dat");
    // std::ofstream coefs_imag;
    // coefs_imag.open("Coefs__interp_imag.dat");

    // double R, I;

    // for ( int k = 0; k < tnew.size(); k++ )
    // {
    //         std::vector<double> t(1,tnew(k));
    //     for ( int i = 0; i < Alfas.cols(); i++ )
    //     {
     
    //         surr_coefs_real[i].evaluate(t, R);
    //         surr_coefs_imag[i].evaluate(t, I);

    //         coefs_real << R << " ";
    //         coefs_imag << I << " ";

    //     } 
    //     coefs_real << std::endl;
    //     coefs_imag << std::endl;
    // }

    // coefs_real.close();
    // coefs_imag.close();

    // std::cout << "Done" << std::endl; 

    std::cout << "Computing Reconstruction... " << std::endl;
    // for ( int i = 0; i < tnew.size(); i++ )
    // {
    //     Rec.col(i) = Reconstruction_Hybrid_DMD ( tnew(i),
    //                                             t_vec,
    //                                             Alfas,
    //                                             Phi,
    //                                             omega,
    //                                             "SCALAR",
    //                                             "LINEAR" );
    // }


    Rec = TimeEvo_DMD ( tnew,
                        dt,
                        alfa,
                        Phi,
                        lam );
    std::cout << "Done" << std::endl;
}


// int Nm;

// if ( r == 0 )
//     Nm = SVHT ( lam_POD, Ns, Np );
// else 
//     Nm = r;


// std::cout << "Writing Modes ..." << std::endl;

// std::ofstream dmd_data;
// dmd_data.open("Modes_DMD.dat");

// for ( int k = 0; k < Np; k++ )
// {

//     for( int j = 0; j < Phi.cols(); j++ ) 
//         dmd_data << Phi(k,j).real() << " ";   
    
// dmd_data << std::endl;

// }

// dmd_data.close();

// std::ofstream dmdi_data;
// dmdi_data.open("ModesI_DMD.dat");

// for ( int k = 0; k < Np; k++ )


std::cout << "Writing Reconstruction ..." << std::endl;

std::ofstream dmd_rec;
dmd_rec.open("Rec_DMD.dat");


for ( int k = 0; k < Npx*Npy; k++ )
{

    for( int j = 0; j < tnew.size(); j++ ) 
        dmd_rec << Rec(k,j).real() << " ";   
    
dmd_rec << std::endl;

}

dmd_rec.close();

std::cout << "Done" << std::endl;

return 0;
}