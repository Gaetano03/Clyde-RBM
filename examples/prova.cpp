#include <iostream>
#include <stdio.h> 
#include <cmath>
#include <vector>
#include "LinearAlgebra/Eigen/Dense"
#include "LinearAlgebra/Eigen/Eigenvalues"


int main()
{

int step = 4;
int j = 2;

std::cout << "(int)1/(int)step : " << 1/step << std::endl;
std::cout << "1.0/(double)step : " << 1.0/(double)step << std::endl;

double aa = 1.0/(double)step;

std::cout << "ans*(int)j : " << aa*j << std::endl;
std::cout << "ans*(double)j : " << aa*(double)j << std::endl;
std::cout << "floor without double : " << std::floor(50/16) << std::endl;
std::cout << "floor with double : " << std::floor((double)50/(double)16) << std::endl;

Eigen::MatrixXd m(3,3);
m << -1.0, 5.0, 3.2,
    9.1, -2.0, 5.5,
    -3.5, -7.1, 4.6;
    // 2.0, 8.0, 6.3;

Eigen::VectorXd v(3);
v << -1.0, 2.5, -9.3;

Eigen::MatrixXd cc = Eigen::MatrixXd::Zero(5,5);

cc = m;

std::cout << "Matrix m:\n" << cc << std::endl; 


std::cout << "abs of vector v" << v.cwiseAbs() << std::endl;

Eigen::BDCSVD<Eigen::MatrixXd> svd( m, Eigen::ComputeThinU | Eigen::ComputeThinV );
Eigen::VectorXd s = svd.singularValues();

Eigen::MatrixXd V = svd.matrixV();
Eigen::MatrixXd U = svd.matrixU();
Eigen::MatrixXd Sigma = Eigen::MatrixXd::Zero(s.size(),s.size());

for ( int i = 0; i < s.size(); i++ )
    Sigma(i,i) = s(i);

std::cout << "Matrix Rec= USV :\n " << U*Sigma*V << std::endl; 
std::cout << "Matrix Rec= USV^T :\n " << U*Sigma*V.transpose() << std::endl; //This is the one that works


std::vector<Eigen::MatrixXd> pul(4, Eigen::MatrixXd::Zero(5,5));


pul[0](1,2) = 9.8;

std::cout << " Element : " << pul[0](1,2) << std::endl;

// Eigen::VectorXd s = Eigen::VectorXd::LinSpaced(10,0.0001,0.0000000001);
// 
// std::cout << " s : " << s << std::endl;
// 
// int N_red1 = 0;
// double eps1 = 0.0;
// double tol = 1e-5;

// while ( eps1 < tol )
// {
    // eps1 = s.tail(N_red1 + 1).sum()/s.sum();
    // std::cout << " eps1 at iteration " << N_red1 << " : " << eps1 << std::endl;
    // N_red1++;

// }
// 
// std::cout << "N_red1 : " << N_red1-1 << std::endl;

Eigen::VectorXd u(10);
u << 2.0, 3.5, 9.0, 5.1, 0.5, 5.1, 3.6, 2.2, 5.1, 9.9;

int pos;
double min_val = u.minCoeff( &pos );

std::cout << " Min val= " << min_val << " at index : " << pos << std::endl;

std::cout << "Sciao sciao " << std::endl;


return 0;

}