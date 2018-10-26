#include <iostream>
#include <stdio.h> 
#include <cmath>

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




return 0;

}