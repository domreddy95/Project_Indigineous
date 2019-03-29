#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <random>
#include <time.h>
#include <fstream>
#include <iomanip>



int main()
{
    double a = 1.0381083198731098371927310937;
    system("mkdir Output");
    system("mkdir Output/Particles");
    system("mkdir Output/Charge");
    system("mkdir Output/ElectricFeild");
    system("mkdir Output/PotentialFeild");

    std::ofstream write;
    write.open ("Output/Particles/phase.dat");
    write << std::setprecision(2) <<  a << " " << "2\n";
    write << a << " " << "2\n";
    write << a << " " << "2\n";
    write << a << " " << "2\n";
    write << a << " " << "2\n";
    write << a << " " << "2\n";
    write.close();
    return 0;
}