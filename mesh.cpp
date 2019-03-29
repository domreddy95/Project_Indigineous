#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <random>

class mesh
{
private:
    /* data */
public:
    double charge;
    double nodal_volume;
    double charge_density;
    double potential;
    std::vector <double> electric_feild;
    std::vector <double> magnetic_feild;
    mesh();
    ~mesh();
};

mesh::mesh()
{
    std::vector <double> electric_feild (3,0);
    std::vector <double> magnetic_feild (3,0);
    this->charge = 0;
    this->nodal_volume = 0;
    this->charge_density = 0;
    this->potential = 0;
    this->electric_feild = electric_feild;
    this->magnetic_feild = magnetic_feild;
}

mesh::~mesh()
{
}

std::vector< std::vector<mesh> > CreateMesh(double min_dim1, double max_dim1, double min_dim2, double max_dim2, int num1, int num2)
{
    std::vector<mesh> temp_array;
    temp_array.reserve(num2+1);
    std::vector< std::vector<mesh> > temp_array2;
    temp_array2.reserve(num1+1);
    for (int i=0;i==num1;i++)
    {
        for (int j=0;j==num2;j++)
        {
            temp_array.emplace_back(mesh());
        }
        temp_array2.emplace_back(temp_array);
    }

    return temp_array2;
}

void calc_nodal_volume(double min_dim1, double max_dim1, double min_dim2, double max_dim2)
{
    
}