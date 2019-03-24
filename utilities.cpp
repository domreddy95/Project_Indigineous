#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <random>

class particle;
class species;// public particle;
class mesh;


class particle
{
public:
    std::vector <double> position;
    std::vector <double> velocity;
    std::vector <double> GetPosition(){return position;}
    void SetPosition(std::vector <double> position){position = position;}
    std::vector <double> GetVelocity(){return velocity;}
    void SetVelocity(std::vector <double> velocity){velocity = velocity;}
    particle();
    particle(double, double);
    ~particle();  
};

//Define default Constructer
particle::particle()
{
    std::vector <double> position (2,0.0);
    std::vector <double> velocity (2,0.0);
    this->position = position;
    this->velocity = velocity;
    //std::cout << " A particle is created" << std::endl;
}

//Parameterized Constructor
particle::particle(double charge, double mass)
{
    std::vector <double> position (2,0.0);
    std::vector <double> velocity (2,0.0);
    std::cout << " A parametrized particle is created" << std::endl;
}

particle::~particle()
{
    //std::cout << " A particle is destroyed";
}

class species: public particle
{
private:
    std::string species_name;
    int species_identifier;
    double charge;
    double mass;

public:
    species();
    species(std::string, int, double, double);
    ~species();
    double GetCharge(){return charge;}
    void SetCharge(double charge){this->charge=charge;}
    double GetMass(){return mass;}
    void SetMass(double mass){this->mass = mass;}
    void print_properties();
};

species::species() : particle()
{
    this->species_name = "no name";
    this->species_identifier = 0;
    this->charge = 0;
    this->mass = 0;
    //std::cout << " A species particle is created" << std::endl;
};


species::species(std::string species_name , int species_identifier, double charge,
double mass) : particle()
{
    this->species_name = species_name;
    this->species_identifier = species_identifier;
    this->charge = charge;
    this->mass = mass;
    //std::cout << " A parametrized species is created" << std::endl;
}

species::~species()
{
    //std::cout << " A species particle is destroyed" << std::endl;
}

void species::print_properties()
{
    using std::cout;
    using std::endl;
    cout << "charge = " << this->charge << endl;
    cout << "mass = " << this->mass << endl;
    cout << "position_x = " << position[0] <<endl;
    cout << "velocity_x = " << velocity[0] <<endl;
}

class mesh
{
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

std::vector <species> create_species (double num_of_particles = 0)
{
    std::vector <species> temp_array;
    temp_array.reserve(num_of_particles);
    for (int i = 0; i<num_of_particles; i++ )
    {
        temp_array.emplace_back(species());
    }
    return temp_array;
}

std::vector <species> load_species (double min_dim1, double max_dim1, double min_dim2, double max_dim2, double density)
{
    double delta_dim1 = (max_dim1-min_dim1);
    double delta_dim2 = (max_dim2-min_dim2);
    int num_of_particles = delta_dim1*((max_dim2*max_dim2)-(min_dim2*min_dim2))*density;
    int N = 1000;
    double dy = (max_dim2-min_dim2)/N;
    std::vector <species> temp_array1 = create_species(num_of_particles);
    int counter1 = 0;
    int counter2 = 0;
    for (int i = 0; i<N; i++)
    {
        double r1 = (i*dy) + min_dim2;
        double r2 = ((i+1)*dy) + min_dim2;
        int num = ((r2*r2) - (r1*r1))*delta_dim1*density;
        counter2 += num;
        for (int j = counter1; j<counter2; j++)
        {
            temp_array1[j].position[0] =  min_dim1 + ((rand() / (RAND_MAX + 1.))*delta_dim1);
            temp_array1[j].position[1] = r1 + ((rand() / (RAND_MAX + 1.))*dy);
        }
        counter1 = counter2;   
    }
    //std::cout << counter << std::endl;
    if (counter2 < num_of_particles)
    {
        int num = num_of_particles - counter2;
        counter2 += num;
        for (int i = counter1; i < counter2; i++)
        {
            temp_array1[i].position[0] =  min_dim1 + ((rand() / (RAND_MAX + 1.))*delta_dim1);
            temp_array1[i].position[1] = (0.70*delta_dim2) + ((rand() / (RAND_MAX + 1.))*(0.29*delta_dim2)); 
        }
    }
    
    return temp_array1;
}

void set_thermal_velocity(std::vector <species> &species_array, double ThermalVelocity = 0, int seed = 4444)
{   
    std::default_random_engine e(seed);
    std::normal_distribution <double> distrN (1.0,1.0);
    int N = species_array.size();
    for (int i=0; i<N; i++)
    {
        species_array[i].velocity[0] = distrN(e)*ThermalVelocity;
        species_array[i].velocity[1] = distrN(e)*ThermalVelocity;
    }
}

void set_beam_velocity(std::vector <species> &species_array, std::vector <double> &BeamVelocity)
{   
    int N = species_array.size();
    for (int i=0; i<N; i++)
    {
        species_array[i].velocity[0] += BeamVelocity[0];
        species_array[i].velocity[1] += BeamVelocity[1];
    }
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