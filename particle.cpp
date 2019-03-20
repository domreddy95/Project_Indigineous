#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>

class particle
{
private:

public:
    std::vector <long double> position;
    std::vector <long double> velocity;
    std::vector <long double> GetPosition(){return position;}
    void SetPosition(std::vector <long double> position){position = position;}
    std::vector <long double> GetVelocity(){return velocity;}
    void SetVelocity(std::vector <long double> velocity){velocity = velocity;}
    particle();
    particle(long double, long double);
    ~particle();  
};

//Define default Constructer
particle::particle()
{
    std::vector <long double> position (2,0.0);
    std::vector <long double> velocity (2,0.0);
    std::cout << " A particle is created" << std::endl;
}

//Parameterized Constructor
particle::particle(long double charge, long double mass)
{
    std::vector <long double> position (2,0.0);
    std::vector <long double> velocity (2,0.0);
    std::cout << " A parametrized particle is created" << std::endl;
}

particle::~particle()
{
    std::cout << " A particle is destroyed";
}



class species: public particle
{
private:
    std::string species_name;
    int species_identifier;
    long double charge;
    long double mass;

public:
    species(): particle(){};
    species(std::string, int, long double, long double);
    ~species();
    long double GetCharge(){return charge;}
    void SetCharge(long double charge){this->charge=charge;}
    long double GetMass(){return mass;}
    void SetMass(long double mass){this->mass = mass;}
    void print_properties();
};

species::species(std::string species_name , int species_identifier, long double charge,
long double mass) : particle()
{
    this->species_name = species_name;
    this->species_identifier = species_identifier;
    this->charge = charge;
    this->mass = mass;
    std::cout << " A parametrized species is created" << std::endl;
}

species::~species()
{
}

void species::print_properties()
{
    using std::cout;
    using std::endl;
    cout << "charge = " << this->charge << endl;
    cout << "mass = " << this->mass << endl;
    cout << "position_x = " << this->position[0]<<endl;
    cout << "velocity_x = " << this->velocity[0]<<endl;
}
int main()
{
    //long double chrg = 1.983742974237498;
    //long double mas = 116.92847928349292;
    species elec();
    particle p;
    p.position[0];
    //ion.print_properties();
    return 0;
}