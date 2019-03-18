#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>

class particle
{
private:
    long double charge;
    long double mass;
    std::string species_name;
    int species_identifier;
public:
    particle(long double arg1, long double mass, std::string species_name, int species_identifier);
    ~particle();
    void particle::print_properties(particle);
};

particle::particle(long double arg1 = 0, long double arg2 = 1, std::string arg3 = "unnamed", int arg4 = 0)
{
    charge = arg1;
    mass = arg2;
    species_name = arg3;
    species_identifier = arg4;
}

particle::~particle()
{
}

void particle::print_properties(particle arg1)
{
    using std::cout;
    using std::endl;
    cout << "charge = " << arg1.charge << endl;
    cout << "mass = " << arg1.mass << endl;
    cout << "species_name = " << arg1.species_name << endl;
    cout << "species_identifier = " << arg1.species_identifier << endl;
}