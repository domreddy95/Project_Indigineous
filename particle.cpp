#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>

class particle
{
private:
    long double charge;
    long double mass;
    long double position[3];
    long double velocity[3];
    std::string species_name;
    int species_identifier;
public:
    particle();
    particle(long double arg1, long double mass, long double position[3], long double velocity[3], std::string species_name, int species_identifier);
    ~particle();
    void particle::print_properties(particle);
    std::vector <particle> particle::load_particles(double min_x,double max_x,double min_y,double max_y,double dx,double dy,double density);
};

particle::particle()
{
    charge = 0;
    mass = 1;
    position[0] = 0;    position[1] = 0;    position[2] = 0;
    velocity[0] = 0;    velocity[1] = 0;    velocity[2] = 0;
    species_name = "unnamed";
    species_identifier = 0;
}

particle::particle(long double arg1, long double arg2,long double arg3[3], long double arg4[3], std::string arg5, int arg6)
{
    charge = arg1;
    mass = arg2;
    position[0] = arg3[0];    position[1] = arg3[1];    position[2] = arg3[2];
    velocity[0] = arg4[0];    velocity[1] = arg4[1];    velocity[2] = arg4[2];
    species_name = arg5;
    species_identifier = arg6;
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

std::vector <particle> particle::load_particles(double min_x,double max_x,double min_y,double max_y,double dx,double dy,double density)
{
    std::vector <particle> V;

}
