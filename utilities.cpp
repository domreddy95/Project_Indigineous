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
//private:
//    int num_dim1;
//    int num_dim2;
//    double h_dim1;
//    double h_dim2;
//    
public:
    double charge;
    double nodal_volume;
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
    std::vector< std::vector<mesh> > temp_array2;
    temp_array2.reserve(num1+2);

    for (int i=0;i<num1+2;i++)
    {
        std::vector<mesh> temp_array;
        temp_array.reserve(num2+2);
        for (int j=0;j<num2+2;j++)
        {
            temp_array.emplace_back(mesh());
        }
        temp_array2.emplace_back(temp_array);
    }
    
    return temp_array2;
}

void calc_nodal_volume(std::vector< std::vector<mesh> > &mesh_array, double min_dim2, double h_dim1, double h_dim2, int num1, int num2)
{
    for (int i=1;i<num1+1;i++)
    {
        for (int j=1;j<num2+1;j++)
        {
            double j_min = (j-1)-0.5;
            double j_max = (j-1)+0.5;
            if (j_min<0){j_min=0;}
            if (j_max>num2){j_max=num2;}
            double factor;
            if (i == 1 || i == num1){factor = 0.5;}
            else{factor = 1;}
            double r_min = min_dim2 + (h_dim2*j_min);
            double r_max = min_dim2 + (h_dim2*j_max);
            mesh_array[i][j].nodal_volume = factor*h_dim1*((r_max*r_max)-(r_min*r_min));
        }
    }
    
}

void charge_weighting(std::vector< std::vector<mesh>> &mesh_arr, std::vector<species> &species_arr,
double min_dim1, double min_dim2, double h_dim1, double h_dim2, int num1, int num2)
{
    for (int i=1;i<num1+1;i++)
    {    
        for (int j=1;j<num2+1;j++)
        {
            mesh_arr[i][j].charge = 0;
        }
        
    }
    
    int n = species_arr.size();
    double particle_charge = species_arr[0].GetCharge();
    for (int i=5000; i<5001; i++)
    {
        double px = species_arr[i].position[0] - min_dim1;
        double py = species_arr[i].position[1] - min_dim2;
        int ix0 = px/h_dim1;
        int iy0 = py/h_dim2;
        int ix1 = ix0 + 1;
        int iy1 = iy0 + 1;
        double hx = ((px/h_dim1)-ix0);
        double hy = ((py/h_dim2)-iy0);
        mesh_arr[ix0][iy0].charge += (1-hx)*(1-hy)*particle_charge;
        mesh_arr[ix1][iy0].charge += hx*(1-hy)*particle_charge;
        mesh_arr[ix0][iy1].charge += (1-hx)*hy*particle_charge;
        mesh_arr[ix1][iy1].charge += hx*hy*particle_charge;
    }
}


void possion_solve(std::vector< std::vector<mesh>> &mesh_arr,double h_dim1, double h_dim2, double epsilon, double w)
{
    int num1 = mesh_arr.size() - 2;
    int num2 = mesh_arr[0].size() - 2;

    for (int j=1;j<num2+1;j++)
    {
        mesh_arr[1][j].potential = 200;       // SET LEFT BOUNDARY CONDITIONS
    }
    
    for (int j=1;j<num2+1;j++)
    {
        mesh_arr[num1][j].potential = -200;    // SET RIGHT BOUNDARY CONDITIONS
    }

    for (int i=0;i<num1+2;i++)
    {
        mesh_arr[i][1].potential = mesh_arr[i][2].potential;
    }

    for (int i=0;i<num1+2;i++)
    {
        mesh_arr[i][num2].potential = mesh_arr[i][num2-1].potential;
    }
    
    for (int i=2;i<num1;i++)
    {    
        for (int j=2;j<num2;j++)
        {
            if ((i+j)%2 == 0)
            {
                double prev = mesh_arr[i][j].potential;
                double t1 = mesh_arr[i][j].charge/mesh_arr[i][j].nodal_volume/epsilon;
                //std::cout << "v = " << t1 << std::endl;
                double t2 = (mesh_arr[i][j+1].potential + mesh_arr[i][j-1].potential)/(h_dim2*h_dim2);
                double t3 = (mesh_arr[i][j+1].potential - mesh_arr[i][j-1].potential)/(2.0*h_dim2*h_dim2*(j-1));
                double t4 = (mesh_arr[i+1][j].potential + mesh_arr[i-1][j].potential)/(h_dim1*h_dim1);
                double t5 = (2.0/(h_dim2*h_dim2)) + (2.0/(h_dim1*h_dim1));
                mesh_arr[i][j].potential = prev + w*(((t1+t2+t4+t3)/t5)-prev);
                
            }
        }
        
    }
    
    for (int i=2;i<num1;i++)
    {    
        for (int j=2;j<num2;j++)
        {
            if ((i+j)%2 != 0)
            {
                double prev = mesh_arr[i][j].potential;
                double t1 = mesh_arr[i][j].charge/mesh_arr[i][j].nodal_volume/epsilon;
                double t2 = (mesh_arr[i][j+1].potential + mesh_arr[i][j-1].potential)/(h_dim2*h_dim2);
                double t3 = (mesh_arr[i][j+1].potential - mesh_arr[i][j-1].potential)/(2.0*h_dim2*h_dim2*(j-1));
                double t4 = (mesh_arr[i+1][j].potential + mesh_arr[i-1][j].potential)/(h_dim1*h_dim1);
                double t5 = (2.0/(h_dim2*h_dim2)) + (2.0/(h_dim1*h_dim1));
                mesh_arr[i][j].potential = prev + w*(((t1+t2+t3+t4)/t5)-prev);
            }
        }
        
    }
    
    for (int i=1;i<num1+1;i++)
    {
        mesh_arr[i][0].potential = 2.0*mesh_arr[i][1].potential - mesh_arr[i][2].potential;
    }

    for (int i=1;i<num1+1;i++)
    {
        mesh_arr[i][num2+1].potential = (2.0*mesh_arr[i][num2].potential) - (mesh_arr[i][num2-1].potential);
    }

    for (int j=0;j<num2+2;j++)
    {
        mesh_arr[num1+1][j].potential = 2.0*mesh_arr[num1][j].potential - mesh_arr[num1-1][j].potential;
    }

    for (int j=0;j<num2+2;j++)
    {
        mesh_arr[0][j].potential = 2.0*mesh_arr[1][j].potential - mesh_arr[2][j].potential;
    }
    
}


//int main()
//{
//    std::vector<std::vector<mesh>> grid = CreateMesh(0,2,0,1,21,11);
//    std::vector <species> electrons = load_species(1,3,1,2,10000);
//    double x = electrons[5000].position[0];
//    double y = electrons[5000].position[1];
//    std::cout << "x = " << x << std::endl;
//    std::cout << "y = " << y << std::endl;
//    charge_weighting(grid,electrons,1,1,0.1,0.1,21,11);
//    calc_nodal_volume(grid,0,0.1,0.1,21,11);
//    for (int i=0;i<50001;i++)
//    {
//        possion_solve(grid,0.1,0.1,1,1.5);
//    }
//    int num1 = 21;
//    int num2 = 11;
//    
//    for (int i=0;i<13;i++)
//    {
//    std::cout << "v = " << grid[0][i].potential << std::endl;
//    }
//    return 0;
//}