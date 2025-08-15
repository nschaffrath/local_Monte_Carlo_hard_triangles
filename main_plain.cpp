#include <iostream>
#include "include/vec.h"
#include "include/particle.h"
#include "include/simulation_volume.h"
#include "include/mc_manager.h"
using namespace std;


int main()
{
    cout << "Setting up the simulation!" << endl;

    int number_of_particles = 2000;
    double radius           = 0.5;
    double packing_fraction = 0.01;
    int iterations          = 1;
    int steps_per_sweep     = 100000;//int(3.5*number_of_particles);
    int number_of_sweeps    = 100;


    cout << "Start of simulation" << endl;
    cout << "Number of particles: " << number_of_particles << endl;

    double boxlength = pow(((number_of_particles*4*M_PI*pow(radius,3))/packing_fraction),1./3.);
    cout << boxlength << endl;
    double array_box[3] ={boxlength, boxlength, boxlength};
    bool pbc[3] = {true, true, true};

    Simulationvolume simvol(3, array_box, pbc);
    MC_Manager manager(1.00, simvol);

//    manager.add_n_particles_with_random_positions(number_of_particles, radius, -1);

    manager.add_single_triangle_with_random_position(1, 10.01);
//    manager.add_single_triangle_with_random_position(1, 0.01);
//    manager.add_single_triangle_with_random_position(1, 0.01);
//    manager.add_single_triangle_with_random_position(1, 0.01);
/*    manager.add_single_triangle_with_random_position(1, 0.01);
    manager.add_single_triangle_with_random_position(1, 0.01);
    manager.add_single_triangle_with_random_position(1, 0.01);
*/
    manager.setup_neighbor_list_3d();

    ofstream outputfile_position;
	outputfile_position.open("test_n_1000.txt");
    manager.print_configuration_to_file(outputfile_position);

//    double new_spring_constant = 0.03;

    cout << "Starting simulation..." << endl;
    for(int k = 0; k < number_of_sweeps; k++)
    {
        for(int i = 0; i < steps_per_sweep; i++)
        {
            manager.mc_step();
        }

   /*     manager.change_spring_constants_of_all_triangles(new_spring_constant);
        if(new_spring_constant < 10)
        {
            cout << "Spring constant: " << new_spring_constant << endl;
            new_spring_constant *= 2;
        }*/
        
        if(k % 1 == 0)
        {
            cout << "Step: " << k << " of " << number_of_sweeps << endl;
            cout << manager.check_validity() << " overlap(s) detected!" << endl;
        }
        manager.print_configuration_to_file(outputfile_position);
    }
    cout << "Simulation done!" << endl;
    manager.remove_all_particles();

    outputfile_position.close();

    return 0;
}
