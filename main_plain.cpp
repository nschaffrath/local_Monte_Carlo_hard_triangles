#include <iostream>
#include "include/vec.h"
#include "include/particle.h"
#include "include/simulation_volume.h"
#include "include/mc_manager.h"
using namespace std;


int main()
{
    cout << "Setting up the simulation!" << endl;

    int number_of_particles = 4000;
    double radius           = 0.5;
    double packing_fraction = 0.002;
    int iterations          = 1;
    int steps_per_sweep     = 1000000;
    int number_of_sweeps    = 750;


    cout << "Start of simulation" << endl;
    cout << "Number of particles: " << number_of_particles << endl;

    double boxlength = pow(((number_of_particles*4*M_PI*pow(radius,3))/packing_fraction),1./3.);
    cout << boxlength << endl;
    double array_box[3] ={boxlength, boxlength, boxlength};
    bool pbc[3] = {true, true, true};

    Simulationvolume simvol(3, array_box, pbc);
    MC_Manager manager(0.20, simvol);

    manager.triangulate_sphere_using_triangles(5, 0.05, 4.0, 10.0);

    vec middle_of_simulation_volume(boxlength/2.0, boxlength/2.0, boxlength/2.0);
    //manager.add_n_particles_with_random_positions(number_of_particles, 0.5, -1);

    for(int i = 0; i < 300; i++)
    {
        manager.add_single_particle_with_given_position(middle_of_simulation_volume, 2.5, -1);
    }
    manager.setup_neighbor_list_3d();
//    manager.print_entire_cell_list();

    int previous_number_of_overlaps = manager.check_validity();
    int current_number_of_overlaps = 0;

    ofstream outputfile_position;
	outputfile_position.open("sphere_radius_0.05_spring_constant_10.0_and_n_300_particles_with_radius_2.5_within.txt");
    manager.print_configuration_to_file(outputfile_position);

    cout << "Starting simulation..." << endl;
    for(int k = 0; k < number_of_sweeps; k++)
    {
        for(int i = 0; i < steps_per_sweep; i++)
        {
            manager.mc_step();
        }        

        if(k % 1 == 0)
        {
            cout << "Step: " << k << " of " << number_of_sweeps << endl;
            current_number_of_overlaps = manager.check_validity();
            cout << current_number_of_overlaps << " overlap(s) detected!" << endl;

            if(current_number_of_overlaps > previous_number_of_overlaps)
            {
                cout << "Error in neighborlist!" << endl;
            }
            previous_number_of_overlaps = current_number_of_overlaps;
        }
        manager.print_configuration_to_file(outputfile_position);
    }

    cout << "Adjust spring constants within a single triangle!" << endl;

    cout << "Simulation done!" << endl;
    manager.remove_all_particles();

    outputfile_position.close();

    return 0;
}
