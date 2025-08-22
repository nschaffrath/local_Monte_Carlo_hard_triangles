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
    int steps_per_sweep     = 10000000;
    int number_of_sweeps    = 750;


    cout << "Start of simulation" << endl;
    cout << "Number of particles: " << number_of_particles << endl;

    double boxlength = pow(((number_of_particles*4*M_PI*pow(radius,3))/packing_fraction),1./3.);
    cout << boxlength << endl;
    double array_box[3] ={boxlength, sqrt(3./4.) * boxlength, boxlength};
    bool pbc[3] = {true, true, true};

    Simulationvolume simvol(3, array_box, pbc);
    MC_Manager manager(0.20, simvol);

    int n = 30;
    double rest_length = boxlength / n;
    manager.triangulate_sphere_using_triangles(4, 0.25, 4.0, 100.0);

    vec middle_of_simulation_volume(boxlength/2.0, boxlength/2.0, boxlength/2.0);
    //manager.add_n_particles_with_random_positions(number_of_particles, 0.5, -1);

    for(int i = 0; i < 30; i++)
    {
        manager.add_single_particle_with_given_position(middle_of_simulation_volume, 4.5, -1);
    }

    vec point_in_simulation_volume(0.75 * boxlength, boxlength/2.0, boxlength/5.0);
    //manager.add_n_particles_with_random_positions(number_of_particles, 0.5, -1);

    for(int i = 0; i < 2500; i++)
    {
        manager.add_single_particle_with_given_position(point_in_simulation_volume, 3., -1);
    }


    manager.triangulate_periodic_membrane(n, n, 0.3 * boxlength, 0.5, rest_length, 5);
    manager.triangulate_periodic_membrane(n, n, 0.7 * boxlength, 0.5, rest_length, 5);

    manager.setup_neighbor_list_3d();

    int previous_number_of_overlaps = manager.check_validity();
    int current_number_of_overlaps = 0;

    ofstream outputfile_position;
	outputfile_position.open("test2.txt");
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
