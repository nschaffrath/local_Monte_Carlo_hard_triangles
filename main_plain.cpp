#include <iostream>
#include "include/vec.h"
#include "include/particle.h"
#include "include/simulation_volume.h"
#include "include/mc_manager.h"
using namespace std;


int main()
{
    cout << "Setting up the simulation!" << endl;

    int number_of_particles = 10000;
    double radius           = 0.5;
    double packing_fraction = 0.01;
    int iterations          = 1;
    int steps_per_sweep     = 500 * number_of_particles;
    int number_of_sweeps    = 1000;


    cout << "Start of simulation" << endl;
    cout << "Number of particles: " << number_of_particles << endl;

    double boxlength = pow(((number_of_particles*4*M_PI*pow(radius,3))/packing_fraction),1./3.);
    cout << boxlength << endl;
    double array_box[3] ={boxlength, boxlength, boxlength};
    bool pbc[3] = {true, true, true};

    // Create points to triangulate a sphere
    int triangulation_depth = 2;
    if(triangulation_depth < 0)
    {
        cout << "Error! Triangulation depth must be larger than one!" << endl;
        return 0;
    } 

    Simulationvolume simvol(3, array_box, pbc);
    MC_Manager manager(0.01, simvol);

    manager.triangulate_sphere_using_triangles(4, 0.50, 4.0, 20.0);

    vec middle_of_simulation_volume(boxlength/2.0, boxlength/2.0, boxlength/2.0);
    /*for(int i = 0; i < 6; i++)
    {
        manager.add_single_particle_with_given_position(middle_of_simulation_volume, 2.5, -1);
    }*/
    manager.setup_neighbor_list_3d();

    ofstream outputfile_position;
	outputfile_position.open("test_n_1000.txt");
    manager.print_configuration_to_file(outputfile_position);

//    manager.refresh_cell_list();

    manager.change_mc_length(0.20);
    cout << "Starting simulation..." << endl;
    for(int k = 0; k < 100; k++)
    {
        for(int i = 0; i < 100 * number_of_particles; i++)
        {
            manager.mc_step();
        }
       
        if(k % 1 == 0)
        {
            cout << "Step: " << k << " of " << number_of_sweeps << endl;
            cout << manager.check_validity() << " overlap(s) detected!" << endl;
        }
        manager.print_configuration_to_file(outputfile_position);
    }

    manager.change_mc_length(0.10);
    cout << "Starting simulation..." << endl;
    for(int k = 0; k < 10; k++)
    {
        for(int i = 0; i < 100 * number_of_particles; i++)
        {
            manager.mc_step();
        }
       
        if(k % 1 == 0)
        {
            cout << "Step: " << k << " of " << number_of_sweeps << endl;
            cout << manager.check_validity() << " overlap(s) detected!" << endl;
        }
        manager.print_configuration_to_file(outputfile_position);
    }

    manager.change_mc_length(0.05);
    cout << "Starting simulation..." << endl;
    for(int k = 0; k < 10; k++)
    {
        for(int i = 0; i < 100 * number_of_particles; i++)
        {
            manager.mc_step();
        }
       
        if(k % 1 == 0)
        {
            cout << "Step: " << k << " of " << number_of_sweeps << endl;
            cout << manager.check_validity() << " overlap(s) detected!" << endl;
        }
        manager.print_configuration_to_file(outputfile_position);
    }

    manager.change_mc_length(0.02);
    cout << "Starting simulation..." << endl;
    for(int k = 0; k < 10; k++)
    {
        for(int i = 0; i < 100 * number_of_particles; i++)
        {
            manager.mc_step();
        }
       
        if(k % 1 == 0)
        {
            cout << "Step: " << k << " of " << number_of_sweeps << endl;
            cout << manager.check_validity() << " overlap(s) detected!" << endl;
        }
        manager.print_configuration_to_file(outputfile_position);
    }



    for(int k = 0; k < number_of_sweeps; k++)
    {
        for(int i = 0; i < steps_per_sweep; i++)
        {
            manager.mc_step();
        }
       
        if(k % 1 == 0)
        {
            cout << "Step: " << k << " of " << number_of_sweeps << endl;
            cout << manager.check_validity() << " overlap(s) detected!" << endl;
        }
        manager.print_configuration_to_file(outputfile_position);
    }

    cout << "Adjust spring constants within a single triangle!" << endl;

    cout << "Simulation done!" << endl;
    manager.remove_all_particles();

    outputfile_position.close();

    return 0;
}
