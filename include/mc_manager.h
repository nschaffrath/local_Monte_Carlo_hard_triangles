#pragma once
#include <vector>
#include <cmath>
#include <random>
#include <fstream>

#include "particle.h"
#include "triangle.h"
#include "neighborlist_3d.h"
using namespace std;

// Create a random device and seed the generator
random_device rd;  // Obtain a random number from hardware
mt19937 gen(rd()); // Seed the generator

class MC_Manager
{
    private:
        double kBT = 1.00;
        // Maximal displacement length per EC move
        double mclength;
        vec mc_direction;

        // Vector that contains all particles within the simulation volume
        vector< Particle > particles;
        
        // Vector that contains all triangles within the simulation volume
        vector< Triangle > triangles;
        int number_of_triangles = 0;
        
        // Simulation volume, in which the particles are moved
        Simulationvolume* container;

        // Neighborlist 3d
        Neighborlist_3d neighborlist;

        // Current particle
        Particle* current_particle;

    public:
        // Constructor of the class        
        MC_Manager(double given_mclength, Simulationvolume & given_container)
            : mclength(given_mclength)
            {
                container = &given_container;
            };

        // Destructor of the class
        ~MC_Manager() {};

        // Add a single particle
        void add_particle(const Particle new_particle)
        {
            particles.push_back(new_particle);
        }
    
        void add_n_particles_with_random_positions(int given_number_of_particles_to_add, double given_radius, int given_triangle_index)
        {
            uniform_real_distribution<double> distribution_x(0.0, container->give_size_of_simulation_volume().x);
            uniform_real_distribution<double> distribution_y(0.0, container->give_size_of_simulation_volume().y);
            uniform_real_distribution<double> distribution_z(0.0, container->give_size_of_simulation_volume().z);


            for(int i = 0; i < given_number_of_particles_to_add; i++)
            {
                vec random_position_vector(distribution_x(gen), distribution_y(gen), distribution_z(gen));
                Particle new_particle(random_position_vector, random_position_vector, given_radius, particles.size(), given_triangle_index);
                particles.push_back(new_particle);
            }
        }

        void add_single_triangle_with_random_position(double given_restlength, double given_spring_constant)
        {
            uniform_real_distribution<double> distribution_x(0.0, container->give_size_of_simulation_volume().x);
            uniform_real_distribution<double> distribution_y(0.0, container->give_size_of_simulation_volume().y);
            uniform_real_distribution<double> distribution_z(0.0, container->give_size_of_simulation_volume().z);

            vec random_position_vector1(distribution_x(gen), distribution_y(gen), distribution_z(gen));
            Particle particle_1(random_position_vector1, random_position_vector1, 0.0, particles.size(), number_of_triangles);
            particles.push_back(particle_1);
            vec random_position_vector2(distribution_x(gen), distribution_y(gen), distribution_z(gen));
            Particle particle_2(random_position_vector2, random_position_vector2, 0.0, particles.size(), number_of_triangles);
            particles.push_back(particle_2);
            vec random_position_vector3(distribution_x(gen), distribution_y(gen), distribution_z(gen));
            Particle particle_3(random_position_vector3, random_position_vector3, 0.0, particles.size(), number_of_triangles);           
            particles.push_back(particle_3);

            Triangle new_triangle(particles.size() - 3, particles.size() - 2, particles.size() - 1, given_restlength, given_spring_constant, number_of_triangles);
            number_of_triangles++;
            triangles.push_back(new_triangle);
        }

        vec create_normalized_3d_random_vector()
        {
            uniform_real_distribution<double> distribution(-1.0, 1.0);

            vec random_vector(distribution(gen),distribution(gen), distribution(gen));
            random_vector.normalize();
            return random_vector;
        }

        // Remove all particles
        void remove_all_particles()
        {
            cout << "Removing all particles..." << endl;
            particles.clear();
            cout << "Done!" << endl;
        }

        // Prints the number of registered particles
        void print_number_of_registered_particles() const
        {
            cout << "Current number of registered particles: " << particles.size() << endl;
        }

        // Prints the number of registered particles
        int give_number_of_registered_particles() const
        {
            return particles.size();
        }

        void print_configuration_to_file(std::ofstream &outputfile)
        {
            outputfile << particles.size() << endl;
            // Specify periodic boundary conditions
            outputfile << "pbc=\"" << 1 << " " << 1 << " " << 0 << "\" ";
            // Specify simulation cell geometry
            outputfile << "Lattice=\"" << container->give_size_of_simulation_volume().x << " 0.0 0.0 0.0 " << container->give_size_of_simulation_volume().y << " 0.0 0.0 0.0 " << 0 << "\" ";
            // Specify input types
            outputfile << "Properties=pos:R:3:radius:R:1:color:R:3 " << endl;

            for (auto element : particles)
            {
                outputfile << element.give_position().x << " " << element.give_position().y << " " << element.give_position().z << " " << 0.5 << " " << 0 << " " << 0 << " " << 1 << endl;
            }
        }

        void setup_neighbor_list_3d()
        {
            double perfect_cell_size = 3.5;

            int intended_number_of_cells_x = int(container->give_size_of_simulation_volume().x / perfect_cell_size);
            int intended_number_of_cells_y = int(container->give_size_of_simulation_volume().y / perfect_cell_size);
            int intended_number_of_cells_z = int(container->give_size_of_simulation_volume().z / perfect_cell_size);
            
            double cellsize_x = container->give_size_of_simulation_volume().x / (1. * intended_number_of_cells_x);
            double cellsize_y = container->give_size_of_simulation_volume().y / (1. * intended_number_of_cells_y);
            double cellsize_z = container->give_size_of_simulation_volume().z / (1. * intended_number_of_cells_z);

            cout << "Cell size: " << cellsize_x << endl;

            neighborlist.setup_cell_size_and_number_of_cells(cellsize_x, cellsize_y, cellsize_z, intended_number_of_cells_x, intended_number_of_cells_y, intended_number_of_cells_z);
            neighborlist.initialize_cell_list(particles);
        }

        vector< Particle* > create_list_with_neighboring_particles(Particle* current_particle)
        {
            int bin_x = int( current_particle->give_old_position().x / neighborlist.size_cell_element_x);
            int bin_y = int( current_particle->give_old_position().y / neighborlist.size_cell_element_y);
            int bin_z = int( current_particle->give_old_position().z / neighborlist.size_cell_element_z);

            vector< Particle* > neighboring_particles;

            // Particle at the boundary (assuming periodic boundary conditions)
            int array_bins_x[3] = {-1, -1, -1};
            if(bin_x == 0)
            {
                array_bins_x[0] = neighborlist.number_of_cells_x - 1;
                array_bins_x[1] = 0;
                array_bins_x[2] = 1;
            }
            if(bin_x == neighborlist.number_of_cells_x - 1)
            {
                array_bins_x[0] = neighborlist.number_of_cells_x - 2;
                array_bins_x[1] = neighborlist.number_of_cells_x - 1;
                array_bins_x[2] = 0;
            }
            if(bin_x > 0 && bin_x < neighborlist.number_of_cells_x - 1)
            {
                array_bins_x[0] = bin_x - 1;
                array_bins_x[1] = bin_x;
                array_bins_x[2] = bin_x + 1;
            }

            int array_bins_y[3];
            if(bin_y == 0)
            {
                array_bins_y[0] = neighborlist.number_of_cells_y - 1;
                array_bins_y[1] = 0;
                array_bins_y[2] = 1;
            }
            if(bin_y == neighborlist.number_of_cells_y - 1)
            {
                array_bins_y[0] = neighborlist.number_of_cells_y - 2;
                array_bins_y[1] = neighborlist.number_of_cells_y - 1;
                array_bins_y[2] = 0;
            }
            if(bin_y > 0 && bin_y < neighborlist.number_of_cells_y - 1)
            {
                array_bins_y[0] = bin_y - 1;
                array_bins_y[1] = bin_y;
                array_bins_y[2] = bin_y + 1;
            }

            int array_bins_z[3];
            if(bin_z == 0)
            {
                array_bins_z[0] = neighborlist.number_of_cells_z - 1;
                array_bins_z[1] = 0;
                array_bins_z[2] = 1;
            }
            if(bin_z == neighborlist.number_of_cells_z - 1)
            {
                array_bins_z[0] = neighborlist.number_of_cells_z - 2;
                array_bins_z[1] = neighborlist.number_of_cells_z - 1;
                array_bins_z[2] = 0;
            }
            if(bin_z > 0 && bin_z < neighborlist.number_of_cells_z - 1)
            {
                array_bins_z[0] = bin_z - 1;
                array_bins_z[1] = bin_z;
                array_bins_z[2] = bin_z + 1;
            }

            for(int i = 0; i <= 2; i++)
            {
                for(int j = 0; j <= 2; j++)
                {
                    for(int k = 0; k <= 2; k++)
                    {
                        for(auto element : neighborlist.cell_list[array_bins_x[i]][array_bins_y[j]][array_bins_z[k]])
                        {
                            if(element != current_particle)
                            {
                                neighboring_particles.push_back(element);
                            }
                        }
                    } 
                }
            }
            return neighboring_particles;
        }

        void change_spring_constants_of_all_triangles(double given_spring_constant)
        {
            for(auto element: triangles)
            {
                element.change_spring_constant(given_spring_constant);
                element.print_spring_constant();
            }
        }

        void mc_step()
        {
            uniform_int_distribution<int> distribution_particle(0, give_number_of_registered_particles() - 1);

            current_particle = &particles[distribution_particle(gen)];
            mc_direction = create_normalized_3d_random_vector();

            vector< Particle* > neighboring_particles = create_list_with_neighboring_particles(current_particle);
            vec proposed_position = container->calculate_size_of_vector_boundary_conform(current_particle->position + mc_direction * mclength);

            bool particle_can_be_moved = true;

            // Check for hard-sphere collisions
            for(auto element : neighboring_particles)
            {
                double distance = (container->calculate_periodic_difference_vector(proposed_position, element->position)).norm();
                if(distance < element->radius + current_particle->radius)
                {
                    particle_can_be_moved = false;
                }
            }

            // Check for triangle-interactions
            double delta_E = 0.0;
            if( current_particle->triangle_index >= 0)
            {
                for(auto element: triangles[current_particle->triangle_index].indices_of_particles_in_triangle)
                {
                    if(current_particle->particle_index != element)
                    {
                        double current_distance_between_particles = (container->calculate_periodic_difference_vector(current_particle->position, particles[element].position)).norm();
                        double current_energy = 0.5 * triangles[current_particle->triangle_index].spring_constant * pow(current_distance_between_particles - triangles[current_particle->triangle_index].rest_length, 2);
                        double proposed_distance_between_particles = (container->calculate_periodic_difference_vector(proposed_position, particles[element].position)).norm();
                        double energy_of_proposed_position = 0.5 * triangles[current_particle->triangle_index].spring_constant * pow(proposed_distance_between_particles - triangles[current_particle->triangle_index].rest_length, 2);
                        delta_E += (energy_of_proposed_position - current_energy);
                    }
                }
            }

            uniform_real_distribution<> udist(0, 1);

            if(udist(gen) > exp(- delta_E / kBT))
            {
                particle_can_be_moved = false;
            }


            if(particle_can_be_moved == true)
            {
                // Update particle coordinates
                this->move_particle(current_particle);
                // Update neighbor list
                if(this->neighborlist.remove_particle_from_cell_list(current_particle) == true)
                {
                    this->neighborlist.add_particle_to_cell_list(current_particle);
                }
            } 
            neighboring_particles.clear();
        }

        void move_particle(Particle* current_particle)
        {
            current_particle->old_position = current_particle->position;
            current_particle->position     = current_particle->position + mc_direction * mclength;
            current_particle->position     = container->calculate_size_of_vector_boundary_conform(current_particle->position);
        }

        int check_validity()
        {
            int number_of_overlaps = 0;
            for(auto particle_1 : particles)
            {
                for(auto particle_2 : particles)
                {
                    double periodic_distance = container->calculate_periodic_difference_vector(particle_1.position, particle_2.position).norm();
                    if(periodic_distance > 0)
                    {
                        if(periodic_distance < particle_1.radius + particle_2.radius)
                        {
                            number_of_overlaps++;
                        }
                    }
                }
            }
            return number_of_overlaps/2;
        }
};
