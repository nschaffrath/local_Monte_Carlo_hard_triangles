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
    

        void change_mc_length(double given_mc_length)
        {
            mclength = given_mc_length;
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


        void add_single_particle_with_given_position(vec given_position, double given_radius, int given_triangle_index)
        {
                Particle new_particle(given_position, given_position, given_radius, particles.size(), given_triangle_index);
                particles.push_back(new_particle);
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


        void add_single_triangle_with_given_positions(int given_index_first_edge, int given_index_second_edge, int given_index_third_edge, double given_restlength, double given_spring_constant)
        {
            Triangle new_triangle(given_index_first_edge, given_index_second_edge, given_index_third_edge, given_restlength, given_spring_constant, number_of_triangles);
            number_of_triangles++;
            triangles.push_back(new_triangle);
        }


        void divide_triangle_into_smaller_triangles_and_add_new_particles(vec given_position_1, vec given_position_2, vec given_position_3, double given_radius, int given_layer)
        {
            if(given_layer > 0)
            {
                vec position_1;
                position_1.x = 0.5 * (given_position_1 + given_position_2).x;
                position_1.y = 0.5 * (given_position_1 + given_position_2).y;
                position_1.z = 0.5 * (given_position_1 + given_position_2).z;

                vec position_2; 
                position_2.x = 0.5 * (given_position_1 + given_position_3).x;
                position_2.y = 0.5 * (given_position_1 + given_position_3).y;
                position_2.z = 0.5 * (given_position_1 + given_position_3).z;
            
                vec position_3;
                position_3.x = 0.5 * (given_position_2 + given_position_3).x;
                position_3.y = 0.5 * (given_position_2 + given_position_3).y;
                position_3.z = 0.5 * (given_position_2 + given_position_3).z;
            
                vec middle_point_1 = container->calculate_size_of_vector_boundary_conform(position_1);
                vec middle_point_2 = container->calculate_size_of_vector_boundary_conform(position_2);
                vec middle_point_3 = container->calculate_size_of_vector_boundary_conform(position_3);  
            
                bool add_p1 = true;
                bool add_p2 = true;
                bool add_p3 = true;
                                
                for(auto element : particles)
                {
                    if(container->calculate_periodic_difference_vector(middle_point_1, element.position).norm() == 0)
                    {
                        add_p1 = false;
                    }
                    if(container->calculate_periodic_difference_vector(middle_point_2, element.position).norm() == 0)
                    {
                        add_p2 = false;
                    }
                    if(container->calculate_periodic_difference_vector(middle_point_3, element.position).norm() == 0)
                    {
                        add_p3 = false;
                    }
                }

                if(add_p1 == true)
                {
                    this->add_single_particle_with_given_position(middle_point_1, given_radius, -1);
                }
                if(add_p2 == true)
                {
                    this->add_single_particle_with_given_position(middle_point_2, given_radius, -1);
                }
                if(add_p3 == true)
                {
                    this->add_single_particle_with_given_position(middle_point_3, given_radius, -1);
                }

                this->divide_triangle_into_smaller_triangles_and_add_new_particles(given_position_1, middle_point_1, middle_point_2, given_radius, given_layer - 1);
                this->divide_triangle_into_smaller_triangles_and_add_new_particles(given_position_2, middle_point_1, middle_point_3, given_radius, given_layer - 1);
                this->divide_triangle_into_smaller_triangles_and_add_new_particles(given_position_3, middle_point_2, middle_point_3, given_radius, given_layer - 1);
                this->divide_triangle_into_smaller_triangles_and_add_new_particles(middle_point_1,   middle_point_2, middle_point_3, given_radius, given_layer - 1);
            }
        }


        void triangulate_sphere_using_triangles(int triangulation_depth, double given_radius, double given_restlength, double given_spring_constant)
        {
            vector< vector< vec > > edge_points_initial_octahedron;
            edge_points_initial_octahedron.resize(8);
            for(int i = 0; i < edge_points_initial_octahedron.size(); i++)
            {
                edge_points_initial_octahedron[i].resize(3);
            }

            vec size_of_simulation_volume = container->give_size_of_simulation_volume();
            double boxlength_x = size_of_simulation_volume.x;
            double boxlength_y = size_of_simulation_volume.y;
            double boxlength_z = size_of_simulation_volume.z;

            vector< vec > edge_positions;

            vec edge_position_1(0.35 * boxlength_x, 0.50 * boxlength_y, 0.50 * boxlength_z);
            vec edge_position_2(0.65 * boxlength_x, 0.50 * boxlength_y, 0.50 * boxlength_z);
            vec edge_position_3(0.50 * boxlength_x, 0.35 * boxlength_y, 0.50 * boxlength_z);
            vec edge_position_4(0.50 * boxlength_x, 0.65 * boxlength_y, 0.50 * boxlength_z);    
            vec edge_position_5(0.50 * boxlength_x, 0.50 * boxlength_y, 0.35 * boxlength_z);
            vec edge_position_6(0.50 * boxlength_x, 0.50 * boxlength_y, 0.65 * boxlength_z);   

            double diameter_octahedron = 0.3 * boxlength_x;

            edge_positions.push_back(edge_position_1);
            edge_positions.push_back(edge_position_2);
            edge_positions.push_back(edge_position_3);
            edge_positions.push_back(edge_position_4);
            edge_positions.push_back(edge_position_5);
            edge_positions.push_back(edge_position_6);

            edge_points_initial_octahedron[0][0] = edge_position_1;
            edge_points_initial_octahedron[0][1] = edge_position_3;
            edge_points_initial_octahedron[0][2] = edge_position_5;

            edge_points_initial_octahedron[1][0] = edge_position_1;
            edge_points_initial_octahedron[1][1] = edge_position_3;
            edge_points_initial_octahedron[1][2] = edge_position_6;


            edge_points_initial_octahedron[2][0] = edge_position_1;
            edge_points_initial_octahedron[2][1] = edge_position_4;
            edge_points_initial_octahedron[2][2] = edge_position_5;

            edge_points_initial_octahedron[3][0] = edge_position_1;
            edge_points_initial_octahedron[3][1] = edge_position_4;
            edge_points_initial_octahedron[3][2] = edge_position_6;


            edge_points_initial_octahedron[4][0] = edge_position_2;
            edge_points_initial_octahedron[4][1] = edge_position_4;
            edge_points_initial_octahedron[4][2] = edge_position_5;

            edge_points_initial_octahedron[5][0] = edge_position_2;
            edge_points_initial_octahedron[5][1] = edge_position_4;
            edge_points_initial_octahedron[5][2] = edge_position_6;


            edge_points_initial_octahedron[6][0] = edge_position_2;
            edge_points_initial_octahedron[6][1] = edge_position_3;
            edge_points_initial_octahedron[6][2] = edge_position_5;

            edge_points_initial_octahedron[7][0] = edge_position_2;
            edge_points_initial_octahedron[7][1] = edge_position_3;
            edge_points_initial_octahedron[7][2] = edge_position_6;
    
            double initial_edge_length = (container->calculate_periodic_difference_vector(edge_position_1, edge_position_3)).norm();
            cout << initial_edge_length << endl;
            initial_edge_length /= pow(2,triangulation_depth);
            cout << initial_edge_length << endl;

            for(int i = 0; i < edge_positions.size(); i++)
            {
                this->add_single_particle_with_given_position(edge_positions[i], given_radius, -1);
            }

            for(int i = 0; i < edge_points_initial_octahedron.size(); i++)
            {
                this->divide_triangle_into_smaller_triangles_and_add_new_particles(edge_points_initial_octahedron[i][0], edge_points_initial_octahedron[i][1], edge_points_initial_octahedron[i][2], given_radius, triangulation_depth);
            }

            cout << "Number of particles: " << particles.size() << endl;
            for(int i = 0; i < particles.size(); i++)
            {
                for(int j = 0; j < particles.size(); j++)
                {
                    if(&particles[i] != &particles[j])
                    {
                        double distance_1_2 = (container->calculate_periodic_difference_vector(particles[i].position, particles[j].position)).norm();
                        if(distance_1_2 < 1.01 * initial_edge_length)
                        {
                            for(int k = 0; k < particles.size(); k++)
                            {
                                if(&particles[i] != &particles[k] && &particles[j] != &particles[k])
                                {
                                    double distance_1_3 = (container->calculate_periodic_difference_vector(particles[i].position, particles[k].position)).norm();
                                    if(distance_1_3 < 1.01 * initial_edge_length)
                                    {
                                        double distance_2_3 = (container->calculate_periodic_difference_vector(particles[j].position, particles[k].position)).norm();
                                        if(distance_2_3 < 1.01 * initial_edge_length)
                                        {
                                            bool triangle_already_registered = false;

                                            for(int l = 0; l < particles[i].triangle_index.size(); l++)
                                            {
                                                for(int m = 0; m < particles[j].triangle_index.size(); m++)
                                                {
                                                    if(particles[i].triangle_index[l] == particles[j].triangle_index[m])
                                                    {
                                                        for(int n = 0; n < particles[k].triangle_index.size(); n++)
                                                        {
                                                            if(particles[i].triangle_index[l] == particles[k].triangle_index[n])
                                                            {
                                                                if(particles[i].triangle_index[l] != -1)
                                                                {
                                                                    triangle_already_registered = true;
                                                                }
                                                            }                                                                        
                                                        }
                                                    }
                                                }
                                            }                                                
                
                                            if(triangle_already_registered == false)
                                            {
                                                this->add_single_triangle_with_given_positions(particles[i].particle_index, particles[j].particle_index, particles[k].particle_index, given_restlength, given_spring_constant);   
                                                particles[i].triangle_index.push_back(number_of_triangles - 1);
                                                particles[j].triangle_index.push_back(number_of_triangles - 1);
                                                particles[k].triangle_index.push_back(number_of_triangles - 1);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            for(int i = 0; i < particles.size(); i++)
            {
                vec middle_of_simulation_volume(boxlength_x/2.0 , boxlength_y/2.0, boxlength_z/2.0);
                vec radial_vector = container->calculate_periodic_difference_vector(middle_of_simulation_volume, particles[i].position);
                radial_vector.normalize();

                radial_vector.x *= diameter_octahedron / 2.0;
                radial_vector.y *= diameter_octahedron / 2.0;
                radial_vector.z *= diameter_octahedron / 2.0;
                particles[i].move_to(middle_of_simulation_volume + radial_vector);
                
            }

            for(auto triangle : triangles)
            {
                particles[triangle.indices_of_particles_in_triangle[0]].triangle_index.push_back(triangle.triangle_index);
                particles[triangle.indices_of_particles_in_triangle[1]].triangle_index.push_back(triangle.triangle_index);
                particles[triangle.indices_of_particles_in_triangle[2]].triangle_index.push_back(triangle.triangle_index);
            }

            for(int i = 0; i < triangles.size(); i++)
            {
                double new_rest_length = (container->calculate_periodic_difference_vector(particles[triangles[i].indices_of_particles_in_triangle[0]].position, particles[triangles[i].indices_of_particles_in_triangle[1]].position)).norm();
                triangles[i].change_rest_length(new_rest_length);
            }
            cout << "A sphere, which consists out of " << triangles.size() << " triangles was successfully constructed!" << endl;
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
                outputfile << element.give_position().x << " " << element.give_position().y << " " << element.give_position().z << " " << element.radius << " " << 0 << " " << 0 << " " << 1 << endl;
            }
        }

        void refresh_cell_list()
        {
            neighborlist.refresh_cell_list(particles);
        }

        void setup_neighbor_list_3d()
        {
            double perfect_cell_size = 5.5;

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


        void mc_warmup_step()
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
            int number_of_triangles_connected_with_particle = current_particle->triangle_index.size();
            for(int i = 0; i < number_of_triangles_connected_with_particle; i++)
            {
                int current_index = current_particle->triangle_index[i];
                if( current_index >= 0)
                {
                    for(auto element: triangles[current_index].indices_of_particles_in_triangle)
                    {
                        if(current_particle->particle_index != element)
                        {
                            double current_distance_between_particles   = (container->calculate_periodic_difference_vector(current_particle->position, particles[element].position)).norm();
                            if(current_distance_between_particles < 1.25 * triangles[current_index].rest_length)
                            {
                                double proposed_distance_between_particles  = (container->calculate_periodic_difference_vector(proposed_position, particles[element].position)).norm();
                                if(proposed_distance_between_particles > 1.25 * triangles[current_index].rest_length)
                                {
                                    particle_can_be_moved = false;
                                }
                            }
                        }
                    }
                }
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
            if(particle_can_be_moved == true)
            {
                int number_of_triangles_connected_with_particle = current_particle->triangle_index.size();
                for(int i = 0; i < number_of_triangles_connected_with_particle; i++)
                {
                    int current_index = current_particle->triangle_index[i];
                    if( current_index >= 0)
                    {
                        for(auto element: triangles[current_index].indices_of_particles_in_triangle)
                        {
                            if(current_particle->particle_index != element)
                            {
                                double current_distance_between_particles = (container->calculate_periodic_difference_vector(current_particle->position, particles[element].position)).norm();
                                double current_energy = 0.5 * triangles[current_index].spring_constant * pow(current_distance_between_particles - triangles[current_index].rest_length, 2);
                                double proposed_distance_between_particles = (container->calculate_periodic_difference_vector(proposed_position, particles[element].position)).norm();
                                double energy_of_proposed_position = 0.5 * triangles[current_index].spring_constant * pow(proposed_distance_between_particles - triangles[current_index].rest_length, 2);
                                delta_E += (energy_of_proposed_position - current_energy);
                            }
                        }
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
