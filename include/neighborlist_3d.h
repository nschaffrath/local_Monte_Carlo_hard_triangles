#pragma once
#include<cmath>
#include<vector>
#include <algorithm>
#include "particle.h"
using namespace std;

// The class represents a cell neighbor list
class Neighborlist_3d
{
    public:
        double size_cell_element_x, size_cell_element_y, size_cell_element_z; 
        int number_of_cells_x, number_of_cells_y, number_of_cells_z;

        vector< vector< vector< vector< Particle* > > > > cell_list;

        // Constructors of the class
        Neighborlist_3d() {};
        Neighborlist_3d(double given_size_cell_element_x, double given_size_cell_element_y, double given_size_cell_element_z,
            int given_number_of_cells_x, int given_number_of_cells_y, int given_number_of_cells_z)
            : size_cell_element_x(given_size_cell_element_x), size_cell_element_y(given_size_cell_element_y), size_cell_element_z(given_size_cell_element_z),
                number_of_cells_x(given_number_of_cells_x), number_of_cells_y(given_number_of_cells_y), number_of_cells_z(given_number_of_cells_z) {};

        // Destructor of the class
        ~Neighborlist_3d() {};

        // Setup of the cell size
        void setup_cell_size_and_number_of_cells(const double given_size_cell_element_x, const double given_size_cell_element_y, 
            const double given_size_cell_element_z, const int given_number_of_cells_x, const int given_number_of_cells_y, const int given_number_of_cells_z)
        {
            size_cell_element_x = given_size_cell_element_x;
            size_cell_element_y = given_size_cell_element_y;
            size_cell_element_z = given_size_cell_element_z;

            number_of_cells_x = given_number_of_cells_x;
            number_of_cells_y = given_number_of_cells_y;
            number_of_cells_z = given_number_of_cells_z;
        }

        void initialize_cell_list(const vector< Particle > & given_particles)
        {
            cell_list.resize(number_of_cells_x);
            for (int i = 0; i < number_of_cells_x; i++) 
            {
                cell_list[i].resize(number_of_cells_y);
            }

            for (int i = 0; i < number_of_cells_x; i++) 
            {
                for (int j = 0; j < number_of_cells_y; j++)
                {
                    cell_list[i][j].resize(number_of_cells_z);
                }
            }

            cout << "Dimensions of the cell list: " << cell_list.size() << " x " << cell_list[0].size() << " x " << cell_list[0][0].size() << endl;
        
//            for(auto element : given_particles)
//            {
//                this->add_particle_to_cell_list(&element);
//            }
        }

        void refresh_cell_list(const vector< Particle > & given_particles)
        {
            cell_list.clear();
            initialize_cell_list(given_particles);
        }

        bool remove_particle_from_cell_list(Particle* given_particle)
        {
            int old_cell_element_x = int(given_particle->give_old_position().x / size_cell_element_x);
            int old_cell_element_y = int(given_particle->give_old_position().y / size_cell_element_y);
            int old_cell_element_z = int(given_particle->give_old_position().z / size_cell_element_z);

            int new_cell_element_x = int(given_particle->give_position().x / size_cell_element_x);
            int new_cell_element_y = int(given_particle->give_position().y / size_cell_element_y);
            int new_cell_element_z = int(given_particle->give_position().z / size_cell_element_z);

            if(old_cell_element_x == new_cell_element_x && old_cell_element_y == new_cell_element_y && old_cell_element_z == new_cell_element_z)
            {
                //cout << "Don't erase!" << endl;
                return false;
            }
            else
            {
                for(auto element = cell_list[old_cell_element_x][old_cell_element_y][old_cell_element_z].begin(); element != cell_list[old_cell_element_x][old_cell_element_y][old_cell_element_z].end(); ++element)
                {
                    if((*element)->particle_index == given_particle->particle_index)
                    {
                        //cout << "Erase" << endl;
                        cell_list[old_cell_element_x][old_cell_element_y][old_cell_element_z].erase(element);
                        return true;
                    }
                }
                cout << "Strange" << endl;
                return false;
            }
        }

        void add_particle_to_cell_list(Particle* given_particle)
        {
            int cell_element_x = int(given_particle->give_position().x / size_cell_element_x);
            int cell_element_y = int(given_particle->give_position().y / size_cell_element_y);
            int cell_element_z = int(given_particle->give_position().z / size_cell_element_z);
            cell_list[cell_element_x][cell_element_y][cell_element_z].push_back(given_particle);
        }       

        void print_entire_cell_list()
        {
            cout << "Particles in cell list!" << endl;
            for(int i = 0; i < number_of_cells_x; i++)
            {
                for(int j = 0; j < number_of_cells_y; j++)
                {
                    for(int k = 0; k < number_of_cells_z; k++)
                    {
                        cout << "Cell: " << i << " " << j << " " << k << " " << cell_list[i][j][k].size() << endl;
                        auto cell = cell_list[i][j][k];
                        for(auto element = cell.begin(); element != cell.end(); ++element)
                        {
                            cout << (*element)->position.x << " " << (*element)->position.y << " " << (*element)->position.z << endl;
                        }
                    }
                }
            }
        }


        void calculate_total_number_of_particles_in_cell_list()
        {
            int counter = 0;
            for(int i = 0; i < number_of_cells_x; i++)
            {
                for(int j = 0; j < number_of_cells_y; j++)
                {
                    for(int k = 0; k < number_of_cells_z; k++)
                    {
                        counter += cell_list[i][j][k].size();
                    }
                }
            }
            cout << "Number of particles in cell list: " << counter << endl;
        }
};