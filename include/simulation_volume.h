#pragma once
#include<cmath>

// The class "Simulationvolume" refers to a two- or three-dimensional container. 
class Simulationvolume
{
    private:
        // Dimensionality of the simulation volume
        int dimensionality;
        // Size of the simulation volume
        double boxlength_x, boxlength_y, boxlength_z;
        // Specification of boundary conditions (true = periodic)
        bool periodic_boundary_conditions[3];

    public:
        // Constructor of the class
        Simulationvolume(int given_dimensionality, double given_size_of_simulationsvolume[3], bool given_periodic_boundary_conditions[3])
            : dimensionality(given_dimensionality), boxlength_x(given_size_of_simulationsvolume[0]), boxlength_y(given_size_of_simulationsvolume[1]), 
                boxlength_z(given_size_of_simulationsvolume[2]) 
        {
            periodic_boundary_conditions[0] = given_periodic_boundary_conditions[0];
            periodic_boundary_conditions[1] = given_periodic_boundary_conditions[1];
            periodic_boundary_conditions[2] = given_periodic_boundary_conditions[2];
        };

        // Destructor of the class
        ~Simulationvolume() {};

        vec give_size_of_simulation_volume()
        {
            vec size_simulation_volume(boxlength_x, boxlength_y, boxlength_z);
            return size_simulation_volume;
        }

        vec calculate_periodic_difference_vector(vec start, vec end)
        {
            vec difference = end - start;
            if(difference.x < -boxlength_x/2.)
            {
                difference.x += boxlength_x;
            }
            if(difference.x >= boxlength_x/2.)
            {
                difference.x -= boxlength_x;
            }
            if(difference.y < -boxlength_y/2.)
            {
                difference.y += boxlength_y;
            }
            if(difference.y >= boxlength_y/2.)
            {
                difference.y -= boxlength_y;
            }
            if(difference.z < -boxlength_z/2.)
            {
                difference.z += boxlength_z;
            }
            if(difference.z >= boxlength_z/2.)
            {
                difference.z -= boxlength_z;
            }
            return difference;
        }

        vec calculate_size_of_vector_boundary_conform(vec given_vector)
        {
            if(given_vector.x < 0)
            {
                given_vector.x += boxlength_x;
            }
            if(given_vector.x >= boxlength_x)
            {
                given_vector.x -= boxlength_x;
            }
            if(given_vector.y < 0)
            {
                given_vector.y += boxlength_y;
            }
            if(given_vector.y >= boxlength_y)
            {
                given_vector.y -= boxlength_y;
            }
            if(given_vector.z < 0)
            {
                given_vector.z += boxlength_z;
            }
            if(given_vector.z >= boxlength_z)
            {
                given_vector.z -= boxlength_z;
            }
            return given_vector;        
        }
};
