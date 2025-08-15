#pragma once
#include<cmath>
#include "vec.h"

class Particle
{   
    public:
        vec old_position;
        vec position;
        double radius;
        int triangle_index;
        int particle_index;
        
        // Constructor of the class
        Particle(vec given_old_position, vec given_position, double given_radius, int given_particle_index, int given_triangle_index)
            : old_position(given_old_position), position(given_position), radius(given_radius), particle_index(given_particle_index),
            triangle_index(given_triangle_index) {};

        // Destructor of the class
        ~Particle() {};

        vec give_position() const
        {
            return position;
        }

        vec give_old_position() const
        {
            return old_position;
        }

        double give_radius() const
        {
            return radius;
        }
};