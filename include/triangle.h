#pragma once
#include<cmath>
#include<vector>
#include "vec.h"
#include "particle.h"
using namespace std;

class Triangle
{   
    public:
        vector< int > indices_of_particles_in_triangle;
        double rest_length;
        double spring_constant;
        int triangle_index;

        // Constructor of the class
        Triangle(int given_index_1, int given_index_2, int given_index_3, double given_restlength, double given_spring_constant, int given_triangle_index):
            rest_length(given_restlength), spring_constant(given_spring_constant), triangle_index(given_triangle_index)
        {
            indices_of_particles_in_triangle.push_back(given_index_1);
            indices_of_particles_in_triangle.push_back(given_index_2);
            indices_of_particles_in_triangle.push_back(given_index_3);
        };
            
        // Destructor of the class
        ~Triangle() {};

        void change_spring_constant(double given_spring_constant)
        {
            spring_constant = given_spring_constant;
        }

        void change_rest_length(double given_rest_length)
        {
            rest_length = given_rest_length;
        }


        void print_spring_constant()
        {
            cout << "Spring constant of triangle: " << spring_constant << endl;
        }

};
