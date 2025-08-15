#pragma once
#include<cmath>

// The class "vec" refers to a standard three-dimensional vector and contains all of the standard arithmetic operations.
class vec
{       
    public:
        double x, y, z;

        // Constructor of the class
        vec(double given_x, double given_y, double given_z) : x(given_x), y(given_y), z(given_z) {}; 
        vec() : x(0.0), y(0.0), z(0.0) {}; 

        // Destructor of the class
        ~vec() {};

        // Overloading the "=" operator
        void operator=(const vec& given_vector)
        {
            x = given_vector.x;
            y = given_vector.y;
            z = given_vector.z;
        }

        // Overloading the "+" operator
        vec operator+(const vec& given_vector) const
        {
            return vec(x + given_vector.x, y + given_vector.y, z + given_vector.z);
        }

        // Overloading the "-" operator
        vec operator-(const vec& given_vector) const
        {
            return vec(x - given_vector.x, y - given_vector.y, z - given_vector.z);
        }

        // Overloading the "*" operator
        double operator*(const vec& given_vector) const
        {
            return x * given_vector.x + y * given_vector.y + z * given_vector.z;
        }

        // Overloading the "*" operator
        vec operator*(const double& scalar) const
        {
            return vec(scalar * x, scalar * y, scalar * z);
        }

        // Calculate the euclidian norm of a vector
        double norm() const
        {
            return sqrt(x * x + y * y + z * z);
        }

        // Normalize a vector
        void normalize()
        {
            double norm_of_vector = this->norm();
            x /= norm_of_vector;
            y /= norm_of_vector;
            z /= norm_of_vector;
        }

        double norm_squared() const
        {
            return x*x + y*y + z*z;
        }

        void reflect_vector(const vec& normal_vector)
        {
            x = x - normal_vector.x * 2 * (x * normal_vector.x);
            y = y - normal_vector.y * 2 * (y * normal_vector.y);
            z = z - normal_vector.z * 2 * (z * normal_vector.z);            
        }
};