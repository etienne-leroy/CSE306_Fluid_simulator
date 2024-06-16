#ifndef VEC3_H
#define VEC3_H

#include <cmath>

class vec3 {
public:
    explicit vec3(double x=0, double y=0, double z=0) {e[0]=x;e[1]=y;e[2]=z;}   // Initialize coordinates

    double length() const {return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];}   
    double length_sqrt() const {return sqrt(length());}     // vector length

    void norm() {
        double n = length_sqrt();
        if(n > 0) {e[0]/=n;e[1]/=n;e[2]/=n;}
    }

    double operator[](int i) const { return e[i]; }
    double &operator[](int i) { return e[i]; }

    double e[3];
};

// Operators
inline vec3 operator-(const vec3 &a, const vec3 &b) {return vec3(a[0]-b[0], a[1]-b[1], a[2]-b[2]);} // substraction of two vectors 
inline vec3 operator+(const vec3 &a, const vec3 &b) {return vec3(a[0]+b[0], a[1]+b[1], a[2]+b[2]);} // addition of two vectors 

inline vec3 operator+=(vec3 &a, const vec3 &b) {    // In-place addition of two vectors 
    a[0]+=b[0]; a[1]+=b[1]; a[2]+=b[2];
    return a;
}

inline vec3 operator/(const vec3 &a, const double b) {return vec3(a[0]/b, a[1]/b, a[2]/b);} // Vector division by scalar

inline vec3 operator*(const vec3 &a, const vec3 &b) {return vec3(a[0]*b[0], a[1]*b[1], a[2]*b[2]);} // Vector multiplication
inline vec3 operator*(const double a, const vec3 &b) {return vec3(a*b[0], a*b[1], a*b[2]);}     // Scalar / Vector 
inline vec3 operator*(const vec3 &a, const double b) {return vec3(a[0]*b, a[1]*b, a[2]*b);}     // Vector / Scalar


inline vec3 cross(const vec3 &a, const vec3 &b) {return vec3(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);}
inline double dot(const vec3 &a, const vec3 &b) {return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];}


#endif // VEC3_H
