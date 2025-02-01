//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>

class Vector {
	public:
    double x, y, z; // coordenadas x,y,z

    // Constructor del vector, parametros por default en cero
    Vector(double x_= 0, double y_= 0, double z_= 0){ x=x_; y=y_; z=z_; }

    // operador para suma y resta de vectores
    Vector operator+(const Vector &b) const { return Vector(x + b.x, y + b.y, z + b.z); }
    Vector operator-(const Vector &b) const { return Vector(x - b.x, y - b.y, z - b.z); }
    // operator multiplicacion vector y escalar
    Vector operator*(double b) const { return Vector(x * b, y * b, z * b); }

    // operator % para producto cruz
    Vector operator%(Vector&b){return Vector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);}

    // producto punto con vector b
    double dot(const Vector &b) const { return x * b.x + y * b.y + z * b.z; }

    // producto elemento a elemento (Hadamard product)
    Vector mult(const Vector &b) const { return Vector(x * b.x, y * b.y, z * b.z); }

    // normalizar vector
    Vector& normalize(){ return *this = *this * (1.0/sqrt(x * x + y * y + z * z)); }
};
typedef Vector Point;
typedef Vector Color;

extern unsigned short seed[3];

#endif //VECTOR_H
