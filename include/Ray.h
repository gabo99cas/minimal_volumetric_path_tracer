//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef RAY_H
#define RAY_H

#include "Vector.h"

class Ray {
    public:
    Point o;
    Vector d; // origen y direcccion del rayo
    Ray(Point o_, Vector d_) : o(o_), d(d_) {} // constructor
};



#endif //RAY_H
