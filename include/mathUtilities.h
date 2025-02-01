//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef MATHUTILITIES_H
#define MATHUTILITIES_H

#include "Vector.h"

inline void coordinateSystem( Vector &n, Vector &s, Vector &t) {
    if (std::abs(n.x) > std::abs(n.y)) {
        double invLen = 1.0 / std::sqrt(n.x * n.x + n.z * n.z);
        t = Vector(n.z * invLen, 0.0f, -n.x * invLen);
    }	else {
        double invLen = 1.0 / std::sqrt(n.y * n.y + n.z * n.z);
        t = Vector(0.0f, n.z * invLen, -n.y * invLen);
    }
    s = t%n;
}

inline void coordinateTraspose(Vector n, Vector &w){ //pasa w direccion local con base a la normal
    Vector wlocal = w;
    Vector s,t;
    coordinateSystem(n, s, t);
    Vector ninv, sinv, tinv;
    sinv = Vector(s.x, t.x, n.x);
    tinv = Vector(s.y, t.y, n.y);
    ninv = Vector(s.z, t.z, n.z);
    w = sinv*wlocal.x+tinv*wlocal.y+ninv*wlocal.z; //vector en local
}


// limita el valor de x a [0,1]
inline double clamp(const double x) {
    if(x < 0.0)
        return 0.0;
    else if(x > 1.0)
        return 1.0;
    return x;
}

// convierte un valor de color en [0,1] a un entero en [0,255]
inline int toDisplayValue(const double x) {
    return int( pow( clamp(x), 1.0/2.2 ) * 255 + .5);
}

#endif //MATHUTILITIES_H
