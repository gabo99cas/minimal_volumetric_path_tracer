//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef PATHTRACINGUTILITIES_H
#define PATHTRACINGUTILITIES_H

#include "Sphere.h"



inline bool intersect(const Ray &r, double &t, int &id) {
    double tmin = __DBL_MAX__; // distancia mínima inicializada al máximo valor
    double tact;              // distancia actual a la esfera
    int contact = 0;          // contador de contactos con esferas

    for (size_t i = 0; i < spheres.size(); i++) { // Iterar sobre las esferas en el vector
        tact = spheres[i].intersect(r);          // Calcular intersección con la esfera actual

        if (tact > 0 && abs(tact) > 0.0001) {    // Verificar intersección válida
            contact++;                           // Incrementar contador de contactos
            if (tact < tmin) {                   // Actualizar la distancia mínima y el índice
                tmin = tact;
                id = i;
            }
        }
    }

    if (contact > 0) {  // Si hay al menos una intersección
        t = tmin;       // Actualizar la distancia mínima de intersección
        return true;
    } else {            // Si no hubo intersecciones
        t = 0;
        return false;
    }
}

//resuelve la visibilidad para dos puntos
bool visibility(Point light, Point x){
    Vector lx = light-x; //vector con direccion a la fuente
    double distance = sqrt(lx.dot(lx)); //distancia a la fuente
    lx.normalize();
    lx=lx*-1; //invierte la direccion
    Ray r2 = Ray(light, lx);
    int id=0;
    double t;
    intersect(r2, t, id);
    if(t>distance || t==0)
    {
        return true; //devolver radiancia
    }
    else return false;
}

//resuelve la visibilidad de forma implicita (no sirve para muestreo de luz)
inline Color rayTracer(Point x, Vector wi, int &sourceid){
    Ray r1 = Ray(x, wi);
    double t;
    int id = 0;
    if(!intersect(r1, t, id)) return Color();
    sourceid = id;
    const Sphere &source = spheres[id];
    return source.radiance;
}

inline double cosinethetaMax(int sourceid, Point x){
    double radio = spheres[sourceid].r;
    Vector cx = spheres[sourceid].p-x;
    double normcx = sqrt(cx.dot(cx));
    cx.normalize();
    double costheta_max = sqrt(1-(radio/normcx)*(radio/normcx));
    return costheta_max;
}

#endif //PATHTRACINGUTILITIES_H
