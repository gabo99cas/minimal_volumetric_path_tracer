//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef SPHERE_H
#define SPHERE_H

#include "Ray.h"
#include <vector>


class Sphere {
public:
    double r;	// radio de la esfera
    Point p;	// posicion
    Color c;	// color
    Color radiance; //radiancia > 0 si es fuente de luz
    int material; //material 0 lambertiano 1 microfacet
    Color eta;
    Color kappa;
    double alpha; //aspereza

    Sphere(double r_, Point p_, Color c_, Color radiance, int material, Color eta, Color kappa, double alpha): r(r_), p(p_), c(c_), radiance(radiance), material(material), eta(eta), kappa(kappa), alpha(alpha){}

    // PROYECTO 1
    // determina si el rayo intersecta a esta esfera
    double intersect(const Ray &ray) const {
        // regresar distancia si hay intersección
        double t1,t2;
        double det=(ray.o-this->p).dot(ray.d)*(ray.o-this->p).dot(ray.d)-(ray.o-this->p).dot(ray.o-this->p)+this->r*this->r;
        if (det<0) return 0.0; // regresar 0.0 si no hay interseccion return 0.0;
        t2 = -(ray.o-this->p).dot(ray.d)+sqrt(det);
        t1 = -(ray.o-this->p).dot(ray.d)-sqrt(det);
        if(t1<0 || abs(t1)<0.0001) return t2; //caso especial donde la primera interseccion es invalida
        return t1; //t1 siempre es la primera interseccion con excepcion del caso especial
        //los numeros negativos se pueden filtrar aqui o en bool intersect()
    }
    //version de intersect para vpt, regresa ambas intersecciones
    void intersectVPT(const Ray &ray, double &t1, double &t2) const {
        // regresar distancia si hay intersección
        double det=(ray.o-this->p).dot(ray.d)*(ray.o-this->p).dot(ray.d)-(ray.o-this->p).dot(ray.o-this->p)+this->r*this->r;
        if (det<0) {t1=0.0; t2=0.0; return;} // regresar 0.0 si no hay interseccion return 0.0;
        t2 = -(ray.o-this->p).dot(ray.d)+sqrt(det);
        t1 = -(ray.o-this->p).dot(ray.d)-sqrt(det);
    }

};

extern std::vector<Sphere> spheres;

#endif //SPHERE_H
