//
// Created by Gabriel Castillo on 29/01/25.
//

#ifndef VPTSAMPLINGFUNCTIONS_H
#define VPTSAMPLINGFUNCTIONS_H

#include "pathTracingUtilities.h"

//muestreo sobre la distancia (para medio homogeneo)
inline double freeFlightSample(double sigma_t){
    //generar un numero aleatorio
    double xi = erand48(seed);
    //calcular la distancia de vuelo libre
    return -log(1-xi)/sigma_t;
}


//pdf del free flight sampling (para medio homogeneo)
inline double freeFlightProb(double sigma_t, double d){
    return sigma_t*exp(sigma_t*d*-1.0);
}
//pdf success = 1 - pdf failure
inline double pdfSuccess(double sigma_t, double tmax){
    return 1.0-exp(-sigma_t*tmax);
}

//pdf failure = transmitance (para medio homogeneo)
inline double pdfFailure(double sigma_t, double tmax){
    return exp(-sigma_t*tmax);
}

//muestreo de fase isotropica
inline Vector isotropicPhaseSample(){
    //generar dos numeros aleatorios
    double xi1 = erand48(seed);
    double xi2 = erand48(seed);
    //calcular el angulo theta
    double theta = acos(1-2*xi1);
    //calcular el angulo phi
    double phi = 2*M_PI*xi2;
    //devolver el vector de la direccion del rayo
    Vector wi = Vector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    wi.normalize();
    return wi;
}

//pdf de la fase isotropica
inline double isotropicPhaseProb(){
    return 1.0/(4.0*M_PI);
}

//equi-angular sampling
inline double equiAngularSample(double D, double thetaA, double thetaB){
    double xi = erand48(seed);
    return D*tan((1-xi)*thetaA + xi*thetaB);
}

//pdf equi-angular sampling
inline double equiAngularProb(double D, double thetaA, double thetaB, double sample_t){
    return D/fabs(thetaB-thetaA)/(sample_t*sample_t+D*D);
}

#endif //VPTSAMPLINGFUNCTIONS_H
