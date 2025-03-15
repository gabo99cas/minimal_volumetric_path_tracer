//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef VOLUMETRICBASICFUNCTIONS_H
#define VOLUMETRICBASICFUNCTIONS_H

#include "mathUtilities.h"
#include "samplingFunctions.h"
#include "pathTracingUtilities.h"
#include "vptSamplingFunctions.h"

//funciones para volumenes
inline double transmitance(Point x1, Point x2, double sigma_t){
    //obten la distancia entre dos puntos
    Vector aux = x2-x1;
    double d2 = (aux).dot(aux);
    double d = sqrt(d2);
    //devuelve la transmitancia
    return exp(sigma_t*d*-1.0);
}

//funcion para calcular la suma de las transmitacias entre dos puntos
//se debe considerar que pueden existir objetos volumetricos en el camino y espacios vacios
//por ahora todos tienen la misma sigma_t, pero se debe tomar de la interseccion
inline double multipleT(Point x1, Point x2 , double sigma_t){
    //las transmitancias se multiplican
    double T = 1;
    double t_cercano;
    int id = 0;
    double t_aux, t_aux2;
    //lanzar un rayo entre x1 y x2
    Vector w = x2-x1;
    w.normalize();
    Ray r = Ray(x1, w);
    //calcula la interseccion con cada esfera volumetrica, si no hay devuelve 1
    //for each volumetric object
    int size = spheres.size();  // Devuelve el número de elementos en el vector
    for(int i = 0; i<size; i++){
        //si es tipo 3
        if(spheres[i].material==3){
            spheres[i].intersectVPT(r,t_aux,t_aux2 ); //tiene que ser una version de intersect que regrese t1 y t2
            if(t_aux2<0){
                T = T*exp(-sigma_t*t_aux);
            }
            //se debe cumplir que t1-t2>0
            if(t_aux2-t_aux>0){
                //calcular la transmitancia entre x1 y x2
                T = T*exp(-sigma_t*(t_aux2-t_aux));
            }

        }

    }

    return T;
}

inline double isotropicPhaseFunction(){
    //devuelve la funcion de fase isotropica
    return 1/(4*M_PI);
}

inline bool intersectVPT(const Ray &r, double &t, int &id) {
    double tmin = __DBL_MAX__;  // Distancia mínima inicializada al valor máximo
    double tact;                // Distancia actual a la esfera
    int contact = 0;            // Contador de contactos con esferas

    for (size_t i = 0; i < spheres.size(); i++) {  // Iterar sobre el vector de esferas
        if (spheres[i].material != 3) {  // Ignorar esferas de tipo 3 (volumétricas)
            tact = spheres[i].intersect(r);  // Calcular la intersección con la esfera
            if (tact > 0 && abs(tact) > 0.0001) {  // Verificar intersección válida
                contact++;  // Se cambia el estado a que el contacto ocurrió
                if (tact < tmin) {  // Si la intersección es más cercana, actualizamos
                    tmin = tact;
                    id = i;
                }
            }
        }
    }

    if (contact > 0) {  // Si hubo algún contacto con una esfera válida
        t = tmin;       // Actualizar la distancia mínima de intersección
        return true;
    } else {            // Si no hubo intersección
        t = 0;
        return false;
    }
}

//visibilidad para vpt
inline bool visibilityVPT(Point light, Point x){
    Vector lx = light-x; //vector con direccion a la fuente
    double distance = sqrt(lx.dot(lx)); //distancia a la fuente
    lx.normalize();
    lx=lx*-1; //invierte la direccion
    Ray r2 = Ray(light, lx);
    int id=0;
    double t;
    intersectVPT(r2, t, id);
    if(t>distance || t==0)
    {
        return true; //devolver radiancia
    }
    else return false;
}

//intersect v2, devulve la primera y segunda interseccion de la esfera mas cercana
inline bool intersectV2(const Ray &r, double &t1, double &t2, int &id) {
    double tmin = __DBL_MAX__;
    double tact1, tact2; //distancia actual a la esfera
    int contact=0;

    int elements = spheres.size();  // Devuelve el número de elementos en el vector

    for(int i=0; i<elements; i++){
        spheres[i].intersectVPT(r, tact1, tact2);
        if(tact1>0 && abs(tact1)>0.0001){ //abs es necesario para eliminar intersecciones con el interior de las esfera dadas por errores de precision
            contact++; //se cambia el estado a que el contacto sucedio
            if(tact1<tmin){
                tmin=tact1;
                t1=tact1;
                t2=tact2;
                id=i;
            }	}	}
    if(contact>0){
        return true;
    }
    else{
        t1=0;
        t2=0;
        return false;
    }
}

//calcular interseccion tomando en cuenta que la escena esta inmersa en un medio
inline bool intersectMedium(const Ray &r, double &tSurface, int &id, double &tMedium, bool &hitMedium, double sigma_t)
{
    double tmin = __DBL_MAX__;
    double tact;
    int contact = 0;
    int size = spheres.size();


    // Regular surface intersection
    for (int i = 0; i < size; i++) {
        tact = spheres[i].intersect(r);
        if (tact > 0 && abs(tact) > 0.0001) {
            contact++;
            if (tact < tmin) {
                tmin = tact;
                id = i;
            }
        }
    }

    tSurface = (contact > 0) ? tmin : __DBL_MAX__;

    // Medium interaction (only check for the enclosing sphere)
    hitMedium = false;
    tMedium = freeFlightSample(sigma_t);  // Sample a distance inside the medium

    // If the medium interaction happens before the surface, register it
    if (tMedium < tSurface) {
        hitMedium = true;
    }

    return contact > 0 || hitMedium;

}

//funcion utilitaria para Sampleo de phase
inline Vector samplePhaseFunction(Vector &wi, const Vector &wo, double &prob) {
    wi = isotropicPhaseSample(); // Pick a random direction
    prob = 1.0 / (4.0 * M_PI); // Isotropic PDF
    return {1, 1, 1}; // No color change
}

//recibe el indice de la fuente de luz y regresa x0, D, thetaA y thetaB
inline bool equiAngularParams(int idsource, Point x, Point &x0, Ray r, double &D, double &thetaA, double &thetaB){
    Point c = spheres[idsource].p;
    //calcular la proyeccion ortogonal de c en el rayo, punto mas cercano a c en el rayo, x0
    x0 = r.o + r.d*((c-r.o).dot(r.d)/r.d.dot(r.d));
    //el parametro x es x_s
    //verificar si el punto x0 esta entre r.o y x
    if((x0-r.o).dot(r.d)<0)
    {
        x0 = r.o;
        //return false;
    }
    if((x0-x).dot(r.d)>0){
        x0 = x;
        //return false;

    }


    //calcular la magnitud de x0-c
    D = sqrt((x0-c).dot(x0-c));
    //theta A y theta B son los angulos de apertura de la fuente de luz, intervalo de integracion, en este caso usaremos r.o y x como puntos de integracion
    //calcular los lados del triangulo rectangulo, comparten D
    const double A = sqrt((x0-r.o).dot(x0-r.o))*-1;
    const double B = sqrt((x-x0).dot(x-x0));
    //calcular el angulo theta A y theta B
    thetaA = atan2(A,D);
    thetaB = atan2(B,D);
    return true;
}

//resuelve la iluminación por single scattering para un punto en particular previamente muestreado
inline Color singleScattering(Point xt, int idsource, double sigma_t, double sigma_s, double transmitanceXT, double probSource)
{
    Color Ld;
    //calcular la dirección a la fuente de luz o light
    //se debe usar angulo sólido o dirección única en caso de fuente puntual

    //PARA FUENTE PUNTUAL
    //calculando un single scattering
    if (spheres[idsource].r == 0) {
        if (const Point light = spheres[idsource].p; visibility(light, xt)) {
            Color Le = spheres[idsource].radiance;
            const double distanceLight = (light - xt).dot(light - xt);
            Le = Le * (1 / distanceLight);
            const Color Ls = Le * transmitance(xt, light, sigma_t) * isotropicPhaseFunction();
            Ld = Ls * transmitanceXT * sigma_s * (1 / probSource);
        }

    }

    //para fuentes no puntuales usar angulo solido
    //todo ES POSIBLE AGREGAR LA RTUTINA EN UNA FUNCIÓN
    //obtener el vector que conecta el centro de la fuente con xt, wc

    Vector wc = spheres[idsource].p-xt;
    //obtener su norma
    double wc_magnitud = sqrt(wc.dot(wc));
    wc = wc*(1/wc_magnitud); //aprovechamos para normalizarlo
    //calcular el angulo del cono de visibilidad
    double costheta_max = sqrt(1-(spheres[idsource].r/wc_magnitud)*(spheres[idsource].r/wc_magnitud));
    //obtener la dirección en el cono, llamémosla wl por ser single scattering
    const Vector wl = solidAngle(wc, costheta_max);
    //calcular la probabilidad para dicha dirección
    const double prob_wl = solidAngleProb(costheta_max);
    //integrar usando la misma expresión de single scattering Ls
    //todo si hay visibilidad
    Color Le = spheres[idsource].radiance;
    Point light = spheres[idsource].p; //TODO recalcularlo en la superficie de la fuente para que el calculo de la transmitancia sea correcto
    //TODO para hacer el recalculo aprovecha que tienes que verificar la visibilidad
    //evaluar la contribución single scattering Ls
    const Color Ls = Le * transmitance(xt, light, sigma_t) * isotropicPhaseFunction();
    //aplicar el resto de componentes para calcular la iluminación directa en el medio
    Ld = Ls * transmitanceXT * sigma_s * (1/prob_wl)* (1 / probSource);

    return Ld;
}

#endif //VOLUMETRICBASICFUNCTIONS_H
