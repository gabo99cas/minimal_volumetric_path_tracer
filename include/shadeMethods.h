//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef SHADEMETHODS_H
#define SHADEMETHODS_H

#include "misSamplingFunctions.h"
#include "mathUtilities.h"
#include "vptShadeMethods.h"

//L de Solid Angle en ray marching
//CUIDADO FUNCION INCOMPLETA
inline Color solidAngleMarching(Vector n, Vector cx, Vector wray, double costheta_max, Point x, int indice, const Sphere &obj, Vector &aux, double alpha){
    Color L , Le;
    Vector wi;
    Color fr;
    Vector wolocal = wray*-1;
    Vector wilocal;
    wi = solidAngle(cx, costheta_max); //wi guarda la direccion muestreada
    aux = wi;
    Ray r1 = Ray(x, wi); //lanzar el rayo desde x con direccion wi
    wilocal = wi;

    coordinateTraspose(n, wilocal);
    coordinateTraspose(n, wolocal);

    wilocal.normalize();//el cambio en el sistema de coordenadas altera la magnitud del vector
    wolocal.normalize();

    Vector wh = (wilocal+wolocal);
    wh.normalize();

    //para aluminio eta =[1.66058, 0.88143, 0.521467] kappa = [9.2282, 6.27077, 4.83803]
    //para oro eta = [0.143245, 0.377423, 1.43919] kappa = [3.98479, 2.3847, 1.60434]
    if(obj.material==0){
        fr = obj.c*(1/M_PI);
    }
    else if(obj.material==2){
        //la evaluacion de un dielectrico suave siempre dara 0 en muestreo de luz, no es factible dar con la direccion correcta
        fr = Color(0,0,0);
    }
    else fr = frMicroFacet(obj.eta, obj.kappa, wilocal, wh, wolocal, alpha, Vector(0,0,1));
    double t;
    int id = 0;
    intersect(r1, t, id);
    const Sphere &source = spheres[id];

    if(indice==id) Le=source.radiance; //nos aseguramos que exista visibilidad entre la fuente y no otras fuentes de luz
    else Le = Color();
    //funcion L
    L = Le.mult(fr)*n.dot(wi)*(1/solidAngleProb(costheta_max));

    return Le;
}

//calcula la fs y la nueva direccion para un camino en Path tracing. Lind
inline Color BDSF(Vector &aux, Vector wray, Vector n,  double &prob, int id){
    Vector wi;
    Color fs1;
    Vector wo = wray*-1;
    if(spheres[id].material==0){
        wi = cosineHemispheric(n);
        fs1 = spheres[id].c*(1/M_PI);
        prob = hemiCosineProb(n.dot(wi));
        aux = wi;
    }
    else if(spheres[id].material == 2){
        Vector wt = refraxDielectric(1.0, 1.5, wo, n);
        wt.normalize();
        double F = fresnelDie(1.0, 1.5, n.dot(wt), n.dot(wo));
        if(erand48(seed)<F){
            //se obtiene reflexion
            wi = reflexDielectric(wo, n);
            wi.normalize();
            fs1=Color(1,1,1)*(1/n.dot(wi))*(F);

            prob = F;
        }
        else{
            wi = wt;
            fs1=Color(1,1,1)*(1/n.dot(wi))*(1-F)*(1.5)*1.5;

            prob = 1-F;

        }
        aux = wi;
    }
    else if(spheres[id].material == 1){
        double alpha = spheres[id].alpha;
        Vector wh = vectorFacet(alpha); //local
        Vector s, t;
        coordinateSystem(n, s, t);
        wh = s*wh.x+t*wh.y+n*wh.z; //global
        wi = wo*(-1)+wh*2*(wh.dot(wo));
        fs1 = frMicroFacet(spheres[id].eta, spheres[id].kappa, wi, wh, wo, alpha, n);
        prob = microFacetProb(wo, wh, alpha, n);
        aux = wi;
    }
    return fs1;
}

//Path tracer explicito iterativo con MIS y ruleta rusa (Version definitiva, las anteriores quedan como legavy)
inline Color iterativePathTracer(Ray r){
    Color Accum = Color();
    Color fs = Vector(1,1,1);
    Color Ld, fs1 = Color();
    Vector wi;
    int bounces = 0;
    double factor = 1;
    double q = 0.4;
    double continueprob = 1.0 - q;
    double prob;
    int i = 0;
    while(true){ //el ciclo se interrumpe en puntos especificos, se puede agregar un contador para limitar los caminos
        double t;
        int id = 0;
        if (!intersect(r, t, id))
            break;	// el rayo no intersecto objeto, return Vector() == negro


        if(spheres[id].radiance.x > 0){ //si se impacta una fuente de luz, regresa radiancia
            if(bounces<1) return spheres[id].radiance; //esto solo aplica para el primer rayo que no hace muestreo explicito
            break;
        }


        Point x = r.o + r.d*t;
        Vector n = (x-spheres[id].p);
        n.normalize();
        //iluminacion directa
        //for each light source
        int size = spheres.size();  // Devuelve el número de elementos en el vector

        for(int light = 0; light<size; light++)
        {
            //resuelve las fuentes puntuales
            if(spheres[light].r==0) Ld = pLight(spheres[id], x, n, r.d, spheres[light].radiance, spheres[light].p, spheres[id].alpha) + Ld;

        }
        Ld = MIS(spheres[id], x, n, r.d, spheres[id].alpha) + Ld; //mis tiene su propio ciclo, por eso esta fuera del for
        //ruleta rusa
        if(erand48(seed) < q)
        {
            //Accum = Accum+ fs.mult(Ld)*factor;
            break;
        }
        //muestreo de BDSF para obtener Lind
        fs1 = BDSF(wi, r.d, n, prob, id);

        Ray newray = Ray(x, wi);
        r.o = newray.o;
        r.d = newray.d;
        double cosine = n.dot(wi);
        Accum = Accum + fs.mult(Ld)*factor;
        fs = fs.mult(fs1);
        factor = factor*cosine*(1/(prob*continueprob));
        bounces++;
        Ld = Color();
    }
    return Accum;

}

#endif //SHADEMETHODS_H
