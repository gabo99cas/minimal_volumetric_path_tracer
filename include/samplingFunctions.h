//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef SAMPLINGFUNCTIONS_H
#define SAMPLINGFUNCTIONS_H

#include "microFacetUtilities.h"
#include "pathTracingUtilities.h"

//devuelve wi con muestreo hemisferico uniforme
inline Vector uniformHemispheric(Vector n){

		double theta = acos(erand48(seed));
		double phi = 2*M_PI*(erand48(seed));

		//es direccion por lo tanto r = 1 (direcciones locales)
		double x1 = sin(theta)*cos(phi);
		double y1 = sin(theta)*sin(phi);
		double z1 = cos(theta);

		Vector s1,t1 = Vector(); //la normal ya esta dada por el punto x
		coordinateSystem(n, s1, t1);
		Vector wi = s1*x1+t1*y1+n*z1;

		wi.normalize();
		return wi;
}

//devuelve wi con muestreo esferico uniforme
inline Vector uniformSpheric(){
		double theta = acos(1-2*erand48(seed));
		double phi = 2*M_PI*erand48(seed);

		//es direccion por lo tanto r = 1 (direcciones locales)
		double x1 = sin(theta)*cos(phi);
		double y1 = sin(theta)*sin(phi);
		double z1 = cos(theta);

		Vector wi = Vector(x1,y1,z1);

		wi.normalize();
		return wi;
}

//muestreo de coseno hemisferico
inline Vector cosineHemispheric(Vector n){ //n es z para la direccion local muestreada
		double theta = acos(sqrt(1-erand48(seed)));
		double phi = 2*M_PI*erand48(seed);

		//es direccion por lo tanto r = 1 (direcciones locales)
		double x1 = sin(theta)*cos(phi);
		double y1 = sin(theta)*sin(phi);
		double z1 = cos(theta);

		Vector s1,t1 = Vector(); //la normal ya esta dada por el punto x
		coordinateSystem(n, s1, t1);
		Vector wi = s1*x1+t1*y1+n*z1;

		wi.normalize();
		return wi;
}

inline Vector solidAngle(Vector wc, double costheta_max){
	double e0 = erand48(seed);
	double theta = acos((1-e0)+e0*costheta_max);
	double phi = 2*M_PI*erand48(seed);

	//convertir el angulo a coordenadas locales
	double x1 = sin(theta)*cos(phi);
	double y1 = sin(theta)*sin(phi);
	double z1 = cos(theta);

	//trasladar a coordenadas globales
	Vector s1,t1 = Vector(); //la normal ya esta dada por el punto x
	coordinateSystem(wc, s1, t1);
	Vector wi = s1*x1+t1*y1+wc*z1;

	wi.normalize();
	return wi;
}

//devuelve la probabilidad para muestreo de angulo solido
inline double solidAngleProb(double costheta_max){
	return 1/(2*M_PI*(1-costheta_max));
}



//devuelve la probabilidad para muestreo de coseno hemisferico
inline double hemiCosineProb(double cosine){
	return cosine*1/M_PI;
}

//muestreo de materiales microfacet (el modelo emplea coordenadas locales)
inline Color microfacet(Point x, Vector wo, Vector wh, Vector n, Sphere &obj, double alpha, int &idsource){
	Color Le, fr;
	Vector s, t;
	Vector nlocal = Vector(0,0,1);
	int sourceid;
	wo = wo*(-1); //se invierte para respetar la convencion de direcciones salientes
	//wo debe pasarse a local para que todos los vectores esten en el mismo marco de referencia n=(0,0,1)
	coordinateTraspose(n, wo);
	wo.normalize();
	Vector wi = wo*(-1)+wh*2*(wh.dot(wo)); //reflexion especular, direccion reflejada
	wi.normalize();

	coordinateSystem(n, s, t);
	Vector wiglobal = s*wi.x+t*wi.y+n*wi.z;
	wiglobal.normalize();
	//lanzar un rayo direccion wi, debe globalizarse para resolver visibilidad
	Le = rayTracer(x, wiglobal, sourceid);
	idsource = sourceid;
	fr = frMicroFacet(obj.eta, obj.kappa, wi, wh, wo, alpha, nlocal);

	return Le.mult(fr)*Vector(0,0,1).dot(wi)*(1/microFacetProb(wo, wh, alpha, nlocal));
}



//L muestreo de area
inline Color areaLight(Point center, double radio, Point x, Color radiance, const Sphere &obj, Vector n, Vector wo, double alpha){
	Vector xl;
	Point light;
	Color le, fr, L;
	Vector aux = uniformSpheric(); //respecto a z
	light = center + aux*radio; //punto muestreado
	xl = x-light;
	if(aux.dot(xl.normalize())<0)
		return 0;
	if(visibility(light, x))
		le = radiance; //radiancia de la fuente de luz
	else
		return 0;
	Vector wilocal = xl*-1;
	Vector wolocal = wo*-1;

	coordinateTraspose(n, wolocal);
	coordinateTraspose(n, wilocal);
	wolocal.normalize();
	wilocal.normalize();

	Vector wh = wilocal+wolocal;
	wh.normalize();

	if(obj.material==0)	fr = obj.c*(1/M_PI);
	else{
		fr = frMicroFacet(obj.eta, obj.kappa, wilocal, wh, wolocal, 0.3, Vector(0,0,1));
	}
	double prob = (light-x).dot(light-x)/(4*M_PI*radio*radio*aux.dot(xl));
	L = le.mult(fr)*n.dot((light-x).normalize())*(1/prob);
	return L;
}

//devuelve color calculado para muestreo de area
inline Color muestreoArea(const Sphere &source, const Sphere &obj, Point x, Vector n, Vector wray, double alpha){
	Color L = areaLight(source.p, source.r, x, source.radiance, obj, n, wray, 0.3);
	return L;
}

//L de Solid Angle
inline Color solidAngle(Vector n, Vector cx, Vector wray, double costheta_max, Point x, int indice, const Sphere &obj, Vector &aux, double alpha){
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

	return L;
}

//muestreo de dielectrico suave
inline Color softDielectric(double etat, double etai, Vector wi, Vector n, Point x, int &idsource){ //wi sera la direccion incidente en este modelo (r.d)
	Color Ld;
	int sourceid;
	//direcciones de transmision para calcular Fresnel
	Vector wt = refraxDielectric(etai, etat, wi, n);
	wt.normalize();
	double F = fresnelDie(etai, etat, n.dot(wt), n.dot(wi));
	//muestreo
	if(erand48(seed)<F){
		//calcula reflexion
		Vector wr = reflexDielectric(wi, n);
		wr.normalize();
		//con reflexion
		//la probailidad de muestreo es F, F esta en el numerador ... 1
		Ld = rayTracer(x, wr, sourceid)*(1/abs(n.dot(wr)));

	}
	else{
		//con transmision
		double ratio = etat/etai;
		//la probabilidad de muestreo es 1-F, 1-F esta en el denominador ... 1
		Ld = rayTracer(x, wt, sourceid)*(1/abs(n.dot(wt)))*ratio*ratio;

	}
	idsource = sourceid;
	return Ld;
}

//muestreo de angulo solido
inline Color muestreoSA(Sphere &source, Point x, int indice, const Sphere &obj, Vector n, Vector wray, Vector &aux, double &costhetaMax, double alpha){
	Color L;
	Vector cx = source.p-x;
	double normcx = sqrt(cx.dot(cx));
	cx = cx*(1/normcx);
	double costheta_max = sqrt(1-(source.r/normcx)*(source.r/normcx));
	costhetaMax = costheta_max;
	L = solidAngle(n, cx, wray, costheta_max, x, indice, obj, aux, alpha);
	return L;
}

//muestreo uniforme
inline Color uniform(Vector n, Point x, Color BDRF, Vector &aux, int &idsource){ //aux obtiene de regreso la direccion muestreada
	int sourceid;
	Color L, Le;
	Vector wi;
	wi = cosineHemispheric(n);
	wi.normalize();
	Le = rayTracer(x, wi, sourceid);
	idsource = sourceid;
	L = L + Le.mult(BDRF*(1/M_PI))*n.dot(wi)*(1/hemiCosineProb(n.dot(wi)));
	aux = wi;
	return L;
}


#endif //SAMPLINGFUNCTIONS_H
