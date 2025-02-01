//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef MICROFACETUTILITIES_H
#define MICROFACETUTILITIES_H

#include "Vector.h"

//calcula el fresnel para solo un espectro
inline double fresnelSpectre(double cosine, double sine, double eta, double kappa){
	double a2b2 = sqrt(((eta)*(eta)-(kappa)*(kappa)-sine*sine)*((eta)*(eta)-(kappa)*(kappa)-sine*sine)+4*(eta)*(eta)*(kappa)*(kappa));
	double a = sqrt(0.5*(a2b2+eta*eta-kappa*kappa-sine*sine));
	double parallel, perpendicular;
	perpendicular = (a2b2+cosine*cosine-2*a*cosine)/(a2b2+cosine*cosine+2*a*cosine);
	parallel = perpendicular*(a2b2*cosine*cosine+sine*sine*sine*sine-2*a*cosine*sine*sine)/(a2b2*cosine*cosine+sine*sine*sine*sine+2*a*cosine*sine*sine);
	return 0.5*(parallel+perpendicular);
}

//calcula la naturaleza de la luz reflejada para modelos con reflexion
inline Color fresnel(double cosine_wh, Vector eta, Vector kappa){ //el argumento es el coseno del angulo con respecto a wh
	//se debe calcular para cada espectro en eta y kappa (x,y,z)
	double sine_wh = sqrt(1-cosine_wh*cosine_wh);
	Color fwh;
	fwh.x = fresnelSpectre(cosine_wh, sine_wh, eta.x, kappa.x);
	fwh.y = fresnelSpectre(cosine_wh, sine_wh, eta.y, kappa.y);
	fwh.z = fresnelSpectre(cosine_wh, sine_wh, eta.z, kappa.z);
	return fwh;
}



//calcula la distribucion NDF, usa el angulo theta_h, cos(theta_h)=wh.dot(n);
inline double NDF(double cosine, double alpha){
	double fac1, fac2, tang, sine;
	if(cosine>=0){
		sine = sqrt(1-cosine*cosine);
		fac1 = M_PI*alpha*alpha*cosine*cosine*cosine*cosine;
		tang = sine/cosine;
		fac2 = exp((-1*tang*tang)/(alpha*alpha));
		return (1/fac1)*fac2;

	}
	else return 0;
}

inline double Gn(Vector n, Vector wv, Vector wh, double alpha){
	double sin = sqrt(1-n.dot(wv)*n.dot(wv));
	double tan = sin/(n.dot(wv)); //tan = sen/cos
	double a = 1/(alpha*tan);
	double num, den;
	if(((wv.dot(wh))/(wv.dot(n)))>0){
		if(a<1.6){
			num = 3.535*a+2.181*a*a;
			den = 1+2.276*a+2.577*a*a;
			return num/den;
		}
		else return 1;
	}
	else return 0;
}

inline double G_smith(Vector n, Vector wi, Vector wo, Vector wh, double alpha){
	double g1, g2;
	g1 = Gn(n, wi, wh, alpha);
	g2 = Gn(n, wo, wh, alpha);
	return g1*g2;
}

//obtener direccion wh
inline Vector vectorFacet(double alpha){
	Vector wh;
	double theta = atan(sqrt(-alpha*alpha*log(1-erand48(seed))));
	double phi = 2*M_PI*erand48(seed);

	double x1 = sin(theta)*cos(phi);
	double y1 = sin(theta)*sin(phi);
	double z1 = cos(theta);

	wh = Vector(x1, y1, z1);
	wh.normalize();

	return wh;
}
//calcula la probabilidad para el modelo microfacet
inline double microFacetProb(Vector wo, Vector wh, double alpha, Vector n){//n debe estar en mismo sistema que wh y wo
	//wo en local
	double num, den;
	num = wh.dot(n);
	den = 4 * abs(wo.dot(wh));
	return NDF(wh.dot(n), alpha)*num/den;
}

//evaluacion de la BDRF para el modelo microfacet
inline Color frMicroFacet(Color eta, Color kappa, Vector wi, Vector wh, Vector wo, double alpha, Vector n){
	//Vector n = (0,0,1);
	//es necesario que los vectores cumplan estar normalizados, ser direcciones salientes y en direccion local
	double den =(4*abs(n.dot(wi))*abs(n.dot(wo)));
	return fresnel(wi.dot(wh), eta, kappa)*NDF(n.dot(wh), alpha)*G_smith(n, wi, wo, wh, alpha)*(1/den);
}

//calcula la naturaleza de la luz reflejada para modelos con reflexion y transmision
/*
etai corresponderia al medio
etat corresponderia al material
*/
inline double fresnelDie(double etai, double etat, double cosinet, double cosinei){
    double parallel, perpendicular;
    parallel = ((etat*cosinei-etai*cosinet)/(etat*cosinei+etai*cosinet))*((etat*cosinei-etai*cosinet)/(etat*cosinei+etai*cosinet));
    perpendicular = ((etai*cosinei-etat*cosinet)/(etai*cosinei+etat*cosinet))*((etai*cosinei-etat*cosinet)/(etai*cosinei+etat*cosinet));
    return 0.5*(parallel+perpendicular);
}



//muestreo dielectrico reflexion
inline Vector reflexDielectric(Vector wi, Vector n){ //wi debe entrar como saliente y no incidente wi = r.d*-1
    //direccion reflejada wr
    return wi*-1+n*(n.dot(wi))*2;
}

//muestreo dielectrico transmision
inline Vector refraxDielectric(double etai, double etat, Vector wi, Vector n){ //wo sustituye a wi dado que la direccion incidente es la de observacion
    //falta convertir a direcciones locales y obtener la wt global de cualquier form
	//direccion refractada wt
	Vector wilocal = wi;
	Vector wtglobal;
	coordinateTraspose(n, wilocal);
	//std::cout<<acos((wilocal).dot(Vector(0,0,1)))<<std::endl;
    double ratio = etat/etai*-1;
    double cosinei = (wi).dot(n);
	double invratio = etai/etat;
    double cosinet = sqrt(1-invratio*invratio*(1-cosinei*cosinei))-1; //cosine de la nueva direccion
    Vector wtlocal = Vector(wilocal.x*ratio, wilocal.y*ratio, cosinet);
	Vector s, t;

	coordinateSystem(n, s, t);
	//Vector global
	wtglobal = s*wtlocal.x+t*wtlocal.y+n*wtlocal.z;
	return wtglobal;
}



#endif //MICROFACETUTILITIES_H
