// rt: un lanzador de rayos minimalista
// g++ -O3 -fopenmp rt.cpp -o rt
#include <chrono>
#include <cmath> //nota: clang si es compatible con cmath (no es necesario)
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <sys/random.h>

//include de clases

#include "Sphere.h"

//includes de funciones inline

#include "mathUtilities.h"
#include "pathTracingUtilities.h"
#include "samplingFunctions.h"
#include "misSamplingFunctions.h"
#include "shadeMethods.h"
#include "volumetricBasicFunctions.h"
#include "vptShadeMethods.h"
#include "vptSamplingFunctions.h"

//prototipos
double transmitance(Point x1, Point x2, double sigma_t);
double multipleT(Point x1, Point x2 , double sigma_t);

// PROYECTO 1
// calcular la intersección del rayo r con todas las esferas
// regresar true si hubo una intersección, falso de otro modo
// almacenar en t la distancia sobre el rayo en que sucede la interseccion
// almacenar en id el indice de spheres[] de la esfera cuya interseccion es mas cercana

//Calcula la escena usando muestreo de angulo sólido
inline Color marchingSA(Sphere &source, Point x, int indice, const Sphere &obj, Vector n, Vector wray, Vector &aux, double &costhetaMax, double alpha){
	Color L; 
	Vector cx = source.p-x;
	double normcx = sqrt(cx.dot(cx));
	cx = cx*(1/normcx);
	double costheta_max = sqrt(1-(source.r/normcx)*(source.r/normcx));
	costhetaMax = costheta_max;
	L = solidAngle(n, cx, wray, costheta_max, x, indice, obj, aux, alpha);
	return L;
}

/*Funciones legacy */

//path tracing explicito recursivo con MIS y ruleta rusa constante
inline Color explicitPathRecursive(const Ray &r, int bounce){
	double t;
	int id = 0;

	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() == negro
  

	if(spheres[id].radiance.x > 0)
		return Color();

	Point x = r.o + r.d*t;

	Vector n = (x-spheres[id].p);  
	n.normalize();
	
	Vector wi, wh;
	Vector wo = r.d * -1;
	Color  fs, Ld, Lind;
	double cosine, prob;
	
	//calcular Ld for each light 
	Ld = MIS(spheres[id], x, n, r.d, 0.001); //multiple importance sampling para iluminacion directa

	//ruleta rusa constante
	
	double q = 0.1;
	double continueprob = 1.0 - q;
	if(erand48(seed) < q)
		return Ld; // no continuar el camino (regresa 0 en la integral)
	//muestreo de BDSF para obtener Lind
	if(spheres[id].material==0){
		wi = cosineHemispheric(n);
		fs = spheres[id].c*(1/M_PI);
		prob = hemiCosineProb(n.dot(wi));
	}
	else{
		double alpha = 0.001;
		Vector wh = vectorFacet(alpha); //local
		Vector s, t;
		coordinateSystem(n, s, t);
		wh = s*wh.x+t*wh.y+n*wh.z; //global
		wi = wo*(-1)+wh*2*(wh.dot(wo));
		fs = frMicroFacet(spheres[id].eta, spheres[id].kappa, wi, wh, wo, alpha, n);
		prob = microFacetProb(wo, wh, alpha, n);
	}

	
	Ray recursiveRay = Ray(x, wi);
	cosine = n.dot(wi);
	bounce++;
	Lind = fs.mult(explicitPathRecursive(recursiveRay, bounce))*abs(cosine)*(1/(prob*continueprob));

	return (Ld + Lind);	
	//return explicitPathRecursive(recursiveRay, bounce, Ld, fs, cosine, prob)
}

inline Color explicitPath(const Ray &r){
	double t;
	int id = 0;

	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() == negro
  
	const Sphere &obj = spheres[id];
	

	if(obj.radiance.x > 0)
		return obj.radiance;
	else return explicitPathRecursive(r, 0);
}

//recursividad de cola, se puede convertir a iteracion facilmente
inline Color tailExplicitPath(const Ray &r, int bounces, Color Accum, Color fs, double factor){
	double t;
	int id = 0;

	if (!intersect(r, t, id))
		return Accum;	// el rayo no intersecto objeto, return Vector() == negro
  

	if(spheres[id].radiance.x > 0)
		return Accum;

	Point x = r.o + r.d*t;

	Vector n = (x-spheres[id].p);  
	n.normalize();
	
	Vector wi, wh;
	Vector wo = r.d * -1;
	Color  fs1, Ld, Lind;
	double cosine, prob;
	
	//calcular Ld for each light 
	Ld = MIS(spheres[id], x, n, r.d, 0.001); //multiple importance sampling para iluminacion directa

	//ruleta rusa constante
	
	double q = 0.1;
	double continueprob = 1.0 - q;
	if(erand48(seed) < q)
		return Accum+fs.mult(Ld)*factor; // no continuar el camino (regresa 0 en la integral)
	//muestreo de BDSF para obtener Lind
	fs1 = BDSF(wi, r.d, n, prob, id);
	Ray recursiveRay = Ray(x, wi);
	cosine = n.dot(wi);
	//se añade el color directo
	Accum = Accum + fs.mult(Ld)*factor;
	fs = fs.mult(fs1);
	return tailExplicitPath(recursiveRay, bounces++, Accum, fs, factor*abs(cosine)*(1/(prob*continueprob)));
}


inline Color implicitPath(const Ray &r, int bounces){
	double t;
	int id = 0;

	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() == negro
  
	const Sphere &obj = spheres[id];
	

	Point x = r.o + r.d*t;

	// determinar la dirección normal en el punto de interseccion
	Vector n = (x-obj.p);  // calcular x con el centro de la esfera como marco de referencia (radio)
	n.normalize();
	
	Vector wi, wh;
	Vector wo = r.d * -1;
	Color value, fs;
	double prob, cosine;

	//value = Le(x) light emitter con los que intersecta x 
	if(obj.radiance.x>0) 
		value =  obj.radiance; 
	
	if(bounces>5)
		return value;
	//sample BDSF
	if(obj.material==0){
		wi = cosineHemispheric(n);
		fs = obj.c*(1/M_PI);
		prob = hemiCosineProb(n.dot(wi));
	}
	else if(obj.material == 2){
		double etat=1.5;
		double etai=1.0;
		Vector wt = refraxDielectric(etai, etat, wo, n);
    	wt.normalize();
    	double F = fresnelDie(etai, etat, n.dot(wt), n.dot(wo));
		if(erand48(seed) < F)
		{
			wi = reflexDielectric(wo, n);
			wi.normalize();
			fs = (1/abs(n.dot(wi)));
			prob = 1.0;
		}	
		else{
			wi = wt;
			double ratio = etat/etai;
			fs = (1/abs(n.dot(wi)))*ratio*ratio;
			prob = 1.0;
		}
	
	}
	else if(obj.material == 1){
		double alpha = 0.3;
		Vector wh = vectorFacet(alpha); //local
		Vector s, t;
		coordinateSystem(n, s, t);
		wh = s*wh.x+t*wh.y+n*wh.z; //global
		wi = wo*(-1)+wh*2*(wh.dot(wo));
		fs = frMicroFacet(obj.eta, obj.kappa, wi, wh, wo, alpha, n);
		prob = microFacetProb(wo, wh, alpha, n);
	}
	wi.normalize();
	cosine = n.dot(wi);
	double q = 0.1;
	double continueprob = 1.0 - q;
	if(erand48(seed) < q)
		return value; // no continuar el camino (regresa 0 en la integral)

	Ray newray = Ray(x, wi);
	bounces++;
	value = (fs*abs(cosine)).mult(implicitPath(newray, bounces))*(1/(prob*continueprob)) + value;
	return value;

}


// Calcula el valor de color para el rayo dado
Color shade(const Ray &r) {
	
	double t;
	int id = 0;
	// determinar que esfera (id) y a que distancia (t) el rayo intersecta
	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() == negro
  
	const Sphere &obj = spheres[id];
	
	// PROYECTO 1
	// determinar coordenadas del punto de interseccion
	const Point x = r.o + r.d*t;

	// determinar la dirección normal en el punto de interseccion
	Vector n = (x-obj.p);  // calcular x con el centro de la esfera como marco de referencia (radio)
	n.normalize();

	if(obj.radiance.x>0) 
		return obj.radiance; //si un rayo impacta una fuente de luz se devuelve la radiancia directamente
	
	//muestreo de fuentes de luz
	
	int size = spheres.size();  // Devuelve el número de elementos en el vector
	
	Color L_total, L = Color();

	for(int i=0; i<size; i++){
		if(spheres[i].r==0){
			L = pLight(obj,x ,n, r.d, spheres[i].radiance, spheres[i].p, 0.0003); //luces puntuales
		}	
		
	}
	L = MIS(spheres[id], x, n, r.d, 0.003)+L;	
	L_total = L_total+L;
	/*
	Muestreo no basado en fuentes de luz (resuelve la escena completa en una solo muestreo, no fuentes puntuales)
	*/	
	//Color L_uniform = uniform(64, n, x, obj.c*(1/M_PI));

	//luces puntuales
	/*
	Color pcolor = Color();
	pcolor = plight(obj, x, n, spheres[7].radiance, spheres[7].p); 
	*/

	return L_total;
}

/*Funciones legacy terminan */

//volumetric path tracer
Color volumetricPathTracer(const Ray &r, double sigma_a, double sigma_s, int profundidad){
	Color Li = Color();
	Color Lo = Color();
	Color Ls = Color();
	double pdf_binary = 0;	
	//obtener la distancia de interseccion
	double t;
	int id = 0;
	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() ==
	Point x = r.o + r.d*t;

	//Transmitancia entre x y r.o
	double T = exp(-(sigma_a+sigma_s)*t); 
	//si intersecta una fuente de luz Lo = radiancia
	if(spheres[id].radiance.x>0) {
		Lo = spheres[id].radiance*T;
		return Lo;

	}
	//si la profundidad llega a 10 se regresa 0
	if(profundidad>=5) return Color();
	//muestrea la distancia de vuelo libre
	double d = freeFlightSample(sigma_a+sigma_s);
	
	//calcular el punto de interseccion xti
	Point x_new = r.o + r.d*d;
	//calcular la transmitancia entre x y xti
	//double T2 = transmitance(x, x_new, sigma_a+sigma_s);
	//calcular la funcion de fase
	double phase = isotropicPhaseFunction();
	//para calcular ls se debe hacer una llamada recursiva a volumetricPathTracer para obtener Li
	//prob = probabilidad del muestreo de fase 
	//obtener la nueva direccion del rayo con muestreo isotropico
	Vector wi_new = isotropicPhaseSample();
	//escoger entre pdf success o pdf failure
	if(d>=t){
		//pdf failure
		pdf_binary = pdfFailure(sigma_a+sigma_s, t);
		return Color();
		
	}
	else{
		//pdf success
		pdf_binary = pdfSuccess(sigma_a+sigma_s, t);
		Li = volumetricPathTracer(Ray(x_new, wi_new), sigma_a, sigma_s, profundidad+1);
	}

	
	double prob = isotropicPhaseProb();
	//Ls = Li*phase*(1/prob); //montecarlo de una sola muestra
	//despejando phase/phase
	Ls = Li;
	//print Li para debug si Li > 0
	
	double sigma_t = sigma_a+sigma_s;

	//montecarlo 
	//Color montecarlo = (Ls*T*sigma_s*(1/freeFlightProb(sigma_a+sigma_s, d)))*(1/pdf_binary);
	Color montecarlo = Ls*(sigma_s/sigma_t); //simplificacion despejando T
	return montecarlo*(1/pdf_binary);	
}

//volumetric path tracer explicito
Color volumetricPathTracerExplicit(const Ray &r, const double sigma_a, const double sigma_s, int profundidad, int idsource){
	Color Li = Color();
	Color Lo = Color();
	Color Ls = Color();
	Color Ld = Color();
	Point x_new = Point();
	Vector wi_new = Vector();
	double continueprob = 0.9;

	//obtener la distancia de interseccion
	double t;
	int id = 0;
	if (!intersect(r, t, id))
		return {};	// el rayo no intersecto objeto, return Vector() ==
	Point x = r.o + r.d*t;
	//Muestreo de fuentes de luz
	//no se debe regresar la radiancia de la fuente de luz dado que ya se muestreo (se puede hacer una excepcion para el primer impacto)
	if(spheres[id].radiance.x>0) {
		if(profundidad==0) return spheres[id].radiance*transmitance(r.o, x, sigma_a+sigma_s);
		else return Color();
	}

	//control de profundidad
	if(profundidad>=5) return Color();
	
	//ruleta rusa
	//if(erand48(seed) < 0.1) return Color();
	
	//muestrea la distancia de vuelo libre, aqui continua como si fuera volumetricPathTracer 
	double d = freeFlightSample(sigma_a+sigma_s);

	double pdf = 0;	
	//escoger entre pdf success o pdf failure
	if(d>=t){
		d = t;
		//pdf failure
		pdf = pdfFailure(sigma_a+sigma_s, t);
		return {};
		
	}
	else{
		x_new = r.o + r.d*d;
		wi_new = isotropicPhaseSample();
		//pdf success
		pdf = 1-pdfFailure(sigma_a+sigma_s, t);
	}
	
	//calcular la transmitancia entre x y xti
	double T = transmitance(r.o, x_new, sigma_a+sigma_s);
	//calcular la funcion de fase
	double phase = isotropicPhaseFunction();
	//para calcular ls se debe hacer una llamada recursiva a volumetricPathTracer para obtener Li
	//prob = probabilidad del muestreo de fase
	//obtener la nueva direccion del rayo con muestreo isotropico
	
	//antes de continuar, muestrear la iluminacion directa en x_new
	//calcular el color con angulo solido
	Vector wc = spheres[idsource].p-x_new;
	double normcx = sqrt(wc.dot(wc)); //magnitud de wc
	wc = wc*(1/normcx); //normalizar wc
	double costheta_max = sqrt(1-(spheres[idsource].r/normcx)*(spheres[idsource].r/normcx));
	//vector wi sampleado con angulo solido
	Vector wi_ld = solidAngle(wc, costheta_max);
	//rayo con origen en x_new y direccion wi
	Ray shadowRay = Ray(x_new, wi_ld);
	//verificar la visibilidad con la fuente de luz
	//si el shadow ray intersecta en el mismo id que la fuente de luz, entonces hay visibilidad
	//en caso contrario no hay visibilidad
	double t2;
	int id2 = 0;
	intersect(shadowRay, t2, id2);	
	//si hay visibilidad, la contribucion es la radiancia de la fuente de luz
	if(id2==idsource){
		//obtener la radiancia de la fuente de luz
		Color Le = spheres[idsource].radiance;
		//calcular Ls
		//montecarlo de una sola muestra
		Color Ls = Le*(phase*transmitance(x_new, spheres[idsource].p, sigma_a+sigma_s))*sigma_s;
		//probabilidad de muestreo
		double prob = solidAngleProb(costheta_max);
		//calcular Ls aproximado
		Color Ls_aprox = Ls*T*(1/prob);
		//acumular Ls aproximado
		Ld = Ls_aprox;
	}
	//en caso de no haber visibilidad
	else{
		Color Ls_aprox = Color();
		Ld = Ls_aprox;
	}

	

	Li = volumetricPathTracerExplicit(Ray(x_new, wi_new), sigma_a, sigma_s, profundidad+1, idsource);


	
	double prob = isotropicPhaseProb();
	//Ls = Li*phase*(1/prob); //montecarlo de una sola muestra
	//despejando phase/phase
	Ls = Li;
	
	double sigma_t = sigma_a+sigma_s;
	//montecarlo
	Color montecarlo = Ls*(sigma_s/sigma_t); //simplificacion despejando T, solo es valido con muestreo free flight
	return (Ld*(1/freeFlightProb(sigma_t,d)) + montecarlo)*(1/pdf);

}


//volumetric path tracer explicito
Color volumetricPathTracerExplicitEquiAngular(const Ray &r, double sigma_a, double sigma_s, int profundidad, int idsource){
	Color Li = Color();
	Color Lo = Color();
	Color Ls = Color();
	Color Ld = Color();
	double continueprob = 0.9;

	//obtener la distancia de interseccion
	double t;
	int id = 0;
	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() ==
	Point x = r.o + r.d*t;
	//Muestreo de fuentes de luz
	//no se debe regresar la radiancia de la fuente de luz dado que ya se muestreo (se puede hacer una excepcion para el primer impacto)
	if(spheres[id].radiance.x>0) {
		if(profundidad==0 )return spheres[id].radiance*transmitance(r.o, x, sigma_a+sigma_s);
		else return Color();
	}
	
	//ruleta rusa
	if(erand48(seed) < 0.1) return Color();

	//muestrea la distancia de vuelo libre, aqui continua como si fuera volumetricPathTracer 
	//double d = freeFlightSample(sigma_a+sigma_s);
	
	//calculos necesarios para equi-angular sampling
	double thetaA = 0;
	double thetaB = 0;
	Point c = spheres[idsource].p;
	//calcular la proyeccion ortogonal de c en el rayo, punto mas cercano a c en el rayo, x0
	Point x0 = r.o + r.d*((c-r.o).dot(r.d)/(r.d.dot(r.d)));
	//verificar si el punto x0 esta entre r.o y x
	if((x0-r.o).dot(r.d)<0) x0=r.o;
	if((x0-x).dot(r.d)>0){
		x0 = x;
	}

	//calcular la magnitud de x0-c
	double D = sqrt((x0-c).dot(x0-c)); 
	//theta A y theta B son los angulos de apertura de la fuente de luz, intervalo de integracion, en este caso usaremos r.o y x como puntos de integracion
	//calcular los lados del triangulo rectangulo, comparten D
	double A = sqrt((x0-r.o).dot(x0-r.o))*-1;
	double B = sqrt((x-x0).dot(x-x0));
	//calcular el angulo theta A y theta B

	//if(((t+A)-B)>0.01) std::cout<<t+A-B<<std::endl;

	thetaA = atan2(A,D);
	thetaB = atan2(B,D);
	//muestrear la distancia de con equi-angular sampling
	double d = equiAngularSample(D, thetaA, thetaB);



	//print d para debug
	//std::cout<<A<<"ww"<<d<<"ww"<<B<<"ww"<<t<<std::endl;

	//calcular el punto de interseccion xti,esto en muestreo de free flight
	//Point x_new = r.o + r.d*d;

	Point x_new = x0 + r.d*d;

	
	
	//calcular la transmitancia entre x y xti
	double T = transmitance(r.o, x_new, sigma_a+sigma_s);
	//calcular la funcion de fase
	double phase = isotropicPhaseFunction();
	//para calcular ls se debe hacer una llamada recursiva a volumetricPathTracer para obtener Li
	//prob = probabilidad del muestreo de fase
	//obtener la nueva direccion del rayo con muestreo isotropico
	Vector wi_new = isotropicPhaseSample();
	//antes de continuar, muestrear la iluminacion directa en x_new
	Vector n = (x_new-spheres[id].p);
	n.normalize();
	//calcular el color con angulo solido
	Vector wc = spheres[idsource].p-x_new;
	double normcx = sqrt(wc.dot(wc));
	wc = wc*(1/normcx);
	double costheta_max = sqrt(1-(spheres[idsource].r/normcx)*(spheres[idsource].r/normcx));
	//vector wi sampleado con angulo solido
	Vector wi_ld = solidAngle(wc, costheta_max);
	//rayo con origen en x_new y direccion wi
	Ray shadowRay = Ray(x_new, wi_ld);
	//verificar la visibilidad con la fuente de luz
	//si el shadow ray intersecta en el mismo id que la fuente de luz, entonces hay visibilidad
	//en caso contrario no hay visibilidad
	double t2;
	int id2 = 0;
	intersect(shadowRay, t2, id2);
	//si hay visibilidad, la contribucion es la radiancia de la fuente de luz
	if(id2==idsource){
		//obtener la radiancia de la fuente de luz
		Color Le = spheres[idsource].radiance;
		//calcular Ls
		//montecarlo de una sola muestra
		Color Ls = Le*(phase*transmitance(x_new, spheres[idsource].p, sigma_a+sigma_s))*sigma_s;
		//probabilidad de muestreo
		double prob = solidAngleProb(costheta_max);
		//calcular Ls aproximado
		Color Ls_aprox = Ls*T*(1/prob);
		//acumular Ls aproximado
		Ld = Ls_aprox;
	}
	//en caso de no haber visibilidad
	else{
		Color Ls_aprox = Color();
		Ld = Ls_aprox;
	}

	
	
	Li = volumetricPathTracerExplicitEquiAngular(Ray(x_new, wi_new), sigma_a, sigma_s, profundidad+1, idsource);
	
	double prob = isotropicPhaseProb();
	//Ls = Li*phase*(1/prob); //montecarlo de una sola muestra
	//despejando phase/phase
	Ls = Li;
	
	double sigma_t = sigma_a+sigma_s;
	//montecarlo

	//montercarlo con muestreo equi-angular
	Color montecarlo = (Ls*T*sigma_s);
	double t_ajustada = 1;
	return (Ld + montecarlo)*(1/equiAngularProb(D, thetaA, thetaB, d))*(1/continueprob);

}

//volumetric path tracer explicito con equi-angular sampling luz puntual
Color volumetricPathTracerExplicit2(const Ray &r, double sigma_a, double sigma_s, int profundidad, int idsource){
	Color Li = Color();
	Color Lo = Color();
	Color Ls = Color();
	Color Ld = Color();
	Point x_new = Point();
	Vector wi_new = Vector();
	double continueprob = 0.9;
	bool continuar = false;
	double pdf = 0;
	double sigma_t = sigma_a+sigma_s;

	//obtener la distancia de interseccion
	double t;
	int id = 0;
	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() ==
	Point x = r.o + r.d*t;


	double q = 1 - continueprob;
	//muestrea la distancia de vuelo libre, aqui continua como si fuera volumetricPathTracer 
	double d = freeFlightSample(sigma_a+sigma_s);
/*	

	//calculos necesarios para equi-angular sampling
	double thetaA = 0;
	double thetaB = 0;
	Point c = spheres[idsource].p;
	//calcular la proyeccion ortogonal de c en el rayo, punto mas cercano a c en el rayo, x0
	Point x0 = r.o + r.d*((c-r.o).dot(r.d)/(r.d.dot(r.d)));
	//verificar si el punto x0 esta entre r.o y x
	if((x0-r.o).dot(r.d)<0) {
		x0 = r.o;
	}
	if((x0-x).dot(r.d)>0){
		x0 = x;
	} 
	


	//calcular la magnitud de x0-c
	double D = sqrt((x0-c).dot(x0-c)); 
	//theta A y theta B son los angulos de apertura de la fuente de luz, intervalo de integracion, en este caso usaremos r.o y x como puntos de integracion
	//calcular los lados del triangulo rectangulo, comparten D
	double A = sqrt((x0-r.o).dot(x0-r.o))*-1;
	double B = sqrt((x-x0).dot(x-x0));
	//calcular el angulo theta A y theta B

	//if(((t+A)-B)>0.01) std::cout<<t+A-B<<std::endl;

	thetaA = atan2(A,D);
	thetaB = atan2(B,D);
	//muestrear la distancia de con equi-angular sampling
	double d = equiAngularSample(D, thetaA, thetaB);
*/

	//print d para debug
	//std::cout<<A<<"ww"<<d<<"ww"<<B<<"ww"<<t<<std::endl;


	
	//escoger entre pdf success o pdf failure
	if(d>=t){
		
		//pdf failure
		pdf = pdfFailure(sigma_t, t);
		d = t;
		return Color();
		
	}
	else{
		//pdf success
		//pdf = sigma_t/(exp(sigma_t*d)-exp(sigma_t*(d-t)));
		pdf = freeFlightProb(sigma_t, d);
		x_new = r.o + r.d*d;
		wi_new = isotropicPhaseSample();
		
		
	}


	

	//x_new = x0 + r.d*d; //para equi-angular sampling

	
	//calcular la transmitancia entre el origen del rayo y el nuevo punto
	double T = transmitance(r.o, x_new, sigma_t);
	//calcular la funcion de fase
	double phase = isotropicPhaseFunction();
	

	
	//verificar visibilidad con luz puntual
	//si hay visibilidad, la contribucion es la radiancia de la fuente de luz
	Point light = spheres[idsource].p;

	if(visibility(light, x_new)){
		//obtener la radiancia de la fuente de luz
		Color Le = spheres[idsource].radiance;
		//calcular ls distancia entre la fuente de luz y x_new
		double distanceLight = (light-x_new).dot(light-x_new);
		//dividir la radiancia por la distancia al cuadrado
		Le = Le*(1/distanceLight);	
		//montecarlo de una sola muestra
		Ls = Le*phase*transmitance(x_new, light, sigma_t);//la probabilidad de muestreo es 1;
		Ld = Ls*T*sigma_s;
	}
	//en caso de no haber visibilidad
	else{
		Ld = Color();
	}

	//if(erand48(seed)<q) return Ld*(1/equiAngularProb(D, thetaA, thetaB, d))*(1/q);
	if(erand48(seed)<q) return Ld*(1/pdf)*(1/q);
	//recursion equi-angular sampling
	//wi_new = isotropicPhaseSample();//obtener la nueva direccion del rayo con muestreo isotropico
	
	
	Li = volumetricPathTracerExplicit2(Ray(x_new, wi_new), sigma_a, sigma_s, profundidad+1, idsource);
	


	//double prob = isotropicPhaseProb();
	//Ls = Li*phase*(1/prob); //multiple scattering
	//despejando phase/phase
	Ls = Li;
	
	//montecarlo
	//Color montecarlo = Ls*(sigma_s/sigma_t); //simplificacion despejando T, solo es valido con muestreo free flight
	//montercarlo con muestreo equi-angular
	Color montecarlo = (Ls*T*sigma_s);
	//return ((Ld + montecarlo)*(1/equiAngularProb(D, thetaA, thetaB, d)))*(1/continueprob);//equiangular sampling 

	return (Ld*(1/pdf) + montecarlo*(1/pdf))*(1/continueprob);

}



int main(int argc, char *argv[]) {
	// semilla para generador de números aleatorios
	while(getentropy(seed,3));

	std::chrono::time_point<std::chrono::system_clock> start, end;

	start = std::chrono::system_clock::now();

	int w = 1024, h = 768; // image resolution
  
	// fija la posicion de la camara y la dirección en que mira
	Ray camera( Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize() );

	// parametros de la camara
	Vector cx = Vector( w * 0.5095 / h, 0., 0.); 
	Vector cy = (cx % camera.d).normalize() * 0.5095;
  
	// auxiliar para valor de pixel y matriz para almacenar la imagen
	Color *pixelColors = new Color[w * h];


	// PROYECTO 1
	// usar openmp para paralelizar el ciclo: cada hilo computara un renglon (ciclo interior),
	#pragma omp parallel for schedule(dynamic, 1)
	for(int y = 0; y < h; y++) 
	{ 	
		// recorre todos los pixeles de la imagen
		fprintf(stderr,"\r%5.2f%%",100.*y/(h-1));
		for(int x = 0; x < w; x++ ) {
			int idx = (h - y - 1) * w + x; // index en 1D para una imagen 2D x,y son invertidos
			Color pixelValue = Color(); // pixelValue en negro por ahora
			
			// para el pixel actual, computar la dirección que un rayo debe tener
			//Vector cameraRayDir = cx * (double(x)/w - .5) + cy * (double(y)/h - .5) + camera.d;
			
			// computar el color del pixel para el punto que intersect*o el rayo desde la camara
			
			//pixelValue = pathTracer( Ray(camera.o, cameraRayDir.normalize()));
		
			//montercarlo path tracing con n rpp (rays per pixel)
			int rpp = atoi(argv[1]); //en esta implementacion tenemos 1 muestra por rayo desde la camara 1 rpp = 1 spp

			for(int i=0; i<rpp; i++){
				Vector cameraRayDir = cx * ((double(x)+erand48(seed)-0.5)/w - .5) + cy * ((double(y)+erand48(seed)-0.5)/h - .5) + camera.d;
			
				// computar el color del pixel para el punto que intersectó el rayo desde la camara

				//pixelValue = explicitPathRecursive2(Ray(camera.o, cameraRayDir.normalize()),0)+ pixelValue;

				//pixelValue = iterativePathTracer(Ray(camera.o, cameraRayDir.normalize()))+ pixelValue;

				//pixelValue = tailExplicitPath(Ray(camera.o, cameraRayDir.normalize()), 0, Color(0,0,0),Vector(1,1,1), 1)+pixelValue;
				//pixelValue = rayMarching2(Ray(camera.o, cameraRayDir.normalize()), 0.001,0.009, 0.1, 7)+pixelValue;

				//pixelValue = rayMarchingGlobal(Ray(camera.o, cameraRayDir.normalize()),0.001,0.0125,10)+pixelValue;
				//pixelValue = volumetricPathTracer(Ray(camera.o, cameraRayDir.normalize()), 0.001, 0.009, 0)+pixelValue;
				//pixelValue = volumetricPathTracerExplicit2(Ray(camera.o, cameraRayDir.normalize()), 0.001, 0.009, 0,5)+pixelValue;
				//pixelValue = volumetricPathTracerExplicitEquiAngular(Ray(camera.o, cameraRayDir.normalize()), 0.001, 0.009, 0,7)+pixelValue;

				//pixelValue = volumetricPathTracer3(Ray(camera.o, cameraRayDir.normalize()), 0.001, 0.009, 0)+pixelValue;
				//pixelValue = volumetricPathTracerRecursive(Ray(camera.o, cameraRayDir.normalize()),0.01,0.009)+pixelValue;
				pixelValue = implicitVPTracerRecursive(Ray(camera.o, cameraRayDir.normalize()),0.01,0.009)+pixelValue;
			}

			pixelValue = pixelValue * (1/static_cast<double>(rpp)); //promedio de color de cada pixel

			// limitar los tres valores de color del pixel a [0,1]
			pixelColors[idx] = Color(clamp(pixelValue.x), clamp(pixelValue.y), clamp(pixelValue.z));
		}
	}
	

	fprintf(stderr,"\n");

	// PROYECTO 1
	// Investigar formato ppm
	FILE *f = fopen("image.ppm", "w");
	// escribe cabecera del archivo ppm, ancho, alto y valor maximo de color
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
	for (int p = 0; p < w * h; p++) 
	{ // escribe todos los valores de los pixeles
    		fprintf(f,"%d %d %d ", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y), 
				toDisplayValue(pixelColors[p].z));
  	}
  	fclose(f);

  	delete[] pixelColors;

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout<< "elapsed time: " << elapsed_seconds.count() << "s\n";

	return 0;
}
