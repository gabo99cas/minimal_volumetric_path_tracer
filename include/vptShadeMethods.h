//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef VPTSHADEMETHODS_H
#define VPTSHADEMETHODS_H

#include "samplingFunctions.h"
#include "volumetricBasicFunctions.h"
#include "vptSamplingFunctions.h"
#include "rayMarchingMethods.h"

//Las dos funciones a continuación no son como tal parte de shadow methods, fue necesario agregarlas aquí por problemas al compilar

//calcula la fs y la nueva direccion para un camino en Path tracing. Lind
inline Color bdsf(Vector &aux, Vector wray, Vector n,  double &prob, int id){
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
			fs1=Color(1,1,1)*(1/n.dot(wi))*F;

			prob = F;
		}
		else{
			wi = wt;
			fs1=Color(1,1,1)*(1/n.dot(wi))*(1-F)*1.5*1.5;

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

//calcula la contribucion de una luz puntual
inline Color pLight(const Sphere &obj, Point x, Vector n, Vector wray, Color I, Point light, double alpha){
    Color Le;
    Color fr;
    if(visibility(light, x))
    {
        Le = I*(1/((light-x).dot(light-x))); //devolver radiancia
    }
    else if(visibilityVPT(light, x))
    {
        Le = I*(1/((light-x).dot(light-x))); //devolver radiancia
        Le = Le*multipleT(x, light,0.05+0.009);
    }
    else Le = Color(); //no hay visibilidad
    //modelo microfacet
    Vector wi = light-x;
    wi.normalize();
    Vector wo = wray*-1;
    coordinateTraspose(n, wo);
    coordinateTraspose(n, wi);
    wi.normalize();
    wo.normalize();
    Vector wh = wi+wo;
    wh.normalize();
    if(obj.material==1){
        fr=frMicroFacet(obj.eta, obj.kappa, wi, wh, wo, alpha, Vector(0,0,1));
    }
    else fr = obj.c*(1/M_PI);
    Color L = Le.mult(fr)*n.dot((light-x).normalize());
    return L;
}

//terminan funciones externas

//empiezan funciones LEGACY

//volumetric path tracer explicito VPT para multiples fuentes de luz
inline Color volumetricPathTracer3(const Ray &r, double sigma_a, double sigma_s, int profundidad){
	Color Li = Color();
	Color Lo = Color();
	Color Lo2 = Color();
	Color Ls = Color();
	Color Ld = Color();
	Color fs1 = Color();
	Point x_new = Point();
	Vector wi_new = Vector();
	double continueprob = 0.9;
	double q = 1 - continueprob;
	double pdf = 0;
	double sigma_t = sigma_a+sigma_s;

	double phase = isotropicPhaseFunction();
	//obtener la distancia de interseccion
	double t;
	int id = 0;
	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() ==
	Point xt = r.o + r.d*t;

	//obtener la contribucion de xt para sumarla
	//esto consiste en calcular la iluminacion directa en xt (sin olvidar la emision)
	//finalmente se multiplica por la transmitancia Tr*Lo

	//normal en xt
	Vector normalXT = xt-spheres[id].p;
	normalXT.normalize();

	//MIS resuelve la iluminacion directa para un punto en particular, ya toma en cuenta la bdsf
	//Lo = MIS(spheres[id], xt, normalXT, r.d, spheres[id].alpha);
	//for each light source
	int number = spheres.size();  // Devuelve el número de elementos en el vector

	//equivalente a MIS pero para fuentes puntuales
	for(int light = 0; light<number; light++)
	{
		//resuelve las fuentes puntuales
		if(spheres[light].r==0) Lo = pLight(spheres[id], xt, normalXT, r.d, spheres[light].radiance, spheres[light].p, spheres[id].alpha) + Lo;

	}

	//calcular la TR para xt

	double Tr1 = transmitance(r.o, xt, sigma_t);

	//hasta aqui Lo solo hace handle de la iluminacion directa (cosa que ya se soluciono en el caso de d>=t)
	//se debe expandir para que continue la recursion en este punto tambien y Lo tenga completa su iluminacion

	//calcular la parte recursiva de Lo
	//bdsf
	Vector wi; //para recuperar la nueva direccion
	double samplingProbability; //para recuperar la probabilidad
	fs1 = bdsf(wi, r.d, normalXT, samplingProbability, id);

	wi.normalize();
	//nuevo rayo con la otra direccion
	Ray newRay = Ray(xt, wi);

	//coseno
	double cosine = normalXT.dot(wi);

	if (profundidad>5) return Lo*Tr1;
	if(erand48(seed)<q)
	{
		return Lo*Tr1*(1/q);
	}
	//llamada recursiva con un nuevo rayo
	Lo2 = volumetricPathTracer3(newRay, sigma_a, sigma_s, profundidad+1);

	//aplicarle la bdsf y los factores al color resultante del segundo paso
	Lo2 = fs1.mult(Lo2)*(1/samplingProbability)*cosine;


	//escoger una fuente de luz para conexion directa
	int idsource = 0;

	//arreglo para las fuentes (por el momemnto estatico)
	int arr[4] = {-1,-1,-1,-1};

	int n = spheres.size();  // Devuelve el número de elementos en el vector
	int j = 0;
	//guardar en las posiciones del arreglo las fuentes de luz
	for(int i = 0; i<n; i++){
		if(spheres[i].radiance.x>0 || spheres[i].radiance.y>0 || spheres[i].radiance.z>0){
			arr[j] = i;
			j++;
		}
	}

	//contar las fuentes

	int count = 0;
	for(int i = 0; i<4; i++){
		if(arr[i]!=-1) count++;
	}

	if(count==0) return Color();
	//para metodo 2 muestrear una fuente de luz

	//probabildad uniforme para cada fuente
	double prob = 1.0/count;

	//muestrear una fuente de luz
	idsource = arr[(int)(erand48(seed)*count)];

/* //en caso de hacer muestreo equi-angular
	double D = 0;
	double thetaA = 0;
	double thetaB = 0;
	Point x0 = Point();
	equiAngularParams(idsource, x, x0, r, D, thetaA, thetaB);
	double d = equiAngularSample(D, thetaA, thetaB);
*/

	//muestrea la distancia de vuelo libre, aqui continua como si fuera volumetricPathTracer
	double d = freeFlightSample(sigma_a+sigma_s);

	//escoger entre pdf success o pdf failure
	if(d>=t){
		//pdf failure
		//pdf = pdfFailure(sigma_t, t); esta parte sera descartada debido a que no se posee un pdf failure para equiangular
		d = t;
		//en caso de fallo ya no es posible evaluar el medio y la recursividad sucederá unicamente en Lo
		return (Lo+Lo2)*Tr1*(1/continueprob);

	}
	else{
		//pdf success
		//si la distancia es menor que la superficie es necesario que la nueva direccion sea muestreada por fase
		//pdf = sigma_t/(exp(sigma_t*d)-exp(sigma_t*(d-t)));
		pdf = freeFlightProb(sigma_t, d);
		x_new = r.o + r.d*d;
		wi_new = isotropicPhaseSample();
	}

/*
	//para equi-angular sampling
	x_new = x0 + r.d*d;
	pdf = equiAngularProb(D, thetaA, thetaB, d);
	wi_new = isotropicPhaseSample();
*/
	//calcular la transmitancia entre el origen del rayo y el nuevo punto
	double T = transmitance(r.o, x_new, sigma_t);

	//calcula la contribucion directa de la fuente de luz
	//determinar si es puntual o esferica
	if(spheres[idsource].r==0){
		//ejecutar rutina de luz puntual
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
		Ld = Ls*T*sigma_s*(1/prob); //1/prob por el evento de muestreo de la fuente
		}
		//en caso de no haber visibilidad
		else	Ld = Color();
	}
	else{
		//ejecutar rutina de angulo solido
	}

	if(erand48(seed)<q)
	{
		return Ld*(1/pdf)*(1/q);
	}

	//recursion
	Li = volumetricPathTracer3(Ray(x_new, wi_new), sigma_a, sigma_s, profundidad+1);

	//probabilidad de muestreo de fase se cancela con phase/phase
	Ls = Li;

	//montecarlo
	Color Lind = (Ls*T*sigma_s);
	Color montecarlo = Ld + Lind;
	return ((Lo+Lo2)*Tr1+Ld*(1/pdf))*(1/continueprob)+Ls;
}

//volumetric path tracer explicito VPT para multiples fuentes de luz, version iterativa
inline Color volumetricPathTracer3alt(const Ray &r, double sigma_a, double sigma_s, int profundidad){
	Color Li = Color();
	Color Lo = Color();
	Color Ls = Color();
	Color Ld = Color();
	Point x_new = Point();
	Vector wi_new = Vector();
	double continueprob = 0.5;
	double pdf = 0;
	double sigma_t = sigma_a+sigma_s;

	double phase = isotropicPhaseFunction();
	//obtener la distancia de interseccion
	double t;
	int id = 0;
	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() ==
	Point x = r.o + r.d*t;

	double q = 1 - continueprob;
	//escoger una fuente de luz para conexion directa
	int idsource = 0;

	//arreglo para las fuentes (por el momento estatico)
	int arr[4] = {-1,-1,-1,-1};

	int n = spheres.size();  // Devuelve el número de elementos en el vector
	int j = 0;
	//guardar en las posiciones del arreglo las fuentes de luz
	for(int i = 0; i<n; i++){
		if(spheres[i].radiance.x>0 || spheres[i].radiance.y>0 || spheres[i].radiance.z>0){
			arr[j] = i;
			j++;
		}
	}

	//contar las fuentes

	int count = 0;
	for(int i = 0; i<4; i++){
		if(arr[i]!=-1) count++;
	}

	if(count==0) return Color();
	//solo se puede hacer free flight sampling
	double d = freeFlightSample(sigma_t);

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

	//calcular la transmitancia entre el origen del rayo y el nuevo punto
	double T = transmitance(r.o, x_new, sigma_t);

	//para metodo 1 iterar sobre arr
	Color accum = Color();
	for(int i = 0; i<count; i++){

		//para cada indice de fuente de luz
		idsource = arr[i];

		//calcula la contribucion directa de la fuente de luz
		//determinar si es puntual o esferica
		if(spheres[idsource].r==0){
			//ejecutar rutina de luz puntual
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
			Ld = Ls*T*sigma_s; //1/prob por el evento de muestreo de la fuente
			}
			//en caso de no haber visibilidad
			else	Ld = Color();
		}
		else{
			//ejecutar rutina de angulo solido
			}
		accum = accum + Ld;

	}

	if(erand48(seed)<q) return accum*(1/pdf)*(1/q);
	//recursion
	Li = volumetricPathTracer3alt(Ray(x_new, wi_new), sigma_a, sigma_s, profundidad+1);

	//probabilidad de muestreo de fase se cancela con phase/phase
	Ls = Li;

	//montecarlo
	Color montecarlo = (Ls*T*sigma_s);
	return (accum*(1/pdf) + montecarlo*(1/pdf))*(1/continueprob);

}

//path tracing explicito recursivo con MIS y ruleta rusa constante, se incluyen objetos volumetricos
inline Color explicitPathRecursive2(const Ray &r, int bounce){
	double sigma_a = 0.05;
	double sigma_s = 0.009;
	double sigma_t = sigma_a+sigma_s;
	double t,t2;
	int id = 0;
	Color Ls = Color();


	Color  fs, Ld, Lind;

	int size = spheres.size();  // Devuelve el número de elementos en el vector

	if (!intersectV2(r, t,t2, id))
		return Color();	// el rayo no intersecto objeto, return Vector() == negro


	if(spheres[id].radiance.x > 0)
		return Color();

	Point x = r.o + r.d*t;

	Point xt;

	if(spheres[id].material==3){
		//haz ray marching hasta terminar la esfera
		int steps = 100;
		//para calcular el paso se debe conocer la distancia al limite de la esfera
		double distance = t2-t;
		double step = distance/steps;
		for (int i=0; i<steps; i++){
			xt = x + r.d*step*i;
			//for each light source
			for(int light = 0; light<size; light++){
				//resuelve las fuentes puntuales
				if(spheres[light].r==0){
					Ls = punctualVolumetric(light,xt,isotropicPhaseFunction(), sigma_t, sigma_s)*step*transmitance(x,xt,sigma_t)+Ls;
				}
			}


		}
		return Ls + explicitPathRecursive2(Ray(xt,r.d),bounce)*transmitance(x,xt, sigma_t);

	}



	Vector n = (x-spheres[id].p);
	n.normalize();


	Vector wi, wh;
	Vector wo = r.d * -1;

	double cosine, prob;

	//calcular Ld for each light
	//for each light source

	for(int light = 0; light<size; light++) {
			//resuelve las fuentes puntuales
			if(spheres[light].r==0) Ld = pLight(spheres[id], x, n, r.d, spheres[light].radiance, spheres[light].p, spheres[id].alpha) + Ld;
	}
	Ld = MIS(spheres[id], x, n, r.d, spheres[id].alpha) + Ld;

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
	Lind = fs.mult(explicitPathRecursive2(recursiveRay, bounce))*abs(cosine)*(1/(prob*continueprob));

	return (Ld + Lind);
	//return explicitPathRecursive(recursiveRay, bounce, Ld, fs, cosine, prob)
}

//implementacions en uso, las anteriores ya se consideran legacy y deben ser borradas proximamente

//version iterativa explicita
inline Color volumetricPathTracerIterative(const Ray &r, double sigma_a, double sigma_s) {
    Color Li, Lo, Ls, Ld, fsActual;
    std::stack<std::tuple<Ray, double, double, int, double, Color>> stack;
    stack.emplace(r, sigma_a, sigma_s, 0, 1.0, Color(1,1,1));
	int number = spheres.size();
	Lo = Color(0, 0, 0);


    while (!stack.empty()) {
        auto [ray, sigma_a, sigma_s, profundidad, factor, pathThroughput] = stack.top();
        stack.pop();


        double sigma_t = sigma_a + sigma_s;
        double continueprob = 0.6;
        double q = 1 - continueprob;

        double t;
        int id = 0;
        if (!intersect(ray, t, id)) {
            continue;
        }
    	//ruleta al inicio de la iteracion para decidir si el camino continua o termina
    	if (erand48(seed) < q) {
    		//Li = Ld + Li;
    		continue;
    	}

    	//determinar el limite del rayo previamente
    	Point xs = ray.o + ray.d * t;
    	Vector normalXS = xs - spheres[id].p;
    	normalXS.normalize();


        int arr[4] = {-1, -1, -1, -1};
        int j = 0, count = 0;
        for (int i = 0; i < number; i++) {
            if (spheres[i].radiance.x > 0 || spheres[i].radiance.y > 0 || spheres[i].radiance.z > 0) {
                arr[j++] = i;
            }
        }

        count = j;

        if (count == 0) {
            continue;
        }

    	//probabildad uniforme para cada fuente
    	double prob = 1.0/count;

        int idsource = arr[static_cast<int>(erand48(seed) * count)];

    	/* Metodo de muestreo */
		//escoger entre uno u otro

		//free flight sampling

/*
    	double d = freeFlightSample(sigma_t);
    	double pdf;
*/

    	//operaciones en caso de equiangular sampling
    	double D = 0;
    	double thetaA = 0;
    	double thetaB = 0;
    	Point x0 = Point();
    	equiAngularParams(idsource, xs, x0, ray,D, thetaA, thetaB);
    	double d = equiAngularSample(D, thetaA, thetaB);
    	double pdf;

    	double TrActual = transmitance(ray.o, xs, sigma_t); //tr de x a xs

    	//evento de fallo
		if (erand48(seed) < TrActual) { //evento para decidir fallo en equiangular
        //if (d>=t) {

        	//d = t;
        	pdf = TrActual;

			//esto solo aplica para el freeflight, equiangular no tiene posibilidad de caer aqui
        	//si falló, la dispersion no ocurrio, el rayo llego a la superficie directamente

        	for (int light = 0; light < number; light++) {
        		if (spheres[light].r == 0) {
        			double Trs = transmitance(xs, spheres[light].p, sigma_t); //es valido en tanto sea puntual
        			Lo = Lo + pLight(spheres[id], xs, normalXS, ray.d, spheres[light].radiance, spheres[light].p, spheres[id].alpha)*Trs;
        			//incluir el efecto de la transmitancia entre xs y p (single scattering)
        		}
        	}




        	Lo = pathThroughput.mult(Lo)*factor*(1/continueprob);//TrActual*(1/pdf) = 1; no tiene caso calcular la transmitancia debido a que al dividirla
        	//entre la probabilidad de fallo que es la misma transmitancia, queda en 1


        	Vector wi;
        	double samplingProbability;
        	fsActual = bdsf(wi, ray.d, normalXS, samplingProbability, id);
        	wi.normalize();

        	double cosine = normalXS.dot(wi);

        	//calcular el path troughput, ahora esto se agrega en el stack directamente para que este disponible en next iteracio
        	//factor = 1/samplingProbability*cosine*factor; pdf no se incluye porque cancela con la transmitancia
        	//agregamos en el stack los valores que corresponden a la siguiente iteracion
        	stack.emplace(Ray(xs, wi), sigma_a, sigma_s, profundidad + 1, factor*1/samplingProbability*cosine*(1/continueprob), pathThroughput.mult(fsActual));

        	Li = Lo + Li;
        }
    	else {
    		//omitir esta parte en equiangular
    		//pdf = freeFlightProb(sigma_t,d)*(1.0-TrActual); //pfallo = Tractual

    		//probabilidad con equiangular
    		pdf = equiAngularProb(D, thetaA, thetaB, d)*(1.0-TrActual);

    		//Point xt = ray.o + ray.d * d;
    		//para equiangular sampling xt tiene como referencia x0
    		Point xt = x0 + ray.d * d;

    		double T = transmitance(ray.o, xt, sigma_t);



    		//calculando un single scattering
    		if (spheres[idsource].r == 0) {
    			Point light = spheres[idsource].p;
    			if (visibility(light, xt)) {
    				Color Le = spheres[idsource].radiance;
    				double distanceLight = (light - xt).dot(light - xt);
    				Le = Le * (1 / distanceLight);
    				Ls =  Le * transmitance(xt, light, sigma_t)* isotropicPhaseFunction(); //leer Ls ray marchind
    				Ld = Ls * T * sigma_s * (1 / prob); //este es el calculo de la luz directa, hasta aqui
    			}

    		}

    		Vector wi_new = isotropicPhaseSample(); //nueva direccion para siguiente rayo
    		//pdf *= isotropicPhaseProb(); la fase cancela a la probabilidad

    		Ld = pathThroughput.mult(Ld) * (1 / pdf) * factor * (1/continueprob);

    		Li = Li + Ld;


    		// -legacy Ld = pathThroughput.mult(Ld) * (1 / pdf) * factor; //aplicar el path throughput
    		//es necesario multiplicar por la pdf inmediata para completar el montecarlo

    		stack.emplace(Ray(xt, wi_new), sigma_a, sigma_s, profundidad + 1, factor*T*sigma_s*(1/continueprob)*(1/pdf), pathThroughput);


    	}
    }

    return Li;
}

//version recursiva explicita
inline Color volumetricPathTracerRecursive(const Ray &r, double sigma_a, double sigma_s){
	int number = spheres.size();
	double sigma_t = sigma_a + sigma_s;
	double continueprob = 0.6;
	double q = 1 - continueprob;

	double t;
	int id = 0;
	if (!intersect(r, t, id)) {
		return {0, 0, 0};
	}

	//RULETA RUSA terminar el camino
	if (erand48(seed) < q) {
		return {0, 0, 0};
	}

	//determinar el limite del rayo previamente
	const Point xs = r.o + r.d * t;
	Vector normalXS = xs - spheres[id].p;
	normalXS.normalize();

	//determinar transmitancia del recorrido
	const double TrActual = transmitance(r.o, xs, sigma_t); //tr de x a xs



	double pSuccess;
	double d;

	/* Metodo de muestreo */
	//escoger entre uno u otro


	/*

	//free flight sampling
	const double d = freeFlightSample(sigma_t);
	pSuccess = freeFlightProb(sigma_t, d)*(1.0-TrActual);
	*/

	 //pre operaciones para seleccionar fuentes
	//guarda en un arreglo los indices de las fuentes de luz
	int arr[4] = {-1, -1, -1, -1};
	int j = 0, count = 0;
	for (int i = 0; i < number; i++) {
		if (spheres[i].radiance.x > 0 || spheres[i].radiance.y > 0 || spheres[i].radiance.z > 0) {
			arr[j++] = i;
		}
	}

	count = j;

	if (count == 0) {
		return {0, 0, 0};
	}
	//probabilidad uniforme para cada fuente
	double prob = 1.0/count;




	int idsource = arr[static_cast<int>(erand48(seed) * count)];


	//operaciones en caso de equiangular sampling
	double D = 0;
	double thetaA = 0;
	double thetaB = 0;
	auto x0 = Point();
	equiAngularParams(idsource, xs, x0, r,D, thetaA, thetaB);
	d = equiAngularSample(D, thetaA, thetaB);
	//probabilidad con equiangular en caso de success
	pSuccess = equiAngularProb(D, thetaA, thetaB, d)*(1.0 - TrActual);

	/*
	else
	{
		//no es posible aplicar equiangular usa free flight para calcular d y la probabilidad
		d = freeFlightSample(sigma_t);
		//probabilidad con freeflight
		pSuccess = freeFlightProb(sigma_t, d)*(1.0-TrActual);

	}*/



	//if(d>t){
	if ( erand48(seed) <= TrActual) {

		auto Ld = Color(0, 0, 0);
		for (int light = 0; light < number; light++) {
			if (spheres[light].r == 0) {
				double Trs = transmitance(xs, spheres[light].p, sigma_t); //es valido en tanto sea puntual
				Ld = Ld + pLight(spheres[id], xs, normalXS, r.d, spheres[light].radiance, spheres[light].p, spheres[id].alpha)*Trs;
				//incluir el efecto de la transmitancia entre xs y p (single scattering)
			}
		}




		Color Lo = Ld*(1/continueprob);//TrActual*(1/pdf) = 1; no tiene caso calcular la transmitancia debido a que al dividirla
		//entre la probabilidad de fallo que es la misma transmitancia, queda en 1


		Vector wi;
		double samplingProbability;
		Color fsActual = bdsf(wi, r.d, normalXS, samplingProbability, id);
		wi.normalize();

		double cosine = normalXS.dot(wi);

		return Lo + fsActual.mult(volumetricPathTracerRecursive(Ray(xs, wi), sigma_a, sigma_s))*(1/continueprob)*cosine*(1/samplingProbability);

	}
	else{

		//const Point xt = r.o + r.d * d; //para freeflight
		//para equiangular sampling xt tiene como referencia x0
		const Point xt = x0 + r.d * d;

		const double T = transmitance(r.o, xt, sigma_t);

		auto Li_parcial = Color(0, 0, 0);


		int lightSource = arr[static_cast<int>(erand48(seed) * count)];


		//calculando un single scattering
		if (spheres[lightSource].r == 0) {
			Point light = spheres[lightSource].p;
			if (visibility(light, xt)) {
				Color Le = spheres[lightSource].radiance;
				const double distanceLight = (light - xt).dot(light - xt);
				Le = Le * (1.0 / distanceLight);
				const Color Ls =  Le * exp(sigma_t*-1.0*sqrt(D*D+d*d));//transmitance(xt, light, sigma_t);
				Li_parcial = Ls * T * sigma_s * isotropicPhaseFunction() * (1/prob) ; //este es el calculo de la luz directa, hasta aqui
			}

		}

		const Vector wi_new = isotropicPhaseSample(); //nueva direccion para siguiente rayo
		//pdf *= isotropicPhaseProb(); la fase cancela a la probabilidad


		return Li_parcial*(1/pSuccess)*(1/continueprob) + volumetricPathTracerRecursive(Ray(xt, wi_new), sigma_a, sigma_s)*sigma_s*T*(1/continueprob)*(1/pSuccess);

	}
}


//version implicita basica con fines de comparacion
inline Color implicitVPTracerRecursive(const Ray &r, double sigma_a, double sigma_s)
{
	int number = spheres.size();
	double sigma_t = sigma_a + sigma_s;
	double continueprob = 0.6;
	double q = 1 - continueprob;

	double TrActual;
	Point xs;

	//RULETA RUSA terminar el camino
	if (erand48(seed) < q) {
		return {0, 0, 0};
	}


	double t;
	int id = 0;
	if (!intersect(r, t, id)) {
		//no regreses cero, propon un limite diferente para x_s
		t = MAXFLOAT; //se ajusta la distancia maxima
		//adicionalmente la probabilida de fallo se debe ajustar a 0 en este caso
		TrActual = 0; //la probabilidad de salir del medio es 0 porque no hay superficie y toda la escena esta inmersa ahi
		//¿que hacemos con id?
		//id se mantiene como 0, los calculos referentes a la superficie no se efectuarán
		//determinar el limite del rayo
		xs = r.o + r.d * t;
	}
	else
	{
		//determinar el limite del rayo
		xs = r.o + r.d * t;
		//determinar transmitancia del recorrido
		TrActual = transmitance(r.o, xs, sigma_t); //tr de x a xs
	}

	Vector normalXS = xs - spheres[id].p;
	normalXS.normalize();


	//pre operaciones para seleccionar fuentes
	//guarda en un arreglo los indices de las fuentes de luz
	int arr[4] = {-1, -1, -1, -1};
	int j = 0, count = 0;
	for (int i = 0; i < number; i++) {
		if (spheres[i].radiance.x > 0 || spheres[i].radiance.y > 0 || spheres[i].radiance.z > 0) {
			arr[j++] = i;
		}
	}

	count = j;

	if (count == 0) {
		return {0, 0, 0};
	}
	//probabilidad uniforme para cada fuente
	//esta probalidad no se debe introducir de este modo, borrar
	double prob = 1.0/count;


	int idsource = arr[static_cast<int>(erand48(seed) * count)];


	//operaciones en caso de equiangular sampling
	double D = 0;
	double thetaA = 0;
	double thetaB = 0;
	auto x0 = Point();
	equiAngularParams(idsource, xs, x0, r,D, thetaA, thetaB);
	const double d = equiAngularSample(D, thetaA, thetaB);

	//probabilidad con equiangular en caso de success  --- alternativa sin MIS
	const double pSuccess = equiAngularProb(D, thetaA, thetaB, d) * (1.0 - TrActual);

	/*
	// --- MIS: calcular el PDF combinado sobre todas las fuentes ---
	// En vez de usar solo el PDF de la fuente idsource, sumamos las contribuciones
	double p_total = 0.0;
	for (int i = 0; i < count; i++) {
		int lightID = arr[i];
		double D_i = 0, thetaA_i = 0, thetaB_i = 0;
		Point x0_i;
		// Se recalculan los parámetros para cada fuente
		equiAngularParams(lightID, xs, x0_i, r, D_i, thetaA_i, thetaB_i);
		p_total += prob * equiAngularProb(D_i, thetaA_i, thetaB_i, d);
	}
	p_total *= 1.0-TrActual;
*/
	if ( erand48(seed) <= TrActual) {


		if (spheres[id].radiance.x>0 || spheres[id].radiance.y>0 || spheres[id].radiance.z>0){
			return spheres[id].radiance;
		}


		Vector wi;
		double samplingProbability;
		const Color fsActual = bdsf(wi, r.d, normalXS, samplingProbability, id);
		wi.normalize();

		const double cosine = normalXS.dot(wi);

		//importante, la transmitancia se omite debido a que es la probabilidad de salir, por tanto T/pFail = 1

		return fsActual.mult(implicitVPTracerRecursive(Ray(xs, wi), sigma_a, sigma_s))*(1/continueprob)*cosine*(1/samplingProbability);

	}
	else{

		//para equiangular sampling xt tiene como referencia x0
		const Point xt = x0 + r.d * d;

		const double T = transmitance(r.o, xt, sigma_t);


		const Vector wi_new = isotropicPhaseSample(); //nueva direccion para siguiente rayo

		return implicitVPTracerRecursive(Ray(xt, wi_new), sigma_a, sigma_s)*sigma_s*T*(1/continueprob)*(1/pSuccess);
	}

}

//version implicita para funcionar de base
inline Color implicitVPTracerRecursiveFree(const Ray &r, double sigma_a, double sigma_s)
{
	double sigma_t = sigma_a + sigma_s;
	double continueprob = 0.6;
	double q = 1 - continueprob;

	Point xs;
	Vector normalXS;
	double TrActual;

	//RULETA RUSA terminar el camino
	if (erand48(seed) < q) {
		return {0, 0, 0};
	}

	double t;
	int id = 0;
	if (!intersect(r, t, id)) {
		t = MAXFLOAT; //equivale definir t como infinito para evitar que se active la condicion de fallo, NO HAY SUPERFICIE
		TrActual = 0; //se debe despreciar la transmitancia, no hay superficie	PROBABILIDAD DE FALLO = 0
		//return {0, 0, 0};
	}
	else
	{
		//determinar el limite del rayo
		xs = r.o + r.d * t;
		normalXS = xs - spheres[id].p;
		normalXS.normalize();

		//determinar transmitancia del recorrido
		TrActual = transmitance(r.o, xs, sigma_t); //tr de x a xs

	}


	//free flight sampling
	const double d = freeFlightSample(sigma_t);
	const double pSuccess = freeFlightProb(sigma_t, d)*(1.0-TrActual);

	if(d>t){ //condicion de fallo

		if (spheres[id].radiance.x>0 || spheres[id].radiance.y>0 || spheres[id].radiance.z>0){
			return spheres[id].radiance;
		}

		Vector wi;
		double samplingProbability;
		Color fsActual = bdsf(wi, r.d, normalXS, samplingProbability, id);
		wi.normalize();

		double cosine = normalXS.dot(wi);

		//importante, la transmitancia se omite debido a que es la probabilidad de salir, por tanto T/pFail = 1

		return fsActual.mult(implicitVPTracerRecursiveFree(Ray(xs, wi), sigma_a, sigma_s))*(1/continueprob)*cosine*(1/samplingProbability);

	}
	else{

		const Point xt = r.o + r.d * d; //para freeflight

		const double T = transmitance(r.o, xt, sigma_t);

		const Vector wi_new = isotropicPhaseSample(); //nueva direccion para siguiente rayo


		return implicitVPTracerRecursiveFree(Ray(xt, wi_new), sigma_a, sigma_s)*sigma_s*T*(1/continueprob)*(1/pSuccess);

	}

}

//version explicita para el VPT con equiangular, es decir, con conexión single scattering
//sustituye a la version anterior volumetricPathTracerRecursive
inline Color explicitVPTracerRecursive(const Ray &r, double sigma_a, double sigma_s, int profundidad)
{

	int number = spheres.size();
	double sigma_t = sigma_a + sigma_s;
	double continueprob = 0.6;
	double q = 1 - continueprob;

	//RULETA RUSA terminar el camino
	if (erand48(seed) < q) {
		return {0, 0, 0};
	}

	double TrActual;
	Point xs;


	double t;
	int id = 0;
	if (!intersect(r, t, id)) {
		//no regreses cero, propon un limite diferente para x_s
		t = MAXFLOAT; //se ajusta la distancia maxima
		//adicionalmente la probabilida de fallo se debe ajustar a 0 en este caso
		TrActual = 0; //la probabilidad de salir del medio es 0 porque no hay superficie y toda la escena esta inmersa ahi
		//determinar el limite del rayo
		xs = r.o + r.d * t;
	}
	else
	{
		//determinar el limite del rayo
		xs = r.o + r.d * t;

		//determinar transmitancia del recorrido
		TrActual = transmitance(r.o, xs, sigma_t); //tr de x a xs
	}



	Vector normalXS = xs - spheres[id].p;
	normalXS.normalize();


	//pre operaciones para seleccionar fuentes
	//guarda en un arreglo los indices de las fuentes de luz
	int arr[4] = {-1, -1, -1, -1};
	int j = 0, count = 0;
	for (int i = 0; i < number; i++) {
		if (spheres[i].radiance.x > 0 || spheres[i].radiance.y > 0 || spheres[i].radiance.z > 0) {
			arr[j++] = i;
		}
	}

	count = j;

	if (count == 0) {
		return {0, 0, 0};
	}
	//probabilidad uniforme para cada fuente
	//esta probalidad no se debe introducir de este modo, borrar
	const double probSource = 1.0/count;


	int idsource = arr[static_cast<int>(erand48(seed) * count)];

	//operaciones en caso de equiangular sampling
	double D = 0;
	double thetaA = 0;
	double thetaB = 0;
	//auto x0 = Point();
	double sample_d;


	/*
	equiAngularParams(idsource, xs, x0, r,D, thetaA, thetaB);

	const double d = equiAngularSample(D, thetaA, thetaB);
*/
	const double d_final = equiAngularParams2(idsource, t, r, D, thetaA, thetaB, sample_d);
	//probabilidad con equiangular en caso de success  --- alternativa sin MIS
	const double pSuccess = equiAngularProb(D, thetaA, thetaB, sample_d) * (1.0 - TrActual);


	if (erand48(seed) <= TrActual) {

		//una versión explicita no debe regresar luz al llegar a la fuente
		if (spheres[id].radiance.x>0 || spheres[id].radiance.y>0 || spheres[id].radiance.z>0){
			//solo regresar radiancia en la primera llamada
			if (profundidad>0)
			{
				return {0, 0, 0};
			}
			return spheres[id].radiance;
		}

		//obtener la luz directa a superfice con plight o MIS - esto equivale a single scattering

		//vamos a aprovechar el muestreo para seleccionar una fuente y solo calcularemos la contribución de idsource / probabilidad
		double Trs = transmitance(xs, spheres[idsource].p, sigma_t); //es valido en tanto sea puntual
		const Color Ld_parcial = pLight(spheres[id], xs, normalXS, r.d, spheres[idsource].radiance, spheres[idsource].p, spheres[id].alpha)*Trs*(1/probSource);

		//para luces esfericas
		//NOTA: esta funcion ya toma en cuenta todas las luces no esfericas, no es para una en particular
		const Color Ld = MISv2(spheres[id], xs, normalXS, r.d, spheres[id].alpha, sigma_t);

		//operaciones para obtener el camino por donde continúa
		Vector wi;
		double samplingProbability;
		const Color fsActual = bdsf(wi, r.d, normalXS, samplingProbability, id);
		wi.normalize();

		const double cosine = normalXS.dot(wi);

		//importante, la transmitancia se omite debido a que es la probabilidad de salir, por tanto T/pFail = 1

		return (Ld_parcial + Ld)*(1/continueprob) + fsActual.mult(explicitVPTracerRecursive(Ray(xs, wi), sigma_a, sigma_s, profundidad+1))*(1/continueprob)*cosine*(1/samplingProbability);

	}
	else{

		//const Point xt = x0 + r.d * d;
		const Point xt = r.o + r.d * d_final; // reparametrizado

		const double T = transmitance(r.o, xt, sigma_t);

		//obtener la iluminación directa aqui haciendo single scattering para luz puntual o de area - angulo solido
		const Color Ld = singleScattering(xt, idsource, sigma_t, sigma_s, T, probSource);

		//operaciones par aobtner la nueva direccion
		const Vector wi_new = isotropicPhaseSample(); //nueva direccion para siguiente rayo


		//sumar Ld transformado aqui
		return Ld*(1/pSuccess)*(1/continueprob) + explicitVPTracerRecursive(Ray(xt, wi_new), sigma_a, sigma_s, profundidad+1)*sigma_s*T*(1/continueprob)*(1/pSuccess);
	}

}

//version explicita para el VPT con equiangular, es decir, con conexión single scattering
//sustituye a la version anterior volumetricPathTracerRecursive
inline Color explicitVPTracerRecursiveFree(const Ray &r, double sigma_a, double sigma_s, int profundidad)
{

	int number = spheres.size();
	double sigma_t = sigma_a + sigma_s;
	double continueprob = 0.6;
	double q = 1 - continueprob;

	//RULETA RUSA terminar el camino
	if (erand48(seed) < q) {
		return {0, 0, 0};
	}


	double t;
	int id = 0;
	if (!intersect(r, t, id)) {
		//no regreses cero, propon un limite diferente para x_s
		t = MAXFLOAT; //se ajusta la distancia maxima

	}
	const Point xs = r.o + r.d * t;
	//const double TrActual = transmitance(r.o, xs, sigma_t); //tr de x a xs
	//NOTA: termina estando implicita

	Vector normalXS = xs - spheres[id].p;
	normalXS.normalize();


	//pre operaciones para seleccionar fuentes
	//guarda en un arreglo los indices de las fuentes de luz
	int arr[4] = {-1, -1, -1, -1};
	int j = 0, count = 0;
	for (int i = 0; i < number; i++) {
		if (spheres[i].radiance.x > 0 || spheres[i].radiance.y > 0 || spheres[i].radiance.z > 0) {
			arr[j++] = i;
		}
	}

	count = j;

	if (count == 0) {
		return {0, 0, 0};
	}
	//probabilidad de seleccionar una fuente de todas, de esta forma es valida y esta normalizada únicamente ...
	// para aplicaciones que traten sobre capturar radiancia JAMAS para acercarse mas o menos con un xt a una fuente
	const double probSource = 1.0/count;


	int idsource = arr[static_cast<int>(erand48(seed) * count)];

	//free flight sampling
	const double d = freeFlightSample(sigma_t);

	if (d>t) {

		//una versión explicita no debe regresar luz al llegar a la fuente
		if (spheres[id].radiance.x>0 || spheres[id].radiance.y>0 || spheres[id].radiance.z>0){
			//solo regresar radiancia en la primera llamada
			if (profundidad>0)
			{
				return {0, 0, 0};
			}
			return spheres[id].radiance;
		}

		//obtener la luz directa a superfice con plight o MIS - esto equivale a single scattering

		//vamos a aprovechar el muestreo para seleccionar una fuente y solo calcularemos la contribución de idsource / probabilidad
		const double Trs = transmitance(xs, spheres[idsource].p, sigma_t); //es valido en tanto sea puntual
		const Color Ld_parcial = pLight(spheres[id], xs, normalXS, r.d, spheres[idsource].radiance, spheres[idsource].p, spheres[id].alpha)*Trs*(1/probSource);

		//para luces esfericas
		//NOTA: esta funcion ya toma en cuenta todas las luces no esfericas, no es para una en particular
		const Color Ld = MISv2(spheres[id], xs, normalXS, r.d, spheres[id].alpha, sigma_t);

		//operaciones para obtener el camino por donde continúa
		Vector wi;
		double samplingProbability;
		const Color fsActual = bdsf(wi, r.d, normalXS, samplingProbability, id);
		wi.normalize();

		const double cosine = normalXS.dot(wi);

		//importante, la transmitancia se omite debido a que es la probabilidad de salir, por tanto T/pFail = 1

		return (Ld_parcial + Ld)*(1/continueprob) + fsActual.mult(explicitVPTracerRecursiveFree(Ray(xs, wi), sigma_a, sigma_s, profundidad+1))*(1/continueprob)*cosine*(1/samplingProbability);

	}
	else{

		//para equiangular sampling xt tiene como referencia x0
		const Point xt = r.o + r.d * d;

		//const double T = transmitance(r.o, xt, sigma_t);
		//NOTA; por simplificaciones propuestas por el PBRT se prescinde de transmitancia de forma explicita

		//obtener la iluminación directa aqui haciendo single scattering para luz puntual o de area - angulo solido
		const Color Ld = freeSingleScattering(xt, idsource, sigma_t, probSource);

		//operaciones par aobtner la nueva direccion
		const Vector wi_new = isotropicPhaseSample(); //nueva direccion para siguiente rayo


		//sumar Ld transformado aqui
		return Ld*(sigma_s/sigma_t)*(1/continueprob) + explicitVPTracerRecursiveFree(Ray(xt, wi_new), sigma_a, sigma_s, profundidad+1)*(sigma_s/sigma_t)*(1/continueprob);
	}

}

inline Color iterativeVPTracerFree(const Ray &r, double sigma_a, double sigma_s) {
    struct StackFrame {
        Ray ray;
        int profundidad;
        Color pathThroughput;
    };

    std::stack<StackFrame> stack;
    stack.push({r, 0, Color(1, 1, 1)});

    Color finalColor(0, 0, 0);
    double sigma_t = sigma_a + sigma_s;
    double continueprob = 0.6;
    double q = 1 - continueprob;

    while (!stack.empty()) {
        auto [currentRay, profundidad, pathThroughput] = stack.top();
        stack.pop();

        if (erand48(seed) < q) continue;

        double t;
        int id = 0;
        if (!intersect(currentRay, t, id)) {
            t = MAXFLOAT;
        }
        Point xs = currentRay.o + currentRay.d * t;
        Vector normalXS = xs - spheres[id].p;
        normalXS.normalize();

        int arr[4] = {-1, -1, -1, -1};
        int count = 0;
        for (int i = 0; i < spheres.size(); i++) {
            if (spheres[i].radiance.x > 0 || spheres[i].radiance.y > 0 || spheres[i].radiance.z > 0) {
                arr[count++] = i;
            }
        }

        if (count == 0) continue;
        double probSource = 1.0 / count;
        int idsource = arr[static_cast<int>(erand48(seed) * count)];

        double d = freeFlightSample(sigma_t);

        if (d > t) {
            if (spheres[id].radiance.x > 0 || spheres[id].radiance.y > 0 || spheres[id].radiance.z > 0) {
                if (profundidad == 0) {
                    finalColor = spheres[id].radiance.mult(pathThroughput);
                }
                continue;
            }

            double Trs = transmitance(xs, spheres[idsource].p, sigma_t);
            Color Ld_parcial = pLight(spheres[id], xs, normalXS, currentRay.d, spheres[idsource].radiance, spheres[idsource].p, spheres[id].alpha) * Trs * (1 / probSource);
            Color Ld = MISv2(spheres[id], xs, normalXS, currentRay.d, spheres[id].alpha, sigma_t);

            Vector wi;
            double samplingProbability;
            Color fsActual = bdsf(wi, currentRay.d, normalXS, samplingProbability, id);
            wi.normalize();
            double cosine = normalXS.dot(wi);

            finalColor = finalColor+  (Ld_parcial + Ld).mult(pathThroughput) * (1 / continueprob);
            stack.push({Ray(xs, wi), profundidad + 1, pathThroughput.mult(fsActual) * (1 / continueprob) * cosine * (1 / samplingProbability)});
        } else {
            Point xt = currentRay.o + currentRay.d * d;
            Color Ld = freeSingleScattering(xt, idsource, sigma_t, probSource);
            Vector wi_new = isotropicPhaseSample();

            finalColor = finalColor + Ld.mult(pathThroughput) * (sigma_s / sigma_t) * (1 / continueprob);
            stack.push({Ray(xt, wi_new), profundidad + 1, pathThroughput * (sigma_s / sigma_t) * (1 / continueprob)});
        }

    }

    return finalColor;
}



//implementa equiangular sampling y free flight al mismo tiempo. El principal objetivo es eliminar los problemas de normalizacion de equiangular sampling
inline Color MISVPTTracerRecursive(const Ray &r, double sigma_a, double sigma_s, int profundidad)
{

	int number = spheres.size();
	double sigma_t = sigma_a + sigma_s;
	double continueprob = 0.6;
	double q = 1 - continueprob;

	//RULETA RUSA terminar el camino
	if (erand48(seed) < q) {
		return {0, 0, 0};
	}

	double TrActual;
	Point xs;


	double t;
	int id = 0;
	if (!intersect(r, t, id)) {
		//no regreses cero, propon un limite diferente para x_s
		t = MAXFLOAT; //se ajusta la distancia maxima
		//adicionalmente la probabilidad de fallo se debe ajustar a 0 en este caso
		TrActual = 0; //la probabilidad de salir del medio es 0 porque no hay superficie y toda la escena esta inmersa ahi
		//determinar el limite del rayo
		xs = r.o + r.d * t;
	}
	else
	{
		//determinar el limite del rayo
		xs = r.o + r.d * t;

		//determinar transmitancia del recorrido
		TrActual = transmitance(r.o, xs, sigma_t); //tr de x a xs
	}



	Vector normalXS = xs - spheres[id].p;
	normalXS.normalize();


	//pre operaciones para seleccionar fuentes
	//guarda en un arreglo los indices de las fuentes de luz
	int arr[4] = {-1, -1, -1, -1};
	int j = 0, count = 0;
	for (int i = 0; i < number; i++) {
		if (spheres[i].radiance.x > 0 || spheres[i].radiance.y > 0 || spheres[i].radiance.z > 0) {
			arr[j++] = i;
		}
	}

	count = j;

	if (count == 0) {
		return {0, 0, 0};
	}
	//probabilidad uniforme para cada fuente
	//esta probalidad no se debe introducir de este modo, borrar
	const double probSource = 1.0/count;


	int idsource = arr[static_cast<int>(erand48(seed) * count)];

	//operaciones en caso de equiangular sampling
	double D = 0;
	double thetaA = 0;
	double thetaB = 0;
	//auto x0 = Point();
	double sample_d;

	//calcular una d de decision con free flight sampling
	//const double d_dec = freeFlightSample(sigma_t);
	//obtener la probabilidad de ir la superficie para fines de normalizacion
	const double psurf = exp(sigma_t*t*-1.0); // es equivalente a la transmitancia TrActual

	const double d_final = equiAngularParams2(idsource, t, r, D, thetaA, thetaB, sample_d);
	//el complemento de psurf es la probabilidad de quedarse en el medio
	const double pdf = equiAngularProb(D, thetaA, thetaB, sample_d) * (1-psurf);


	if (erand48(seed) < psurf) { //se sustituye la condición impuesta por una organica


		if (spheres[id].radiance.x>0 || spheres[id].radiance.y>0 || spheres[id].radiance.z>0){

			if (profundidad>0)
			{
				return {0, 0, 0};
			}
			//solo existe la posibilidad de regresar la radiancia en la primera llamada
			//debido a que lleva una direccion fija que no esta sesgada hacia la fuente
			return spheres[id].radiance;
		}

		//obtener la luz directa a superfice con plight o MIS - esto equivale a single scattering

		//vamos a aprovechar el muestreo para seleccionar una fuente y solo calcularemos la contribución de idsource / probabilidad
		double Trs = transmitance(xs, spheres[idsource].p, sigma_t); //es valido en tanto sea puntual
		const Color Ld_parcial = pLight(spheres[id], xs, normalXS, r.d, spheres[idsource].radiance, spheres[idsource].p, spheres[id].alpha)*Trs*(1/probSource);

		//para luces esfericas
		//NOTA: esta funcion ya toma en cuenta todas las luces no esfericas, no es para una en particular
		const Color Ld = MISv2(spheres[id], xs, normalXS, r.d, spheres[id].alpha, sigma_t);

		//operaciones para obtener el camino por donde continúa
		Vector wi;
		double samplingProbability;
		const Color fsActual = bdsf(wi, r.d, normalXS, samplingProbability, id);
		wi.normalize();

		const double cosine = normalXS.dot(wi);

		//importante, la transmitancia se omite debido a que es la probabilidad de salir, por tanto T/pFail = 1

		return (Ld_parcial + Ld)*(1/continueprob) + fsActual.mult(MISVPTTracerRecursive(Ray(xs, wi), sigma_a, sigma_s, profundidad+1))*(1/continueprob)*cosine*(1/samplingProbability);

	}
	else{

		//const Point xt = x0 + r.d * d;
		const Point xt = r.o + r.d * d_final; // reparametrizado, es decir, va incluida la delta

		const double T = transmitance(r.o, xt, sigma_t);

		//obtener la iluminación directa aqui haciendo single scattering para luz puntual o de area - angulo solido
		const Color Ld = singleScattering(xt, idsource, sigma_t, sigma_s, T, probSource);

		//operaciones par aobtner la nueva direccion
		const Vector wi_new = isotropicPhaseSample(); //nueva direccion para siguiente rayo


		//sumar Ld transformado aqui
		return Ld*(1/pdf)*(1/continueprob) + MISVPTTracerRecursive(Ray(xt, wi_new), sigma_a, sigma_s, profundidad+1)*sigma_s*T*(1/continueprob)*(1/pdf);
	}

}

#endif //VPTSHADEMETHODS_H
