//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef VPTSHADEMETHODS_H
#define VPTSHADEMETHODS_H

#include "samplingFunctions.h"
#include "volumetricBasicFunctions.h"
#include "vptSamplingFunctions.h"
#include "rayMarchingMethods.h"

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

inline Color volumetricPathTracerIterative(const Ray &r, double sigma_a, double sigma_s) {
    Color Li, Lo, Ls, Ld, fsActual;
	double scatterFactor = 1.0;
    std::stack<std::tuple<Ray, double, double, int, double, Color>> stack;
    stack.emplace(r, sigma_a, sigma_s, 0, 1.0, Color(1,1,1));

    while (!stack.empty()) {
        auto [ray, sigma_a, sigma_s, profundidad, factor, pathThroughput] = stack.top();
        stack.pop();

        double sigma_t = sigma_a + sigma_s;
        double continueprob = 0.7;
        double q = 1 - continueprob;

        double t;
        int id = 0;
        if (!intersect(ray, t, id)) {
            break;
        }

        Point xs = ray.o + ray.d * t;
        Vector normalXS = xs - spheres[id].p;
        normalXS.normalize();

        Lo = Color(0, 0, 0);
        int number = spheres.size();
        for (int light = 0; light < number; light++) {
            if (spheres[light].r == 0) {
                Lo = Lo + pLight(spheres[id], xs, normalXS, ray.d, spheres[light].radiance, spheres[light].p, spheres[id].alpha);
            }
        }

        double TrActual = transmitance(ray.o, xs, sigma_t);
    	

    	Lo = pathThroughput.mult(Lo)*factor*TrActual;

    	if (erand48(seed) < q) {
    		Li = Li + Lo*(1/q);
    		continue; //regresa negro
    	}



        Vector wi;
        double samplingProbability;
        fsActual = bdsf(wi, ray.d, normalXS, samplingProbability, id);
        wi.normalize();

        double cosine = normalXS.dot(wi);



        Ray newRay = Ray(xs, wi);
    	//calcular el path troughput, ahora esto se agrega en el stack directamente para que este disponible en next iteracion
    	//Accum = fs1.mult(Accum);
    	//factor = 1/samplingProbability*cosine*factor;
    	//agregamos en el stack los valores que corresponden a la siguiente iteracion
        stack.emplace(newRay, sigma_a, sigma_s, profundidad + 1, factor*1/samplingProbability*cosine*(1/continueprob)*TrActual, pathThroughput.mult(fsActual));

    	Li = Lo + Li;
    	//actualizando el factor despues para que la probabilidad de continuar tenga efecto hasta la siguiente iteracion
    	//esto no es definitivo y es para coincidir con el iterativePathTracer original
    	//factor = factor*(1/continueprob);



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

        int idsource = arr[(int)(erand48(seed) * count)];

    	/* Metodo de muestreo */
		//escoger entre uno u otro

		//free flight sampling

    	double d = freeFlightSample(sigma_t);
    	double pdf;


		/*
    	//operaciones en caso de equiangular sampling
    	double D = 0;
    	double thetaA = 0;
    	double thetaB = 0;
    	Point x0 = Point();
    	equiAngularParams(idsource, xs, x0, ray,D, thetaA, thetaB);
    	double d = equiAngularSample(D, thetaA, thetaB);
    	double pdf = equiAngularProb(D, thetaA, thetaB, d);
		*/

        if (d >= t) {
        	d = t;
        	pdf = freeFlightProb(sigma_t, d)/pdfFailure(sigma_t, t);
        	//esto solo aplica para el freeflight, equiangular no tiene posibilidad de caer aqui
        	//si falló, saltar la parte de scattering
			continue; //al hacer continue aqui, le estas dando oportunidad al nuevo rayo
        	//return Li; //este cambio concide con la implementacion recrusiva pero lo correcto es hacer continue
        	//podriamos continuar en la version recursiva haciendo un salto
        }
    	else {
    		//omitir esta parte en equiangular
    		pdf = freeFlightProb(sigma_t,d)/pdfSuccess(sigma_t, t);
    	}


        Point xt = ray.o + ray.d * d;
    	//para equiangular sampling xt tiene como referencia x0
    	//Point xt = x0 + ray.d * d;
        Vector wi_new = isotropicPhaseSample();
        double T = transmitance(ray.o, xt, sigma_t);

        if (spheres[idsource].r == 0) {
            Point light = spheres[idsource].p;
            if (visibility(light, xt)) {
                Color Le = spheres[idsource].radiance;
                double distanceLight = (light - xt).dot(light - xt);
                Le = Le * (1 / distanceLight);
                Ls = Le * isotropicPhaseFunction() * transmitance(xt, light, sigma_t);
                Ld = Ls * T * sigma_s * (1 / prob); //este es el calculo de la luz directa, hasta aqui
            }

        }

    	Ld = pathThroughput.mult(Ld) * (1 / pdf) * factor;

        if (erand48(seed) < q) {
        	Li = Li + Ld * (1/q);
        	continue;

        }

    	//Ld = pathThroughput.mult(Ld) * (1 / pdf) * factor; //aplicar el path throughput
    	//es necesario multiplicar por la pdf inmediata para completar el montecarlo

        stack.emplace(Ray(xt, wi_new), sigma_a, sigma_s, profundidad + 1, factor*(1/pdf)*T*sigma_s*(1/continueprob), pathThroughput);

    	//scatterFactor *=  pdf * T * sigma_s; se acumula directamente en el stack para que este disponible hasta la siguiente llamada (iteracion)


		Li = Li + Ld;
    }

    return Li;
}



#endif //VPTSHADEMETHODS_H
