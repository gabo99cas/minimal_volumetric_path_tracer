//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef RAYMARCHINGMETHODS_H
#define RAYMARCHINGMETHODS_H

#include "volumetricBasicFunctions.h"
#include "samplingFunctions.h"

//luz puntual single scattering v1 para ray marching
Color punctualVolumetric(int idsource, Point x, double phase, double sigma_t, double sigma_s){
	Color Ls = Color();
	Color Ld = Color();
	Point light = spheres[idsource].p;
	if(visibilityVPT(light, x)){
		//obtener la radiancia de la fuente de luz
		Color Le = spheres[idsource].radiance;
		//calcular ls distancia entre la fuente de luz y x_new
		double distanceLight = (light-x).dot(light-x);
		//dividir la radiancia por la distancia al cuadrado
		Le = Le*(1/distanceLight);
		//montecarlo de una sola muestra
		Ls = Le*phase*multipleT(x, light, sigma_t);//la probabilidad de muestreo es 1;
		Ld = Ls*sigma_s;
	}
	//en caso de no haber visibilidad
	else	Ld = Color();
	return Ld;

}

//funcion para muestreo explicito en ray marching con iluminacion global
Color rayMarching(const Ray &r, double sigma_t, double sigma_s, double steps, Point &x_new, int &idsource){

	//obtener la distancia de interseccion
	double t;
	int id = 0;
	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() ==
	idsource = id;
	Point x = r.o + r.d*t;
	//se guarda el punto de interseccion
	x_new = x;
	Color Li = Color();
	Color Lo = Color();
	//si el rayo intersecta una fuente de luz, Lo = radiancia * transmitancia
	if(spheres[id].radiance.x>0) {
		//como la fuente ya fue muestreada, no se toma en cuenta en el ray marching
		return Color();
	}
	//divide la distancia entre steps para obtener la distancia de cada segmento
	double step = t/steps;
	//para cada segmento se calcula la transmitancia y la funcion de fase
	for(int i = 0; i<steps; i++){
		Point xt = r.o + r.d*step*i;
		//calcular la transmitancia
		double T = transmitance(x, xt, sigma_t);
		//calcular la funcion de fase
		double phase = isotropicPhaseFunction();
		//calcular Ls aproximado
		//lanzar un shadow rays con angulo solido
		//vector wc de xt a la fuente
		Vector wc = spheres[5].p-xt;
		double normcx = sqrt(wc.dot(wc));
		wc = wc*(1/normcx);
		double costheta_max = sqrt(1-(spheres[5].r/normcx)*(spheres[5].r/normcx));
		//vector wi sampleado con angulo solido
		Vector wi = solidAngle(wc, costheta_max);
		//rayo con origen en xt y direccion wi
		Ray shadowRay = Ray(xt, wi);
		//verificar la visibilidad con la fuente de luz
		//si el shadow ray intersecta en el mismo id que la fuente de luz, entonces hay visibilidad
		//en caso contrario no hay visibilidad
		double t2;
		int id2 = 0;
		intersect(shadowRay, t2, id2);
		//si hay visibilidad, la contribucion es la radiancia de la fuente de luz

		if(id2==5){//se tiene que modificar para que tome en cuenta varias fuentes (workaround)
			//obtener la radiancia de la fuente de luz
			Color Le = spheres[5].radiance;
			//calcular Ls
			//montecarlo de una sola muestra
			Color Ls = Le*(phase*transmitance(xt, spheres[5].p, sigma_t));
			//probabilidad de muestreo
			double prob = solidAngleProb(costheta_max);
			//calcular Ls aproximado
			Color Ls_aprox = Ls*(T*1/prob)*sigma_s*step;
			//acumular Ls aproximado
			Li = Li + Ls_aprox;
		}
		//en caso de no haber visibilidad
		else{
			Color Ls_aprox = Color();
			Li = Li + Ls_aprox;
		}

	}
		return Li;


}

//ray marching para volumenes, segmento o paso variable
Color rayMarchingGlobal(const Ray &r, double sigma_a,double sigma_s, double  segmentos){

	double sigma_t = sigma_a+sigma_s;

	//obtener la distancia de interseccion
	double t;
	int id = 0;
	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() ==
	Point x = r.o + r.d*t;
	Color Li = Color();
	Color Lo = Color();
	//si el rayo intersecta una fuente de luz, Lo = radiancia * transmitancia
	if(spheres[id].radiance.x>0) {
		return spheres[id].radiance*transmitance(r.o, x, sigma_t);
		//Lo = spheres[id].radiance*transmitance(r.o, x, sigma_t);

	}


	/*
	//calcular el color con angulo solido
	Vector aux = Vector();
	double costheta_max;
	Color Ld = muestreoSA(spheres[5], x, 5, spheres[id], n, r.d, aux, costheta_max, 0.001);
	//sumar Ld * transmitancia a Lo, solo si LO es 0
	if(Lo.x==0) Lo = Lo + Ld*transmitance(r.o, x, sigma_t);
	//si se desea aplicar multiple scattering, se debe aplicar ray marching despues de muestrear la direccion de angulo solido
	//sustituir en muestreo SA el montecarlo normal por ray marching
	//Le = emision de la fuente de luz * transmitancia + suma de reiman ... Le = Li
	//es decir luz entrante en el punto x = Li
	*/

	//obtener Le correspondiente a la fuente de luz
	//se debe iterar esta parte para incluir mas rebotes
	Color Le = Color();
	Color fs = Color(1,1,1);
	double factor = 1;
	Color Ld = Color();
	for(int i=0; i<10; i++){


	Color fr = spheres[id].c*(1/M_PI); //se debe calcular la fr de acuerdo al material

	Vector n = (x-spheres[id].p);
	n.normalize();
	//calcular el color con angulo solido
	Vector wc = spheres[5].p-x;
	double normcx = sqrt(wc.dot(wc));
	wc = wc*(1/normcx);
	double costheta_max = sqrt(1-(spheres[5].r/normcx)*(spheres[5].r/normcx));
	//obtener vector wi
	Vector wi = solidAngle(wc, costheta_max);
	//rayo con origen en x y direccion wi
	Ray shadowRay = Ray(x, wi);
	//verificar la visibilidad con la fuente de luz
	//si el shadow ray intersecta en el mismo id que la fuente de luz, entonces hay visibilidad
	//en caso contrario no hay visibilidad
	double t_aux;
	int id_aux = 0;
	intersect(shadowRay, t_aux, id_aux);
	Point x2 = shadowRay.o + shadowRay.d*t_aux;
	//si hay visibilidad, la contribucion es la radiancia de la fuente de luz
	if(id_aux==5){
		//Le = emision de la fuente de luz
		Le = spheres[5].radiance*transmitance(x, spheres[5].p, sigma_t);

	//std::cout<<Ld.x<<" "<<Ld.y<<" "<<Ld.z<<std::endl;
	//Ld unicamente es la emision de la fuente de luz con direccion a x
	//integrar con montecarlo de angulo solido
	//fr = albedo

	Ld = Le.mult(fr)*(1/solidAngleProb(costheta_max))*n.dot(wi);
	}


	//muestrea una dirección en base al material albedo
	Vector wray = cosineHemispheric(n);
	//probabilidad coseno
	double prob = hemiCosineProb(n.dot(wray));
	//rayo con origen en x y direccion wi
	Ray newray = Ray(x, wray);
	Point x_new = Point();
	//ray marching
	Color Lm = rayMarching(newray, sigma_t, sigma_s, segmentos, x_new, id); //se actualiza la id del impacto

	//acumular Lm
	Ld=Ld + Lm.mult(fr)*n.dot(wray)*(1/prob);
	//acumular Ld a Lo
	Lo = Lo + Ld.mult(fs)*transmitance(r.o, x, sigma_t)*factor;
	//si Lm es 0, se regresa lo que se acumulo hasta el momento
	if(Lm.x==0 && Lm.y==0 && Lm.z==0) return Lo;
	//guardar la fs en el acumulador
	fs = fs.mult(fr);	 //acumula las fr
	factor = factor*n.dot(wray)*(1/(prob)); //acumula la probabilidad de muestreo
	//actualiza la x
	x = x_new;
	}

	//se divide t entre segmentos para obtener la distancia de cada segmento
	double step = t/segmentos;
	//para cada segmento se calcula la transmitancia y la funcion de fase
	for(int i = 0; i<segmentos; i++){
		Point xt = r.o + r.d*step*i;
		//calcular la transmitancia
		double T = transmitance(x, xt, sigma_t);
		//calcular la funcion de fase
		double phase = isotropicPhaseFunction();
		//calcular Ls aproximado
		//lanzar un shadow rays con angulo solido
		//vector wc de xt a la fuente
		Vector wc = spheres[5].p-xt;
		double normcx = sqrt(wc.dot(wc));
		wc = wc*(1/normcx);
		double costheta_max = sqrt(1-(spheres[5].r/normcx)*(spheres[5].r/normcx));
		//vector wi sampleado con angulo solido
		Vector wi = solidAngle(wc, costheta_max);
		//rayo con origen en xt y direccion wi
		Ray shadowRay = Ray(xt, wi);
		//verificar la visibilidad con la fuente de luz
		//si el shadow ray intersecta en el mismo id que la fuente de luz, entonces hay visibilidad
		//en caso contrario no hay visibilidad
		double t2;
		int id2 = 0;
		intersect(shadowRay, t2, id2);
		//si hay visibilidad, la contribucion es la radiancia de la fuente de luz

		if(id2==5){//se tiene que modificar para que tome en cuenta varias fuentes (workaround)
			//obtener la radiancia de la fuente de luz
			Color Le = spheres[5].radiance;
			//calcular Ls
			//montecarlo de una sola muestra
			Color Ls = Le*(phase*transmitance(xt, spheres[5].p, sigma_t));
			//probabilidad de muestreo
			double prob = solidAngleProb(costheta_max);
			//calcular Ls aproximado
			Color Ls_aprox = Ls*(T*1/prob)*sigma_s*step;
			//acumular Ls aproximado
			Li = Li + Ls_aprox;
		}
		//en caso de no haber visibilidad, la contribucion es 0
		else{
			Color Ls_aprox = Color();
			Li = Li + Ls_aprox;
		}



}
	return Li+Lo;
}




//ray marching de paso constante
Color rayMarching2(const Ray &r, double sigma_a, double sigma_s, double step, int idsource){

	//obtener la distancia de interseccion
	double t;
	int id = 0;
	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() ==
	Point x = r.o + r.d*t;
	Color Li = Color();
	Color Lo = Color();
	//si el rayo intersecta una fuente de luz, Lo = radiancia * transmitancia
	if(spheres[id].radiance.x>0) {
		//return spheres[id].radiance*transmitance(r.o, x, sigma_t);
		Lo = spheres[id].radiance*transmitance(r.o, x, sigma_a+sigma_s);
	}
	//distancia entre x y la fuente de luz divida entre segmento para obtener el nimero de segmentos
	double steps = t/step;
	//ciclo para cada segmento
	for(int i = 0; i<steps; i++){
		Point xt = r.o + r.d*step*i;
		//calcular la transmitancia
		double T = transmitance(x, xt, sigma_a+sigma_s);
		//calcular la funcion de fase
		double phase = isotropicPhaseFunction();

		//calcular Ls aproximado
		//lanzar un shadow rays con angulo solido
		//vector wc de xt a la fuente
		Vector wc = spheres[idsource].p-xt;
		double normcx = sqrt(wc.dot(wc));
		wc = wc*(1/normcx);
		double costheta_max = sqrt(1-(spheres[idsource].r/normcx)*(spheres[idsource].r/normcx));
		//vector wi sampleado con angulo solido
		Vector wi = solidAngle(wc, costheta_max);
		//rayo con origen en xt y direccion wi
		Ray shadowRay = Ray(xt, wi);
		//verificar la visibilidad con la fuente de luz
		//si el shadow ray intersecta en el mismo id que la fuente de luz, entonces hay visibilidad
		//en caso contrario no hay visibilidad
		double t2;
		int id2 = 0;
		intersect(shadowRay, t2, id2);
		//si hay visibilidad, la contribucion es la radiancia de la fuente de luz

		if(id2==idsource){//se tiene que modificar para que tome en cuenta varias fuentes (workaround)
			//obtener la radiancia de la fuente de luz
			Color Le = spheres[idsource].radiance;
			//calcular Ls
			//montecarlo de una sola muestra
			Color Ls = Le*(phase*transmitance(xt, spheres[idsource].p, sigma_a+sigma_s));
			//probabilidad de muestreo
			double prob = solidAngleProb(costheta_max);
			//calcular Ls aproximado
			Color Ls_aprox = Ls*(T*1/prob)*sigma_s*step;
			//acumular Ls aproximado
			Li = Li + Ls_aprox;
		}
		//en caso de no haber visibilidad, la contribucion es 0
		else{
			Color Ls_aprox = Color();
			Li = Li + Ls_aprox;
		}

}
	return Li+Lo;
}

//ray marching de paso constante luz puntual
Color rayMarching3(const Ray &r, double sigma_a, double sigma_s, double step, int idsource){

	//obtener la distancia de interseccion
	double t;
	int id = 0;
	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() ==
	Point x = r.o + r.d*t;
	Color Li = Color();
	Color Lo = Color();
	/*
	//si el rayo intersecta una fuente de luz, Lo = radiancia * transmitancia
	if(spheres[id].radiance.x>0) {
		//return spheres[id].radiance*transmitance(r.o, x, sigma_t);
		Lo = spheres[id].radiance*transmitance(r.o, x, sigma_a+sigma_s);
	}
	*/
	//distancia entre x y la fuente de luz divida entre segmento para obtener el nimero de segmentos
	double steps = t/step;
	//ciclo para cada segmento
	for(int i = 0; i<steps; i++){
		Point xt = r.o + r.d*step*i;
		//calcular la transmitancia
		double T = transmitance(x, xt, sigma_a+sigma_s);
		//calcular la funcion de fase
		double phase = isotropicPhaseFunction();
		//calcular Ls aproximado
		//lanzar un shadow rays a la fuente puntual
		//vector wc de xt a la fuente
		Vector wc = spheres[idsource].p-xt;
		//obtener la distancia entre xt y la fuente
		double normwc = wc.dot(wc);
		//verificar la visibilidad con la fuente de luz

		if(visibility(spheres[idsource].p,xt)){//se tiene que modificar para que tome en cuenta varias fuentes (workaround)
			//obtener la radiancia de la fuente de luz
			Color Le = spheres[idsource].radiance*(1/normwc);
			//calcular Ls
			//montecarlo de una sola muestra
			Color Ls = Le*(phase*transmitance(xt, spheres[idsource].p, sigma_a+sigma_s));

			//calcular Ls aproximado
			Color Ls_aprox = Ls*(T)*sigma_s*step;
			//acumular Ls aproximado
			Li = Li + Ls_aprox;
		}
		//en caso de no haber visibilidad, la contribucion es 0
		else{
			Color Ls_aprox = Color();
			Li = Li + Ls_aprox;
		}

}
	return Li;
}

#endif //RAYMARCHINGMETHODS_H
