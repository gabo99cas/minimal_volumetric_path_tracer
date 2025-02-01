//
// Created by Gabriel Castillo on 28/01/25.
//

#ifndef MISSAMPLINGFUNCTIONS_H
#define MISSAMPLINGFUNCTIONS_H

#include "samplingFunctions.h"

//calcula la heuristica de potencia B=2
inline double powerHeuristics(double fpdf, double gpdf){
	double f2 = fpdf*fpdf;
	double g2 = gpdf*gpdf;
	return f2/(f2+g2);
}

//devuelve la integral para multiple importance sampling
inline Color MIS(Sphere &obj, Point x, Vector n, Vector wray, double alpha){ //for each ligth source
	Color f, g, montecarlo = Color();
	Vector wiLight, wiBDRF, wo;
	double wf, wg, fpdf, gpdf, costhetaMax;
	int sourceid, sourceid2=0;
	wo = wray*-1;

//se muestrea each light source
	const int size = spheres.size();  // Devuelve el número de elementos en el vector
	for(int light = 0; light<size; light++) {
		if(spheres[light].r>0 && spheres[light].radiance.x>0){
			f = muestreoSA(spheres[light], x, light, obj, n, wray, wiLight, costhetaMax, alpha);
			fpdf = solidAngleProb(costhetaMax);
			if(obj.material==0) gpdf = hemiCosineProb(n.dot(wiLight));
			else if(obj.material==2){
			//calcular el fresnel para la direccion muestreada por luz prob = F en reflexion
				Vector wt = refraxDielectric(1.0, 1.5, wo, n);
				wt.normalize();
				gpdf = fresnelDie(1.0, 1.5, n.dot(wt), n.dot(wo));
				if(erand48(seed)>gpdf){
				gpdf=1-gpdf;
				}
	}
			else{
				Vector wh = wiLight+wo; //global
				wh.normalize();
				gpdf = microFacetProb(wo, wh, alpha, n);
	}
			wf = powerHeuristics(fpdf, gpdf);
			montecarlo = montecarlo + f*wf;

		}

	}

//se muestrea BDSF (es agnostico a la fuente, puede darle a cualquiera)
	if(obj.material==0){
		g = uniform(n, x, obj.c, wiBDRF, sourceid);
		gpdf = hemiCosineProb(n.dot(wiBDRF));
		//para obtener el costhetaMax apropiado, considerar la fuente a la que la direccion apunta
		if(g.x>0 && g.y>0 && g.z>0){ costhetaMax = cosinethetaMax(sourceid, x);
			fpdf = solidAngleProb(costhetaMax);
			wg = powerHeuristics(gpdf, fpdf);
		}
		else wg = 0;


	}
	else if(obj.material==2){
		//calcula iluminacion para dielectrico suave
		g = softDielectric(1.5, 1.0, wo, n, x, sourceid);
		if(g.x>0 && g.y>0 && g.z>0){ costhetaMax = cosinethetaMax(sourceid, x);
		fpdf = solidAngleProb(costhetaMax);
		wg = powerHeuristics(gpdf, fpdf);
		}
		else wg = 0;
	}
	else{
		Vector wh = vectorFacet(alpha);
		coordinateTraspose(n, wo);
		wo.normalize();
		g = microfacet(x, wray, wh, n, obj, alpha, sourceid2);
		gpdf = microFacetProb(wo, wh, alpha, Vector(0,0,1));//local
		//para obtener el costhetaMax apropiado, considerar la fuente a la que la direccion apunta
		if(g.x>0)	costhetaMax = cosinethetaMax(sourceid2, x);
		fpdf = solidAngleProb(costhetaMax);
		wg = powerHeuristics(gpdf, fpdf);
	}

	//MIS

	montecarlo = montecarlo + g*wg;

	return montecarlo;
}

#endif //MISSAMPLINGFUNCTIONS_H
