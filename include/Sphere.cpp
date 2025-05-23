//
// Created by Gabriel Castillo on 28/01/25.
//

#include "Sphere.h"

std::vector<Sphere> spheres = {

	//Escena: radio, posicion, color, radiancia, material, eta, kappa, alpha
	//paredes
	Sphere(1e5,  Point(-1e5 - 49, 0, 0),  Color(.5, .5, .5),    Color (0, 0, 0),0, Color(), Color(),0), // left wall
  	Sphere(1e5,  Point(1e5 + 49, 0, 0),  Color(.0, .0, .5),    Color (0, 0, 0),0, Color(), Color(),0), // right wall
  	Sphere(1e5,  Point(0, 0, -1e5 - 81.6),  Color(.5, .5, .5),    Color (0, 0, 0),0, Color(), Color(),0), // back wall
  	Sphere(1e5,  Point(0, -1e5 - 40.8, 0),  Color(.5, .5, .5),    Color (0, 0, 0),0, Color(), Color(),0), // floor
  	Sphere(1e5,  Point(0, 1e5 + 40.8, 0),  Color(.5, .5, .5),     Color (0, 0, 0),0, Color(), Color(),0), // ceiling

	Sphere(16.5, Point(-23, -24.3, -34.6), Color(), Color(), 1, Color(1.66058, 0.88143, 0.521467), Color(9.2282, 6.27077, 4.83803), 0.09), // esfera abajo-izq
	Sphere(16.5,  Point(23, -24.3, -3.6),  Color(.0, .0, .9),    Color (0, 0, 0),0, Color(), Color(),0), // bottom right sphere
	//fuentes
	Sphere(2,  Point(0, 24.3, -35),  Color(),    Color (100, 100, 0),0, Color(), Color(),0),// light source
	Sphere(0,  Point(-23, 24.3, 0),  Color(),    Color (6000, 0, 0),0, Color(), Color(),0), // light source
	Sphere(2,  Point(23, 24.3, 35),  Color(),    Color (75, 75, 60),0, Color(), Color(),0) // light source




 /*
 //ESCENA DOS - COMPARAR LOS EFECTOS DE SUBIR Y BAJAR LOS VALORES DE SIGMA
	//paredes
	Sphere(1e5,  Point(-1e5 - 49, 0, 0),  Color(), Color(), 1, Color(1.66058, 0.88143, 0.521467), Color(9.2282, 6.27077, 4.83803), 0.03), // left wall
  	Sphere(1e5,  Point(1e5 + 49, 0, 0),  Color(), Color(), 1, Color(1.66058, 0.88143, 0.521467), Color(9.2282, 6.27077, 4.83803), 0.03), // right wall
  	Sphere(1e5,  Point(0, 0, -1e5 - 81.6),  Color(.25, .75, .25),    Color (0, 0, 0),0, Color(), Color(),0), // back wall
  	Sphere(1e5,  Point(0, -1e5 - 40.8, 0),  Color(.25, .75, .75),    Color (0, 0, 0),0, Color(), Color(),0), // floor
  	Sphere(1e5,  Point(0, 1e5 + 40.8, 0),  Color(.75, .75, .25),     Color (0, 0, 0),0, Color(), Color(),0), // ceiling

	//esferas dentro
	Sphere(16.5, Point(-23, -24.3, -34.6), Color(.75, .75, .25), Color(), 0, Color(), Color(), 0), // esfera abajo-izq
	Sphere(16.5,  Point(23, -24.3, -3.6),  Color(.4, .3, .2),    Color (0, 0, 0),0, Color(), Color(),0), // bottom right sphere


	//Sphere(0,  Point(24, 24.3, -50),  Color(),   Color(0, 16000, 16000),0, Color(), Color(),0), // light source 1
	Sphere(0,  Point(14, -24.3, -35),  Color(),   Color(2000, 2000, 3000),0, Color(), Color(),0), // light source 2
	//Sphere(0,  Point(-23, -24.3, 10),  Color(),   Color(4000, 4000, 2000),0, Color(), Color(),0), // light source 2

	//Sphere(0, Point(0, -24.3, -30), Color(),   Color(8000, 8000, 0),0, Color(), Color(),0) //light source 3

*/
/*
	//ESCENA 3 - QUE SUCEDE CUANDO LA FUENTE ESTA MAS CERCA DE LA CÁMARA
	Sphere(30, Point(0, 11.2, 165), Color(.0, .25, .75), Color(), 0, Color(), Color(), 0), // esfera abajo-izq
	//Sphere(16.5, Point(23, -24.3, -3.6),   Color(), Color(), 1, Color(0.143245, 0.377423, 1.43919), Color(3.98479, 2.3847, 1.60434), 0.3), // esfera abajo-der
	//Sphere(100, Point(0, -24.3, -200),   Color(), Color(), 1, Color(0.143245, 0.377423, 1.43919), Color(3.98479, 2.3847, 1.60434), 0.02), // esfera abajo-der

	//esferas al rededor
	Sphere(16.5, Point(0, -10, 200), Color(.75, .75, .75), Color(), 0, Color(), Color(), 0), // esfera abajo-izq


	Sphere(0,  Point(0, 11.2, 204),  Color(),   Color(400, 400,400),0, Color(), Color(),0), // light source 1 izq arriba

	//Sphere(0,  Point(-24, 10, -34.6),  Color(),   Color(2000, 5000, 1000),0, Color(), Color(),0), // light source 2 der arriba
	//Sphere(0, Point(0, -24.3, -30), Color(),   Color(4000, 8000, 4000),0, Color(), Color(),0) //light source 3 abajo al centro

*/
	/*
	//fuentes de area que tienden a puntuales
	Sphere(1e5,  Point(-1e5 - 49, 0, 0),  Color(.75, .25, .25),    Color (0, 0, 0),0, Color(), Color(),0), // left wall
	Sphere(1e5,  Point(1e5 + 49, 0, 0),  Color(.25, .25, .75),    Color (0, 0, 0),0, Color(), Color(),0), // right wall
	Sphere(1e5,  Point(0, 0, -1e5 - 81.6),  Color(.25, .75, .25),    Color (0, 0, 0),0, Color(), Color(),0), // back wall
	Sphere(1e5,  Point(0, -1e5 - 40.8, 0),  Color(.25, .75, .75),    Color (0, 0, 0),0, Color(), Color(),0), // floorSphere(1e5,  Point(0, 1e5 + 40.8, 0),  Color(.75, .75, .25),     Color (0, 0, 0),0, Color(), Color(),0), // ceiling
	//materiales
	Sphere(16.5, Point(-23, -24.3, -34.6), Color(), Color(), 1, Color(1.66058, 0.88143, 0.521467), Color(9.2282, 6.27077, 4.83803), 0.03), // esfera abajo-izq

	Sphere(12,  Point(24, 24.3, -50),  Color(),   Color(0, 800, 800),0, Color(), Color(),0), // light source 1
	//Sphere(6,  Point(-24, 24.3, -50),  Color(),   Color(800, 0, 400),0, Color(), Color(),0), // light source 2
	//Sphere(6, Point(0, -24.3, -30), Color(),   Color(800, 800, 0),0, Color(), Color(),0) //light source 3

	*/
	/*
	 *Escena  1 primitive infinite - NO MODIFICAR
	//escena especial, no paredes, unicamente fuentes y esferas flotando en el espacio
	Sphere(16.5, Point(-23, -24.3, -34.6), Color(), Color(), 1, Color(1.66058, 0.88143, 0.521467), Color(9.2282, 6.27077, 4.83803), 0.03), // esfera abajo-izq
	Sphere(16.5, Point(23, -24.3, -3.6),   Color(), Color(), 1, Color(0.143245, 0.377423, 1.43919), Color(3.98479, 2.3847, 1.60434), 0.3), // esfera abajo-der
	Sphere(100, Point(0, -24.3, -200),   Color(), Color(), 1, Color(0.143245, 0.377423, 1.43919), Color(3.98479, 2.3847, 1.60434), 0.02), // esfera abajo-der


	Sphere(0,  Point(24, 24.3, -3.6),  Color(),   Color(2000, 2000, 2000),0, Color(), Color(),0), // light source 1 izq arriba
	Sphere(0,  Point(-24, 10, -34.6),  Color(),   Color(2000, 5000, 1000),0, Color(), Color(),0), // light source 2 der arriba
	Sphere(0, Point(0, -24.3, -30), Color(),   Color(4000, 8000, 4000),0, Color(), Color(),0) //light source 3 abajo al centro
*/

	/*
	//Escena: radio, posicion, color
	Sphere(1e5,  Point(-1e5 - 49, 0, 0),  Color(.5, .5, .5), Color(), 0, Color(), Color(), 0), // pared izq
	Sphere(1e5,  Point(1e5 + 49, 0, 0),    Color(.5, .5, .5), Color(), 0, Color(), Color(), 0), // pared der
	Sphere(1e5,  Point(0, 0, -1e5 - 81.6), Color(.5, .5, .5), Color(), 0, Color(), Color(), 0), // pared detras
	Sphere(1e5,  Point(0, -1e5 - 40.8, 0),Color(.5, .5, .5), Color(), 0, Color(), Color(), 0), // suelo
	Sphere(1e5,  Point(0, 1e5 + 40.8, 0), Color(.5, .5, .5), Color(), 0, Color(), Color(), 0), // techo
	//Sphere(0, Point(0,24.3,-10), Color(1, 1, 1), Color(8000,8000,0),0, Color(), Color(),0), // esfera arriba
	//Sphere(16.5, Point(-23, -24.3, -34.6), Color(), Color(), 0, Color(), Color(), 0),
	//Sphere(16.5, Point(-23, -24.3, -34.6), Color(), Color(), 1, Color(1.66058, 0.88143, 0.521467), Color(9.2282, 6.27077, 4.83803), 0.03), // esfera abajo-izq
	//Sphere(16.5, Point(23, -24.3, -3.6),   Color(), Color(), 1, Color(0.143245, 0.377423, 1.43919), Color(3.98479, 2.3847, 1.60434), 0.3), // esfera abajo-der
	//Sphere(16.5, Point(-23, -24.3, -20.6), Color(.75, .75, .25), Color(), 3, Color(), Color(), 0),
	Sphere(16.5, Point(23, -24.3, -3.6), Color(.50, .50, 0), Color(), 0, Color(), Color(), 0),
	Sphere(0, Point(-23, 0, -10.6), Color(1, 1, 1), Color(6000,6000,6000),0, Color(), Color(), 0),
	Sphere(0, Point(23, 24.3, -50), Color(1, 1, 1), Color(4000,4000,4000),0, Color(), Color(), 0)
*/
};