//
// Created by Gabriel Castillo on 28/01/25.
//

#include "Sphere.h"

std::vector<Sphere> spheres = {
	/*
	//Escena: radio, posicion, color, radiancia, material, eta, kappa, alpha
	Sphere(1e5,  Point(-1e5 - 49, 0, 0),  Color(.75, .25, .25),    Color (0, 0, 0),0, Color(), Color(),0), // left wall
  	Sphere(1e5,  Point(1e5 + 49, 0, 0),  Color(.25, .25, .75),    Color (0, 0, 0),0, Color(), Color(),0), // right wall
  	Sphere(1e5,  Point(0, 0, -1e5 - 81.6),  Color(.25, .75, .25),    Color (0, 0, 0),0, Color(), Color(),0), // back wall
  	Sphere(1e5,  Point(0, -1e5 - 40.8, 0),  Color(.25, .75, .75),    Color (0, 0, 0),0, Color(), Color(),0), // floor
  	Sphere(1e5,  Point(0, 1e5 + 40.8, 0),  Color(.75, .75, .25),     Color (0, 0, 0),0, Color(), Color(),0), // ceiling
  	Sphere(16.5,  Point(-23, -24.3, -34.6),  Color(.2, .3, .4),    Color (0, 0, 0),0, Color(), Color(),0), // bottom left sphere
  	Sphere(16.5,  Point(23, -24.3, -3.6),  Color(.4, .3, .2),    Color (0, 0, 0),0, Color(), Color(),0), // bottom right sphere
  	Sphere(5,  Point(14, -24.3, -35),  Color(),    Color (12, 12, 12),0, Color(), Color(),0) // light source
*/

/*
 //escena para vpt
	Sphere(1e5,  Point(-1e5 - 49, 0, 0),  Color(.75, .25, .25),    Color (0, 0, 0),0, Color(), Color(),0), // left wall
  	Sphere(1e5,  Point(1e5 + 49, 0, 0),  Color(.25, .25, .75),    Color (0, 0, 0),0, Color(), Color(),0), // right wall
  	Sphere(1e5,  Point(0, 0, -1e5 - 81.6),  Color(.25, .75, .25),    Color (0, 0, 0),0, Color(), Color(),0), // back wall
  	Sphere(1e5,  Point(0, -1e5 - 40.8, 0),  Color(.25, .75, .75),    Color (0, 0, 0),0, Color(), Color(),0), // floor
  	Sphere(1e5,  Point(0, 1e5 + 40.8, 0),  Color(.75, .75, .25),     Color (0, 0, 0),0, Color(), Color(),0), // ceiling
	//pruebas
	Sphere(16.5, Point(-23, -24.3, -34.6), Color(), Color(), 1, Color(1.66058, 0.88143, 0.521467), Color(9.2282, 6.27077, 4.83803), 0.03), // esfera abajo-izq

	Sphere(0,  Point(24, 24.3, -50),  Color(),   Color(0, 16000, 16000),0, Color(), Color(),0), // light source 1
	Sphere(0,  Point(-24, 24.3, -50),  Color(),   Color(8000, 0, 8000),0, Color(), Color(),0), // light source 2
	Sphere(0, Point(0, -24.3, -30), Color(),   Color(8000, 8000, 0),0, Color(), Color(),0) //light source 3
*/


/*
	//escena para vpt paredes negras solo fuentes
	Sphere(1e5,  Point(-1e5 - 49, 0, 0),  Color(),    Color (0, 0, 0),0, Color(), Color(),0), // left wall
	Sphere(1e5,  Point(1e5 + 49, 0, 0),  Color(),    Color (0, 0, 0),0, Color(), Color(),0), // right wall
	Sphere(1e5,  Point(0, 0, -1e5 - 81.6),  Color(),    Color (0, 0, 0),0, Color(), Color(),0), // back wall
	Sphere(1e5,  Point(0, -1e5 - 40.8, 0),  Color(),    Color (0, 0, 0),0, Color(), Color(),0), // floor
	Sphere(1e5,  Point(0, 1e5 + 40.8, 0),  Color(),     Color (0, 0, 0),0, Color(), Color(),0), // ceiling
	//pruebas con un material microfacet
	//Sphere(16.5, Point(-23, -24.3, -34.6), Color(), Color(), 1, Color(1.66058, 0.88143, 0.521467), Color(9.2282, 6.27077, 4.83803), 0.03), // esfera abajo-izq

	Sphere(0,  Point(24, 24.3, -50),  Color(),   Color(0, 8000, 8000),0, Color(), Color(),0), // light source 1
	//Sphere(0,  Point(-24, 24.3, -50),  Color(),   Color(8000, 0, 8000),0, Color(), Color(),0), // light source 2
	//Sphere(0, Point(0, -24.3, -30), Color(),   Color(8000, 8000, 0),0, Color(), Color(),0) //light source 3
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
	//escena especial, no paredes, unicamente fuentes

	Sphere(12,  Point(24, 24.3, -50),  Color(),   Color(0, 800, 800),0, Color(), Color(),0), // light source 1
	//Sphere(6,  Point(-24, 24.3, -50),  Color(),   Color(400, 0, 800),0, Color(), Color(),0), // light source 2
	//Sphere(6, Point(0, -24.3, -30), Color(),   Color(800, 800, 0),0, Color(), Color(),0) //light source 3
/*
	//Escena: radio, posicion, color
	Sphere(1e5,  Point(-1e5 - 49, 0, 0),  Color(.5, .5, .5), Color(), 0, Color(), Color(), 0), // pared izq
	Sphere(1e5,  Point(1e5 + 49, 0, 0),    Color(.5, .5, .5), Color(), 0, Color(), Color(), 0), // pared der
	Sphere(1e5,  Point(0, 0, -1e5 - 81.6), Color(.5, .5, .5), Color(), 0, Color(), Color(), 0), // pared detras
	Sphere(1e5,  Point(0, -1e5 - 40.8, 0),Color(.5, .5, .5), Color(), 0, Color(), Color(), 0), // suelo
	Sphere(1e5,  Point(0, 1e5 + 40.8, 0), Color(.5, .5, .5), Color(), 0, Color(), Color(), 0), // techo
	Sphere(0, Point(0,24.3,-10), Color(1, 1, 1), Color(8000,8000,0),0, Color(), Color(),0), // esfera arriba
	//Sphere(16.5, Point(-23, -24.3, -34.6), Color(), Color(), 0, Color(), Color(), 0),
	//Sphere(16.5, Point(-23, -24.3, -34.6), Color(), Color(), 1, Color(1.66058, 0.88143, 0.521467), Color(9.2282, 6.27077, 4.83803), 0.03), // esfera abajo-izq
	//Sphere(16.5, Point(23, -24.3, -3.6),   Color(), Color(), 1, Color(0.143245, 0.377423, 1.43919), Color(3.98479, 2.3847, 1.60434), 0.3), // esfera abajo-der
	//Sphere(16.5, Point(-23, -24.3, -20.6), Color(.75, .75, .25), Color(), 3, Color(), Color(), 0),
	Sphere(16.5, Point(23, -24.3, -3.6), Color(.70, .3, 0), Color(), 0, Color(), Color(), 0),
	Sphere(0, Point(-23, 0, -10.6), Color(1, 1, 1), Color(0,8000,8000),0, Color(), Color(), 0),
	Sphere(0, Point(23, 24.3, -50), Color(1, 1, 1), Color(8000,0,8000),0, Color(), Color(), 0)
*/
};