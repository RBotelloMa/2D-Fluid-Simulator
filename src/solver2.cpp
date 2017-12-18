#include "solver.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

void Solver::Init(unsigned N, float dt, float diff, float visc, bool Tadvec, bool Tadvec2, bool Tproject, bool TWalls1, bool TWalls2)
{
	this->dt = dt;
	this->diff = diff;
	this->visc = visc;
	this->N = N;
	this->Tadvec = Tadvec;
	this->Tadvec2 = Tadvec2;
	this->Tproject = Tproject;
	this->TWalls1 = TWalls1;
	this->TWalls2 = TWalls2;
}


void Solver::Tadvecf() { Tadvec = !(Tadvec); }
void Solver::Tadvec2f() { Tadvec2 = !(Tadvec2); }
void Solver::Tprojectf() { Tproject = !(Tproject); };
void Solver::TWalls1f() { TWalls1 = !(TWalls1); };
void Solver::TWalls2f() { TWalls2 = !(TWalls2); };

/*
----------------------------------------------------------------------
free/clear/allocate simulation data
----------------------------------------------------------------------
*/
void Solver::FreeData(void)
{
	//DONE: Libera los buffers de memoria.
	free(u_prev);
	free(v_prev);
	free(dens_prev);
	free(u);
	free(v);
	free(dens);

}

void Solver::ClearData(void)
{
	//DONE: Borra todo el contenido de los buffers
	for (int i = 0; i < (N+2)*(N+2); i++) {
		u_prev[i] = 0;
		v_prev[i] = 0;
		dens_prev[i] = 0;
		u[i] = 0;
		v[i] = 0;
		dens[i] = 0;
		
	}
}

bool Solver::AllocateData(void)
{
	//DONE:
	//Reservamos memoria, en caso de fallo devlvemos false.
	int _size = (2 + N)*(2 + N)*sizeof(float);

	float * u_prev2;
	float *v_prev2;
	float *dens_prev2;
	float *u2;
	float *v2;
	float *dens2;

	u_prev2 = (float*)malloc(_size);
	v_prev2 = (float*)malloc(_size);
	dens_prev2 = (float*)malloc(_size);
	u2 = (float*)malloc(_size);
	v2 = (float*)malloc(_size);
	dens2 = (float*)malloc(_size);

	this->u_prev = u_prev2;
	this->v_prev = v_prev2;
	this->dens_prev = dens_prev2;
	this->u = u2;
	this->v = v2;
	this->dens = dens2;


	if ((u_prev == NULL)||(v_prev == NULL)||(dens_prev == NULL)||(u == NULL)||(v == NULL)||(dens == NULL)) exit(1);
	
	ClearData();
	return true;
}

void Solver::ClearPrevData() 
{
	//DONE: Borra el contenido de los buffers _prev
	for (int i = 0; i < (N + 2)*(N + 2); i++) {
		u_prev[i] = 0;
		v_prev[i] = 0;
		dens_prev[i] = 0;
	};
}

void Solver::AddDensity(unsigned x, unsigned y, float source)
{
	//DONE: Añade el valor de source al array de densidades. Sería interesante usar la macro: XY_TO_ARRAY
	dens_prev[XY_TO_ARRAY(x, y)] += source;



}

void Solver::AddVelocity(unsigned x, unsigned y, float forceX, float forceY)
{
	//DONE: Añade el valor de fuerza a sus respectivos arrays. Sería interesante usar la macro: XY_TO_ARRAY
	u_prev[XY_TO_ARRAY(x, y)] += forceX;
	v_prev[XY_TO_ARRAY(x, y)] += forceY;
}

void Solver::Solve()
{
	VelStep();
	DensStep();
	/*float sum0 = 0;
	float sum1 = 0;
	float sum2 = 0;
	int i;
	int j;
	FOR_EACH_CELL
		sum0 += dens[XY_TO_ARRAY(i, j)];
	sum1 += u[XY_TO_ARRAY(i, j)];
	sum2 += v[XY_TO_ARRAY(i, j)];
	END_FOR
	printf("%f, %f, %f \n", sum0, sum1, sum2);*/
	
	
}

void Solver::DensStep()
{
	AddSource(dens, dens_prev);			//Adding input density (dens_prev) to final density (dens).
	SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	Diffuse(0, dens, dens_prev);		//Writing result in dens because we made the swap before. bi = dens_prev. The initial trash in dens matrix, doesnt matter, because it converges anyways.
	//SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	if (Tadvec) Advect(0, dens, dens_prev, u, v);	//Advect phase, result in dens.
	//SWAP(dens_prev, dens);
}

void Solver::VelStep()
{
	AddSource(u, u_prev);
	AddSource(v, v_prev);
	SWAP (u_prev,u)			
	SWAP (v_prev, v)
	Diffuse(1, u, u_prev);  
	Diffuse(2, v, v_prev); 
	if (Tproject) Project(u, v, u_prev, v_prev);		//Mass conserving.
	//SWAP (u_prev,u)			
	//SWAP (v_prev,v)
	if (Tadvec2) Advect(1, u, u_prev, u, v);
	if (Tadvec2) Advect(2, v, v_prev, u, v);
	if (Tproject) Project(u, v, u_prev, v_prev);		//Mass conserving.
}

void Solver::AddSource(float * base, float * source)
{
	//DONE: Teniendo en cuenta dt (Delta Time), incrementar el array base con nuestro source. Esto sirve tanto para añadir las nuevas densidades como las nuevas fuerzas.
	int _size = (2 + N)*(2 + N);
	for (int i = 0; i < _size; i++) {
		base[i] += source[i];
	}
}


void Solver::SetBounds(int b, float * x)
{
	/*DONE:
	Input b: 0, 1 or 2.
	0: borders = same value than the inner value.
	1: x axis borders inverted, y axis equal.
	2: y axis borders inverted, x axis equal.
	Corner values allways are mean value between associated edges.
	*/
	if (TWalls1) {
		for (int i = 1; i <= N; i++) {
			x[XY_TO_ARRAY(i, 0)] = x[XY_TO_ARRAY(i, N)];
			x[XY_TO_ARRAY(i, N + 1)] = x[XY_TO_ARRAY(i, 1)];
			x[XY_TO_ARRAY(N + 1, i)] = x[XY_TO_ARRAY( 1, i)];
			x[XY_TO_ARRAY(0, i)] = x[XY_TO_ARRAY(N, i)];
		}
		x[XY_TO_ARRAY(N + 1, N + 1)] = x[XY_TO_ARRAY(1, 1)];
		x[XY_TO_ARRAY(0, 0)] = x[XY_TO_ARRAY(N,N)];
		x[XY_TO_ARRAY(0, N + 1)] = x[XY_TO_ARRAY(N, 1)];
		x[XY_TO_ARRAY(N + 1, 0)] = x[XY_TO_ARRAY(1, N)];
		
	}
	else if (b == 1) {
		for (unsigned i = 1; i <= N; i++) {
			x[XY_TO_ARRAY(i, 0)] = x[XY_TO_ARRAY(i, 1)];//y normal x opuesto
			x[XY_TO_ARRAY(i, N + 1)] = x[XY_TO_ARRAY(i, N)];
			x[XY_TO_ARRAY(N + 1, i)] = -x[XY_TO_ARRAY(N, i)];
			x[XY_TO_ARRAY(0, i)] = -x[XY_TO_ARRAY(1, i)];
		}
		x[XY_TO_ARRAY(0, 0)] = (x[XY_TO_ARRAY(1, 0)] + x[XY_TO_ARRAY(0, 1)]) / 2;
		x[XY_TO_ARRAY(0, N + 1)] = (x[XY_TO_ARRAY(1, N + 1)] + x[XY_TO_ARRAY(0, N)]) / 2;
		x[XY_TO_ARRAY(N + 1, 0)] = (x[XY_TO_ARRAY(N, 0)] + x[XY_TO_ARRAY(N + 1, 1)]) / 2;
		x[XY_TO_ARRAY(N + 1, N + 1)] = (x[XY_TO_ARRAY(N, N + 1)] + x[XY_TO_ARRAY(N + 1, N)]) / 2;
	}
	else if (b == 2) {
		for (unsigned i = 1; i <= N; i++) {
			x[XY_TO_ARRAY(i, 0)] = -x[XY_TO_ARRAY(i, 1)];//y opuesto x normal
			x[XY_TO_ARRAY(i, N + 1)] = -x[XY_TO_ARRAY(i, N)];
			x[XY_TO_ARRAY(N + 1, i)] = x[XY_TO_ARRAY(N, i)];
			x[XY_TO_ARRAY(0, i)] = x[XY_TO_ARRAY(1, i)];
		}
		x[XY_TO_ARRAY(0, 0)] = (x[XY_TO_ARRAY(1, 0)] + x[XY_TO_ARRAY(0, 1)]) / 2;
		x[XY_TO_ARRAY(0, N + 1)] = (x[XY_TO_ARRAY(1, N + 1)] + x[XY_TO_ARRAY(0, N)]) / 2;
		x[XY_TO_ARRAY(N + 1, 0)] = (x[XY_TO_ARRAY(N, 0)] + x[XY_TO_ARRAY(N + 1, 1)]) / 2;
		x[XY_TO_ARRAY(N + 1, N + 1)] = (x[XY_TO_ARRAY(N, N + 1)] + x[XY_TO_ARRAY(N + 1, N)]) / 2;
	}
	else if (b == 0) {
		for (int i = 1; i <= N; i++) {
			x[XY_TO_ARRAY(i, 0)] = x[XY_TO_ARRAY(i, 1)];
			x[XY_TO_ARRAY(i, N + 1)] = x[XY_TO_ARRAY(i, N)];
			x[XY_TO_ARRAY(N + 1, i)] = x[XY_TO_ARRAY(N, i)];
			x[XY_TO_ARRAY(0, i)] = x[XY_TO_ARRAY(1, i)];
		}
		x[XY_TO_ARRAY(0, 0)] = (x[XY_TO_ARRAY(1, 0)] + x[XY_TO_ARRAY(0, 1)]) / 2;
		x[XY_TO_ARRAY(0, N + 1)] = (x[XY_TO_ARRAY(1, N + 1)] + x[XY_TO_ARRAY(0, N)]) / 2;
		x[XY_TO_ARRAY(N + 1, 0)] = (x[XY_TO_ARRAY(N, 0)] + x[XY_TO_ARRAY(N + 1, 1)]) / 2;
		x[XY_TO_ARRAY(N + 1, N + 1)] = (x[XY_TO_ARRAY(N, N + 1)] + x[XY_TO_ARRAY(N + 1, N)]) / 2;
	}
	else exit(1);
}


/*
https://www.youtube.com/watch?v=62_RUX_hrT4
https://es.wikipedia.org/wiki/M%C3%A9todo_de_Gauss-Seidel <- Solución de valores independientes.
Despreciando posibles valores de x no contiguos, se simplifica mucho. Mirar diapositivas y la solución de Gauss Seidel de términos independientes.
Gauss Seidel -> Matrix x and x0
*/

void Solver::LinSolve(int b, float * x, float * x0, float aij, float aii)
{
	//DONE: Se recomienda usar FOR_EACH_CELL, END_FOR y XY_TO_ARRAY.
	//En paso inicial hacer x=x0; 
	int i;
	int j;
	for (i = 0; i <= N + 1; i++) {
		for (j = 0; j <= N + 1; j++)
		{
			x[(N + 2)*i + j] = x0[(N + 2)*i + j];
		}
	}
	for (int r = 0; r < 10; r++) {
		FOR_EACH_CELL
			x[XY_TO_ARRAY2(j, i)] = (1 / aii)*(-aij*(x[XY_TO_ARRAY2(j - 1, i)] + x[XY_TO_ARRAY2(j, i - 1)]) - aij*(x[XY_TO_ARRAY2(j, i + 1)] + x[XY_TO_ARRAY2(j + 1, i)]) + x0[XY_TO_ARRAY2(j, i)]);
		END_FOR	
	}
}

/*
Nuestra función de difusión solo debe resolver el sistema de ecuaciones simplificado a las celdas contiguas de la casilla que queremos resolver,
por lo que solo con la entrada de dos valores, debemos poder obtener el resultado.
*/
void Solver::Diffuse(int b, float * x, float * x0)
{
//DONE: Solo necesitaremos pasar dos parámetros a nuestro resolutor de sistemas de ecuaciones de Gauss Seidel. Calculamos dichos valores y llamamos a la resolución del sistema.
	float aij = - diff* (N)*(N)*dt;
	float aii = 1 + diff*(N)*(N)* 4 *dt;
	LinSolve(b, x, x0, aij, aii);
	SetBounds(b, x);
}

/*
d is overwrited with the initial d0 data and affected by the u & v vectorfield.
Hay que tener en cuenta que el centro de las casillas representa la posición entera dentro de la casilla, por lo que los bordes estan
en las posiciones x,5.
*/
//
void Solver::Advect(int b, float * d, float * d0, float * u, float * v)
{
//TODO: Se aplica el campo vectorial realizando una interploación lineal entre las 4 casillas más cercanas donde caiga el nuevo valor.
	int i, j, p, q, s1, s2;
	float x, y, a, c, f1, f2, f3, f4, sum0;
	for (j = 0; j < (2 + N)*(2 + N); j++) {
		d0[j] = 0;
	}
	if (!TWalls2) {
			for (j = 1; j < N + 1; j++) {
				for (i = 1; i < N + 1; i++)
				{
					x = (float)i + u[XY_TO_ARRAY(j, i)] * dt;
					y = (float)j + v[XY_TO_ARRAY(j, i)] * dt;
					p = round(x);
					q = round(y);
		
					if (!(p<2 || q<2 || p>N - 1 || q>N - 1)) {
						a = p - x;
						c = q - y;
						s1 = ((a > 0) - (a < 0));
						s2 = ((c > 0) - (c < 0));
						f1 = 0.1f*fabs(a)*fabs(c);
						f2 = 0.1f*fabs(1 - a)*fabs(c);
						f3 = 0.1f*fabs(a)*fabs(1 - c);
						f4 = 0.1f*fabs(1 - a)*fabs(1 - c);
		
						sum0 = f1*d[XY_TO_ARRAY(q, p)];
						d0[XY_TO_ARRAY(j, i)] += sum0;
						d0[XY_TO_ARRAY(q, p)] -= sum0;
		
						sum0 = f2*d[XY_TO_ARRAY(q, p + s1)];
						d0[XY_TO_ARRAY(j, i)] += sum0;
						d0[XY_TO_ARRAY(q, p + s1)] -= sum0;
		
						sum0 = f3*d[XY_TO_ARRAY(q + s2, p)];
						d0[XY_TO_ARRAY(j, i)] += sum0;
						d0[XY_TO_ARRAY(q + s2, p)] -= sum0;
		
						sum0 = f4*d[XY_TO_ARRAY(q + s2, p + s1)];
						d0[XY_TO_ARRAY(j, i)] += sum0;
						d0[XY_TO_ARRAY(q + s2, p + s1)] -= sum0;	
					}
				}
			}
			for (j = 0; j < (2 + N)*(2 + N); j++) {
						d[j] += d0[j];
				}
		}

	else{
		for (j = 1; j < N + 1; j++) {
			for (i = 1; i < N + 1; i++)
			{

				x = (float)i - u[XY_TO_ARRAY2(j, i)] *N/N* dt;
				y = (float)j - v[XY_TO_ARRAY2(j, i)] *N/N* dt;
				p = round(x);
				q = round(y);
				a = -p + x;
				c = -q + y;


				if (!(p<1 || q<1 || p>N || q>N)) {

					s1 = ((a > 0) - (a < 0));
					s2 = ((c > 0) - (c < 0));
					f1 = 0.2f*fabs(a)*fabs(c);
					f2 = 0.2f*(1 - fabs(a))*fabs(c);
					f3 = 0.2f*fabs(a)*(1 - fabs(c));
					f4 = 0.2f*(1 - fabs(a))*(1 - fabs(c));
					s1 = p + s1;
					s2 = q + s2;
					if (s1 >= (N + 1)) { s1 = 1; }
					if (s2 >= (N + 1)) { s2 = 1; }
					if (s1 <= 0) { s1 = N; }
					if (s2 <= 0) { s2 = N; }

					sum0 = f1*d[XY_TO_ARRAY2(q, p)];
					d0[XY_TO_ARRAY2(j, i)] += sum0;
					d0[XY_TO_ARRAY2(q, p)] -= sum0;

					sum0 = f2*d[XY_TO_ARRAY2(q, s1)];
					d0[XY_TO_ARRAY2(j, i)] += sum0;
					d0[XY_TO_ARRAY2(q, s1)] -= sum0;

					sum0 = f3*d[XY_TO_ARRAY2(s2, p)];
					d0[XY_TO_ARRAY2(j, i)] += sum0;
					d0[XY_TO_ARRAY2(s2, p)] -= sum0;

					sum0 = f4*d[XY_TO_ARRAY2(s2, s1)];
					d0[XY_TO_ARRAY2(j, i)] += sum0;
					d0[XY_TO_ARRAY2(s2, s1)] -= sum0;

				}
			}
		}

		for (j = 0; j < (2 + N)*(2 + N); j++) {
			d[j] += d0[j];
		}
	}
}

//Advec no cíclico

/*
Se encarga de estabilizar el fluido y hacerlo conservativo de masa. Se usa solo en las matrices de velocidades.
No necesaria implementación por su complejidad.
*/
void Solver::Project(float * u, float * v, float * p, float * div)
{
	int i, j;

	FOR_EACH_CELL
		div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
		p[XY_TO_ARRAY(i, j)] = 0;
	END_FOR
	SetBounds(0, div);
	SetBounds(0, p);

	LinSolve(0, p, div, 1, 4);

	//Aproximamos: Laplaciano de q a su gradiente.
	FOR_EACH_CELL
		u[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i + 1, j)] - p[XY_TO_ARRAY(i - 1, j)]);
		v[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i, j + 1)] - p[XY_TO_ARRAY(i, j - 1)]);
	END_FOR
	SetBounds(1, u);
	SetBounds(2, v);
}