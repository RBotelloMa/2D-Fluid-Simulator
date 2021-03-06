#ifndef _SOLVER_H_
#define _SOLVER_H_

#define XY_TO_ARRAY(i,j) ((i)+(N+2)*(j))
#define XY_TO_ARRAY2(j,i) ((i)+(N+2)*(j))
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define FOR_EACH_CELL2 for ( i=2 ; i<N ; i++ ) { for ( j=2 ; j<N ; j++ ) {
#define END_FOR }}
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define MODY(x) (((x-1+N)%N)+1)

class Solver
{
	float dt, diff, visc;
	int N;
	float * u_prev, *v_prev, *dens_prev;
	bool Tadvec, Tadvec2, Tproject, TWalls1, TWalls2;

public:
	float * u, * v, * dens;
	void Init(unsigned N, float dt, float diff, float visc, bool Tadvec, bool Tadvec2, bool Tproject, bool TWalls1, bool TWalls2);
	void FreeData(void);
	void ClearData(void);
	bool AllocateData(void);
	void ClearPrevData(void);
	void AddDensity(unsigned i, unsigned j, float source);
	void AddVelocity(unsigned i, unsigned j, float forceX, float forceY);
	void Solve(void);
	void Tadvecf();
	void Tadvec2f();
	void Tprojectf();
	void TWalls1f();
	void TWalls2f();


private:
	void DensStep(void);
	void VelStep(void);

	void AddSource(float * x, float * s);
	void SetBounds(int b, float * x);
	void LinSolve( int b, float * x, float * x0, float a, float c);
	void Diffuse(int b, float * x, float * x0);
	void Advect(int b, float * d, float * d0, float * u, float * v);
	void Project(float * u, float * v, float * p, float * div);
};

#endif
