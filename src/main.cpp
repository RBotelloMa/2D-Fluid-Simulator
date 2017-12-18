#include <stdlib.h>
#include <stdio.h>
#include <glut.h>
#include <cmath> 
#include "solver.h"

#define LIGHTCG(color) if (Toggel2==0) color = sqrtf(sqrtf(color*color*color));
#define LIGHTCG2(color) if (Toggel2==1) color = sqrtf(color);


/* global variables */
static int N;
static float dt, diff, visc;
static float force, source;
static bool dvel;

static int win_id;
static int win_x, win_y;	//Window Size.
static int mouse_down[3];	//Mouse button states.
static int omx, omy, mx, my;

static Solver solver;

const float OMG = 2*3.1415/3;
const float cOMG = cosf(OMG);
const float sOMG = sinf(OMG);
bool Toggel = true;
int Toggel2 = 0;

bool Tadvec=true, Tadvec2 = true, Tproject=true, TWalls1= false, TWalls2= false;



/*
----------------------------------------------------------------------
OpenGL specific drawing routines
----------------------------------------------------------------------
*/
static void PreDisplay(void)
{
	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
}

static void PostDisplay(void)
{
	glutSwapBuffers();
}

static void DrawVelocity(void)
{
	//DONE? NOPE
	glBegin(GL_LINES);
	float color1, color2,color3,color4,color5,color6;
	float prop = (1 / (float)(N + 2));
	float posx, posy,norm;
	
	for (unsigned j = 0; j <= N + 1; j++) {
		for (unsigned i = 0; i <= N + 1; i++) {
			posx = solver.u[XY_TO_ARRAY(i, j)];
			posy = solver.v[XY_TO_ARRAY(i, j)];
			norm = 1;// fmax(fabs(posx), fabs(posy));
			color1 = fmax(fmin(posx, 1),0);
			color2 = fmax(fmin((cOMG*posx-sOMG*posy) ,1),0);
			color3 = fmax(fmin((cOMG*posx + sOMG*posy),1),0);
			color4 = fmax(fmin(-posx , 1),0)/2;
			color5 = fmax(fmin(-(cOMG*(posx) - sOMG*posy) , 1),0)/2;
			color6 = fmax(fmin(-(cOMG*posx + sOMG*posy) , 1),0)/2;
			glColor3f(color1+color5+color6, color2+color4+color6,color3+color4+color5);
			glVertex2f((i + 0.5 - posx*norm*0.45)*prop, (j + 0.5 - posy*norm*0.45)*prop);
			glVertex2f((i + 0.5 + posx*norm*0.45)*prop, (j + 0.5 + posy*norm*0.45)*prop);
			
		}
	}
	glEnd();
}


static void DrawDensity(void)
{
	
	if (!Toggel){
		glBegin(GL_QUADS);
		float color;
		float prop = (1 / (float)(N + 2));
		for (int j = 0; j <= N + 1; j++) {
			for (int i = 0; i <= N + 1; i++) {
				color = fmin(solver.dens[XY_TO_ARRAY(i, j)], 1);
				LIGHTCG(color);
				LIGHTCG2(color);
				glColor3f(color, color, color);
				glVertex2f((i)*prop, (j)*prop);
				glVertex2f((i)*prop, (j + 1.0f)*prop);
				glVertex2f((i + 1.0f)*prop, (j + 1.0f)*prop);
				glVertex2f((i + 1.0f)*prop, (j)*prop);

			}
		}
		glEnd();
	}
	else {
		glBegin(GL_QUADS);
		float color, color1, color2, color3, color4, color5, color6, color7, color8, color9;
		float prop = (1 / (float)(N + 2));
		for (int j = 1; j < N + 1; j++) {
			for (int i = 1; i < N + 1; i++) {
				color1 = (fmin(solver.dens[XY_TO_ARRAY(i-1, j-1)], 1));
				color2 = (fmin(solver.dens[XY_TO_ARRAY(i-1, j)], 1));
				color3 = (fmin(solver.dens[XY_TO_ARRAY(i-1, j+1)], 1));
				color4 = (fmin(solver.dens[XY_TO_ARRAY(i, j - 1)], 1));
				color5 = (fmin(solver.dens[XY_TO_ARRAY(i, j)], 1));
				color6 = (fmin(solver.dens[XY_TO_ARRAY(i, j + 1)], 1));
				color7 = (fmin(solver.dens[XY_TO_ARRAY(i+1, j-1)], 1));
				color8 = (fmin(solver.dens[XY_TO_ARRAY(i+1, j)], 1));
				color9 = (fmin(solver.dens[XY_TO_ARRAY(i+1, j+1)], 1));
				
				LIGHTCG(color1) LIGHTCG(color2) LIGHTCG(color3) LIGHTCG(color4) 
				LIGHTCG(color5) LIGHTCG(color6) LIGHTCG(color7) LIGHTCG(color8) LIGHTCG(color9)

				LIGHTCG2(color1) LIGHTCG2(color2) LIGHTCG2(color3) LIGHTCG2(color4)
				LIGHTCG2(color5) LIGHTCG2(color6) LIGHTCG2(color7) LIGHTCG2(color8) LIGHTCG2(color9)

				color = (color1+2*color2 + 2*color4 + 4 * color5) / 9;
				glColor3f(color, color, color);


				glVertex2f((i)*prop, (j)*prop);
				glVertex2f((i)*prop, (j + 0.3333f)*prop);
				glVertex2f((i + 0.3333f)*prop, (j + 0.3333)*prop);
				glVertex2f((i + 0.3333f)*prop, (j)*prop);


				color = (3*color4 + 6 * color5) / 9;
				glColor3f(color, color, color);

				glVertex2f((i + 0.3333f)*prop, (j)*prop);
				glVertex2f((i + 0.3333f)*prop, (j + 0.3333f)*prop);
				glVertex2f((i + 0.6666f)*prop, (j + 0.3333f)*prop);
				glVertex2f((i + 0.6666f)*prop, (j)*prop);

				color = (color7 + 2 * color4 + 2 * color8 + 4 * color5) / 9;
				glColor3f(color, color, color);

				glVertex2f((i + 0.6666f)*prop, (j)*prop);
				glVertex2f((i + 0.6666f)*prop, (j + 0.3333f)*prop);
				glVertex2f((i + 1.0f)*prop, (j + 0.3333f)*prop);
				glVertex2f((i + 1.0f)*prop, (j)*prop);


				color = (3*color2 + 6 * color5) / 9;
				glColor3f(color, color, color);

				glVertex2f((i)*prop, (j + 0.3333f)*prop);
				glVertex2f((i)*prop, (j + 0.6666f)*prop);
				glVertex2f((i + 0.3333f)*prop, (j + 0.6666f)*prop);
				glVertex2f((i + 0.3333f)*prop, (j + 0.3333f)*prop);

				color = color5;
				glColor3f(color, color, color);

				glVertex2f((i + 0.3333f)*prop, (j + 0.3333f)*prop);
				glVertex2f((i + 0.3333f)*prop, (j + 0.6666f)*prop);
				glVertex2f((i + 0.6666f)*prop, (j + 0.6666f)*prop);
				glVertex2f((i + 0.6666f)*prop, (j + 0.3333f)*prop);


				color = (3*color8 + 6 * color5) / 9;
				glColor3f(color, color, color);

				glVertex2f((i + 0.6666f)*prop, (j + 0.3333f)*prop);
				glVertex2f((i + 0.6666f)*prop, (j + 0.6666f)*prop);
				glVertex2f((i + 1.0f)*prop, (j + 0.6666f)*prop);
				glVertex2f((i + 1.0f)*prop, (j + 0.3333f)*prop);

				color = (color3 + 2 * color6 + 2 * color2 + 4 * color5) / 9;
				glColor3f(color, color, color);

				glVertex2f((i)*prop, (j + 0.6666f)*prop);
				glVertex2f((i)*prop, (j + 1.0f)*prop);
				glVertex2f((i + 0.3333f)*prop, (j + 1.0f)*prop);
				glVertex2f((i + 0.3333f)*prop, (j + 0.6666f)*prop);

				color = (3*color6 + 6 * color5) / 9;
				glColor3f(color, color, color);

				glVertex2f((i + 0.3333f)*prop, (j + 0.6666f)*prop);
				glVertex2f((i + 0.3333f)*prop, (j + 1.0f)*prop);
				glVertex2f((i + 0.6666f)*prop, (j + 1.0f)*prop);
				glVertex2f((i + 0.6666f)*prop, (j + 0.6666f)*prop);

				color = (color9 + 2 * color6 + 2 * color8 + 4 * color5) / 9;
				glColor3f(color, color, color);

				glVertex2f((i + 0.6666f)*prop, (j + 0.6666f)*prop);
				glVertex2f((i + 0.6666f)*prop, (j + 1.0f)*prop);
				glVertex2f((i + 1.0f)*prop, (j + 1.0f)*prop);
				glVertex2f((i + 1.0f)*prop, (j + 0.6666f)*prop);



			}
		}
		glEnd();
	}
}

/*
----------------------------------------------------------------------
relates mouse movements to forces and sources
----------------------------------------------------------------------
*/
static void AddInteractionFromUI()
{
	int i, j;

	if (!mouse_down[0] && !mouse_down[2]) return;

	i = (int)((mx / (float)win_x)*N + 1);
	j = (int)(((win_y - my) / (float)win_y)*N + 1);

	if (i<1 || i>N || j<1 || j>N) return;

	if (mouse_down[GLUT_LEFT_BUTTON]) {
		solver.AddVelocity(i, j, force * (mx - omx), force * (omy - my));
	}

	if (mouse_down[GLUT_RIGHT_BUTTON]) {
		solver.AddDensity(i, j, source);
	}
	
	omx = mx;
	omy = my;

	return;
}

/*
----------------------------------------------------------------------
GLUT callback routines
----------------------------------------------------------------------
*/

static void KeyFunc(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'c':
	case 'C':
		solver.ClearData();
		break;

	case 'q':
	case 'Q':
		solver.FreeData();
		exit(0);
		break;
	case 'v':
	case 'V':
		dvel = !dvel;
		break;
	case 't':
	case 'T':
		Toggel = !Toggel;
		break;
	case 'k':
	case 'K':
		Toggel2 += 1;
		Toggel2 = Toggel2 % 3;
		break;
	case 'W':
	case 'w':
		solver.TWalls1f();
		TWalls1 = !TWalls1;
		printf("Wall1 turned %i\n", TWalls1);
		break;
	case 'u':
	case 'U':
		solver.TWalls2f();
		TWalls2 = !TWalls2;
		printf("Wall2 for advec turned %i \n", TWalls2);
		break;
	case 'a':
	case 'A':
		solver.Tadvecf();
		Tadvec = !Tadvec;
		printf("Advect turned %i\n", Tadvec);
		break;
	case 's':
	case 'S':
		solver.Tadvec2f();
		Tadvec2 = !Tadvec2;
		printf("Advect for vectors turned %i\n", Tadvec2);
		break;

	case 'p':
	case 'P':
		solver.Tprojectf();
		Tproject = !Tproject;
		printf("Project turned %i \n", Tproject);
		break;


	}
}

static void MouseFunc(int button, int state, int x, int y)
{
	omx = mx = x;
	omx = my = y;

	mouse_down[button] = state == GLUT_DOWN;
}

static void MotionFunc(int x, int y)
{
	mx = x;
	my = y;
}

static void ReshapeFunc(int width, int height)
{
	glutSetWindow(win_id);
	glutReshapeWindow(width, height);

	win_x = width;
	win_y = height;
}

static void IdleFunc(void)
{
	solver.ClearPrevData(); //Clean last step forces
	AddInteractionFromUI();	//Add Forces and Densities

	solver.Solve();			//Calculate the next step

	glutSetWindow(win_id);
	glutPostRedisplay();
}

static void DisplayFunc(void)
{
	PreDisplay();

	if (dvel) DrawVelocity();
	else		DrawDensity();

	PostDisplay();
}


/*
----------------------------------------------------------------------
open_glut_window --- open a glut compatible window and set callbacks
----------------------------------------------------------------------
*/

static void OpenGlutWindow(void)
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("Alias | wavefront");

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();

	PreDisplay();
	
	glutKeyboardFunc(KeyFunc);
	glutMouseFunc(MouseFunc);
	glutMotionFunc(MotionFunc);
	glutReshapeFunc(ReshapeFunc);
	glutIdleFunc(IdleFunc);
	glutDisplayFunc(DisplayFunc);
}


/*
----------------------------------------------------------------------
main --- main routine
----------------------------------------------------------------------
*/

int main(int argc, char ** argv)
{
	glutInit(&argc, argv);

	if (argc != 1 && argc != 6) {
		fprintf(stderr, "usage : %s N dt diff visc force source\n", argv[0]);
		fprintf(stderr, "where:\n"); \
		fprintf(stderr, "\t N      : grid resolution\n");
		fprintf(stderr, "\t dt     : time step\n");
		fprintf(stderr, "\t diff   : diffusion rate of the density\n");
		fprintf(stderr, "\t visc   : viscosity of the fluid\n");
		fprintf(stderr, "\t force  : scales the mouse movement that generate a force\n");
		fprintf(stderr, "\t source : amount of density that will be deposited\n");
		exit(1);
	}

	if (argc == 1) {
		N = 101;
		dt = 0.1f;
		diff = 0.0001f;
		visc = 0.1f;
		force = 5.0f;
		source = 10.0f;
		fprintf(stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n",
			N, dt, diff, visc, force, source);
	}
	else {
		N = atoi(argv[1]);
		dt = atof(argv[2]);
		diff = atof(argv[3]);
		visc = atof(argv[4]);
		force = atof(argv[5]);
		source = atof(argv[6]);
	}

	printf("\n\nHow to use this demo:\n\n");
	printf("\t Add densities with the right mouse button\n");
	printf("\t Add velocities with the left mouse button and dragging the mouse\n");
	printf("\t Toggle density/velocity display with the 'v' key\n");
	printf("\t Clear the simulation by pressing the 'c' key\n");
	printf("\t Turn on/offa filter pressing the 't' key\n");
	printf("\t Change the bright by pressing the 'k' key\n");
	printf("\t Change the typo of boundaries between non-cyclic/cyclic by pressing the 'w' key\n");
	printf("\t Turn on/off the advection/vector_advection/project by pressing the a/s/p keys\n");
	printf("\t Press the 'u' key to change how advec works \n\t\t The effect obtained is not realistic but is cool \n\t\t If boundaries are cyclic and project on, cool explosion effects may happen\n ");
	printf("\t Quit by pressing the 'q' key\n");

	dvel = false;
	
	solver.Init(N, dt, diff, visc, Tadvec, Tadvec2, Tproject, TWalls1, TWalls2);

	if (!solver.AllocateData()) {
		exit(1);
	}
	
	win_x = 512;
	win_y = 512;
	OpenGlutWindow();

	glutMainLoop();

	exit(0);
}