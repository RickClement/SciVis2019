// Usage: Drag with the mouse to add smoke to the fluid. This will also move a "rotor" that disturbs
//        the velocity field at the mouse location. Press the indicated keys to change options
//--------------------------------------------------------------------------------------------------

#include <rfftw.h>              //the numerical simulation FFTW library
#include <stdio.h>              //for printing the help text
#include <math.h>               //for various math functions
#include <GL/glut.h>            //the GLUT graphics library
#include <algorithm>
#include "Quad.h"
#include "glui.h"




//--- SIMULATION PARAMETERS ------------------------------------------------------------------------
const int DIM = 50;				//size of simulation grid
double dt = 0.4;				//simulation time step
float visc = 0.001;				//fluid viscosity
fftw_real *vx, *vy;             //(vx,vy)   = velocity field at the current moment
fftw_real *vx0, *vy0;           //(vx0,vy0) = velocity field at the previous moment
fftw_real *fx, *fy;	            //(fx,fy)   = user-controlled simulation forces, steered with the mouse
fftw_real *rho, *rho0;			//smoke density at the current (rho) and previous (rho0) moment
rfftwnd_plan plan_rc, plan_cr;  //simulation domain discretization


//--- VISUALIZATION PARAMETERS ---------------------------------------------------------------------
int   winWidth, winHeight;      //size of the graphics window, in pixels
int   color_dir = 0;            //use direction color-coding or not
float vec_scale = 1000;			//scaling of hedgehogs
int   draw_smoke = 0;           //draw the smoke or not
int   draw_vecs = 1;            //draw the vector field or not
const int COLOR_BLACKWHITE=0;   //different types of color mapping: black-and-white, rainbow, heat
const int COLOR_RAINBOW=1;
const int COLOR_HEAT=2;
const int VIS_DENSITY=0;
const int VIS_VELOCITY=1;
const int VIS_FORCE=2;
int   scalar_col = 0;           //method for scalar coloring
int   vis = 0;					//toggles between visualisation datasets
int   frozen = 0;               //toggles on/off the animation
int   numcols = 128;			//parameterises the number of colours in the colourmap
int   scalclam = 0;				//toggles between colormap scaling or clamping. 0 means scaling, 1 means clamping
float minValueData =  999;
float maxValueData = -999;

GLUI *glui;

//------ SIMULATION CODE STARTS HERE -----------------------------------------------------------------


//init_simulation: Initialize simulation data structures as a function of the grid size 'n'.
//                 Although the simulation takes place on a 2D grid, we allocate all data structures as 1D arrays,
//                 for compatibility with the FFTW numerical library.
void init_simulation(int n)
{
	int i; size_t dim;

	dim     = n * 2*(n/2+1)*sizeof(fftw_real);        //Allocate data structures
	vx       = (fftw_real*) malloc(dim);
	vy       = (fftw_real*) malloc(dim);
	vx0      = (fftw_real*) malloc(dim);
	vy0      = (fftw_real*) malloc(dim);
	dim     = n * n * sizeof(fftw_real);
	fx      = (fftw_real*) malloc(dim);
	fy      = (fftw_real*) malloc(dim);
	rho     = (fftw_real*) malloc(dim);
	rho0    = (fftw_real*) malloc(dim);
	plan_rc = rfftw2d_create_plan(n, n, FFTW_REAL_TO_COMPLEX, FFTW_IN_PLACE);
	plan_cr = rfftw2d_create_plan(n, n, FFTW_COMPLEX_TO_REAL, FFTW_IN_PLACE);

	for (i = 0; i < n * n; i++)                      //Initialize data structures to 0
	{ vx[i] = vy[i] = vx0[i] = vy0[i] = fx[i] = fy[i] = rho[i] = rho0[i] = 0.0f; }
}

//FFT: Execute the Fast Fourier Transform on the dataset 'vx'.
//     'direction' indicates if we do the direct (1) or inverse (-1) Fourier Transform
void FFT(int direction,void* vx)
{
	if(direction==1) rfftwnd_one_real_to_complex(plan_rc,(fftw_real*)vx,(fftw_complex*)vx);
	else             rfftwnd_one_complex_to_real(plan_cr,(fftw_complex*)vx,(fftw_real*)vx);
}

int clamp(float x)
{ return ((x)>=0.0?((int)(x)):(-((int)(1-(x))))); }

float max(float x, float y)
{ return x < y ? x : y; }

//solve: Solve (compute) one step of the fluid flow simulation
void solve(int n, fftw_real* vx, fftw_real* vy, fftw_real* vx0, fftw_real* vy0, fftw_real visc, fftw_real dt)
{
	fftw_real x, y, x0, y0, f, r, U[2], V[2], s, t;
	int i, j, i0, j0, i1, j1;

	for (i=0;i<n*n;i++)
	{ vx[i] += dt*vx0[i]; vx0[i] = vx[i]; vy[i] += dt*vy0[i]; vy0[i] = vy[i]; }

	for ( x=0.5f/n,i=0 ; i<n ; i++,x+=1.0f/n )
	   for ( y=0.5f/n,j=0 ; j<n ; j++,y+=1.0f/n )
	   {
	      x0 = n*(x-dt*vx0[i+n*j])-0.5f;
	      y0 = n*(y-dt*vy0[i+n*j])-0.5f;
	      i0 = clamp(x0); s = x0-i0;
	      i0 = (n+(i0%n))%n;
	      i1 = (i0+1)%n;
	      j0 = clamp(y0); t = y0-j0;
	      j0 = (n+(j0%n))%n;
	      j1 = (j0+1)%n;
	      vx[i+n*j] = (1-s)*((1-t)*vx0[i0+n*j0]+t*vx0[i0+n*j1])+s*((1-t)*vx0[i1+n*j0]+t*vx0[i1+n*j1]);
	      vy[i+n*j] = (1-s)*((1-t)*vy0[i0+n*j0]+t*vy0[i0+n*j1])+s*((1-t)*vy0[i1+n*j0]+t*vy0[i1+n*j1]);
	   }

	for(i=0; i<n; i++)
	  for(j=0; j<n; j++)
	  {  vx0[i+(n+2)*j] = vx[i+n*j]; vy0[i+(n+2)*j] = vy[i+n*j]; }

	FFT(1,vx0);
	FFT(1,vy0);

	for (i=0;i<=n;i+=2)
	{
	   x = 0.5f*i;
	   for (j=0;j<n;j++)
	   {
	      y = j<=n/2 ? (fftw_real)j : (fftw_real)j-n;
	      r = x*x+y*y;
	      if ( r==0.0f ) continue;
	      f = (fftw_real)exp(-r*dt*visc);
	      U[0] = vx0[i  +(n+2)*j]; V[0] = vy0[i  +(n+2)*j];
	      U[1] = vx0[i+1+(n+2)*j]; V[1] = vy0[i+1+(n+2)*j];

	      vx0[i  +(n+2)*j] = f*((1-x*x/r)*U[0]     -x*y/r *V[0]);
	      vx0[i+1+(n+2)*j] = f*((1-x*x/r)*U[1]     -x*y/r *V[1]);
	      vy0[i+  (n+2)*j] = f*(  -y*x/r *U[0] + (1-y*y/r)*V[0]);
	      vy0[i+1+(n+2)*j] = f*(  -y*x/r *U[1] + (1-y*y/r)*V[1]);
	   }
	}

	FFT(-1,vx0);
	FFT(-1,vy0);

	f = 1.0/(n*n);
 	for (i=0;i<n;i++)
	   for (j=0;j<n;j++)
	   { vx[i+n*j] = f*vx0[i+(n+2)*j]; vy[i+n*j] = f*vy0[i+(n+2)*j]; }
}


// diffuse_matter: This function diffuses matter that has been placed in the velocity field. It's almost identical to the
// velocity diffusion step in the function above. The input matter densities are in rho0 and the result is written into rho.
void diffuse_matter(int n, fftw_real *vx, fftw_real *vy, fftw_real *rho, fftw_real *rho0, fftw_real dt)
{
	fftw_real x, y, x0, y0, s, t;
	int i, j, i0, j0, i1, j1;

	for ( x=0.5f/n,i=0 ; i<n ; i++,x+=1.0f/n )
		for ( y=0.5f/n,j=0 ; j<n ; j++,y+=1.0f/n )
		{
			x0 = n*(x-dt*vx[i+n*j])-0.5f;
			y0 = n*(y-dt*vy[i+n*j])-0.5f;
			i0 = clamp(x0);
			s = x0-i0;
			i0 = (n+(i0%n))%n;
			i1 = (i0+1)%n;
			j0 = clamp(y0);
			t = y0-j0;
			j0 = (n+(j0%n))%n;
			j1 = (j0+1)%n;
			rho[i+n*j] = (1-s)*((1-t)*rho0[i0+n*j0]+t*rho0[i0+n*j1])+s*((1-t)*rho0[i1+n*j0]+t*rho0[i1+n*j1]);
		}
}

//set_forces: copy user-controlled forces to the force vectors that are sent to the solver.
//            Also dampen forces and matter density to get a stable simulation.
void set_forces(void)
{
	int i;
	for (i = 0; i < DIM * DIM; i++)
	{
        rho0[i]  = 0.995 * rho[i];
        fx[i] *= 0.85;
        fy[i] *= 0.85;
        vx0[i]    = fx[i];
        vy0[i]    = fy[i];
	}
}


//do_one_simulation_step: Do one complete cycle of the simulation:
//      - set_forces:
//      - solve:            read forces from the user
//      - diffuse_matter:   compute a new set of velocities
//      - gluPostRedisplay: draw a new visualization frame
void do_one_simulation_step(void)
{
	if (!frozen)
	{
	  set_forces();
	  solve(DIM, vx, vy, vx0, vy0, visc, dt);
	  diffuse_matter(DIM, vx, vy, rho, rho0, dt);
	  glutPostRedisplay();
	}
}


//------ VISUALIZATION CODE STARTS HERE -----------------------------------------------------------------


//rainbow: Implements a color palette, mapping the scalar 'value' to a rainbow color RGB
void rainbow(float value,float* R,float* G,float* B)
{
   const float dx=0.8;
   value = (6-2*dx)*value+dx;
   *R = max(1.0,(float)((3-fabs(value-4)-fabs(value-5))/2.0));
   *G = max(1.0,(float)((4-fabs(value-2)-fabs(value-4))/2.0));
   *B = max(1.0,(float)((3-fabs(value-1)-fabs(value-2))/2.0));
}

//heatmap: Implements a red->yellow->white heat inspired colour palette.
void heatmap(float value, float* R, float* G, float* B)
{
   const float dx=0.8;
   value = (6-2*dx)*value+dx;
   *R = max(1.0,(float)((6-fabs(value-3)-fabs(value-4))/2.0));
   *G = max(1.0,(float)((4-fabs(value-4)-fabs(value-5))/2.0));
   *B = max(1.0,(float)((3-fabs(value-5)-fabs(value-5))/2.0));
}

//set_colormap: Sets three different types of colormaps
void set_colormap(float vy)
{
   float R,G,B;

   vy *= numcols;
   vy = (int)(vy);
   vy/= numcols;

   if(scalclam) {
	   // Clamp the data between the given min and max values
	   if(vy < minValueData) {
		   vy = minValueData;
	   } else if (vy > maxValueData) {
		   vy = maxValueData;
	   }
   }

   // Normalise vy to the interval 0-1.
   float interval = maxValueData - minValueData;
   vy = (vy - minValueData) / interval;

   if (scalar_col==COLOR_BLACKWHITE) {
       R = G = B = vy;
   }
   else if (scalar_col==COLOR_RAINBOW) {
	   rainbow(vy,&R,&G,&B);
   }
   else if (scalar_col==COLOR_HEAT) {
       heatmap(vy,&R,&G,&B);
   }

   glColor3f(R,G,B);
}


//direction_to_color: Set the current color by mapping a direction vector (x,y), using
//                    the color mapping method 'method'. If method==1, map the vector direction
//                    using a rainbow colormap. If method==0, simply use the white color
void direction_to_color(float x, float y, int method)
{
	float r,g,b,f;
	if (method)
	{
	  f = atan2(y,x) / 3.1415927 + 1;
	  r = f;
	  if(r > 1) r = 2 - r;
	  g = f + .66667;
      if(g > 2) g -= 2;
	  if(g > 1) g = 2 - g;
	  b = f + 2 * .66667;
	  if(b > 2) b -= 2;
	  if(b > 1) b = 2 - b;
	}
	else
	{ r = g = b = 1; }
	glColor3f(r,g,b);
}

//drawColorLegend: Draws the colour legend using the current colour map.
void drawColorLegend(){
    /*
    Quad q;
    glColor3f(1,0,0);
    q.addPoint(5,5,0);
    q.addPoint(5,30,0);
    q.addPoint(winWidth-5,30,0);
    q.addPoint(winWidth-5,5,0);
    q.draw();
    //printf("Ik heb een quad gemaakt!\n");
    */
    
    float interval;
    float value = 0;
    float width = 0;
    
	glBegin(GL_QUAD_STRIP);
	interval = 1 / (float) numcols;
	minValueData = 0;
	maxValueData = 1;
    while(value <= 1){
		set_colormap(value);
		width = value * winWidth;
		glVertex2i(width,15);
		glVertex2i(width,0);
		glVertex2i(width+interval,0);
		glVertex2i(width+interval,15);
		
		value += interval;
	}
	glEnd();
	
}

//Probably something wrong with it: array size seems to remain at 1. I don't get this data structure...
int getSizeArray(fftw_real *array){
    /**
    int size = sizeof(array)/sizeof(array[0]);
    return size;
     **/
     /**
     int size = 0;
     fftw_real j = 1;
    while (j != 0){
        j = array[size];
        size++;
    }
    return size;
    **/
    return DIM * DIM;
}

void updateMinValue(fftw_real *data) {

	if(scalclam == 0) { //only when we are scaling

		int arraySize = getSizeArray(data);
		fftw_real minValue = 9999;

		for (int i = 0; i < arraySize; i++) {
			if (data[i] < minValue) {
				minValue = data[i];
			}
		}
		minValueData = minValue;
		glui->sync_live();

	}

}

void updateMaxValue(fftw_real *data) {

	if(scalclam == 0) { //only when we are scaling

		int arraySize = getSizeArray(data);
		fftw_real maxValue = 0;

		for (int i = 0; i < arraySize; i++) {
			if (data[i] > maxValue) {
				maxValue = data[i];
			}
		}
		maxValueData = maxValue;
		glui->sync_live();
	}


}

void updateMinMaxValues(fftw_real *data){
	updateMinValue(data);
	updateMaxValue(data);
}

void updateMinMaxValues(){

    if (scalclam == 0){ // When scaling, use the values from the data set as out min and max

        switch(vis) {
            case VIS_DENSITY:
                updateMinMaxValues(rho);

                /**
            case VIS_VELOCITY:
                return sqrt((pow(vx[idx], 2) + pow(vy[idx], 2)));
            case VIS_FORCE:
                return sqrt((pow(fx[idx], 2) + pow(fy[idx], 2)));
                 **/
            default:
                break;
        }
    }
    else{ // When clamping, use the values as set in the clamp

    }
}


fftw_real getVariable(int idx){
	switch(vis){
			//case VIS_DENSITY: updateMinMaxValues(rho); return rho[idx];
            case VIS_DENSITY:  return rho[idx];
			case VIS_VELOCITY: return sqrt((pow(vx[idx],2)+pow(vy[idx],2)));
			case VIS_FORCE:    return sqrt((pow(fx[idx],2)+pow(fy[idx],2)));
		}
}




//visualize: This is the main visualization function
void visualize(void)
{
	int        i, j, idx; double px,py;
	fftw_real  wn = (fftw_real)winWidth / (fftw_real)(DIM + 1);   // Grid cell width
	fftw_real  hn = (fftw_real)winHeight / (fftw_real)(DIM + 1);  // Grid cell heigh
	
	fftw_real variable;

	if (draw_smoke)
	{
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	for (j = 0; j < DIM - 1; j++)			//draw smoke
	{

	    glBegin(GL_TRIANGLE_STRIP);

		i = 0;
		px = wn + (fftw_real)i * wn;
		py = hn + (fftw_real)j * hn;
		idx = (j * DIM) + i;
		
		variable = getVariable(idx);
		glColor3f(variable,variable,variable);
		glVertex2f(px,py);

		for (i = 0; i < DIM - 1; i++)
		{
			px = wn + (fftw_real)i * wn;
			py = hn + (fftw_real)(j + 1) * hn;
			idx = ((j + 1) * DIM) + i;
			variable = getVariable(idx);
			set_colormap(variable);
			glVertex2f(px, py);
			px = wn + (fftw_real)(i + 1) * wn;
			py = hn + (fftw_real)j * hn;
			idx = (j * DIM) + (i + 1);
			variable = getVariable(idx);
			set_colormap(variable);
			glVertex2f(px, py);
		}

        updateMinMaxValues();

		px = wn + (fftw_real)(DIM - 1) * wn;
		py = hn + (fftw_real)(j + 1) * hn;
		idx = ((j + 1) * DIM) + (DIM - 1);
		variable = getVariable(idx);
		set_colormap(variable);
		glVertex2f(px, py);
		glEnd();
	}
	}

	if (draw_vecs)
	{
	  glBegin(GL_LINES);				//draw velocities
	  for (i = 0; i < DIM; i++)
	    for (j = 0; j < DIM; j++)
	    {
		  idx = (j * DIM) + i;
		  direction_to_color(vx[idx],vy[idx],color_dir);
		  glVertex2f(wn + (fftw_real)i * wn, hn + (fftw_real)j * hn);
		  glVertex2f((wn + (fftw_real)i * wn) + vec_scale * vx[idx], (hn + (fftw_real)j * hn) + vec_scale * vy[idx]);
	    }
	  glEnd();
	}
	drawColorLegend();
}


//------ INTERACTION CODE STARTS HERE -----------------------------------------------------------------

//display: Handle window redrawing events. Simply delegates to visualize().
void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	visualize();
	glFlush();
	glutSwapBuffers();
}

//reshape: Handle window resizing (reshaping) events
void reshape(int w, int h)
{
 	glViewport(0.0f, 0.0f, (GLfloat)w, (GLfloat)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, (GLdouble)w, 0.0, (GLdouble)h);
	winWidth = w; winHeight = h;
}

//keyboard: Handle key presses
void keyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 't': dt -= 0.001; glui->sync_live(); break;
        case 'T': dt += 0.001; glui->sync_live(); break;
        case 'c': color_dir = 1 - color_dir; glui->sync_live(); break;
        case 'S': vec_scale *= 1.2; glui->sync_live(); break;
        case 's': vec_scale *= 0.8; glui->sync_live(); break;
        case 'V': visc *= 5; glui->sync_live(); break;
        case 'v': visc *= 0.2; glui->sync_live(); break;
        case 'x': draw_smoke = 1 - draw_smoke;
            if (draw_smoke==0) draw_vecs = 1; glui->sync_live(); break;
        case 'y': draw_vecs = 1 - draw_vecs;
            if (draw_vecs==0) draw_smoke = 1; glui->sync_live(); break;
        case 'm': scalar_col = (scalar_col + 1) % 3; glui->sync_live(); break;
        case 'a': frozen = 1-frozen; glui->sync_live(); break;
        case 'n': scalclam = 1 - scalclam; glui->sync_live(); break;
        case 'b': vis = (vis + 1) % 3; minValueData = 999; maxValueData = -999; glui->sync_live(); break;
        case '+': if(numcols < 255) numcols += 1; glui->sync_live(); break;
        case '-': if(numcols > 1)   numcols -= 1; glui->sync_live(); break;
        case 'q': exit(0);
    }
}



// drag: When the user drags with the mouse, add a force that corresponds to the direction of the mouse
//       cursor movement. Also inject some new matter into the field at the mouse location.
void drag(int mx, int my)
{
	int xi,yi,X,Y; double  dx, dy, len;
	static int lmx=0,lmy=0;				//remembers last mouse location

	// Compute the array index that corresponds to the cursor location
	xi = (int)clamp((double)(DIM + 1) * ((double)mx / (double)winWidth));
	yi = (int)clamp((double)(DIM + 1) * ((double)(winHeight - my) / (double)winHeight));

	X = xi; Y = yi;

	if (X > (DIM - 1))  X = DIM - 1; if (Y > (DIM - 1))  Y = DIM - 1;
	if (X < 0) X = 0; if (Y < 0) Y = 0;

	// Add force at the cursor location
	my = winHeight - my;
	dx = mx - lmx; dy = my - lmy;
	len = sqrt(dx * dx + dy * dy);
	if (len != 0.0) {  dx *= 0.1 / len; dy *= 0.1 / len; }
	fx[Y * DIM + X] += dx;
	fy[Y * DIM + X] += dy;
	rho[Y * DIM + X] = 10.0f;
	lmx = mx; lmy = my;
}

void GLUI_interface(GLUI *glui){

    GLUI_Checkbox   *checkboxScaling, *checkboxDirectionColoring, *checkboxDrawMatter, *checkboxHedgehogs;
    GLUI_Spinner    *numberOfColours;
    GLUI_RadioGroup *radio, *radio2;
    GLUI_EditText   *edittext, *edittext2;

    //GLUI *glui = GLUI_Master.create_glui( "Options", 0, 400, 50 ); /* name, flags, x, and y */
    new GLUI_StaticText( glui, "Additional options visualization" );
    new GLUI_Separator( glui );
    //checkbox = new GLUI_Checkbox( glui, "Cool color button!", &scalar_col, 1 );
    numberOfColours  = new GLUI_Spinner( glui, "Segments colormap:", &numcols, 2);
    numberOfColours->set_int_limits( 3, 256 );

	//new GLUI_StaticText( glui, minValueData);
    //new GLUI_StaticText(glui, )

    checkboxScaling = new GLUI_Checkbox( glui, "Use scaling (rather than clamping) (n)", &scalclam, 1 );
    checkboxDirectionColoring = new GLUI_Checkbox( glui, "Toggle direction coloring (c)", &color_dir, 1 );
    checkboxDrawMatter = new GLUI_Checkbox( glui, "Toggle drawing matter (x)", &draw_smoke, 1 );
    checkboxHedgehogs = new GLUI_Checkbox( glui, "Toggle drawing hedgehogs (y)", &draw_vecs, 1 );

	GLUI_Panel *obj_panel = new GLUI_Panel( glui, "Colormap type (m)" );
	radio = new GLUI_RadioGroup( obj_panel, &scalar_col);
	new GLUI_RadioButton( radio, "Black/white" );
	new GLUI_RadioButton( radio, "Rainbow" );
	new GLUI_RadioButton( radio, "Heatmap" );

    GLUI_Panel *obj_panel2 = new GLUI_Panel( glui, "Data to visualize (b)" );
    radio2 = new GLUI_RadioGroup( obj_panel2, &vis);
    new GLUI_RadioButton( radio2, "Fluid density" );
    new GLUI_RadioButton( radio2, "Fluid velocity magnitude" );
    new GLUI_RadioButton( radio2, "Force field magnitude" );

    edittext = new GLUI_EditText( glui, "Min value data points:", &minValueData, 3);
    edittext2 = new GLUI_EditText( glui, "Max value data points:", &maxValueData, 3);

//    glui->sync_live();

    /**
    edittext = new GLUI_EditText( glui, "Text:", text, 3, control_cb );
    GLUI_Panel *obj_panel = new GLUI_Panel( glui, "Object Type" );
    radio = new GLUI_RadioGroup( obj_panel,&obj,4,control_cb );
    new GLUI_RadioButton( radio, "Sphere" );
    new GLUI_RadioButton( radio, "Torus" );
    new GLUI_RadioButton( radio, "Teapot" );
    new GLUI_Button( glui, "Quit", 0,(GLUI_Update_CB)exit );

    glui->set_main_gfx_window( main_window );

    // We register the idle callback with GLUI, *not* with GLUT
    //GLUI_Master.set_glutIdleFunc( myGlutIdle );
    GLUI_Master.set_glutIdleFunc( NULL );
    **/
}


//main: The main program
int main(int argc, char **argv)
{
	printf("Fluid Flow Simulation and Visualization\n");
	printf("=======================================\n");
	printf("Click and drag the mouse to steer the flow!\n");
	printf("T/t:   increase/decrease simulation timestep\n");
	printf("S/s:   increase/decrease hedgehog scaling\n");
	printf("c:     toggle direction coloring on/off\n");
	printf("V/v:   increase decrease fluid viscosity\n");
	printf("x:     toggle drawing matter on/off\n");
	printf("y:     toggle drawing hedgehogs on/off\n");
	printf("m:     toggle thru scalar coloring\n");
	printf("a:     toggle the animation on/off\n");
	printf("n	   toggle between scaling and clamping\n");
	printf("b:     toggle through visualisations\n");
	printf("+      increase colours in the colourmap\n");
	printf("-      decrease colours in the colourmap\n");
	printf("q:     quit\n\n");


	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(500,500);
	glutCreateWindow("Real-time smoke simulation and visualization");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(do_one_simulation_step);
    glutKeyboardFunc(keyboard);

    glutMotionFunc(drag);

    init_simulation(DIM);    //initialize the simulation data structures

    glui = GLUI_Master.create_glui( "Options", 0, 400, 50 );
    GLUI_interface(glui);
    glutMainLoop();            //calls do_one_simulation_step, keyboard, display, drag, reshape
    return 0;
}
