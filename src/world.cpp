#include <string>
#include <stdio.h>
#include <math>
#include <stdlib.h>

#define numPoints 10
#define pi 3.141592653589793238462643383279 

struct point 
{
   double x;
   double y;
   double z;
};

struct world
{
  char integrator[10]; // "RK4" or "Euler"
  double dt; // timestep, e.g.. 0.001
  double dt_damp; //damping loop timestep
  int n; // display only every nth timestep
  double kElastic; // Hook's elasticity coefficient for all springs except collision springs
  double dElastic; // Damping coefficient for all springs except collision springs
  double kStretch; // Spring  constant for stretch spring
  double dStretch; // Damping constant for stretch spring
  double kCollision; // Hook's elasticity coefficient for collision springs
  double dCollision; // Damping coefficient collision springs
  double mass; // mass of each of the 512 control points, mass assumed to be equal for every control point
  int incPlanePresent; // Is the inclined plane present? 1 = YES, 0 = NO
  double a,b,c,d; // inclined plane has equation a * x + b * y + c * z + d = 0; if no inclined plane, these four fields are not used
  int resolution; // resolution for the 3d grid specifying the external force field; value of 0 means that there is no force field
  struct point * forceField; // pointer to the array of values of the force field
  struct point p[numPoints]; // position of the 512 control points
  struct point v[numPoints]; // velocities of the 512 control points
};


/* writes the world parameters to a world file on disk*/
/* fileName = string containing the name of the output world file, ex: xxx.w */
void writeWorld(const char * fileName, struct world * hairs)
{
  int i,j,k;
  FILE * file;
  
  file = fopen(fileName, "w");
  if (file == NULL) {
    printf ("can't open file\n");
    exit(1);
  }

  /* write integrator algorithm */ 
  fprintf(file,"%s\n",hairs->integrator);

  /* write timestep */
  fprintf(file,"%lf %d\n",hairs->dt,hairs->n);
  fprintf(file,"%lf %d\n",hairs->dt_damp,hairs->n);

  /* write physical parameters */
  fprintf(file, "%lf %lf %lf %lf %lf %lf\n", 
    hairs->kElastic, hairs->dElastic, hairs->dStretch, hairs->kStretch, hairs->kCollision, hairs->dCollision);

  /* write mass */
  fprintf(file, "%lf\n", hairs->mass);

  /* write info about the plane */
  fprintf(file, "%d\n", hairs->incPlanePresent);
  if (hairs->incPlanePresent == 1)
    fprintf(file, "%lf %lf %lf %lf\n", hairs->a, hairs->b, hairs->c, hairs->d);

  /* write info about the force field */
  fprintf(file, "%d\n", hairs->resolution);
  if (hairs->resolution != 0)
    for (i=0; i<= hairs->resolution-1; i++)
      for (j=0; j<= hairs->resolution-1; j++)
        for (k=0; k<= hairs->resolution-1; k++)
          fprintf(file, "%lf %lf %lf\n", 
             hairs->forceField[i * hairs->resolution * hairs->resolution + j * hairs->resolution + k].x, 
             hairs->forceField[i * hairs->resolution * hairs->resolution + j * hairs->resolution + k].y, 
             hairs->forceField[i * hairs->resolution * hairs->resolution + j * hairs->resolution + k].z);
  

  /* write initial point positions */
  for (i = 0; i < numPoints; i++)
  {
    fprintf(file, "%lf %lf %lf\n", 
      hairs->p[i].x, hairs->p[i].y, hairs->p[i].z);
  }
      
  /* write initial point velocities */
  for (i = 0; i < numPoints; i++)
  {
    fprintf(file, "%lf %lf %lf\n", 
      hairs->v[i].x, hairs->v[i].y, hairs->v[i].z);
  }

  fclose(file);
  
  return;
}

/* modify main to create your own world */
int main()
{
  struct world hairs;
  int i,j,k;
  double x,y,z;

  // set the integrator and the physical parameters
  // the values below are EXAMPLES, to be modified by you as needed
  strcpy(hairs.integrator,"Euler");
  hairs.dt=0.00005000;
  hairs.dt_damp = 0.00006000; 
  hairs.n=1;
  hairs.kElastic=800;
  hairs.dElastic=40;
  hairs.kStretch = 4000;
  hairs.dStretch = 450; 
  hairs.kCollision = 1400.0;
  hairs.dCollision=0.25;
  hairs.mass= 5; // / (double)512;

  // set the inclined plane (not used in this assignment; ignore)
  hairs.incPlanePresent=1;
  hairs.a=-1;
  hairs.b=1;
  hairs.c=1;
  hairs.d=2;

  // set the external force field
  hairs.resolution=3;   // FIX: CHANGED FROM 30 
  hairs.forceField = 
    (struct point *)malloc(hairs.resolution*hairs.resolution*hairs.resolution*sizeof(struct point));
  for (i=0; i<= hairs.resolution-1; i++)
    for (j=0; j<= hairs.resolution-1; j++)
      for (k=0; k<= hairs.resolution-1; k++)
      {
        // set the force at node i,j,k
        // actual space location = x,y,z
        x = -2 + 4*(1.0 * i / (hairs.resolution-1));
        y = -2 + 4*(1.0 * j / (hairs.resolution-1));
        z = -2 + 4*(1.0 * k / (hairs.resolution-1));

        hairs.forceField[i * hairs.resolution * hairs.resolution 
          + j * hairs.resolution + k].x = 0; 
        hairs.forceField[i * hairs.resolution * hairs.resolution 
          + j * hairs.resolution + k].y = 0;
        hairs.forceField[i * hairs.resolution * hairs.resolution 
          + j * hairs.resolution + k].z = 0;
      }

  double radius = 0.25;
  double c = 0.05;
  double radians = 0.0;

  // set the positions of control points
  for (i=0; i< numPoints; i++)
  {
    radians = -((i / 5.0) * 360.0 * pi/180.0);
    hairs.p[i].x = radius * cos(radians);
    hairs.p[i].y = radius * sin(radians);
    hairs.p[i].z = c * radians + 1.5;
  }

  // set the velocities of control points
  for (i=0; i< numPoints; i++)
  {
    hairs.v[i].x=0.0;
    hairs.v[i].y=0.0;
    hairs.v[i].z=0.0;
  }

  // write the hairs variable out to file on disk
  // change hairs.w to whatever you need
  writeWorld("world/hair.w",&hairs);

  return 0;
}


