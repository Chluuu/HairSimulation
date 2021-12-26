#define _CRT_SECURE_NO_WARNINGS

#include "hair.h"
#include "input.h"
#include "simu.h"

/* Write a screenshot, in the PPM format, to the specified filename, in PPM format */
void saveScreenshot(int windowWidth, int windowHeight, char *filename)
{
  if (filename == NULL)
    return;

  // Allocate a picture buffer 
  Pic * in = pic_alloc(windowWidth, windowHeight, 3, NULL);

  printf("File to save to: %s\n", filename);

  for (int i=windowHeight-1; i>=0; i--) 
  {
    glReadPixels(0, windowHeight-i-1, windowWidth, 1, GL_RGB, GL_UNSIGNED_BYTE,
      &in->pix[i*in->nx*in->bpp]);
  }

  if (ppm_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);
}

/* converts mouse drags into information about rotation/translation/scaling */
void mouseMotionDrag(int x, int y)
{
  int vMouseDelta[2] = {x-g_vMousePos[0], y-g_vMousePos[1]};

  if (g_iRightMouseButton) // handle camera rotations
  {
    Phi += vMouseDelta[0] * 0.01;
    Theta += vMouseDelta[1] * 0.01;
    
    if (Phi>2*pi)
      Phi -= 2*pi;
    
    if (Phi<0)
      Phi += 2*pi;
    
    if (Theta>pi / 2 - 0.01) // dont let the point enter the north pole
      Theta = pi / 2 - 0.01;
    
    if (Theta<- pi / 2 + 0.01)
      Theta = -pi / 2 + 0.01;
    
    g_vMousePos[0] = x;
    g_vMousePos[1] = y;
  }
}

void mouseMotion (int x, int y)
{
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

void mouseButton(int button, int state, int x, int y)
{
  switch (button)
  {
    case GLUT_LEFT_BUTTON:
      g_iLeftMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_MIDDLE_BUTTON:
      g_iMiddleMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_RIGHT_BUTTON:
      g_iRightMouseButton = (state==GLUT_DOWN);
      break;
  }
 
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

// gets called whenever a key is pressed
void keyboardFunc (unsigned char key, int x, int y)
{
  switch (key)
  {
    case 27:
      exit(0);
      break;

    case 'e':
      Theta = pi / 6;
      Phi = pi / 6;
      viewingMode = 0;
      break;

    case 'v':
      viewingMode = 1 - viewingMode;
      break;

    case 'h':
      shear = 1 - shear;
      break;

    case 's':
      structural = 1 - structural;
      break;

    case 'b':
      bend = 1 - bend;
      break;

    case 'p':
      pause = 1 - pause;
      break;

    case 'z':
      R -= 0.2;
      if (R < 0.2)
        R = 0.2;
      break;

    case 'x':
      R += 0.2;
      break;

    case ' ':
      saveScreenToFile = 1 - saveScreenToFile;
      break;
    
    // add keys for wind 
    case 'i':
      up = 1 - up; 
    case 'k':
      down = 1 - down;
    case 'j':
      left = 1 - left;
    case 'l':
      right = 1 - right; 
  }
}

/* reads the world parameters from a world file */
/* fileName = string containing the name of the world file, ex: hair1.w */
/* function fills the structure 'hair' with parameters read from file */
/* structure 'hair' will typically be declared (probably statically, not on the heap)
   by the caller function */
/* function aborts the program if can't access the file */
void readWorld (char * fileName, struct world * hair, double offsetX, double offsetY)
{
  int i,j,k;
  FILE * file;
  
 file= fopen(fileName, "r");
  if (file == NULL) {
    printf ("can't open file\n");
    exit(1);
  }
 
/* 

  File should first contain a line specifying the integrator (EULER or RK4).
  Example: EULER
  
  Then, follows one line specifying the size of the timestep for the integrator, and
  an integer parameter n specifying  that every nth timestep will actually be drawn
  (the other steps will only be used for internal calculation)
  
  Example: 0.001 5
  Now, timestep equals 0.001. Every fifth time point will actually be drawn,
  i.e. frame1 <--> t = 0
  frame2 <--> t = 0.005
  frame3 <--> t = 0.010
  frame4 <--> t = 0.015
  ...
  
  Then, there should be two lines for physical parameters and external acceleration.
  Format is:
    kElastic dElastic kCollision dCollision
    mass
  Here
    kElastic = elastic coefficient of the spring (same for all springs except collision springs)
    dElastic = damping coefficient of the spring (same for all springs except collision springs)
    kCollision = elastic coefficient of collision springs (same for all collision springs)
    dCollision = damping coefficient of collision springs (same for all collision springs)
    mass = mass in kilograms for each of mass points 
    (mass assumed to be the same for all the points; total mass of the hair = numPoints * mass)
  
  Example:
    10000 25 10000 15
    0.002
  
  Then, there should be one or two lines for the inclined plane, with the obvious syntax. 
  If there is no inclined plane, there should be only one line with a 0 value. There
  is no line for the coefficient. Otherwise, there are two lines, first one containing 1,
  and the second one containing the coefficients.
  Note: there is no inclined plane in this assignment (always 0).
  Example:
    1
    0.31 -0.78 0.5 5.39
  
  Next is the forceField block, first with the resolution and then the data, one point per row.
  Example:
    30
    <here 30 * 30 * 30 = 27 000 lines follow, each containing 3 real numbers>
  
  After this, there should be 1024 lines, each containing three floating-point numbers.
  The first 512 lines correspond to initial point locations.
  The last 512 lines correspond to initial point velocities.
  
  There should no blank lines anywhere in the file.

*/
       
  /* read integrator algorithm */ 
  fscanf(file,"%s\n",&hair->integrator);

  /* read timestep size and render */
  fscanf(file,"%lf %d\n",&hair->dt,&hair->n);
  fscanf(file,"%lf %d\n",&hair->dt_damp,&hair->n);

  /* read physical parameters */
  fscanf(file, "%lf %lf %lf %lf %lf %lf\n", 
    &hair->kElastic, &hair->dElastic, &hair->kStretch, &hair->dStretch, &hair->kCollision, &hair->dCollision);

  /* read mass of each of the 512 points */
  fscanf(file, "%lf\n", &hair->mass);

  /* read info about the plane */
  fscanf(file, "%d\n", &hair->incPlanePresent);
  if (hair->incPlanePresent == 1)
    fscanf(file, "%lf %lf %lf %lf\n", &hair->a, &hair->b, &hair->c, &hair->d);

  /* read info about the force field */
  fscanf(file, "%d\n", &hair->resolution);
  hair->forceField = 
    (struct point *)malloc(hair->resolution*hair->resolution*hair->resolution*sizeof(struct point));
  if (hair->resolution != 0)
    for (i=0; i<= hair->resolution-1; i++)
      for (j=0; j<= hair->resolution-1; j++)
        for (k=0; k<= hair->resolution-1; k++) {
          fscanf(file, "%lf %lf %lf\n", 
              &hair->forceField[i * hair->resolution * hair->resolution + j * hair->resolution + k].x, 
              &hair->forceField[i * hair->resolution * hair->resolution + j * hair->resolution + k].y, 
              &hair->forceField[i * hair->resolution * hair->resolution + j * hair->resolution + k].z);
        }
          
  
  /* read initial point positions */
  for (i= 0; i < numPoints; i++)
  {
			fscanf(file, "%lf %lf %lf\n",&hair->p[i].x, &hair->p[i].y, &hair->p[i].z);
            hair->p[i].y += offsetY;
            hair->p[i].x += offsetX;
			hair->p0[i].x = hair->p[i].x;
			hair->p0[i].y = hair->p[i].y;
			hair->p0[i].z = hair->p[i].z;
  }

  smoothing(hair);  // FIXME
  // TODO: might need to initialize frames
  getFrames(hair, hair->F0); 
  getInitialReferenceVectors(hair);
      
  /* read initial point velocities */
  for (i= 0; i < numPoints; i++)
  {
			fscanf(file, "%lf %lf %lf\n",&hair->v[i].x, &hair->v[i].y, &hair->v[i].z);
  }

  fclose(file);
  
  return;
}

/* writes the world parameters to a world file on disk*/
/* fileName = string containing the name of the output world file, ex: hair1.w */
/* function creates the output world file and then fills it corresponding to the contents
   of structure 'hair' */
/* function aborts the program if can't access the file */
void writeWorld (char * fileName, struct world * hair)
{
  int i,j,k;
  FILE * file;
  
  file = fopen(fileName, "w");
  if (file == NULL) {
    printf ("can't open file\n");
    exit(1);
  }

  /* write integrator algorithm */ 
  fprintf(file,"%s\n",hair->integrator);

  /* write timestep */
  fprintf(file,"%lf %d\n",hair->dt,hair->n);
  fprintf(file,"%lf %d\n",hair->dt_damp,hair->n);

  /* write physical parameters */
  fprintf(file, "%lf %lf %lf %lf\n", 
    hair->kElastic, hair->dElastic, hair->kCollision, hair->dCollision);

  fprintf(file, "%lf %lf %lf %lf\n", 
    hair->kStretch, hair->dStretch, hair->kCollision, hair->dCollision);

  /* write mass */
  fprintf(file, "%lf\n", 
    hair->mass);

  /* write info about the plane */
  fprintf(file, "%d\n", hair->incPlanePresent);
  if (hair->incPlanePresent == 1)
    fprintf(file, "%lf %lf %lf %lf\n", hair->a, hair->b, hair->c, hair->d);

  /* write info about the force field */
  fprintf(file, "%d\n", hair->resolution);
  if (hair->resolution != 0)
    for (i=0; i<= hair->resolution-1; i++)
      for (j=0; j<= hair->resolution-1; j++)
        for (k=0; k<= hair->resolution-1; k++)
          fprintf(file, "%lf %lf %lf\n", 
             hair->forceField[i * hair->resolution * hair->resolution + j * hair->resolution + k].x, 
             hair->forceField[i * hair->resolution * hair->resolution + j * hair->resolution + k].y, 
             hair->forceField[i * hair->resolution * hair->resolution + j * hair->resolution + k].z);
  


 /* write initial point positions */
  for (i = 0; i < numPoints ; i++)
  {
        fprintf(file, "%lf %lf %lf\n", 
          hair->p[i].x, hair->p[i].y, hair->p[i].z);
  }
      
  /* write initial point velocities */
  for (i = 0; i < numPoints; i++)
  {
        fprintf(file, "%lf %lf %lf\n", 
          hair->v[i].x, hair->v[i].y, hair->v[i].z);
  }

  fclose(file);
  
  return;
}

