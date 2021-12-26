#ifndef _HAIR_H_
#define _HAIR_H_

#include "GLheader.h"
#include "pic.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <map>

#define pi 3.141592653589793238462643383279 
#define numPoints 10
#define strandOffset 0.5
#define numHairs 9
#define strandOffsetMulti 0.2

// camera angles
extern double Theta;
extern double Phi;
extern double R;

// number of images saved to disk so far
extern int sprite;

// mouse control
extern int g_vMousePos[2];
extern int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;

struct point 
{
   double x;
   double y;
   double z;
};

// these variables control what is displayed on the screen
extern int shear, bend, structural, pause, viewingMode, saveScreenToFile;
extern int up, down, left, right; 

struct world
{
  char integrator[10]; // "RK4" or "Euler"
  double dt; // timestep, e.g.. 0.001
  double dt_damp; // timestep for damping loop
  int n; // display only every nth timepoint
  double kElastic; // Hook's elasticity coefficient for all springs except collision springs
  double dElastic; // Damping coefficient for all springs except collision springs
  double kStretch; // Spring  constant for stretch spring
  double dStretch; // Damping constant for stretch spring
  double kCollision; // Hook's elasticity coefficient for collision springs
  double dCollision; // Damping coefficient collision springs
  double mass; // mass of each of control points, mass assumed to be equal for every control point
  int incPlanePresent; // Is the inclined plane present? 1 = YES, 0 = NO (always NO in this assignment)
  double a, b, c, d; // inclined plane has equation a * x + b * y + c * z + d = 0; if no inclined plane, these four fields are not used
  int resolution; // resolution for the 3d grid specifying the external force field; value of 0 means that there is no force field
  struct point * forceField;  // pointer to the array of values of the force field
  struct point p[numPoints];  // position of control points
  struct point v[numPoints];  // velocities of control points
  struct point f[numPoints];  // forces for each of control points
  struct point p0[numPoints]; // rest positions of control points
  struct point p_smooth[numPoints]; // smoothed positions of control points
  double F[numPoints][3][3];  // local frames; F[i] is 3x3 matrix holding the axes for the ith frame
  double F0[numPoints][3][3]; // initial local frames
  struct point t[numPoints];  // reference vectors
  struct point t0[numPoints]; // initial reference vectors
  bool pushedByWind; // if true, then gets affected by wind
  std::map<int, int> collisions; // maps index of a node to displacement
                                          // from a node on the other strand that
                                          // it collides with
};

extern struct world hair;

// computes crossproduct of three vectors, which are given as points
// struct point vector1, vector2, dest
// result goes into dest
#define CROSSPRODUCTp(vector1,vector2,dest)\
  CROSSPRODUCT( (vector1).x, (vector1).y, (vector1).z,\
                (vector2).x, (vector2).y, (vector2).z,\
                (dest).x, (dest).y, (dest).z )

// computes crossproduct of three vectors, which are specified by floating-point coordinates
// double x1,y1,z1,x2,y2,z2,x,y,z
// result goes into x,y,z
#define CROSSPRODUCT(x1,y1,z1,x2,y2,z2,x,y,z)\
\
  x = (y1) * (z2) - (y2) * (z1);\
  y = (x2) * (z1) - (x1) * (z2);\
  z = (x1) * (y2) - (x2) * (y1)

// normalizes vector dest
// struct point dest
// result returned in dest
// must declare a double variable called 'length' somewhere inside the scope of the NORMALIZE macrp
// macro will change that variable
#define pNORMALIZE(dest)\
\
  length = sqrt((dest).x * (dest).x + (dest).y * (dest).y + (dest).z * (dest).z);\
  (dest).x /= length;\
  (dest).y /= length;\
  (dest).z /= length;

// copies vector source to vector dest
// struct point source,dest
#define pCPY(source,dest)\
\
  (dest).x = (source).x;\
  (dest).y = (source).y;\
  (dest).z = (source).z;
  
// assigns values x,y,z to point vector dest
// struct point dest
// double x,y,z
#define pMAKE(x,y,z,dest)\
\
  (dest).(x) = (x);\
  (dest).(y) = (y);\
  (dest).(z) = (z);

// sums points src1 and src2 to dest
// struct point src1,src2,dest
#define pSUM(src1,src2,dest)\
\
  (dest).x = (src1).x + (src2).x;\
  (dest).y = (src1).y + (src2).y;\
  (dest).z = (src1).z + (src2).z;

// dest = src2 - src1
// struct point src1,src2,dest
#define pDIFFERENCE(src1,src2,dest)\
\
  (dest).x = (src1).x - (src2).x;\
  (dest).y = (src1).y - (src2).y;\
  (dest).z = (src1).z - (src2).z;

// mulitplies components of point src by scalar and returns the result in dest
// struct point src,dest
// double scalar
#define pMULTIPLY(src,scalar,dest)\
\
  (dest).x = (src).x * (scalar);\
  (dest).y = (src).y * (scalar);\
  (dest).z = (src).z * (scalar);

// Performs dot product of src1 and src2 and returns as scalar in dest
// struct point src1, src2
#define pDOT(src1, src2, dest)\
\
  (dest) = ((src1).x * (src2).x) + ((src1).y * (src2).y) + ((src1).z * (src2).z);

// Calculates length of vector src and stores as scalar in dest
#define pLENGTH(src, dest)\
\
  (dest) = sqrt((src).x * (src).x + (src).y * (src).y + (src).z * (src).z);

#endif