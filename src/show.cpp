#include "hair.h"
#include "show.h"

void showHair(struct world * hair)
{
  int i, ip; 

  #define PROCESS_NEIGHBOUR(di) \
    ip=i+(di);\
    if\
    (!( (ip>=numPoints) || (ip<0) ) && !((i>=numPoints) || (i<0))) \
    {\
      glVertex3f(hair->p[i].x,hair->p[i].y,hair->p[i].z);\
      glVertex3f(hair->p[ip].x,hair->p[ip].y,hair->p[ip].z);\
    }\

  // render wireframe as hair 
  glLineWidth(1);
	glPointSize(5);
	glDisable(GL_LIGHTING);

	for (int i=0; i<numPoints; i++)
	{
		glBegin(GL_POINTS); // draw point
    // printf("%f\n", hair->forceField->y );
		glColor4f(hair->forceField->x / 10,hair->forceField->y/ 10 ,hair->forceField->z /10,0);  
		glVertex3f(hair->p[i].x,hair->p[i].y,hair->p[i].z);        
		glEnd();

		glBegin(GL_LINES);      
		// structural
		if (structural == 1)
		{
		  glColor4f(0,0,1,1);
		  PROCESS_NEIGHBOUR(1);
		}
                     
		glEnd();
    
		// glEnable(GL_LIGHTING);
    }
  glFrontFace(GL_CCW);
}

void showBoundingBox()
{
  int i,j;

  glColor4f(0.6,0.6,0.6,0);

  glBegin(GL_LINES);

  // front face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(i,-2,-2);
    glVertex3f(i,-2,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(-2,-2,j);
    glVertex3f(2,-2,j);
  }

  // back face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(i,2,-2);
    glVertex3f(i,2,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(-2,2,j);
    glVertex3f(2,2,j);
  }

  // left face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(-2,i,-2);
    glVertex3f(-2,i,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(-2,-2,j);
    glVertex3f(-2,2,j);
  }

  // right face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(2,i,-2);
    glVertex3f(2,i,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(2,-2,j);
    glVertex3f(2,2,j);
  }
  
  glEnd();

  return;
}
