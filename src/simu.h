#ifndef _SIMU_H_
#define _SIMU_H_

void computeAcceleration(struct world * hair, struct world *collider, struct point a[numPoints]);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the hair structure accordingly
void Euler(struct world * hair, struct world * collider);

void stretchSpringForce(struct world *hair);
void smoothing(struct world *hair);
void getFrames(struct world *hair, double frames[numPoints][3][3]);
void transpose(double m[3][3], double mT[3][3]);
void frameTimesVector(double frame[3][3], struct point& vec, struct point& reference);
void getInitialReferenceVectors(struct world *hair);
void getReferenceVectors(struct world *hair);
void bendSpringForce(struct world *hair);
void gravity(struct world *hair);
void wind(struct world *hair); 
void collisions(struct world *hair, struct world *collider);
void dampingForces(struct world *hair, struct point a[numPoints]);
void stretchDampingForce(struct world *hair); 
void bendDampingForce(struct world *hair);
void integrateForces(struct world* hair);
void coreSpringForce(struct world* hair);
#endif
