#include "hair.h"
#include "simu.h"




void smoothing(struct world *hair)
{
     struct point displacement;
     double len, avg_rest_len, total_len;
     total_len = 0.0; 
     for (int i = 0; i < (numPoints - 1); i++)
     {
         pDIFFERENCE(hair->p0[i + 1], hair->p0[i], displacement);
         pLENGTH(displacement, len);
         total_len += len;
     }

     avg_rest_len = total_len / (numPoints - 1);  

     double alpha = 0.5; // smoothing amount
     double beta = fmin(1.0, 1.0 - exp(-1.0 * avg_rest_len / alpha));

     struct point d[numPoints - 1];
     struct point d_prior, d_2prior, to_next;
     struct point d_init;
     pDIFFERENCE(hair->p[1], hair->p[0], d_init);

     for (int i = 0; i < (numPoints - 1); i++)
     {
         d_prior  = ((i - 1) < 0) ? d_init : d[i - 1];
         d_2prior = ((i - 2) < 0) ? d_init : d[i - 2];

         pDIFFERENCE(hair->p[i + 1], hair->p[i], to_next);

         // d_prior now holds first term
         pMULTIPLY(d_prior, 2.0 * (1 - beta), d_prior);
         // d_2prior now holds second term
         pMULTIPLY(d_2prior, -1.0 * (1.0 - beta) * (1.0 - beta), d_2prior);
         // to_next now holds third term
         pMULTIPLY(to_next, beta * beta, to_next);

         pSUM(d_prior, d_2prior, d[i]);
         pSUM(d[i], to_next, d[i]);
     }
    
     for (int i = 0; i < numPoints; i++) {
          if (i == 0) {
            //hair->p[i] = hair->p0[i];
              hair->p_smooth[i] = hair->p0[i];
          } else {
            //pSUM(hair->p[i - 1], d[i - 1], hair->p[i]);
              pSUM(hair->p_smooth[i - 1], d[i - 1], hair->p_smooth[i]);
          }
     }
}

void getFrames(struct world* hair, double frames[numPoints][3][3]){
    struct point up, start, end, aim, cross;
    //initial up
    up.x = 0.0;
    up.y = 0.0;
    up.z = 1.0;
    for (int i = 0; i < (numPoints - 1); i++) {
        pCPY(hair->p_smooth[i],     start); 
        pCPY(hair->p_smooth[i + 1], end);  

        pDIFFERENCE(end, start, aim); 
        double length; 
        pNORMALIZE(aim);

        CROSSPRODUCTp(aim, up, cross);

        pNORMALIZE(cross);

        CROSSPRODUCTp(cross, aim, up);
        pNORMALIZE(up);


        // 1st column is aim; 2nd up; 3rd cross
        frames[i][0][0] = aim.x;
        frames[i][1][0] = aim.y;
        frames[i][2][0] = aim.z;

        frames[i][0][1] = up.x;
        frames[i][1][1] = up.y;
        frames[i][2][1] = up.z;

        frames[i][0][2] = cross.x;
        frames[i][1][2] = cross.y;
        frames[i][2][2] = cross.z;
    }

}

// Takes in matrix m, transposes it, and stores in mT
void transpose(double m[3][3], double mT[3][3]) { // TODO: reference syntax ??

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++) {
            mT[j][i] = m[i][j];
        }
    }
}

void frameTimesVector(double frame[3][3], struct point& vec, struct point& reference) {
     double ref[3];
     ref[0] = 0.0;
     ref[1] = 0.0; 
     ref[2] = 0.0; 
     
     double v[3];
     v[0] = vec.x;
     v[1] = vec.y;
     v[2] = vec.z;

    for (int i = 0; i < 3; i++)
     {
         for (int j = 0; j < 3; j++)
         {
             ref[i] += frame[i][j] * v[j];
         }
     }

     reference.x = ref[0];
     reference.y = ref[1];
     reference.z = ref[2];
}

void getInitialReferenceVectors(struct world* hair){
    struct point edge;
    struct point zero; 
    zero.x = 0.0;
    zero.y = 0.0;
    zero.z = 0.0;
    hair->t0[0] = zero; 
    double Ft[3][3];
    for(int i=1; i< numPoints-1; i++){
        transpose(hair->F0[i - 1], Ft);
        pDIFFERENCE(hair->p0[i + 1], hair->p0[i], edge); 
        frameTimesVector(hair->F0[i - 1], edge, hair->t0[i]);
    }
}

void getReferenceVectors(struct world *hair) {
    struct point zero; 
    zero.x = 0.0;
    zero.y = 0.0;
    zero.z = 0.0;
    pCPY(hair->t[0], zero); 
    pCPY(hair->t[numPoints - 1], zero);
    double Ft[3][3];
    for (int i = 1; i < (numPoints - 1); i++)
    {
        transpose(hair->F[i - 1], Ft);
        frameTimesVector(hair->F[i - 1], hair->t0[i], hair->t[i]);
    }
}

void stretchSpringForce(struct world *hair) {
    // For each point mass, add stretch spring force to its force accumulator and
    // to that of the next point (equal and opposite)
    for (int i = 0; i < (numPoints - 1); i++) {
        // Damped stretch spring equation from paper
        struct point edge;
        struct point rest_edge;
        double length; // length of edge (will be assigned value by pNORMALIZE)
        double curr_length;
        double rest_length;
        // struct point delta_v;
        // double dot; // dot product of delta_v and edge unit vector
        struct point spring_term;
        // struct point damp_term;
        // struct point stretch_force;
        struct point neighbor_force; // equal and opposite of stretch_force, to be
                                     // applied to neighbor                                                                                                                                                                                          
        pDIFFERENCE(hair->p[i+1], hair->p[i], edge);
        pDIFFERENCE(hair->p0[i+1], hair->p0[i], rest_edge);

        pNORMALIZE(edge); // edge now holds edge unit vector
        curr_length = length;
        pNORMALIZE(rest_edge);
        rest_length = length;
        pMULTIPLY(edge, (hair->kStretch * (curr_length - rest_length)), spring_term);
        pMULTIPLY(spring_term, -1.0, neighbor_force);

        pSUM(spring_term, hair->f[i], hair->f[i]);     
        pSUM(neighbor_force, hair->f[i + 1], hair->f[i + 1]);
    }
 }

void bendSpringForce(struct world *hair) {
  // Control bending similar to elastic rods model 
    for (int i = 0; i < (numPoints - 1); i++) {
        // Damped stretch spring equation from paper
        struct point edge;
        double length; // length of edge (will be assigned value by pNORMALIZE)
        struct point delta_v;
        double dot; // dot product of delta_v and edge unit vector
        struct point spring_term;
        
        struct point neighbor_force; 

        pDIFFERENCE(hair->p[i+1], hair->p[i], edge);

        pDIFFERENCE(edge, hair->t[i], spring_term)

        pMULTIPLY(spring_term, hair->kElastic, spring_term);

        pSUM(spring_term, hair->f[i], hair->f[i]);

        pMULTIPLY(spring_term, -1.0, neighbor_force);
        pSUM(neighbor_force, hair->f[i+1], hair->f[i+1]);
    }
}

void integrateForces(struct world* hair) {
    struct point a[numPoints];
    hair->v[0].x = 0.0;
    hair->v[0].y = 0.0;
    hair->v[0].z = 0.0;
    for (int i = 1; i < numPoints; i++) {
         a[i].x = hair->f[i].x / hair->mass;
         a[i].y = hair->f[i].y / hair->mass;
         a[i].z = hair->f[i].z / hair->mass;
         hair->v[i].x += hair->dt * a[i].x;
         hair->v[i].y += hair->dt * a[i].y;
         hair->v[i].z += hair->dt * a[i].z;

         hair->f[i].x = 0;
         hair->f[i].y = 0;
         hair->f[i].z = 0;
     }
}

void dampingForces(struct world *hair, struct point a[numPoints]){
    // First point is stationary
    hair->v[0].x = 0.0;
    hair->v[0].y = 0.0;
    hair->v[0].z = 0.0;

    for (int i = 1; i < numPoints; i++){
        stretchDampingForce(hair); 
        bendDampingForce(hair);
    }
 }

void stretchDampingForce(struct world *hair){
    for (int i = 0; i < (numPoints - 1); i++) {
        // Damped stretch spring equation from paper
        struct point edge;
        struct point rest_edge;
        double length; // length of edge (will be assigned value by pNORMALIZE)
        double curr_length;
        double rest_length;
        struct point delta_v;
        double dot; // dot product of delta_v and edge unit vector
        // struct point spring_term;
        struct point damp_term;
        // struct point stretch_force;
        struct point neighbor_force; // equal and opposite of stretch_force, to be
                                     // applied to neighbor

        pDIFFERENCE(hair->p[i+1], hair->p[i], edge);
        // pDIFFERENCE(hair->p0[i+1], hair->p0[i], rest_edge);

        pNORMALIZE(edge); // edge now holds edge unit vector
        curr_length = length;
        // pNORMALIZE(rest_edge);
        // rest_length = length;
        pDIFFERENCE(hair->v[i + 1], hair->v[i], delta_v);
        pDOT(delta_v, edge, dot);
        pMULTIPLY(edge, (hair->dStretch * dot), damp_term);
        // pMULTIPLY(edge, (hair->kStretch * (rest_length - curr_length)), spring_term);      
        // pSUM(spring_term, damp_term, stretch_force);
        pMULTIPLY(damp_term, -1.0, neighbor_force);

        pSUM(damp_term, hair->f[i], hair->f[i]);     
        pSUM(neighbor_force, hair->f[i + 1], hair->f[i + 1]);
    }
}

void bendDampingForce(struct world *hair){
    for (int i = 0; i < (numPoints - 1); i++) {
        // Damped stretch spring equation from paper
        struct point og_edge; 
        struct point edge;
        double length; // length of edge (will be assigned value by pNORMALIZE)
        struct point delta_v;
        double dot; // dot product of delta_v and edge unit vector
        // struct point spring_term;
        struct point damp_term;
        // struct point bending_force;
        
        struct point neighbor_force; // equal and opposite of stretch_force, to be
                                     // applied to neighbor

        pDIFFERENCE(hair->p[i+1], hair->p[i], og_edge);
        edge = og_edge; 
        pNORMALIZE(edge); // edge now holds edge unit vector
        pDIFFERENCE(hair->v[i + 1], hair->v[i], delta_v);
        pDOT(delta_v, edge, dot);

        pMULTIPLY(edge, dot, damp_term);   
        pDIFFERENCE(delta_v, damp_term, damp_term);
        pMULTIPLY(damp_term, hair->dElastic, damp_term); 

        pSUM(damp_term, hair->f[i], hair->f[i]);

        pMULTIPLY(damp_term, -1.0, neighbor_force);
        pSUM(neighbor_force, hair->f[i+1], hair->f[i+1]);
    }
}

void gravity(struct world *hair) {
     struct point g;
     g.x = 0.0;
     g.y = 0.0;
     g.z = -9.81;

     struct point w; 
     for (int i = 0; i < numPoints; i++)
     {
        pMULTIPLY(g, hair->mass, w);
        pSUM(hair->f[i], w, hair->f[i]); 
     }
}

 void wind(struct world *hair) {
    for (int i = 1; i < numPoints; i++)
    {
      if(up == 1){
            hair->f[i].z += 1000.0; 
        }
      else if (down == 1){
            hair->f[i].z -= 1000.0;
        }
      else if(left == 1){
            hair->f[i].y -= 1000.0; 
        }
      else if(right == 1){
            hair->f[i].y += 1000.0; 
        }
    }
 }

 void collisions(struct world* hair, struct world *collider) {
     // Treat collision between two nodes as undamped stretch spring
     // (But only apply collision force to THIS strand's node! Each strand
     // will be detecting collisions for themselves)
     for (std::pair<int, int> collision : hair->collisions) {
         struct point edge;
         struct point rest_edge;
         double length; // length of edge (will be assigned value by pNORMALIZE)
         double curr_length;
         double rest_length;

         struct point spring_term;

         pDIFFERENCE(collider->p[collision.second],  hair->p[collision.first],  edge);
         pDIFFERENCE(collider->p0[collision.second], hair->p0[collision.first], rest_edge);

         pNORMALIZE(edge); // edge now holds edge unit vector
         curr_length = length;
         pNORMALIZE(rest_edge);
         rest_length = length;

         pMULTIPLY(edge, (hair->kCollision * (rest_length - curr_length)), spring_term);
         pSUM(spring_term, hair->f[collision.first], hair->f[collision.first]);
     }
}

 void coreSpringForce(struct world* hair) {

}

void computeAcceleration(struct world * hair, struct world *collider, struct point a[numPoints])
 {

    //clear f
    for (int i = 0; i < numPoints; i++) {
        hair->f[i].x = 0;
        hair->f[i].y = 0;
        hair->f[i].z = 0;
    }

   // INTEGRATE internal hair forces
   stretchSpringForce(hair);
   smoothing(hair);
   getFrames(hair, hair->F); // Frames (hair->local_frames) are computed along smoothed representation of hair curve (to reduce sensitivity of frame to changes in hair positions)
   getReferenceVectors(hair); // Uses local frames to calculate reference vectors t
   bendSpringForce(hair); // Uses reference vectors to calculate bending spring force and apply to each point and their neighbor
   coreSpringForce(hair);

   //integrateForces(hair);

   // INTEGRATE external forces
   gravity(hair);
   if (hair->pushedByWind) {
       wind(hair);
   }
   collisions(hair, collider);
   //integrateForces(hair);

   for (int i = 0; i < 7; i++) {
       dampingForces(hair, a);
       integrateForces(hair);
   }
 }

  /* performs one step of Inmplecit Euler Integration */
void Euler(struct world * hair, struct world * collider)
{
    int i,j,k;
    struct point a[numPoints];

    // Detect collisions
    struct point displacement;
    double distance;
    hair->collisions.clear();
    for (i = 0; i < numPoints; i++) {
        for (j = 0; j < numPoints; j++) {
            pDIFFERENCE(collider->p[j], hair->p[i], displacement);
            pLENGTH(displacement, distance);
            // This is a collision! Need to push hair ID (NOT COLLIDER) back!
            if (distance < 0.1) {
                hair->collisions.insert({i, j});
            }
        }
    }

    // Force loop
    for (i = 0; i < 10; i++) {
        computeAcceleration(hair, collider, a);
    }
    
    for (i = 1; i < numPoints; i++) {
        hair->p[i].x += hair->dt * hair->v[i].x;
        hair->p[i].y += hair->dt * hair->v[i].y;
        hair->p[i].z += hair->dt * hair->v[i].z;
    }
  }


