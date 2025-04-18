/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include <iostream>

// find neighbour points for structural
void findNeighbours(struct world * jello, int i, int j, int k, struct point neighbours[32], struct point Vneighbours[32]){
  int ip, jp, kp;
  int structural_neighbours[32][3] = { // total 6+20+6=32 neighbours
    // structual
    {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
    {-1, 0, 0}, {0, -1, 0}, {0, 0, -1},
    // shear
    {1, 1, 0}, {-1, 1, 0}, {-1, -1, 0}, {1, -1, 0},
    {0, 1, 1}, {0, -1, 1}, {0, -1, -1}, {0, 1, -1},
    {1, 0, 1}, {-1, 0, 1}, {-1, 0, -1}, {1, 0, -1},
    {1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {1, -1, 1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, -1}, {1, -1, -1},
    // bend
    {2, 0, 0}, {0, 2, 0}, {0, 0, 2},
    {-2, 0, 0}, {0, -2, 0}, {0, 0, -2}
  };
  for (int n=0; n<32; n++){
    ip = i + structural_neighbours[n][0];
    jp = j + structural_neighbours[n][1];
    kp = k + structural_neighbours[n][2];
    if (!(ip>7 || ip<0 || jp>7 || jp<0 || kp>7 || kp<0)){
      neighbours[n] = jello->p[ip][jp][kp];
      Vneighbours[n] = jello->v[ip][jp][kp];
    } else { // when neighbour does not exist
      neighbours[n] = jello->p[i][j][k];
      Vneighbours[n] = jello->v[i][j][k];
    }
  }
}

// compute Hook Force for point A
struct point computeHook(struct world * jello, int i, int j, int k){
  struct point a = jello->p[i][j][k]; // A
  struct point b; // B
  struct point F; F.x = 0; F.y = 0; F.z = 0; // Hook force
  struct point L; // L = A - B
  double length; // length of L
  double kHook = jello->kElastic;
  struct point neighbours[32];
  struct point Vneighbours[32];
  findNeighbours(jello, i, j, k, neighbours, Vneighbours);

  // compute hook force for each neighbour of A
  for (int n=0; n<32; n++){
    b = neighbours[n];
    pDIFFERENCE(a, b, L); // compute L
    // normalize L
    if (L.x == 0 && L.y == 0 && L.z ==0){
      continue;
    } else {
      pNORMALIZE(L);
    }

    // find spring rest length
    double R;
    if (n < 6) {
      R = 1.0 / 7.0;
    } else if (n < 18){
      R = sqrt(2.0) / 7.0;
    } else if (n < 26){
      R = sqrt(3.0) / 7.0;
    } else {
      R = 2.0 / 7.0;
    }

    struct point F1;
    pMULTIPLY(L, (-kHook * (length - R)), F1); // calculate hook
    pSUM(F1, F, F); // sum all hook forces
  }
  return F;
}

// compute Dampling Force for point A
struct point computeDampling(struct world * jello, int i, int j, int k){
  struct point a = jello->p[i][j][k]; // A
  struct point b; // B
  struct point F; F.x=0; F.y=0; F.z=0;// Dampling force
  struct point L; // L = A - B
  struct point Lnorm; // L / |L|
  double length; // |L|
  double kDampling = jello->dElastic;

  // find neighbour
  struct point neighbours[32];
  struct point Vneighbours[32];
  findNeighbours(jello, i, j, k, neighbours, Vneighbours);

  // velocity
  struct point Va = jello->v[i][j][k];
  struct point Vb;
  struct point Vdiff;

  for (int n=0; n<32; n++){
    b = neighbours[n];
    Vb = Vneighbours[n];
    pDIFFERENCE(Va, Vb, Vdiff); // compute Va - Vb
    pDIFFERENCE(a, b, L); // compute L = a - b
    
    // skip if neighbour does not exist
    if (L.x == 0 && L.y == 0 && L.z ==0){
      continue;
    }

    pCPY(L, Lnorm) 
    pNORMALIZE(Lnorm); // compute |L|
    double dotproduct = -kDampling * ((Vdiff.x * L.x + Vdiff.y * L.y + Vdiff.z * L.z) / length);
    struct point F1;
    pMULTIPLY(Lnorm, dotproduct, F1);
    pSUM(F1, F, F);
  }
  return F;
}

struct point computeCollisionDetection(struct world * jello, int i, int j, int k){
  struct point a = jello->p[i][j][k];
  struct point v = jello->v[i][j][k];
  double kCollision = jello->kCollision;
  double dCollision = jello->dCollision;
  double penetration;
  struct point normal;
  struct point Fhook, Fdamping, Ftotal;
  double length;

  // compute Hook Force
  if (a.x <= -2){
    pMAKE(1, 0, 0, normal);
    penetration = -2 - a.x;
  } else if (a.x >= 2) {
    pMAKE(-1, 0, 0, normal);
    penetration = a.x - 2;
  } else if (a.y <= -2) {
    pMAKE(0, 1, 0, normal);
    penetration = -2 - a.y;
  } else if (a.y >= 2) {
    pMAKE(0, -1, 0, normal);
    penetration = a.y - 2;
  } else if (a.z <= -2) {
    pMAKE(0, 0, 1, normal);
    penetration = -2 - a.z;
  } else if (a.z >= 2) {
    pMAKE(0, 0, -1, normal);
    penetration = a.z - 2;
  }
  pMULTIPLY(normal, (kCollision * penetration), Fhook);

  // compute Damping Force
  pMULTIPLY(normal, (dCollision * (v.x * normal.x + v.y * normal.y + v.z * normal.z)), Fdamping);

  // compute total collision detection force
  pSUM(Fhook, Fdamping, Ftotal);
  return Ftotal;
}

struct point computePlaneCollision(struct world * jello, int i, int j, int k)
{
    struct point p = jello->p[i][j][k];
    struct point v = jello->v[i][j][k];

    double A = jello->a, B = jello->b, C = jello->c, D = jello->d;
    // plane eqn: A*x + B*y + C*z + D = 0
    double val = A*p.x + B*p.y + C*p.z + D;

    // if val >= 0, no collision
    if (val >= 0) {
      struct point zero={0,0,0}; 
      return zero;
    }

    // otherwise, we have negative => behind plane => collision
    double length = sqrt(A*A + B*B + C*C);
    double penetration = -(val)/length;  // negative val => positive pen
    struct point normal = {A/length, B/length, C/length};

    // collision spring
    struct point Fhook;
    pMULTIPLY(normal, jello->kCollision * penetration, Fhook);

    // collision damping
    double normalVel = (v.x * normal.x + v.y * normal.y + v.z * normal.z);
    struct point Fdamp;
    pMULTIPLY(normal, jello->dCollision*normalVel, Fdamp);

    // sum
    struct point Ftotal;
    pSUM(Fhook, Fdamp, Ftotal);
    return Ftotal;
}

struct point computeSphereCollision(struct world * jello, int i, int j, int k) {
  struct point Fsphere;
  Fsphere.x = 0.0; Fsphere.y = 0.0; Fsphere.z = 0.0;
  
  double kCol = jello->kCollision;
  double dCol = jello->dCollision;
  double px = jello->p[i][j][k].x;
  double py = jello->p[i][j][k].y;
  double pz = jello->p[i][j][k].z;
  double vx = jello->v[i][j][k].x;
  double vy = jello->v[i][j][k].y;
  double vz = jello->v[i][j][k].z;
  
  // distance from the center
  double d = sqrt(pow((px - 1.0),2) + pow((py - 1.5),2) + pow((pz - 1.5),2));
  
  if (d < 1e-6) return Fsphere;

  if (d < 0.4) {
    // normal
    struct point normal;
    normal.x = (px - 1.0) / d;
    normal.y = (py - 1.5) / d;
    normal.z = (pz - 1.5) / d;

    // Hook force
    struct point Fhook;
    pMULTIPLY(normal, kCol*(0.5 - d), Fhook);

    // Damping force
    // velocity dot normal
    double normalVel = (vx * normal.x + vy * normal.y + vz * normal.z);
    struct point Fdamp;
    pMULTIPLY(normal, dCol * normalVel, Fdamp);

    // total
    pSUM(Fhook, Fdamp, Fsphere);
  }
  
  return Fsphere;
}

// detect collision detection and compute force if collide
bool detectCollisionDetection(struct world * jello, int i, int j, int k){
  struct point a = jello->p[i][j][k];
  if (a.x <= -2 || a.x >= 2 
  || a.y <= -2 || a.y >= 2
  || a.z <= -2 || a.z >= 2){
    return true;
  }
  return false;
}

// I used Trilinear Interpolation, equations are referenced from https://spie.org/samples/PM159.pdf
struct point computeForceField(struct world * jello, int i, int j, int k){
  struct point a = jello->p[i][j][k];
  int n = jello->resolution;
  struct point * forceField = jello->forceField;

  double x, y, z; // actual space location
  x = a.x;
  y = a.y;
  z = a.z;
  if (x < -2.0) { x = -2.0; } else if (x > 2.0) { x = 2.0; } // avoid out of bounding box
  if (y < -2.0) { y = -2.0; } else if (y > 2.0) { y = 2.0; }
  if (z < -2.0) { z = -2.0; } else if (z > 2.0) { z = 2.0; }

  double u, v, w; // index
  u = ((x + 2.0) * (n - 1.0)) / 4.0;
  v = ((y + 2.0) * (n - 1.0)) / 4.0;
  w = ((z + 2.0) * (n - 1.0)) / 4.0;
  int u0 = floor(u);
  int u1 = ceil(u);
  int v0 = floor(v);
  int v1 = ceil(v);
  int w0 = floor(w);
  int w1 = ceil(w);
  
  struct point p000 = jello->forceField[u0*n*n + v0*n + w0];
  struct point p100 = jello->forceField[u1*n*n + v0*n + w0];
  struct point p110 = jello->forceField[u1*n*n + v1*n + w0];
  struct point p010 = jello->forceField[u0*n*n + v1*n + w0];
  struct point p001 = jello->forceField[u0*n*n + v0*n + w1];
  struct point p011 = jello->forceField[u0*n*n + v1*n + w1];
  struct point p111 = jello->forceField[u1*n*n + v1*n + w1];
  struct point p101 = jello->forceField[u1*n*n + v0*n + w1];

  struct point c0 = p000;
  struct point c1; pDIFFERENCE(p100, p000, c1);
  struct point c2; pDIFFERENCE(p010, p000, c2);
  struct point c3; pDIFFERENCE(p001, p000, c3);
  struct point c4;
  c4.x = p110.x - p010.x - p100.x + p000.x;
  c4.y = p110.y - p010.y - p100.y + p000.y;
  c4.z = p110.z - p010.z - p100.z + p000.z;
  struct point c5;
  c5.x = p011.x - p001.x - p010.x + p000.x;
  c5.y = p011.y - p001.y - p010.y + p000.y;
  c5.z = p011.z - p001.z - p010.z + p000.z;
  struct point c6;
  c6.x = p101.x - p001.x - p100.x + p000.x;
  c6.y = p101.y - p001.y - p100.y + p000.y;
  c6.z = p101.z - p001.z - p100.z + p000.z;
  struct point c7;
  c7.x = p111.x - p011.x - p101.x - p110.x + p100.x + p001.x + p010.x - p000.x;
  c7.y = p111.y - p011.y - p101.y - p110.y + p100.y + p001.y + p010.y - p000.y;
  c7.z = p111.z - p011.z - p101.z - p110.z + p100.z + p001.z + p010.z - p000.z;

  double u_delta, v_delta, w_delta; 
  if (u1 == u0) {
    u_delta = 0.0;
  } else {
    u_delta = (u - u0) / (u1 - u0);
  }
  if (v1 == v0) {
    v_delta = 0.0;
  } else {
    v_delta = (v - v0) / (v1 - v0);
  }
  if (w1 == w0) {
    w_delta = 0.0;
  } else {
    w_delta = (w - w0) / (w1 - w0);
  }
  
  struct point F;
  F.x = c0.x + c1.x * u_delta + c2.x * v_delta + c3.x * w_delta + c4.x * u_delta * v_delta + c5.x * v_delta * w_delta + c6.x * w_delta * u_delta + c7.x * u_delta * v_delta * w_delta;
  F.y = c0.y + c1.y * u_delta + c2.y * v_delta + c3.y * w_delta + c4.y * u_delta * v_delta + c5.y * v_delta * w_delta + c6.y * w_delta * u_delta + c7.y * u_delta * v_delta * w_delta;
  F.z = c0.z + c1.z * u_delta + c2.z * v_delta + c3.z * w_delta + c4.z * u_delta * v_delta + c5.z * v_delta * w_delta + c6.z * w_delta * u_delta + c7.z * u_delta * v_delta * w_delta;

  return F;
}

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  int i,j,k;
  // notice we are missing Fforcefield right now !!!
  double m = jello->mass;
  struct point acc;
  struct point Ftotal, Fhook, Fdampling, Ffield, Fcollision;
  for (i=0; i<=7; i++){
    for (j=0; j<=7; j++){
      for (k=0; k<=7; k++){
        Fhook = computeHook(jello, i, j, k);
        Fdampling = computeDampling(jello, i, j, k);
        pSUM(Fhook, Fdampling, Ftotal);

        Ffield = computeForceField(jello, i, j, k);
        pSUM(Ffield, Ftotal, Ftotal);

        // compute collision detection
        if (detectCollisionDetection(jello, i, j, k)) {
          Fcollision = computeCollisionDetection(jello, i, j, k);
          pSUM(Fcollision, Ftotal, Ftotal);
        }

        if (jello->incPlanePresent == 1) {
          struct point Fplane = computePlaneCollision(jello, i, j, k);
          pSUM(Fplane, Ftotal, Ftotal);
        }

        if (spherePresent == 1) {
          struct point Fsphere = computeSphereCollision(jello, i, j , k);
          pSUM(Fsphere, Ftotal, Ftotal);
        }

        if (grav == 1) {
          struct point Fgrav;
          pMAKE(0.0, 0.0, -0.981, Fgrav);
          pSUM(Fgrav, Ftotal, Ftotal);
        }

        pMULTIPLY(Ftotal, 1/m, acc);
        a[i][j][k] = acc;
      }
    }
  }
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{ 
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
