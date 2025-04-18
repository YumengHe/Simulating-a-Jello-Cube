/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include "IPC.h"
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <cassert>
#include <algorithm>

// Eigen
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Sparse"
#include "eigen-3.4.0/Eigen/SparseCholesky"

using namespace std;
using namespace Eigen;

static const struct point gravity = {0.0, 0.0, -9.81};
static const double dhat = 0.01;
static const double kappa = 1e5;
static const double BMIN = -2.0;
static const double BMAX =  2.0;


// math operations
struct point vec_minus_vec(struct point a, struct point b) {
  struct point result;
  result.x = a.x - b.x;
  result.y = a.y - b.y;
  result.z = a.z - b.z;
  return result;
}

struct point vec_plus_vec(struct point a, struct point b) {
  struct point result;
  result.x = a.x + b.x;
  result.y = a.y + b.y;
  result.z = a.z + b.z;
  return result;
}

struct point vec_time_scale(struct point a, double b) {
  struct point result;
  result.x = a.x * b;
  result.y = a.y * b;
  result.z = a.z * b;
  return result;
}

double vec_dot_vec(struct point a, struct point b) {
  double result;
  result = a.x * b.x + a.y * b.y + a.z * b.z;
  return result;
}

// find index
int nodeIndex(int i, int j, int k) {
  return (i * 8 * 8) + (j * 8) + k;
}

int indexI(int x) {
  return int(x / 64);
}

int indexJ(int x) {
  return int((x / 8) % 8);
}

int indexK(int x) {
  return int(x % 8);
}

// helper functions
void copyJelloPositions(world *jello, vector<point> &x) {
  x.resize(512);
  for(int i = 0; i < 8; i++){
    for(int j = 0; j < 8; j++){
      for(int k = 0; k < 8; k++){
        x[nodeIndex(i,j,k)] = jello->p[i][j][k];
      }
    }
  }
}

// compute edges
vector< array<int, 2> > generateJelloEdges() {
  vector< array<int, 2> > edges;
  edges.reserve(6000);

  // ------ structural edges ------
  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      for(int k = 0; k < 8; k++) {
        int idx = nodeIndex(i, j, k);
        // along x
        if(i + 1 < 8) {
          array<int,2> edgeTmp = {idx, nodeIndex(i+1, j, k)};
          edges.push_back(edgeTmp);
        }
        // along y
        if(j + 1 < 8) {
          array<int,2> edgeTmp = {idx, nodeIndex(i, j+1, k)};
          edges.push_back(edgeTmp);
        }
        // along z
        if(k + 1 < 8) {
          array<int,2> edgeTmp = {idx, nodeIndex(i, j, k+1)};
          edges.push_back(edgeTmp);
        }
      }
    }
  }

  // ------ shear edges ------
  int shearOffsets[20][3] = {
    {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1,-1, 0},
    {0, 1, 1}, {0, 1,-1}, {0,-1, 1}, {0, -1, -1},
    {1, 0, 1}, {1, 0,-1}, {-1,0, 1}, {-1,0,-1},
    {1, 1, 1}, {1, 1,-1}, {1,-1, 1}, {-1,1, 1}, { -1,-1,1}, { -1,1,-1}, { 1,-1,-1}, { -1,-1,-1}
  };

  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      for(int k = 0; k < 8; k++) {
        int idx = nodeIndex(i, j, k);
        // Try each offset in the shearOffsets list
        for(auto &ofs : shearOffsets) {
            int ni = i + ofs[0]; // neighbor i
            int nj = j + ofs[1]; // neighbor j
            int nk = k + ofs[2]; // neighbor k

          // check boundaries
          if(ni >= 0 && ni < 8 &&
              nj >= 0 && nj < 8 &&
              nk >= 0 && nk < 8) {
            // To avoid duplicating edges in both directions,
            // only push if (ni, nj, nk) > (i, j, k) in flattened order.
            // That is, we only add edge if new_idx > idx.
            int new_idx = nodeIndex(ni, nj, nk);
            if(new_idx > idx)
            {
              array<int,2> edgeTmp = {idx, new_idx};
              edges.push_back(edgeTmp);
            }
          }
        }
      }
    }
  }

  // ------ bend edges ------
  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      for(int k = 0; k < 8; k++) {
        int idx = nodeIndex(i, j, k);
        // i+2
        if(i + 2 < 8) {
          array<int,2> edgeTmp = {idx, nodeIndex(i+2, j, k)};
          edges.push_back(edgeTmp);
        }
        // j+2
        if(j + 2 < 8) {
          array<int,2> edgeTmp = {idx, nodeIndex(i, j+2, k)};
          edges.push_back(edgeTmp);
        }
        // k+2
        if(k + 2 < 8) {
          array<int,2> edgeTmp = {idx, nodeIndex(i, j, k+2)};
          edges.push_back(edgeTmp);
        }
      }
    }
  }

  return edges;
}

// Inertia Term
double valInertia(world *jello, vector<point> &xTilde, double m) {
  double sum = 0.0;
  for (int i = 0; i < 512; i++) {
    point diff = vec_minus_vec(jello->p[indexI(i)][indexJ(i)][indexK(i)], xTilde[i]);
    sum += 0.5 * m * vec_dot_vec(diff, diff);
  }
  return sum;
}

double temp_valInertia(vector<point> &x, vector<point> &xTilde, double m)
{
  double sum = 0.0;
  for (int i = 0; i < 512; i++) {
    point diff = vec_minus_vec(x[i], xTilde[i]);
    sum += 0.5 * m * vec_dot_vec(diff, diff);
  }
  return sum;
}

vector<point> gradInertia(world *jello, vector<point> &xTilde, double m) {
  vector<point> g(512);
  for (int i = 0; i < 512; i++) {
    point diff = vec_minus_vec(jello->p[indexI(i)][indexJ(i)][indexK(i)], xTilde[i]);
    g[i] = vec_time_scale(diff, m);
  }
  return g;
}

IJV hessInertia(world *jello, vector<point> &xTilde, double m) {
  IJV h;
  h.I.resize(3 * 512);
  h.J.resize(3 * 512);
  h.V.resize(3 * 512);

  for (int i = 0; i < 512; i++) {
    for (int d = 0; d < 3; d++) {
      h.I[i * 3 + d] = i * 3 + d;
      h.J[i * 3 + d] = i * 3 + d;
      h.V[i * 3 + d] = m;
    }
  }
  return h;
}

// Mass-Spring Potential Energy
double valSpring(world *jello, vector< array<int,2> > &e, vector<double> &l2, double k)
{
  double sum = 0.0;
  for(int i = 0; i < e.size(); i++) {
    point diff = vec_minus_vec(jello->p[indexI(e[i][0])][indexJ(e[i][0])][indexK(e[i][0])], jello->p[indexI(e[i][1])][indexJ(e[i][1])][indexK(e[i][1])]);
    sum += l2[i] * 0.5 * k * pow((vec_dot_vec(diff, diff) / l2[i] - 1), 2);
  }
  return sum;
}

double temp_valSpring(vector<point> &x, vector< array<int,2> > &e, vector<double> &l2, double k)
{
  double sum = 0.0;
  for(int i = 0; i < e.size(); i++) {
    point diff = vec_minus_vec(x[e[i][0]], x[e[i][1]]);
    sum += l2[i] * 0.5 * k * pow((vec_dot_vec(diff, diff) / l2[i] - 1), 2);
  }
  return sum;
}

vector<point> gradSpring(world *jello, vector< array<int,2> > &e, vector<double> &l2, double k) {
  vector<point> g(512);
    
  // initialize
  for (int i = 0; i < 512; i++){
    pMAKE(0.0, 0.0, 0.0, g[i]);
  }

  for(int i = 0; i < e.size(); i++) {
    point diff = vec_minus_vec(jello->p[indexI(e[i][0])][indexJ(e[i][0])][indexK(e[i][0])], jello->p[indexI(e[i][1])][indexJ(e[i][1])][indexK(e[i][1])]);
    point g_diff = vec_time_scale(diff, 2.0 * k * (vec_dot_vec(diff, diff) / l2[i] - 1));
    g[e[i][0]] = vec_plus_vec(g[e[i][0]], g_diff);
    g[e[i][1]] = vec_minus_vec(g[e[i][1]], g_diff);
  }

  return g;
}

IJV hessSpring(world *jello, vector< array<int,2> > &e, vector<double> &l2, double k) {
  // Each edge contributes 6x6 = 36 entries to the global Hessian
  IJV h;
  h.I.resize(e.size() * 36);
  h.J.resize(e.size() * 36);
  h.V.resize(e.size() * 36);

  for(int i = 0; i < e.size(); i++) {
    // compute diff
    point diff = vec_minus_vec(jello->p[indexI(e[i][0])][indexJ(e[i][0])][indexK(e[i][0])], jello->p[indexI(e[i][1])][indexJ(e[i][1])][indexK(e[i][1])]);
    
    // compute H_diff
    double factor = 2.0 * k / l2[i];
    double Hdiff[9];
    Hdiff[0] = 2.0 * diff.x * diff.x + vec_dot_vec(diff, diff) - l2[i];
    Hdiff[1] = 2.0 * diff.x * diff.y;
    Hdiff[2] = 2.0 * diff.x * diff.z;
    Hdiff[3] = 2.0 * diff.y * diff.x;
    Hdiff[4] = 2.0 * diff.y * diff.y + vec_dot_vec(diff, diff) - l2[i];
    Hdiff[5] = 2.0 * diff.y * diff.z;
    Hdiff[6] = 2.0 * diff.z * diff.x;
    Hdiff[7] = 2.0 * diff.z * diff.y;
    Hdiff[8] = 2.0 * diff.z * diff.z + vec_dot_vec(diff, diff) - l2[i];
    for(int r=0; r<9; r++) {
      Hdiff[r] *= factor;
    }

    // 6x6 block = [[ Hdiff, -Hdiff ],
    //              [ -Hdiff, Hdiff  ]]
    int baseIdx = i * 36; // each edge uses 36 entries
    int n0 = e[i][0];
    int n1 = e[i][1];
    // For local row R in [0..5], local col C in [0..5]:
    //   global row = (R<3 ? n0 : n1)*3 + (R%3)
    //   global col = (C<3 ? n0 : n1)*3 + (C%3)
    //   sign = +1 if (R<3 and C<3) or (R>=3 and C>=3), else -1
    for(int R=0; R<6; R++) {
      for(int C=0; C<6; C++) {
        int nodePartR = (R < 3) ? n0 : n1;
        int nodePartC = (C < 3) ? n0 : n1;

        int coordR = R % 3;
        int coordC = C % 3;

        double sign = ((R<3 && C<3) || (R>=3 && C>=3)) ? 1.0 : -1.0;
        double valRC = sign * Hdiff[coordR*3 + coordC];

        int globalRow = nodePartR * 3 + coordR;
        int globalCol = nodePartC*3 + coordC;

        int offset = baseIdx + (R*6 + C);
        h.I[offset] = globalRow;
        h.J[offset] = globalCol;
        h.V[offset] = valRC;
      }
    }
  }
  return h;
}

// Gravity Energy
double valGravity(world *jello, double m) {
  double sum = 0.0;
  for (int i = 0; i < 512; i++) {
    sum += -m * vec_dot_vec(jello->p[indexI(i)][indexJ(i)][indexK(i)], gravity);
  }
  return sum;
}

double temp_valGravity(vector<point> &x, double m) {
  double sum = 0.0;
  for (int i = 0; i < 512; i++) {
    sum += -m * vec_dot_vec(x[i], gravity);
  }
  return sum;
}

vector<point> gradGravity(world *jello, double m) {
  vector<point> g(512);
  for (int i = 0; i < 512; i++) {
    g[i].x = -m * gravity.x;
    g[i].y = -m * gravity.y;
    g[i].z = -m * gravity.z;
  }
  return g;
}

// Barrier Energy
double barrierVal(world *jello, double y_ground, double contact_area) {
  double sumVal = 0.0;
  for(size_t i=0; i < 512; i++)
  { 
    double px = jello->p[indexI(i)][indexJ(i)][indexK(i)].x;
    double py = jello->p[indexI(i)][indexJ(i)][indexK(i)].y;
    double pz = jello->p[indexI(i)][indexJ(i)][indexK(i)].z;
    double d[6];
    d[0] = (BMAX - px);
    d[1] = (px - BMIN);
    d[2] = (BMAX - py);
    d[3] = (py - BMIN);
    d[4] = (BMAX - pz);
    d[5] = (pz - BMIN);
    for(int k=0; k<6; k++) {
      if(d[k] < dhat)
      {
          double s = d[k] / dhat;
          sumVal += contact_area * dhat * (0.5*kappa) * (s - 1.0)*log(s);
      }
    }
  }
  return sumVal;
}

double temp_barrierVal(vector<point> &x, double y_ground, double contact_area) {
  double sumVal = 0.0;
  for(size_t i=0; i < 512; i++)
  { 
    double px = x[i].x;
    double py = x[i].y;
    double pz = x[i].z;
    double d[6];
    d[0] = (BMAX - px);
    d[1] = (px - BMIN);
    d[2] = (BMAX - py);
    d[3] = (py - BMIN);
    d[4] = (BMAX - pz);
    d[5] = (pz - BMIN);
    for(int k=0; k<6; k++) {
      if(d[k] < dhat)
      {
          double s = d[k] / dhat;
          sumVal += contact_area * dhat * (0.5*kappa) * (s - 1.0)*log(s);
      }
    }
  }
  return sumVal;
}

vector<point> barrierGrad(world *jello, double y_ground, double contact_area) {
  vector<point> g(512);
  for(size_t i=0; i<512; i++){
      g[i].x = 0.0; g[i].y = 0.0; g[i].z = 0.0;
  }

  for(size_t i=0; i < 512; i++)
  { 
    double px = jello->p[indexI(i)][indexJ(i)][indexK(i)].x;
    double py = jello->p[indexI(i)][indexJ(i)][indexK(i)].y;
    double pz = jello->p[indexI(i)][indexJ(i)][indexK(i)].z;
    double d[6];
    // plane x=+2 => d[0], normal is -x
    d[0] = BMAX - px;
    // plane x=-2 => d[1], normal is +x
    d[1] = px - BMIN;

    // plane y=+2 => d[2], normal is -y
    d[2] = BMAX - py;
    // plane y=-2 => d[3], normal is +y
    d[3] = py - BMIN;

    // plane z=+2 => d[4], normal is -z
    d[4] = BMAX - pz;
    // plane z=-2 => d[5], normal is +z
    d[5] = pz - BMIN;
    for(int k=0; k<6; k++)
        {
            if(d[k] < dhat)
            {
                double s = d[k]/dhat;
                // from the 2D code: 
                // grad[i][normal_dir] = contact_area[i]*dhat*(kappa/2)*( log(s)/dhat + (s-1)/d )
                // We must figure out the sign to get the correct direction.

                // Let's define the partial formula as:
                double val = contact_area*dhat*(0.5*kappa)*( std::log(s)/dhat + (s-1.0)/d[k] );

                // Then apply it in the direction of the plane's normal
                // plane 0 => normal = -X => g[i].x -= val
                // plane 1 => normal = +X => g[i].x += val
                // plane 2 => normal = -Y => g[i].y -= val
                // plane 3 => normal = +Y => g[i].y += val
                // plane 4 => normal = -Z => g[i].z -= val
                // plane 5 => normal = +Z => g[i].z += val
                switch(k)
                {
                    case 0: g[i].x -= val; break; // x=+2
                    case 1: g[i].x += val; break; // x=-2
                    case 2: g[i].y -= val; break; // y=+2
                    case 3: g[i].y += val; break; // y=-2
                    case 4: g[i].z -= val; break; // z=+2
                    case 5: g[i].z += val; break; // z=-2
                }
            }
        }
  }
  return g;
}

IJV barrierHess(world *jello, double y_ground, double contact_area) {
  IJV out;
  out.I.resize(512);
  out.J.resize(512);
  out.V.resize(512);

  for(int i=0; i < 512; i++)
  { 
    double px = jello->p[indexI(i)][indexJ(i)][indexK(i)].x;
        double py = jello->p[indexI(i)][indexJ(i)][indexK(i)].y;
        double pz = jello->p[indexI(i)][indexJ(i)][indexK(i)].z;

        double d[6];
        d[0] = BMAX - px; // plane x=+2 => normal -X
        d[1] = px - BMIN; // plane x=-2 => normal +X
        d[2] = BMAX - py; // plane y=+2 => normal -Y
        d[3] = py - BMIN; // plane y=-2 => normal +Y
        d[4] = BMAX - pz; // plane z=+2 => normal -Z
        d[5] = pz - BMIN; // plane z=-2 => normal +Z

        for(int k=0; k<6; k++)
        {
            if(d[k]<dhat)
            {
                // from 2D code: H(i) = contact_area[i]*dhat*kappa / (2*d*d*dhat)*(d + dhat)
                // We'll store it in the correct dof (3*i+0, or +1, or +2).
                double val = contact_area*dhat*kappa/(2.0*d[k]*d[k]*dhat)*(d[k]+dhat);

                // which dof?
                int dof = -1;
                // plane 0 => x dimension => dof=3*i+0
                // plane 1 => x dimension => dof=3*i+0
                // plane 2 => y dimension => dof=3*i+1
                // plane 3 => y dimension => dof=3*i+1
                // plane 4 => z dimension => dof=3*i+2
                // plane 5 => z dimension => dof=3*i+2
                int coord = 0; // x=0, y=1, z=2
                if(k==0 || k==1) coord=0;
                else if(k==2 || k==3) coord=1;
                else if(k==4 || k==5) coord=2;

                int globalRow = 3*static_cast<int>(i)+coord;
                int globalCol = globalRow;

                // We place it as a diagonal entry (globalRow, globalCol)
                out.I.push_back(globalRow);
                out.J.push_back(globalCol);
                out.V.push_back(val);
            }
        }
  }

  return out;
}

double init_step_size(world *jello, double y_ground, vector<point> &p) {
  double alpha = 1.0; // default
  for(size_t i=0; i<512; i++) {
    double py = p[i].y;
    if(py < 0.0 && jello->p[indexI(i)][indexJ(i)][indexK(i)].y > y_ground)
    {
      double candidate = 0.9 * (y_ground - jello->p[indexI(i)][indexJ(i)][indexK(i)].y) / py;
      if(candidate < alpha && candidate > 0.0)
        alpha = candidate;
    }
  }
  return alpha;
}

// Optimization Time Integrator
double IP_val(world *jello, vector<point> &x_tilde, double m, vector< array<int,2> > &e, vector<double> &l2, double k, double h, double y_ground, double contact_area) {
  double result;
  result = valInertia(jello, x_tilde, m) + h * h * (valSpring(jello, e, l2, k) + barrierVal(jello, y_ground, contact_area));

  if (grav == 1) {
    result += (h * h * valGravity(jello, m));
  }
  return result;
}

double temp_IP_val(vector<point> &x, vector<point> &x_tilde, double m, vector< array<int,2> > &e, vector<double> &l2, double k, double h, double y_ground, double contact_area) {
  double result;
  result = temp_valInertia(x, x_tilde, m) + h * h * (temp_valSpring(x, e, l2, k) + temp_barrierVal(x, y_ground, contact_area));

  if (grav == 1) {
    result += (h * h * temp_valGravity(x, m));
  }
  return result;
}

vector<point> IP_grad(world *jello, vector<point> &x_tilde, double m, vector< array<int,2> > &e, vector<double> &l2, double k, double h, double y_ground, double contact_area) {
  vector<point> gI = gradInertia(jello, x_tilde, m);
  vector<point> gS = gradSpring(jello, e, l2, k);
  vector<point> gB = barrierGrad(jello, y_ground, contact_area);
  
  for(int i=0; i<gI.size(); i++) {
    gI[i].x += (h * h * (gS[i].x + gB[i].x));
    gI[i].y += (h * h * (gS[i].y + gB[i].y));
    gI[i].z += (h * h * (gS[i].z + gB[i].z));
  }

  if (grav == 1) {
    vector<point> gG = gradGravity(jello, m);
    for(int i=0; i<gI.size(); i++) {
      gI[i].x += (h * h * gG[i].x);
      gI[i].y += (h * h * gG[i].y);
      gI[i].z += (h * h *gG[i].z);
    }
  }

  return gI;
}

IJV IP_hess(world *jello, vector<point> &x_tilde, double m, vector< array<int,2> > &e, vector<double> &l2, double k, double h, double y_ground, double contact_area) {
  IJV hIn = hessInertia(jello, x_tilde, m);
  IJV hSp = hessSpring(jello, e, l2, k);
  IJV Hb = barrierHess(jello, y_ground, contact_area);

  for(size_t i=0; i < hSp.V.size(); i++){
      hSp.V[i] *= (h * h);
  }
  for(size_t i=0; i < Hb.V.size(); i++){
      Hb.V[i] *= (h * h);
  }
  
  // combine hIn and hSp
  IJV combined;
  combined.I.reserve(hIn.I.size() + hSp.I.size() + Hb.I.size());
  combined.J.reserve(hIn.J.size() + hSp.J.size() + Hb.J.size());
  combined.V.reserve(hIn.V.size() + hSp.V.size() + Hb.V.size());

  // inertia parts
  combined.I.insert(combined.I.end(), hIn.I.begin(), hIn.I.end());
  combined.J.insert(combined.J.end(), hIn.J.begin(), hIn.J.end());
  combined.V.insert(combined.V.end(), hIn.V.begin(), hIn.V.end());

  // spring parts
  combined.I.insert(combined.I.end(), hSp.I.begin(), hSp.I.end());
  combined.J.insert(combined.J.end(), hSp.J.begin(), hSp.J.end());
  combined.V.insert(combined.V.end(), hSp.V.begin(), hSp.V.end());

  // barrier parts
  combined.I.insert(combined.I.end(), Hb.I.begin(), Hb.I.end());
  combined.J.insert(combined.J.end(), Hb.J.begin(), Hb.J.end());
  combined.V.insert(combined.V.end(), Hb.V.begin(), Hb.V.end());

  return combined;
}

vector<double> solveSparseSystemEigen(IJV &A, vector<double> &rhs, int n)
{
  // convert IJV triplets to Eigen triplets
  vector< Triplet<double> > triplets;
  triplets.reserve(A.I.size());
  for(size_t k = 0; k < A.I.size(); k++) {
    triplets.push_back(Triplet<double>(A.I[k], A.J[k], A.V[k]));
  }

  // 2) Build the sparse matrix
  SparseMatrix<double> mat(n, n);
  mat.setFromTriplets(triplets.begin(), triplets.end());

  // 3) Build the Eigen vector for RHS
  VectorXd b(n);
  for(int i = 0; i < n; i++) {
    b[i] = rhs[i];
  }

  // 4) Choose and run a sparse direct solver.
  //    For a typical (semi)definite Hessian in a physics problem,
  //    SimplicialLDLT or SimplicialLLT is common. If not sure, you
  //    can use SparseLU or SparseQR, but they might be slower.
  SimplicialLDLT< SparseMatrix<double> > solver;
  solver.compute(mat);
  if(solver.info() != Success) {
    return vector<double>(n, 0.0);
  }

  VectorXd x = solver.solve(b);
  if(solver.info() != Success) {
    return std::vector<double>(n, 0.0);
  }

  // 5) Convert x back to std::vector<double>
  vector<double> xSol(n);
  for(int i = 0; i < n; i++) {
    xSol[i] = x[i];
  }
  return xSol;
}

vector<point> search_dir(world *jello, vector<point> &x_tilde, double m, vector< array<int,2> > &e, vector<double> &l2, double k, double h, double y_ground, double contact_area) {
  IJV A = IP_hess(jello, x_tilde, m, e, l2, k, h, y_ground, contact_area);
  vector<point> g = IP_grad(jello, x_tilde, m, e, l2, k, h, y_ground, contact_area);

  // Flatten g into -g for the RHS
  int N = g.size();
  vector<double> rhs(3*N);
  for(int i=0; i < N; i++) {
    rhs[3*i + 0] = -g[i].x;
    rhs[3*i + 1] = -g[i].y;
    rhs[3*i + 2] = -g[i].z;
  }

  // Solve the system: A * p = rhs
  vector<double> p_sol = solveSparseSystemEigen(A, rhs, 3*N);
  // Unflatten p_sol into a vector<point>
  vector<point> p(N);
  for(int i=0; i < N; i++){
    p[i].x = p_sol[3*i + 0];
    p[i].y = p_sol[3*i + 1];
    p[i].z = p_sol[3*i + 2];
  }
  return p;
}

// Simulation
double normInfOverH(vector<point> &p, double h) {
  double maxVal = 0.0;
  for(int i = 0; i < p.size(); i++) {
    double absx = fabs(p[i].x);
    double absy = fabs(p[i].y);
    double absz = fabs(p[i].z);
    if (absx > maxVal) {
      maxVal = absx;
    }
    if (absy > maxVal) {
      maxVal = absy;
    }
    if (absz > maxVal) {
      maxVal = absz;
    }
  }
  return maxVal / h;
};

void stepForwardImplicitEuler(world *jello, vector< array<int,2> > &edges, double m, vector<double> &l2, double k, double h, double tol, vector<point> &x_old, double y_ground, double contact_area) {

  // x_tilde
  vector<point> x_tilde(512);
  for(size_t i=0; i < 512; i++) {
    x_tilde[i].x = jello->p[indexI(i)][indexJ(i)][indexK(i)].x + h * jello->v[indexI(i)][indexJ(i)][indexK(i)].x;
    x_tilde[i].y = jello->p[indexI(i)][indexJ(i)][indexK(i)].y + h * jello->v[indexI(i)][indexJ(i)][indexK(i)].y;
    x_tilde[i].z = jello->p[indexI(i)][indexJ(i)][indexK(i)].z + h * jello->v[indexI(i)][indexJ(i)][indexK(i)].z;
  }

  // newton loop
  int iter = 0;
  double E_last = IP_val(jello, x_tilde, m, edges, l2, k, h, y_ground, contact_area);
  vector<point> p = search_dir(jello, x_tilde, m, edges, l2, k, h, y_ground, contact_area);

  while(normInfOverH(p, h) > tol && iter < 30) {
    
    // line search
    double alpha = init_step_size(jello, y_ground, p);
    while(true) {
      vector<point> tempX(512);
      for(size_t i=0; i < 512; i++){
        tempX[i].x = jello->p[indexI(i)][indexJ(i)][indexK(i)].x + alpha * p[i].x;
        tempX[i].y = jello->p[indexI(i)][indexJ(i)][indexK(i)].y + alpha * p[i].y;
        tempX[i].z = jello->p[indexI(i)][indexJ(i)][indexK(i)].z + alpha * p[i].z;
      }
      double E_try = temp_IP_val(tempX, x_tilde, m, edges, l2, k, h, y_ground, contact_area);
      if(E_try <= E_last){
        // accept
        break;
      }
      alpha *= 0.5;
      if(alpha < 1e-12) {
        break;
      }
    }

    // update x
    for(int i=0; i < 512; i++){
      jello->p[indexI(i)][indexJ(i)][indexK(i)].x += alpha * p[i].x;
      jello->p[indexI(i)][indexJ(i)][indexK(i)].y += alpha * p[i].y;
      jello->p[indexI(i)][indexJ(i)][indexK(i)].z += alpha * p[i].z;
    }
    E_last = IP_val(jello, x_tilde, m, edges, l2, k, h, y_ground, contact_area);
    p = search_dir(jello, x_tilde, m, edges, l2, k, h, y_ground, contact_area);
    iter++;
  }

  // v = (x - x_n) / h
  for(int i=0; i < 512; i++) {
    jello->v[indexI(i)][indexJ(i)][indexK(i)].x = (jello->p[indexI(i)][indexJ(i)][indexK(i)].x - x_old[i].x) / h;
    jello->v[indexI(i)][indexJ(i)][indexK(i)].y = (jello->p[indexI(i)][indexJ(i)][indexK(i)].y - x_old[i].y) / h;
    jello->v[indexI(i)][indexJ(i)][indexK(i)].z = (jello->p[indexI(i)][indexJ(i)][indexK(i)].z - x_old[i].z) / h;
  }

  // write results back to jello->p and jello->v
  // copyVectorsToJello(jello, x, v);
}

void simulate(world *jello) {
  auto edges = generateJelloEdges();
  double tolerance = 1e-2;
  double timeStep = jello->dt;
  double m = jello->mass;
  double k = jello->kElastic;

  vector<point> xFlat;
  copyJelloPositions(jello, xFlat);

  vector<double> restLenSquared(edges.size());
  for (int i = 0; i < edges.size(); i++) {
    int idxA = edges[i][0];
    int idxB = edges[i][1];
    point diff;
    diff.x = xFlat[idxA].x - xFlat[idxB].x;
    diff.y = xFlat[idxA].y - xFlat[idxB].y;
    diff.z = xFlat[idxA].z - xFlat[idxB].z;

    restLenSquared[i] = pDOT(diff, diff);
  }

  // barrier
  double y_ground = -2;
  double contact_area = 1.0 / 7.0;

  stepForwardImplicitEuler(jello, edges, m, restLenSquared, k, timeStep, tolerance, xFlat, y_ground, contact_area);
}
