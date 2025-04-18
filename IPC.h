#ifndef _IPC_H_
#define _IPC_H_

#include <vector>
#include <array>
#include "jello.h"  // struct world, struct point, etc.

/**
 * A simple struct to represent IJV (Index,Index,Value) for a sparse matrix.
 * You use it when building the Hessian, etc.
 */
struct IJV
{
    std::vector<int>    I; // row indices
    std::vector<int>    J; // column indices
    std::vector<double> V; // values
};

void copyJelloPositions(world *jello, std::vector<point> &x);

std::vector< std::array<int,2> > generateJelloEdges();

void simulate(world *jello);

void stepForwardImplicitEuler(
    world *jello,
    std::vector< std::array<int,2> > &edges,
    double m, 
    std::vector<double> &l2,
    double k,
    double h,
    double tol,
    std::vector<point> &x_old,
    double y_ground,
    double contact_area);

double init_step_size(world *jello, double y_ground, std::vector<point> &p);

std::vector<double> solveSparseSystemEigen(IJV &A, std::vector<double> &rhs, int n);

#endif  // _IPC_H_