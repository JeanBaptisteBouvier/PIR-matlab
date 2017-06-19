#ifndef CR3BP_DERIVATIVES_H_INCLUDED
#define CR3BP_DERIVATIVES_H_INCLUDED
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//GSL
#include "../lib/gsl/gsl_matrix.h"
#include "../lib/gsl/gsl_vector.h"
#include "../lib/gsl/gsl_blas.h"


int cr3bp_derivatives_rel_48 (double t, const double y[], double f[], void *params);
int cr3bp_derivatives_rel_12 (double t, const double y[], double f[], void *params);

int bcp_derivatives_rel_12 (double t, const double y[], double f[], void *params);
int bcp_derivatives_rel_48 (double t, const double y[], double f[], void *params);

#endif // CR3BP_DERIVATIVES_H_INCLUDED
