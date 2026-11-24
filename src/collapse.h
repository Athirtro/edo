#ifndef EDO_COLLAPSE_H
#define EDO_COLLAPSE_H

#include "model.h"

typedef struct {
    int collapsed;
    double t_colap;
    double W_at_event;
    double max_conservation_error;
    double final_W;
    double final_T;
    double final_G;
    long steps;
} SimResult;

/* Nivel de agua que define el evento de desierto. */
double collapse_threshold_water(const Params *p);

/* Integracion RK4 adaptativa sobre S = W_tot - W; el evento se localiza
 * por biseccion de subpasos RK4. */
SimResult simulate_until_collapse(const Params *p, double t_max, double tol);

/* Estimacion por cuadratura bajo la aproximacion cuasi-estacionaria de G. */
double quadrature_t_colap(const Params *p);

#endif /* EDO_COLLAPSE_H */
