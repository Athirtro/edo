#ifndef EDO_ANALYSIS_H
#define EDO_ANALYSIS_H

#include "model.h"

typedef struct {
    double W;
    double T;
    double G;
} Equilibrium;

typedef enum {
    STABLE_NODE,
    UNSTABLE_NODE,
    FOCUS_STABLE,
    FOCUS_UNSTABLE,
    SADDLE,
    INDETERMINATE
} Stability;

/* Devuelve 1 si existe la rama de persistencia, 0 en caso contrario.
 * El argumento `persistence` sólo es significativo cuando devuelve 1;
 * en caso contrario se rellena con ceros. */
int find_equilibria(const Params *p, Equilibrium *desert,
                    Equilibrium *persistence);

/* Valor de r en el que la rama de persistencia colisiona con el desierto. */
double bifurcation_r(const Params *p);

/* Jacobiano del sistema reducido (T,G). */
void jacobian_reduced(const Params *p, double T, double G, double J[2][2]);

/* Clasificacion por traza y determinante. */
Stability classify_equilibrium(double J[2][2]);

/* Umbral convencional: f_desert * W0. */
double conventional_threshold_water(const Params *p);

/* 1 si W* de persistencia esta por debajo del umbral convencional. */
int persistence_below_threshold(const Params *p);

#endif /* EDO_ANALYSIS_H */
