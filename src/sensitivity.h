#ifndef EDO_SENSITIVITY_H
#define EDO_SENSITIVITY_H

#include <stdint.h>

#include "model.h"
#include "parameters.h"

/* Tiempo de simulacion usado en el analisis global; si no hay colapso se
 * censura t_colap a este valor. */
#define SENS_T_MAX 20000.0
/* Tolerancia del integrador en todos los calculos de sensibilidad. */
#define SENS_TOL 1e-6
/* Paso en coordenada normalizada del cribado Morris. */
#define SENS_MORRIS_DELTA 0.25

typedef enum {
    SENS_R,
    SENS_K,
    SENS_MU,
    SENS_DELTA,
    SENS_ALPHA,
    SENS_BETA,
    SENS_ETA,
    SENS_XI
} FreeParam;

#define SENS_FREE_COUNT 8

typedef struct {
    double mu_star;  /* media de |efecto elemental| */
    double sigma;    /* desviacion del efecto elemental */
    int n;           /* efectos validos acumulados */
} MorrisEffect;

const char *free_param_symbol(FreeParam which);
const char *free_param_name(FreeParam which);
void free_param_range(FreeParam which, double *min, double *max);
double free_param_get(const Params *p, FreeParam which);
void free_param_set(Params *p, FreeParam which, double value);

/* Derivada central de t_colap respecto al parametro; NaN si el paso sale
 * de los rangos o viola las restricciones de validez. */
double local_derivative(const Params *base, FreeParam which, double h);
double local_elasticity(const Params *base, FreeParam which, double h);

/* Cribado global: trayectorias Morris de dos niveles con semilla fija. */
void morris_screen(const Params *base, int trajectories, uint32_t seed,
                   MorrisEffect effects[SENS_FREE_COUNT]);
/* Orden de parametros de mayor a menor mu_star. */
void morris_rank(const MorrisEffect effects[SENS_FREE_COUNT],
                 FreeParam order[SENS_FREE_COUNT]);

/* Simetria exacta de la conservacion: alpha y beta escalados, m al inverso. */
void conservation_scaling(const Params *base, double lambda, Params *out);

#endif /* EDO_SENSITIVITY_H */
