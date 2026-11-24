#ifndef EDO_VISUALIZATION_H
#define EDO_VISUALIZATION_H

#include "analysis.h"
#include "model.h"

/* Numeros maximos de muestras por figura. */
#define SERIES_MAX 8192
#define PHASE_MAX 8192
#define HEAT_R 20
#define HEAT_K 20

typedef struct {
    double t[SERIES_MAX];
    double W[SERIES_MAX];
    double T[SERIES_MAX];
    double G[SERIES_MAX];
    int n;
} Series;

typedef struct {
    double T[PHASE_MAX];
    double G[PHASE_MAX];
    int n;
} Curve;

typedef struct {
    Curve trajectory;
    Curve separatrix;
    Equilibrium desert;
    Equilibrium persist;
} PhaseData;

typedef struct {
    int nr;
    int nk;
    double r[HEAT_R];
    double K[HEAT_K];
    double t[HEAT_R][HEAT_K];
} Heatmap;

/* Generacion de datos. Devuelven 0 en caso de exito. */
int series_build(const Params *p, Series *s);
int phase_build(const Params *p, PhaseData *d);
int heatmap_build(const Params *p, Heatmap *h);

/* t_colap de una celda; devuelve 0.0 si el indice esta fuera de rejilla. */
double heatmap_at(const Heatmap *h, int i, int j);

/* Escritura de SVG (SVG 1.1 basico, sin dependencias). Devuelven 0 si OK. */
int write_series_svg(const Series *s, const char *path);
int write_phase_svg(const PhaseData *d, const char *path);
int write_heatmap_svg(const Heatmap *h, const char *path);

#endif /* EDO_VISUALIZATION_H */
