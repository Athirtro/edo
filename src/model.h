#ifndef EDO_MODEL_H
#define EDO_MODEL_H

/*
 * Modelo agua-trucha-gusano (bloque 1: formulación).
 * Unidades:
 *   W       km^3 de agua liquida disponible
 *   T, G    toneladas de biomasa
 *   tiempo  anios
 */
typedef struct {
    double r;        /* tasa maxima de crecimiento de truchas, 1/anio */
    double K;        /* constante de media saturacion en W, km^3 */
    double mu;       /* tasa de maduracion trucha -> gusano, 1/anio */
    double delta;    /* tasa de mortalidad del gusano, 1/anio */
    double alpha;    /* agua retenida por tonelada de trucha, km^3/t */
    double beta;     /* agua retenida por tonelada de gusano, km^3/t */
    double eta;      /* rendimiento gusano/trucha, t/t */
    double xi;       /* rendimiento trucha/gusano, t/t */
    double W0;       /* agua liquida inicial, km^3 */
    double m_trucha; /* masa de una trucha, t */
    double f_desert; /* fraccion umbral de W0, adimensional */
} Params;

typedef struct {
    double W;
    double T;
    double G;
} State;

/* Conversiones de volumen de agua a masa y viceversa. */
double ton_to_km3(double tons);
double km3_to_ton(double km3);

/* Verifica positividad y las dos condiciones de invariancia. */
int params_valid(const Params *p);

/* Respuesta saturante de Michaelis-Menten: r*W/(K+W). */
double growth(double W, const Params *p);

/* Estado inicial: W = W0, T = m_trucha (una trucha), G = 0. */
void initial_state(State *s, const Params *p);

/* Cantidad conservada: W + alpha*T + beta*G. */
double total_water(const Params *p, const State *s);

/* Sistema completo de tres EDO. */
void rhs(const Params *p, const State *s, double *dW, double *dT, double *dG);

/* Sistema reducido a (T, G) usando la conservacion del agua. */
void rhs_reduced(const Params *p, double T, double G, double *dT, double *dG);

/* Reconstruye W a partir de (T, G) y el agua total inicial. */
void state_from_reduced(const Params *p, double T, double G, State *s);

#endif /* EDO_MODEL_H */
