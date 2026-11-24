#include "model.h"

#include <math.h>

/* 1 km^3 de agua = 10^9 toneladas. */
static const double TON_PER_KM3 = 1.0e9;

double ton_to_km3(double tons)
{
    return tons / TON_PER_KM3;
}

double km3_to_ton(double km3)
{
    return km3 * TON_PER_KM3;
}

int params_valid(const Params *p)
{
    if (p->r <= 0.0 || p->K <= 0.0 || p->mu <= 0.0 || p->delta <= 0.0 ||
        p->alpha <= 0.0 || p->beta <= 0.0 || p->eta <= 0.0 ||
        p->xi <= 0.0 || p->W0 <= 0.0 || p->m_trucha <= 0.0 ||
        p->f_desert <= 0.0 || p->f_desert >= 1.0) {
        return 0;
    }
    /* Evita que el campo empuje W por debajo de cero en las transiciones. */
    if (p->alpha < p->beta * p->eta) {
        return 0;
    }
    if (p->beta < p->alpha * p->xi) {
        return 0;
    }
    return 1;
}

double growth(double W, const Params *p)
{
    if (W <= 0.0) {
        return 0.0;
    }
    return p->r * W / (p->K + W);
}

void initial_state(State *s, const Params *p)
{
    s->W = p->W0;
    s->T = p->m_trucha;
    s->G = 0.0;
}

double total_water(const Params *p, const State *s)
{
    return s->W + p->alpha * s->T + p->beta * s->G;
}

void rhs(const Params *p, const State *s, double *dW, double *dT, double *dG)
{
    double g = growth(s->W, p);
    double maturation = p->mu * s->T;
    double mortality = p->delta * s->G;

    /* El agua retenida en cada transicion se ajusta para conservar el total. */
    *dW = -p->alpha * g * s->T
          + (p->alpha - p->beta * p->eta) * maturation
          + (p->beta - p->alpha * p->xi) * mortality;

    *dT = g * s->T - maturation + p->xi * mortality;
    *dG = p->eta * maturation - mortality;
}

void state_from_reduced(const Params *p, double T, double G, State *s)
{
    double W_tot = p->W0 + p->alpha * p->m_trucha; /* G(0) = 0 */
    s->W = W_tot - p->alpha * T - p->beta * G;
    s->T = T;
    s->G = G;
}

void rhs_reduced(const Params *p, double T, double G, double *dT, double *dG)
{
    State s;
    double dW;

    state_from_reduced(p, T, G, &s);
    rhs(p, &s, &dW, dT, dG);
    (void)dW;
}
