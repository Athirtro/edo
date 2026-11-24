#include "analysis.h"

int find_equilibria(const Params *p, Equilibrium *desert,
                    Equilibrium *persistence)
{
    double W_tot = p->W0 + p->alpha * p->m_trucha;
    double c = p->mu * (1.0 - p->xi * p->eta);
    double W_star, D;

    desert->W = W_tot;
    desert->T = 0.0;
    desert->G = 0.0;
    persistence->W = 0.0;
    persistence->T = 0.0;
    persistence->G = 0.0;

    if (c < 0.0) {
        return 0;
    }
    if (c == 0.0) {
        W_star = 0.0;
    } else {
        if (p->r <= c) {
            return 0;
        }
        W_star = p->K * c / (p->r - c);
    }
    if (W_star >= W_tot) {
        return 0;
    }

    D = p->alpha + p->beta * p->eta * p->mu / p->delta;
    persistence->W = W_star;
    persistence->T = (W_tot - W_star) / D;
    persistence->G = p->eta * p->mu * persistence->T / p->delta;
    return 1;
}

double bifurcation_r(const Params *p)
{
    double W_tot = p->W0 + p->alpha * p->m_trucha;
    return p->mu * (1.0 - p->xi * p->eta) * (1.0 + p->K / W_tot);
}

void jacobian_reduced(const Params *p, double T, double G, double J[2][2])
{
    double W = p->W0 + p->alpha * p->m_trucha - p->alpha * T - p->beta * G;
    double gp = p->r * p->K / ((p->K + W) * (p->K + W));

    J[0][0] = p->r * W / (p->K + W) - p->mu - p->alpha * T * gp;
    J[0][1] = -p->beta * T * gp + p->xi * p->delta;
    J[1][0] = p->eta * p->mu;
    J[1][1] = -p->delta;
}

Stability classify_equilibrium(double J[2][2])
{
    double tr = J[0][0] + J[1][1];
    double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    double disc = tr * tr - 4.0 * det;

    if (det < 0.0) {
        return SADDLE;
    }
    if (det == 0.0 || tr == 0.0) {
        return INDETERMINATE;
    }
    if (tr < 0.0) {
        return disc >= 0.0 ? STABLE_NODE : FOCUS_STABLE;
    }
    return disc >= 0.0 ? UNSTABLE_NODE : FOCUS_UNSTABLE;
}

double conventional_threshold_water(const Params *p)
{
    return p->f_desert * p->W0;
}

int persistence_below_threshold(const Params *p)
{
    Equilibrium desert, persist;

    if (!find_equilibria(p, &desert, &persist)) {
        return 0;
    }
    return persist.W < conventional_threshold_water(p) ? 1 : 0;
}
