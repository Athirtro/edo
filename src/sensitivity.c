#include "sensitivity.h"

#include <math.h>
#include <string.h>

#include "collapse.h"

static const char *SYMBOLS[SENS_FREE_COUNT] = {
    "r", "K", "mu", "delta", "alpha", "beta", "eta", "xi"
};

static const char *NAMES[SENS_FREE_COUNT] = {
    "tasa maxima de crecimiento",
    "agua de media saturacion",
    "tasa de maduracion",
    "tasa de mortalidad del gusano",
    "agua retenida por tonelada de trucha",
    "agua retenida por tonelada de gusano",
    "rendimiento gusano/trucha",
    "rendimiento trucha/gusano"
};

static const ParamMeta *meta_of(FreeParam which)
{
    size_t i;

    if ((int)which < 0 || (int)which >= SENS_FREE_COUNT) {
        return NULL;
    }
    for (i = 0; i < PARAMETERS_COUNT; i++) {
        if (strcmp(PARAMETERS[i].symbol, SYMBOLS[which]) == 0) {
            return &PARAMETERS[i];
        }
    }
    return NULL;
}

const char *free_param_symbol(FreeParam which)
{
    const ParamMeta *m = meta_of(which);
    return m != NULL ? m->symbol : NULL;
}

const char *free_param_name(FreeParam which)
{
    if ((int)which < 0 || (int)which >= SENS_FREE_COUNT) {
        return NULL;
    }
    return NAMES[which];
}

void free_param_range(FreeParam which, double *min, double *max)
{
    const ParamMeta *m = meta_of(which);

    if (m != NULL) {
        if (min != NULL) {
            *min = m->min;
        }
        if (max != NULL) {
            *max = m->max;
        }
    } else {
        if (min != NULL) {
            *min = NAN;
        }
        if (max != NULL) {
            *max = NAN;
        }
    }
}

double free_param_get(const Params *p, FreeParam which)
{
    switch (which) {
    case SENS_R:     return p->r;
    case SENS_K:     return p->K;
    case SENS_MU:    return p->mu;
    case SENS_DELTA: return p->delta;
    case SENS_ALPHA: return p->alpha;
    case SENS_BETA:  return p->beta;
    case SENS_ETA:   return p->eta;
    case SENS_XI:    return p->xi;
    }
    return NAN;
}

void free_param_set(Params *p, FreeParam which, double value)
{
    switch (which) {
    case SENS_R:     p->r = value; break;
    case SENS_K:     p->K = value; break;
    case SENS_MU:    p->mu = value; break;
    case SENS_DELTA: p->delta = value; break;
    case SENS_ALPHA: p->alpha = value; break;
    case SENS_BETA:  p->beta = value; break;
    case SENS_ETA:   p->eta = value; break;
    case SENS_XI:    p->xi = value; break;
    }
}

static double collapse_time(const Params *p)
{
    SimResult r = simulate_until_collapse(p, SENS_T_MAX, SENS_TOL);

    if (r.collapsed && isfinite(r.t_colap)) {
        return r.t_colap;
    }
    return SENS_T_MAX;
}

double local_derivative(const Params *base, FreeParam which, double h)
{
    Params plus, minus;
    const ParamMeta *m = meta_of(which);
    double v, lo, hi, tp, tm;

    if (m == NULL || !(h > 0.0)) {
        return NAN;
    }
    v = free_param_get(base, which);
    lo = fmax(m->min, v - h);
    hi = fmin(m->max, v + h);
    if (!(hi > lo)) {
        return NAN;
    }

    plus = *base;
    minus = *base;
    free_param_set(&plus, which, hi);
    free_param_set(&minus, which, lo);
    if (!params_valid(&plus) || !params_valid(&minus)) {
        return NAN;
    }

    tp = collapse_time(&plus);
    tm = collapse_time(&minus);
    return (tp - tm) / (hi - lo);
}

double local_elasticity(const Params *base, FreeParam which, double h)
{
    double d = local_derivative(base, which, h);
    double t = collapse_time(base);
    double v = free_param_get(base, which);

    if (!isfinite(d) || !(t > 0.0) || v == 0.0) {
        return NAN;
    }
    return (v / t) * d;
}

static uint32_t next_u32(uint32_t *state)
{
    *state = 1664525u * (*state) + 1013904223u;
    return *state;
}

static double unit_random(uint32_t *state)
{
    return (double)(next_u32(state) >> 8) / 16777216.0;
}

static void set_coordinate(Params *p, FreeParam which, double x)
{
    double lo, hi;

    free_param_range(which, &lo, &hi);
    if (x < 0.0) {
        x = 0.0;
    }
    if (x > 1.0) {
        x = 1.0;
    }
    free_param_set(p, which, lo + x * (hi - lo));
}

void morris_screen(const Params *base, int trajectories, uint32_t seed,
                   MorrisEffect effects[SENS_FREE_COUNT])
{
    double sum[SENS_FREE_COUNT] = {0.0};
    double sum_abs[SENS_FREE_COUNT] = {0.0};
    double sum2[SENS_FREE_COUNT] = {0.0};
    int n[SENS_FREE_COUNT] = {0};
    uint32_t rng = seed;
    int tr, attempt, i, j;

    for (tr = 0; tr < trajectories; tr++) {
        Params p;
        int order[SENS_FREE_COUNT];
        int sign[SENS_FREE_COUNT];
        double x[SENS_FREE_COUNT];
        int valid = 0;

        for (attempt = 0; attempt < 100 && !valid; attempt++) {
            for (i = 0; i < SENS_FREE_COUNT; i++) {
                double u = unit_random(&rng);
                order[i] = i;
                sign[i] = unit_random(&rng) > 0.5 ? 1 : -1;
                x[i] = sign[i] > 0 ? (1.0 - SENS_MORRIS_DELTA) * u :
                       SENS_MORRIS_DELTA +
                       (1.0 - SENS_MORRIS_DELTA) * u;
            }
            for (i = SENS_FREE_COUNT - 1; i > 0; i--) {
                int other = (int)(unit_random(&rng) * (double)(i + 1));
                int tmp = order[i];
                order[i] = order[other];
                order[other] = tmp;
            }

            p = *base;
            for (i = 0; i < SENS_FREE_COUNT; i++) {
                set_coordinate(&p, (FreeParam)order[i], x[order[i]]);
            }
            if (!params_valid(&p)) {
                continue;
            }
            valid = 1;

            for (i = 0; i < SENS_FREE_COUNT; i++) {
                FreeParam q = (FreeParam)order[i];
                Params trial;
                double x_old = x[q];
                double x_new = x_old + SENS_MORRIS_DELTA * (double)sign[q];
                double dx, t0, t1, ee;

                if (x_new < 0.0 || x_new > 1.0) {
                    continue;
                }
                trial = p;
                set_coordinate(&trial, q, x_new);
                dx = x_new - x_old;
                if (!params_valid(&trial)) {
                    double x_alt = x_old - SENS_MORRIS_DELTA * (double)sign[q];
                    if (x_alt < 0.0 || x_alt > 1.0) {
                        continue;
                    }
                    set_coordinate(&trial, q, x_alt);
                    dx = x_alt - x_old;
                    if (!params_valid(&trial)) {
                        continue;
                    }
                    x_new = x_alt;
                }

                t0 = collapse_time(&p);
                t1 = collapse_time(&trial);
                ee = (t1 - t0) / dx;
                sum[q] += ee;
                sum_abs[q] += fabs(ee);
                sum2[q] += ee * ee;
                n[q]++;

                p = trial;
                x[q] = x_new;
            }
        }
    }

    for (j = 0; j < SENS_FREE_COUNT; j++) {
        double mean, var;

        effects[j].n = n[j];
        if (n[j] == 0) {
            effects[j].mu_star = 0.0;
            effects[j].sigma = 0.0;
            continue;
        }
        effects[j].mu_star = sum_abs[j] / n[j];
        mean = sum[j] / n[j];
        var = sum2[j] / n[j] - mean * mean;
        effects[j].sigma = sqrt(fmax(0.0, var));
    }
}

void morris_rank(const MorrisEffect effects[SENS_FREE_COUNT],
                 FreeParam order[SENS_FREE_COUNT])
{
    int i, j;

    for (i = 0; i < SENS_FREE_COUNT; i++) {
        order[i] = (FreeParam)i;
    }
    for (i = 1; i < SENS_FREE_COUNT; i++) {
        FreeParam current = order[i];
        j = i - 1;
        while (j >= 0 &&
               (effects[current].mu_star > effects[order[j]].mu_star ||
                (effects[current].mu_star == effects[order[j]].mu_star &&
                 effects[current].sigma > effects[order[j]].sigma))) {
            order[j + 1] = order[j];
            j--;
        }
        order[j + 1] = current;
    }
}

void conservation_scaling(const Params *base, double lambda, Params *out)
{
    *out = *base;
    if (!(lambda > 0.0)) {
        return;
    }
    out->alpha = base->alpha * lambda;
    out->beta = base->beta * lambda;
    out->m_trucha = base->m_trucha / lambda;
}
