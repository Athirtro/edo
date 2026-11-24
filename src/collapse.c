#include "collapse.h"

#include <math.h>

#include "analysis.h"

#define H_MIN 1.0e-6
#define H_MAX 50.0

/* Variable de integracion: S = W_tot - W = alpha*T + beta*G.
 * Al inicio S = alpha*m_trucha ~ 5e-11 km^3, mientras que W ~ 1.36e9:
 * integrar W directamente haria que la tolerancia ignorara la cantidad
 * pequena que esta creciendo. */
typedef struct {
    double S;
    double T;
    double G;
} SState;

double collapse_threshold_water(const Params *p)
{
    return conventional_threshold_water(p);
}

static void rhs_s(const Params *p, const SState *y, double W_tot,
                  double *dS, double *dT, double *dG)
{
    State w;
    double dW;

    w.W = W_tot - y->S;
    w.T = y->T;
    w.G = y->G;
    rhs(p, &w, &dW, dT, dG);
    *dS = -dW;
}

static void rk4_step(const Params *p, const SState *y, double W_tot,
                     double h, SState *out)
{
    SState tmp;
    double a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3;

    rhs_s(p, y, W_tot, &a1, &a2, &a3);
    tmp.S = y->S + 0.5 * h * a1;
    tmp.T = y->T + 0.5 * h * a2;
    tmp.G = y->G + 0.5 * h * a3;
    rhs_s(p, &tmp, W_tot, &b1, &b2, &b3);

    tmp.S = y->S + 0.5 * h * b1;
    tmp.T = y->T + 0.5 * h * b2;
    tmp.G = y->G + 0.5 * h * b3;
    rhs_s(p, &tmp, W_tot, &c1, &c2, &c3);

    tmp.S = y->S + h * c1;
    tmp.T = y->T + h * c2;
    tmp.G = y->G + h * c3;
    rhs_s(p, &tmp, W_tot, &d1, &d2, &d3);

    out->S = y->S + (h / 6.0) * (a1 + 2.0 * b1 + 2.0 * c1 + d1);
    out->T = y->T + (h / 6.0) * (a2 + 2.0 * b2 + 2.0 * c2 + d2);
    out->G = y->G + (h / 6.0) * (a3 + 2.0 * b3 + 2.0 * c3 + d3);
}

/* Operador aceptado por el integrador: dos medios pasos corregidos por
 * Richardson. Se usa tanto para avanzar como para localizar el evento. */
static void rk4_step_extrapolated(const Params *p, const SState *y,
                                  double W_tot, double h, SState *out)
{
    SState full, half, both;

    rk4_step(p, y, W_tot, h, &full);
    rk4_step(p, y, W_tot, 0.5 * h, &half);
    rk4_step(p, &half, W_tot, 0.5 * h, &both);

    out->S = both.S + (both.S - full.S) / 15.0;
    out->T = both.T + (both.T - full.T) / 15.0;
    out->G = both.G + (both.G - full.G) / 15.0;
}

static double scaled_error(double a, double b, double tol)
{
    double scale = 1.0e-30 + tol * fmax(fabs(a), fabs(b));
    return fabs(b - a) / scale;
}

SimResult simulate_until_collapse(const Params *p, double t_max, double tol)
{
    SimResult r;
    SState y, y_half, y_full, y_both, y_event;
    State w;
    double t = 0.0, h = 1.0;
    double W_event = collapse_threshold_water(p);
    double W_tot = p->W0 + p->alpha * p->m_trucha;
    double S_event = W_tot - W_event;
    double err, ce, factor;

    r.collapsed = 0;
    r.t_colap = INFINITY;
    r.W_at_event = 0.0;
    r.max_conservation_error = 0.0;
    r.steps = 0;

    y.S = p->alpha * p->m_trucha;
    y.T = p->m_trucha;
    y.G = 0.0;
    r.final_W = W_tot - y.S;
    r.final_T = y.T;
    r.final_G = y.G;

    if (y.S >= S_event) {
        r.collapsed = 1;
        r.t_colap = 0.0;
        r.W_at_event = W_tot - y.S;
        return r;
    }

    while (t < t_max) {
        double t0 = t;
        if (h > t_max - t) {
            h = t_max - t;
        }

        rk4_step(p, &y, W_tot, h, &y_full);
        rk4_step(p, &y, W_tot, 0.5 * h, &y_half);
        rk4_step(p, &y_half, W_tot, 0.5 * h, &y_both);

        /* Diferencia entre un paso completo y dos medios pasos; el factor
         * 15 es el de Richardson para RK4. */
        err = fmax(scaled_error(y_full.S, y_both.S, tol),
                   fmax(scaled_error(y_full.T, y_both.T, tol),
                        scaled_error(y_full.G, y_both.G, tol))) / 15.0;

        if (err <= 1.0 || h <= H_MIN) {
            t += h;
            /* RK4 tiene error O(h^5): y_exacta ~= y_both + (y_both-y_full)/15. */
            y_both.S += (y_both.S - y_full.S) / 15.0;
            y_both.T += (y_both.T - y_full.T) / 15.0;
            y_both.G += (y_both.G - y_full.G) / 15.0;

            w.W = W_tot - y_both.S;
            w.T = y_both.T;
            w.G = y_both.G;
            ce = fabs(total_water(p, &w) - W_tot) / W_tot;
            if (ce > r.max_conservation_error) {
                r.max_conservation_error = ce;
            }

            if (y.S < S_event && y_both.S >= S_event) {
                /* Biseccion con el mismo operador que se acepta: asi el
                 * cruce no depende de donde caiga dentro del paso. */
                double lo = 0.0, hi = h, mid;
                int k;
                for (k = 0; k < 50; k++) {
                    mid = 0.5 * (lo + hi);
                    rk4_step_extrapolated(p, &y, W_tot, mid, &y_event);
                    if (y_event.S < S_event) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                mid = hi;
                rk4_step_extrapolated(p, &y, W_tot, mid, &y_event);
                w.W = W_tot - y_event.S;
                w.T = y_event.T;
                w.G = y_event.G;
                ce = fabs(total_water(p, &w) - W_tot) / W_tot;
                if (ce > r.max_conservation_error) {
                    r.max_conservation_error = ce;
                }
                r.collapsed = 1;
                r.t_colap = t0 + mid;
                r.W_at_event = W_tot - y_event.S;
                y = y_event;
            }

            if (!r.collapsed) {
                y = y_both;
            }
            r.steps++;

            factor = err > 0.0 ? 0.9 * pow(err, -0.2) : 5.0;
            if (factor < 0.2) {
                factor = 0.2;
            }
            if (factor > 5.0) {
                factor = 5.0;
            }
            h = fmin(h * factor, H_MAX);

            if (r.collapsed) {
                break;
            }
        } else {
            h = fmax(h * 0.5, H_MIN);
        }
    }

    r.final_W = W_tot - y.S;
    r.final_T = y.T;
    r.final_G = y.G;
    return r;
}

/* Con G ~= eta*mu*T/delta, el sistema se reduce a una EDO separable en W:
 *   dW/dt = -A(W) (W_tot - W),
 * con A(W) = [alpha*g(W) - mu(alpha-beta*eta) - eta*mu(beta-alpha*xi)] / D
 * y D = alpha + beta*eta*mu/delta. Integrando en u = ln(W_tot - W):
 *   t = int du / A(W_tot - exp(u)).
 */
double quadrature_t_colap(const Params *p)
{
    double W_event = collapse_threshold_water(p);
    double W_tot = p->W0 + p->alpha * p->m_trucha;
    double D = p->alpha + p->beta * p->eta * p->mu / p->delta;
    double released = p->mu * (p->alpha - p->beta * p->eta) +
                      p->eta * p->mu * (p->beta - p->alpha * p->xi);
    /* W_tot - W(0) = alpha*m_trucha se pierde al restar en coma flotante;
     * se calcula directamente. */
    double S0 = p->alpha * p->m_trucha;
    double S1 = W_tot - W_event;
    double u0, u1, du, sum = 0.0;
    int N = 2000, i;

    if (W_event >= p->W0 || S0 <= 0.0 || S1 <= 0.0) {
        return 0.0;
    }

    u0 = log(S0);
    u1 = log(S1);
    du = (u1 - u0) / N;

    for (i = 0; i <= N; i++) {
        double u = u0 + i * du;
        double W = W_tot - exp(u);
        double A = (p->alpha * growth(W, p) - released) / D;
        double f = 1.0 / A;
        double weight = (i == 0 || i == N) ? 1.0 :
                        (i % 2 == 0 ? 2.0 : 4.0);
        if (!isfinite(f) || f <= 0.0) {
            return INFINITY;
        }
        sum += weight * f;
    }

    return (du / 3.0) * sum;
}
