/*
 * Referencia independiente para t_colap: RK4 de paso fijo sobre (W,T,G).
 * Ejecutar con: make reference
 * No comparte el integrador adaptativo de src/collapse.c.
 */
#include <stdio.h>

#include "collapse.h"
#include "model.h"
#include "parameters.h"

static void rk4_step(const Params *p, const State *y, double h, State *out)
{
    State tmp;
    double a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3;

    rhs(p, y, &a1, &a2, &a3);

    tmp.W = y->W + 0.5 * h * a1;
    tmp.T = y->T + 0.5 * h * a2;
    tmp.G = y->G + 0.5 * h * a3;
    rhs(p, &tmp, &b1, &b2, &b3);

    tmp.W = y->W + 0.5 * h * b1;
    tmp.T = y->T + 0.5 * h * b2;
    tmp.G = y->G + 0.5 * h * b3;
    rhs(p, &tmp, &c1, &c2, &c3);

    tmp.W = y->W + h * c1;
    tmp.T = y->T + h * c2;
    tmp.G = y->G + h * c3;
    rhs(p, &tmp, &d1, &d2, &d3);

    out->W = y->W + (h / 6.0) * (a1 + 2.0 * b1 + 2.0 * c1 + d1);
    out->T = y->T + (h / 6.0) * (a2 + 2.0 * b2 + 2.0 * c2 + d2);
    out->G = y->G + (h / 6.0) * (a3 + 2.0 * b3 + 2.0 * c3 + d3);
}

int main(void)
{
    Params p;
    State y, next, event;
    double h = 0.001, t = 0.0, lo, hi, mid;
    double W_event;
    int k;

    params_default(&p);
    initial_state(&y, &p);
    W_event = collapse_threshold_water(&p);

    while (t < 20000.0) {
        rk4_step(&p, &y, h, &next);
        if (y.W > W_event && next.W <= W_event) {
            lo = 0.0;
            hi = h;
            for (k = 0; k < 50; k++) {
                mid = 0.5 * (lo + hi);
                rk4_step(&p, &y, mid, &event);
                if (event.W > W_event) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            printf("RK4 paso fijo h=%.3g: t_colap = %.9f anios\n",
                   h, t + hi);
            return 0;
        }
        y = next;
        t += h;
    }

    printf("No se detecto colapso antes de 20000 anios\n");
    return 1;
}
