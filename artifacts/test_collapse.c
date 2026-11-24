/*
 * Pruebas del bloque 4: tiempo de colapso.
 * Ejecutar con: make test
 *
 * Fase TDD: este archivo se escribio antes que src/collapse.c.
 */
#include <math.h>
#include <stdio.h>

#include "analysis.h"
#include "collapse.h"
#include "parameters.h"

static int failures = 0;

#define CHECK(cond, msg)                                                        \
    do {                                                                        \
        if (!(cond)) {                                                          \
            fprintf(stderr, "FAIL  %s  [%s:%d]\n", (msg), __FILE__, __LINE__);  \
            failures++;                                                         \
        }                                                                       \
    } while (0)

static int close_rel(double a, double b, double tol)
{
    double scale = fmax(fabs(a), fabs(b));
    if (scale < 1e-30) {
        scale = 1.0;
    }
    return fabs(a - b) <= tol * scale;
}

static void test_event_definition(void)
{
    Params p;

    params_default(&p);
    CHECK(close_rel(collapse_threshold_water(&p), p.f_desert * p.W0, 1e-12),
          "el evento ocurre cuando W = f_desert * W0");
    CHECK(close_rel(collapse_threshold_water(&p), 1.815477e8, 1e-6),
          "el umbral base coincide con la cifra documentada");
}

static void test_quadrature(void)
{
    Params p;
    double tq;

    params_default(&p);
    tq = quadrature_t_colap(&p);
    CHECK(isfinite(tq) && tq > 0.0, "la cuadratura da un tiempo finito positivo");
    CHECK(close_rel(tq, 4366.03092047, 1e-7),
          "la cuadratura reproduce la referencia Simpson con N=2000");
}

static void test_numerical_collapse(void)
{
    Params p;
    SimResult r;
    double W_tot;

    params_default(&p);
    r = simulate_until_collapse(&p, 20000.0, 1e-6);
    W_tot = p.W0 + p.alpha * p.m_trucha;

    CHECK(r.collapsed == 1, "el conjunto base colapsa antes de 20000 anios");
    CHECK(fabs(r.t_colap - 4024.278215) < 0.002,
          "t_colap coincide con RK4 de paso fijo h=0.001");
    CHECK(close_rel(r.W_at_event, collapse_threshold_water(&p), 1e-6),
          "el estado del evento esta en el umbral");
    CHECK(close_rel(r.final_W, collapse_threshold_water(&p), 1e-6),
          "el estado final esta en el umbral, no en un paso posterior");
    CHECK(close_rel(r.final_W + p.alpha * r.final_T + p.beta * r.final_G,
                    W_tot, 1e-12),
          "el estado final conserva W_tot");
    CHECK(r.max_conservation_error < 1e-13,
          "el integrador conserva el agua total en toda la trayectoria");
    CHECK(r.steps > 0 && r.steps < 10000,
          "el integrador usa un numero razonable de pasos");
}

static void test_agreement_between_methods(void)
{
    Params p;
    SimResult r;
    double tq, rel;

    params_default(&p);
    r = simulate_until_collapse(&p, 20000.0, 1e-6);
    tq = quadrature_t_colap(&p);
    rel = fabs(tq - r.t_colap) / r.t_colap;

    CHECK(rel < 0.10,
          "la aproximacion cuasi-estacionaria coincide con la simulacion (<10%)");
}

static void test_tolerance_convergence(void)
{
    Params p;
    SimResult coarse, fine;
    double ref = 4024.278215;

    params_default(&p);
    coarse = simulate_until_collapse(&p, 20000.0, 1e-6);
    fine = simulate_until_collapse(&p, 20000.0, 1e-7);

    CHECK(coarse.collapsed == 1 && fine.collapsed == 1,
          "ambas tolerancias detectan el colapso");
    CHECK(fabs(fine.t_colap - coarse.t_colap) < 0.005,
          "reducir la tolerancia x10 cambia t_colap menos de 0.005 anios");
    CHECK(fabs(fine.t_colap - ref) < fabs(coarse.t_colap - ref),
          "la tolerancia mas finita se acerca a la referencia");
}

static void test_event_operator_consistency(void)
{
    Params p;
    SimResult r;
    double W_event;

    /* Umbral elegido para que el cruce caiga muy cerca del final de un paso
     * grande: detecta si la biseccion usa un operador distinto del aceptado. */
    params_default(&p);
    p.f_desert = 0.11344477252853814;
    W_event = collapse_threshold_water(&p);
    r = simulate_until_collapse(&p, 20000.0, 1e-6);

    CHECK(r.collapsed == 1, "el umbral alternativo tambien se alcanza");
    CHECK(close_rel(r.W_at_event, W_event, 1e-8),
          "el evento se localiza con el operador aceptado");
    CHECK(close_rel(r.final_W, W_event, 1e-8),
          "el estado devuelto esta en el umbral, no pasado el paso");
}

static void test_no_collapse_when_desert_is_stable(void)
{
    Params p;
    SimResult r;
    double W_tot;

    params_default(&p);
    p.r = 0.9 * bifurcation_r(&p);
    r = simulate_until_collapse(&p, 20000.0, 1e-6);
    W_tot = p.W0 + p.alpha * p.m_trucha;

    CHECK(r.collapsed == 0, "si el desierto es estable no hay colapso");
    CHECK(r.final_T < p.m_trucha, "la trucha inicial decae");
    CHECK(close_rel(r.final_W + p.alpha * r.final_T + p.beta * r.final_G,
                    W_tot, 1e-12),
          "el estado final conserva W_tot sin colapso");
    CHECK(r.max_conservation_error < 1e-13,
          "la conservacion tambien se cumple sin colapso");
}

int main(void)
{
    test_event_definition();
    test_quadrature();
    test_numerical_collapse();
    test_agreement_between_methods();
    test_tolerance_convergence();
    test_event_operator_consistency();
    test_no_collapse_when_desert_is_stable();

    if (failures == 0) {
        printf("OK: todas las pruebas del bloque 4 pasaron\n");
        return 0;
    }
    printf("%d prueba(s) fallaron\n", failures);
    return 1;
}
