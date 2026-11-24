/*
 * Pruebas del bloque 5: sensibilidad.
 * Ejecutar con: make test
 *
 * Fase TDD: este archivo se escribio antes que src/sensitivity.c.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "collapse.h"
#include "parameters.h"
#include "sensitivity.h"

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

static double find_range(const char *symbol, int is_min)
{
    size_t i;
    for (i = 0; i < PARAMETERS_COUNT; i++) {
        if (strcmp(PARAMETERS[i].symbol, symbol) == 0) {
            return is_min ? PARAMETERS[i].min : PARAMETERS[i].max;
        }
    }
    return NAN;
}

static double collapse_time(const Params *p)
{
    SimResult r = simulate_until_collapse(p, SENS_T_MAX, SENS_TOL);
    return r.collapsed ? r.t_colap : SENS_T_MAX;
}

static void test_free_parameter_table(void)
{
    Params p, copy;
    int i;

    params_default(&p);
    CHECK(SENS_FREE_COUNT == 8,
          "se analizan 8 de los 9 supuestos; m fija la escala de biomasa");

    for (i = 0; i < SENS_FREE_COUNT; i++) {
        double lo, hi, value;
        const char *symbol = free_param_symbol((FreeParam)i);

        CHECK(symbol != NULL && free_param_name((FreeParam)i) != NULL,
              "cada parametro libre tiene simbolo y nombre");
        free_param_range((FreeParam)i, &lo, &hi);
        CHECK(close_rel(lo, find_range(symbol, 1), 1e-12) &&
              close_rel(hi, find_range(symbol, 0), 1e-12),
              "los rangos vienen de la tabla de parametros");
        value = free_param_get(&p, (FreeParam)i);
        CHECK(value > lo && value < hi, "el valor base esta dentro del rango");

        copy = p;
        free_param_set(&copy, (FreeParam)i, value);
        CHECK(close_rel(free_param_get(&copy, (FreeParam)i), value, 1e-15),
              "set/get conserva el valor del parametro");
    }
}

static void test_local_derivative(void)
{
    Params p, plus, minus;
    double h, d, expected, t0, tp, tm, elasticity;
    double er, ek, emu;

    params_default(&p);
    h = 0.01 * p.r;
    d = local_derivative(&p, SENS_R, h);

    plus = p;
    minus = p;
    /* Diferencia externa con paso mitad del usado por la implementacion. */
    plus.r += 0.5 * h;
    minus.r -= 0.5 * h;
    t0 = collapse_time(&p);
    tp = collapse_time(&plus);
    tm = collapse_time(&minus);
    expected = (tp - tm) / h;

    CHECK(isfinite(d), "la derivada local es finita");
    CHECK(close_rel(d, expected, 0.02),
          "la derivada de r coincide con una diferencia central externa");
    CHECK(d < 0.0, "aumentar r acorta el colapso");
    CHECK(t0 > 0.0, "el tiempo base es finito");

    elasticity = (p.r / t0) * d;
    CHECK(close_rel(elasticity, (p.r / t0) * expected, 0.02),
          "la elasticidad es (p/t) veces la derivada");

    er = fabs(elasticity);
    ek = fabs(local_elasticity(&p, SENS_K, 0.01 * p.K));
    emu = fabs(local_elasticity(&p, SENS_MU, 0.01 * p.mu));
    CHECK(er > ek && ek > emu,
          "localmente dominan r, despues K y despues mu");
}

static void test_morris_screening(void)
{
    Params p;
    MorrisEffect e[SENS_FREE_COUNT];
    FreeParam order[SENS_FREE_COUNT];
    int i;

    params_default(&p);
    morris_screen(&p, 100, 12345u, e);

    for (i = 0; i < SENS_FREE_COUNT; i++) {
        CHECK(isfinite(e[i].mu_star) && e[i].mu_star > 0.0,
              "todo parametro tiene efecto elementary finito y positivo");
        CHECK(isfinite(e[i].sigma) && e[i].sigma >= 0.0,
              "la desviacion del efecto es finita");
        CHECK(e[i].n > 0, "hay muestras validas para cada parametro");
    }

    morris_rank(e, order);
    CHECK(order[0] == SENS_R && order[1] == SENS_K && order[2] == SENS_MU,
          "el cribado global rankea r, K y mu como dominantes");
    CHECK(e[SENS_R].mu_star > e[SENS_BETA].mu_star,
          "el efecto de r supera claramente al de beta");
    CHECK(order[0] != order[7], "el ranking no devuelve un orden degenerado");
}

static void test_morris_determinism(void)
{
    Params p;
    MorrisEffect a[SENS_FREE_COUNT], b[SENS_FREE_COUNT];
    int i;

    params_default(&p);
    morris_screen(&p, 12, 777u, a);
    morris_screen(&p, 12, 777u, b);
    for (i = 0; i < SENS_FREE_COUNT; i++) {
        CHECK(a[i].mu_star == b[i].mu_star &&
              a[i].sigma == b[i].sigma && a[i].n == b[i].n,
              "la misma semilla reproduce exactamente el cribado");
    }

    morris_screen(&p, 12, 778u, b);
    CHECK(a[SENS_R].mu_star != b[SENS_R].mu_star,
          "semillas distintas dan trayectorias distintas");
}

static void test_derivative_limits(void)
{
    Params p;
    double lo;

    params_default(&p);
    CHECK(isnan(local_derivative(&p, SENS_R, 0.0)),
          "paso cero da NaN");
    CHECK(isnan(local_derivative(&p, (FreeParam)99, 1e-6)),
          "parametro desconocido da NaN");
    CHECK(isnan(local_derivative(&p, SENS_ALPHA, 1e6)) &&
          isnan(local_derivative(&p, SENS_BETA, 1e6)),
          "pasos que violan las restricciones dan NaN");

    free_param_range(SENS_R, &lo, NULL);
    CHECK(isfinite(local_derivative(&p, SENS_R, 1e6)),
          "un paso grande se recorta a los extremos del rango");
}

static void test_conservation_nonidentifiability(void)
{
    Params p, scaled, broken;
    SimResult base, equivalent, changed;
    double lambdas[] = {0.25, 0.5, 2.0, 5.0};
    int i;

    params_default(&p);
    base = simulate_until_collapse(&p, SENS_T_MAX, SENS_TOL);

    for (i = 0; i < 4; i++) {
        conservation_scaling(&p, lambdas[i], &scaled);
        CHECK(params_valid(&scaled),
              "la transformacion de escala conserva las restricciones");
        CHECK(close_rel(scaled.alpha * scaled.m_trucha,
                        p.alpha * p.m_trucha, 1e-12) &&
              close_rel(scaled.W0 + scaled.alpha * scaled.m_trucha,
                        p.W0 + p.alpha * p.m_trucha, 1e-12),
              "el agua total inicial no cambia con la escala de biomasa");

        equivalent = simulate_until_collapse(&scaled, SENS_T_MAX, SENS_TOL);
        CHECK(close_rel(equivalent.t_colap, base.t_colap, 1e-6) &&
              close_rel(equivalent.final_W, base.final_W, 1e-9),
              "alpha, beta y m escalados dan el mismo evento");
    }

    broken = p;
    broken.alpha *= 2.0;
    broken.beta *= 2.0;
    changed = simulate_until_collapse(&broken, SENS_T_MAX, SENS_TOL);
    CHECK(fabs(changed.t_colap - base.t_colap) > 0.01 * base.t_colap,
          "escalar alpha y beta sin reescalar la biomasa cambia t_colap >1%");
}

int main(void)
{
    test_free_parameter_table();
    test_local_derivative();
    test_morris_screening();
    test_morris_determinism();
    test_derivative_limits();
    test_conservation_nonidentifiability();

    if (failures == 0) {
        printf("OK: todas las pruebas del bloque 5 pasaron\n");
        return 0;
    }
    printf("%d prueba(s) fallaron\n", failures);
    return 1;
}
