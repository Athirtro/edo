/*
 * Pruebas del bloque 3: analisis cualitativo.
 * Ejecutar con: make test
 *
 * Fase TDD: este archivo se escribio antes que src/analysis.c.
 */
#include <math.h>
#include <stdio.h>

#include "analysis.h"
#include "parameters.h"

static int failures = 0;

#define CHECK(cond, msg)                                                        \
    do {                                                                        \
        if (!(cond)) {                                                          \
            fprintf(stderr, "FAIL  %s  [%s:%d]\n", (msg), __FILE__, __LINE__);  \
            failures++;                                                         \
        }                                                                       \
    } while (0)

static int close_to(double a, double b, double tol)
{
    double scale = fmax(1.0, fmax(fabs(a), fabs(b)));
    return fabs(a - b) <= tol * scale;
}

static int close_rel(double a, double b, double tol)
{
    double scale = fmax(fabs(a), fabs(b));
    if (scale < 1e-30) {
        scale = 1.0;
    }
    return fabs(a - b) <= tol * scale;
}

static double residual_T(const Params *p, double T, double G)
{
    double dT, dG;
    rhs_reduced(p, T, G, &dT, &dG);
    (void)dG;
    return dT;
}

static double residual_G(const Params *p, double T, double G)
{
    double dT, dG;
    rhs_reduced(p, T, G, &dT, &dG);
    (void)dT;
    return dG;
}

static void test_desert_equilibrium(void)
{
    Params p;
    Equilibrium desert, persist;

    params_default(&p);
    (void)find_equilibria(&p, &desert, &persist);

    CHECK(desert.T == 0.0 && desert.G == 0.0, "desierto: T=G=0");
    CHECK(close_to(desert.W, p.W0 + p.alpha * p.m_trucha, 1e-12),
          "desierto: W = W_tot");
    CHECK(fabs(residual_T(&p, desert.T, desert.G)) < 1e-15,
          "desierto: dT/dt = 0");
    CHECK(fabs(residual_G(&p, desert.T, desert.G)) < 1e-15,
          "desierto: dG/dt = 0");
}

static void test_persistence_equilibrium(void)
{
    Params p;
    Equilibrium desert, persist;
    int exists;
    double c;

    params_default(&p);
    exists = find_equilibria(&p, &desert, &persist);
    CHECK(exists == 1, "con los parametros base existe persistencia");

    c = p.mu * (1.0 - p.xi * p.eta);
    CHECK(close_rel(growth(persist.W, &p), c, 1e-10),
          "W de persistencia satisface g(W)=mu(1-xi*eta)");
    CHECK(persist.T > 0.0 && persist.G > 0.0,
          "la rama de persistencia tiene biomasa positiva");
    CHECK(close_to(persist.G, p.eta * p.mu * persist.T / p.delta, 1e-10),
          "G = eta*mu*T/delta en equilibrio");
    {
        double scale = fabs(growth(persist.W, &p) * persist.T) +
                       p.mu * persist.T + p.xi * p.delta * persist.G;
        CHECK(fabs(residual_T(&p, persist.T, persist.G)) / scale < 1e-12,
              "persistencia: dT/dt = 0 (escala relativa)");
    }
    CHECK(fabs(residual_G(&p, persist.T, persist.G)) < 1e-12,
          "persistencia: dG/dt = 0");
    CHECK(close_to(persist.W + p.alpha * persist.T + p.beta * persist.G,
                   p.W0 + p.alpha * p.m_trucha, 1e-10),
          "persistencia conserva el agua total");
}

static void test_bifurcation_threshold(void)
{
    Params p;
    Equilibrium desert, persist;
    double rcrit, J[2][2], det, c;

    params_default(&p);
    rcrit = bifurcation_r(&p);

    /* En r = r_crit, g(W_tot) = c y el desierto deja de ser hiperbolico. */
    p.r = rcrit;
    c = p.mu * (1.0 - p.xi * p.eta);
    CHECK(close_rel(growth(p.W0 + p.alpha * p.m_trucha, &p), c, 1e-10),
          "en r_crit, g(W_tot) = c");
    jacobian_reduced(&p, 0.0, 0.0, J);
    det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    CHECK(fabs(det) < 1e-12, "en r_crit, det J del desierto es cero");

    p.r = 0.99 * rcrit;
    CHECK(find_equilibria(&p, &desert, &persist) == 0,
          "por debajo de r_crit no hay persistencia");

    p.r = 1.01 * rcrit;
    CHECK(find_equilibria(&p, &desert, &persist) == 1 && persist.T > 0.0,
          "por encima de r_crit aparece la rama de persistencia");
}

static void test_jacobian_finite_differences(void)
{
    Params p;
    Equilibrium desert, persist;
    double J[2][2], Jfd[2][2];
    double T, G;
    double hT, hG;
    double dT, dG;

    params_default(&p);
    (void)find_equilibria(&p, &desert, &persist);
    T = persist.T;
    G = persist.G;
    hT = 1.0e13;
    hG = 1.0e13;
    jacobian_reduced(&p, T, G, J);

    rhs_reduced(&p, T + hT, G, &dT, &dG);
    Jfd[0][0] = dT;
    rhs_reduced(&p, T - hT, G, &dT, &dG);
    Jfd[0][0] = (Jfd[0][0] - dT) / (2.0 * hT);

    rhs_reduced(&p, T, G + hG, &dT, &dG);
    Jfd[0][1] = dT;
    rhs_reduced(&p, T, G - hG, &dT, &dG);
    Jfd[0][1] = (Jfd[0][1] - dT) / (2.0 * hG);

    rhs_reduced(&p, T + hT, G, &dT, &dG);
    Jfd[1][0] = dG;
    rhs_reduced(&p, T - hT, G, &dT, &dG);
    Jfd[1][0] = (Jfd[1][0] - dG) / (2.0 * hT);

    rhs_reduced(&p, T, G + hG, &dT, &dG);
    Jfd[1][1] = dG;
    rhs_reduced(&p, T, G - hG, &dT, &dG);
    Jfd[1][1] = (Jfd[1][1] - dG) / (2.0 * hG);

    CHECK(close_rel(J[0][0], Jfd[0][0], 1e-4), "J[0][0] analitico vs FD");
    CHECK(close_rel(J[0][1], Jfd[0][1], 1e-4), "J[0][1] analitico vs FD");
    CHECK(close_rel(J[1][0], Jfd[1][0], 1e-4), "J[1][0] analitico vs FD");
    CHECK(close_rel(J[1][1], Jfd[1][1], 1e-4), "J[1][1] analitico vs FD");
}

static void test_stability_classes(void)
{
    Params p;
    Equilibrium desert, persist;
    double J[2][2], tr, det;

    params_default(&p);
    (void)find_equilibria(&p, &desert, &persist);

    jacobian_reduced(&p, desert.T, desert.G, J);
    tr = J[0][0] + J[1][1];
    det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    CHECK(det < 0.0, "el desierto es silla con los parametros base");
    CHECK(classify_equilibrium(J) == SADDLE, "clasificacion: silla");

    jacobian_reduced(&p, persist.T, persist.G, J);
    tr = J[0][0] + J[1][1];
    det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    CHECK(tr < 0.0 && det > 0.0, "la persistencia es estable");
    CHECK(classify_equilibrium(J) == STABLE_NODE ||
          classify_equilibrium(J) == FOCUS_STABLE,
          "clasificacion: estable");
}

static void test_threshold_comparison(void)
{
    Params p;
    double w_thr;

    params_default(&p);
    w_thr = conventional_threshold_water(&p);
    CHECK(close_to(w_thr, p.f_desert * p.W0, 1e-12),
          "umbral convencional = f_desert * W0");
    CHECK(persistence_below_threshold(&p) == 1,
          "W* de persistencia queda por debajo del umbral convencional");
}

static void test_center_classification(void)
{
    double J[2][2] = {{0.0, 1.0}, {-1.0, 0.0}};

    CHECK(classify_equilibrium(J) == INDETERMINATE,
          "traza cero no se clasifica como foco inestable");
}

int main(void)
{
    test_desert_equilibrium();
    test_persistence_equilibrium();
    test_bifurcation_threshold();
    test_jacobian_finite_differences();
    test_stability_classes();
    test_center_classification();
    test_threshold_comparison();

    if (failures == 0) {
        printf("OK: todas las pruebas del bloque 3 pasaron\n");
        return 0;
    }
    printf("%d prueba(s) fallaron\n", failures);
    return 1;
}
