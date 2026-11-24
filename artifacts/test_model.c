/*
 * Pruebas del bloque 1: formulación del modelo.
 * Ejecutar con: make test
 *
 * Fase TDD: este archivo se escribió antes que src/model.c.
 */
#include <math.h>
#include <stdio.h>

#include "model.h"

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

static Params sample_params(void)
{
    Params p;
    p.r = 0.02;             /* 1/año */
    p.K = 1.0e9;            /* km³ */
    p.mu = 1.0 / 4000.0;    /* 1/año */
    p.delta = 1.0 / 500.0;  /* 1/año */
    p.alpha = 1.0e-9;       /* km³ por tonelada de trucha */
    p.beta = 5.0e-10;       /* km³ por tonelada de gusano */
    p.eta = 1.5;            /* t gusano / t trucha */
    p.xi = 0.3;             /* t trucha / t gusano */
    p.W0 = 1.3e9;           /* km³ */
    p.m_trucha = 0.05;      /* toneladas (una trucha de 50 kg) */
    return p;
}

/* Tarea 1: unidades y factores de conversión */
static void test_units(void)
{
    CHECK(close_to(ton_to_km3(1.0e9), 1.0, 1e-12),
          "10^9 toneladas de agua equivalen a 1 km^3");
    CHECK(close_to(km3_to_ton(1.0), 1.0e9, 1e-12),
          "1 km^3 equivale a 10^9 toneladas de agua");
    CHECK(close_to(ton_to_km3(0.05), 5.0e-11, 1e-12),
          "0.05 t de agua equivalen a 5e-11 km^3");
}

/* Tarea 1: restricciones de los parámetros */
static void test_params_valid(void)
{
    Params p = sample_params();
    CHECK(params_valid(&p) == 1, "el conjunto base de parámetros es válido");

    Params q = p;
    q.eta = (p.alpha / p.beta) * 1.01;  /* rompe alpha >= beta*eta */
    CHECK(params_valid(&q) == 0, "se rechaza eta demasiado grande");

    q = p;
    q.xi = (p.beta / p.alpha) * 1.01;   /* rompe beta >= alpha*xi */
    CHECK(params_valid(&q) == 0, "se rechaza xi demasiado grande");

    q = p;
    q.r = -1.0;
    CHECK(params_valid(&q) == 0, "se rechaza una tasa de crecimiento negativa");

    q = p;
    q.K = 0.0;
    CHECK(params_valid(&q) == 0, "se rechaza una constante de saturación nula");
}

/* Tarea 2: la respuesta de crecimiento es saturante tipo Michaelis-Menten */
static void test_growth(void)
{
    Params p = sample_params();

    CHECK(close_to(growth(0.0, &p), 0.0, 1e-15),
          "sin agua no hay crecimiento");
    CHECK(close_to(growth(p.K, &p), p.r / 2.0, 1e-12),
          "en W = K el crecimiento es la mitad del máximo");
    CHECK(growth(1.0e6, &p) < growth(1.0e8, &p),
          "el crecimiento crece con W");
    CHECK(growth(1.0e30, &p) > 0.99 * p.r && growth(1.0e30, &p) <= p.r,
          "el crecimiento satura en r");
}

/* Tarea 3: condiciones iniciales */
static void test_initial_state(void)
{
    Params p = sample_params();
    State s;

    initial_state(&s, &p);

    CHECK(close_to(s.W, p.W0, 1e-12), "W(0) es el agua líquida inicial");
    CHECK(close_to(s.T, p.m_trucha, 1e-12), "T(0) es una trucha en toneladas");
    CHECK(s.G == 0.0, "G(0) es cero");
    CHECK(close_to(total_water(&p, &s), p.W0 + p.alpha * p.m_trucha, 1e-12),
          "el agua total inicial incluye la retenida por la trucha");
}

/* Tarea 2: conservación de la masa de agua */
static void test_water_conservation(void)
{
    Params p = sample_params();
    State s;
    double dW, dT, dG;

    s.W = 5.0e8;
    s.T = 1000.0;
    s.G = 200.0;
    rhs(&p, &s, &dW, &dT, &dG);

    CHECK(fabs(dW + p.alpha * dT + p.beta * dG) < 1e-18,
          "d/dt (W + alpha*T + beta*G) = 0");
}

/* Tarea 4: invariancia del dominio y cota superior de W */
static void test_invariance(void)
{
    Params p = sample_params();
    State s;
    double dW, dT, dG;

    s.W = 0.0; s.T = 1000.0; s.G = 200.0;
    rhs(&p, &s, &dW, &dT, &dG);
    CHECK(dW >= -1e-18, "en W = 0 el campo no apunta a W negativo");

    s.W = 1.0e8; s.T = 0.0; s.G = 200.0;
    rhs(&p, &s, &dW, &dT, &dG);
    CHECK(dT >= -1e-18, "en T = 0 el campo no apunta a T negativo");

    s.W = 1.0e8; s.T = 1000.0; s.G = 0.0;
    rhs(&p, &s, &dW, &dT, &dG);
    CHECK(dG >= -1e-18, "en G = 0 el campo no apunta a G negativo");

    s.W = p.W0 + p.alpha * p.m_trucha; s.T = 0.0; s.G = 0.0;
    CHECK(s.W <= p.W0 + p.alpha * p.m_trucha,
          "W está acotado por el agua total cuando T,G >= 0");
}

/* Tarea 5: la reducción a 2D reproduce el sistema completo */
static void test_reduction(void)
{
    Params p = sample_params();
    State s;
    double T = 1000.0, G = 200.0;
    double dW, dT_full, dG_full, dT_red, dG_red;

    state_from_reduced(&p, T, G, &s);

    CHECK(close_to(s.W, p.W0 + p.alpha * p.m_trucha - p.alpha * T - p.beta * G, 1e-12),
          "la reconstrucción de W usa la cantidad conservada");
    CHECK(close_to(total_water(&p, &s), p.W0 + p.alpha * p.m_trucha, 1e-12),
          "la reconstrucción preserva el agua total");

    rhs(&p, &s, &dW, &dT_full, &dG_full);
    rhs_reduced(&p, T, G, &dT_red, &dG_red);

    CHECK(close_to(dT_full, dT_red, 1e-12), "dT/dt coincide en ambos sistemas");
    CHECK(close_to(dG_full, dG_red, 1e-12), "dG/dt coincide en ambos sistemas");
}

int main(void)
{
    test_units();
    test_params_valid();
    test_growth();
    test_initial_state();
    test_water_conservation();
    test_invariance();
    test_reduction();

    if (failures == 0) {
        printf("OK: todas las pruebas del bloque 1 pasaron\n");
        return 0;
    }
    printf("%d prueba(s) fallaron\n", failures);
    return 1;
}
