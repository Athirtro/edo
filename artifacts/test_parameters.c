/*
 * Pruebas del bloque 2: datos y parametros.
 * Ejecutar con: make test
 *
 * Fase TDD: este archivo se escribio antes que src/parameters.c.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

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

static const ParamMeta *find_param(const char *symbol)
{
    size_t i;
    for (i = 0; i < PARAMETERS_COUNT; i++) {
        if (strcmp(PARAMETERS[i].symbol, symbol) == 0) {
            return &PARAMETERS[i];
        }
    }
    return NULL;
}

/* Tareas 8 y 9: toda la tabla tiene metadatos completos y rangos coherentes */
static void test_metadata_complete(void)
{
    size_t i;

    CHECK(PARAMETERS_COUNT >= 12, "la tabla incluye al menos 12 parametros");
    for (i = 0; i < PARAMETERS_COUNT; i++) {
        const ParamMeta *m = &PARAMETERS[i];
        CHECK(m->symbol != NULL && m->symbol[0] != '\0', "simbolo no vacio");
        CHECK(m->unit != NULL && m->unit[0] != '\0', "unidad no vacia");
        CHECK(m->source != NULL && m->source[0] != '\0', "fuente no vacia");
        CHECK(m->note != NULL && m->note[0] != '\0', "nota no vacia");
        CHECK(m->value > 0.0, "valor positivo");
        CHECK(m->min > 0.0 && m->min < m->value, "minimo menor que el valor");
        CHECK(m->max > m->value, "maximo mayor que el valor");
        CHECK(m->origin == PARAM_DATA || m->origin == PARAM_ASSUMPTION ||
              m->origin == PARAM_CONVENTION, "origen enumerado");
    }
}

/* Tarea 8: simbolos requeridos por el modelo y por el umbral */
static void test_required_symbols(void)
{
    static const char *symbols[] = {
        "W0", "A", "f_desert", "r", "K", "mu", "delta",
        "alpha", "beta", "eta", "xi", "m_trucha"
    };
    size_t i;

    for (i = 0; i < sizeof(symbols) / sizeof(symbols[0]); i++) {
        CHECK(find_param(symbols[i]) != NULL, "parametro presente en la tabla");
    }
}

/* Tarea 6: datos reales citados */
static void test_data_values(void)
{
    const ParamMeta *w0 = find_param("W0");
    const ParamMeta *area = find_param("A");

    CHECK(w0 != NULL && close_to(w0->value, 1.36160761e9, 1e-9),
          "W0 es la suma de reservas liquidas USGS");
    CHECK(w0 != NULL && w0->origin == PARAM_DATA,
          "W0 se marca como dato");
    CHECK(w0 != NULL && strstr(w0->source, "usgs.gov") != NULL,
          "W0 cita USGS");

    CHECK(area != NULL && close_to(area->value, 5.10064472e8, 1e-9),
          "el area superficial es la de NASA/JPL");
    CHECK(area != NULL && area->origin == PARAM_DATA,
          "el area se marca como dato");
}

/* Tarea 7: convencion del umbral a partir de humedad del suelo */
static void test_threshold_convention(void)
{
    const ParamMeta *f = find_param("f_desert");

    CHECK(f != NULL && f->origin == PARAM_CONVENTION,
          "el umbral se declara como convencion");
    CHECK(f != NULL && f->value > 0.0 && f->value < 1.0,
          "el umbral es una fraccion");
    CHECK(f != NULL && close_to(f->value, 0.02 / 0.15, 1e-9),
          "f_desert = theta_Sahara / theta_referencia");
    CHECK(f != NULL && strstr(f->source, "doi.org") != NULL,
          "la convencion cita fuentes con DOI");
}

/* Tarea 9: supuestos ficticios registrados como parametros libres */
static void test_assumptions_registered(void)
{
    static const char *symbols[] = {
        "r", "K", "mu", "delta", "alpha", "beta", "eta", "xi", "m_trucha"
    };
    size_t i;
    const ParamMeta *mu;

    for (i = 0; i < sizeof(symbols) / sizeof(symbols[0]); i++) {
        const ParamMeta *m = find_param(symbols[i]);
        CHECK(m != NULL && m->origin == PARAM_ASSUMPTION,
              "parametro registrado como supuesto");
    }

    mu = find_param("mu");
    CHECK(mu != NULL && close_to(mu->value, 1.0 / 4000.0, 1e-12),
          "la maduracion por defecto es 1/4000 anios");
}

/* Los valores por defecto alimentan al modelo y respetan sus restricciones */
static void test_defaults(void)
{
    Params p;
    const ParamMeta *w0 = find_param("W0");
    const ParamMeta *f = find_param("f_desert");

    params_default(&p);

    CHECK(params_valid(&p) == 1, "los valores por defecto son validos");
    CHECK(w0 != NULL && close_to(p.W0, w0->value, 1e-12),
          "params_default toma W0 de la tabla");
    CHECK(f != NULL && close_to(p.f_desert, f->value, 1e-12),
          "params_default toma f_desert de la tabla");
    CHECK(close_to(p.mu, 1.0 / 4000.0, 1e-12),
          "params_default usa maduracion de 4000 anios");
    CHECK(p.alpha >= p.beta * p.eta && p.beta >= p.alpha * p.xi,
          "los supuestos cumplen la invariancia del bloque 1");
}

int main(void)
{
    test_metadata_complete();
    test_required_symbols();
    test_data_values();
    test_threshold_convention();
    test_assumptions_registered();
    test_defaults();

    if (failures == 0) {
        printf("OK: todas las pruebas del bloque 2 pasaron\n");
        return 0;
    }
    printf("%d prueba(s) fallaron\n", failures);
    return 1;
}
