/*
 * Pruebas del bloque 6: visualizacion y entrega.
 * Ejecutar con: make test
 *
 * Fase TDD: este archivo se escribio antes que src/visualization.c.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "analysis.h"
#include "collapse.h"
#include "parameters.h"
#include "visualization.h"

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

static int file_contains(const char *path, const char *needle)
{
    char buffer[65536];
    size_t n;
    FILE *f = fopen(path, "rb");

    if (f == NULL) {
        return 0;
    }
    n = fread(buffer, 1, sizeof(buffer) - 1, f);
    buffer[n] = '\0';
    fclose(f);
    return strstr(buffer, needle) != NULL;
}

static int xml_well_formed(const char *path)
{
    char buffer[1 << 20];
    char stack[64][64];
    size_t n, pos = 0;
    int depth = 0;
    FILE *f = fopen(path, "rb");

    if (f == NULL) {
        return 0;
    }
    n = fread(buffer, 1, sizeof(buffer) - 1, f);
    buffer[n] = '\0';
    fclose(f);

    while (pos < n) {
        char *open = strchr(buffer + pos, '<');
        char *close;
        char *end;
        char name[64];
        size_t len, k;
        int self_closing;

        if (open == NULL) {
            break;
        }
        close = strchr(open, '>');
        if (close == NULL) {
            return 0;
        }
        end = close - 1;
        while (end > open && (*end == ' ' || *end == '\t')) {
            end--;
        }
        self_closing = (*end == '/');
        if (open[1] == '?' || open[1] == '!') {
            pos = (size_t)(close - buffer) + 1;
            continue;
        }

        end = open + 1;
        if (*end == '/') {
            end++;
        }
        while (*end == ' ' || *end == '\t') {
            end++;
        }
        len = 0;
        while (end + len < close && end[len] > ' ' && end[len] != '/' &&
               len + 1 < sizeof(name)) {
            name[len] = end[len];
            len++;
        }
        name[len] = '\0';
        if (len == 0) {
            return 0;
        }

        if (open[1] == '/') {
            if (depth == 0 || strcmp(stack[depth - 1], name) != 0) {
                return 0;
            }
            depth--;
        } else if (!self_closing) {
            if (depth >= (int)(sizeof(stack) / sizeof(stack[0]))) {
                return 0;
            }
            for (k = 0; k <= len; k++) {
                stack[depth][k] = name[k];
            }
            depth++;
        }
        pos = (size_t)(close - buffer) + 1;
    }
    return depth == 0;
}

static int file_is_svg(const char *path, const char *label)
{
    return file_contains(path, "<svg") &&
           file_contains(path, "</svg>") &&
           file_contains(path, label) &&
           !file_contains(path, "nan") &&
           !file_contains(path, "inf") &&
           xml_well_formed(path);
}

static void test_series(void)
{
    Params p;
    Series s;
    double W_tot, event;
    int i, ok;

    params_default(&p);
    ok = series_build(&p, &s);
    CHECK(ok == 0 && s.n > 100 && s.n <= SERIES_MAX,
          "la serie tiene cientos de muestras");
    if (!ok) {
        return;
    }

    W_tot = p.W0 + p.alpha * p.m_trucha;
    event = collapse_threshold_water(&p);
    CHECK(s.t[0] == 0.0 && close_rel(s.W[0], p.W0, 1e-12) &&
          close_rel(s.T[0], p.m_trucha, 1e-12) && s.G[0] == 0.0,
          "la serie empieza en la condicion inicial");
    CHECK(close_rel(s.W[s.n - 1], event, 1e-6),
          "la serie termina en el evento de desierto");
    CHECK(fabs(s.t[s.n - 1] - 4024.278215360) < 0.01,
          "la serie reproduce la referencia RK4 de t_colap");

    for (i = 1; i < s.n; i++) {
        CHECK(s.t[i] > s.t[i - 1], "el tiempo de la serie crece");
        if (s.W[i] > s.W[i - 1] + 1e-3) {
            CHECK(0, "W no crece a lo largo de la serie");
            break;
        }
    }
    for (i = 0; i < s.n; i += 37) {
        CHECK(close_rel(s.W[i] + p.alpha * s.T[i] + p.beta * s.G[i],
                        W_tot, 1e-12),
              "la serie conserva el agua total");
    }

    CHECK(file_is_svg("figures/series.svg", "t_colap"),
          "series.svg es un SVG con la marca del evento");
}

static void test_phase_portrait(void)
{
    Params p;
    PhaseData d;
    Equilibrium desert, persist;
    double Tx, Gx, W_tot;
    int i, ok;

    params_default(&p);
    (void)find_equilibria(&p, &desert, &persist);
    ok = phase_build(&p, &d);
    CHECK(ok == 0, "se construye el retrato de fase");
    if (!ok) {
        return;
    }

    CHECK(d.desert.T == 0.0 && d.desert.G == 0.0,
          "el desierto se marca en el origen del sistema reducido");
    CHECK(close_rel(d.persist.T, persist.T, 1e-9) &&
          close_rel(d.persist.G, persist.G, 1e-9),
          "la persistencia se marca en su equilibrio");
    CHECK(d.trajectory.n > 100 && d.separatrix.n > 20,
          "hay trayectoria y separatriz suficientes");

    Tx = d.separatrix.T[d.separatrix.n - 1] / persist.T;
    Gx = d.separatrix.G[d.separatrix.n - 1] / persist.G;
    CHECK(Tx > 0.1 && Gx > 0.1,
          "la separatriz se aleja del desierto hasta la zona de interes");
    CHECK(d.separatrix.G[0] <= 1.1 * persist.G,
          "la separatriz esta recortada al area de la trayectoria");

    W_tot = p.W0 + p.alpha * p.m_trucha;
    for (i = 0; i < d.separatrix.n; i += 7) {
        CHECK(close_rel(p.alpha * d.separatrix.T[i] +
                        p.beta * d.separatrix.G[i],
                        W_tot - collapse_threshold_water(&p), 1e-12),
              "la separatriz satisface W = W_event");
        CHECK(d.separatrix.T[i] >= 0.0 && d.separatrix.G[i] >= 0.0,
              "la separatriz se mantiene en el cuadrante fisico");
    }
    Gx = 0.0;
    for (i = 0; i < d.trajectory.n; i++) {
        Gx = fmax(Gx, d.trajectory.G[i]);
    }
    CHECK(Gx > 0.1 * persist.G,
          "la trayectoria ocupa una fraccion visible del retrato");

    CHECK(file_is_svg("figures/retrato_fase.svg", "separatriz"),
          "retrato_fase.svg es un SVG con separatriz y equilibrios");
}

static void test_heatmap(void)
{
    Params p;
    Heatmap h;
    int i, j, ok;

    params_default(&p);
    ok = heatmap_build(&p, &h);
    CHECK(ok == 0 && h.nr == HEAT_R && h.nk == HEAT_K,
          "el mapa de calor tiene la rejilla prevista");
    if (!ok) {
        return;
    }

    for (i = 0; i < h.nr; i++) {
        for (j = 0; j < h.nk; j++) {
            double t = heatmap_at(&h, i, j);
            CHECK(isfinite(t) && t > 0.0 && t < 2.0e5,
                  "t_colap de la celda es finito y no es el centinela");
            if (i + 1 < h.nr) {
                CHECK(t > heatmap_at(&h, i + 1, j),
                      "aumentar r acorta el colapso entre celdas vecinas");
            }
            if (j + 1 < h.nk) {
                CHECK(t < heatmap_at(&h, i, j + 1),
                      "aumentar K alarga el colapso entre celdas vecinas");
            }
        }
    }
    for (j = 0; j < h.nk; j++) {
        CHECK(heatmap_at(&h, 0, j) > heatmap_at(&h, h.nr - 1, j),
              "aumentar r acorta el colapso en toda la rejilla");
    }
    for (i = 0; i < h.nr; i++) {
        CHECK(heatmap_at(&h, i, 0) < heatmap_at(&h, i, h.nk - 1),
              "aumentar K alarga el colapso en toda la rejilla");
    }
    {
        Params cell = p;
        SimResult rk;
        int mid_i = HEAT_R / 2, mid_j = HEAT_K / 2;

        cell.r = h.r[mid_i];
        cell.K = h.K[mid_j];
        rk = simulate_until_collapse(&cell, 2.0e5, 1.0e-6);
        CHECK(rk.collapsed &&
              fabs(heatmap_at(&h, mid_i, mid_j) - rk.t_colap) / rk.t_colap < 0.12,
              "la cuadratura del mapa concuerda con RK4 en una celda central");
    }

    CHECK(file_is_svg("figures/heatmap.svg", "t_colap"),
          "heatmap.svg es un SVG con escala de t_colap");
}

static void test_null_inputs(void)
{
    Params p;
    Series s;
    PhaseData d;
    Heatmap h;

    params_default(&p);
    CHECK(series_build(NULL, &s) != 0 && series_build(&p, NULL) != 0,
          "series_build rechaza argumentos nulos");
    CHECK(phase_build(NULL, &d) != 0 && phase_build(&p, NULL) != 0,
          "phase_build rechaza argumentos nulos");
    CHECK(heatmap_build(NULL, &h) != 0 && heatmap_build(&p, NULL) != 0,
          "heatmap_build rechaza argumentos nulos");
}

static void test_deliverables(void)
{
    CHECK(file_contains("README.md", "## Bloque 6"),
          "README contiene el bloque 6");
    CHECK(file_contains("docs/informe.md", "limitaciones"),
          "el informe explica las limitaciones");
    CHECK(file_contains("docs/informe.md", "4024"),
          "el informe da la cifra central de t_colap");
    CHECK(file_contains("docs/exposicion.md", "figures/"),
          "la exposicion referencia las figuras generadas");
}

int main(void)
{
    test_series();
    test_phase_portrait();
    test_heatmap();
    test_null_inputs();
    test_deliverables();

    if (failures == 0) {
        printf("OK: todas las pruebas del bloque 6 pasaron\n");
        return 0;
    }
    printf("%d prueba(s) fallaron\n", failures);
    return 1;
}
