/*
 * Visualizacion del bloque 6: datos y figuras SVG.
 *
 * Las series se integran con RK4 de paso fijo y el cruce del umbral se
 * localiza por biseccion con el mismo operador RK4. Se trabaja en variables
 * reducidas (T,G) y se reconstruye W desde la cantidad conservada, de modo
 * que W + alpha*T + beta*G se conserve por construccion.
 */
#include "visualization.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "collapse.h"
#include "parameters.h"
#include "svg.h"

/* Paso objetivo: ~1000 muestras para el conjunto base (t_colap ~ 4024 anios).
 * Se estima con la cuadratura para no agotar SERIES_MAX en celdas lentas. */
#define SERIES_STEP 4.0
#define SERIES_STEP_MAX 50.0
#define T_MAX_SEARCH 2.0e5

/* Zoom de r donde el evento existe para todo K de la tabla. Con K = 1.4e9,
 * hace falta r > c*K/(K+W_event) ~ 1.2e-3 para que g(W_event) > c. */
#define HEAT_R_MIN 0.002
#define HEAT_R_MAX 0.05

static void rk4_reduced(const Params *p, double T, double G, double h,
                        double *T_out, double *G_out)
{
    double a1, a2, b1, b2, c1, c2, d1, d2;

    rhs_reduced(p, T, G, &a1, &a2);
    rhs_reduced(p, T + 0.5 * h * a1, G + 0.5 * h * a2, &b1, &b2);
    rhs_reduced(p, T + 0.5 * h * b1, G + 0.5 * h * b2, &c1, &c2);
    rhs_reduced(p, T + h * c1, G + h * c2, &d1, &d2);

    *T_out = T + (h / 6.0) * (a1 + 2.0 * b1 + 2.0 * c1 + d1);
    *G_out = G + (h / 6.0) * (a2 + 2.0 * b2 + 2.0 * c2 + d2);
}

static void append_sample(Series *s, const Params *p, double t, double T,
                          double G)
{
    double W_tot = p->W0 + p->alpha * p->m_trucha;
    int n = s->n;

    if (n >= SERIES_MAX) {
        return;
    }
    s->t[n] = t;
    s->T[n] = T;
    s->G[n] = G;
    s->W[n] = W_tot - p->alpha * T - p->beta * G;
    s->n = n + 1;
}

int series_build(const Params *p, Series *s)
{
    double W_event, W_tot, S_event;
    double T, G;
    double t = 0.0, h = SERIES_STEP;
    double t_est, T_next, G_next, S_next;

    if (p == NULL || s == NULL || !params_valid(p)) {
        return 1;
    }
    W_event = collapse_threshold_water(p);
    W_tot = p->W0 + p->alpha * p->m_trucha;
    S_event = W_tot - W_event;
    T = p->m_trucha;
    G = 0.0;
    t_est = quadrature_t_colap(p);
    if (isfinite(t_est) && t_est > 0.0) {
        h = fmax(1.0, fmin(SERIES_STEP_MAX, t_est / 1000.0));
    }
    s->n = 0;
    append_sample(s, p, 0.0, T, G);

    if (p->alpha * p->m_trucha >= S_event) {
        return 0;
    }

    while (t < T_MAX_SEARCH) {
        double f_here = 0.0, f_hi = 1.0;

        if (h > T_MAX_SEARCH - t) {
            h = T_MAX_SEARCH - t;
        }
        rk4_reduced(p, T, G, h, &T_next, &G_next);
        S_next = p->alpha * T_next + p->beta * G_next;

        if (S_next >= S_event) {
            double T_event = T_next, G_event = G_next;
            int k;

            /* Misma relacion de Richardson que en el bloque 4: se refina el
             * subpaso que cruza para que el ultimo punto caiga en el evento. */
            for (k = 0; k < 60; k++) {
                double f_mid = 0.5 * (f_here + f_hi);
                double T_mid, G_mid, S_mid;

                rk4_reduced(p, T, G, f_mid * h, &T_mid, &G_mid);
                S_mid = p->alpha * T_mid + p->beta * G_mid;
                if (S_mid < S_event) {
                    f_here = f_mid;
                } else {
                    f_hi = f_mid;
                    T_event = T_mid;
                    G_event = G_mid;
                }
            }
            f_here = f_hi;
            rk4_reduced(p, T, G, f_here * h, &T_event, &G_event);
            if (t + f_here * h <= t) {
                return 1;
            }
            append_sample(s, p, t + f_here * h, T_event, G_event);
            return 0;
        }

        t += h;
        T = T_next;
        G = G_next;
        append_sample(s, p, t, T, G);
        if (s->n >= SERIES_MAX) {
            return 1;
        }
    }
    return 1;
}

int phase_build(const Params *p, PhaseData *d)
{
    Series s;
    Equilibrium desert, persist;
    double S_event, i;
    int k, n_sep = 64;

    if (p == NULL || d == NULL || !params_valid(p)) {
        return 1;
    }
    if (series_build(p, &s) != 0) {
        return 1;
    }

    memset(d, 0, sizeof(*d));
    d->trajectory.n = s.n;
    for (k = 0; k < s.n && k < PHASE_MAX; k++) {
        d->trajectory.T[k] = s.T[k];
        d->trajectory.G[k] = s.G[k];
    }
    d->trajectory.n = (s.n < PHASE_MAX) ? s.n : PHASE_MAX;

    (void)find_equilibria(p, &desert, &persist);
    d->desert.T = desert.T;
    d->desert.G = desert.G;
    d->desert.W = desert.W;
    d->persist.T = persist.T;
    d->persist.G = persist.G;
    d->persist.W = persist.W;

    /*
     * El desierto es silla, pero su variedad estable sale del cuadrante
     * fisico (el autovector estable tiene G < 0): no hay frontera de cuenca
     * interior. La "separatriz" dibujada es la recta operativa que define el
     * evento W = W_event, es decir, alpha*T + beta*G = W_tot - W_event.
     */
    S_event = (p->W0 + p->alpha * p->m_trucha) - collapse_threshold_water(p);
    if (persist.T > 0.0 && persist.G > 0.0) {
        /* Se recorta la recta a G/G* en [0.12, 1.05] para que la trayectoria
         * no quede pegada al eje y el extremo siga dentro del cuadrante. */
        double G_hi = 1.05 * persist.G;
        double G_lo = 0.12 * persist.G;

        for (k = 0; k < n_sep && k < PHASE_MAX; k++) {
            i = (double)k / (double)(n_sep - 1);
            double G_sep = G_hi + i * (G_lo - G_hi);
            double T_sep = (S_event - p->beta * G_sep) / p->alpha;

            d->separatrix.T[k] = T_sep;
            d->separatrix.G[k] = G_sep;
        }
        d->separatrix.n = n_sep;
    }
    return 0;
}

static int parameter_range(const char *symbol, double *lo, double *hi)
{
    size_t i;

    for (i = 0; i < PARAMETERS_COUNT; i++) {
        if (strcmp(PARAMETERS[i].symbol, symbol) == 0) {
            *lo = PARAMETERS[i].min;
            *hi = PARAMETERS[i].max;
            return 0;
        }
    }
    return 1;
}

int heatmap_build(const Params *p, Heatmap *h)
{
    Params cell;
    double K_lo, K_hi;
    int i, j;

    if (p == NULL || h == NULL || !params_valid(p)) {
        return 1;
    }
    if (parameter_range("K", &K_lo, &K_hi) != 0) {
        return 1;
    }
    /* Se usa el mismo conjunto base y solo se varia r y K. */
    cell = *p;
    h->nr = HEAT_R;
    h->nk = HEAT_K;

    for (i = 0; i < HEAT_R; i++) {
        h->r[i] = HEAT_R_MIN +
                  (HEAT_R_MAX - HEAT_R_MIN) * i / (double)(HEAT_R - 1);
    }
    for (j = 0; j < HEAT_K; j++) {
        h->K[j] = K_lo + (K_hi - K_lo) * j / (double)(HEAT_K - 1);
    }
    for (i = 0; i < HEAT_R; i++) {
        for (j = 0; j < HEAT_K; j++) {
            cell.r = h->r[i];
            cell.K = h->K[j];
            h->t[i][j] = quadrature_t_colap(&cell);
            if (!isfinite(h->t[i][j]) || h->t[i][j] <= 0.0) {
                h->t[i][j] = T_MAX_SEARCH;
            }
        }
    }
    return 0;
}

double heatmap_at(const Heatmap *h, int i, int j)
{
    if (h == NULL || i < 0 || i >= h->nr || j < 0 || j >= h->nk) {
        return 0.0;
    }
    return h->t[i][j];
}

int write_series_svg(const Series *s, const char *path)
{
    const double x0 = 80.0, y0 = 55.0, w = 760.0, hh = 370.0;
    double px[SERIES_MAX], py[SERIES_MAX];
    double t_end, W_max, T_max, G_max, W_ref, T_ref, G_ref;
    char label[128];
    Svg svg;
    int i;

    if (s == NULL || path == NULL || s->n < 2) {
        return 1;
    }
    t_end = s->t[s->n - 1];
    W_max = s->W[0];
    T_max = s->T[0];
    G_max = s->G[0];
    for (i = 0; i < s->n; i++) {
        T_max = fmax(T_max, s->T[i]);
        G_max = fmax(G_max, s->G[i]);
    }
    if (!(t_end > 0.0) || !(W_max > 0.0) || !(T_max > 0.0) ||
        !(G_max > 0.0)) {
        return 1;
    }
    if (svg_open(&svg, path, 900.0, 490.0,
                 "Series temporales W(t), T(t), G(t)") != 0) {
        return 1;
    }

    svg_rect(&svg, x0, y0, w, hh, "#fbfbfb");
    svg_text(&svg, x0 + w / 2.0, 28.0, "Series temporales y evento de desierto",
             18.0, "middle", "#111111");
    svg_line(&svg, x0, y0 + hh, x0 + w, y0 + hh, "#333333", 1.5);
    svg_line(&svg, x0, y0, x0, y0 + hh, "#333333", 1.5);

    for (i = 0; i < s->n; i++) {
        px[i] = x0 + w * s->t[i] / t_end;
    }
    W_ref = W_max;
    T_ref = T_max;
    G_ref = G_max;
    for (i = 0; i < s->n; i++) {
        py[i] = y0 + hh * (1.0 - s->W[i] / W_ref);
    }
    svg_polyline(&svg, px, py, s->n, "#1f77b4", 2.5);
    for (i = 0; i < s->n; i++) {
        py[i] = y0 + hh * (1.0 - s->T[i] / T_ref);
    }
    svg_polyline(&svg, px, py, s->n, "#d62728", 2.0);
    for (i = 0; i < s->n; i++) {
        py[i] = y0 + hh * (1.0 - s->G[i] / G_ref);
    }
    svg_polyline(&svg, px, py, s->n, "#2ca02c", 2.0);

    svg_line(&svg, x0 + w, y0, x0 + w, y0 + hh, "#7a3db8", 1.5);
    snprintf(label, sizeof(label), "t_colap = %.2f anios", t_end);
    svg_text(&svg, x0 + w - 6.0, y0 + 16.0, label, 13.0, "end", "#7a3db8");
    svg_text(&svg, x0 + w / 2.0, y0 + hh + 34.0, "tiempo (anios)", 14.0,
             "middle", "#111111");
    svg_text(&svg, x0 - 12.0, y0 + hh + 4.0, "0", 12.0, "end", "#111111");
    svg_text(&svg, x0 - 12.0, y0 + 10.0, "max", 12.0, "end", "#111111");

    svg_line(&svg, x0 + 12.0, y0 + 18.0, x0 + 42.0, y0 + 18.0, "#1f77b4", 2.5);
    svg_text(&svg, x0 + 50.0, y0 + 22.0, "W / W0", 13.0, "start", "#111111");
    svg_line(&svg, x0 + 150.0, y0 + 18.0, x0 + 180.0, y0 + 18.0, "#d62728", 2.0);
    svg_text(&svg, x0 + 188.0, y0 + 22.0, "T / max", 13.0, "start", "#111111");
    svg_line(&svg, x0 + 300.0, y0 + 18.0, x0 + 330.0, y0 + 18.0, "#2ca02c", 2.0);
    svg_text(&svg, x0 + 338.0, y0 + 22.0, "G / max", 13.0, "start", "#111111");

    svg_close(&svg);
    return 0;
}

int write_phase_svg(const PhaseData *d, const char *path)
{
    const double x0 = 90.0, y0 = 55.0, w = 700.0, hh = 430.0;
    const double margin = 1.08;
    double px[PHASE_MAX], py[PHASE_MAX];
    double max_T, max_G;
    char label[96];
    Svg svg;
    int i;

    if (d == NULL || path == NULL || d->trajectory.n < 2) {
        return 1;
    }
    max_T = d->persist.T;
    max_G = d->persist.G;
    for (i = 0; i < d->trajectory.n; i++) {
        max_T = fmax(max_T, d->trajectory.T[i]);
        max_G = fmax(max_G, d->trajectory.G[i]);
    }
    /* La escala la fijan los equilibrios y la trayectoria; la separatriz se
     * recorta en phase_build para caber en este mismo rectangulo. */
    max_T *= margin;
    max_G *= margin;
    if (!(max_T > 0.0) || !(max_G > 0.0)) {
        return 1;
    }
    if (svg_open(&svg, path, 860.0, 540.0,
                 "Retrato de fase reducido (T,G)") != 0) {
        return 1;
    }

    svg_rect(&svg, x0, y0, w, hh, "#fbfbfb");
    svg_text(&svg, x0 + w / 2.0, 28.0, "Retrato de fase: trayectoria y separatriz operativa",
             18.0, "middle", "#111111");
    svg_line(&svg, x0, y0 + hh, x0 + w, y0 + hh, "#333333", 1.5);
    svg_line(&svg, x0, y0, x0, y0 + hh, "#333333", 1.5);

    for (i = 0; i < d->separatrix.n; i++) {
        px[i] = x0 + w * d->separatrix.T[i] / max_T;
        py[i] = y0 + hh * (1.0 - d->separatrix.G[i] / max_G);
    }
    svg_polyline(&svg, px, py, d->separatrix.n, "#7a3db8", 2.5);
    for (i = 0; i < d->trajectory.n; i++) {
        px[i] = x0 + w * d->trajectory.T[i] / max_T;
        py[i] = y0 + hh * (1.0 - d->trajectory.G[i] / max_G);
    }
    svg_polyline(&svg, px, py, d->trajectory.n, "#1f77b4", 2.5);

    svg_circle(&svg, x0, y0 + hh, 5.0, "#d62728");
    svg_text(&svg, x0 + 10.0, y0 + hh - 8.0, "desierto", 13.0, "start", "#d62728");
    if (d->persist.T > 0.0 && d->persist.G > 0.0) {
        double xp = x0 + w * d->persist.T / max_T;
        double yp = y0 + hh * (1.0 - d->persist.G / max_G);
        svg_circle(&svg, xp, yp, 6.0, "#d62728");
        svg_text(&svg, xp - 10.0, yp - 12.0, "persistencia", 13.0, "end", "#d62728");
    }
    svg_text(&svg, x0 + w / 2.0, y0 + hh + 34.0, "T / max_T (escala lineal)", 14.0,
             "middle", "#111111");
    svg_text(&svg, x0 - 16.0, y0 + 14.0, "G / max_G", 13.0, "end", "#111111");
    snprintf(label, sizeof(label), "max_T = %.2e t;  max_G = %.2e t",
             max_T, max_G);
    svg_text(&svg, x0 + w, y0 + hh + 47.0, label, 11.0, "end", "#555555");
    svg_line(&svg, x0 + 14.0, y0 + 22.0, x0 + 44.0, y0 + 22.0, "#1f77b4", 2.5);
    svg_text(&svg, x0 + 52.0, y0 + 26.0, "trayectoria", 13.0, "start", "#111111");
    svg_line(&svg, x0 + 180.0, y0 + 22.0, x0 + 210.0, y0 + 22.0, "#7a3db8", 2.5);
    svg_text(&svg, x0 + 218.0, y0 + 26.0, "separatriz (W = W_event)", 13.0,
             "start", "#111111");

    svg_close(&svg);
    return 0;
}

static void heat_color(double q, char out[8])
{
    int red, green, blue;

    if (!isfinite(q) || q < 0.0) {
        q = 0.0;
    }
    if (q > 1.0) {
        q = 1.0;
    }
    red = (int)(255.0 * q);
    green = (int)(90.0 * (1.0 - fabs(q - 0.5) * 2.0));
    blue = (int)(255.0 * (1.0 - q));
    snprintf(out, 8, "#%02x%02x%02x", red, green, blue);
}

int write_heatmap_svg(const Heatmap *h, const char *path)
{
    const double x0 = 95.0, y0 = 60.0, w = 560.0, hh = 440.0;
    double lo = INFINITY, hi = -INFINITY;
    char color[8], label[96];
    Svg svg;
    int i, j;

    if (h == NULL || path == NULL || h->nr <= 0 || h->nk <= 0) {
        return 1;
    }
    for (i = 0; i < h->nr; i++) {
        for (j = 0; j < h->nk; j++) {
            double v = heatmap_at(h, i, j);
            if (isfinite(v) && v > 0.0) {
                lo = fmin(lo, log10(v));
                hi = fmax(hi, log10(v));
            }
        }
    }
    if (!isfinite(lo) || !isfinite(hi) || hi <= lo) {
        return 1;
    }
    if (svg_open(&svg, path, 820.0, 585.0,
                 "Mapa de calor de t_colap cuasi-estacionario frente a r y K") != 0) {
        return 1;
    }

    svg_text(&svg, 410.0, 28.0,
             "Mapa de calor: t_colap cuasi-estacionario (anios) vs r y K", 18.0,
             "middle", "#111111");
    for (i = 0; i < h->nr; i++) {
        for (j = 0; j < h->nk; j++) {
            double q = (log10(heatmap_at(h, i, j)) - lo) / (hi - lo);
            double cx = x0 + w * j / (double)h->nk;
            double cy = y0 + hh * i / (double)h->nr;
            heat_color(q, color);
            svg_rect(&svg, cx, cy, w / (double)h->nk, hh / (double)h->nr + 0.2,
                     color);
        }
    }
    svg_rect(&svg, x0, y0, w, hh, "none");
    svg_text(&svg, x0 + w / 2.0, y0 + hh + 35.0, "K (km^3)", 14.0, "middle",
             "#111111");
    svg_text(&svg, x0 - 12.0, y0 + 14.0, "r", 14.0, "end", "#111111");
    snprintf(label, sizeof(label), "%.3f", h->r[0]);
    svg_text(&svg, x0 - 12.0, y0 + 10.0, label, 11.0, "end", "#555555");
    snprintf(label, sizeof(label), "%.3f", h->r[h->nr - 1]);
    svg_text(&svg, x0 - 12.0, y0 + hh, label, 11.0, "end", "#555555");
    snprintf(label, sizeof(label), "%.2e", h->K[0]);
    svg_text(&svg, x0, y0 + hh + 18.0, label, 11.0, "start", "#555555");
    snprintf(label, sizeof(label), "%.2e", h->K[h->nk - 1]);
    svg_text(&svg, x0 + w, y0 + hh + 18.0, label, 11.0, "end", "#555555");

    /* Escala de color: se dibuja con la misma funcion que las celdas. */
    for (i = 0; i < 40; i++) {
        double q = i / 39.0;
        heat_color(q, color);
        svg_rect(&svg, x0 + w + 34.0, y0 + hh - (i + 1) * hh / 40.0,
                 22.0, hh / 40.0 + 0.2, color);
    }
    snprintf(label, sizeof(label), "%.0f", pow(10.0, hi));
    svg_text(&svg, x0 + w + 62.0, y0 + 12.0, label, 11.0, "start", "#111111");
    snprintf(label, sizeof(label), "%.0f", pow(10.0, lo));
    svg_text(&svg, x0 + w + 62.0, y0 + hh, label, 11.0, "start", "#111111");
    svg_text(&svg, x0 + w + 34.0, y0 - 12.0, "t_colap", 12.0, "start", "#111111");
    svg_text(&svg, 410.0, 565.0,
             "escala logaritmica; cada celda usa el par (r,K) de su fila y columna",
             12.0, "middle", "#555555");

    svg_close(&svg);
    return 0;
}
