/*
 * Genera las tres figuras SVG del bloque 6 con los parametros base.
 * No requiere Python, gnuplot ni bibliotecas externas.
 */
#include <stdio.h>

#include "parameters.h"
#include "visualization.h"

int main(void)
{
    Params p;
    Series series;
    PhaseData phase;
    Heatmap heatmap;
    int rc = 0;

    params_default(&p);

    if (series_build(&p, &series) != 0) {
        fprintf(stderr, "No se pudo construir la serie temporal\n");
        return 1;
    }
    if (phase_build(&p, &phase) != 0) {
        fprintf(stderr, "No se pudo construir el retrato de fase\n");
        return 1;
    }
    if (heatmap_build(&p, &heatmap) != 0) {
        fprintf(stderr, "No se pudo construir el mapa de calor\n");
        return 1;
    }

    rc |= write_series_svg(&series, "figures/series.svg");
    rc |= write_phase_svg(&phase, "figures/retrato_fase.svg");
    rc |= write_heatmap_svg(&heatmap, "figures/heatmap.svg");

    if (rc != 0) {
        fprintf(stderr, "Fallo al escribir alguna figura\n");
        return 1;
    }
    printf("Figuras escritas: figures/series.svg, figures/retrato_fase.svg, "
           "figures/heatmap.svg\n");
    return 0;
}
