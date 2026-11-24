#ifndef EDO_SVG_H
#define EDO_SVG_H

#include <stdio.h>

/* Escritor minimalista de SVG 1.1: suficiente para las figuras del proyecto
 * y sin dependencias externas. */
typedef struct {
    FILE *f;
    double width;
    double height;
} Svg;

int svg_open(Svg *s, const char *path, double width, double height,
             const char *title);
void svg_close(Svg *s);

void svg_rect(Svg *s, double x, double y, double w, double h,
              const char *fill);
void svg_line(Svg *s, double x1, double y1, double x2, double y2,
              const char *color, double stroke);
void svg_polyline(Svg *s, const double *xs, const double *ys, int n,
                  const char *color, double stroke);
void svg_circle(Svg *s, double x, double y, double r, const char *fill);
void svg_text(Svg *s, double x, double y, const char *text, double size,
              const char *anchor, const char *fill);

#endif /* EDO_SVG_H */
