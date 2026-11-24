#include "svg.h"

#include <math.h>

static void svg_escape(FILE *f, const char *text)
{
    const unsigned char *p;

    for (p = (const unsigned char *)text; *p != '\0'; p++) {
        switch (*p) {
        case '<':  fputs("&lt;", f); break;
        case '>':  fputs("&gt;", f); break;
        case '&':  fputs("&amp;", f); break;
        case '"':  fputs("&quot;", f); break;
        case '\'': fputs("&apos;", f); break;
        default:   fputc(*p, f); break;
        }
    }
}

int svg_open(Svg *s, const char *path, double width, double height,
             const char *title)
{
    s->f = fopen(path, "wb");
    if (s->f == NULL) {
        return 1;
    }
    s->width = width;
    s->height = height;

    fputs("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n", s->f);
    fprintf(s->f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" "
                  "height=\"%.0f\" viewBox=\"0 0 %.0f %.0f\">\n",
            width, height, width, height);
    fputs("<title>", s->f);
    svg_escape(s->f, title);
    fputs("</title>\n", s->f);
    fprintf(s->f, "<rect x=\"0\" y=\"0\" width=\"%.0f\" height=\"%.0f\" "
                  "fill=\"white\"/>\n", width, height);
    return 0;
}

void svg_close(Svg *s)
{
    if (s->f != NULL) {
        fputs("</svg>\n", s->f);
        fclose(s->f);
        s->f = NULL;
    }
}

void svg_rect(Svg *s, double x, double y, double w, double h,
              const char *fill)
{
    fprintf(s->f, "<rect x=\"%.4f\" y=\"%.4f\" width=\"%.4f\" "
                  "height=\"%.4f\" fill=\"%s\"/>\n",
            x, y, w, h, fill);
}

void svg_line(Svg *s, double x1, double y1, double x2, double y2,
              const char *color, double stroke)
{
    fprintf(s->f, "<line x1=\"%.4f\" y1=\"%.4f\" x2=\"%.4f\" y2=\"%.4f\" "
                  "stroke=\"%s\" stroke-width=\"%.3f\"/>\n",
            x1, y1, x2, y2, color, stroke);
}

void svg_polyline(Svg *s, const double *xs, const double *ys, int n,
                  const char *color, double stroke)
{
    int i;

    if (n <= 0) {
        return;
    }
    fprintf(s->f, "<polyline fill=\"none\" stroke=\"%s\" "
                  "stroke-width=\"%.3f\" points=\"", color, stroke);
    for (i = 0; i < n; i++) {
        double x = isfinite(xs[i]) ? xs[i] : 0.0;
        double y = isfinite(ys[i]) ? ys[i] : 0.0;
        fprintf(s->f, "%.3f,%.3f%s", x, y, i + 1 < n ? " " : "");
    }
    fputs("\"/>\n", s->f);
}

void svg_circle(Svg *s, double x, double y, double r, const char *fill)
{
    fprintf(s->f, "<circle cx=\"%.4f\" cy=\"%.4f\" r=\"%.4f\" "
                  "fill=\"%s\"/>\n", x, y, r, fill);
}

void svg_text(Svg *s, double x, double y, const char *text, double size,
              const char *anchor, const char *fill)
{
    fprintf(s->f, "<text x=\"%.4f\" y=\"%.4f\" font-size=\"%.2f\" "
                  "text-anchor=\"%s\" fill=\"%s\">",
            x, y, size, anchor, fill);
    svg_escape(s->f, text);
    fputs("</text>\n", s->f);
}
