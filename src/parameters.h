#ifndef EDO_PARAMETERS_H
#define EDO_PARAMETERS_H

#include <stddef.h>

#include "model.h"

/* Procedencia de cada valor de la tabla. */
typedef enum {
    PARAM_DATA,        /* medición o compilación con fuente citada */
    PARAM_ASSUMPTION,  /* supuesto ficticio declarado */
    PARAM_CONVENTION   /* convención derivada de datos, no dato directo */
} ParamOrigin;

typedef struct {
    const char *symbol;
    const char *unit;
    double value;
    double min;
    double max;
    ParamOrigin origin;
    const char *source;
    const char *note;
} ParamMeta;

/* Tabla única de verdad para el README y las pruebas. */
extern const ParamMeta PARAMETERS[];
extern const size_t PARAMETERS_COUNT;

/* Valores por defecto del modelo, tomados de PARAMETERS. */
void params_default(Params *p);

#endif /* EDO_PARAMETERS_H */
