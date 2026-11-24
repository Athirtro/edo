#include "parameters.h"

#include <math.h>
#include <string.h>

/*
 * Fuentes consultadas en el bloque 2:
 *  - USGS Water Science School (tabla de Shiklomanov/Gleick 1993):
 *    https://www.usgs.gov/water-science-school/science/how-much-water-there-earth
 *  - NASA/JPL, Compare Earth and the Moon:
 *    https://science.nasa.gov/moon/by-the-numbers/
 *  - Shao, Fraedrich & Ishizuka (2021), Boundary-Layer Meteorology:
 *    https://doi.org/10.1007/s10546-020-00596-9
 *  - Pan et al. (2019), Water 11(5):883:
 *    https://doi.org/10.3390/w11050883
 */

#define THETA_SAHARA 0.02
#define THETA_REFERENCE 0.15

const ParamMeta PARAMETERS[] = {
    {
        "W0", "km^3", 1.36160761e9, 1.33e9, 1.39e9, PARAM_DATA,
        "USGS Water Science School, https://www.usgs.gov/water-science-school/science/how-much-water-there-earth",
        "Suma de reservas liquidas: oceanos 1.338e9 + aguas subterraneas 2.34e7 + lagos, humedales, rios, humedad del suelo y agua biologica; excluye hielo y vapor."
    },
    {
        "A", "km^2", 5.10064472e8, 5.099e8, 5.101e8, PARAM_DATA,
        "NASA/JPL, https://science.nasa.gov/moon/by-the-numbers/",
        "Area superficial total de la Tierra; contexto para convertir volumenes en altura de capa."
    },
    {
        "f_desert", "adimensional", THETA_SAHARA / THETA_REFERENCE,
        0.05, 0.25, PARAM_CONVENTION,
        "Shao et al. 2021 https://doi.org/10.1007/s10546-020-00596-9 ; Pan et al. 2019 https://doi.org/10.3390/w11050883",
        "Convencion: theta_Sahara / theta_referencia = 0.02 / 0.15. No es un volumen literal; una Tierra totalmente seca seria mas arida que el Sahara."
    },
    {
        "r", "1/anio", 0.02, 0.001, 0.2, PARAM_ASSUMPTION,
        "supuesto del modelo (Herbert no da cifras)",
        "Tasa maxima de crecimiento de las truchas."
    },
    {
        "K", "km^3", 1.0e9, 1.0e8, 1.4e9, PARAM_ASSUMPTION,
        "supuesto del modelo",
        "Agua de media saturacion; del orden de W0."
    },
    {
        "mu", "1/anio", 1.0 / 4000.0, 1.0 / 10000.0, 1.0 / 1000.0,
        PARAM_ASSUMPTION, "supuesto del modelo (~4000 anios de maduracion)",
        "Tasa media de maduracion trucha -> gusano."
    },
    {
        "delta", "1/anio", 1.0 / 500.0, 1.0 / 5000.0, 1.0 / 100.0,
        PARAM_ASSUMPTION, "supuesto del modelo",
        "Tasa de mortalidad del gusano; tiempo medio de vida ~500 anios."
    },
    {
        "alpha", "km^3/t", 1.0e-9, 1.0e-10, 1.0e-8, PARAM_ASSUMPTION,
        "supuesto del modelo",
        "Agua retenida por tonelada de trucha (biomasa + reservas); 1e-9 km^3 = 1 m^3."
    },
    {
        "beta", "km^3/t", 5.0e-10, 5.0e-11, 5.0e-9, PARAM_ASSUMPTION,
        "supuesto del modelo",
        "Agua retenida por tonelada de gusano; debe satisfacer beta >= alpha*xi."
    },
    {
        "eta", "t/t", 1.5, 0.5, 2.0, PARAM_ASSUMPTION,
        "supuesto del modelo",
        "Biomasa de gusano producida por tonelada de trucha que madura; eta <= alpha/beta."
    },
    {
        "xi", "t/t", 0.3, 0.05, 0.5, PARAM_ASSUMPTION,
        "supuesto del modelo",
        "Biomasa de trucha liberada por tonelada de gusano que muere; xi <= beta/alpha."
    },
    {
        "m_trucha", "t", 0.05, 0.01, 0.2, PARAM_ASSUMPTION,
        "supuesto del modelo (una trucha de 50 kg)",
        "Masa de la unica trucha introducida en t(0)."
    }
};

const size_t PARAMETERS_COUNT = sizeof(PARAMETERS) / sizeof(PARAMETERS[0]);

static double value_of(const char *symbol)
{
    size_t i;
    for (i = 0; i < PARAMETERS_COUNT; i++) {
        if (strcmp(PARAMETERS[i].symbol, symbol) == 0) {
            return PARAMETERS[i].value;
        }
    }
    return NAN;
}

void params_default(Params *p)
{
    p->r = value_of("r");
    p->K = value_of("K");
    p->mu = value_of("mu");
    p->delta = value_of("delta");
    p->alpha = value_of("alpha");
    p->beta = value_of("beta");
    p->eta = value_of("eta");
    p->xi = value_of("xi");
    p->W0 = value_of("W0");
    p->m_trucha = value_of("m_trucha");
    p->f_desert = value_of("f_desert");
}
