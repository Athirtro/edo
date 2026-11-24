# Extracción de cifras (bloque 2)

Notas de trabajo sobre cómo se obtuvieron los valores que aparecen en el README. No sustituyen a las referencias, sólo las desglosan.

## [1] USGS — volumen de agua

Página: https://www.usgs.gov/water-science-school/science/how-much-water-there-earth

La tabla «Where is Earth's water located?» atribuye las cifras a Igor Shiklomanov, en Peter H. Gleick (ed.), *Water in Crisis: A Guide to the World's Fresh Water Resources* (Oxford University Press, 1993). Los valores en km³ son:

- Océanos, mares y bahías: 1 338 000 000
- Casquetes, glaciares y nieve permanente: 24 064 000 → **excluido** (sólido)
- Aguas subterráneas: 23 400 000 → **incluido**
  - dulces 10 530 000, salinas 12 870 000
- Humedad del suelo: 16 500 → **incluido**
- Hielo y permafrost: 300 000 → **excluido** (sólido)
- Lagos: 176 400 → **incluido** (dulces 91 000, salinos 85 400)
- Vapor atmosférico: 12 900 → **excluido** (gas)
- Humedales: 11 470 → **incluido**
- Ríos: 2 120 → **incluido**
- Agua biológica: 1 120 → **incluido**

Suma de líquidos: 1 361 607 610 km³. El total de **toda** el agua (incluido hielo y vapor) es 1 386 000 000 km³; no se usa ese valor porque $W$ representa agua líquida disponible.

La NASA reproduce la misma tabla en https://science.nasa.gov/earth/earth-observatory/the-water-cycle/ con total 1 385 984 ×10³ km³, coherente dentro del redondeo.

## [2] NASA/JPL — área superficial

Página: https://science.nasa.gov/moon/by-the-numbers/

La tabla «Earth and the Moon, by the numbers» da el área superficial terrestre como 510 064 472 km² (radio medio 6371,00 km). El valor es compatible con $4\pi R^2$ usando el radio volumen medio.

## [3] y [4] — humedad del suelo y umbral

- Shao, Fraedrich & Ishizuka (2021), https://doi.org/10.1007/s10546-020-00596-9 : en suelos hiperáridos el contenido residual de aire seco es $\theta_r \approx 0{,}02$ m³/m³ para arena (parámetros de Brooks–Corey) y $\theta_r \approx 0{,}045$ m³/m³ para van Genuchten. Se toma 0,02 como referencia de arena seca. El mismo artículo documenta valores medidos tan bajos como 0,005 m³/m³ en el desierto de Taklamakan, de modo que 0,02 no es un mínimo absoluto.
- Pan et al. (2019), https://doi.org/10.3390/w11050883 : con ESA CCI SM v4.2 (capa de 10 cm) las medias hemisféricas mensuales van de 0,05 a 0,25 m³/m³, y el mapa global va de 0 a 0,41 m³/m³. El Sahara y la península arábiga aparecen como las zonas más secas; trabajos derivados sitúan el Sahara por debajo de 0,05 m³/m³. Se adopta $\theta_{\rm ref}=0{,}15$ m³/m³ como valor medio representativo.
- Bekin et al. (2025), https://doi.org/10.1029/2025JG009193 : mediciones in situ en suelos del Sahara (sur de Marruecos) muestran incrementos nocturnos de 0,7–1 % en $\theta_v$ por adsorción de vapor, confirmando que el contenido base es extremadamente bajo y variable.

**Conversión a fracción:** $f_{\rm desert} = 0{,}02/0{,}15 = 0{,}1333$. Es una razón entre contenidos volumétricos (adimensional), no un cociente entre volúmenes globales. Se usa para convertir el umbral físico intuitivo «Sahara» en el parámetro adimensional que necesita la EDO. Una Tierra totalmente desértica tendría $W$ muy por debajo de este valor, por lo que el criterio es conservador.

## Supuestos ficticios

Herbert no proporciona cifras, de modo que $r, K, \mu, \delta, \alpha, \beta, \eta, \xi$ y $m_{\rm trucha}$ son supuestos declarados. El único anclaje narrativo es la maduración de ~4000 años, que se traduce en $\mu = 1/4000$ años⁻¹ como tasa media. El análisis de sensibilidad posterior tratará todos ellos como parámetros libres dentro de los rangos de la tabla.
