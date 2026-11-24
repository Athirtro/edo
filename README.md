# ¿Cuánto tardaría una trucha de arena en desertificar la Tierra?

Modelo conceptual de ecuaciones diferenciales ordinarias (EDO) inspirado en el ciclo de agua, truchas y gusanos de *Dune*. **No es una réplica del canon ni una predicción física**: es un laboratorio para estudiar cómo un sistema con un recurso limitado puede transitar hacia el colapso. La pregunta que guía el proyecto es:

> ¿Cómo depende el tiempo de desertificación de los parámetros biológicos y ambientales de un sistema agua-trucha-gusano cuando se introduce una única trucha, y existen umbrales a partir de los cuales el sistema colapsa de forma irreversible?

Este documento es el entregable principal del proyecto. El código fuente está en `src/`, las pruebas en `artifacts/` y las notas de diseño en `docs/`.

## Bloque 1 — Formulación del modelo

### Variables de estado

| Símbolo | Significado | Unidad |
| --- | --- | --- |
| $W(t)$ | agua líquida disponible (superficie y subsuelo) | km³ |
| $T(t)$ | biomasa de truchas de arena | toneladas |
| $G(t)$ | biomasa de gusanos adultos | toneladas |
| $t$ | tiempo | años |

### Unidades y conversiones

Un km³ de agua equivale a $10^9$ toneladas, por lo que $1\ \mathrm{t} = 10^{-9}\ \mathrm{km^3}$. Las funciones `ton_to_km3` y `km3_to_ton` de `src/model.c` encapsulan esta conversión.

Además de las variables de estado se definen cuatro factores que cierran el balance:

| Símbolo | Significado | Unidad |
| --- | --- | --- |
| $\alpha$ | agua retenida por tonelada de trucha (biomasa + reservas asociadas) | km³/t |
| $\beta$ | agua retenida por tonelada de gusano | km³/t |
| $\eta$ | rendimiento de gusano por tonelada de trucha que madura | t/t |
| $\xi$ | rendimiento de trucha por tonelada de gusano que muere | t/t |

La masa de una trucha, $m_{\rm trucha}$ (toneladas), convierte la condición inicial «una trucha» en biomasa.

### Ecuaciones

La tasa de crecimiento de las truchas es saturante, de tipo Michaelis-Menten:

$$g(W) = \frac{rW}{K+W},$$

con $r>0$ la tasa máxima (años⁻¹) y $K>0$ el agua de media saturación (km³). El sistema completo es:

$$
\begin{aligned}
\dot W &= -\alpha\, g(W)\, T
       + (\alpha - \beta\eta)\,\mu T
       + (\beta - \alpha\xi)\,\delta G,\\[4pt]
\dot T &= g(W)\,T - \mu T + \xi\,\delta G,\\[4pt]
\dot G &= \eta\,\mu T - \delta G.
\end{aligned}
$$

Los parámetros $\mu$ y $\delta$ son las tasas de maduración y mortalidad (años⁻¹). Interpretación de los términos:

- $-\alpha g(W)T$: el agua que consumen las truchas al crecer se descuenta de $W$ y pasa a la biomasa y a las reservas subterráneas asociadas.
- $(\alpha-\beta\eta)\mu T$: al madurar, el agua retenida por la biomasa que sale de $T$ ($\alpha$ por tonelada) se reorganiza en la biomasa que entra en $G$ ($\beta\eta$ por tonelada). La diferencia regresa a $W$ o se extrae de él.
- $(\beta-\alpha\xi)\delta G$: lo análogo cuando un gusano muere y libera truchas.
- $-\mu T$ y $+\eta\mu T$: transferencia continua trucha → gusano.
- $+\xi\delta G$ y $-\delta G$: transferencia continua gusano → trucha.

La liberación de truchas se modela como una **tasa promedio continua**, no como un evento discreto; así el sistema sigue siendo una EDO autónoma. La variante discreta queda reservada para una validación posterior (véase `docs/decisiones.md`).

### Conservación del agua

La cantidad conservada no es $W$ por sí solo, sino el agua líquida más el agua retenida en la biomasa:

$$W_{\rm tot} = W + \alpha T + \beta G.$$

En efecto, al derivar y sustituir las tres EDO los términos se cancelan:

$$
\begin{aligned}
\frac{dW_{\rm tot}}{dt}
&= \dot W + \alpha\dot T + \beta\dot G\\
&= \bigl[-\alpha gT + (\alpha-\beta\eta)\mu T + (\beta-\alpha\xi)\delta G\bigr]
   + \alpha\bigl[gT-\mu T+\xi\delta G\bigr]
   + \beta\bigl[\eta\mu T-\delta G\bigr] = 0.
\end{aligned}
$$

Esta identidad está verificada por la prueba `test_water_conservation`.

### Condiciones iniciales

$$W(0) = W_0,\qquad T(0) = m_{\rm trucha},\qquad G(0) = 0.$$

Por tanto, el agua total del sistema es $W_{\rm tot} = W_0 + \alpha\,m_{\rm trucha}$. Los valores numéricos de $W_0$ y $m_{\rm trucha}$ se documentarán en el siguiente bloque, junto con el resto de parámetros; las pruebas actuales usan valores arbitrarios porque sólo verifican propiedades estructurales.

### Dominio invariante y restricciones

Se exigen todos los parámetros estrictamente positivos y dos condiciones adicionales:

$$\alpha \ge \beta\eta,\qquad \beta \ge \alpha\xi,$$

que implican $\eta\xi \le 1$. Su interpretación: ninguna transición entre truchas y gusanos puede requerir más agua de la que libera. Con estas condiciones el dominio

$$\mathcal{D} = \{(W,T,G) : W\ge 0,\ T\ge 0,\ G\ge 0,\ W+\alpha T+\beta G = W_{\rm tot}\}$$

es positivamente invariante:

- Si $W=0$, entonces $g(0)=0$ y $\dot W = (\alpha-\beta\eta)\mu T + (\beta-\alpha\xi)\delta G \ge 0$.
- Si $T=0$, entonces $\dot T = \xi\delta G \ge 0$.
- Si $G=0$, entonces $\dot G = \eta\mu T \ge 0$.
- Como $\alpha,\beta>0$ y $T,G\ge 0$, se tiene $W = W_{\rm tot}-\alpha T-\beta G \le W_{\rm tot}$.

Estas cuatro afirmaciones están cubiertas por `test_invariance` y por la validación `params_valid`.

### Reducción a dos dimensiones

Como $W_{\rm tot}$ es constante, se puede eliminar $W$ usando $W = W_{\rm tot}-\alpha T-\beta G$ y estudiar el sistema en el plano $(T,G)$:

$$
\begin{aligned}
\dot T &= g\!\left(W_{\rm tot}-\alpha T-\beta G\right)T - \mu T + \xi\delta G,\\
\dot G &= \eta\mu T - \delta G,
\end{aligned}
$$

con restricción $\alpha T+\beta G \le W_{\rm tot}$, $T\ge 0$, $G\ge 0$. Este sistema reducido es el que se linealizará en el siguiente bloque: hacerlo directamente sobre las tres variables originales daría un autovalor cero espurio debido a la conservación. La prueba `test_reduction` comprueba que `rhs_reduced` coincide con `rhs` tras reconstruir $W$.

## Bloque 2 — Datos y parámetros

### Volumen de agua líquida y área superficial

La fuente es la compilación de la USGS (tabla de Shiklomanov/Gleick, 1993), citada en [1]. Se suman únicamente los reservorios **líquidos**:

| Reservorio | Volumen (km³) |
| --- | ---: |
| Océanos, mares y bahías | 1 338 000 000 |
| Aguas subterráneas (dulces y salinas) | 23 400 000 |
| Lagos dulces | 91 000 |
| Lagos salinos | 85 400 |
| Humedales | 11 470 |
| Humedad del suelo | 16 500 |
| Ríos | 2 120 |
| Agua biológica | 1 120 |
| **Total, $W_0$** | **1 361 607 610** |

Se excluyen casquetes y glaciares (24 064 000 km³), permafrost (300 000 km³) y vapor atmosférico (12 900 km³) por no ser agua líquida disponible. El área superficial de la Tierra es $A = 5{,}10064472\times10^8$ km², según NASA/JPL [2]. Con estos dos datos, la capa equivalente de agua líquida es $W_0/A \approx 2{,}67$ km.

### Umbral de desierto

La humedad volumétrica del suelo en el Sahara es del orden de $\theta_{\rm Sahara} \approx 0{,}02$ m³/m³: es el valor residual de aire seco para arena usado en modelos de superficie [3], y los productos satelitales sitúan el Sahara por debajo de 0,05 m³/m³ [4]. Como referencia húmeda se toma $\theta_{\rm ref} \approx 0{,}15$ m³/m³, valor representativo del rango hemisférico 0,05–0,25 m³/m³ observado por ESA CCI [4].

La convención adoptada es:

$$f_{\rm desert} = \frac{\theta_{\rm Sahara}}{\theta_{\rm ref}} = \frac{0{,}02}{0{,}15} \approx 0{,}133.$$

Es decir, se declara el desierto cuando $W < f_{\rm desert} W_0 \approx 1{,}815\times10^8$ km³: el agua restante, *en relación con la inicial*, es comparable a la humedad del Sahara *en relación con un suelo húmedo de referencia*. **No es un cálculo de volumen global literal**: una Tierra totalmente seca sería más árida que el Sahara, por lo que este umbral es una convención deliberadamente anticipada. Su rango exploratorio es $[0{,}05,\,0{,}25]$.

### Tabla de parámetros

Los símbolos, unidades, valores, rangos y orígenes están también codificados en `src/parameters.c`, de modo que el README y las pruebas comparten una única fuente de verdad.

| Símbolo | Unidad | Valor | Rango | Origen |
| --- | --- | ---: | --- | --- |
| $W_0$ | km³ | $1{,}36160761\times10^9$ | $[1{,}33\times10^9,\,1{,}39\times10^9]$ | Dato [1] |
| $A$ | km² | $5{,}10064472\times10^8$ | $[5{,}099\times10^8,\,5{,}101\times10^8]$ | Dato [2] |
| $f_{\rm desert}$ | adim. | $0{,}1333$ | $[0{,}05,\,0{,}25]$ | Convención [3, 4] |
| $r$ | 1/año | $0{,}02$ | $[0{,}001,\,0{,}2]$ | Supuesto |
| $K$ | km³ | $1{,}0\times10^9$ | $[1{,}0\times10^8,\,1{,}4\times10^9]$ | Supuesto |
| $\mu$ | 1/año | $2{,}5\times10^{-4}$ | $[1\times10^{-4},\,1\times10^{-3}]$ | Supuesto (~4000 años) |
| $\delta$ | 1/año | $2{,}0\times10^{-3}$ | $[2\times10^{-4},\,1\times10^{-2}]$ | Supuesto |
| $\alpha$ | km³/t | $1{,}0\times10^{-9}$ | $[1\times10^{-10},\,1\times10^{-8}]$ | Supuesto |
| $\beta$ | km³/t | $5{,}0\times10^{-10}$ | $[5\times10^{-11},\,5\times10^{-9}]$ | Supuesto |
| $\eta$ | t/t | $1{,}5$ | $[0{,}5,\,2{,}0]$ | Supuesto |
| $\xi$ | t/t | $0{,}3$ | $[0{,}05,\,0{,}5]$ | Supuesto |
| $m_{\rm trucha}$ | t | $0{,}05$ | $[0{,}01,\,0{,}2]$ | Supuesto |

Los rangos de $\eta$ y $\xi$ respetan las restricciones del bloque 1, $\eta \le \alpha/\beta = 2$ y $\xi \le \beta/\alpha = 0{,}5$. Todos los supuestos son parámetros libres: el análisis posterior estudiará cómo cambia el resultado al variarlos, no buscará un único valor «correcto».

**Referencias**

1. U.S. Geological Survey, *How Much Water is There on Earth?* — https://www.usgs.gov/water-science-school/science/how-much-water-there-earth
2. NASA/JPL, *Compare Earth and the Moon* — https://science.nasa.gov/moon/by-the-numbers/
3. Shao, Y., Fraedrich, K. & Ishizuka, M. (2021), *Modelling Soil Moisture in Hyper-Arid Conditions*, Boundary-Layer Meteorology — https://doi.org/10.1007/s10546-020-00596-9
4. Pan, N. et al. (2019), *Global Surface Soil Moisture Dynamics in 1979–2016 Observed from ESA CCI SM Dataset*, Water 11(5):883 — https://doi.org/10.3390/w11050883

La extracción detallada de las cifras está en `docs/fuentes.md`.

## Bloque 3 — Análisis cualitativo

Trabajamos sobre el sistema reducido $(T,G)$, con $W = W_{\rm tot}-\alpha T-\beta G$. Escribiendo $c=\mu(1-\xi\eta)$, los equilibrios son:

- **Desierto:** $(T,G)=(0,0)$, con $W=W_{\rm tot}$. Existe siempre.
- **Persistencia:** si $c>0$ y $r>c$, se requiere $g(W^*)=c$, de donde

$$W^* = \frac{Kc}{r-c},\qquad
T^* = \frac{W_{\rm tot}-W^*}{\alpha+\beta\eta\mu/\delta},\qquad
G^* = \frac{\eta\mu}{\delta}T^*.$$

La rama existe sólo si $W^*<W_{\rm tot}$. El valor crítico en que colisiona con el desierto es

$$r_{\rm crit} = c\left(1+\frac{K}{W_{\rm tot}}\right).$$

Es una **bifurcación de transcrítica**: para $r<r_{\rm crit}$ el desierto es el único equilibrio; para $r>r_{\rm crit}$ aparece el estado de persistencia.

El jacobiano reducido, con $g'(W)=rK/(K+W)^2$, es

$$J(T,G) = \begin{pmatrix}
g(W)-\mu-\alpha T g'(W) & -\beta T g'(W)+\xi\delta\\
\eta\mu & -\delta
\end{pmatrix}.$$

En el desierto, $\det J = -\delta\,[g(W_{\rm tot})-c]$: con los parámetros base $g(W_{\rm tot})\approx0{,}58\,r \gg c$, luego el desierto es **silla** (invasible por una sola trucha). En la persistencia, traza negativa y determinante positivo: es un **nodo estable**.

Con los valores del bloque 2: $c=1{,}375\times10^{-4}$ años⁻¹, $r_{\rm crit}\approx2{,}3848\times10^{-4}$ años⁻¹, $W^*\approx6{,}92\times10^{6}$ km³ y $f_{\rm desert}W_0\approx1{,}815\times10^{8}$ km³.

**Comparación de umbrales.** El umbral dinámico ($r=r_{\rm crit}$) y el convencional ($W=f_{\rm desert}W_0$) responden preguntas distintas. El primero decide *si existe* un estado no desértico; el segundo decide *cuándo se cruza* una convención de aridificación. En el conjunto base $W^* < f_{\rm desert}W_0$: el sistema se dirige hacia un equilibrio de persistencia, pero éste está por debajo del umbral convencional, así que el tiempo de cruce se alcanzará antes de estabilizarse. La prueba `test_threshold_comparison` verifica esta relación.

## Cómo ejecutar las pruebas

```bash
make test
```

Requiere `gcc` y `make`. Los arneses están en `artifacts/test_model.c`, `artifacts/test_parameters.c` y `artifacts/test_equilibria.c`; no usan bibliotecas externas. El flujo fue **TDD** en todos los bloques: primero se escribieron las pruebas y se comprobó que fallaban por ausencia del módulo correspondiente (fase roja), y después se implementó hasta hacerlas pasar (fase verde).

| Prueba | Propiedad que especifica |
| --- | --- |
| `test_units` | 1 km³ = 10⁹ t de agua |
| `test_params_valid` | positividad y condiciones $\alpha\ge\beta\eta$, $\beta\ge\alpha\xi$ |
| `test_growth` | $g(0)=0$, $g(K)=r/2$, monótona y saturada en $r$ |
| `test_initial_state` | una trucha inicial y agua total correcta |
| `test_water_conservation` | $\dot W_{\rm tot}=0$ |
| `test_invariance` | no negatividad y cota superior de $W$ |
| `test_reduction` | equivalencia del sistema 2D y el 3D |
| `test_metadata_complete` | metadatos y rangos válidos en toda la tabla |
| `test_required_symbols` | presencia de los 12 parámetros requeridos |
| `test_data_values` | $W_0$ y $A$ coinciden con las fuentes citadas |
| `test_threshold_convention` | $f_{\rm desert}=\theta_{\rm Sahara}/\theta_{\rm ref}$ |
| `test_assumptions_registered` | supuestos marcados como tales; $\mu=1/4000$ |
| `test_defaults` | `params_default` alimenta al modelo y cumple la invariancia |
| `test_desert_equilibrium` | el desierto existe con $W=W_{\rm tot}$ y residuos nulos |
| `test_persistence_equilibrium` | $g(W^*)=c$, $G^*=\eta\mu T^*/\delta$, conservación |
| `test_bifurcation_threshold` | $r_{\rm crit}$ y existencia de la rama a cada lado |
| `test_jacobian_finite_differences` | jacobiano analítico vs. diferencias finitas |
| `test_stability_classes` | desierto silla, persistencia estable |
| `test_threshold_comparison` | umbral dinámico vs. convencional |

## Estado y siguientes bloques

- [x] Bloque 1 — Formulación del modelo (unidades, EDO, condiciones iniciales, invariancia y reducción).
- [x] Bloque 2 — Datos y parámetros (fuentes citadas, umbral convencional y tabla de parámetros).
- [x] Bloque 3 — Equilibrios, linealización y umbrales de bifurcación.
- [ ] Bloque 4 — Tiempo de colapso.
- [ ] Bloque 5 — Sensibilidad.
- [ ] Bloque 6 — Visualización y entrega.
