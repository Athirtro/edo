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

## Cómo ejecutar las pruebas

```bash
make test
```

Requiere `gcc` y `make`. El arnés de pruebas está en `artifacts/test_model.c`; no usa bibliotecas externas. El flujo fue **TDD**: primero se escribieron las pruebas y se comprobó que fallaban por ausencia de `src/model.c` (fase roja), y después se implementó el modelo hasta hacerlas pasar (fase verde).

| Prueba | Propiedad que especifica |
| --- | --- |
| `test_units` | 1 km³ = 10⁹ t de agua |
| `test_params_valid` | positividad y condiciones $\alpha\ge\beta\eta$, $\beta\ge\alpha\xi$ |
| `test_growth` | $g(0)=0$, $g(K)=r/2$, monótona y saturada en $r$ |
| `test_initial_state` | una trucha inicial y agua total correcta |
| `test_water_conservation` | $\dot W_{\rm tot}=0$ |
| `test_invariance` | no negatividad y cota superior de $W$ |
| `test_reduction` | equivalencia del sistema 2D y el 3D |

## Estado y siguientes bloques

- [x] Bloque 1 — Formulación del modelo (unidades, EDO, condiciones iniciales, invariancia y reducción).
- [ ] Bloque 2 — Búsqueda y citación de datos reales y construcción de la tabla de parámetros.
- [ ] Bloque 3 — Equilibrios, linealización y umbrales de bifurcación.
- [ ] Bloques restantes — Tiempo de colapso, sensibilidad, visualización e informe, según `plan-trucha-arena.md`.
