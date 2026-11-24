# Informe final — Desertificación por una trucha de arena

## Pregunta y respuesta condicional

¿Cuánto tardaría una trucha de arena en desertificar la Tierra bajo el modelo agua–trucha–gusano propuesto? Con los parámetros base (todos los biológicos son supuestos explícitos), el evento convencional $W=f_{\rm desert}W_0$ se alcanza en $t_{\rm colap}\approx 4024$ años. La cuadratura cuasi-estacionaria da 4366 años y la referencia independiente con RK4 de paso fijo $h=0{,}001$, 4024,2782 años. Por tanto, la conjetura de 8000 años queda descartada **para este conjunto de supuestos**, pero no en sentido absoluto: $r$ y $K$ dominan la respuesta y sus rangos son epistémicos.

## Modelo

Se conserva el agua total $W_{\rm tot}=W+\alpha T+\beta G$. La dinámica reducida en $(T,G)$ evita el autovalor cero espurio del sistema 3D. Los parámetros de datos son $W_0$ y $A$; el umbral $f_{\rm desert}=0{,}02/0{,}15$ es una convención; los nueve restantes son supuestos. El dominio físico es invariante si $\alpha\ge\beta\eta$ y $\beta\ge\alpha\xi$.

## Resultados por bloque

1. **Formulación:** conversiones, balance de masa, invariancia y reducción a 2D verificadas.
2. **Datos y parámetros:** fuentes citadas, tabla única en `src/parameters.c`.
3. **Análisis cualitativo:** desierto silla, persistencia nodo estable, $r_{\rm crit}\approx2{,}38\times10^{-4}$ años⁻¹. El umbral dinámico y el convencional responden preguntas distintas.
4. **Tiempo de colapso:** evento por bisección con el operador RK4 aceptado; $\max|\Delta W_{\rm tot}|/W_{\rm tot}\sim5\times10^{-16}$.
5. **Sensibilidad:** elasticidades $E_r=-1{,}020$, $E_K=+0{,}443$; Morris (semilla 12345) rankea $r>K>\mu>\alpha>\xi>\eta>\delta>\beta$. La simetría $\alpha,\beta,m_{\rm trucha}$ hace que $m_{\rm trucha}$ no sea identificable por separado desde $t_{\rm colap}$.
6. **Visualización:** figuras SVG generadas desde C puro, sin dependencias externas. `series.svg` marca el evento; `retrato_fase.svg` muestra trayectoria, equilibrios y la separatriz operativa; `heatmap.svg` explora $r$ y $K$ mediante la cuadratura cuasi-estacionaria (4366 años en el conjunto base, frente a 4024,2793 de RK4). La tabla completa de parámetros y las pruebas están en el README.

## Separatriz: precisión importante

El desierto es un punto de silla, pero su autovector estable tiene $G<0$ y por tanto la variedad estable sale del cuadrante físico: **no existe una frontera de cuenca interna** que separe cuencas dentro del dominio permitido. En el retrato se dibuja como “separatriz operativa” la recta $W=W_{\rm event}$, equivalente a $\alpha T+\beta G=W_{\rm tot}-W_{\rm event}$. Esa recta separa el estado declarado desierto por convención del resto del dominio; no debe interpretarse como frontera de cuenca dinámica.

## Limitaciones

Esta sección resume las limitaciones del modelo y del análisis; ninguna de ellas se oculta en las cifras centrales.

- No es una simulación del canon de *Dune*: las tasas biológicas y los rendimientos son supuestos sin calibración empírica.
- El umbral 0,133 es convencional: una Tierra totalmente seca sería más árida que el Sahara.
- La liberación de truchas es continua; un evento discreto de aparición de gusanos cambiaría la dinámica.
- La cuadratura supone $G\approx(\eta\mu/\delta)T$; por eso se contrasta contra RK4.
- Morris rankea pero no descompone varianza; no hay inferencia bayesiana ni ajuste a datos.
- $\alpha$, $\beta$ y $m_{\rm trucha}$ son no identificables sin una escala externa de biomasa.
- El mapa de calor usa un zoom $r\in[0{,}002,0{,}05]$ donde el evento existe para todo $K$ de la tabla; fuera de esa zona el modelo puede no cruzar el umbral.

## Reproducción

```bash
make test       # incluye generacion de figuras y todas las pruebas
make reference  # referencia RK4 de paso fijo
make figures    # solo regenera figures/*.svg
```

Los SVG son texto legible y se pueden inspeccionar o versionar directamente. `make clean` borra binarios, no las figuras. No se necesita Python, gnuplot ni matplotlib.
