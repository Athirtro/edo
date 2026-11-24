# Exposición (8–10 min)

## Guion y figuras

1. **Pregunta (30 s).** Una trucha, agua global y un ciclo trucha→gusano→trucha. ¿Cuándo se declara desierto y qué controla ese tiempo?
2. **Modelo (1 min).** $\dot W,\dot T,\dot G$ con $W_{\rm tot}=W+\alpha T+\beta G$. Reducimos a $(T,G)$ para linealizar sin autovalores espurios.
3. **Parámetros (1 min).** Dos datos con fuente, un umbral convencional y nueve supuestos. Enseñar la tabla del README; remarcar que el resultado es condicional.
4. **Análisis (1,5 min).** Desierto silla, persistencia estable, $r_{\rm crit}\approx2{,}38\times10^{-4}$ años⁻¹. El umbral dinámico y el convencional no son lo mismo.
5. **Colapso (2 min).** Usar `figures/series.svg`: $W$ cae, $T$ y $G$ crecen, evento en $t_{\rm colap}\approx4024$ años. Comparar 4366 (cuadratura) y 4024,2782 (referencia).
6. **Retrato de fase (1,5 min).** Usar `figures/retrato_fase.svg`: trayectoria desde $(m_{\rm trucha},0)$ hacia persistencia, con ejes escalados por `max_T` y `max_G` (se indican las cifras en el propio SVG). La línea morada es la separatriz operativa $W=W_{\rm event}$; aclarar que la variedad estable del desierto sale del cuadrante físico, así que no hay frontera de cuenca interna.
7. **Sensibilidad (1,5 min).** Usar `figures/heatmap.svg` (cuadratura cuasi-estacionaria): aumentar $r$ acorta el colapso; aumentar $K$ lo alarga. Morris confirma $r>K>\mu$; $\alpha,\beta,m_{\rm trucha}$ no son identificables por separado.
8. **Cierre (30 s).** Respuesta: ~4024 años **bajo estos supuestos**, no 8000. Mensaje condicional: sin escala de biomasa y sin calibración, el número es una demostración del método, no una predicción.

## Notas para preguntas

- «¿Por qué no 8000?». El parámetro $r=0{,}02$ años⁻¹ y $K=10^9$ km³ dan ~4024; duplicar $r$ lo reduce aproximadamente a la mitad ($E_r\approx-1$).
- «¿Es el Sahara un volumen?». No: es una fracción $0{,}02/0{,}15$ aplicada a $W_0$.
- «¿Hay colapso irreversible?». El desierto es inestable ante una trucha; lo irreversible aquí es el cruce del convenio de aridificación, no una transición de cuenca.
- «¿Cómo se generaron las figuras?». Con `make figures`: C puro escribe SVG, sin dependencias externas.

Todas las figuras citadas están en `figures/` y se regeneran con `make figures`.
