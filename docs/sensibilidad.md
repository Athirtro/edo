# Sensibilidad (bloque 5)

## Conjunto analizado

Se consideran ocho de los nueve parámetros con `PARAM_ASSUMPTION` en `src/parameters.c`: $r,K,\mu,\delta,\alpha,\beta,\eta,\xi$. El noveno, $m_{\rm trucha}$, se mantiene fijo como escala de biomasa por la no identificabilidad demostrada más abajo. Se excluyen además $W_0$ y $A$ (datos) y $f_{\rm desert}$ (convención que define el evento, no la dinámica). Los rangos se leen directamente de `PARAMETERS`, de modo que README, pruebas y simulaciones comparten una única tabla.

Todas las simulaciones usan `SENS_TOL = 1e-6` y `SENS_T_MAX = 20000` años. Si no se alcanza el umbral, se censura $t_{\rm colap}$ a $SENS_T_MAX$ para que los efectos elementales sean finitos.

## Sensibilidad local

Para el parámetro $p$, con $h$ pequeño y contenido en su rango:

$$\frac{\partial t_{\rm colap}}{\partial p}
\approx \frac{t_{\rm colap}(p+h)-t_{\rm colap}(p-h)}{2h},
\qquad
E_p=\frac{p}{t_{\rm colap}}\frac{\partial t_{\rm colap}}{\partial p}.$$

La elasticidad $E_p$ es adimensional: $E_r=-1{,}020$ significa que un 1 % más de $r$ reduce $t_{\rm colap}$ aproximadamente un 1,02 %. Las derivadas de la tabla se calculan con $h=0{,}01\,p$. La prueba usa ese mismo $h$ para `local_derivative` y lo contrasta con una diferencia externa de paso $h/2=0{,}005\,p$, de modo que la comprobación no repite exactamente el cálculo interno.

## Cribado global de Morris

Se trabaja en coordenadas normalizadas $x_p\in[0,1]$, con $p=p_{\min}+x_p(p_{\max}-p_{\min})$. En cada trayectoria:

1. Se elige al azar un signo $s_p\in\{+1,-1\}$ y un punto inicial $x_p\in[0,3/4]$ si $s_p=+1$, o $x_p\in[1/4,1]$ si $s_p=-1$. Así el sentido opuesto suele existir dentro del dominio; si cae fuera, ese efecto se descarta.
2. Se baraja el orden de los ocho parámetros.
3. Para cada parámetro $p$ en ese orden, se evalúa el efecto elemental

$${\rm EE}_p=\frac{t_{\rm colap}(x_p+s_p\Delta)-t_{\rm colap}(x_p)}{s_p\Delta},
\qquad \Delta=\tfrac14,$$

manteniendo los demás parámetros en los extremos ya fijados por la trayectoria. Si el nuevo punto viola $\beta\ge\alpha\xi$ o $\alpha\ge\beta\eta$, se prueba el sentido opuesto; si ambos son inválidos, se descarta ese efecto.

Se usan 100 trayectorias y semilla fija `12345` (LCG `1664525`, `1013904223`), para que el resultado sea reproducible sin dependencias externas. Se resumen

$$\mu^*_p=\frac{1}{n_p}\sum |{\rm EE}_p|,
\qquad
\sigma_p=\sqrt{\frac{1}{n_p}\sum({\rm EE}_p-\overline{\rm EE}_p)^2}.$$

$\mu^*$ mide influencia absoluta; $\sigma$ indica no linealidad o interacciones. $n_p$ se publica junto a los resultados porque las restricciones descartan efectos y los puestos bajos del ranking dependen de ese submuestreo. La función `morris_screen` es determinista y se ejecuta desde las pruebas, sin artefactos adicionales.

## No identificabilidad inducida por la conservación

Sea $\lambda\in[0{,}25,\,5]$ (intersección de los rangos admisibles de $m$, $\alpha$ y $\beta$) y la transformación

$$\alpha'=\lambda\alpha,\quad \beta'=\lambda\beta,\quad
m'=\frac{m}{\lambda},\quad T'=\frac{T}{\lambda},\quad G'=\frac{G}{\lambda}.$$

El agua total se conserva porque $\alpha'T'+\beta'G'=\alpha T+\beta G$, luego $W'=W$. Las ecuaciones de biomasa escalan de forma coherente:

$$\dot T'=g(W)T'-\mu T'+\xi\delta G',\qquad
\dot G'=\eta\mu T'-\delta G',$$

ya que $\dot T'=\dot T/\lambda$ y $\dot G'=\dot G/\lambda$. Por tanto $W(t)$ y el instante del evento son exactamente los mismos. Ésta es una no identificabilidad estructural: observando sólo $W(t)$ o $t_{\rm colap}$ no se pueden separar $\alpha$, $\beta$ y $m$; hace falta una escala absoluta de biomasa. Como $m=0{,}05$ t es un supuesto del modelo, la cifra de 4024 años es condicional a esa elección.

## Resultados

Con 100 trayectorias y $\Delta=1/4$, los efectos $\mu^*$ (años) son: $r=8559$ ($n=100$), $K=1202$ ($100$), $\mu=172$ ($100$), $\alpha=58{,}5$ ($89$), $\xi=51{,}6$ ($99$), $\eta=15{,}4$ ($97$), $\delta=5{,}34$ ($100$) y $\beta=2{,}65$ ($84$).

El ranking global es $r > K > \mu > \alpha > \xi > \eta > \delta > \beta$. Los efectos elementales de $r$ tienen $\sigma$ comparable a $\mu^*$, señal de que el sistema responde de forma muy no lineal dentro del rango de $r$: en algunos puntos $r$ está por debajo de $r_{\rm crit}$ y no hay colapso (caso censurado), mientras que en otros el colapso es mucho más rápido. El cribado no sustituye a una descomposición de la varianza tipo Sobol; identifica qué parámetros merecen un estudio más detallado.
