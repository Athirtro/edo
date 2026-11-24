# Derivaciones del análisis cualitativo (bloque 3)

Desarrollo algebraico que respalda las fórmulas del README.

## Equilibrios

Sistema reducido, con $W=W_{\rm tot}-\alpha T-\beta G$:

$$\dot T = g(W)T-\mu T+\xi\delta G,\qquad \dot G=\eta\mu T-\delta G.$$

De $\dot G=0$ se obtiene $G=(\eta\mu/\delta)T$. Sustituyendo en $\dot T=0$:

$$\bigl[g(W)-\mu+\xi\eta\mu\bigr]T = \bigl[g(W)-c\bigr]T = 0,\qquad c\equiv\mu(1-\xi\eta).$$

- $T=0$ da el equilibrio desértico $(0,0)$, con $W=W_{\rm tot}$.
- $T>0$ exige $g(W^*)=c$. Si $c\le 0$ no hay solución con $W>0$ (salvo el caso límite $c=0$, donde $W^*=0$). Si $c>0$ y $r>c$:

$$W^*=\frac{Kc}{r-c}.$$

La condición $W^*<W_{\rm tot}$ equivale a $r>c(1+K/W_{\rm tot})\equiv r_{\rm crit}$. Entonces

$$T^*=\frac{W_{\rm tot}-W^*}{\alpha+\beta\eta\mu/\delta},\qquad G^*=\frac{\eta\mu}{\delta}T^*.$$

El denominador surge de $\alpha T+\beta G=T(\alpha+\beta\eta\mu/\delta)=W_{\rm tot}-W^*$.

En $r=r_{\rm crit}$ se tiene $W^*=W_{\rm tot}$ y por tanto $T^*=G^*=0$: la rama de persistencia nace del equilibrio desértico. Es una transcrítica: cruza un autovalor por cero, lo que se detecta porque $\det J_{\rm desierto}=0$.

## Jacobiano

Con $g'(W)=rK/(K+W)^2$ y $\partial W/\partial T=-\alpha$, $\partial W/\partial G=-\beta$:

$$\frac{\partial \dot T}{\partial T}=g(W)-\mu-\alpha T g'(W),\qquad
\frac{\partial \dot T}{\partial G}=-\beta T g'(W)+\xi\delta,$$

$$\frac{\partial \dot G}{\partial T}=\eta\mu,\qquad
\frac{\partial \dot G}{\partial G}=-\delta.$$

En el desierto ($T=G=0$, $W=W_{\rm tot}$):

$$J_{\rm desierto}=\begin{pmatrix}g(W_{\rm tot})-\mu & \xi\delta\\ \eta\mu & -\delta\end{pmatrix},$$

$$\operatorname{tr}J=g(W_{\rm tot})-\mu-\delta,\qquad
\det J=-\delta\bigl[g(W_{\rm tot})-c\bigr].$$

Por tanto el desierto es silla si $g(W_{\rm tot})>c$, estable si $g(W_{\rm tot})<c$, y no hiperbólico en la igualdad. En $r=r_{\rm crit}$ la identidad es exacta: $g(W_{\rm tot})=c(1+K/W_{\rm tot})W_{\rm tot}/(K+W_{\rm tot})=c$.

## Verificación numérica

`test_jacobian_finite_differences` compara cada entrada del jacobiano analítico con diferencias centrales del sistema reducido en $(T,G)=(10^6,10^5)$. `test_persistence_equilibrium` comprueba los residuos en forma relativa porque los términos individuales son del orden de $10^{14}$ y se cancelan; una tolerancia absoluta no tendría sentido.

## Comparación de umbrales

- Umbral **dinámico**: $r=r_{\rm crit}$; responde a «¿existe un estado no desértico?».
- Umbral **convencional**: $W=f_{\rm desert}W_0$; responde a «¿cuándo se declara desierto?».

Con los parámetros base, $W^*\approx6{,}92\times10^6$ km³ mientras que $f_{\rm desert}W_0\approx1{,}815\times10^8$ km³. El sistema converge hacia un equilibrio con agua por debajo del umbral convencional, de modo que el cruce del umbral ocurre durante la transitoria.
