# Tiempo de colapso (bloque 4)

## Detalle de la cuadratura

La aproximación cuasi-estacionaria sustituye $G(t)$ por su valor de equilibrio instantáneo $G=(\eta\mu/\delta)T$. Bajo esa hipótesis:

$$\alpha T+\beta G = T\left(\alpha+\frac{\beta\eta\mu}{\delta}\right) \equiv D T = W_{\rm tot}-W,$$

y por tanto $T=(W_{\rm tot}-W)/D$. Sustituyendo en $\dot W$:

$$\dot W = -\alpha g(W)T+(\alpha-\beta\eta)\mu T+(\beta-\alpha\xi)\delta\frac{\eta\mu}{\delta}T
       = -\Bigl[\alpha g(W)-\mu(\alpha-\beta\eta)-\eta\mu(\beta-\alpha\xi)\Bigr]\frac{W_{\rm tot}-W}{D}.$$

Definiendo $A(W)$ como el corchete dividido por $D$, se obtiene $\dot W=-A(W)(W_{\rm tot}-W)$.

### Por qué integrar en logaritmos

$W_{\rm tot}-W(0)=\alpha m_{\rm trucha}=5\times10^{-11}$ km³, mientras que $W_{\rm tot}-W_{\rm event}\approx1{,}18\times10^{9}$ km³: 20 órdenes de magnitud. Una cuadratura de Simpson uniforme en $W$ concentraría casi todos los nodos en la región incorrecta. Con el cambio $u=\ln(W_{\rm tot}-W)$:

$$dt = -\frac{dW}{A(W)(W_{\rm tot}-W)} = \frac{du}{A(W_{\rm tot}-e^u)},$$

y el integrando es suave (del orden de $1/A$, entre 96 y 373 años). Se usa Simpson con 2000 intervalos, que da $t_{\rm colap}^{\rm quad}=4366{,}0309$ años (valor con $N=2000$; el límite convergido es $4366{,}0300$).

### Escala numérica

Al sumar $W_0\approx1{,}36\times10^9$ y $\alpha m_{\rm trucha}=5\times10^{-11}$ en doble precisión, el segundo término se pierde por completo ($\epsilon_{\rm máq}W_0\approx3\times10^{-7}$). Por eso el código calcula $S_0=\alpha m_{\rm trucha}$ directamente en lugar de restar $W_{\rm tot}-W_0$. Es un detalle numérico, no una corrección del modelo.

## Integrador RK4 adaptativo

Se integra $S=W_{\rm tot}-W$, con $\dot S=-\dot W$ y $\dot T,\dot G$ dados por el modelo. Integrar $W$ directamente sería equivalente en precisión infinita, pero en coma flotante la tolerancia se mediría contra $W\sim10^9$ e ignoraría la variable pequeña que está creciendo.

- En cada paso se calculan un paso completo $y^{(1)}$ y dos medios pasos $y^{(2)}$. El estimador es $|y^{(1)}-y^{(2)}|/(10^{-30}+{\rm tol}\max(|y^{(1)}|,|y^{(2)}|))/15$; el factor $15$ viene de Richardson para RK4.
- Si el error normalizado es $\le1$ se acepta $y^{(2)}+ (y^{(2)}-y^{(1)})/15$ (extrapolación de Richardson); si no, se reduce el paso a la mitad. El factor de crecimiento tras aceptar es $0{,}9\,{\rm err}^{-0{,}2}$, acotado entre 0,2 y 5, con paso máximo 50 años.
- El evento se detecta cuando $S$ cruza $S_{\rm event}=W_{\rm tot}-W_{\rm event}$. Se localiza por bisección de 50 iteraciones usando **el mismo operador extrapolado** que se acepta en el paso; se devuelve el estado en el extremo superior del intervalo de bisección, garantizando que el cruce ya ocurrió. Usar un RK4 simple en la bisección sesgaría $t_{\rm colap}$ y podría dejar el estado final por encima del umbral.
- La conservación se verifica con $W+\alpha T+\beta G = W_{\rm tot}$, no con la identidad algebraica $S=\alpha T+\beta G$. RK4 preserva esa combinación lineal salvo redondeo.

### Convergencia

| `tol` | pasos | $t_{\rm colap}$ (años) |
| ---: | ---: | ---: |
| $10^{-3}$ | 83 | 4024,3132 |
| $10^{-6}$ | 172 | 4024,2793 |
| $10^{-7}$ | 274 | 4024,2783 |
| $10^{-9}$ | 694 | 4024,2782 |
| $10^{-10}$ | 1104 | 4024,2782 |

La referencia independiente, RK4 de paso fijo $h=0{,}001$ sobre $(W,T,G)$, da $4024{,}278215360$ años y se reproduce con `make reference` (`artifacts/reference_rk4.c`). La versión adaptativa se acerca a ese valor de forma monótona al reducir la tolerancia; de $10^{-6}$ a $10^{-7}$ el cambio es de $9{,}4\times10^{-4}$ años.

## Resultados

| Método | $t_{\rm colap}$ (años) |
| --- | ---: |
| Cuadratura cuasi-estacionaria | 4366,0309 |
| RK4 adaptativo | 4024,2793 |
| RK4 paso fijo (referencia) | 4024,2782 |

Diferencia relativa cuadratura/RK4: 8,49 %. Error máximo de conservación: $\sim5\times10^{-16}$ relativo. La conjetura de 8000 años no se reproduce con los parámetros base; el modelo da aproximadamente la mitad.
