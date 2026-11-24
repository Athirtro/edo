# Decisiones de modelado (bloque 1)

Notas de contexto que complementan el README. No son necesarias para entender el modelo, pero explican por qué se eligió una formulación y no otra.

## 1. Liberación de truchas como tasa continua

La propuesta original dejaba abierta la liberación de truchas por muerte del gusano: «tasa promedio continua o evento discreto». Se eligió la **tasa continua** $\xi\delta G$ para mantener una EDO autónoma de tres variables y poder aplicar linealización y análisis de bifurcación con las herramientas estándar. El ciclo real es discreto (milenios de latencia y explosiones reproductivas); la tasa continua es su promedio. La variante discreta (modelo híbrido o sistema con retardo) queda como comprobación de robustez en bloques posteriores.

## 2. Factores de conversión separados

No se asumió que truchas y gusanos retuvieran la misma agua por tonelada ni que las conversiones conservaran biomasa. En su lugar se introdujeron $\alpha,\beta,\eta,\xi$ como parámetros libres. El balance de agua en cada transición se escribe de modo que el agua que deja una forma reaparezca en la otra o vuelva a $W$:

- maduración: sale $\alpha$ por tonelada de trucha, entra $\beta\eta$ por tonelada de gusano;
- mortalidad: sale $\beta$ por tonelada de gusano, entra $\alpha\xi$ por tonelada de trucha.

Esto evita crear o destruir agua implícitamente, que era el punto débil del borrador original.

## 3. Condiciones de invariancia como supuesto declarado

Las desigualdades $\alpha\ge\beta\eta$ y $\beta\ge\alpha\xi$ se imponen para que el campo vectorial no saque el sistema del dominio físico cuando $W=0$. No provienen del canon de *Dune*: son una convención matemática necesaria. Si en bloques posteriores se exploran parámetros que las violan, debe interpretarse como un modelo diferente (con recarga externa o agua fósil), no como una corrección del actual.

## 4. Reducir antes de linealizar

Con $W_{\rm tot}$ conservada, el jacobiano del sistema 3D tiene un autovalor cero que no aporta información sobre estabilidad. Por eso el análisis cualitativo se hará sobre el sistema 2D $(T,G)$; la prueba `test_reduction` garantiza que la reducción no altera las ecuaciones de $T$ y $G$.

## 5. Pruebas como especificación ejecutable

El bloque 1 se desarrolló con TDD. Las pruebas no comprueban valores numéricos de los parámetros (aún no están documentados), sino propiedades estructurales: conservación, no negatividad, saturación, equivalencia de la reducción y correcta conversión de unidades. Así, cuando el bloque 2 asigne magnitudes reales a los parámetros, las pruebas seguirán siendo válidas.
