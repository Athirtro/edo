# Propuesta de proyecto EDO: ¿Cuánto tardaría una trucha de arena en desertificar la Tierra?

## Pregunta de investigación

**¿Cómo depende el tiempo de desertificación de la Tierra de los parámetros biológicos y ambientales de un sistema agua-trucha-gusano inspirado en Dune, cuando se introduce una única trucha de arena, y existen umbrales a partir de los cuales el sistema colapsa de forma irreversible?**

## Motivación

El ciclo de vida del gusano de arena en la saga Dune es un mecanismo ecológico con retroalimentación: las truchas de arena retienen el agua del planeta y, al madurar, dan lugar a los gusanos. Este proyecto toma ese mecanismo como un laboratorio conceptual para estudiar, con ecuaciones diferenciales ordinarias, cómo un sistema con un recurso limitado puede transitar hacia un colapso. Es un modelo conceptual inspirado en la novela, no una réplica del canon ni una predicción física.

## Objetivo general

Modelar con EDO la dinámica del agua líquida disponible, la biomasa de truchas y la biomasa de gusanos tras la introducción de una única trucha, y analizar de qué parámetros depende el tiempo hasta alcanzar el estado desértico.

## Objetivos específicos

- Formular un sistema de tres EDO acopladas con conservación de la masa de agua.
- Identificar los puntos de equilibrio y analizar su estabilidad mediante linealización (jacobiano).
- Definir un umbral de desierto y estimar el tiempo de cruce para distintos valores de los parámetros.
- Realizar un análisis de sensibilidad para determinar qué parámetros dominan el tiempo de colapso.
- (Extra) Simular numéricamente el sistema y visualizar las trayectorias.

## Modelo propuesto (borrador)

### Variables de estado

- **W(t):** agua líquida disponible (superficie y subsuelo), en km³.
- **T(t):** biomasa de truchas de arena, en toneladas.
- **G(t):** biomasa de gusanos adultos, en toneladas.

### Estructura de las ecuaciones

- **Crecimiento de truchas:** dT/dt crece proporcionalmente a T y a una función saturante de W (tipo Michaelis-Menten), menos la tasa de maduración hacia gusano.
- **Retención de agua:** el agua consumida por las truchas se descuenta de W y se acumula en su biomasa y en reservas subterráneas (conservación de masa de agua).
- **Maduración:** una tasa lenta de T hacia G, que permite representar el retardo de milenios dentro de una EDO.
- **Gusanos:** dG/dt aumenta por la maduración y disminuye por mortalidad; al morir, liberan truchas (término a decidir: tasa promedio continua o evento discreto).

### Umbral de desierto

Se define como el instante en que W cae por debajo de una fracción de su valor inicial. Esta fracción se calibrará con un dato de referencia (por ejemplo, humedad del suelo del Sahara) que se buscará en fuentes confiables. Debe justificarse como convención, ya que una Tierra totalmente seca sería más árida que el Sahara.

## Datos y parámetros

- **Datos reales a buscar:** volumen de agua líquida de la Tierra, superficie del planeta, referencia de humedad del Sahara.
- **Parámetros ficticios (supuestos declarados):** masa de una trucha, truchas por gusano, consumo de agua por unidad de biomasa, tiempo de maduración (supuesto de aprox. 4.000 años), mortalidad del gusano.
- Los valores ficticios se tratan como parámetros libres, y el estudio se centra en cómo cambia el resultado al variarlos.

## Limitaciones declaradas

- El ciclo real es discreto (largos periodos y explosiones reproductivas); una EDO continua lo aproxima con tasas promedio.
- Herbert no proporciona cifras; los parámetros son supuestos.
- El resultado es condicional a los supuestos: ilustra un mecanismo, no demuestra nada sobre la Tierra real.
- La hipótesis de colapso hacia los 8.000 años es una conjetura inicial que el modelo puede refutar; un resultado distinto es igualmente válido.

## Relación con Ingeniería de Sistemas y Computación

El sistema es un lazo de retroalimentación con un recurso limitado. Su análisis combina modelado matemático, estudio cualitativo de estabilidad, simulación numérica y análisis de sensibilidad, que son herramientas centrales del análisis de sistemas. La interpretación conecta con la relación entre eficiencia, fragilidad y colapso en sistemas complejos.

## Plan tentativo (dos semanas)

- **Semana 1:** formulación del modelo, búsqueda de datos, equilibrios y linealización.
- **Semana 2:** simulación numérica (extra), análisis de sensibilidad, interpretación y preparación de la exposición.
