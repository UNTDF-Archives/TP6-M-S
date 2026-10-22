# Guía Teórica Completa: Bioimpedancia Eléctrica y Modelo de Cole
## Trabajo Práctico N° 6 — Modelos y Simulación (2026)

Este documento explica de forma clara, intuitiva y paso a paso cada concepto, fórmula, parámetro biológico y gráfico utilizado en el TP.

---

## 1. ¿Qué es la Bioimpedancia Eléctrica (EIS)?

La **impedancia ($Z$)** es la oposición total que presenta un material o tejido al paso de una corriente eléctrica alterna ($AC$). A diferencia de un circuito puramente resistivo (como una estufa), un tejido biológico contiene:
1. **Líquidos con iones** (agua con sales dentro y fuera de las células), que se comportan como **resistencias puras** ($R$).
2. **Membranas celulares lipídicas** (bicapas de grasa muy delgadas que aíslan el interior del exterior), que actúan como pequeños **capacitores** ($C$).

Al inyectar una señal eléctrica senoidal de frecuencia variable $f$ (en Hertz, Hz) y medir la respuesta, obtenemos una impedancia compleja:

$$Z(f) = R(f) + j \cdot X(f)$$

- **Parte Real $R(f)$**: Resistencia óhmica (disipación de energía por movimiento de iones).
- **Parte Imaginaria $X(f)$**: Reactancia capacitiva (almacenamiento de carga en las membranas celulares). En medios biológicos capacitivos, esta reactancia es negativa: $X(f) < 0$.
- **$j = \sqrt{-1}$**: La unidad imaginaria.

---

## 2. La Ecuación de Cole (Punto 1 del TP)

En 1941, los hermanos Kenneth y Robert Cole demostraron que los tejidos vivos no se comportan como un capacitor ideal (relajación de Debye pura), sino que exhiben una dispersión dieléctrica más ancha debido a que las células tienen distintas formas, orientaciones y tamaños. Propusieron la fórmula:

$$Z(\omega) = R_\infty + \frac{R_0 - R_\infty}{1 + (j \omega \tau)^\alpha}$$

Donde:
- $\omega = 2\pi f$ es la **frecuencia angular** (rad/s).
- $\tau = \frac{1}{2\pi f_c}$ es la **constante de tiempo de relajación** (segundos).

### ¿Qué significa cada parámetro en el cuerpo humano o animal?

| Parámetro | Nombre | Unidad | Significado Biológico y Físico Intuitivo |
| :--- | :--- | :--- | :--- |
| **$R_0$** | Resistencia a frecuencia cero ($f \to 0$, DC) | Ohms ($\Omega$) | A muy baja frecuencia, las membranas celulares actúan como aislantes perfectos (los capacitores se abren). La corriente eléctrica **no puede entrar a las células** y se ve forzada a viajar exclusivamente por el líquido extracelular (por fuera). Por ende, $R_0$ es alta. |
| **$R_\infty$** | Resistencia a frecuencia infinita ($f \to \infty$) | Ohms ($\Omega$) | A frecuencias muy altas, las membranas celulares ya no oponen resistencia capacitiva ($X_C \to 0$, actúan como un cortocircuito). La corriente atraviesa libremente **tanto el interior celular (intracelular) como el exterior (extracelular)** en paralelo. Por ende, $R_\infty$ es la resistencia mínima del tejido. |
| **$\Delta R = R_0 - R_\infty$** | Amplitud de dispersión | Ohms ($\Omega$) | Representa la contribución resistiva neta ligada a la presencia de las membranas celulares intactas. Si las células mueren o se rompen (lisis), $\Delta R \to 0$. |
| **$f_c$** | Frecuencia característica | Hertz (Hz / kHz) | Es la frecuencia exacta donde la reactancia capacitiva de las células alcanza su valor máximo. Marca el centro de la transición entre baja y alta frecuencia. |
| **$\tau$** | Tiempo de relajación | Microsegundos ($\mu s$) | El tiempo promedio que tardan las cargas en redistribuirse sobre las membranas celulares ante un cambio del campo eléctrico. Se calcula como $\tau = \frac{1}{2\pi f_c}$. |
| **$\alpha$** | Exponente de dispersión | Adimensional ($0 < \alpha \le 1$) | Mide la **heterogeneidad** del tejido. Si todas las células fuesen esferas idénticas del mismo tamaño, $\alpha = 1$ (modelo ideal de Debye, semicírculo perfecto). En tejidos reales con células alargadas, vasos sanguíneos y fibras, $\alpha$ suele oscilar entre $0.20$ y $0.80$, lo que aplana el semicírculo y deprime su centro bajo el eje real. |

---

## 3. Para qué sirve cada Gráfico

La bioimpedancia se analiza principalmente con dos tipos de diagramas:

### A. Diagrama de Nyquist / Cole-Cole ($-Im(Z)$ vs $Re(Z)$)
- **Eje Horizontal ($X$)**: Resistencia $Re(Z)$ en Ohms.
- **Eje Vertical ($Y$)**: Reactancia opuesta $-Im(Z)$ en Ohms (se grafica con signo negativo para que el semicírculo quede hacia arriba, en el primer cuadrante).
- **Cómo leerlo**:
  - Al extremo derecho (baja frecuencia, $f \to 0$), la curva toca el eje $X$ en el valor **$R_0$**.
  - Al extremo izquierdo (alta frecuencia, $f \to \infty$), la curva toca el eje $X$ en el valor **$R_\infty$**.
  - El punto más alto del arco (el vértice) ocurre exactamente a la frecuencia característica **$f_c$**.
  - Cuanto más chico sea $\alpha$, más achatado o aplastado es el arco.

### B. Diagrama de Bode (Magnitud y Fase vs Frecuencia)
- **Gráfico 1: Magnitud $|Z|$ vs $\log_{10}(f)$**: Muestra cómo la oposición total del tejido disminuye progresivamente desde $R_0$ (plataforma a baja frecuencia) hasta $R_\infty$ (plataforma a alta frecuencia). La caída más pronunciada se ubica en $f_c$.
- **Gráfico 2: Fase $\theta$ vs $\log_{10}(f)$**: Muestra el desfase angular entre el voltaje y la corriente. En DC es $0^\circ$ (resistivo puro), alcanza un mínimo negativo (desfasaje capacitivo) cerca de $f_c$, y vuelve a acercarse a $0^\circ$ a muy alta frecuencia.

---

## 4. Análisis del Punto 2: Comparación de los 5 Tejidos (Tabla 1)

El TP solicita comparar los siguientes 5 sistemas biológicos:

1. **Músculo Longitudinal ($R_0=385\,\Omega, R_\infty=70\,\Omega, f_c=37\,\text{kHz}, \alpha=0.27$)**:
   - Las fibras musculares conducen la corriente a lo largo de su eje con menor resistencia que de forma transversal.
2. **Músculo Transversal ($R_0=730\,\Omega, R_\infty=75\,\Omega, f_c=27\,\text{kHz}, \alpha=0.24$)**:
   - Al medir transversalmente a las fibras, la corriente choca contra muchas más capas de membrana por milímetro; por eso $R_0$ ($730\,\Omega$) casi duplica al músculo longitudinal.
3. **Hígado ($R_0=540\,\Omega, R_\infty=100\,\Omega, f_c=72\,\text{kHz}, \alpha=0.45$)**:
   - Tejido parenquimatoso con hepatocitos ordenados en lóbulos. Posee un $\alpha$ más alto ($0.45$), lo que genera arcos de Cole con menos depresión geométrica.
4. **Pulmón ($R_0=440\,\Omega, R_\infty=70\,\Omega, f_c=100\,\text{kHz}, \alpha=0.59$)**:
   - Gran presencia de alvéolos llenos de aire (medio no conductor) y vasculatura capilar muy densa. Posee el $\alpha$ más alto de la tabla ($0.59$).
5. **Bazo ($R_0=550\,\Omega, R_\infty=60\,\Omega, f_c=373\,\text{kHz}, \alpha=0.36$)**:
   - Pulpa roja y blanca con abundancia de eritrocitos y linfocitos pequeños. Su frecuencia característica es notablemente más alta ($373\,\text{kHz}$), desplazando su pico de relajación hacia la derecha en el espectro.

---

## 5. Análisis del Punto 3: Ajuste de Curvas con Dataset Experimental (Dell'Osa et al.)

### ¿Qué se hace en este punto?
En una situación real de laboratorio o clínica médica, **no conocemos $R_0$, $R_\infty$, $f_c$ ni $\alpha$ de antemano**. Lo único que entrega el equipo de medición (ej. Solartron, Zurich Instruments, Impedimed) es una tabla con:
- Una lista de frecuencias aplicadas $f_i$ (ej. desde $1\,\text{Hz}$ hasta $65\,\text{kHz}$).
- La parte real medida $Re_{exp}(f_i)$.
- La parte imaginaria medida $Im_{exp}(f_i)$.

El objetivo es encontrar matemáticamente cuáles son los valores de $R_0, R_\infty, f_c, \alpha$ que mejor reproducen esos puntos medidos.

### Método Matemático Utilizado: Mínimos Cuadrados No Lineales Ponderados
Se utiliza la función de optimización de SciPy (`scipy.optimize.least_squares` con algoritmo *Trust Region Reflective*), que minimiza la suma cuadrática de las diferencias relativas normalizadas:

$$\min_{R_0, R_\infty, f_c, \alpha} \sum_{i=1}^{N} \left[ \left(\frac{Re_{teo}(f_i) - Re_{exp}(f_i)}{|Z_{exp}(f_i)|}\right)^2 + \left(\frac{Im_{teo}(f_i) - Im_{exp}(f_i)}{|Z_{exp}(f_i)|}\right)^2 \right]$$

### Cómo juzgar la calidad del ajuste:
- **Coeficiente de Determinación ($R^2$)**:
  - $R^2 \ge 0.99$: Ajuste excelente. El modelo de Cole describe con fidelidad casi absoluta la respuesta del tejido o fantoma.
  - $R^2 \approx 0.95 - 0.98$: Buen ajuste; pequeñas dispersiones por efectos de contacto de electrodos en bajas frecuencias o capacitancias parásitas a muy altas frecuencias.
- **Gráfico de Residuos**: Muestra el error punto a punto ($Z_{exp} - Z_{teo}$). Si los residuos están distribuidos aleatoriamente cerca de cero sin patrones sistemáticos evidentes, el ajuste es óptimo.

---

## 6. El Dataset de Zenodo (Dell'Osa et al.): ¿Qué es, para qué sirve y cómo razonarlo?

### A. ¿De qué se trata este dataset?
El archivo `registros_EIS.mat` proviene de la publicación abierta en Zenodo:
> *Dell'Osa, A. H., Bravo, N. Y., Villanueva Jousset, T. (2025). Single dispersion EIS measurements: a compilation.*

Es una compilación histórica y multicéntrica de **1457 ensayos de Espectroscopía de Bioimpedancia** realizados entre 1992 y 2018 con equipos de alta precisión (Solartron, Zurich Instruments HF2IS, Impedimed SFB7, Xitron 4001B).

Contiene mediciones sobre:
- **Tejido humano in-vivo** (mediciones de cuerpo entero, extremidades, composición corporal).
- **Tejidos animales y ex-vivo** (ej. huesos, músculo cárnico).
- **Muestras vegetales** (frutas, tubérculos para estudiar membranas vegetales).
- **Fantomas de laboratorio** (geles de agar-agar que simulan propiedades eléctricas conocidas y controladas).

### B. ¿Cómo está estructurado internamente cada registro?
Cada uno de los 1457 ensayos contiene una submatriz llamada `eis` con 5 columnas:
1. **Columna 1 — Frecuencia ($f$)**: El barrido de frecuencias inyectado en Hz (por ejemplo, desde $1\,\text{Hz}$ hasta $65\,\text{kHz}$, o hasta $1.2\,\text{MHz}$).
2. **Columna 2 — Resistencia Real ($Re(Z)$)**: Cuántos Ohms resistivos opone la muestra a esa frecuencia.
3. **Columna 3 — Reactancia Imaginaria ($Im(Z)$)**: La componente capacitiva en Ohms.
4. **Columna 4 — Módulo ($|Z| = \sqrt{Re^2 + Im^2}$)**: La amplitud total de la impedancia en Ohms.
5. **Columna 5 — Fase ($\theta = \arctan(Im/Re)$)**: El desfasaje angular en grados.

### C. ¿Por qué lo estamos utilizando en el TP? (El propósito didáctico y profesional)
En los Puntos 1 y 2 del TP trabajamos en sentido **directo (simulación hacia adelante)**:
> *"Tengo los parámetros teóricos ($R_0, R_\infty, f_c, \alpha$) $\longrightarrow$ genero la curva perfecta de $Z(f)$."*

Pero en la práctica biomédica y científica real ocurre el **problema inverso (identificación de sistemas / ajuste de curvas)**:
> *"Tengo sólo los números medidos con ruido e imperfecciones por un aparato $\longrightarrow$ necesito deducir cuáles son los parámetros biológicos de ese paciente o tejido."*

### D. ¿Cómo razonar los resultados al mirar una muestra?
Cuando seleccionas un registro en la pestaña del Punto 3 y el algoritmo encuentra los valores:

1. **Si $R_0$ da alto y $R_\infty$ da bajo**:
   - Significa que las membranas celulares están **intactas y sanas**. La diferencia $\Delta R = R_0 - R_\infty$ es grande porque hay una pared celular sólida que bloquea la corriente a baja frecuencia y la deja pasar a alta frecuencia.
2. **Si $\Delta R \approx 0$ o $R_0 \approx R_\infty$**:
   - Indica ausencia de membranas o tejido no viable (células rotas, o solución salina sin células).
3. **El valor de $f_c$**:
   - Dice a qué velocidad promedio responden las células. Si la muestra es un fantoma o tejido denso, $f_c$ puede estar en torno a los $15-50\,\text{kHz}$.
4. **El valor de $\alpha$**:
   - Si $\alpha \approx 0.90 - 0.95$, la muestra es sumamente homogénea (por ejemplo, un gel de agar-agar bien mezclado o suspensión uniforme).
   - Si $\alpha \approx 0.40 - 0.60$, la muestra tiene muchas estructuras distintas mezcladas (células de diferente forma, líquido intersticial y fibras conectivas).
