# Informe de verificación — TP6 Modelo de Cole

Fecha: 2026-10-27. Reproducible con `python informe/verificacion.py` (desde la raíz del repo).

## Veredicto

El código es correcto: fórmula, Tabla 1, lectura del dataset y ajuste verificados. Hay 3 observaciones menores (ninguna rompe resultados).

## 1. Enunciado vs código

| Ítem | Enunciado | Código | Estado |
|---|---|---|---|
| Ecuación | `Z = R∞ + (R0−R∞)/(1 − (jωτ)^α)` | `... /(1 + (jωτ)^α)` | Difiere a propósito (ver abajo) |
| `τ = 1/(2π fc)`, `ω = 2πf` | sí | sí (`cole_engine.py:83-84`) | OK |
| Tabla 1 (5 tejidos × 4 parámetros) | — | `cole_engine.py:12-48` | OK, valores idénticos |
| Punto 2: todo en un mismo espectro | — | sección 2 (Nyquist + Bode superpuestos) | OK |
| Punto 3: ajuste con datos Dell'Osa | — | sección 3 y 4 | OK |

**Signo `−` del PDF:** con `1 − (jωτ)^α`, para músculo longitudinal `Re(Z)` llega a −115 Ω e `Im(Z)` es positivo (no es una impedancia pasiva). Se usa `+`, que es el modelo de Cole estándar. El código usó `+` desde el primer commit (`234ca08`). Conviene aclararlo en el informe del TP o consultarlo a la cátedra. El PDF actual es idéntico al del commit inicial; no hay versión anterior para comparar.

## 2. Fórmula (Punto 1 y 2)

Para los 5 tejidos de la Tabla 1:

- `Z(f→0) ≈ R0` y `Z(f→∞) ≈ R∞`, con `Im(Z)` siempre negativo (arco capacitivo). OK.
- `Z(fc)` coincide exacto con la forma cerrada `R∞ + ΔR/(1+e^{jαπ/2})` (error 0.0).
- El máximo de `−Im(Z)` cae en `fc` (p. ej. 36.9998 kHz para 37 kHz; 373.01 kHz para 373 kHz).
- Altura del arco numérica = analítica `ΔR·sin(απ/2)/(2(1+cos(απ/2)))`: 33.91, 62.47, 81.16, 92.44 y 71.18 Ω (coinciden a 1e-9).
- Fase mínima: −9.9°, −11.2°, −16.7°, −24.2° y −16.6°. Crece con α como se espera.

## 3. Dataset (`registros_EIS.mat`)

- 1457 registros leídos, ninguno descartado por el cargador.
- Rangos de origen coherentes con `Data detail.txt`: origen 1 = índices 1310–1456 (147 registros, Dell'Osa).
- Los 5 registros del Punto 4 (índices 1311, 1312, 1350, 1400, 1450) son Dell'Osa, fantoma agar-agar, 50 puntos, 1 Hz–65 kHz. OK.
- Columnas 4 (módulo) y 5 (fase) coinciden con las columnas 2 y 3, salvo 14 registros de origen 5 (Bertemes Filho, índices 1260–1273), donde la fase no coincide con `atan2(Im, Re)`. La app no usa esas columnas (recalcula módulo y fase), así que no afecta.
- En Dell'Osa, `Im(Z)` es negativo casi siempre. En algunos registros pasa a positivo en las frecuencias más altas (artefacto inductivo del equipo; ejemplo: 1350 a 65 kHz, `Im = +1.77 Ω`). El modelo de Cole no lo representa; es normal.

## 4. Ajuste (Punto 3)

Prueba A — **recuperación sintética** (parámetros de la Tabla 1 + ruido de 0.2 Ω, 60 puntos de 100 Hz a 10 MHz):

| Tejido | R0 ajustado/real | R∞ | fc [kHz] | α |
|---|---|---|---|---|
| Músculo Long. | 384.8 / 385 | 69.99 / 70 | 37.21 / 37 | 0.270 / 0.27 |
| Músculo Trans. | 729.7 / 730 | 75.05 / 75 | 27.11 / 27 | 0.240 / 0.24 |
| Hígado | 539.9 / 540 | 100.01 / 100 | 72.15 / 72 | 0.450 / 0.45 |
| Pulmón | 439.9 / 440 | 70.01 / 70 | 100.20 / 100 | 0.590 / 0.59 |
| Bazo | 549.8 / 550 | 60.35 / 60 | 372.9 / 373 | 0.361 / 0.36 |

Error menor a 1 % en todo. El ajuste es correcto.

Prueba B — **ajuste independiente** (`scipy.curve_fit`, 9 puntos iniciales distintos, residuo complejo sin ponderar) sobre los 147 registros de Dell'Osa:

- 147 de 147 ajustes convergen, ninguno con α pegado a un límite.
- R²: mínimo 0.9546, mediana 0.9980, máximo 0.9998.
- Diferencia relativa mediana app vs ajuste independiente: R0 1.1 %, R∞ 0.5 %, fc 4.0 %, α 1.0 %. Percentil 95: 9.8 %, 1.6 %, 25 %, 6.0 %. La diferencia viene de la ponderación distinta (la app pondera por módulo); es esperable.

Los 5 registros del Punto 4 (app):

| Registro | R0 | R∞ | fc [kHz] | α | R² |
|---|---|---|---|---|---|
| 1312 | 153.12 | 68.18 | 0.016 | 0.677 | 0.9981 |
| 1313 | 148.93 | 67.98 | 0.018 | 0.693 | 0.9979 |
| 1351 | 404.57 | 68.14 | 0.001 | 0.732 | 0.9546 |
| 1401 | 71.19 | 58.14 | 0.047 | 0.590 | 0.9994 |
| 1451 | 605.63 | 73.88 | 0.001 | 0.811 | 0.9911 |

### Modelo y método de ajuste (especificación exacta)

**Modelo ajustado:** el modelo de Cole de una sola dispersión, de 4 parámetros `θ = (R0, R∞, fc, α)`, en forma de impedancia compleja:

`Z(ω; θ) = R∞ + (R0 − R∞) / (1 + (jωτ)^α)`, con `ω = 2πf` y `τ = 1/(2π·fc)`.

No es un modelo de circuito equivalente con elementos separados. Es la ecuación del enunciado, con el signo `+` en el denominador.

**Método:** mínimos cuadrados no lineales (`scipy.optimize.least_squares`, algoritmo Trust Region Reflective, `method="trf"`, hasta 5000 evaluaciones, derivadas por diferencias finitas). Es un ajuste simultáneo sobre la parte real y la imaginaria.

**Función que se minimiza** (ponderación "por módulo", la que usa el Punto 4):

`S(θ) = Σᵢ [ (Re Zmedᵢ − Re Z(ωᵢ;θ))² + (Im Zmedᵢ − Im Z(ωᵢ;θ))² ] / wᵢ²`, con `wᵢ = |Zmedᵢ| + 10⁻⁴`.

Cada residuo se divide por el módulo medido, por lo que es un error relativo. La app ofrece otra opción ("por componentes independientes"), donde `Re` se divide por `|Re|+1` e `Im` por `|Im|+1`.

**Restricciones (límites):** `R0` e `R∞` positivos; `0.05 ≤ α ≤ 1`; `fc` entre `0.01·f_min` y `100·f_max` de cada registro.

**Punto de partida:** `R0` = percentil 95 de `Re`; `R∞` = percentil 5 de `Re`; `fc` = frecuencia del pico de `−Im` (o la mediana de las frecuencias si el pico queda en un extremo); `α` = 0.65.

**Métricas reportadas:** R² sobre `Re` e `Im` juntos (`r2_total`) y por separado, y RMSE de cada parte. Código: `cole_engine.py:150-276`.

**Ajuste de contraste:** mismo modelo, pero con `scipy.optimize.curve_fit` (también `trf`), residuos sin ponderar en ohmios, `fc` parametrizado como `log10(fc)` (límite inferior 1 Hz) y 9 puntos de partida.

## 5. Ejecución de la app

Se ejecutaron las 4 secciones de `app.py` con `streamlit.testing.AppTest`: 0 excepciones y 0 errores en todas. Las métricas de la sección 3 coinciden con el motor (R² 0.9976 en el primer registro). El servidor `http://10.7.226.199:8501/` responde HTTP 200.

**Capturas:** no se pudieron sacar. Chromium headless captura solo el esqueleto de carga de Streamlit (necesita websocket y JS), y las secciones 2 a 4 requieren clicks. No hay herramienta de automatización de navegador cargada en esta sesión. Los gráficos se validaron por los datos que los generan (secciones 2 y 4), no visualmente. Conviene mirar los gráficos una vez a mano.

## 6. Observaciones (corregidas el 2026-10-27)

1. **Registros 1351 y 1451 (y probablemente otros de Dell'Osa con fc muy baja):** el ajuste da `fc ≈ 1 Hz`, en el borde inferior de la banda medida (1 Hz–65 kHz). El pico de dispersión queda casi fuera de los datos, así que `R0` y `fc` salen mal determinados. El ajuste independiente da `R0 = 323` contra `405` de la app para 1351, y `fc` de 1.97 contra 0.97 Hz. **Corregido:** la sección 3 muestra una advertencia cuando `fc` queda en el borde de la banda medida. Falta mencionarlo en el informe del TP.
2. **Formato de `fc` en la UI:** `fc` se muestra en kHz con 3 decimales; para fc ≈ 1 Hz muestra `0.001 kHz` (`app.py:469`, `app.py:646`) y τ sale como 164 138 µs. Es correcto pero poco legible. **Corregido:** `fc` se muestra en Hz (métrica, tabla del Punto 4 y etiquetas del gráfico de Nyquist).
3. **Gráfico de residuos (`app.py:596-617`):** el título dice "Residuos Relativos" pero se grafica el residuo absoluto de `Re(Z)` en Ω, y no se muestra el residuo de `Im(Z)`. **Corregido:** el gráfico se titula "Residuos [Ω]" y muestra `Re(Z)` e `Im(Z)`.

Cambios de código: `app.py` y `cole_engine.py` (verificado con `AppTest`: 0 excepciones; el registro 1351 dispara la advertencia). Archivos agregados: `informe/INFORME_VERIFICACION.md` y `informe/verificacion.py`.

## 7. Comparación de modelos: Cole vs Cole + inductancia

**Hipótesis:** el modelo de Cole puro no reproduce la reactancia positiva a frecuencias altas (efecto de cables y equipo). Agregar una inductancia en serie debería bajar el error.

**Modelo alternativo:** `Z(ω) = R∞ + (R0 − R∞)/(1 + (jωτ)^α) + jωL`, con `L ≥ 0`. Tiene 5 parámetros (`R0, R∞, fc, α, L`) contra 4 del Cole puro. El Cole puro es el caso `L = 0`. Código: `cole_engine.py` (`cole_inductive_impedance`, `compare_cole_models`).

**Método:** ambos modelos se ajustan al mismo registro con el mismo método (mínimos cuadrados `trf`) y la misma ponderación por módulo. El modelo con `L` arranca desde el resultado del Cole puro. Se comparan con **AIC** (`n·ln(SSE/n) + 2k`, donde `k` es el número de parámetros; menor es mejor). Es una medida que penaliza los parámetros extra. `ΔAIC = AIC(Cole) − AIC(Cole+L)`; valores positivos favorecen a Cole+L. Se toma ΔAIC > 10 como mejora clara y ΔAIC ≤ 2 como diferencia no relevante.

**Resultados sobre los 147 registros de Dell'Osa:**

| Resultado | Registros |
|---|---|
| Cole + L claramente mejor (ΔAIC > 10) | 2 |
| Mejora moderada (2 < ΔAIC ≤ 10) | 16 |
| Sin diferencia relevante (ΔAIC ≤ 2) | 129 |

- Reducción mediana del RMSE: 0.20 % (2.652 Ω → 2.642 Ω).
- `L` estimada: mediana 2.6 µH, máximo 5.7 µH. En 45 registros el ajuste dejó `L = 0`.
- Los parámetros de Cole casi no se mueven al agregar `L`. Ejemplo, registro 1312: `R0` 153.12 → 153.36, `fc` 16.46 → 16.35 Hz.

| Registro (Punto 4) | ΔAIC | L [µH] |
|---|---|---|
| 1312 | 2.5 | 3.75 |
| 1313 | 1.5 | 3.50 |
| 1351 | −1.9 | 4.51 |
| 1401 | −2.0 | 0.00 |
| 1451 | −1.1 | 2.91 |

**Conclusión:** la hipótesis se confirma solo en parte. Agregar `L` mejora el ajuste en 18 de 147 registros, y de forma clara en 2. Para la gran mayoría no justifica el parámetro extra, y la reducción de error es de 0.2 %. **El modelo de Cole puro es el que mejor equilibra ajuste y simplicidad**, y es el que pide el enunciado.

**Por qué ayuda tan poco (explicación probable, no demostrada):** la reactancia positiva aparece solo en los últimos puntos (cerca de 65 kHz), donde `|Z|` ronda 70 Ω y la reactancia vale unos 2 Ω. Pesa poco frente al resto del espectro. Los registros que mejoran tienen algo más de puntos con reactancia positiva (mediana 8 y 6, contra 5 en el resto), pero la correlación entre ΔAIC y el tamaño de ese efecto es casi nula (0.02). El error que queda (RMSE mediano de unos 2.6 Ω) parece venir de otra fuente. No lo investigué. Una candidata es la asimetría del arco (modelo Cole-Davidson), que no probé.

**En la app:** nueva sección "5. Comparación de Modelos (Cole vs Cole + L)", con la tabla de parámetros y AIC de ambos modelos para el registro elegido, Nyquist y fase superpuestos (datos, Cole, Cole+L), un veredicto por registro y el resumen de los 147 con histograma de ΔAIC. Verificada con `AppTest`: 0 excepciones.

