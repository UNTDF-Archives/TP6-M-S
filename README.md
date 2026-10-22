# TP6 - Modelos y Simulación (2026)
## Bioimpedancia Eléctrica y Modelo de Cole

Aplicación interactiva y motor de cálculo numérico desarrollados para resolver íntegramente las actividades del Trabajo Práctico N° 6.

---

### Estructura del Proyecto

- [cole_engine.py](file:///home/lk/Proyectos/TP6-MS/cole_engine.py): Motor matemático puro con la ecuación de Cole, cálculo de componentes complejas ($Re$, $Im$, $|Z|$, $\theta$), carga estructurada de los 1457 registros de `registros_EIS.mat` y ajuste no lineal ponderado por mínimos cuadrados (`scipy.optimize.least_squares`).
- [app.py](file:///home/lk/Proyectos/TP6-MS/app.py): Interfaz web en tiempo real desarrollada con **Streamlit + Plotly**.
- [Data/](file:///home/lk/Proyectos/TP6-MS/Data): Dataset original de Dell'Osa et al. (Zenodo DOI: 10.5281/zenodo.15795175).

---



### Inicializacion

# 1. Crear entorno virtual
```bash
python3 -m venv .venv

# 2. Activar entorno
source .venv/bin/activate
## (o en fish)
source .venv/bin/activate.fish


# 3. Instalar librerías del proyecto
pip install streamlit plotly pandas scipy numpy

# 4. Iniciar aplicación
streamlit run app.py
```

### Ejecución de la Aplicación

Para iniciar el servidor interactivo:

```bash
source .venv/bin/activate
streamlit run app.py

## O iniciar directamente con

.venv/bin/streamlit run app.py
.venv/bin/streamlit run app.py --server.address 127.0.0.1

```

Luego abra en el navegador: [http://localhost:8501](http://localhost:8501)

---

### Módulos Implementados

1. **Punto 1: Ecuación y Explorador Cole**:
   - $Z(\omega) = R_\infty + \frac{R_0 - R_\infty}{1 + (j \omega \tau)^\alpha}$.
   - Controles dinámicos (sliders en tiempo real) para $R_0$, $R_\infty$, $f_c$, $\alpha$ y rango frecuencial.
   - Diagrama de Nyquist / Cole-Cole interactivo con escala isométrica 1:1, punto crítico $f_c$, y diagramas de Bode (magnitud $|Z|$ y fase $\theta$).

2. **Punto 2: Simulación de los 5 Tejidos (Tabla 1)**:
   - Tabla interactiva editable para Músculo Longitudinal, Músculo Transversal, Hígado, Pulmón y Bazo.
   - Superposición en un mismo espectro de frecuencias (diagrama Cole-Cole y Bode multicurva).
   - Paleta de color científica de alto contraste (sin colores saturados genéricos de IA).

3. **Punto 3: Ajuste de Curvas con Dataset Zenodo (Dell'Osa et al.)**:
   - Selector filtrable entre los 1457 registros experimentales de la base abierta.
   - Algoritmo de optimización *Trust Region Reflective* con estimaciones automáticas de parámetros y soporte de convenciones de signo.
   - Métricas de bondad de ajuste: $R^2$ total, $R^2_{Re}$, $R^2_{Im}$, RMSE y gráfico de residuos.

4. **Punto 4: Síntesis Global y Comparador**:
   - Tabla de parámetros característicos extraídos para múltiples muestras biológicas y fantomas.
   - Interpretación biofísica de las variables de Cole.
