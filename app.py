"""
Interactive Web Application for Bioimpedance Spectroscopy & Cole Model (TP6 - MyS 2026).
Built with Streamlit, Plotly and SciPy.
"""

import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
from typing import Dict, List

import cole_engine as ce

# Page configuration
st.set_page_config(
    page_title="Bioimpedancia | Modelo de Cole (TP6)",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom Styling for polished, high-contrast, theme-aware UI
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');
    
    html, body, [class*="css"] {
        font-family: 'Inter', sans-serif;
    }
    
    /* Theme-adaptive metric card */
    div[data-testid="stMetric"] {
        background-color: var(--secondary-background-color, rgba(30, 41, 59, 0.7));
        border: 1px solid rgba(148, 163, 184, 0.25);
        padding: 14px 18px;
        border-radius: 12px;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.08);
    }
    
    div[data-testid="stMetric"] label {
        color: var(--text-color, #94A3B8) !important;
        font-weight: 500;
        font-size: 0.85rem;
    }

    div[data-testid="stMetric"] div[data-testid="stMetricValue"] {
        color: var(--text-color, #F8FAFC) !important;
        font-weight: 700;
    }
    
    .highlight-card {
        background: linear-gradient(135deg, rgba(15, 23, 42, 0.95) 0%, rgba(30, 41, 59, 0.95) 100%);
        border: 1px solid rgba(56, 189, 248, 0.3);
        color: #F8FAFC;
        padding: 20px 24px;
        border-radius: 14px;
        margin-bottom: 24px;
        box-shadow: 0 4px 16px rgba(0, 0, 0, 0.25);
    }
    
    .badge {
        display: inline-block;
        padding: 4px 12px;
        font-size: 0.82rem;
        font-weight: 600;
        border-radius: 9999px;
        background: rgba(56, 189, 248, 0.15);
        color: #38BDF8;
        border: 1px solid rgba(56, 189, 248, 0.3);
        margin-right: 8px;
    }
</style>
""", unsafe_allow_html=True)


# --- CACHED DATA LOADING ---
@st.cache_data(show_spinner=False)
def get_cached_records():
    return ce.load_dataset("Data/registros_EIS.mat")


records = get_cached_records()

# --- SIDEBAR NAVIGATION ---
st.sidebar.markdown("## Modelos y Simulación")
st.sidebar.caption("Trabajo Práctico N° 6 - Modelo de Cole")

section = st.sidebar.radio(
    "Navegación del TP:",
    [
        "1. Ecuación y Explorador Cole",
        "2. Simulación de los 5 Tejidos (Tabla 1)",
        "3. Ajuste de Curvas con Dataset Zenodo",
        "4. Comparador Global y Reporte",
        "5. Comparación de Modelos (Cole vs Cole + L)"
    ],
    index=0
)

st.sidebar.markdown("---")
theme_choice = st.sidebar.selectbox("Tema visual de gráficos:", ["Dark Modern (Sleek)", "Light Minimal (Paper)"], index=0)
plot_template = "plotly_dark" if "Dark" in theme_choice else "plotly_white"
bg_color = "#0B0F19" if "Dark" in theme_choice else "#FFFFFF"
grid_color = "#1E293B" if "Dark" in theme_choice else "#E2E8F0"

# --- SECTION 1: ECUACIÓN Y EXPLORADOR DINÁMICO ---
if section == "1. Ecuación y Explorador Cole":
    st.title("Punto 1: Modelo Teórico de Cole")
    
    with st.container(border=True):
        st.markdown("<h3 style='margin:0 0 8px 0; color:#38BDF8;'>Ecuación de Cole (1941)</h3>", unsafe_allow_html=True)
        st.markdown("Describe la dispersión dieléctrica en medios biológicos considerando la relajación no ideal de las membranas celulares:")
        st.latex(r"Z(\omega) = R_\infty + \frac{R_0 - R_\infty}{1 + (j \omega \tau)^\alpha}")
        
        c_sub1, c_sub2, c_sub3 = st.columns(3)
        with c_sub1:
            st.latex(r"\omega = 2\pi f")
        with c_sub2:
            st.latex(r"\tau = \frac{1}{2\pi f_c}")
        with c_sub3:
            st.latex(r"\alpha \in (0, 1]")

    col_ctrl, col_plot = st.columns([1, 2])
    
    with col_ctrl:
        st.subheader("Parámetros Dinámicos")
        preset = st.selectbox(
            "Cargar Preset de Tejido:",
            ["Personalizado"] + list(ce.TABLE_1_TISSUES.keys()),
            index=1
        )
        
        if preset in ce.TABLE_1_TISSUES:
            p_data = ce.TABLE_1_TISSUES[preset]
            init_R0 = float(p_data["R0"])
            init_Rinf = float(p_data["Rinf"])
            init_fc = float(p_data["fc"] / 1000.0)
            init_alpha = float(p_data["alpha"])
        else:
            init_R0 = 500.0
            init_Rinf = 80.0
            init_fc = 50.0
            init_alpha = 0.50

        r0_val = st.slider("R₀ (Resistencia a f → 0) [Ω]:", 10.0, 2000.0, init_R0, 5.0)
        rinf_val = st.slider("R∞ (Resistencia a f → ∞) [Ω]:", 5.0, min(r0_val - 1.0, 1000.0), min(init_Rinf, r0_val - 1.0), 1.0)
        fc_val_khz = st.slider("Frecuencia característica fc [kHz]:", 0.1, 1000.0, init_fc, 0.5)
        alpha_val = st.slider("Exponente de dispersión α:", 0.05, 1.00, init_alpha, 0.01)

        st.markdown("#### Rango Frecuencial")
        f_min_exp = st.slider("Frecuencia Mínima log₁₀(Hz):", 0.0, 3.0, 1.0, 0.5)
        f_max_exp = st.slider("Frecuencia Máxima log₁₀(Hz):", 5.0, 9.0, 7.0, 0.5)
        points = st.slider("Puntos de simulación:", 100, 1000, 500, 50)
        
        fc_hz = fc_val_khz * 1000.0
        tau_us = (1.0 / (2.0 * np.pi * fc_hz)) * 1e6
        
        st.markdown("---")
        st.markdown(f"**Constante de relajación τ:** `{tau_us:.2f} µs`")
        st.markdown(f"**ΔR (R₀ - R∞):** `{r0_val - rinf_val:.1f} Ω`")

    with col_plot:
        f_vec = np.logspace(f_min_exp, f_max_exp, points)
        comp = ce.compute_cole_components(f_vec, r0_val, rinf_val, fc_hz, alpha_val)
        
        # Nyquist Plot
        fig_nyq = go.Figure()
        
        # Arc trace
        fig_nyq.add_trace(go.Scatter(
            x=comp["re"],
            y=comp["minus_im"],
            mode="lines",
            line=dict(color="#38BDF8", width=3.5),
            name="Z(ω) Modelo Cole",
            customdata=np.stack([f_vec, comp["modulus"], comp["phase_deg"]], axis=-1),
            hovertemplate="<b>Re(Z):</b> %{x:.2f} Ω<br><b>-Im(Z):</b> %{y:.2f} Ω<br><b>f:</b> %{customdata[0]:.1f} Hz<br><b>|Z|:</b> %{customdata[1]:.2f} Ω<br><b>Fase:</b> %{customdata[2]:.2f}°<extra></extra>"
        ))
        
        # Characteristic point at fc
        comp_fc = ce.compute_cole_components(np.array([fc_hz]), r0_val, rinf_val, fc_hz, alpha_val)
        fig_nyq.add_trace(go.Scatter(
            x=[comp_fc["re"][0]],
            y=[comp_fc["minus_im"][0]],
            mode="markers+text",
            marker=dict(color="#F43F5E", size=11, symbol="diamond"),
            name=f"Vértice en fc = {fc_val_khz:.1f} kHz",
            text=[f"fc={fc_val_khz:.1f}kHz"],
            textposition="top center",
            hovertemplate="<b>Punto Crítico en fc</b><br>Re: %{x:.2f} Ω<br>-Im: %{y:.2f} Ω<extra></extra>"
        ))
        
        # Theoretical R0 and Rinf markers on real axis
        fig_nyq.add_trace(go.Scatter(
            x=[rinf_val, r0_val],
            y=[0, 0],
            mode="markers+text",
            marker=dict(color="#10B981", size=8, symbol="cross"),
            text=["R∞", "R₀"],
            textposition="bottom center",
            name="Límites asíntoticos"
        ))

        fig_nyq.update_layout(
            title="Diagrama de Nyquist / Cole-Cole (-Im(Z) vs Re(Z))",
            xaxis_title="Resistencia Re(Z) [Ω]",
            yaxis_title="Reactancia -Im(Z) [Ω]",
            template=plot_template,
            plot_bgcolor=bg_color,
            paper_bgcolor=bg_color,
            hovermode="closest",
            height=460,
            yaxis=dict(scaleanchor="x", scaleratio=1, gridcolor=grid_color),
            xaxis=dict(gridcolor=grid_color)
        )
        st.plotly_chart(fig_nyq, width='stretch')

        # Bode Plots (Magnitud y Fase)
        fig_bode = make_subplots(
            rows=2, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.1,
            subplot_titles=["Magnitud |Z| vs Frecuencia", "Ángulo de Fase [°] vs Frecuencia"]
        )

        fig_bode.add_trace(
            go.Scatter(
                x=f_vec, y=comp["modulus"],
                mode="lines", line=dict(color="#06B6D4", width=2.5),
                name="|Z|",
                hovertemplate="f: %{x:.1f} Hz<br>|Z|: %{y:.2f} Ω<extra></extra>"
            ),
            row=1, col=1
        )

        fig_bode.add_trace(
            go.Scatter(
                x=f_vec, y=comp["phase_deg"],
                mode="lines", line=dict(color="#A855F7", width=2.5),
                name="Fase [°]",
                hovertemplate="f: %{x:.1f} Hz<br>Fase: %{y:.2f}°<extra></extra>"
            ),
            row=2, col=1
        )

        fig_bode.update_xaxes(type="log", title_text="Frecuencia [Hz]", row=2, col=1, gridcolor=grid_color)
        fig_bode.update_xaxes(type="log", row=1, col=1, gridcolor=grid_color)
        fig_bode.update_yaxes(title_text="|Z| [Ω]", row=1, col=1, gridcolor=grid_color)
        fig_bode.update_yaxes(title_text="Fase [°]", row=2, col=1, gridcolor=grid_color)
        fig_bode.update_layout(
            template=plot_template,
            plot_bgcolor=bg_color,
            paper_bgcolor=bg_color,
            height=480,
            showlegend=False
        )
        st.plotly_chart(fig_bode, width='stretch')


# --- SECTION 2: SIMULACIÓN DE LOS 5 TEJIDOS (TABLA 1) ---
elif section == "2. Simulación de los 5 Tejidos (Tabla 1)":
    st.title("Punto 2: Espectro Comparativo de 5 Tejidos Biológicos")
    st.markdown("Comparación simultánea en un mismo espectro frecuencial conforme a los valores característicos de la **Tabla 1** del enunciado.")

    # Interactive Table Editor
    st.subheader("Parámetros de la Tabla 1 (Modificables en Tiempo Real)")
    
    table_data = []
    for name, vals in ce.TABLE_1_TISSUES.items():
        table_data.append({
            "Tejido": name,
            "R0 [Ω]": vals["R0"],
            "Rinf [Ω]": vals["Rinf"],
            "fc [kHz]": vals["fc"] / 1000.0,
            "alpha": vals["alpha"],
            "Color": vals["color"]
        })
    df_tissues = pd.DataFrame(table_data)

    edited_df = st.data_editor(
        df_tissues,
        num_rows="dynamic",
        width="stretch",
        key="tissue_editor"
    )

    # Frequency Range selection
    col_f1, col_f2 = st.columns(2)
    with col_f1:
        f_min_table = st.select_slider("Frecuencia Mínima de Simulación:", options=[1, 10, 100, 1000, 10000], value=10)
    with col_f2:
        f_max_table = st.select_slider("Frecuencia Máxima de Simulación:", options=[100000, 500000, 1000000, 10000000, 100000000], value=10000000)

    f_sim = np.logspace(np.log10(f_min_table), np.log10(f_max_table), 600)

    # Nyquist Multi-curve
    fig_multi_nyq = go.Figure()
    
    # Bode Multi-curve
    fig_multi_bode = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.1,
        subplot_titles=["Magnitud de Impedancia |Z| vs Frecuencia", "Fase [°] vs Frecuencia"]
    )

    for _, row in edited_df.iterrows():
        t_name = str(row["Tejido"])
        r0 = float(row["R0 [Ω]"])
        rinf = float(row["Rinf [Ω]"])
        fc_hz = float(row["fc [kHz]"]) * 1000.0
        alpha = float(row["alpha"])
        c_code = str(row["Color"])

        comp_t = ce.compute_cole_components(f_sim, r0, rinf, fc_hz, alpha)

        # Trace Nyquist
        fig_multi_nyq.add_trace(go.Scatter(
            x=comp_t["re"],
            y=comp_t["minus_im"],
            mode="lines",
            line=dict(color=c_code, width=3),
            name=t_name,
            customdata=np.stack([f_sim, comp_t["modulus"], comp_t["phase_deg"]], axis=-1),
            hovertemplate=f"<b>{t_name}</b><br>Re: %{{x:.1f}} Ω<br>-Im: %{{y:.1f}} Ω<br>f: %{{customdata[0]:.0f}} Hz<extra></extra>"
        ))

        # Peak marker at fc
        comp_fc_t = ce.compute_cole_components(np.array([fc_hz]), r0, rinf, fc_hz, alpha)
        fig_multi_nyq.add_trace(go.Scatter(
            x=[comp_fc_t["re"][0]],
            y=[comp_fc_t["minus_im"][0]],
            mode="markers",
            marker=dict(color=c_code, size=9, symbol="diamond"),
            showlegend=False,
            name=f"{t_name} (fc)",
            hovertemplate=f"<b>{t_name} en fc</b> ({row['fc [kHz]']} kHz)<br>Re: %{{x:.1f}} Ω<br>-Im: %{{y:.1f}} Ω<extra></extra>"
        ))

        # Bode traces
        fig_multi_bode.add_trace(
            go.Scatter(
                x=f_sim, y=comp_t["modulus"],
                mode="lines", line=dict(color=c_code, width=2.5),
                name=t_name,
                hovertemplate=f"<b>{t_name}</b><br>f: %{{x:.0f}} Hz<br>|Z|: %{{y:.1f}} Ω<extra></extra>"
            ),
            row=1, col=1
        )
        fig_multi_bode.add_trace(
            go.Scatter(
                x=f_sim, y=comp_t["phase_deg"],
                mode="lines", line=dict(color=c_code, width=2.5),
                showlegend=False,
                name=t_name,
                hovertemplate=f"<b>{t_name}</b><br>f: %{{x:.0f}} Hz<br>Fase: %{{y:.2f}}°<extra></extra>"
            ),
            row=2, col=1
        )

    fig_multi_nyq.update_layout(
        title="Superposición de Tejidos: Diagrama Cole-Cole (-Im(Z) vs Re(Z))",
        xaxis_title="Re(Z) [Ω]",
        yaxis_title="-Im(Z) [Ω]",
        template=plot_template,
        plot_bgcolor=bg_color,
        paper_bgcolor=bg_color,
        height=520,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        xaxis=dict(gridcolor=grid_color),
        yaxis=dict(gridcolor=grid_color)
    )
    st.plotly_chart(fig_multi_nyq, width='stretch')

    fig_multi_bode.update_xaxes(type="log", title_text="Frecuencia [Hz]", row=2, col=1, gridcolor=grid_color)
    fig_multi_bode.update_xaxes(type="log", row=1, col=1, gridcolor=grid_color)
    fig_multi_bode.update_yaxes(title_text="|Z| [Ω]", row=1, col=1, gridcolor=grid_color)
    fig_multi_bode.update_yaxes(title_text="Fase [°]", row=2, col=1, gridcolor=grid_color)
    fig_multi_bode.update_layout(
        title="Espectro Frecuencial Comparativo (Bode Magnitud y Fase)",
        template=plot_template,
        plot_bgcolor=bg_color,
        paper_bgcolor=bg_color,
        height=560,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )
    st.plotly_chart(fig_multi_bode, width='stretch')


# --- SECTION 3: AJUSTE DE CURVAS DATASET ZENODO ---
elif section == "3. Ajuste de Curvas con Dataset Zenodo":
    st.title("Punto 3: Ajuste de Curvas con Dataset Experimental (Dell'Osa et al.)")
    st.markdown("Selección de registros reales del repositorio Zenodo (`registros_EIS.mat`), ajuste no lineal por mínimos cuadrados y extracción de parámetros característicos.")

    col_sel1, col_sel2, col_sel3 = st.columns([1.5, 1.5, 1])
    
    with col_sel1:
        origin_filter = st.selectbox(
            "Filtrar por Origen de Datos:",
            ["Todos"] + list(ce.ORIGIN_LABELS.values()),
            index=1  # Default to Dell'Osa et al.
        )
    with col_sel2:
        sample_filter = st.selectbox(
            "Filtrar por Tipo de Muestra:",
            ["Todas"] + list(ce.SAMPLE_LABELS.values()),
            index=0
        )

    # Filter records
    filtered_records = []
    for r in records:
        if origin_filter != "Todos" and r["origin_label"] != origin_filter:
            continue
        if sample_filter != "Todas" and r["sample_label"] != sample_filter:
            continue
        filtered_records.append(r)

    with col_sel3:
        st.metric("Registros encontrados", len(filtered_records))

    if not filtered_records:
        st.warning("No se encontraron registros con los filtros seleccionados.")
    else:
        # Selection of individual record
        options_map = {
            f"ID {r['id']} | {r['origin_label'][:18]} | {r['sample_label']} ({r['points']} pts)": r
            for r in filtered_records[:300]
        }
        chosen_label = st.selectbox("Seleccione el Registro Experimental a Ajustar:", list(options_map.keys()), index=0)
        rec = options_map[chosen_label]

        # Extract EIS data
        eis_data = rec["eis"]
        f_raw = eis_data[:, 0]
        re_raw = eis_data[:, 1]
        im_raw = eis_data[:, 2]

        col_fit_opts, col_fit_res = st.columns([1, 1])
        with col_fit_opts:
            st.markdown("#### Ajuste de Parámetros")
            weight_type = st.radio("Ponderación del residuo:", ["Por módulo |Z| (Recomendado)", "Por componentes independientes"], index=0)
            weight_param = "modulus" if "módulo" in weight_type else "relative"

            with st.expander("Modificar Estimaciones Iniciales (Opcional)"):
                custom_r0 = st.number_input("R₀ inicial [Ω]:", value=float(np.max(re_raw)), step=5.0)
                custom_rinf = st.number_input("R∞ inicial [Ω]:", value=max(float(np.min(re_raw)), 1.0), step=5.0)
                custom_fc = st.number_input("fc inicial [Hz]:", value=float(f_raw[len(f_raw)//2]), step=100.0)
                custom_alpha = st.slider("α inicial:", 0.1, 1.0, 0.65, 0.05)
                use_custom_p0 = st.checkbox("Forzar estimación inicial personalizada", value=False)

            p0_input = [custom_r0, custom_rinf, custom_fc, custom_alpha] if use_custom_p0 else None

        # Execute fitting
        with st.spinner("Ajustando modelo de Cole a los datos experimentales..."):
            try:
                fit_res = ce.fit_cole_model(f_raw, re_raw, im_raw, p0=p0_input, weights=weight_param)
                params = fit_res["params"]
                metrics = fit_res["metrics"]
            except Exception as e:
                st.error(f"Error durante el ajuste: {str(e)}")
                fit_res = None

        if fit_res:
            with col_fit_res:
                st.markdown("#### Resultados del Ajuste")
                m1, m2, m3, m4 = st.columns(4)
                m1.metric("R₀ [Ω]", f"{params['R0']:.2f}")
                m2.metric("R∞ [Ω]", f"{params['Rinf']:.2f}")
                m3.metric("fc [Hz]", f"{params['fc']:.1f}")
                m4.metric("α", f"{params['alpha']:.3f}")

                m5, m6, m7 = st.columns(3)
                m5.metric("R² Total", f"{metrics['r2_total']:.4f}")
                m6.metric("τ de relajación", f"{params['tau_us']:.2f} µs")
                m7.metric("ΔR (R₀ - R∞)", f"{params['delta_R']:.2f} Ω")

            f_lo, f_hi = fit_res["f_meas"].min(), fit_res["f_meas"].max()
            if not (2 * f_lo <= params["fc"] <= f_hi / 2):
                st.warning(f"fc = {params['fc']:.1f} Hz está en el borde de la banda medida ({f_lo:.0f} Hz – {f_hi:.0f} Hz): R₀ y fc quedan poco determinados.")

            # Dense simulation for fitted model
            f_dense = np.logspace(np.log10(fit_res["f_meas"].min() * 0.5), np.log10(fit_res["f_meas"].max() * 2.0), 500)
            Z_dense = ce.cole_impedance(f_dense, params["R0"], params["Rinf"], params["fc"], params["alpha"])
            re_dense = np.real(Z_dense)
            im_dense = np.imag(Z_dense)
            minus_im_dense = -im_dense

            # Plots: Nyquist Fit & Residuals / Bode
            tab_nyq, tab_bode, tab_res = st.tabs(["Diagrama de Nyquist / Cole", "Diagrama de Bode", "Residuos y Datos"])

            with tab_nyq:
                fig_fit_nyq = go.Figure()

                # Experimental points
                fig_fit_nyq.add_trace(go.Scatter(
                    x=fit_res["re_meas"],
                    y=fit_res["minus_im_meas"],
                    mode="markers",
                    marker=dict(color="#38BDF8", size=8, symbol="circle-open", line=dict(width=2)),
                    name="Datos Experimentales (Zenodo)",
                    customdata=fit_res["f_meas"],
                    hovertemplate="<b>Experimental</b><br>Re: %{x:.2f} Ω<br>-Im: %{y:.2f} Ω<br>f: %{customdata:.1f} Hz<extra></extra>"
                ))

                # Fitted Arc
                fig_fit_nyq.add_trace(go.Scatter(
                    x=re_dense,
                    y=minus_im_dense,
                    mode="lines",
                    line=dict(color="#F43F5E", width=3),
                    name="Curva Ajustada Modelo Cole",
                    customdata=f_dense,
                    hovertemplate="<b>Ajuste Cole</b><br>Re: %{x:.2f} Ω<br>-Im: %{y:.2f} Ω<br>f: %{customdata:.1f} Hz<extra></extra>"
                ))

                # Highlight fc
                Z_fc = ce.cole_impedance(np.array([params["fc"]]), params["R0"], params["Rinf"], params["fc"], params["alpha"])
                fig_fit_nyq.add_trace(go.Scatter(
                    x=[float(np.real(Z_fc)[0])],
                    y=[float(-np.imag(Z_fc)[0])],
                    mode="markers+text",
                    marker=dict(color="#FBBF24", size=10, symbol="star"),
                    name=f"fc = {params['fc']:.1f} Hz",
                    text=[f"fc={params['fc']:.1f}Hz"],
                    textposition="top center"
                ))

                fig_fit_nyq.update_layout(
                    title=f"Ajuste en Diagrama de Nyquist: {chosen_label}",
                    xaxis_title="Re(Z) [Ω]",
                    yaxis_title="-Im(Z) [Ω]",
                    template=plot_template,
                    plot_bgcolor=bg_color,
                    paper_bgcolor=bg_color,
                    height=500,
                    yaxis=dict(scaleanchor="x", scaleratio=1, gridcolor=grid_color),
                    xaxis=dict(gridcolor=grid_color)
                )
                st.plotly_chart(fig_fit_nyq, width='stretch')

            with tab_bode:
                fig_fit_bode = make_subplots(
                    rows=2, cols=1,
                    shared_xaxes=True,
                    vertical_spacing=0.1,
                    subplot_titles=["Magnitud |Z| [Ω]", "Ángulo de Fase [°]"]
                )

                # Experimental
                mod_meas = np.sqrt(fit_res["re_meas"]**2 + fit_res["im_meas"]**2)
                phase_meas = np.rad2deg(np.arctan2(fit_res["im_meas"], fit_res["re_meas"]))

                # Fitted
                mod_dense = np.abs(Z_dense)
                phase_dense = np.angle(Z_dense, deg=True)

                fig_fit_bode.add_trace(
                    go.Scatter(x=fit_res["f_meas"], y=mod_meas, mode="markers", marker=dict(color="#38BDF8", size=6), name="|Z| Exp"),
                    row=1, col=1
                )
                fig_fit_bode.add_trace(
                    go.Scatter(x=f_dense, y=mod_dense, mode="lines", line=dict(color="#F43F5E", width=2.5), name="|Z| Cole"),
                    row=1, col=1
                )

                fig_fit_bode.add_trace(
                    go.Scatter(x=fit_res["f_meas"], y=phase_meas, mode="markers", marker=dict(color="#38BDF8", size=6), name="Fase Exp", showlegend=False),
                    row=2, col=1
                )
                fig_fit_bode.add_trace(
                    go.Scatter(x=f_dense, y=phase_dense, mode="lines", line=dict(color="#F43F5E", width=2.5), name="Fase Cole", showlegend=False),
                    row=2, col=1
                )

                fig_fit_bode.update_xaxes(type="log", title_text="Frecuencia [Hz]", row=2, col=1, gridcolor=grid_color)
                fig_fit_bode.update_xaxes(type="log", row=1, col=1, gridcolor=grid_color)
                fig_fit_bode.update_yaxes(title_text="|Z| [Ω]", row=1, col=1, gridcolor=grid_color)
                fig_fit_bode.update_yaxes(title_text="Fase [°]", row=2, col=1, gridcolor=grid_color)
                fig_fit_bode.update_layout(
                    template=plot_template,
                    plot_bgcolor=bg_color,
                    paper_bgcolor=bg_color,
                    height=520
                )
                st.plotly_chart(fig_fit_bode, width='stretch')

            with tab_res:
                col_df, col_res_plt = st.columns([1, 1])
                with col_df:
                    st.markdown("##### Tabla de Datos")
                    df_out = pd.DataFrame({
                        "Frecuencia [Hz]": fit_res["f_meas"],
                        "Re(Z) Exp": fit_res["re_meas"],
                        "Im(Z) Exp": fit_res["im_meas"],
                        "|Z| Exp": mod_meas,
                    })
                    st.dataframe(df_out, height=350, width='stretch')

                with col_res_plt:
                    st.markdown("##### Residuos [Ω]")
                    # Calculate residuals at experimental frequencies (dataset sign convention)
                    Z_fit_at_f = ce.cole_impedance(fit_res["f_meas"], params["R0"], params["Rinf"], params["fc"], params["alpha"])
                    sign_im = 1.0 if fit_res["sign_im_convention"] == "standard_negative" else -1.0
                    res_re = fit_res["re_meas"] - np.real(Z_fit_at_f)
                    res_im = fit_res["im_meas"] - np.imag(Z_fit_at_f) * sign_im
                    
                    fig_res = go.Figure()
                    fig_res.add_trace(go.Bar(
                        x=fit_res["f_meas"],
                        y=res_re,
                        marker_color="#F43F5E",
                        name="Residual Re(Z)"
                    ))
                    fig_res.add_trace(go.Bar(
                        x=fit_res["f_meas"],
                        y=res_im,
                        marker_color="#38BDF8",
                        name="Residual Im(Z)"
                    ))
                    fig_res.update_layout(barmode="group")
                    fig_res.update_xaxes(type="log", title_text="Frecuencia [Hz]", gridcolor=grid_color)
                    fig_res.update_yaxes(title_text="Residual [Ω]", gridcolor=grid_color)
                    fig_res.update_layout(
                        template=plot_template,
                        plot_bgcolor=bg_color,
                        paper_bgcolor=bg_color,
                        height=350,
                        margin=dict(t=30, b=30)
                    )
                    st.plotly_chart(fig_res, width='stretch')


# --- SECTION 4: COMPARADOR GLOBAL Y REPORTE ---
elif section == "4. Comparador Global y Reporte":
    st.title("Punto 4: Síntesis Global y Reporte de Parámetros")
    st.markdown("Integración de los resultados de simulación teórica (Punto 2) con ajustes experimentales (Punto 3).")

    # Run fit across 3 representative datasets from Dell'Osa
    st.subheader("Ajuste Automatizado de Múltiples Muestras de Dell'Osa et al.")

    sample_indices = [1311, 1312, 1350, 1400, 1450]
    multi_results = []
    
    with st.spinner("Procesando conjunto de muestras seleccionadas..."):
        for s_idx in sample_indices:
            if s_idx < len(records):
                r_item = records[s_idx]
                e = r_item["eis"]
                f_m = e[:, 0]
                re_m = e[:, 1]
                im_m = e[:, 2]
                res_f = ce.fit_cole_model(f_m, re_m, im_m)
                multi_results.append({
                    "ID Muestra": r_item["id"],
                    "Tipo": r_item["sample_label"],
                    "Puntos": r_item["points"],
                    "R₀ [Ω]": f"{res_f['params']['R0']:.2f}",
                    "R∞ [Ω]": f"{res_f['params']['Rinf']:.2f}",
                    "fc [Hz]": f"{res_f['params']['fc']:.1f}",
                    "α": f"{res_f['params']['alpha']:.3f}",
                    "τ [µs]": f"{res_f['params']['tau_us']:.2f}",
                    "R² Total": f"{res_f['metrics']['r2_total']:.4f}",
                })
    
    st.dataframe(pd.DataFrame(multi_results), width='stretch')

    st.markdown("---")
    st.subheader("Resumen Teórico del Trabajo Práctico")
    st.markdown("""
    1. **Ecuación del Modelo de Cole**: Modela la dispersión $\\alpha$ típica de suspensión celular y tejidos conectivos, donde la constante $\\alpha < 1$ refleja la heterogeneidad de tamaños celulares y constantes de tiempo de relajación.
    2. **Significado Físico de Parámetros**:
       - **$R_0$**: Resistencia a bajas frecuencias; la corriente eléctrica no penetra las membranas lipídicas celulares (aislantes a DC), fluyendo casi exclusivamente por el medio extracelular (ECM).
       - **$R_\\infty$**: Resistencia a muy altas frecuencias; la reactancia capacitiva de las bicapas lipídicas tiende a cero ($X_C \\to 0$), permitiendo que la corriente fluya a través del medio intra y extracelular en paralelo.
       - **$f_c$**: Frecuencia de corte o pico de relajación dieléctrica, donde la reactancia capacitiva alcanza su máximo relativo.
       - **$\\alpha$**: Parámetro adimensional que cuantifica el aplanamiento del semicírculo de Cole ($0 < \\alpha \\le 1$). Si $\\alpha = 1$, el modelo se reduce a la relajación ideal de Debye.
    """)


# --- SECTION 5: COMPARACIÓN DE MODELOS (COLE vs COLE + INDUCTANCIA) ---
elif section == "5. Comparación de Modelos (Cole vs Cole + L)":
    st.title("Punto 5: ¿Mejora el ajuste agregar una inductancia?")
    st.markdown(
        "Algunos registros de Dell'Osa muestran reactancia positiva a frecuencias altas (efecto de cables y equipo). "
        "El modelo de Cole puro no lo representa. Se compara con **Cole + inductancia en serie**: "
        "`Z = Z_Cole + jωL` (5 parámetros). Ambos se ajustan con la misma ponderación por módulo."
    )

    @st.cache_data(show_spinner=False)
    def compare_all_dellosa():
        out = []
        for r in records:
            if r["origin"] != "1":
                continue
            e = r["eis"]
            c = ce.compare_cole_models(e[:, 0], e[:, 1], e[:, 2])
            out.append({
                "ID": r["id"],
                "ΔAIC": c["delta_aic"],
                "RMSE Cole [Ω]": c["cole"]["stats"]["rmse"],
                "RMSE Cole+L [Ω]": c["ind"]["stats"]["rmse"],
                "L [µH]": c["ind"]["params"]["L_uH"],
            })
        return pd.DataFrame(out)

    dell = [r for r in records if r["origin"] == "1"]
    labels = {f"ID {r['id']} ({r['sample_label']})": r for r in dell}
    default_idx = next((i for i, r in enumerate(dell) if r["id"] == "1312"), 0)
    rec5 = labels[st.selectbox("Registro:", list(labels.keys()), index=default_idx)]
    e5 = rec5["eis"]
    cmp = ce.compare_cole_models(e5[:, 0], e5[:, 1], e5[:, 2])
    f5, re5 = cmp["f_meas"], cmp["re_meas"]
    im5 = cmp["im_meas"] * cmp["sign_im"]  # convención teórica (negativa)
    pc, pi_ = cmp["cole"]["params"], cmp["ind"]["params"]
    sc, si = cmp["cole"]["stats"], cmp["ind"]["stats"]

    st.dataframe(pd.DataFrame({
        "Modelo": ["Cole (4 parámetros)", "Cole + L (5 parámetros)"],
        "R₀ [Ω]": [f"{pc['R0']:.2f}", f"{pi_['R0']:.2f}"],
        "R∞ [Ω]": [f"{pc['Rinf']:.2f}", f"{pi_['Rinf']:.2f}"],
        "fc [Hz]": [f"{pc['fc']:.1f}", f"{pi_['fc']:.1f}"],
        "α": [f"{pc['alpha']:.3f}", f"{pi_['alpha']:.3f}"],
        "L [µH]": ["—", f"{pi_['L_uH']:.2f}"],
        "RMSE [Ω]": [f"{sc['rmse']:.3f}", f"{si['rmse']:.3f}"],
        "AIC": [f"{sc['aic']:.1f}", f"{si['aic']:.1f}"],
        "BIC": [f"{sc['bic']:.1f}", f"{si['bic']:.1f}"],
    }), width="stretch", hide_index=True)

    d = cmp["delta_aic"]
    if d > 10:
        st.success(f"ΔAIC = {d:.1f}: Cole + L ajusta claramente mejor en este registro.")
    elif d > 2:
        st.info(f"ΔAIC = {d:.1f}: Cole + L ajusta algo mejor, evidencia moderada.")
    else:
        st.warning(f"ΔAIC = {d:.1f}: no hay mejora que justifique el parámetro extra. Se prefiere Cole.")

    f_d = np.logspace(np.log10(f5.min()), np.log10(f5.max()), 400)
    Zc_d = ce.cole_impedance(f_d, pc["R0"], pc["Rinf"], pc["fc"], pc["alpha"])
    Zi_d = ce.cole_inductive_impedance(f_d, pi_["R0"], pi_["Rinf"], pi_["fc"], pi_["alpha"], pi_["L_uH"] * 1e-6)

    fig5 = make_subplots(rows=1, cols=2, subplot_titles=["Nyquist: -Im(Z) vs Re(Z)", "Fase [°] vs Frecuencia"], horizontal_spacing=0.1)
    fig5.add_trace(go.Scatter(x=re5, y=-im5, mode="markers", marker=dict(color="#38BDF8", size=7), name="Datos"), row=1, col=1)
    fig5.add_trace(go.Scatter(x=np.real(Zc_d), y=-np.imag(Zc_d), mode="lines", line=dict(color="#F43F5E", width=2.5), name="Cole"), row=1, col=1)
    fig5.add_trace(go.Scatter(x=np.real(Zi_d), y=-np.imag(Zi_d), mode="lines", line=dict(color="#FBBF24", width=2.5, dash="dash"), name="Cole + L"), row=1, col=1)
    fig5.add_trace(go.Scatter(x=f5, y=np.rad2deg(np.arctan2(im5, re5)), mode="markers", marker=dict(color="#38BDF8", size=7), name="Datos", showlegend=False), row=1, col=2)
    fig5.add_trace(go.Scatter(x=f_d, y=np.angle(Zc_d, deg=True), mode="lines", line=dict(color="#F43F5E", width=2.5), name="Cole", showlegend=False), row=1, col=2)
    fig5.add_trace(go.Scatter(x=f_d, y=np.angle(Zi_d, deg=True), mode="lines", line=dict(color="#FBBF24", width=2.5, dash="dash"), name="Cole + L", showlegend=False), row=1, col=2)
    fig5.update_xaxes(title_text="Re(Z) [Ω]", row=1, col=1, gridcolor=grid_color)
    fig5.update_xaxes(type="log", title_text="Frecuencia [Hz]", row=1, col=2, gridcolor=grid_color)
    fig5.update_yaxes(title_text="-Im(Z) [Ω]", row=1, col=1, gridcolor=grid_color)
    fig5.update_yaxes(title_text="Fase [°]", row=1, col=2, gridcolor=grid_color)
    fig5.update_layout(template=plot_template, plot_bgcolor=bg_color, paper_bgcolor=bg_color, height=440)
    st.plotly_chart(fig5, width="stretch")

    st.subheader("Los 147 registros de Dell'Osa")
    with st.spinner("Ajustando ambos modelos en todos los registros..."):
        allc = compare_all_dellosa()
    n = len(allc)
    n_clear = int((allc["ΔAIC"] > 10).sum())
    n_mod = int(((allc["ΔAIC"] > 2) & (allc["ΔAIC"] <= 10)).sum())
    n_same = n - n_clear - n_mod
    red = float(np.median(1 - allc["RMSE Cole+L [Ω]"] / allc["RMSE Cole [Ω]"])) * 100
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Cole + L claramente mejor (ΔAIC>10)", f"{n_clear} / {n}")
    c2.metric("Mejora moderada (2<ΔAIC≤10)", f"{n_mod} / {n}")
    c3.metric("Sin diferencia relevante", f"{n_same} / {n}")
    c4.metric("Reducción mediana del RMSE", f"{red:.2f} %")
    fig_h = go.Figure(go.Histogram(x=allc["ΔAIC"], nbinsx=30, marker_color="#38BDF8"))
    fig_h.add_vline(x=2, line_dash="dash", line_color="#FBBF24")
    fig_h.add_vline(x=10, line_dash="dash", line_color="#F43F5E")
    fig_h.update_layout(title="Distribución de ΔAIC = AIC(Cole) − AIC(Cole+L)  (>0 favorece Cole+L)",
                        xaxis_title="ΔAIC", yaxis_title="Registros", template=plot_template,
                        plot_bgcolor=bg_color, paper_bgcolor=bg_color, height=340)
    st.plotly_chart(fig_h, width="stretch")

    st.markdown(f"""
**Conclusión.** Para la gran mayoría de los registros ({n_same} de {n}) agregar la inductancia no mejora el ajuste lo suficiente para justificar un parámetro más.
La reducción mediana del error es de solo {red:.2f} %.
Por eso **el modelo de Cole puro es el que mejor equilibra ajuste y simplicidad** (criterio AIC: a igual calidad, gana el modelo con menos parámetros).

**Por qué la inductancia ayuda tan poco.** La reactancia positiva aparece solo en los últimos puntos (cerca de 65 kHz), donde `|Z|` ronda 70 Ω y la reactancia vale ~2 Ω: pesa poco frente al resto del espectro.
La inductancia estimada es de unos pocos µH y su efecto `ωL` solo se nota ahí.
Los {n_clear + n_mod} registros con mejora tienen en promedio algo más de puntos con reactancia positiva, pero la relación es débil: el error que queda parece venir de otra fuente que esta inductancia no explica.
""")
