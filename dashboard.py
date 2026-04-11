import pandas as pd
import numpy as np
from dash import Dash, dcc, html, Input, Output
import plotly.graph_objs as go

# ==============================
# LOAD
# ==============================
CSV_PATH = "sim_export_bodies.csv"

df = pd.read_csv(CSV_PATH, comment="#")
df = df.sort_values(["t", "body_id"])

# ==============================
# PRECOMPUTE
# ==============================
# Energia
energy = df.groupby("t").agg({
    "ke": "sum",
    "pe": "sum"
}).reset_index()

energy["E"] = energy["ke"] + energy["pe"]
E0 = energy["E"].iloc[0]
energy["dE"] = (energy["E"] - E0) / abs(E0)

# Momento angular
df["Lz"] = df["x"] * df["vy"] - df["y"] * df["vx"]
Lz = df.groupby("t")["Lz"].sum().reset_index()

# Drift orbital
orbital = df[["t", "body_id", "orb_a"]].copy()
a0 = orbital.groupby("body_id")["orb_a"].first()
orbital["a0"] = orbital["body_id"].map(a0)
orbital["da"] = (orbital["orb_a"] - orbital["a0"]) / orbital["a0"]

body_ids = df["body_id"].unique()

# ==============================
# APP
# ==============================
app = Dash(__name__)

app.layout = html.Div([
    html.H2("Gravity Simulation Dashboard"),

    # seleção de corpos
    dcc.Dropdown(
        id="body-select",
        options=[{"label": f"Body {b}", "value": int(b)} for b in body_ids],
        value=list(body_ids[:5]),
        multi=True
    ),

    # slider temporal
    dcc.Slider(
        id="time-slider",
        min=df["t"].min(),
        max=df["t"].max(),
        step=None,
        value=df["t"].min(),
        marks={float(t): f"{t:.2f}" for t in df["t"].unique()[::len(df["t"].unique())//10 or 1]}
    ),

    dcc.Graph(id="energy-plot"),
    dcc.Graph(id="drift-plot"),
    dcc.Graph(id="lz-plot"),
    dcc.Graph(id="orbital-plot"),
    dcc.Graph(id="position-plot"),
])

# ==============================
# CALLBACK
# ==============================
@app.callback(
    Output("energy-plot", "figure"),
    Output("drift-plot", "figure"),
    Output("lz-plot", "figure"),
    Output("orbital-plot", "figure"),
    Output("position-plot", "figure"),
    Input("body-select", "value"),
    Input("time-slider", "value"),
)
def update(selected_bodies, t_selected):

    # --------------------------
    # ENERGIA
    # --------------------------
    energy_fig = go.Figure()
    energy_fig.add_trace(go.Scatter(x=energy["t"], y=energy["E"], name="E"))
    energy_fig.add_trace(go.Scatter(x=energy["t"], y=energy["ke"], name="K"))
    energy_fig.add_trace(go.Scatter(x=energy["t"], y=energy["pe"], name="U"))

    # --------------------------
    # DRIFT
    # --------------------------
    drift_fig = go.Figure()
    drift_fig.add_trace(go.Scatter(
        x=energy["t"],
        y=np.abs(energy["dE"]),
        name="|ΔE/E₀|"
    ))

    drift_fig.update_yaxes(type="log")

    # --------------------------
    # Lz
    # --------------------------
    lz_fig = go.Figure()
    lz_fig.add_trace(go.Scatter(x=Lz["t"], y=Lz["Lz"], name="Lz"))

    # --------------------------
    # ORBITAL DRIFT
    # --------------------------
    orbital_fig = go.Figure()

    for bid in selected_bodies:
        b = orbital[orbital["body_id"] == bid]
        orbital_fig.add_trace(go.Scatter(
            x=b["t"], y=b["da"], name=f"Body {bid}"
        ))

    # --------------------------
    # POSIÇÃO (snapshot)
    # --------------------------
    frame = df[df["t"] == t_selected]

    pos_fig = go.Figure()

    pos_fig.add_trace(go.Scatter(
        x=frame["x"],
        y=frame["y"],
        mode="markers",
        marker=dict(size=4),
        name="Bodies"
    ))

    pos_fig.update_layout(
        title=f"Positions at t={t_selected:.3f}",
        xaxis_title="x",
        yaxis_title="y"
    )

    return energy_fig, drift_fig, lz_fig, orbital_fig, pos_fig


# ==============================
# RUN
# ==============================
if __name__ == "__main__":
    app.run(debug=True)