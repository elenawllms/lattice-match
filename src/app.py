# LATTICE MATCH PROGRAM
# Falson Lab, Caltech
# Summer 2024
# Elena Williams / elwilliams@hmc.edu / sewillia@caltech.edu

# IMPORTS ----
from pathlib import Path

from dash import Dash, html, dcc, callback, Output, Input, dash_table, State
import pandas as pd
import numpy as np
from plotly import graph_objects as go
import dash_bootstrap_components as dbc

DATA_DIR = Path(__file__).resolve().parent / "assets" / "data"

# GLOBAL VARIABLES ----
ELEMENTS = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
CRYSTAL_SYSTEMS = ["Cubic", "Orthorhombic", "Monoclinic", "Tetragonal", "Hexagonal"]
POINT_GROUPS = ['6/mmm', 'm-3m', '-42m', '6/m', '4/mmm', '6mm', '-43m', '2/m',
       'm-3', 'mmm', 'mm2', '2', 'm', '-4', '422', '222', '432', '23',
       '4/m', '622', '-6m2', '6', '4mm', '-6', '4', '-1', '-3m', '1',
       '3m']

# MOVE THESE TO A SEPARATE FILE

# Resolution of the (a, b) grid used for the Voronoi / mismatch-heatmap solve.
# The figure trace coordinates must be built with this same value.
RESOLUTION = 100

# Rows returned to the matches DataTable. It paginates 10 at a time, so shipping
# the full sorted catalogue (up to 7,201 rows) to the browser is pure waste.
MAX_TABLE_ROWS = 200

# Above this many superlattices the Voronoi solve is skipped.
MAX_VORONOI_SUPERLATTICES = 1000

def mismatch(sub, film):
    return (film - sub) / sub

def costFunction2d(a_mismatch, b_mismatch, mcia, single_double=0.5, large_superlattice=0.5):
    return (np.abs(a_mismatch)**single_double + np.abs(b_mismatch)**single_double)**(1/single_double) * (mcia**(1-large_superlattice))

def costFunction1d(a_mismatch, mcia, large_superlattice=0.5):
    return np.abs(a_mismatch) * (mcia**(1-large_superlattice))

# SUBSTRATE / FILM DATA ----
sublattices_2d = pd.read_csv(DATA_DIR / "sublattices_2d.csv")
sublattices_1d = pd.read_csv(DATA_DIR / "sublattices_1d.csv")

substrates = pd.concat([sublattices_2d["substrate"], sublattices_1d["substrate"]]).unique()

films_2d = pd.read_csv(DATA_DIR / "stable_films_2d.csv")
films_1d = pd.read_csv(DATA_DIR / "stable_films_1d.csv")


# BASIC METHODS ----

def toggle_modal(open_clicks, close_clicks, is_open):
    # Check if the open or close button was clicked
    if open_clicks or close_clicks:
        return not is_open
    return is_open

# HTML ----

header = html.Div([
    html.H1("Lattice Matching to Commercially Available Substrates"),
    html.H3(["Designed by the ", html.A("Falson Lab at Caltech", href="https://epitaxy.caltech.edu/")]),
    html.P([
        "This program is designed to help you find commercially available substrates that match the lattice parameters of your material. ", 
        html.A(" Learn more about the tool.", id="learn-more", n_clicks=0, className="modal-link"),
        dbc.Modal([
            dbc.ModalHeader(dbc.ModalTitle("Information")),
            dbc.ModalBody([
                html.P("This tool finds periodic substrate-film matches for simple crystal geometries (square/rectangular planes such as cubic (001) and triangular planes, as in cubic (111) or hexagonal (0001)). Enter lattice parameters below the table to find the top, matching substrates of those selected at left. Additionally, for small sets of substrates, the website will map the best match at every point in the 2D lattice parameter space (Voronoi diagram) and show the quality of the match (Mismatch Heatmap)."),
                html.H3("Select materials"),
                html.P(["Choose substrates and films in the Select Data tab. Film data comes from the ", html.A("Materials Project", href="https://next-gen.materialsproject.org/"), ", and can be filtered to obtain a class of films of interest. The substrates are all commercially available from Crystec, and empirical structural constants are used. The site precomputes reasonable superlattices which a film might match to; thus, a particular substrate face might match to about 10-15 different film geometries at different angles. The substrates and films are shown in the Display tab, the substrate superlattices in color and the films as black points. Hovering over each substrate, you can find the MCIA and the geometry of each superlattice."]),
                html.H3("Find matches"),
                html.P("Click on any of the film dots to update the table in the center, or manually enter lattice parameters. The table will then provide optimal matches. Change the sliders at left to customize the parameters."),
                html.H3("User-defined parameters"),
                html.P("The lattice mismatch formula is a measure of strain, a simple percentage difference between the matching lengths of substrate and film. However, minimizing the minimal coincident area (MCIA) also improves the quality of the match, so we offer a slider to allow the user to penalize high coincident areas; sliding the scale all the way to 'Yes' fully eliminates this penalty. Similarly, for square/rectangular faces, the user can favor single-axis matches (where one axis is very well aligned) or double-axis matches (where both axes are fairly well aligned). This penalty affects how the program calculates the best match in the table and on the Voronoi/Mismatch plots. In the Display plot, a larger substrate dot corresponds to lower MCIA, and thus a better match, all else equal."),
                html.H3("Voronoi diagram & mismatch heatmap"),
                html.P("The Voronoi diagram displays the best matching substrate at every point in 2D parameter space, and the mismatch heatmap displays the quality of that best match. However, this is an inefficient calculation, so it is only performed on substrate selections with under 1000 superlattices; even so, the tool can still be slow. This many-to-many matching can provide insight beyond the best way to grow a single film; it can visually describe how to grow many films, find the most valuable substrates for a class of new materials, and identify which sets of substrates might be most valuable to grow new films with a wide range of lattice parameters.")
            ]),
            dbc.ModalFooter(
                dbc.Button("Close", id="learn-more-close", className="ml-auto", n_clicks=0)
            ),
        ], id="methodology-modal", is_open=False, size="lg")
    ])
], className="large-container p-3", id="header")


substrates_modal = dbc.Modal([
    dbc.ModalHeader(dbc.ModalTitle("Select Substrates")),
    dbc.ModalBody([
        dcc.Dropdown(substrates, id="substrate-select", multi=True,
        placeholder="Select substrates..."),
    ]),
    dbc.ModalFooter(
        dbc.Button("Update substrates", id="close-substrates-modal", className="ml-auto", n_clicks=0)
    ),
], id="substrates-modal", is_open=False, size="lg")

@callback(
    Output("films-list", "children"),
    Input("database-pull-status", "children")
)
def update_readout(data):
    return data

@callback(
    Output("database-pull-status", "children"),
    Output("film-select", "options"),
    Output("film-select", "value"),
    Output("selected-films-2d", "data"),
    Output("selected-films-1d", "data"),
    Input("pull-materials", "n_clicks"),
    State("must-include-elements", "value"),
    State("can-include-elements", "value"),
    State("exclude-elements", "value"),
    State("num-elements", "value"),
    State("crystal-systems", "value"),
    State("point-groups", "value"),
    State("film-select", "value"),
    State("cubic-faces", "value"),
    State("tetragonal-faces", "value"),
    State("orthorhombic-faces", "value"),
    State("monoclinic-faces", "value"),
    State("hexagonal-faces", "value"),
)
def pull_materials(n_clicks, must_include, can_include, exclude, num_elements, crystal_systems, point_groups, current_films, cubic_faces, tetragonal_faces, orthorhombic_faces, monoclinic_faces, hexagonal_faces):
    if n_clicks == 0: return "No films selected", [], [], [], []

    try:
        num_elements = None if num_elements in (None, "") else int(num_elements)
    except (TypeError, ValueError):
        return "Number of elements must be a whole number", [], [], [], []

    num_elem_condition_2d = num_elements is None or films_2d["num_elements"] == num_elements
    num_elem_condition_1d = num_elements is None or films_1d["num_elements"] == num_elements
    # NOTE: these guards test falsiness, not `is None`. Clearing a multi-select
    # dropdown yields [], not None, and `all(e in [] for e in ...)` is False for
    # every row -- which silently made every search return nothing.
    crystal_sys_condition_2d = not crystal_systems or films_2d["crystal_system"].isin(crystal_systems)
    crystal_sys_condition_1d = not crystal_systems or films_1d["crystal_system"].isin(crystal_systems)
    point_group_condition_2d = not point_groups or films_2d["point_group"].isin(point_groups)
    point_group_condition_1d = not point_groups or films_1d["point_group"].isin(point_groups)
    must_include_condition_2d = not must_include or films_2d["elements"].apply(lambda x: all([e in x for e in must_include]))
    must_include_condition_1d = not must_include or films_1d["elements"].apply(lambda x: all([e in x for e in must_include]))
    can_include_condition_2d = not can_include or films_2d["elements"].apply(lambda x: all([e in can_include for e in x.split(", ")]))
    can_include_condition_1d = not can_include or films_1d["elements"].apply(lambda x: all([e in can_include for e in x.split(", ")]))
    exclude_condition_2d = not exclude or films_2d["elements"].apply(lambda x: not any([e in x for e in exclude]))
    exclude_condition_1d = not exclude or films_1d["elements"].apply(lambda x: not any([e in x for e in exclude]))
    
    cubic_condition_2d = (films_2d["plane"].isin(cubic_faces)) | ((films_2d["crystal_system"] != "Cubic"))
    tetragonal_condition_2d = (films_2d["plane"].isin(tetragonal_faces)) | ((films_2d["crystal_system"] != "Tetragonal"))
    orthorhombic_condition_2d = (films_2d["plane"].isin(orthorhombic_faces)) | ((films_2d["crystal_system"] != "Orthorhombic"))
    monoclinic_condition_2d = (films_2d["plane"].isin(monoclinic_faces)) | ((films_2d["crystal_system"] != "Monoclinic"))
    hexagonal_condition_2d = (films_2d["plane"].isin(hexagonal_faces)) | ((films_2d["crystal_system"] != "Hexagonal"))
    
    cubic_condition_1d = (films_1d["plane"].isin(cubic_faces)) | ((films_1d["crystal_system"] != "Cubic"))
    tetragonal_condition_1d = (films_1d["plane"].isin(tetragonal_faces)) | ((films_1d["crystal_system"] != "Tetragonal"))
    orthorhombic_condition_1d = (films_1d["plane"].isin(orthorhombic_faces)) | ((films_1d["crystal_system"] != "Orthorhombic"))
    monoclinic_condition_1d = (films_1d["plane"].isin(monoclinic_faces)) | ((films_1d["crystal_system"] != "Monoclinic"))
    hexagonal_condition_1d = (films_1d["plane"].isin(hexagonal_faces)) | ((films_1d["crystal_system"] != "Hexagonal"))
    
    
    # print(np.sum(num_elem_condition_2d), np.sum(crystal_sys_condition_2d), np.sum(point_group_condition_2d), np.sum(must_include_condition_2d), np.sum(can_include_condition_2d), np.sum(exclude_condition_2d))

    dff_2d = films_2d[num_elem_condition_2d & crystal_sys_condition_2d & point_group_condition_2d & must_include_condition_2d & can_include_condition_2d & exclude_condition_2d & cubic_condition_2d & tetragonal_condition_2d & orthorhombic_condition_2d & monoclinic_condition_2d & hexagonal_condition_2d]
    dff_1d = films_1d[num_elem_condition_1d & crystal_sys_condition_1d & point_group_condition_1d & must_include_condition_1d & can_include_condition_1d & exclude_condition_1d & cubic_condition_1d & tetragonal_condition_1d & orthorhombic_condition_1d & monoclinic_condition_1d & hexagonal_condition_1d]
    new_films = pd.concat([dff_2d["name"] + " " +  dff_2d["crystal_system"], dff_1d["name"] + " " +  dff_1d["crystal_system"]]).unique()
    
    if len(new_films) == 0:
        return "Search returned no films", [], [], [], []
    elif len(new_films) > 1000:
        return "Too many films selected", [], [], [], []
    else:
        options = list(new_films)
        return f"{len(new_films)} films selected", options, options, dff_2d.to_dict("records"), dff_1d.to_dict("records")

films_modal = dbc.Modal([
    dbc.ModalHeader(dbc.ModalTitle("Select Films")),
    dbc.ModalBody(html.Div([
        html.Div(
            [
                html.P("Filter films by the parameters below; then click 'Pull materials'. The dropdown menu allows further modification. Films with more than 3 elements or with lattice parameters exceeding 16Å are excluded by default to increase performance. You may not select more than 1,000 film planes at a time. Allow a few seconds for this page to load."),
                # Options are populated by `pull_materials`. They were previously
                # seeded with all 78,417 film names -- 5.89 MB of JSON embedded in
                # the initial page payload, mounted into React-Select before the
                # user clicked anything.
                dcc.Dropdown(options=[], multi=True, placeholder="Select films...", id="film-select"),
                html.P("No films selected", id="database-pull-status"),
                dbc.Button("Pull materials", id="pull-materials", className="btn btn-primary", n_clicks=0, style={"margin-bottom": "0.5em"}),
            ], id="pull-info"
        ),
        html.Div([
        html.Fieldset(children=[
                html.Legend("Filters"),
                dbc.Row([
                dbc.Col([
                    html.Div([
                        dbc.Label("Must Include"),
                        dcc.Dropdown(options=ELEMENTS, multi=True, style={"width": "100%"}, id="must-include-elements"),
                    ], className="border-div"),
                    html.Div([
                        dbc.Label("Can Include"),
                        dcc.Dropdown(options=ELEMENTS, multi=True, style={"width": "100%"}, id="can-include-elements"),
                    ], className="border-div"),
                    html.Div([
                        dbc.Label("Exclude Elements"),
                        dcc.Dropdown(options=ELEMENTS, multi=True, style={"width": "100%"}, id="exclude-elements"),
                    ], className="border-div"),
                    
                ]),
                dbc.Col([
                    html.Div([
                        dbc.Label("Number of Elements"),
                        dbc.Input(type="numeric", style={"width": "100%"}, id="num-elements"),
                    ], className="border-div"),
                    html.Div([
                        dbc.Label("Filter Crystal Systems"),
                        dcc.Dropdown(options=CRYSTAL_SYSTEMS, multi=True, id="crystal-systems"),
                    ], className="border-div"),
                    html.Div([
                        dbc.Label("Filter Point Groups"),
                        dcc.Dropdown(options=POINT_GROUPS, multi=True, id="point-groups"),
                    ], className="border-div")
                    # html.Div([
                    #     dbc.Label("Filter by Band Gap"),
                    #     dcc.Dropdown(options=POINT_GROUPS, multi=True),
                    # ], className="border-div"),

                ])
            ])
        ]),
        html.Fieldset(children=[
            html.Legend("Lattice Planes"),
            html.Div([
                dbc.Label("Cubic"),
                dbc.Checklist(
                    options=[
                        {"label": "(100)", "value": "(100)"},
                        {"label": "(110)", "value": "(110)"},
                        {"label": "(111)", "value": "(111)"}
                    ],
                    value=["(100)", "(110)", "(111)"],
                    id="cubic-faces",
                    inline=True,
                ),
            ], className="border-div"),
            html.Div([
                dbc.Label("Tetragonal"),
                dbc.Checklist(
                    options=[
                        {"label": "(001)", "value": "(001)"},
                        {"label": "(100)", "value": "(100)"},
                        {"label": "(101)", "value": "(101)"},
                        {"label": "(110)", "value": "(110)"}
                    ],
                    value=["(001)", "(100)", "(101)", "(110)"],
                    id="tetragonal-faces",
                    inline=True,
                ),
            ], className="border-div"),
            html.Div([
                dbc.Label("Orthorhombic"),
                dbc.Checklist(
                    options=[
                        {"label": "(001)", "value": "(001)"},
                        {"label": "(010)", "value": "(010)"},
                        {"label": "(100)", "value": "(100)"},
                        {"label": "(011)", "value": "(011)"},
                        {"label": "(101)", "value": "(101)"},
                        {"label": "(110)", "value": "(110)"}
                    ],
                    value=["(001)", "(010)", "(100)", "(011)", "(101)", "(110)"],
                    id="orthorhombic-faces",
                    inline=True,
                ),
            ], className="border-div"),
            html.Div([
                dbc.Label("Monoclinic"),
                dbc.Checklist(
                    options=[
                        {"label": "(001)", "value": "(001)"},
                        {"label": "(100)", "value": "(100)"},
                        {"label": "(1-10)", "value": "(1-10)"},
                        {"label": "(110)", "value": "(110)"},
                    ],
                    value=["(001)", "(100)", "(1-10)", "(110)"],
                    id="monoclinic-faces",
                    inline=True,
                ),
            ], className="border-div"),
            html.Div([
                dbc.Label("Hexagonal"),
                dbc.Checklist(
                    options=[
                        {"label": "C (0001)", "value": "(0001)"},
                        {"label": "M (1-100)", "value": "(1-100)"},
                        {"label": "R (1-102)", "value": "(1-102)"},
                        {"label": "A (11-20)", "value": "(11-20)"},
                    ],
                    value=["(0001)", "(1-100)", "(1-102)", "(11-20)"],
                    id="hexagonal-faces",
                    inline=True,
                ),
            ], className="border-div"),
        ])
        ])])
    ),
    dbc.ModalFooter(
        dbc.Button("Update films", id="close-films-modal", className="ml-auto", n_clicks=0)
    ),
], id="films-modal", is_open=False, size="lg")

select_data = html.Div([
    html.H2("Select Data"),
    html.H3("Substrates"),
    html.P(["Currently selected ", html.A("(modify)", n_clicks=0, className="modal-link")], id="open-substrates-modal"), substrates_modal,
    html.Div(id="substrates-list", className="p-2"),
    html.H3("Films"),
    html.P(["Currently selected ", html.A("(modify)", n_clicks=0, className="modal-link")], id="open-films-modal"), films_modal,
    html.Div(id="films-list", className="p-2"),
    html.H3("Single vs. Double-Axis Matching"),
    dcc.Slider(id="single-double-slider", min=0.05, max=1, step=0.05, value=0.5, marks={0.05: "Single", 1: "Double"}),
    html.H3("Allow large superlattices"),
    dcc.Slider(id="large-superlattice-slider", min=0.05, max=1, step=0.05, value=0.5, marks={0.05: "No", 1: "Yes"}),
], className="large-container p-3")

find_matches = html.Div([
    html.H2("Find Matches"),
    dash_table.DataTable(
    style_data = {
        'font-family': 'sans-serif',
        'text-align': 'center'
    },
    style_header = {
        'font-family': 'sans-serif',
        'text-align': 'center'
    },
    id="matches-table", page_size=10),
    html.H3("Manually enter lattice parameters"),
    html.P(["Currently selected: ", html.Span("Manual Entry", id="selected-film-status")]),
    html.Div(
        [
            dbc.RadioItems(
                options=[
                    {"label": "Square/Rectangle", "value": 2},
                    {"label": "Triangle", "value": 1},
                ],
                value=2,
                id="manual-entry-dimension",
                inline=True,
            ),
        ]
    ),
    html.Div([
        "a (Å): ",
        dcc.Input(id='a-input', value=5, type='number', n_blur=0)
    ]),
    html.Div([
        "b (Å): ",
        dcc.Input(id='b-input', value=5, type='number', n_blur=0)
    ], id="manual-entry-b"),

], className="large-container p-3", id="find-matches")

@callback(
    Output("selected-film-status", "children"),
    Input("a-input", "n_blur"),
    Input("b-input", "n_blur"),
)
def update_selected_film_status(a, b):
    return "Manual Entry"

display = dbc.Stack([
    html.H2("Display"),
    html.Div([
        dbc.RadioItems(options=["Scatter", "Voronoi", "Mismatch Heatmap"], inline=True, value="Scatter", id="display-mode"),
    ]),
    html.P("Voronoi diagrams and heatmaps do not display for substrates with more than 1000 superlattices.", id="display-warning"),
    html.H3("Square/Rectangular Lattices"),
    dcc.Graph(id='2d-figure', className="graph"),
    html.H3("Triangular Lattices"),
    dcc.Graph(id='1d-figure', className="graph"),
], className="large-container p-3")

@callback(
    Output("display-warning", "style"),
    Input("2d-superlattices", "data"),
)
def update_display_warning(data):
    # Must use the same threshold and comparison as the guard in
    # update_best_2d_matches. These were previously two separate literals with
    # opposite comparisons, so at exactly 1000 the warning showed while the
    # solve still ran.
    if not data or len(data) <= MAX_VORONOI_SUPERLATTICES:
        return {"display": "none"}
    return {}

# LAYOUT ----
def getLayout():
    return dbc.Container([dbc.Row([
        dbc.Col([
            dbc.Row(dbc.Col(header)),
            dbc.Row([
                dbc.Col(select_data, sm=5),
                dbc.Col(find_matches, sm=7)
            ])
        ], md=7),
        dbc.Col([display], md=5)
    ]),
        dcc.Store(id="selected-substrates"),
        dcc.Store(id="selected-films-2d"),
        dcc.Store(id="selected-films-1d"),
        dcc.Store(id="2d-superlattices"),
        dcc.Store(id="1d-superlattices"),
        dcc.Store(id="best-2d-matches"),
        dcc.Store(id="2d-window", data={"xaxis.range[0]": 4, "xaxis.range[1]": 10, "yaxis.range[0]": 4, "yaxis.range[1]": 10})
    ], 
    fluid=True,
    id='main-container')

# SERVE APP ----
# compress=True needs Flask-Compress installed (see requirements.txt). Dash
# defaults it to False, so callback payloads were being sent uncompressed --
# the figure traces in particular are highly compressible JSON.
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], compress=True)
server = app.server
app.title = "Lattice Matcher | Falson Lab"
app.layout = getLayout()

# CALLBACKS ----

# Modal toggles
app.callback(
    Output("methodology-modal", "is_open"),
    [Input("learn-more", "n_clicks"), Input("learn-more-close", "n_clicks")],
    State("methodology-modal", "is_open")
)(toggle_modal)

app.callback(
    Output("substrates-modal", "is_open"),
    [Input("open-substrates-modal", "n_clicks"), Input("close-substrates-modal", "n_clicks")],
    State("substrates-modal", "is_open")
)(toggle_modal)

app.callback(
    Output("films-modal", "is_open"),
    [Input("open-films-modal", "n_clicks"), Input("close-films-modal", "n_clicks")],
    State("films-modal", "is_open")
)(toggle_modal)


# Hide b input if triangle selected
@app.callback(
    Output("manual-entry-b", "style"),
    [Input("manual-entry-dimension", "value")]
)
def update_manual_entry_b(dimension):
    if dimension == 1:
        return {"display": "none"}
    return {}

# Update selected substrates
@app.callback(
    Output("selected-substrates", "data"),
    Input("substrates-modal", "is_open"),
    State("substrate-select", "value"),
    State("selected-substrates", "data")
)
def update_selected_substrate_data(_, substrates, currentData):
    if substrates is None:
        return currentData
    return substrates

# Update substrate and film readout
@app.callback(
    Output("substrates-list", "children"),
    Input("selected-substrates", "data")
)
def update_readout(data):
    if not data:
        return "None selected"
    if len(data) == len(substrates):
        return "All substrates selected"
    if len(data) > 15:
        return (", ".join(data[:15]) + " +" + str(len(data) - 15) + " more")
    return ", ".join(data)

@app.callback(
    Output("2d-superlattices", "data"),
    Output("1d-superlattices", "data"),
    Input("selected-substrates", "data")
)
def update_superlattices(data):
    if not data:
        return sublattices_2d.to_dict("records"), sublattices_1d.to_dict("records")
    return sublattices_2d[sublattices_2d["substrate"].isin(data)].to_dict("records"), sublattices_1d[sublattices_1d["substrate"].isin(data)].to_dict("records")

# Update table

@app.callback(
    Output("matches-table", "data"),
    Input("2d-superlattices", "data"),
    Input("1d-superlattices", "data"),
    Input("manual-entry-dimension", "value"),
    Input("a-input", "value"),
    Input("b-input", "value"),
    Input("single-double-slider", "value"),
    Input("large-superlattice-slider", "value")
)
def update_table(data2d, data1d, dimension, a, b, single_double, large_superlattice):
    if a == None or b == None: return []
    if data2d == None or data1d == None: return []
    data = data1d if dimension == 1 else data2d
    if dimension == 2:
        rows = [
            {
                "substrate": row["substrate"],
                "dimensions": row["dimensions"], 
                "a_mismatch": mismatch(row["a"], a),
                "b_mismatch": mismatch(row["b"], b),
                "mcia": row["mcia"]
            } for row in data
        ]
        rows.sort(key=lambda x: costFunction2d(x["a_mismatch"], x["b_mismatch"], x["mcia"], single_double, large_superlattice))
    else:
        rows = [
            {
                "substrate": row["substrate"], 
                "dimensions": row["dimensions"], 
                "a_mismatch": mismatch(row["a"], a),
                "mcia": row["mcia"]
            } for row in data
        ]
        rows.sort(key=lambda x: costFunction1d(x["a_mismatch"], x["mcia"], large_superlattice))
    # The table paginates 10 rows at a time; returning all ~7,201 shipped ~1.8 MB
    # of JSON to render 10 of them.
    return rows[:MAX_TABLE_ROWS]


# Hide b-mismatch column if triangle selected
@app.callback(
    Output("matches-table", "columns"),
    Input("manual-entry-dimension", "value")
)
def update_table_columns(dimension):
    if dimension == 1:
        return [
            {"name": "Substrate", "id": "substrate"},
            {"name": "Dimensions", "id": "dimensions"},
            {"name": "a mismatch", "id": "a_mismatch", "type": "numeric", "format": dash_table.FormatTemplate.percentage(2)},
            {"name": "MCIA", "id": "mcia", "type": "numeric", 
            "format": dash_table.Format.Format(precision=2, scheme=dash_table.Format.Scheme.decimal_integer)}
        ]
    return [
        {"name": "Substrate", "id": "substrate"},
        {"name": "Dimensions", "id": "dimensions"},
        {"name": "a mismatch", "id": "a_mismatch", "type": "numeric", "format": dash_table.FormatTemplate.percentage(2)},
        {"name": "b mismatch", "id": "b_mismatch", "type": "numeric", "format": dash_table.FormatTemplate.percentage(2)},
        {"name": "MCIA", "id": "mcia", "type": "numeric", "format": dash_table.Format.Format(precision=2, scheme=dash_table.Format.Scheme.decimal_integer)}
    ]

@app.callback(
    Output("best-2d-matches", "data"),
    Input("2d-superlattices", "data"),
    Input("single-double-slider", "value"),
    Input("large-superlattice-slider", "value"),
    Input("display-mode", "value"),
)
def update_best_2d_matches(superlattices, single_double, large_superlattice, display_mode):
    # This is by far the most expensive callback in the app, and its result is
    # only ever read by the Voronoi and Heatmap branches of update_2d_figure.
    # Without this guard it ran on every substrate change and every slider
    # release while the user sat on the default "Scatter" view, and the result
    # was discarded.
    if display_mode not in ("Voronoi", "Mismatch Heatmap"): return []
    if not superlattices: return []
    if len(superlattices) > MAX_VORONOI_SUPERLATTICES: return []

    limits = [3, 15, 3, 15]

    X, Y = np.meshgrid(np.linspace(limits[0], limits[1], RESOLUTION), np.linspace(limits[2], limits[3], RESOLUTION))
    sp = pd.DataFrame(superlattices)

    # Broadcast the grid against the whole catalogue instead of looping.
    # This previously used np.vectorize -- a Python for-loop over all 10,000
    # grid points, each running ~8 pandas Series ops and computing the cost
    # function twice per point.
    #
    # Done in row chunks: a full (RESOLUTION, RESOLUTION, N) float64 array is
    # ~80 MB at N=1000, and costFunction2d allocates several temporaries of
    # that size. Chunking bounds peak memory to a few MB, which matters on a
    # 512 MB instance.
    sub_a = sp["a"].to_numpy()[None, None, :]
    sub_b = sp["b"].to_numpy()[None, None, :]
    sub_mcia = sp["mcia"].to_numpy()[None, None, :]

    argmin = np.empty((RESOLUTION, RESOLUTION), dtype=np.intp)
    minCost = np.empty((RESOLUTION, RESOLUTION), dtype=float)

    CHUNK = 10
    for lo in range(0, RESOLUTION, CHUNK):
        hi = min(lo + CHUNK, RESOLUTION)
        # Argument order preserved from the original (grid value passed as
        # `sub`). This disagrees with update_table, which passes the substrate
        # as `sub`; reconciling the two is a Stage 1 change, made under test.
        cost = costFunction2d(
            mismatch(X[lo:hi, :, None], sub_a),
            mismatch(Y[lo:hi, :, None], sub_b),
            sub_mcia,
            single_double,
            large_superlattice,
        )
        idx = np.argmin(cost, axis=2)
        argmin[lo:hi] = idx
        minCost[lo:hi] = np.take_along_axis(cost, idx[:, :, None], axis=2).squeeze(2)

    # Positional lookup. The original indexed pandas Series with a positional
    # argmin, which only worked because the frame happened to have a RangeIndex.
    R = sp["R"].to_numpy()[argmin]
    G = sp["G"].to_numpy()[argmin]
    B = sp["B"].to_numpy()[argmin]

    # Calculate x0, y0 and dx, dy
    dx = (limits[1] - limits[0]) / RESOLUTION
    dy = (limits[3] - limits[2]) / RESOLUTION
    x0 = limits[0] + dx / 2
    y0 = limits[2] + dy / 2
    x1 = limits[1] - dx / 2
    y1 = limits[3] - dy / 2
    
    best_match_data = {"R": R, "G": G, "B": B, "minCost": minCost, "x0": x0, "y0": y0, "dx": dx, "dy": dy, "x1": x1, "y1": y1}
    # print(X, Y)
    
    return best_match_data

# Update 2D figure
@app.callback(
    Output("2d-figure", "figure"),
    Input("2d-superlattices", "data"),
    Input("selected-films-2d", "data"),
    # Input("single-double-slider", "value"),
    Input("large-superlattice-slider", "value"),
    Input("best-2d-matches", "data"),
    Input("display-mode", "value")
)
def update_2d_figure(superlattices, films, large_superlattice, best2d, display_mode):
    fig = go.Figure()
    fig.update_layout(
        # autosize=False, 
        #               width=650, height=650,
                    # title="Ideal Substrates for Square/Rectangular Lattices",
                    xaxis_title="a (Å)",
                    yaxis_title="b (Å)",
                    margin=dict(l=70, r=20, t=20, b=70),)
    fig.update_layout(
        plot_bgcolor='white'
    )
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey',
        rangemode="nonnegative"
    )
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey',
        rangemode="nonnegative"
    )

    fig.update_xaxes(range=[4, 10])
    fig.update_yaxes(range=[4, 10])
    
    sp = pd.DataFrame(superlattices)
    sp_labels = sp["substrate"] + " " + sp["dimensions"] + "<br>MCIA: " + sp["mcia"].astype(int).astype(str) + " sq. Å"
    
    fig.add_trace(go.Scattergl(
        x=sp.a,
        y=sp.b,
        mode='markers',
        marker=dict(color = sp.color_column,
                    size = 10 * (100 / sp.mcia)**(1-large_superlattice) if display_mode == "Scatter" else 7,
                   line=dict(width=1, color='Black')
                   ,symbol="square"),
        hovertemplate="%{text}<br>a: %{x:.4f} Å<br>b: %{y:.4f} Å<br><extra></extra>",
        text=sp_labels,
        customdata=[sp.mcia],
        showlegend=False,
        name=""
    ))
    
    flm = pd.DataFrame(films)
    if films:
        flm_labels = flm["formula"] + " " + flm["plane"] + " " + flm["crystal_system"]
        fig.add_trace(go.Scattergl(
            x = flm["a"],
            y = flm["b"],
            mode='markers',
            marker=dict(color='black', size=7, symbol="circle", line=dict(width=2, color='black')),
            hovertemplate="%{text}<br>a: %{x:.4f} Å<br>b: %{y:.4f} Å<br><extra></extra>",
            text=flm_labels,
            showlegend=False,
            name=""
        ))
        
    if not best2d: return fig
    
    if display_mode == "Voronoi":
        r, g, b = best2d["R"], best2d["G"], best2d["B"]
        x0, y0, dx, dy, x1, y1 = best2d["x0"], best2d["y0"], best2d["dx"], best2d["dy"], best2d["x1"], best2d["y1"]
        rgb = np.array([r, g, b]).transpose(2, 1, 0)
        fig.add_trace(go.Image(
            z=rgb*255,
            dx=dx,
            dy=dy,
            x0=x0,
            y0=y0,
            hoverinfo = "skip",
    #         colorscale="Temps",
            opacity=0.4,
            
        ))
        # fig.update_yaxes(range=[y0, y1])
        # fig.update_xaxes(range=[x0, x1])
    elif display_mode == "Mismatch Heatmap":
        fig.add_trace(go.Heatmap(
            z=best2d["minCost"],
            # These must match the z grid. They were hardcoded to 150 while
            # RESOLUTION was 100, so the heatmap was drawn against the wrong
            # coordinates and misregistered with the scatter beneath it.
            x=np.linspace(best2d["x0"], best2d["x1"], RESOLUTION),
            y=np.linspace(best2d["y0"], best2d["y1"], RESOLUTION),
            colorscale="Temps",
            colorbar=dict(title="Cost"),
            hoverinfo="skip",
            showscale=True
        ))
    return fig

@callback(
    # Output('click-data', 'children'),
    Output("selected-film-status", "children", allow_duplicate=True),
    Output('a-input', 'value'),
    Output('b-input', 'value'),
    Output('manual-entry-dimension', 'value'),
    Input('2d-figure', 'clickData'),
    State('a-input', 'value'),
    State('b-input', 'value'),
    State('manual-entry-dimension', 'value'),
    State('selected-film-status', 'children'),
    prevent_initial_call='initial_duplicate'
)
def display_click_data_2d(clickData, current_a, current_b, current_dimension, current_status):
    if clickData is None: return current_status, current_a, current_b, current_dimension
    newStatus = clickData["points"][0].get("text")
    if not newStatus: return current_status, current_a, current_b, current_dimension
    if not clickData["points"][0].get("curveNumber") == 1: return current_status, current_a, current_b, current_dimension
    new_a = clickData["points"][0]["x"]
    new_b = clickData["points"][0]["y"]
    return newStatus, new_a, new_b, 2

@callback(
    Output("selected-film-status", "children", allow_duplicate=True),
    Output('a-input', 'value', allow_duplicate=True),
    Output('manual-entry-dimension', 'value', allow_duplicate=True),
    Input('1d-figure', 'clickData'),
    State('a-input', 'value'),
    State('manual-entry-dimension', 'value'),
    State('selected-film-status', 'children'),
    prevent_initial_call='initial_duplicate'
)
def display_click_data_1d(clickData, current_a, current_dimension, current_status):
    if clickData is None: return current_status, current_a, current_dimension
    newStatus = clickData["points"][0].get("text")
    if not newStatus: return current_status, current_a, current_dimension
    if not clickData["points"][0].get("curveNumber") == 1: return current_status, current_a, current_dimension
    return newStatus, clickData["points"][0]["x"], 1


@app.callback(
    Output("1d-figure", "figure"),
    Input("1d-superlattices", "data"),
    Input("selected-films-1d", "data"),
    Input("large-superlattice-slider", "value"),
    Input("display-mode", "value")
)
def update_1d_figure(substrates, films, large_superlattice, display_mode):
    fig = go.Figure()
    fig.update_layout(
        height=200,
        xaxis_title="a (Å)",
        plot_bgcolor='white',
                    margin=dict(l=70, r=20, t=20, b=70),)
    fig.update_yaxes(
        fixedrange=True,
        range=[-1, 1],
        zeroline=True,
        zerolinewidth=2,
        zerolinecolor="black",
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey',
    )
    fig.update_xaxes(
        range=[4, 10],
    )
    
    
    if not substrates: return fig
    
    
    sp = pd.DataFrame(substrates)
    flm = pd.DataFrame(films)
        
    fig.add_trace(go.Scattergl(
        x=sp["a"],
        y=np.zeros(len(sp)),
        mode='markers',
        marker=dict(color = sp.color_column,
                    size = 10 * (100 / sp["mcia"])**(1-large_superlattice) if display_mode == "Scatter" else 7,
                   line=dict(width=1, color='Black')
                   ,symbol="square"),
        hovertemplate="%{text}<br>a: %{x:.4f} Å<br><extra></extra>",
        text=sp["substrate"] + " " + sp["dimensions"] + "<br>MCIA: " + sp["mcia"].astype(int).astype(str) + " sq. Å",
        showlegend=False,
        name=""
    ))
    
    if not films: return fig
    fig.add_trace(go.Scattergl(
        x = flm["a"],
        y = np.zeros(len(flm)),
        mode='markers',
        marker=dict(color='black', size=7, symbol="circle", line=dict(width=2, color='black')),
        hovertemplate="%{text}<br>a: %{x:.4f} Å<br><extra></extra>",
        text=flm["formula"] + " " + flm["plane"] + " " + flm["crystal_system"],
        showlegend=False,
        name=""
    ))
    return fig


# Only runs for local development. Under gunicorn this module is imported to
# resolve `app:server`, and without the __main__ guard that import called
# run_server(), which blocks forever in the Werkzeug development server -- so
# the worker never reached gunicorn's accept loop.
if __name__ == "__main__":
    app.run(debug=True)



