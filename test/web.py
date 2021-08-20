import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_table
import plotly.graph_objs as go
from datetime import datetime as dt

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
from collections import OrderedDict

COLORS = [
    {'background': 'white', 'text': '#3CCC3C'},
    {'background': 'white', 'text': '#3CCC3C'},
    {'background': 'white', 'text': '#3CCC3C'},
    {'background': 'white', 'text': '#4FF4C4C'}
]

field_color = '#5F5F5F5'
field_pad = 10
margin = 10
headline_color = '#F1A44F'


def is_numeric(value) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def generate_table1(dataframe, max_rows=10) -> html:
    """
    :return: Return html converted from pandas dataframe
    """

    return html.Table([html.Tr([html.Th(col) for col in dataframe.columns])] + [
        html.Tr([html.Td(dataframe.iloc[i][col]) for col in dataframe.columns]) for i in
        range(min(len(dataframe), max_rows))])


def cell_style(value, min_value, max_value) -> dict:
    style = {}
    if is_numeric(value):
        if value <= 0.25:
            style = {
                'backgroundColor': COLORS[0]['background'],
                'color': COLORS[0]['text']
            }
        elif value <= 0.5:
            style = {
                'backgroundColor': COLORS[1]['background'],
                'color': COLORS[1]['text']
            }
        elif value <= 1:
            style = {
                'backgroundColor': COLORS[2]['background'],
                'color': COLORS[2]['text']
            }
        elif value > 1:
            style = {
                'backgroundColor': COLORS[3]['background'],
                'color': COLORS[3]['text']
            }
        return style


def generate_table(dataframe, max_rows=100):
    max_value = dataframe.max(numeric_only=True).max()
    min_value = dataframe.min(numeric_only=True).min()
    rows = []
    for i in range(min(len(dataframe), max_rows)):
        row = []
        for col in dataframe.columns:
            value = dataframe.iloc[i]['UR[-]']
            style = cell_style(value, min_value, max_value)
            row.append(html.Td(dataframe.iloc[i][col], style=style))
        rows.append(html.Tr(row))

    return html.Table([html.Tr([html.Th(col) for col in dataframe.columns])] + rows,
                      style={'color': 'black', 'background': '#F8F8F8', 'fontSize': 12, 'height': 2, 'line-height': 2})


app = dash.Dash('Section Analysis')

app.layout = html.Div([
    html.H1('Capacity surface visualization', style={'color': headline_color, 'fontSize': 30}),

    html.Div(
        [
            html.H2("Row with columns"),
            dbc.Row(dbc.Col(html.Div("A single column"))),
            dbc.Row(
                [
                    dbc.Col(html.Div(
                        dbc.InputGroup(
                            [dbc.InputGroupAddon("Ec", addon_type="append"),
                             dbc.Input(placeholder="200000"),
                             dbc.InputGroupAddon("MPa", addon_type="append"), ]
                        ),
                    ), ),
                    dbc.Col(html.Div(
                        dbc.InputGroup(
                            [dbc.InputGroupAddon("Ec", addon_type="append"),
                             dbc.Input(placeholder="200000"),
                             dbc.InputGroupAddon("MPa", addon_type="append"), ]
                        ),
                    ), ),
                    dbc.Col(html.Div(
                        dbc.InputGroup(
                            [dbc.InputGroupAddon("Ec", addon_type="append"),
                             dbc.Input(placeholder="200000"),
                             dbc.InputGroupAddon("MPa", addon_type="append"), ]
                        ),
                    ), ),
                ]
            ),
        ]
    ),

    html.Div(
        [
            dbc.Button('Primary', color='primary', outline=True, className='mr-1')
        ]
    ),

    html.Div(
        className='row',
        children=[
            html.Div(
                className='col s2 m2 l2',
                children=[
                    html.Div([
                        # Div for specifying the cross section vertices
                        html.Label('Cross section vertices', style={
                            'color': headline_color, 'fontSize': 16}),
                        dash_table.DataTable(
                            id='section-vertices',
                            columns=(
                                [{'id': p, 'name': p}
                                 for p in ['x[mm]', 'y[mm]']]
                            ),
                            data=[{'x[mm]': '200', 'y[mm]': '200'},
                                  {'x[mm]': '-200', 'y[mm]': '200'},
                                  {'x[mm]': '-200', 'y[mm]': '-200'},
                                  {'x[mm]': '200', 'y[mm]': '-200'},
                                  {'x[mm]': '', 'y[mm]': ''},
                                  {'x[mm]': '', 'y[mm]': ''},
                                  {'x[mm]': '', 'y[mm]': ''},
                                  {'x[mm]': '', 'y[mm]': ''},
                                  {'x[mm]': '', 'y[mm]': ''},
                                  {'x[mm]': '', 'y[mm]': ''}],
                            editable=True,
                            style_table={
                                # 'maxHeight': '500',
                                # 'overflowY': 'scroll'
                            },
                        ),
                    ], ),
                ], style={'backgroundColor': 'white', 'width': '47%', 'padding': field_pad, 'border-radius': 5,
                          'margin': 0, 'float': 'left'},
            ),
            html.Div(
                className='col s2 m2 l2',
                children=[
                    html.Div([
                        # Div for specifying rebar locations
                        html.Label('Rebar coordinates', style={
                            'color': headline_color, 'fontSize': 16}),
                        dash_table.DataTable(
                            id='rebar-locations',
                            columns=(
                                [{'id': p, 'name': p} for p in ['xs[mm]', 'ys[mm]']]
                            ),
                            data=[{'xs[mm]': '140', 'ys[mm]': '140'},
                                  {'xs[mm]': '-140', 'ys[mm]': '140'},
                                  {'xs[mm]': '-140', 'ys[mm]': '-140'},
                                  {'xs[mm]': '140', 'ys[mm]': '-140'},
                                  {'xs[mm]': '140', 'ys[mm]': '0'},
                                  {'xs[mm]': '-140', 'ys[mm]': '0'},
                                  {'xs[mm]': '0', 'ys[mm]': '140'},
                                  {'xs[mm]': '0', 'ys[mm]': '-140'},
                                  {'xs[mm]': '', 'ys[mm]': ''},
                                  {'xs[mm]': '', 'ys[mm]': ''},
                                  ],
                            editable=True
                        ),
                    ], ),
                ], style={'backgroundColor': 'white', 'width': '47%', 'padding': field_pad, 'border-radius': 5,
                          'margin': 0, 'float': 'right'}
            ),
        ], style={'backgroundColor': '#FFFFFF', 'padding': 0, 'border-radius': 5},
    ),

    # Div for alerts
    html.Div([dbc.Alert('Some rebars are located outside the cross section!', color='primary')]),
    html.Div([dbc.Alert('Polygon is self-intersecting!', color='danger')]),

    # Div for graphs
    html.Div(
        className='row',
        children=[
            html.Div(
                className='col s2 m2 l2',
                children=[
                    html.Div([
                        dcc.Graph(
                            id='section-plot',
                        ),
                    ], ),
                ], style={'backgroundColor': field_color, 'width': '47%', 'padding': field_pad, 'border-radius': 5,
                          'margin': margin, 'float': 'left'},
            ),
            html.Div(
                className='col s2 m2 l2',
                children=[
                    html.Div([
                        dcc.Graph(
                            id='capacity-surface',
                        ),
                    ], ),
                ], style={'backgroundColor': field_color, 'width': '47%', 'padding': field_pad, 'border-radius': 5,
                          'margin': margin, 'float': 'right'}
            ),
        ], style={'backgroundColor': '#FFFFFF', 'padding': 0, 'border-radius': 5},
    ),

    # Div for load combinations and results
    html.Div(
        className='row',
        children=[
            html.Div(
                className='col s2 m2 l2',
                children=[
                    html.Label('Load combinations', style={'color': headline_color, 'fontSize': 16}),
                    dash_table.DataTable(
                        id='load-combs',
                        columns=(
                            [{'id': p, 'name': p}
                             for p in ['P[kN]', 'Mx[kNm]', 'My[kNm]']]
                        ),
                        data=[{'P[kN]': '200', 'Mx[kNm]': '200', 'My[kNm]': '200'},
                              {'P[kN]': '-200', 'Mx[kNm]': '250', 'My[kNm]': '150'},
                              {'P[kN]': '-500', 'Mx[kNm]': '-300', 'My[kNm]': '75'},
                              {'P[kN]': '300', 'Mx[kNm]': '-350', 'My[kNm]': '200'}],
                        editable=True
                    ),
                ], style={'backgroundColor': 'white', 'width': '47%', 'padding': field_pad, 'border-radius': 5,
                          'margin': 0, 'float': 'left'},
            ),
            html.Div(
                className='col s2 m2 l2',
                children=[
                    html.Label('Result table', style={'color': headline_color, 'fontSize': 16}),
                    html.Div(
                        # Div for result table
                        id='result-table',
                        # Table is returned here from a callback
                    ),
                ], style={'backgroundColor': 'white', 'width': '47%', 'padding': field_pad, 'border-radius': 5,
                          'margin': 0, 'float': 'right'}
            ),
        ], style={'backgroundColor': '#FFFFFF', 'padding': 0, 'border-radius': 5},
    ),

    # Hidden div for storing the computed capacity surface so it can be shared by many callbacks
    # without computing it over and over each time
    html.Div(id='capacity-surface-results', style={'display': 'none'})

], className='container', style={'width': '95%'})


@app.callback(
    Output('capacity-surface-results', 'children'),
    [Input('section-vertices', 'data'),
     Input('secetion-vertices', 'columns'),
     Input('rebar-locations', 'data'),
     Input('rebar-locations', 'columns'), ])
def calc_surf_and_store(xy, sv_col, xsys, rebar_col):
    df_sv = pd.DataFrame(xy, columns=[c['name'] for c in sv_col])
    df_rebars = pd.DataFrame(xsys, columns=[c['name'] for c in rebar_col])

    x = [float(c) for c in list(df_sv['x[mm]']) if c]
    y = [float(c) for c in list(df_sv['y[mm]']) if c]

    xr = [float(c) for c in list(df_rebars['xs[mm]']) if c]
    yr = [float(c) for c in list(df_rebars['ys[mm]']) if c]

    eps_cu = 0.0035
    fck = 25
    gamma_c = 1.0
    Es = 200 * 10 ** 3
    fcd = fck / gamma_c
    fyk = 500
    gamma_s = 1.0
    fyd = fyk / gamma_s
    As = np.pi() * 25 ** 2 / 4

    
