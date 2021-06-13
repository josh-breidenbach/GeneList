import plotly.express as px
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_table
import dash_bootstrap_components as dbc
import re
import numpy as np
import pandas as pd
import dash
from dash import no_update

columnstyle = [
    {
        'if': {'column_id': 'tax_id'},
        'padding': '5px',
        'minWidth': '100px', 'width': '100px', 'maxWidth': '100px',
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
        'textAlign': 'left',
    },
    {
        'if': {'column_id': 'GeneID'},
        'padding': '5px',
        'minWidth': '100px', 'width': '100px', 'maxWidth': '100px',
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
        'textAlign': 'left',
    },
    {
        'if': {'column_id': 'Symbol'},
        'padding': '5px',
        'minWidth': '100px', 'width': '100px', 'maxWidth': '100px',
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
        'textAlign': 'left',
    },
    {
        'if': {'column_id': 'Synonyms'},
        'padding': '5px',
        'minWidth': '200px', 'width': '200px', 'maxWidth': '200px',
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
    },
    {
        'if': {'column_id': 'description'},
        'padding': '5px',
        'minWidth': '300px', 'width': '300px', 'maxWidth': '300px',
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
    },
    {
        'if': {'column_id': 'type_of_gene'},
        'padding': '5px',
        'minWidth': '200px', 'width': '200px', 'maxWidth': '200px',
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
    },
    {
        'if': {'column_id': 'Symbol_from_nomenclature_authority'},
        'padding': '5px',
        'minWidth': '300px', 'width': '300px', 'maxWidth': '300px',
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
    },
    {
        'if': {'column_id': 'Full_name_from_nomenclature_authority'},
        'padding': '5px',
        'minWidth': '300px', 'width': '300px', 'maxWidth': '300px',
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
    },
    {
        'if': {'column_id': 'Other_designations'},
        'padding': '5px',
        'minWidth': '250px', 'width': '250px', 'maxWidth': '250px',
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
    },
]


class GeneList:
    def __init__(self):
        # Read in original data

        # Original_Data = pd.read_csv(r"C:\Users\Josh\Google Drive\Projects\GeneLists\Homo_sapiens.gene_info.csv")
        # Original_Data = pd.read_csv(r"E:\Google Drive\Projects\GeneLists\Homo_sapiens.gene_info.csv")
        Original_Data = pd.read_csv(
            r"C:\Users\Josh\Documents\PythonProjects\PythonAnywhere\GeneListProject\assets\Homo_sapiens.gene_info.csv")

        columns_to_keep = Original_Data[[
            '#tax_id',
            'GeneID',
            'Symbol',
            'Synonyms',
            'description',
            'type_of_gene',
            'Symbol_from_nomenclature_authority',
            'Full_name_from_nomenclature_authority',
            'Other_designations'
        ]]
        self.slim_df = columns_to_keep.copy()
        self.slim_df['GeneID'] = self.slim_df['GeneID'].apply(str)
        self.matches_df = pd.DataFrame()
        self.alias_dict = {}
        self.no_match = []

    def search(self, search):
        self.search = search
        self.matchNumber = 0
        self.aliasNumber = 0
        self.noMatchNumber = 0
        print("Searched terms:")
        self.search = list(
            dict.fromkeys(self.search))  # Utilizes a transformation to a dictionary and back to remove duplicates
        for term in self.search:

            upper = term.upper()
            # print(upper)

            # print (term)
            if term in list(self.slim_df.Symbol):
                term_found_df = self.slim_df.loc[self.slim_df['Symbol'] == term]
                self.matches_df = pd.concat([self.matches_df, term_found_df]).drop_duplicates().reset_index(drop=True)
                self.matchNumber += 1
            elif term in list(self.slim_df.GeneID):
                term_found_df = self.slim_df.loc[self.slim_df['GeneID'] == term]
                self.matches_df = pd.concat([self.matches_df, term_found_df]).drop_duplicates().reset_index(drop=True)
                self.matchNumber += 1
            elif True in ((pd.Series(list(self.slim_df.Synonyms))).str.contains(term).values):
                alias_df = self.slim_df[self.slim_df['Synonyms'].str.contains(term)]
                self.alias_dict[term] = alias_df  # Dictionaries don't allow redundant key entries
                self.aliasNumber += 1
            elif upper in list(self.slim_df.Symbol):
                alias_df = self.slim_df.loc[self.slim_df['Symbol'] == upper]
                self.alias_dict[term] = alias_df  # Dictionaries don't allow redundant key entries
                self.aliasNumber += 1
            elif True in ((pd.Series(list(self.slim_df.Synonyms))).str.contains(upper).values):
                alias_df = self.slim_df[self.slim_df['Synonyms'].str.contains(upper)]
                self.alias_dict[term] = alias_df  # Dictionaries don't allow redundant key entries
                self.aliasNumber += 1
            else:
                if term not in self.no_match:
                    self.no_match.append(term)
                    self.noMatchNumber += 1

        if self.matchNumber > 0:
            print('\n', self.matchNumber, "Genes matched exactly")
            # print (self.matches_df)
        if self.aliasNumber > 0:
            print('\n', self.aliasNumber, "Genes with suggested matches")
            # print (self.alias_dict)
        if self.noMatchNumber > 0:
            print('\n', self.noMatchNumber, "Genes with No Match")
            # print (self.no_match)

    def clear(self):
        self.matches_df = pd.DataFrame()
        self.alias_dict = {}
        self.no_match = []


# Start empty variables
chosen_row = ''
inputlist = []

Original_df = pd.read_csv(r"C:\Users\Josh\Documents\PythonProjects\PythonAnywhere\GeneListProject\assets\Homo_sapiens.gene_info.csv")

bs = 'https://cdn.jsdelivr.net/npm/bootstrap@4.6.0/dist/css/bootstrap.min.css'

app = dash.Dash(__name__, external_stylesheets=[bs],
                meta_tags=[{'name': 'viewport', 'content': 'width=device-width, initial-scale=1.0'}]
                )

app.layout = dbc.Container([
    dbc.Row(dbc.Col(html.H1("Gene List",
                            className='text-primary, mb-4'),
                    width={'size': True, 'offset': 0},
                    ),
            ),

    dbc.Row([
        dbc.Col(dbc.Card([
            dbc.CardHeader(html.H2("Search")),
            dbc.CardBody(
                [
                    html.H6("Input gene (symbol or ID) or list of genes:", className="card-subtitle, mb-2"),
                    dcc.Textarea(
                        id="userlist",
                        className="mb-2",
                        placeholder="Gene or List",
                    ),
                    dbc.Button('Search', id='submit-search', n_clicks=0, color="primary", className="mb-1"),
                    html.Div(id='textOutput', style={'whiteSpace': 'pre-line'}),
                ]
            ),
        ],
            className="mb-2",
            style={"width": "17rem"},
            color="light",
        ),
            #  width={'size':2, 'offset':0},
        ),

        dbc.Col(dbc.Card([
            dbc.CardHeader(html.H3("Suggested Matches")),
            dbc.CardBody(
                [
                    html.H6("Some genes may have suggested matches.", className="card-subtitle"),
                    html.H6("Select a gene to see suggestions.", className="card-subtitle, mb-2"),
                    dcc.Dropdown(
                        id='aliasfound-dropdown',
                        value='',
                        placeholder="Select a gene",
                        style={"width": "10rem"},
                        className="mb-2",
                    ),
                    html.Div(id='dd-output-container'),
                    html.Div([
                        dash_table.DataTable(
                            id='aliasDataTable',
                            row_selectable="single",
                            selected_rows=[],
                            page_action="native",
                            page_current=0,
                            page_size=3,
                            style_cell_conditional=columnstyle,
                            style_header={
                                'backgroundColor': 'white',
                                'fontWeight': 'bold'
                            },
                            style_table={
                                'overflowX': 'auto'
                            },
                            style_as_list_view=True,
                            derived_virtual_selected_rows=[]
                        ),
                        html.Div(id='aliasDataTable-container'),
                    ], className="mb-2"),
                    dbc.Button('Accept Suggestion', id='submit-val', n_clicks=0, style=dict(display='none'),
                               color="primary", className="mr-1"),
                ]
            ),
            # className="mb-2",
            # style={"width": "18rem"},
        ],
            color="light",
            className="mb-2",
        ),
            width={'size': 9},
        )
    ], justify="start"),

    dbc.Row(dbc.Col(dbc.Card([
        dbc.CardHeader(html.H2("Final List")),
        dbc.CardBody(
            [
                # html.H6("Some genes may have suggested matches", className="card-subtitle, mb-2"),
                html.Div([
                    dash_table.DataTable(
                        id='Final_List_datatable',
                        columns=[
                            {"name": "#tax_id", "id": "#tax_id", "hideable": "last"},
                            {"name": "GeneID", "id": "GeneID", "hideable": "last"},
                            {"name": "Symbol", "id": "Symbol", "hideable": "last"},
                            {"name": "Synonyms", "id": "Synonyms", "hideable": "last"},
                            {"name": "description", "id": "description", "hideable": "last"},
                            {"name": "type_of_gene", "id": "type_of_gene", "hideable": "last"},
                            {"name": "Symbol_from_nomenclature_authority", "id": "Symbol_from_nomenclature_authority",
                             "hideable": "last"},
                            {"name": "Full_name_from_nomenclature_authority",
                             "id": "Full_name_from_nomenclature_authority", "hideable": "last"},
                            {"name": "Other_designations", "id": "Other_designations", "hideable": "last"},
                        ],
                        filter_action="native",
                        sort_action="native",
                        sort_mode="multi",
                        page_action="native",
                        page_current=0,
                        page_size=10,
                        style_table={
                            'overflowX': 'auto'
                        },
                        style_as_list_view=True,
                        style_cell_conditional=columnstyle,
                        style_header={
                            'backgroundColor': 'white',
                            'fontWeight': 'bold'
                        },
                        export_format='xlsx',
                        export_headers='display',
                    ),
                    html.Div(id='Final_List_datatable-container'),
                ]),
            ]
        ),
    ],
        color="light",
        className="mb-2",
    ),
        # style={"width": "18rem"},
    )
    ),

], fluid=True)


# Define callback to show alias data table when a gene is selected in the alias dropdown

@app.callback(
    [Output(component_id='aliasDataTable', component_property='data'),
     Output(component_id='aliasDataTable', component_property='columns'),
     Output(component_id='submit-val', component_property='style')],
    [Input(component_id="aliasfound-dropdown", component_property="value")],
    prevent_initial_call=True
)
# Respond to selections made in the alias dropdown
def aliasDropdown_selected(aliasDropdown_value):
    aliasdrop = GeneList()
    selectedAlias = str(aliasDropdown_value)
    Listof_selectedAlias = [selectedAlias]
    aliasdrop.search(Listof_selectedAlias)
    selectedAlias_df = aliasdrop.alias_dict[selectedAlias]
    aliasDataTableColumns = [{"name": i, "id": i} for i in selectedAlias_df.columns]
    styleout = dict()

    return selectedAlias_df.to_dict('records'), aliasDataTableColumns, styleout


# Define callback to add the selected suggestion to the final list when the "Accept Suggestion" button is pushed.

@app.callback(
    [Output(component_id="userlist", component_property="value"),
     Output(component_id="submit-search", component_property="n_clicks")],
    [Input(component_id="submit-val", component_property="n_clicks")],
    [State(component_id="aliasfound-dropdown", component_property="value"),
     State(component_id="aliasDataTable", component_property="derived_virtual_selected_rows")],
    prevent_initial_call=True
)
def aliasRowSelected(nclick, aliasDropdown_value2, chosen_row):
    aliasselect = GeneList()
    selectedAlias = str(aliasDropdown_value2)
    Listof_selectedAlias = [selectedAlias]
    aliasselect.search(Listof_selectedAlias)
    selectedAlias_df = aliasselect.alias_dict[selectedAlias]

    selected_suggestion = selectedAlias_df.iloc[chosen_row, 2]
    x = []
    x = list(selected_suggestion)
    x = x[0]
    print(x)
    if x != []:
        return x, 0
    else:
        return dash.no_update, dash.no_update


# Define callback to search the gene or gene list input by user.

@app.callback(
    [Output(component_id='textOutput', component_property='children'),
     Output(component_id='Final_List_datatable', component_property='data'),
     Output(component_id='aliasfound-dropdown', component_property='options')],
    [Input(component_id="submit-search", component_property="n_clicks")],
    [State(component_id="userlist", component_property="value")],
    prevent_initial_call=True  # Prevents the initial call so that there are no issues with test.matches_df being empty
)
# Respond to changes made in the userlist input
# Format this userlist input
def update_list(nclick, userlist_value):
    print(userlist_value)
    # print(nclick)
    if (userlist_value != None) and (userlist_value != ''):
        if ' ' in userlist_value:
            userlist_value_split = re.split('[\s,]', userlist_value)
        elif ',' in userlist_value:
            userlist_value_split = re.split('[\s,]', userlist_value)
        else:
            userlist_value_split = userlist_value.split()
        userlist_value_split = [i for i in userlist_value_split if i != ',']
        userlist_value_split = [i for i in userlist_value_split if i != '']

        test = GeneList()

        inputlist.extend(userlist_value_split)  # Adds user input to search list (inputlist)
        test.search(inputlist)  # This inputlist is fed into the search function as the "search" list

        textoutput = '{} genes matched exactly \n {} genes matched to a synonym \n {} genes did not match'.format(
            test.matchNumber, test.aliasNumber, test.noMatchNumber)

        # aliasKeyList = list(test.alias_dict.keys()) # Creates list of keys in alias dictionary
        aliasFoundOptionsOutput = []
        for i in test.alias_dict:
            aliasFoundOptionsOutput.append({'label': i, 'value': i})
        if not aliasFoundOptionsOutput:  # If aliasFoundOptionsOutput is empty
            aliasFoundOptionsOutput = [{'label': "None", 'value': ''}]
        ##generate_table(test.matches_df)
        return textoutput, test.matches_df.to_dict(
            'records'), aliasFoundOptionsOutput  # Equal number of returns for each "Output", in order

    # textoutput is the number of matched genes etc.
    # test.matches_df.to_dict('records') is the data for the dashtable
    # aliasFoundOptionsOutput provides the options in the aliasfound-dropdown

    # Run app and display result inline in the notebook


if __name__ == '__main__':
    app.run_server(debug=False, port=8070)
