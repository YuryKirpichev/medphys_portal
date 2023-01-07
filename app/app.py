import io
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import os
import base64
from flask import Flask

import dash
from dash.dependencies import Input, Output, State
from dash import dcc
from dash import html

import dash_bootstrap_components as dbc

from dash import dash_table

import pydicom as dcm
from math import sqrt

from pylinac import PicketFence, Starshot, LasVegas
from pylinac.picketfence import  MLC
from plotly.tools import mpl_to_plotly
from datetime import datetime

import numpy as np
import plotly.graph_objects as go
import matplotlib.cm as cm
from matplotlib.colors import Normalize

from datetime import datetime

import os
import glob

trs_ptw_31016 = dict({
    10:   1,
    8:    1, 
    6:    1,
    4:    1.001,
    3:    1.001,
    2.5:  1.001,
    2:    1.004,
    1.5:  1.013,
    1.2:  1.025,
    1.0:  1.039
    })

from scipy.interpolate import interp1d
import numpy as np

x = list(trs_ptw_31016.keys())
y = list(trs_ptw_31016.values())
f_ptw_31016 = interp1d(x, y, 'cubic', fill_value='extrapolate')


def clean_temp_folders():

    files_dont_del = ['.gitkeep']
    os.listdir('temp')
    for f in os.listdir('temp'):
        if f not in files_dont_del:
            os.remove('temp/'+f)
        
    os.listdir('assets')
    for f in os.listdir('assets'):
        if f not in files_dont_del:
            os.remove('assets/'+f)


def plot_plotly_from_dicom(image):

    # image dimensions (pixels)
    n1 = 1280 # height
    n2 = 1280 # width
    # Generate an image starting from a numerical function
    x, y = np.mgrid[-3:3:n1*1j, -3:3:n2*1j]
    f = np.array(image)

    # use matplotlib to assign a colormap to the computed values
    norm = Normalize(f.min(), f.max())
    norm_f = norm(f)
    # generate the image:
    # img has a shape of (150, 100, 4).
    # The four channels are R, G, B, alpha
    # All values will be between 0 and 1
    img = cm.viridis(norm_f)

    # convert the image to values between 0 and 255
    # this is required by go.Image
    img = (img * 255).astype(int)

    # Create the image
    fig = go.Figure(data=[
            go.Image(
                # Note that you can move the image around the screen
                # by setting appropriate values to x0, y0, dx, dy
                x0=x.min(),
                y0=y.min(),
                dx=(x.max() - x.min()) / n2,
                dy=(y.max() - y.min()) / n1,
                z=f
            )
        ],
        layout={
            # set equal aspect ratio and axis labels
            "yaxis": {"scaleanchor": "x", "title": "y"},
            "xaxis": {"title": "x"}
        }
    )
    return fig


def number_of_beams_calculation(file):
    return(int(file.FractionGroupSequence[0].NumberOfBeams))

def select_mlc(control_point):
    for i in range(len(control_point)):
        if control_point[i].RTBeamLimitingDeviceType == 'MLCX':
            mlc_id = i
    return mlc_id

def calculafe_filed_size_varian(file, beam_number, control_point, number_of_leafs,  **kwargs):
    
    leafs = kwargs.pop("leafs", None)
    size = pd.DataFrame()
    
    if leafs == None:
        leafs = range(number_of_leafs)
        print('The machine has '+ str(number_of_leafs)+ ' leaf pairs')
    #field size
    for c in (range(control_point)):
        
        
        if hasattr(file.BeamSequence[beam_number].ControlPointSequence[c], 'BeamLimitingDevicePositionSequence'):
            mlc_id = select_mlc(file.BeamSequence[beam_number].ControlPointSequence[c].BeamLimitingDevicePositionSequence)
            size.loc[c, 'field_size'] = 0
            for p in leafs:
                leaf_upper_position = file.BeamSequence[beam_number].BeamLimitingDeviceSequence[2].LeafPositionBoundaries[p]
                leaf_lower_position = file.BeamSequence[beam_number].BeamLimitingDeviceSequence[2].LeafPositionBoundaries[p+1]
                
                leaf_left_position = file.BeamSequence[beam_number].ControlPointSequence[c].BeamLimitingDevicePositionSequence[mlc_id].LeafJawPositions[p]
                leaf_right_position = file.BeamSequence[beam_number].ControlPointSequence[c].BeamLimitingDevicePositionSequence[mlc_id].LeafJawPositions[p+number_of_leafs]
                
                
                size_upper_lower = leaf_lower_position - leaf_upper_position
                size_left_right = leaf_right_position - leaf_left_position
                
                #if size_left_right != 0:
                    #print('size_upper_lower' + str(size_upper_lower))
                    #print('size_left_right' + str(size_left_right))
        
                size.loc[c, 'field_size'] = size.loc[c, 'field_size'] + size_upper_lower*size_left_right
            
        elif len(range(control_point)) == 2:
            size.loc[c, 'field_size'] = 0
            for p in leafs:
                leaf_upper_position = file.BeamSequence[beam_number].BeamLimitingDeviceSequence[2].LeafPositionBoundaries[p]
                leaf_lower_position = file.BeamSequence[beam_number].BeamLimitingDeviceSequence[2].LeafPositionBoundaries[p+1]
                leaf_right_position = file.BeamSequence[beam_number].ControlPointSequence[c-1].BeamLimitingDevicePositionSequence[2].LeafJawPositions[p+number_of_leafs]
                leaf_left_position = file.BeamSequence[beam_number].ControlPointSequence[c-1].BeamLimitingDevicePositionSequence[2].LeafJawPositions[p]
    
                size_upper_lower = leaf_lower_position - leaf_upper_position
                size_left_right = leaf_right_position - leaf_left_position
        
                size.loc[c, 'field_size'] = size.loc[c, 'field_size'] + size_upper_lower*size_left_right
        else:
            size.loc[c, 'field_size'] = 'error'
            print('ERROR')
            
                
    
    #field size
    for c in (range(control_point)):
        if c == 0:
            size.loc[c, 'weigh'] = file.BeamSequence[beam_number].ControlPointSequence[c].CumulativeMetersetWeight
        else:
            size.loc[c, 'weigh'] = file.BeamSequence[beam_number].ControlPointSequence[c].CumulativeMetersetWeight - file.BeamSequence[beam_number].ControlPointSequence[c-1].CumulativeMetersetWeight
    
    #calculate mean filed size between control points
    for s in range(len(size)):
        if s == 0:
            size.loc[s, 'mean_size'] = 0
        else:
            size.loc[s, 'mean_size'] = (size.loc[s, 'field_size'] + size.loc[s - 1, 'field_size'])/2
    
            
    return size
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# creating a web server instance
server = Flask(__name__)

# creating app
app = dash.Dash(__name__,
                server=server,
                external_stylesheets=external_stylesheets,
                title='Small fileds calculator', )

style_passed={
    'margin-bottom': 5,
    'backgroundColor': '#90EE90', #green,
    'padding': 10}

style_failed={
    'margin-bottom': 5,
    'backgroundColor': '#F08080', #red
    'padding': 10}


style_warning={
    'margin-bottom': 5,
    'backgroundColor': '#FFFACD', #yellow
    'padding': 10}


style_text = {
    'width': '200px',
    'height': '30px',
    'font-size': 14,
    'textAlign': 'right',
    'marginTop': 0,
    'display': 'inline-block',}

style = {
    'width': '200px',
    'height': '30px',
    'font-size': 14,
    'marginTop': 5,
    'textAlign': 'right',
    'display': 'inline-block',}

style_drop_down = {
    'width': '200px',
    'height': '30px',
    'font-size': 14,
    'marginTop': 0,
    'vertical-align': 'top',
    'textAlign': 'right',
    'display': 'inline-block',}

style_column = {
    'verticalAlign': 'top',
    'display': 'inline-block',}

test_types = ['PicketFence', 'Star', 'WL', 'CatPhan', 'EffectiveSmallFields']

test_options_dict = dict({
    'EffectiveSmallFields': html.Div([
                                html.Strong('Chamber: ', style=style_text),
                                dcc.Dropdown(
                                    id='chamber-id',
                                    options=['PinPoint 3D 31016'],
                                    value='PinPoint 3D 31016',
                                    style=style_drop_down),
                                html.Strong('First leaf: ', style=style_text, id='show-dots-text'),
                                dcc.Input(
                                    id='first_leaf',
                                    type="number",
                                    value=0,
                                    style=style),
                                html.Strong('Last leaf: ', id='show-cutoff-text', style=style_text),
                                dcc.Input(
                                    id='last_leaf',
                                    type="number",
                                    value=60,
                                    style=style)]),
    'Star': html.Div([
        html.Strong('Radius:', style=style_text),
        dcc.Input(id='star_radius', 
                  type="number",
                  value= 0.9, 
                  min = 0.1,
                  max = 1, 
                  step=0.01,
                  style=style_drop_down),
        html.Strong('Tolerance:', style=style_text),
        dcc.Input(id='star_tolerance', 
                  type="number",
                  value= 0.5, 
                  min = 0.1,
                  max = 1, 
                  step=0.01,
                  style=style_drop_down)
    ],  
        style={'display': 'inline-block', 'marginTop': 0, 'height': '30px',}),
    'PicketFence': html.Div([
        html.Strong('Tollerance:', style=style_text),
        dcc.Input(id='pf_tollerance', 
                  type="number",
                  value= 0.2, 
                  min = 0.1,
                  max = 0.5, 
                  step=0.01,
                  style=style_drop_down)
    ],  
        style={'display': 'inline-block',  'marginTop': 0, 'height': '30px',}),
    'WL': html.Div([html.P('WL test under development')]),
    'CatPhan': html.Div([html.P('CatPhan test under development')])
})

#function to analyse picket fence test
def analyze_picket_fence(files, list_of_names, implementation, tollerance_pf):
    output = []
    clean_temp_folders()
    for l in range(len(files)):
                f = files[l]
                file_id = list_of_names[l]
                file_name = str('./temp/'+str(l)+'_'+str(f.Modality)+'.dcm')
                f.save_as(file_name)
                pf = PicketFence(file_name)
                pf.analyze(tolerance=tollerance_pf)
                if pf.passed:
                    style_f = style_passed
                else:
                    style_f = style_failed
                output.append(html.Div(dbc.Row(dbc.Col(html.P('File {}'.format(file_id))),
                                               style = style_f)))
                image_name = str(str(l) + '_'+ str(implementation)+str(hash(datetime.now())) + '.png')
                pf.save_analyzed_image('./assets/' +image_name)
                
                #creating the report
                pf.publish_pdf('./assets/{}.pdf'.format(image_name))
                #creating a output string
                s = ''
                for picket in pf.pickets: s = s + str(' {:2.2f}'.format(picket.dist2cax))
                
                report_name = '{}.pdf'.format(image_name)
                report_path = './assets/{}.pdf'.format(image_name)
                results_f = html.Div([
                    html.P('{}% Passed'.format(pf.percent_passing)),
                    html.P('Study Date Time: {:%Y-%m-%d %H:%M}'.format(datetime.strptime(f.ContentDate+f.ContentTime, '%Y%m%d%H%M%S.%f'))),
                    html.P('Machine: {}'.format(f.RadiationMachineName)),
                    html.P('Gantry angle: {:3.3f}'.format(f.GantryAngle)),
                    html.P('Collimator angle: {:3.2f}'.format(f.BeamLimitingDeviceAngle)),
                    html.P('Median error: {:10.3f}'.format(pf.abs_median_error)),
                    html.P('Mean picket spacing: {:4.2f}'.format(pf.mean_picket_spacing)),
                    html.P('Picket offsets from CAX (mm): {}'.format(s)),
                    html.P('Max Error: {:2.4}mm on Picket: {}, Leaf: {}'.format(pf.max_error, pf.max_error_picket,pf.max_error_leaf)),
                    html.A('Download Report', download=report_name, href=report_path),
                ], style = {'verticalAlign': 'top', 'margin-top':50})
                
                output.append(html.Div(dbc.Row([
                    dbc.Col(html.Img(src=app.get_asset_url(image_name),)
                                   , style = style_column), 
                    dbc.Col(results_f, style = style_column),
                ], style = {'display': 'inline-block'}), style = {'display': 'inline-block'}))
    return output

#function to analyse star test
def analyze_star(files, list_of_names, implementation, radius, tolerance):
    output = []
    ff_files = []
    
    name_left = str('star'+str(hash(datetime.now()))+'-image_r.png')
    name_right = str('star'+str(hash(datetime.now()))+'-image_l.png')
    
    clean_temp_folders()
    
    for l in range(len(files)):
            f = files[l]
            file_id = list_of_names[l]
            file_name = str('./temp/'+str(l)+'_'+str(f.Modality)+'.dcm')
            f.save_as(file_name)
            ff_files.append(file_name)
    
    mystar = Starshot.from_multiple_images(ff_files)
    try:
        mystar.analyze(radius=radius, tolerance=tolerance)
    except Exception as e:
        return html.P('Star Test analysis failed')
    
    mystar.save_analyzed_subimage('./assets/'+name_right, dpi=200, bbox_inches=None,)
    mystar.save_analyzed_subimage('./assets/'+name_left, subimage = 'out', dpi=200, bbox_inches=None,)

    if mystar.passed:
        style_f = style_passed
    else:
        style_f = style_failed
    output.append(html.Div(dbc.Row(dbc.Col(
        html.P('File analyzed sucsesfully. Passed: {}. Circle diametr {:0.3}'.format(mystar.passed, mystar.results_data().circle_diameter_mm)), 
                                           style = style_f))))
    
    output.append(html.Div(dbc.Row([
        dbc.Col(html.Div(
            html.Img(
                src=app.get_asset_url(name_left), 
                style = {'display': 'inline-block',  'width': '90%', "border":"2px black solid"}
            ), style = {'width':500}), style = {'display': 'inline-block',"border":"2px black solid"}),
        dbc.Col(html.Div(
            html.Img(
                src=app.get_asset_url(name_right), 
                 style = {'display': 'inline-block',  'width': '90%', "border":"2px black solid"}
            ), style = {'width':500}), style = {'display': 'inline-block', "border":"2px black solid"}),
    ], style = {'display': 'inline-block',"border":"2px black solid"}), style = {'display': 'inline-block', "border":"2px black solid"}))
    
    
    return output

#function to analyse effective field size
def effective_field_size_calculation(files, list_of_names, chamber, leaf_min, leaf_max):
    for i in range(len(files)):
        file = files[i]
        filename = list_of_names[i]
        text = []
        leafs = range(leaf_min, leaf_max)
        file_type = file.Modality
        if file_type == 'RTPLAN':
            try:
                plan_name = file.RTPlanLabel
            except Exception:
                plan_name = '-' 
        
            try:
                plan_patient_name = file.PatientName 
            except Exception:
                plan_patient_name ='-'
    
            try:    
                plan_date = file.InstanceCreationDate  
            except Exception:
                plan_date = '-'
        
            try:
                plan_time = file.InstanceCreationTime
            except Exception:
                plan_time = '-'
    
            try:
                plan_vendor = file.Manufacturer
            except Exception:    
                plan_vendor = '-'
    
            if plan_vendor == 'Varian Medical Systems':
                number_of_beams = number_of_beams_calculation(file)
                for n in range(number_of_beams):
                    plan_machime_name = file.BeamSequence[n].TreatmentMachineName  
                    plan_beam_name = file.BeamSequence[n].BeamName
                    plan_beam_type = file.BeamSequence[n].BeamType
                    #print()
            
                    if hasattr(file.FractionGroupSequence[0].ReferencedBeamSequence[n], 'BeamDose'):
                        text.append(str('%s beam name: %s ' %(plan_beam_type, plan_beam_name)))
                        text.append(html.Br())
                
                        plan_beam_dose = file.FractionGroupSequence[0].ReferencedBeamSequence[n].BeamDose
                        plan_beam_meterset = file.FractionGroupSequence[0].ReferencedBeamSequence[n].BeamMeterset
            
                        number_of_leafs = file.BeamSequence[n].BeamLimitingDeviceSequence[2].NumberOfLeafJawPairs

                        number_of_control_points = file.BeamSequence[n].NumberOfControlPoints
                        size_control_point = calculafe_filed_size_varian(file, n, number_of_control_points, number_of_leafs, leafs =leafs)
                
                         
                        #calculate weighted field size
                        size_control_point['weighted_size'] = size_control_point['mean_size']*size_control_point['weigh']
                
                        mean_field_size = size_control_point['weighted_size'].sum()
                        text.append(str('mean field size corrected for the weights is: %s, mm^2' % mean_field_size))
                        text.append(html.Br())
                        text.append(str('effective square filed size is %s, cm' % sqrt(mean_field_size/100)))
                        text.append(html.Br())
                        text.append(str('correction factor for PTW PinPoint 3D 31016 – %s ' %f_ptw_31016(sqrt(mean_field_size/100))))
                        text.append(html.Br())
                
                    else:
                        text.append(str('%s beam name: %s has no dose' %(plan_beam_type, plan_beam_name)))
                        text.append(html.Br())
                    
    export_div = html.Div([
        html.P('The %s was successfully uploaded ' %filename),
        html.P('Modality is %s.' %file.Modality),
        html.P('Patient id is %s.' %file.PatientID),
        html.P('Plan  %s has beams :%s ' %(plan_name, number_of_beams)),
        html.P(text),
        html.Hr(),  # horizontal line
        ])
    
    return export_div


# creating the app layout
app.layout = html.Div([

    dcc.Store(id='memory-output'),
    html.Div([
        dbc.Row([
            dbc.Col('COL1'
            ),
        ], style={'margin-bottom': 5, 'margin-left': 100}),
        html.Br(),

        html.H1('MedPhys portal',
                style={
                    'display': 'inline-block',
                    'padding': 10,
                    'textAlign': 'right',
                    'margin-top': 20,
                    'margin-left': 100,
                    'color': 'white'},
                id='title'),
    ], style={'backgroundColor': '#061c42', 'padding': 10}),


    html.Br(),
    # creating line with parameters
    dbc.Row([
        html.Div([
            html.Strong('Implementation: ', style=style_text),
            dcc.Dropdown(
                id='implementation',
                options=test_types,
                value = 'PicketFence',
                style=style_drop_down
            ),
            
        ], style={'display': 'inline-block'}),
        html.Div(id='test_options', style={'display': 'inline-block', 'padding': 10})

    ],
        style={'height': '30px', 'marginTop': 0,'padding': 10, 'width': '98%','vertical-align': 'top',}),
    html.Br(style={'height': 2}),
    # uploading field
    
    dcc.Upload(
        id='upload-data',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select Files')
        ]),
        style={
            'width': '97%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        # Allow multiple files to be uploaded
        multiple=True),
    # place for the output
    
    html.Div(id='output-data-upload-pf'),
    html.Div(id='output-data-upload-star'),
    html.Div(id='output-data-upload-eff-fs'),
])

def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    file = dcm.dcmread(io.BytesIO(decoded), force=True)
    return file

#callback to select implementation options:
@app.callback(Output('test_options','children'),
             [Input('implementation', 'value')])
def change_implementation(implementation_slt):
    return test_options_dict[implementation_slt]

#effective field size
@app.callback(Output('output-data-upload-eff-fs', 'children'),
            [Input('implementation','value'),
             Input('first_leaf','value'),
             Input('last_leaf','value'),
             Input('upload-data','contents')],
            State('upload-data','filename'),
            State('upload-data','last_modified'))
def update_output_fs(implementation, leaf_first, leaf_last,list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        output = []
        try:
            files = [
                parse_contents(c, n, d) for c, n, d in
                zip(list_of_contents, list_of_names, list_of_dates)]
        except Exception as e:
            print(e)
            return html.P('Error in uploading file ' + str(e))
        
        amount_of_files = len(files)
        output.append(html.P('%s files was uploaded successfully' %amount_of_files))

        if implementation == 'EffectiveSmallFields':
            chamber = 'PinPoint'

            output.append(effective_field_size_calculation(files, list_of_names, chamber, leaf_first, leaf_last))
        
        return output

#picket fence test
@app.callback(Output('output-data-upload-pf', 'children'),
              [Input('implementation', 'value'),
               Input('pf_tollerance', 'value'),
               Input('upload-data', 'contents'),],
              State('upload-data', 'filename'),
              State('upload-data', 'last_modified'))
def update_output_pf(implementation, tollerance_pf, list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        output = []
        try:
            files = [
                parse_contents(c, n, d) for c, n, d in
                zip(list_of_contents, list_of_names, list_of_dates)]
        except Exception as e:
            print(e)
            return html.P('Error in uploading file ' + str(e))
        amount_of_files = len(files)
        output.append(html.P('%s files was uploaded successfully' %amount_of_files))
        
        # Picket fence test analyzation
        if implementation == 'PicketFence':
            output.append(analyze_picket_fence(files, list_of_names, implementation, tollerance_pf))
        return output

    
#star test
@app.callback(Output('output-data-upload-star', 'children'),
              [Input('implementation', 'value'),
               Input('star_radius', 'value'),
               Input('star_tolerance', 'value'),
               Input('upload-data', 'contents'),],
              State('upload-data', 'filename'),
              State('upload-data', 'last_modified'))
def update_output_star(implementation, radius, tolerance, list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        output = []
        try:
            files = [
                parse_contents(c, n, d) for c, n, d in
                zip(list_of_contents, list_of_names, list_of_dates)]
        except Exception as e:
            print(e)
            return html.P('Error in uploading file ' + str(e))
        amount_of_files = len(files)
        output.append(html.P('%s files was uploaded successfully' %amount_of_files))
        
        # Picket fence test analyzation
        if implementation == 'Star':
            output = analyze_star(files, list_of_names, implementation, radius, tolerance)
        return output

if __name__ == '__main__':
    app.run_server(debug=False, host='0.0.0.0', port=8055)
