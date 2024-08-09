#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
POSSUM Status Monitor

This tool collects various functions for updating and checking the status
of POSSUM data, including observations, tiles, and pipeline products.

At present, all of this information is being kept in a central Google Sheet.
An API key is needed to access Google sheets through the Python module. This
needs to be provided to the script through a file. Users can try to set up
their own OAuth token if they just need read access; for write access you will
either need to be given edit permission, or have a bot API token; ask Cameron
about this to get access.

A command line interface will be included, in order to allow automatic updates.

Sheet URL: 
https://docs.google.com/spreadsheets/d/1sWCtxSSzTwjYjhxr1_KVLWG2AnrHwSJf_RWQow7wbH0

@author: cvaneck
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.markers import MarkerStyle
import astropy.table as at
import astropy.coordinates as ac
import astropy.units as u
import healpy as hp
import gspread
import pyvo as vo
import datetime
import cartopy.crs as ccrs
from time import sleep
import psycopg2
import json


def open_sheet(token_file):
    """Opens the POSSUM status Google Sheet as a Python variable.
    Requires a token file to authorize the Google Sheet API access.
    """
    
    gc = gspread.service_account(filename=token_file)
    ps= gc.open_by_url('https://docs.google.com/spreadsheets/d/1sWCtxSSzTwjYjhxr1_KVLWG2AnrHwSJf_RWQow7wbH0')
    
    return ps

def _get_sheet(ps,survey,data):
    """Open a particular sheet within the survey master Google sheet file, 
    based on which survey, band, and data type (observations or tiles)
    are specified.
    Inputs:
        ps: master sheet variable
        survey (str): survey to process: '1'/'2' for full survey band 1/2, 
                        'p1'/'p2' for pilot survey.
        data: 'Tiles' or 'Observations'
    """
    if survey[0].lower() == 'p':
        worksheet_name = f'Pilot {data} - Band {survey[1]}'
    else:
        worksheet_name = f'Survey {data} - Band {survey[0]}'
    return ps.worksheet(worksheet_name)

    
def _access_database(auth_file):
    """
    Access the AusSRC database. Takes in a JSON filename (containing the login
    information) and returns a connection object.
    """
    with open(auth_file,'r') as f:
        db_details = json.load(f)

    connection = psycopg2.connect(host=db_details['host'],
                                  user=db_details['user'],
                                  password=db_details['password'],
                                  database='possum')
    return connection

def verify_sheet(ps):
    """Verifies that the Google Sheet file matches script expectations (in 
    terms of what sheets are present and the column header). If the sheet
    changes (e.g., adding new columns), update all the code that interacts
    with that sheet then update this function to match the new state."""
    
    #Check names of sheets:
    sheets=[x.title for x in ps.worksheets()]
    assert 'Pilot Observations - Band 1' in sheets
    assert 'Pilot Observations - Band 2' in sheets
    assert 'Pilot Tiles - Band 1' in sheets
    assert 'Pilot Tiles - Band 2' in sheets
    assert 'Survey Observations - Band 1' in sheets
    assert 'Survey Tiles - Band 1' in sheets
    

    #Check headers of individual sheets. If someone adds new columns, the code
    #needs to be updated.
    #Check all observation sheets:
    for obs in [x for x in sheets if 'Observations' in x]:
        sheet=ps.worksheet(obs)
        header = sheet.row_values(1)
        try:
            assert header == ['name', 'ra', 'dec', 'ra_deg', 'dec_deg', 'gl', 'gb',
                          'rotation', 'duration', 'centrefreq', 'bandwidth',
                          'footprint', 'associated_tiles', 'obs_start', 'obs_end',
                          'sbid', 'processed', 'validated', 'aus_src','single_SB_1D_pipeline']
        except:
            raise Exception(f'{obs} header had changed.')

    for tiles in [x for x in sheets if 'Tiles' in x]:
        sheet = ps.worksheet(tiles)
        header = sheet.row_values(1)
        try:
            assert header == ['tile_id', 'ra', 'dec', 'ra_deg', 'dec_deg', 'gl', 'gb',
                          'aus_src','1d_pipeline_main','1d_pipeline_borders','3d_pipeline',
                          'field1', 'field2', 'field3', 'field4', 'field5',
                          'field6', 'field7', 'field8', 'field9', 'neighbour_SW',
                          'neighbour_W', 'neighbour_NW', 'neighbour_N',
                          'neighbour_NE', 'neighbour_E', 'neighbour_SE',
                          'neighbour_S']
        except:
            raise Exception(f'{tiles} header had changed.')







def tile_observation_maps(map_file):
    """Create dictionaries that contain the mappings connecting observations to
    tiles, and vice versa.
    Require a input file that contains the mapping (such as generated by 
    Lerato's scripts), with the tile numbers in the first column (titled 'PIXEL'')
    and the field/observation names in the following columns (titled 'SB#', with
    # from 1 to 9). 
    Make sure this file also contains all the tiles that have only a single
    contributing observation!
    
    Returns tile_dict (tile number as keys, list of contributing fields as values)
    and field_dict (field names as keys, list of associated tiles as values)."""
    mapping = at.Table.read(map_file)

    tile_dict={}
    for row in mapping:
        tile_dict[row['PIXEL']] = [ x for x in list(row.as_void())[1:] if x != '0' ] 
    tile_dict=dict(sorted(tile_dict.items())) #Sort, for convenience

    field_dict={}
    for tile in tile_dict:
        field_list=tile_dict[tile]
        for field in field_list:
            if field not in field_dict:
                field_dict[field] = [tile]
            else:
                field_dict[field].append(tile)

    number_of_observations=[len(tile_dict[x]) for x in tile_dict]
    if min(number_of_observations) < 1:
        raise Exception("All tiles have 2+ observations. Are the single-observations tiles missing?")
    
    return tile_dict, field_dict


def query_CASDA_TAP(sb):
    """Gets information from CASDA for a supplied SB number.
    Returns 4 element tuple: observation start and end dates (as strings), the deposited date, and validated date.
    The validated date may be None if the SB has not been validated yet."""
    session=vo.tap.TAPService('https://casda.csiro.au/casda_vo_tools/tap')
    result=session.search(f"SELECT TOP 100 sbid,obs_start,obs_end,event_date,event_type,obs_program,project_code FROM casda.observation_event WHERE sbid = '{sb}' AND (project_code = 'AS103' OR project_code = 'AS203' )").to_table()
    if ((result['event_type'] == 'DEPOSITED') | (result['event_type'] == 'UPDATED')).sum() >= 1:
        start=result['obs_start'][0]
        end=result['obs_end'][0]
        deposit_date=result[(result['event_type'] == 'DEPOSITED') | (result['event_type'] == 'UPDATED')]['event_date'][0]
    else:
        start=None
        end=None
        deposit_date=None
    if (result['event_type'] == 'RELEASED').sum() > 1:
        raise Exception(f"Wrong number of RELEASED events on SB{sb}? What the heck?")
#Search for rejections first: CASDA records both a VALIDATED and REJECTED event for rejected data.
    elif (result['event_type'] == 'REJECTED').sum() == 1: 
        valid_date='REJECTED - '+result[result['event_type'] == 'REJECTED']['event_date'][0]
    elif (result['event_type'] == 'RELEASED').sum() == 1:
        valid_date=result[result['event_type'] == 'VALIDATED']['event_date'][0]
    else:
        valid_date=None

#    print('Successful query.')
    return start,end,deposit_date,valid_date



def update_observed_fields(ps,db_auth_file):
    """Update the information on which fields have been observed, using
    Vanessa Moss' EMU and WALLABY spreadsheets plus CASDA as the information 
    sources. Takes the POSSUM master sheet variable and the AusSRC database
    authenatication file as inputs.
    """

    #Get as-is POSSUM status
    sheet = ps.worksheet('Survey Observations - Band 1')
    current_data=sheet.get_values()
    sleep(1)
    current_data=at.Table(np.array(current_data)[1:],names=current_data[0])

    #Get current EMU status.
    obs_sh = ps.client.open_by_url("https://docs.google.com/spreadsheets/d/1HjSWDvwknndCi5PxTP4aDHYAdQiSliRYFVwZGf8YZJw")
    emu_obs=obs_sh.sheet1.get_values()
    sleep(1)
    emu_obs=at.Table(np.array(emu_obs)[1:],names=emu_obs[0])

    #Update any new observations:


    #Catch observations that have had SB# added since last update.
    new_obs = np.where((current_data['sbid'] == '')  & (emu_obs['sbid'] != '0'))[0]
    #Update rows for all new observations (new SBs)
    for i in new_obs:
        sb=emu_obs[i]['sbid']
        sheet.update(f'P{i+2}',sb)
        sleep(1)


    #Catch observations that have had SB# changed (reobservations) since last update.
    reobs = np.where((current_data['sbid'] != emu_obs['sbid']) & 
                       (emu_obs['sbid'] != '0') & (current_data['sbid'] != ''))[0]
    # #Update rows for all new observations -- removing old obs dates and validation
    for i in reobs:
        sb=emu_obs[i]['sbid']
        sheet.update(f'N{i+2}',[['','',sb,'','']])
        sleep(1)

    cancelled = np.where((current_data['sbid'] != emu_obs['sbid']) & 
                       ((emu_obs['sbid'] == '0') & (current_data['sbid'] != '') & (current_data['processed'] == '')))[0]
    for i in cancelled:
        sb=emu_obs[i]['sbid']
        sheet.update(f'N{i+2}',[['','','','','']])
        sleep(1)


    #Check recent observations for processing
    awaiting_processing = np.where((current_data['sbid'] != '') & (current_data['obs_start'] == ''))[0]
    for i in awaiting_processing:
        sb=emu_obs[i]['sbid']
        obs_start,obs_end,deposit_date=query_CASDA_TAP(sb)[0:3]
        if obs_start is not None:
            sheet.update(f'N{i+2}',[[obs_start,obs_end,sb,deposit_date[0:10]]])
            sleep(1)

        
    #Check/update validation status for previously unvalidated observations
    awaiting_validation=np.where((current_data['processed'] != '') & (current_data['validated'] == ''))[0]
    reject_sheet=ps.worksheet('Rejected Observations')


    for i in awaiting_validation:
        sb=current_data[i]['sbid']
        valid_date=query_CASDA_TAP(sb)[3]
        if valid_date is not None:
            if 'REJECTED' in valid_date:
                sheet.update(f'R{i+2}',valid_date)
                sheet.format(f'R{i+2}',{'backgroundColorStyle':{"rgbColor": {"red":1,"green": 0.,"blue": 0.}}})
                sleep(2)
                
                row=sheet.get_values(f'A{i+2}:Z{i+2}')
                testcol=reject_sheet.get_values("A:A")
                sleep(2)
                newrow_num=len(testcol)+1
                reject_sheet.update(f'A{newrow_num}:Z{newrow_num}',row)
                sleep(1)


            else:
                sheet.update(f'R{i+2}',valid_date[0:10])
                sheet.format(f'R{i+2}',{'backgroundColorStyle':{"rgbColor": {"red":0.,"green": 0.,"blue": 1.}}})
                sleep(2)


    #Check for AusSRC processing status
    connection = _access_database(db_auth_file)
    cursor = connection.cursor()
    try:
        cursor.execute("""SELECT name,sbid,cube_update FROM possum.observation 
        WHERE cube_state = 'COMPLETED' AND band = '1'  """)
        db_data=cursor.fetchall()
    finally:
        connection.rollback()
    db_data=at.Table(np.array(db_data),names=['name','sbid','cube_update'])
    db_data['sbid'] = [ x.split('-')[1] for x in db_data['sbid']] #get rid of first part of  'ASKAP-<sbid>'

    for i,sb in enumerate(db_data['sbid']): #Iterate through all processed SBs to see if status is new.
        sb_row = np.where(current_data['sbid'] == sb)[0]
        if sb_row.size != 1:
            print(f"For some reason there's either zero or more than 1 sb={sb} rows in the status monitor?")
            continue
        if (current_data[sb_row[0]]['validated'] != '') and (current_data[sb_row[0]]['validated'].startswith('20')):
            if (current_data[sb_row[0]]['aus_src'] == '') or (current_data[sb_row[0]]['aus_src'] != db_data[i]['cube_update'].strftime("%Y-%m-%d")):
                #Update cell if blank or has different timestamp
                #print(f'S{sb_row[0]+2}',db_data[i]['cube_update'].strftime("%Y-%m-%d"))
                sheet.update(f'S{sb_row[0]+2}',db_data[i]['cube_update'].strftime("%Y-%m-%d"))
                sheet.format(f'S{sb_row[0]+2}',{'backgroundColorStyle':{"rgbColor": {"red":0,"green": 1.,"blue": 0.}}})
                sleep(2)



    ##WALLABY
    sheet = ps.worksheet('Survey Observations - Band 2')
    current_data=sheet.get_values()
    sleep(1)
    current_data=at.Table(np.array(current_data)[1:],names=current_data[0])

    
    obs_sh = ps.client.open_by_url("https://docs.google.com/spreadsheets/d/1w9hn7b1q8QgeMU7D-4ivkMeBzf-mzvVSfB9ABANd78U/edit#gid=0")
    emu_obs=obs_sh.sheet1.get_values()
    emu_obs=at.Table(np.array(emu_obs)[1:],names=emu_obs[0])


    #Catch observations that have had SB# added since last update.
    new_obs = np.where((current_data['sbid'] == '') & (emu_obs['sbid'] != '0'))[0]
    #Update rows for all new observations (new SBs)
    for i in new_obs:
        sb=emu_obs[i]['sbid']
        sheet.update(f'P{i+2}',sb)
        sleep(1)


    #Catch observations that have had SB# changed (reobservations) since last update.
    reobs = np.where((current_data['sbid'] != emu_obs['sbid']) & 
                       (emu_obs['sbid'] != '0') & (current_data['sbid'] != ''))[0]
    # #Update rows for all new observations -- removing old obs dates and validation
    for i in reobs:
        sb=emu_obs[i]['sbid']
        sheet.update(f'N{i+2}',[['','',sb,'','']])
        sleep(1)


    #Check recent observations for processing
    awaiting_processing = np.where((current_data['sbid'] != '') & (current_data['obs_start'] == ''))[0]
    for i in awaiting_processing:
        sb=emu_obs[i]['sbid']
        obs_start,obs_end,deposit_date=query_CASDA_TAP(sb)[0:3]
        if obs_start is not None:
            sheet.update(f'N{i+2}',[[obs_start,obs_end,sb,deposit_date[0:10]]])
            sleep(1)

        
    #Check/update validation status for previously unvalidated observations
    awaiting_validation=np.where((current_data['processed'] != '') & (current_data['validated'] == ''))[0]
    for i in awaiting_validation:
        sb=current_data[i]['sbid']
        valid_date=query_CASDA_TAP(sb)[3]
        if valid_date is not None:
            if 'REJECTED' in valid_date:
                sheet.update(f'R{i+2}',valid_date)
                sheet.format(f'R{i+2}',{'backgroundColorStyle':{"rgbColor": {"red":1,"green": 0.,"blue": 0.}}})
                sleep(2)
                
                row=sheet.get_values(f'A{i+2}:Z{i+2}')
                testcol=reject_sheet.get_values("A:A")
                sleep(2)
                newrow_num=len(testcol)+1
                reject_sheet.update(f'A{newrow_num}:Z{newrow_num}',row)
                sleep(1)

            else:
                sheet.update(f'R{i+2}',valid_date[0:10])
                sheet.format(f'R{i+2}',{'backgroundColorStyle':{"rgbColor": {"red":0.,"green": 0.,"blue": 1.}}})
                sleep(2)

    try:
        cursor.execute("""SELECT name,sbid,cube_update FROM possum.observation 
        WHERE cube_state = 'COMPLETED' AND band = '2'  """)
        db_data=cursor.fetchall()
    finally:
        connection.rollback()
    connection.close()
    db_data=at.Table(np.array(db_data),names=['name','sbid','cube_update'])
    db_data['sbid'] = [ x.split('-')[1] for x in db_data['sbid']] #get rid of first part of  'ASKAP-<sbid>'

    for i,sb in enumerate(db_data['sbid']): #Iterate through all processed SBs to see if status is new.
        sb_row = np.where(current_data['sbid'] == sb)[0]
        if sb_row.size != 1:
            print(f"For some reason there's either zero or more than 1 sb={sb} rows in the status monitor?")
            continue
        if (current_data[sb_row[0]]['validated'] != '') and (current_data[sb_row[0]]['validated'].startswith('20')):
            if (current_data[sb_row[0]]['aus_src'] == '') or (current_data[sb_row[0]]['aus_src'] != db_data[i]['cube_update'].strftime("%Y-%m-%d")):
                #Update cell if blank or has different timestamp
                #print(f'S{sb_row[0]+2}',db_data[i]['cube_update'].strftime("%Y-%m-%d"))
                sheet.update(f'S{sb_row[0]+2}',db_data[i]['cube_update'].strftime("%Y-%m-%d"))
                sheet.format(f'S{sb_row[0]+2}',{'backgroundColorStyle':{"rgbColor": {"red":0,"green": 1.,"blue": 0.}}})
                sleep(2)






def update_aussrc_field_processed(ps,survey,sb):
    """Update an observation to indicate it has been tiled by the AusSRC pipeline.
    Notes the date in the observation sheet, and color-codes the associated cells
    in the tile sheet.
    Inputs:
        ps: master sheet variable
        survey (str): survey to process: '1'/'2' for full survey band 1/2, 
                        'p1'/'p2' for pilot survey.
        sb (int or str): SBID number for observation"""
    
    obs_sheet = _get_sheet(ps,survey,'Observations')
    current_data=obs_sheet.get_values()
    sleep(1)
    current_data=at.Table(np.array(current_data)[1:],names=current_data[0])
    
    #Find row to update; assign today's date to aussrc column.
    w=np.where(current_data['sbid'] == str(sb))[0]
    if w.size != 1:
        raise Exception('Could not uniquely find SB row. Either missing or duplicated?')
    today=datetime.date.today().isoformat()
    obs_sheet.update(f'S{int(w)+2}',today) 
    obs_sheet.format(f'S{int(w)+2}',{'backgroundColorStyle':{"rgbColor": {"red":0.,"green": 1.,"blue": 0.}}})
    sleep(2)
    
    #Get tiles that should have been produced:
    field_name=current_data[w]['name'][0]
    affected_tiles = [ int(x) for x in current_data[w]['associated_tiles'][0].split(',')]

    tile_sheet = ps.worksheet('Survey Tiles - Band 1')
    tile_data=tile_sheet.get_values()
    sleep(1)
    tile_data=at.Table(np.array(tile_data)[1:],names=tile_data[0])
    
    for tile in affected_tiles:
        w = np.where(tile_data['tile_id'] == str(tile))[0]
        if w.size != 1:
            raise Exception('Could not uniquely find tile row. Either missing or duplicated?')
        for i in range(1,9):
            if tile_data[w][f'field{i}'] == field_name[4:]:
                break
        tile_sheet.format(chr(ord('K')+i)+str(w[0]+2),{'backgroundColorStyle':{"rgbColor": {"red":0,"green": 1,"blue": 50/255.}}})
        sleep(1)


def update_aussrc_tile_processed(ps,survey,tile):
    """Update a tile to indicate it has been mosaicked and uploaded by the AusSRC pipeline.
    Notes the date in the tile sheet, and color-codes the appropriate cells for the neighbouring tiles.
    Inputs:
        ps: master sheet variable
        survey (str): survey to process: '1'/'2' for full survey band 1/2, 
                    'p1'/'p2' for pilot survey.
        tile: tile id number"""

    tile_sheet = _get_sheet(ps,survey,'Tiles')
    tile_data=tile_sheet.get_values()
    sleep(1)
    tile_data=at.Table(np.array(tile_data)[1:],names=tile_data[0])

    #Update status cell with date:
    w = np.where(tile_data['tile_id'] == str(tile))[0]
    if w.size != 1:
        raise Exception('Could not uniquely find tile row. Either missing or duplicated?')
    today=datetime.date.today().isoformat()
    tile_sheet.update(f'H{int(w)+2}',today) 
    sleep(1)

    #Small function to get the right column name, accounting for the wrapping from Z to AA.
    def get_column(start,i):
        x=ord(start)+(i+4)%8
        if x < 91:
            return chr(x)
        else:
            return 'A'+chr(x-26)
    
    #Color this tile in rows of neighbours.
    neighbours = [ int(x) if x != '' else 0 for x in tile_data[w[0]][tile_data.colnames[20:]].values() ]
    for i in range(len(neighbours)):
        if neighbours[i] != 0:
            w = np.where(tile_data['tile_id'] == str(neighbours[i]))[0]
            if w.size != 1:
                raise Exception('Could not uniquely find tile row. Either missing or duplicated?')
            tile_sheet.format(get_column('U',i)+str(w[0]+2),{'backgroundColorStyle':{"rgbColor": {"red":0,"green": 1,"blue": 50/255.}}})
            sleep(1)


def create_plots(ps,survey,basename):
    """Creates plots showing the survey status, and saves them to disk.
    Makes the plots for both the tiles and observations. Currently band 1 only.
    Inputs:
        ps: the master Google Sheet variable
        basename (str): path+base filename for plots. Each plot will append
        its specific name + .pdf.
    """
    
    #First the tile plots. Get the data from the relevant sheet.
    tile_sheet=_get_sheet(ps,survey,'Tiles')
    tile_info=tile_sheet.get_values()
    sleep(1)
    tile_info=at.Table(np.array(tile_info)[1:],names=tile_info[0])

    #Since the tiles are HEALPix pixels, using the HealPy plotting functionality.
    #This requires making a HEALPix array (a list of pixel values).
    tile_map=np.zeros(hp.nside2npix(32),dtype='int')
    tile_map[tile_info['tile_id'].astype('int')] = 1 #Set 1 = inside survey area.
    tile_map[tile_info['tile_id'].astype('int')[tile_info['aus_src'] != '']] = 2
    
    pipeline=np.vstack(([x.startswith('20') for x in tile_info['1d_pipeline_main']],
                        [x.startswith('20') for x in tile_info['1d_pipeline_borders']],
                        [x.startswith('20') for x in tile_info['3d_pipeline']]))
    
    tile_map[tile_info['tile_id'].astype('int')[[x.startswith('20') for x in tile_info['1d_pipeline_main'] ] ]] = 3
    tile_map[tile_info['tile_id'].astype('int')[[x.startswith('20') for x in tile_info['3d_pipeline']]]] = 4
    tile_map[tile_info['tile_id'].astype('int')[ np.array([x.startswith('20') for x in tile_info['1d_pipeline_main']]) & np.array([x.startswith('20') for x in tile_info['3d_pipeline']])]] = 5
    tile_map[tile_info['tile_id'].astype('int')[np.all(pipeline,axis=0)]] = 6
    
    
    col_dict={0:'black',
              1:"grey",
              2:"orange",
              3:"tab:blue",
              4:"tab:red",
              5:"darkviolet",
              6:"green"}
    labels = np.array(['Outside Survey',
                       f'Incomplete ({(tile_map==1).sum()})',
                       f'Mosaicked ({(tile_map==2).sum()})',
                       f'1D processed (core) ({(tile_map==3).sum()})',
                       f'3D processed ({(tile_map==4).sum()})',
                       f'1D(core)+\n3D processed ({(tile_map==5).sum()})',
                       f'Fully Processed ({(tile_map==6).sum()})',
                      ])
    cm = mpl.colors.ListedColormap([col_dict[x] for x in col_dict.keys()])
    norm_bins = np.sort([*col_dict.keys()]) + 0.5
    norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)
    norm = mpl.colors.BoundaryNorm(norm_bins, len(labels), clip=True)
    fmt = mpl.ticker.FuncFormatter(lambda x, pos: labels[norm(x)])
    diff = norm_bins[1:] - norm_bins[:-1]
    tickz = norm_bins[:-1] + diff / 2


    data=hp.mollview(tile_map,title='Survey tile status',cmap=cm,min=0,max=4,return_projected_map=True,xsize=10000)
    plt.figure(figsize=(12,6))
    ax=plt.subplot(projection=ccrs.Mollweide())
    im=ax.imshow(data,transform=ccrs.Mollweide(),cmap=cm,norm=norm,extent=(-18040095.696147293, 18040095.696147293,-9020047.848073646, 9020047.848073646),origin='lower')
    ax.gridlines(draw_labels=False,xlocs=[0,45,90,135,180,-45,-90,-135,-179.99],ylocs=[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])
    plt.colorbar(im,format=fmt, ticks=tickz,shrink=0.7)
    for ra_grid in [0,3,6,9,12]:
        ax.text(-15*ra_grid,15,str(int(ra_grid))+'h',{'ha':'right','va':'top','color':'w'},transform=ccrs.PlateCarree())
    for ra_grid in [12.01,15,18,21]:
        ax.text(-15*ra_grid,15,str(int(ra_grid))+'h',{'ha':'left','va':'top','color':'w'},transform=ccrs.PlateCarree())
    for dec_grid in [75,60,45,30,15,0,-15,-30,-45,-60,-75]:
        ax.text(0,dec_grid,str(dec_grid)+'°',{'ha':'left','va':'bottom','color':'w'},transform=ccrs.PlateCarree())
    xmax=plt.xlim()[1]
    ymax=plt.ylim()[1]
    plt.text(xmax*-0.99,ymax*-0.95,'Equatorial\ncoordinates',{'size':14,'weight':'bold','ha':'left'})
    plt.text(xmax*-0.99,ymax*0.95,'Tiles',{'size':14,'weight':'bold','ha':'left'})

    #Draw lines for HPX projection boundaries
    ax.plot(-np.repeat(0,10),np.linspace(41.8103,-41.8103,10),'r:',transform=ccrs.PlateCarree(),label='HPX projection\nboundaries')
    for ra in [0,90,179.99,180.01,270]:
        ax.plot(-np.repeat(ra,10),np.linspace(-41.8103,-90,10),'r:',transform=ccrs.PlateCarree())
        ax.plot(-np.repeat(ra,10),np.linspace(41.8103,90,10),'r:',transform=ccrs.PlateCarree())
    plt.legend(loc='upper right')


    plt.savefig(basename+'tiles_equatorial.png',bbox_inches='tight',dpi=300)
    plt.close()

    data=hp.mollview(tile_map,title='Survey tile status',coord='CG',cbar=False,return_projected_map=True,xsize=10000)
    #hp.graticule(dmer=45,dpar=15)
    plt.figure(figsize=(12,6))
    ax=plt.subplot(projection=ccrs.Mollweide())
    im=ax.imshow(data,transform=ccrs.Mollweide(),cmap=cm,norm=norm,extent=(-18040095.696147293, 18040095.696147293,-9020047.848073646, 9020047.848073646),origin='lower')
    ax.gridlines(draw_labels=False,xlocs=[0,45,90,135,180,-45,-90,-135,-179.99],ylocs=[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])
    plt.colorbar(im,format=fmt, ticks=tickz,shrink=0.7)
    for ra_grid in [0,45,90,135,180]:
        ax.text(-1*ra_grid,0,str(int(ra_grid))+'°',{'ha':'right','va':'top','color':'w'},transform=ccrs.PlateCarree())
    for ra_grid in [180.01,225,270,315]:
        ax.text(-1*ra_grid,0,str(int(ra_grid))+'°',{'ha':'left','va':'top','color':'w'},transform=ccrs.PlateCarree())
    for dec_grid in [75,60,45,30,15,0,-15,-30,-45,-60,-75]:
        ax.text(0,dec_grid,str(dec_grid)+'°',{'ha':'left','va':'bottom','color':'w'},transform=ccrs.PlateCarree())
    xmax=plt.xlim()[1]
    ymax=plt.ylim()[1]
    plt.text(xmax*-0.99,ymax*-0.95,'Galactic\ncoordinates',{'size':14,'weight':'bold','ha':'left'})
    plt.text(xmax*-0.99,ymax*0.95,'Tiles',{'size':14,'weight':'bold','ha':'left'})


    #Draw lines for HPX projection boundaries
    coords=ac.SkyCoord(ra=np.repeat(0,20),dec=np.linspace(90,-90,20),frame='icrs',unit='deg')
    ax.plot([-x if x < 180 else 360-x for x in coords.galactic.l.deg],coords.galactic.b.deg,'r:',transform=ccrs.PlateCarree(),label='HPX projection\nboundaries')
    for ra in [90,180,270]:
        coords=ac.SkyCoord(ra=np.repeat(ra,10),dec=np.linspace(-41.8103,-90,10),frame='icrs',unit='deg')
        ax.plot(-coords.galactic.l.deg,coords.galactic.b.deg,'r:',transform=ccrs.PlateCarree())
        coords=ac.SkyCoord(ra=np.repeat(ra,10),dec=np.linspace(41.8103,90,10),frame='icrs',unit='deg')
        ax.plot(-coords.galactic.l.deg,coords.galactic.b.deg,'r:',transform=ccrs.PlateCarree())
    plt.legend(loc='lower right')


    plt.savefig(basename+'tiles_galactic.png',bbox_inches='tight',dpi=300)
    plt.close()



    #Now the observations.
    obs_sheet=_get_sheet(ps,survey,'Observations')
    obs_info=obs_sheet.get_values()
    sleep(1)
    obs_info=at.Table(np.array(obs_info)[1:],names=obs_info[0])

    #Get the coordinates of each observation.
    #Since the high-dec fields are doubled up, add a Dec offset so both are visible.
    ra=obs_info['ra_deg'].data.astype('float')
    dec=obs_info['dec_deg'].data.astype('float')
    name_mod = [ x[-1] for x in obs_info['name']]

    A_tiles=np.array(name_mod) == 'A'
    B_tiles=np.array(name_mod) == 'B'
    single_tiles=(np.array(name_mod) != 'A') & (np.array(name_mod) != 'B')



    status=np.array([ 1 if x != '' else 0 for x in obs_info['sbid']])
    status[obs_info['processed'] != ''] = 2
    status[obs_info['validated'] != ''] = 4
    status[np.char.startswith(obs_info['validated'],'REJECTED')] = 3
    status[obs_info['aus_src'] != ''] = 5

    #
    col_dict={0:"grey",
              1: "yellow",
              2:"orange",
              3:'tomato',
              4:"blue",
              5:"green"}
    labels = np.array([f'Unobserved ({(status==0).sum()})',
                       f'Observed ({(status==1).sum()})',
                       f'Observatory\nProcessed ({(status==2).sum()})',
                       f'Rejected ({(status==3).sum()})',
                       f'Released ({(status==4).sum()})',
                       f'AusSRC\nProcessed ({(status==5).sum()})',
                      ])
    cm = mpl.colors.ListedColormap([col_dict[x] for x in col_dict.keys()])
    norm_bins = np.sort([*col_dict.keys()]) + 0.5
    norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)
    norm = mpl.colors.BoundaryNorm(norm_bins, len(labels), clip=True)
    fmt = mpl.ticker.FuncFormatter(lambda x, pos: labels[norm(x)])
    diff = norm_bins[1:] - norm_bins[:-1]
    tickz = norm_bins[:-1] + diff / 2


    plt.figure(figsize=(12,6))
    ax=plt.subplot(projection=ccrs.Mollweide())
    ax.gridlines(draw_labels=False,xlocs=[0,45,90,135,180,-45,-90,-135,-179.99],ylocs=[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])
    ax.set_global()

    im=ax.scatter(-1*ra[A_tiles],dec[A_tiles],
                  c=status[A_tiles],transform=ccrs.PlateCarree(),
                  cmap=cm,norm=norm,marker=MarkerStyle("o", fillstyle="top"),edgecolor='k',linewidth=0.2)
    im=ax.scatter(-1*ra[B_tiles],dec[B_tiles],
                  c=status[B_tiles],transform=ccrs.PlateCarree(),
                  cmap=cm,norm=norm,marker=MarkerStyle("o", fillstyle="bottom"),edgecolor='k',linewidth=0.2)
    im=ax.scatter(-1*ra[single_tiles],dec[single_tiles],
                  c=status[single_tiles],transform=ccrs.PlateCarree(),
                  cmap=cm,norm=norm,marker='o',edgecolor='k',linewidth=0.2)


    for ra_grid in [0,3,6,9,12]:
        ax.text(-15*ra_grid,15,str(int(ra_grid))+'h',{'ha':'right','va':'top'},transform=ccrs.PlateCarree())
    for ra_grid in [12.01,15,18,21]:
        ax.text(-15*ra_grid,15,str(int(ra_grid))+'h',{'ha':'left','va':'top'},transform=ccrs.PlateCarree())
    for dec_grid in [75,60,45,30,15,0,-15,-30,-45,-60,-75]:
        ax.text(0,dec_grid,str(dec_grid)+'°',{'ha':'left','weight':'bold','va':'bottom'},transform=ccrs.PlateCarree())
    xmax=plt.xlim()[1]
    ymax=plt.ylim()[1]
    plt.text(xmax*-0.99,ymax*-0.95,'Equatorial\ncoordinates',{'size':14,'weight':'bold','ha':'left'})
    plt.text(xmax*-0.99,ymax*0.95,'Observations',{'size':14,'weight':'bold','ha':'left'})
    plt.colorbar(im,format=fmt, ticks=tickz,shrink=0.7)

    plt.savefig(basename+'observations_equatorial.png',bbox_inches='tight',dpi=300)
    plt.close()

    coords=ac.SkyCoord(ra,dec,unit='deg')
    plt.figure(figsize=(12,6))
    ax=plt.subplot(projection=ccrs.Mollweide())
    ax.gridlines(draw_labels=False,xlocs=[0,45,90,135,180,-45,-90,-135,-179.99],ylocs=[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])
    ax.set_global()

    im=ax.scatter(-1*coords.galactic.l.deg[A_tiles],coords.galactic.b.deg[A_tiles],
                  c=status[A_tiles],transform=ccrs.PlateCarree(),
                  cmap=cm,norm=norm,marker=MarkerStyle("o", fillstyle="top"),edgecolor='k',linewidth=0.2)
    im=ax.scatter(-1*coords.galactic.l.deg[B_tiles],coords.galactic.b.deg[B_tiles],
                  c=status[B_tiles],transform=ccrs.PlateCarree(),
                  cmap=cm,norm=norm,marker=MarkerStyle("o", fillstyle="bottom"),edgecolor='k',linewidth=0.2)
    im=ax.scatter(-1*coords.galactic.l.deg[single_tiles],coords.galactic.b.deg[single_tiles],
                  c=status[single_tiles],transform=ccrs.PlateCarree(),
                  cmap=cm,norm=norm,edgecolor='k',linewidth=0.2)


    for ra_grid in [0,45,90,135,180]:
        ax.text(-1*ra_grid,0,str(int(ra_grid))+'°',{'ha':'right','va':'top'},transform=ccrs.PlateCarree())
    for ra_grid in [180.01,225,270,315]:
        ax.text(-1*ra_grid,0,str(int(ra_grid))+'°',{'ha':'left','va':'top'},transform=ccrs.PlateCarree())
    for dec_grid in [75,60,45,30,15,0,-15,-30,-45,-60,-75]:
        ax.text(0,dec_grid,str(dec_grid)+'°',{'ha':'left','va':'bottom'},transform=ccrs.PlateCarree())
    
    xmax=plt.xlim()[1]
    ymax=plt.ylim()[1]
    plt.text(xmax*-0.99,ymax*-0.95,'Galactic\ncoordinates',{'size':14,'weight':'bold','ha':'left'})
    plt.text(xmax*-0.99,ymax*0.95,'Observations',{'size':14,'weight':'bold','ha':'left'})

    plt.colorbar(im,format=fmt, ticks=tickz,shrink=0.7)

    plt.savefig(basename+'observations_galactic.png',bbox_inches='tight',dpi=300)
    plt.close()



def get_obs_state(obs_info,i):
    """Get the state of an observation, from a previously grabbed observation 
    sheet. Return a short string describing the state.
    i can be a single value (for one field), or a list. If a list, the return
    value will be the most advanced state of all the fields selected.
    """
    if type(i) == int:
        i=[i]
    
    status='planned'
    if ((obs_info[i]['processed'] != '')  & ['REJECTED' not in x for x in obs_info[i]['validated']]).any():
        status='observed'
    if ((obs_info[i]['validated'] != '') & ['REJECTED' not in x for x in obs_info[i]['validated']]).any():
        status='validated'
    if (obs_info[i]['aus_src'] != '').any():
        status='aus_src'
    return status


def aladin_webpage(ps,survey,outfile):
    """Generage an html file that contains an Aladin Lite applet with overlays
    that show the POSSUM status (both observations and tiles).
    """
    #Get field list.
    obs_sheet=_get_sheet(ps,survey,'Observations')
    obs_info=obs_sheet.get_values()
    sleep(1)
    obs_info=at.Table(np.array(obs_info)[1:],names=obs_info[0])
    
    #Constants relating to field size and orientation.
    D = 3.151*u.deg
    PA1 = 51.79*u.deg
    PA2 = 128.21*u.deg
    
    #Compute the corner locations for all the fields.
    #Make list of all tile corners, and list of markers for all field centers.
    field_markers = []
    field_corners = []
    for i in range(len(obs_info)):
        name1 = obs_info[i]['name']
        rot1 = (45+float(obs_info[i]['rotation'] ))*u.deg 
        if name1.startswith('EMU') and name1.endswith('B'): #Skip the 'B's, they have the same sky positions.
            continue
        elif name1.endswith('A'): #Add both A and B SBs to the 'A's
            name1+='/B'
            sbid = obs_info[i]['sbid'] + ' ' + obs_info[i+1]['sbid']
            state=get_obs_state(obs_info, [i,i+1])
        else:
            sbid = obs_info[i]['sbid']
            state=get_obs_state(obs_info, i)
        if sbid == ' ':
            sbid = ''
    
        sc=ac.SkyCoord(obs_info[i]['ra_deg'],obs_info[i]['dec_deg'],unit='deg',equinox='J2000')
        
        field_marker = [name1, sc.ra.deg, sc.dec.deg, rot1.value, sbid]
        field_markers.append( field_marker )
    
        tc1 = sc.directional_offset_by(position_angle = PA1+rot1, separation=D)
        tc2 = sc.directional_offset_by(position_angle = PA2+rot1, separation=D)
        tc3 = sc.directional_offset_by(position_angle = -PA2+rot1, separation=D)
        tc4 = sc.directional_offset_by(position_angle = -PA1+rot1, separation=D)
        field_corners.append( [name1,tc1.ra.deg,tc1.dec.deg,tc2.ra.deg,tc2.dec.deg,tc3.ra.deg,tc3.dec.deg,tc4.ra.deg,tc4.dec.deg,state] )
    
    # Make all tile markers
    field_markers_aladin = []
    for (name,ra,dec,pa,sbid) in field_markers:
        field_markers_aladin.append( f"field_cat.addSources([A.marker({ra},{dec},{{popupTitle: '{name}', popupDesc: '<em>RA:</em> {ra:.3f}<br/><em>Dec:</em> {dec:.3f}<br/><em>PA:</em> {pa:.1f}<br/>(J2000)<br/><em>SBID:</em> {sbid}'}})]);\n" ) 
    
    #Make all field boxes
    field_corners_aladin = []
    field_corner_line =  'field_overlay.addFootprints([A.polygon([{0}])]);\n'
    for field in field_corners:
        corner1 = ', '.join( [ '[{0},{1}]'.format(field[i*2+1],field[i*2+2]) for i in range(4) ] )
        field_corners_aladin.append( field_corner_line.format( corner1 ) ) 
    
    #Observed fields only
    observed_corners_aladin = []
    field_corner_line =  'observed_field_overlay.addFootprints([A.polygon([{0}])]);\n'
    for field in field_corners:
        if field[-1] == 'observed': #Catch only processed fields
            corner1 = ', '.join( [ '[{0},{1}]'.format(field[i*2+1],field[i*2+2]) for i in range(4) ] )
            observed_corners_aladin.append( field_corner_line.format( corner1 ) ) 
    
    #Validated fields only
    validated_corners_aladin = []
    field_corner_line =  'validated_field_overlay.addFootprints([A.polygon([{0}])]);\n'
    for field in field_corners:
        if field[-1] == 'validated': #Catch only processed fields
            corner1 = ', '.join( [ '[{0},{1}]'.format(field[i*2+1],field[i*2+2]) for i in range(4) ] )
            validated_corners_aladin.append( field_corner_line.format( corner1 ) ) 
    
    #Pre-processed fields only
    aussrc_obs_corners_aladin = []
    field_corner_line =  'processed_field_overlay.addFootprints([A.polygon([{0}])]);\n'
    for field in field_corners:
        if field[-1] == 'aus_src': #Catch only processed fields
            corner1 = ', '.join( [ '[{0},{1}]'.format(field[i*2+1],field[i*2+2]) for i in range(4) ] )
            aussrc_obs_corners_aladin.append( field_corner_line.format( corner1 ) ) 
    
    
    
    #Now process the tiles...
    tile_sheet=_get_sheet(ps,survey,'Tiles')
    tile_info=tile_sheet.get_values()
    tile_info=at.Table(np.array(tile_info)[1:],names=tile_info[0])
    
    #Make markers:
    tile_markers_aladin = []
    for row in tile_info:
        if row['1d_pipeline_main'].startswith('20') and '-' in row['1d_pipeline_borders'] and row['3d_pipeline'].startswith('20'):
            status = 'Fully Processed'
        elif row['1d_pipeline_main'].startswith('20') or row['1d_pipeline_borders'].startswith('20') or row['3d_pipeline'].startswith('20'):
            status = 'Partially Processed'
        elif row['aus_src'] != '':
            status = 'Mosiacked'
        else:
            status = 'Incomplete'
        tile_markers_aladin.append( f"tile_cat.addSources([A.marker({row['ra_deg']},{row['dec_deg']},{{popupTitle: 'Tile# {row['tile_id']}', popupDesc: '<em>RA:</em> {row['ra']}<br/><em>Dec:</em> {row['dec']}<br/><em>Status:</em> {status}<br/>'}})]);\n" ) 
    
    
    
    #All tiles boxes:
    all_tile_corners_aladin = []
    for row in tile_info:
        coords=np.array(hp.vec2ang(hp.boundaries(32,int(row['tile_id']),step=1).T,lonlat=True)).T.tolist()
        corner = ','.join([f'[{ra},{dec}]' for ra,dec in iter(coords)])
        all_tile_corners_aladin.append('all_tile_overlay.addFootprints([A.polygon([{0}])]);\n'.format(corner))
        
    #Mosaicked/AusSRC processed tiles:
    mosaicked_tile_corners_aladin = []
    for row in tile_info:
        if row['aus_src'] == '': continue
        coords=np.array(hp.vec2ang(hp.boundaries(32,int(row['tile_id']),step=1).T,lonlat=True)).T.tolist()
        corner = ','.join([f'[{ra},{dec}]' for ra,dec in iter(coords)])
        mosaicked_tile_corners_aladin.append('mosaicked_tile_overlay.addFootprints([A.polygon([{0}])]);\n'.format(corner))
    
    #Partially processed tiles:
    partial_tile_corners_aladin = []
    for row in tile_info:
        if row['1d_pipeline_main'].startswith('20') and '-' in row['1d_pipeline_borders'] and row['3d_pipeline'].startswith('20'): continue
        if not (row['1d_pipeline_main'].startswith('20') or row['1d_pipeline_borders'].startswith('20') or row['3d_pipeline'].startswith('20')): continue
        coords=np.array(hp.vec2ang(hp.boundaries(32,int(row['tile_id']),step=1).T,lonlat=True)).T.tolist()
        corner = ','.join([f'[{ra},{dec}]' for ra,dec in iter(coords)])
        partial_tile_corners_aladin.append('partial_tile_overlay.addFootprints([A.polygon([{0}])]);\n'.format(corner))
        
    #Fully processed tiles:
    completed_tile_corners_aladin = []
    for row in tile_info:
        if not (row['1d_pipeline_main'].startswith('20') and '-' in row['1d_pipeline_borders'] and row['3d_pipeline'].startswith('20')): continue
        coords=np.array(hp.vec2ang(hp.boundaries(32,int(row['tile_id']),step=1).T,lonlat=True)).T.tolist()
        corner = ','.join([f'[{ra},{dec}]' for ra,dec in iter(coords)])
        completed_tile_corners_aladin.append('completed_tile_overlay.addFootprints([A.polygon([{0}])]);\n'.format(corner))


    if survey[0].lower() == 'p':
        page_name = f'Pilot surveys - Band {survey[1]}'
    else:
        page_name = f'Full survey - Band {survey[0]}'


    html_header = '''
    <!DOCTYPE html>
    <html>
    <head>
        <!-- Mandatory when setting up Aladin Lite v3 for a smartphones/tablet usage -->
        <meta name="viewport" content="width=device-width, height=device-height, initial-scale=1.0, user-scalable=no">
    </head>
    <body>
        <h1> POSSUM Status Monitor </h1>
        <h2> ''' + page_name + ''' </h2>
        <br>
        Last updated: ''' + str(datetime.date.today()) + '''
        <br>
        Both observation and tile overlays available; click the 'Stack' button
        to access tile overlays.
    <!-- Aladin Lite has a dependency on the jQuery library -->
    <script src="https://code.jquery.com/jquery-1.10.1.min.js"></script>

    <!-- Aladin Lite container at requested dimensions -->
    <div id="aladin-lite-div" style="width:700px;height:600px;"></div>

    <!-- Aladin Lite JS code -->
    <script type="text/javascript" src="https://aladin.cds.unistra.fr/AladinLite/api/v3/latest/aladin.js" charset="utf-8"></script>

    <!-- Creation of Aladin Lite instance with initial parameters -->

    <script type="text/javascript">
        let aladin;
        A.init.then(() => {
            aladin = A.aladin('#aladin-lite-div', {fov: 240, projection: "SIN", target:'90 -90', cooFrame: 'equatorial', showReticle:false, showCooGridControl: true, showSimbadPointerControl: true, showCooGrid: false});
            let chipass = aladin.createImageSurvey("CHIPASS", "ov-gso/P/CHIPASS", "http://cade.irap.omp.eu/documents/Ancillary/4Aladin/CHIPASS/", "galactic", 3, {imgFormat: 'png'}); // setting a custom HiPS
            let emu = aladin.createImageSurvey('EMU', 'ASKAP-EMU Full Survey', 'http://www.aoc.nrao.edu/~jmarvil/emuHiPS', 'equatorial', 3, {imgFormat: 'png'}); // setting a custom HiPS
            aladin.setImageSurvey(aladin.createImageSurvey("RACS", "RACS-low-I", "https://casda.csiro.au/hips/RACS/low/I", "equatorial", 8, {imgFormat: 'png'})); // setting a custom HiPS
            let racs = aladin.getBaseImageLayer();
            racs.setColormap('grayscale',{stretch: 'linear', reversed: false});
            aladin.setOverlayImageLayer(chipass,'CHIPASS');
            let chimap = aladin.getOverlayImageLayer('CHIPASS');
            chimap.setOpacity(0.5);
            chimap.setColormap('grayscale',{stretch: 'linear', reversed: false});
            chimap.setBlendingConfig(true);
            racs.setBlendingConfig(true);

            var field_overlay = A.graphicOverlay({color: 'grey', lineWidth: 1, name:'Planned Fields'});
            aladin.addOverlay(field_overlay);
            var observed_field_overlay = A.graphicOverlay({color: 'orange', lineWidth: 2, name:'Observed fields'});
            aladin.addOverlay(observed_field_overlay);
            var validated_field_overlay = A.graphicOverlay({color: 'blue', lineWidth: 2, name:'Released fields'});
            aladin.addOverlay(validated_field_overlay);
            var processed_field_overlay = A.graphicOverlay({color: 'green', lineWidth: 2, name:'Pre-processed fields'});
            aladin.addOverlay(processed_field_overlay);
            var field_cat = A.catalog({name: 'Field markers', sourceSize: 10, color: '#eee'});
            aladin.addCatalog(field_cat);
            field_cat.hide();

            var tile_cat = A.catalog({name: 'Tile markers', sourceSize: 10, color: '#eee'});
            aladin.addCatalog(tile_cat);
            tile_cat.hide();

            var all_tile_overlay = A.graphicOverlay({color: 'grey', lineWidth: 1, name:'Planned Tiles'});
            aladin.addOverlay(all_tile_overlay);
            all_tile_overlay.hide();
            var mosaicked_tile_overlay = A.graphicOverlay({color: 'orange', lineWidth: 2, name:'Mosaicked Tiles'});
            aladin.addOverlay(mosaicked_tile_overlay);
            mosaicked_tile_overlay.hide();
            var partial_tile_overlay = A.graphicOverlay({color: 'blue', lineWidth: 2, name:'Partially Processed Tiles'});
            aladin.addOverlay(partial_tile_overlay);
            partial_tile_overlay.hide();
            var completed_tile_overlay = A.graphicOverlay({color: 'green', lineWidth: 2, name:'Fully Processed Tiles'});
            aladin.addOverlay(completed_tile_overlay);
            completed_tile_overlay.hide();


    '''

    footer3 = '''
    });

</script>

        <h3> Built by Cameron Van Eck. Based on the EMU Survey Progress page by Josh Marvil </h3> 


</body>
</html>
    '''


    with open(outfile,'w') as out1:
        out1.writelines(html_header)
        out1.writelines( field_corners_aladin )
        out1.writelines( observed_corners_aladin )
        out1.writelines( validated_corners_aladin )
        out1.writelines( aussrc_obs_corners_aladin )
        out1.writelines( field_markers_aladin )
        out1.writelines(all_tile_corners_aladin)
        out1.writelines(mosaicked_tile_corners_aladin)
        out1.writelines(partial_tile_corners_aladin)
        out1.writelines(completed_tile_corners_aladin)
        out1.writelines(tile_markers_aladin)
        out1.writelines(footer3)




def set_update_date_on_main_page(webpage_file):
    """Writes today's date on the appropriate line of the main status page.
    This is the only information that (currently) gets written to the page html;
    everything else is just updating the plots and Aladin sub-pages.
    Overwrite the file in-place. Takes the page filepath as input.
    """

    with open(webpage_file, 'r') as f:    # pass an appropriate path of the required file
        lines = f.readlines()

        lines[7] = f"Last updated: {str(datetime.date.today())}\n"
        with open(webpage_file, 'w') as f:
            f.writelines(lines)
    


def create_overlap_plots(ps, basename):
    """Creates the plots showing the overlap regions between the two bands.
    Currently, it only tracks 3 data states, on a per-tile basis:
        - planned but not yet observed.
        - Observed (by at least one SB)
        - Aus-SRC processed/tiled (final tile)
    It tracks this for each band, then creates 3 levels for each state:
        - true for one band; not part of 2nd band footprint
        - true for one band, but not yet for 2nd band
        - true for both bands
    
    This implementation is a bit optimistic: it uses the observation-tile mapping,
    so if an observation even slightly touches a tile it gets counted. So it's
    possible for one corner of a tile to be covered in one band and the opposite
    corner in the other band, which would be counted as 'overlap' since both
    bands have that tile.
    """
    #Get all tiles planned for survey:
    tile_sheet=_get_sheet(ps,'1','Tiles')
    band1_tiles=tile_sheet.get_values()
    sleep(1)
    band1_tiles=at.Table(np.array(band1_tiles)[1:],names=band1_tiles[0])
    tile_sheet=_get_sheet(ps,'2','Tiles')
    band2_tiles=tile_sheet.get_values()
    sleep(1)
    band2_tiles=at.Table(np.array(band2_tiles)[1:],names=band2_tiles[0])

    #Get all observations planned for survey:
    obs_sheet=_get_sheet(ps,'1','Observations')
    band1_obs=obs_sheet.get_values()
    sleep(1)
    band1_obs=at.Table(np.array(band1_obs)[1:],names=band1_obs[0])
    obs_sheet=_get_sheet(ps,'2','Observations')
    band2_obs=obs_sheet.get_values()
    sleep(1)
    band2_obs=at.Table(np.array(band2_obs)[1:],names=band2_obs[0])

    #Initialize blank map. Find all tiles in survey area.
    tile_map=np.zeros(hp.nside2npix(32),dtype='int') #Initialize map. 0 = not in survey

    all_tiles=np.unique(np.concatenate((band1_tiles['tile_id'].data,band2_tiles['tile_id'].data))).astype(int)
    tile_map[all_tiles] = 1  # 1 = in survey area in at least one band.

    two_band_tiles=np.array(list(set(band1_tiles['tile_id'].data.astype(int)) & set(band2_tiles['tile_id'].data.astype(int))))
    tile_map[two_band_tiles] = 2 # 2 = in survey area of both bands


    #Find all the tiles which have at least one observation validated.
    observed_tiles_b1=np.zeros(hp.nside2npix(32),dtype='int')
    observed_tiles_b2=np.zeros(hp.nside2npix(32),dtype='int')
    
    for tile in band1_tiles['tile_id'].data.astype(int):
        x=band1_tiles[band1_tiles['tile_id'] == str(tile)]['field1','field2','field3',
                                                       'field4','field5','field6','field7','field8']
        fields=[ y for y in np.lib.recfunctions.structured_to_unstructured(x.as_array())[0] if y != '' ]
        for field in fields:
            i=[ i for i,x in enumerate(band1_obs['name']) if field in x]
            if (band1_obs[i]['validated'][0] != '') & ('REJECTED' not in band1_obs[i]['validated'][0]):
                observed_tiles_b1[tile] = 1
                
    for tile in band2_tiles['tile_id'].data.astype(int):
        x=band2_tiles[band2_tiles['tile_id'] == str(tile)]['field1','field2','field3',
                                                       'field4','field5','field6','field7','field8']
        fields=[ y for y in np.lib.recfunctions.structured_to_unstructured(x.as_array())[0] if y != '' ]
        for field in fields:
            i=[ i for i,x in enumerate(band2_obs['name']) if field in x]
            if (band2_obs[i]['validated'][0] != '') & ('REJECTED' not in band2_obs[i]['validated'][0]):
                observed_tiles_b2[tile] = 2
                
    observed_tiles=observed_tiles_b1+observed_tiles_b2

    #Assign tiles based on if they're observed at least once per band:
    tile_map[(tile_map == 2) & (observed_tiles > 0) & (observed_tiles < 3)] = 3 #3 = Partially observed (1 of 2 bands)
    tile_map[(tile_map == 1) & (observed_tiles > 0) & (observed_tiles < 3)] = 4 #4 = Fully observed (1 band)
    tile_map[observed_tiles == 3] = 5 # 5 = Fully Observed (both bands)


    #Find all tiles that have their final version produced by AusSRC pipeline:
    mosaicked_tiles=np.zeros(hp.nside2npix(32),dtype='int')
    for tile in band1_tiles['tile_id'].data.astype(int):
        if band1_tiles[band1_tiles['tile_id'] == str(tile)]['aus_src'] != '':
            mosaicked_tiles[tile] += 1
    
    for tile in band2_tiles['tile_id'].data.astype(int):
        if band2_tiles[band2_tiles['tile_id'] == str(tile)]['aus_src'] != '':
            mosaicked_tiles[tile] += 2        

    #Assign tiles based on AusSRC completion (Note this isn't fully tested!)
    tile_map[((tile_map == 3) | (tile_map == 5) ) & (mosaicked_tiles > 0) & (mosaicked_tiles < 3)] = 6 #6 = One of 2 tiles produced
    tile_map[(tile_map == 4) & (mosaicked_tiles > 0) & (mosaicked_tiles < 3)] = 7 #7 = 1 of 1 tiles produced
    tile_map[mosaicked_tiles == 3] = 8 #8 = Both tiles produced


    col_dict={0:'black',
              1:"darkgrey",
              2:"lightgrey",
              3:"orange",
              4:"goldenrod",
              5:"yellow",
              6:"blue",
              7:"purple",
              8:"magenta"
                }
    labels = np.array(['Outside Survey',
                       f'Planned (1 band only) ({(tile_map==1).sum()})',
                       f'Planned (both bands) ({(tile_map==2).sum()})',
                       f'Partly observed (1 of 2 bands) ({(tile_map==3).sum()})',
                       f'Fully observed (1 band) ({(tile_map==4).sum()})',
                       f'Fully observed (both bands) ({(tile_map==5).sum()})',
                       f'Complete tile (1 of 2 bands) ({(tile_map==6).sum()})',
                       f'Complete tile (1 band) ({(tile_map==7).sum()})',
                       f'Complete tile (both bands) ({(tile_map==8).sum()})',
                      ])
    cm = mpl.colors.ListedColormap([col_dict[x] for x in col_dict.keys()])
    norm_bins = np.sort([*col_dict.keys()]) + 0.5
    norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)
    norm = mpl.colors.BoundaryNorm(norm_bins, len(labels), clip=True)
    fmt = mpl.ticker.FuncFormatter(lambda x, pos: labels[norm(x)])
    diff = norm_bins[1:] - norm_bins[:-1]
    tickz = norm_bins[:-1] + diff / 2
    
    
    
    data=hp.mollview(tile_map,title='Survey tile status',cmap=cm,min=0,max=8,return_projected_map=True,xsize=10000)
    plt.figure(figsize=(12,6))
    ax=plt.subplot(projection=ccrs.Mollweide())
    im=ax.imshow(data,transform=ccrs.Mollweide(),cmap=cm,norm=norm,extent=(-18040095.696147293, 18040095.696147293,-9020047.848073646, 9020047.848073646),origin='lower')
    ax.gridlines(draw_labels=False,xlocs=[0,45,90,135,180,-45,-90,-135,-179.99],ylocs=[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])
    plt.colorbar(im,format=fmt, ticks=tickz,shrink=0.7)
    for ra_grid in [0,3,6,9,12]:
        ax.text(-15*ra_grid,15,str(int(ra_grid))+'h',{'ha':'right','va':'top','color':'w'},transform=ccrs.PlateCarree())
    for ra_grid in [12.01,15,18,21]:
        ax.text(-15*ra_grid,15,str(int(ra_grid))+'h',{'ha':'left','va':'top','color':'w'},transform=ccrs.PlateCarree())
    for dec_grid in [75,60,45,30,15,0,-15,-30,-45,-60,-75]:
        ax.text(0,dec_grid,str(dec_grid)+'°',{'ha':'left','va':'bottom','color':'w'},transform=ccrs.PlateCarree())
    xmax=plt.xlim()[1]
    ymax=plt.ylim()[1]
    plt.text(xmax*-0.99,ymax*-0.95,'Equatorial\ncoordinates',{'size':14,'weight':'bold','ha':'left'})
    plt.text(xmax*-0.99,ymax*0.95,'Overlap status',{'size':14,'weight':'bold','ha':'left'})
    
    plt.savefig(basename+'overlap_equatorial.png',bbox_inches='tight',dpi=300)
    plt.close()

    
    data=hp.mollview(tile_map,title='Survey tile status',coord='CG',cbar=False,return_projected_map=True,xsize=10000)
    #hp.graticule(dmer=45,dpar=15)
    plt.figure(figsize=(12,6))
    ax=plt.subplot(projection=ccrs.Mollweide())
    im=ax.imshow(data,transform=ccrs.Mollweide(),cmap=cm,norm=norm,extent=(-18040095.696147293, 18040095.696147293,-9020047.848073646, 9020047.848073646),origin='lower')
    ax.gridlines(draw_labels=False,xlocs=[0,45,90,135,180,-45,-90,-135,-179.99],ylocs=[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90])
    plt.colorbar(im,format=fmt, ticks=tickz,shrink=0.7)
    for ra_grid in [0,45,90,135,180]:
        ax.text(-1*ra_grid,0,str(int(ra_grid))+'°',{'ha':'right','va':'top','color':'w'},transform=ccrs.PlateCarree())
    for ra_grid in [180.01,225,270,315]:
        ax.text(-1*ra_grid,0,str(int(ra_grid))+'°',{'ha':'left','va':'top','color':'w'},transform=ccrs.PlateCarree())
    for dec_grid in [75,60,45,30,15,0,-15,-30,-45,-60,-75]:
        ax.text(0,dec_grid,str(dec_grid)+'°',{'ha':'left','va':'bottom','color':'w'},transform=ccrs.PlateCarree())
    xmax=plt.xlim()[1]
    ymax=plt.ylim()[1]
    plt.text(xmax*-0.99,ymax*-0.95,'Galactic\ncoordinates',{'size':14,'weight':'bold','ha':'left'})
    plt.text(xmax*-0.99,ymax*0.95,'Overlap status',{'size':14,'weight':'bold','ha':'left'})

    plt.savefig(basename+'overlap_galactic.png',bbox_inches='tight',dpi=300)
    plt.close()
    
    
    #Generate and Aladin Lite view showing overlap:
    
    #Constants relating to field size and orientation.
    D = 3.151*u.deg
    PA1 = 51.79*u.deg
    PA2 = 128.21*u.deg
    
    #Compute the corner locations for all the fields.
    #Make list of all tile corners, and list of markers for all field centers.
    field_markers = []
    field_corners = []
    
    #Make band 1 boxes/markers:
    band1_validated=band1_obs[(np.char.startswith(band1_obs['validated'],'20')) & (band1_obs['validated'] != '')]
    for obs in band1_validated:
        name1 = obs['name']
        rot1 = (45+float(obs['rotation'] ))*u.deg 
        sbid = obs['sbid']
    
        sc=ac.SkyCoord(obs['ra_deg'],obs['dec_deg'],unit='deg',equinox='J2000')
    
        field_marker = [name1, sc.ra.deg, sc.dec.deg, rot1.value, sbid]
        field_markers.append( field_marker )
    
        tc1 = sc.directional_offset_by(position_angle = PA1+rot1, separation=D)
        tc2 = sc.directional_offset_by(position_angle = PA2+rot1, separation=D)
        tc3 = sc.directional_offset_by(position_angle = -PA2+rot1, separation=D)
        tc4 = sc.directional_offset_by(position_angle = -PA1+rot1, separation=D)
        field_corners.append( [name1,tc1.ra.deg,tc1.dec.deg,tc2.ra.deg,tc2.dec.deg,tc3.ra.deg,tc3.dec.deg,tc4.ra.deg,tc4.dec.deg] )
    
    band1_markers_aladin = []
    for (name,ra,dec,pa,sbid) in field_markers:
        band1_markers_aladin.append( f"band1_cat.addSources([A.marker({ra},{dec},{{popupTitle: '{name}', popupDesc: '<em>RA:</em> {ra:.3f}<br/><em>Dec:</em> {dec:.3f}<br/><em>PA:</em> {pa:.1f}<br/>(J2000)'}})]);\n" ) 
    
    band1_corners_aladin = []
    band1_corner_line =  'band1_overlay.addFootprints([A.polygon([{0}])]);\n'
    for field in field_corners:
        corner1 = ', '.join( [ '[{0},{1}]'.format(field[i*2+1],field[i*2+2]) for i in range(4) ] )
        band1_corners_aladin.append( band1_corner_line.format( corner1 ) ) 
    
    #Make band 2 boxes/markers:
    field_markers = []
    field_corners = []
    band2_validated=band2_obs[(np.char.startswith(band2_obs['validated'],'20')) & (band2_obs['validated'] != '')]
    for obs in band2_validated:
        name1 = obs['name']
        rot1 = (45+float(obs['rotation'] ))*u.deg 
        sbid = obs['sbid']
    
        sc=ac.SkyCoord(obs['ra_deg'],obs['dec_deg'],unit='deg',equinox='J2000')
    
        field_marker = [name1, sc.ra.deg, sc.dec.deg, rot1.value, sbid]
        field_markers.append( field_marker )
    
        tc1 = sc.directional_offset_by(position_angle = PA1+rot1, separation=D)
        tc2 = sc.directional_offset_by(position_angle = PA2+rot1, separation=D)
        tc3 = sc.directional_offset_by(position_angle = -PA2+rot1, separation=D)
        tc4 = sc.directional_offset_by(position_angle = -PA1+rot1, separation=D)
        field_corners.append( [name1,tc1.ra.deg,tc1.dec.deg,tc2.ra.deg,tc2.dec.deg,tc3.ra.deg,tc3.dec.deg,tc4.ra.deg,tc4.dec.deg] )
    
    # Make all tile markers
    band2_markers_aladin = []
    for (name,ra,dec,pa,sbid) in field_markers:
        band2_markers_aladin.append( f"band2_cat.addSources([A.marker({ra},{dec},{{popupTitle: '{name}', popupDesc: '<em>RA:</em> {ra:.3f}<br/><em>Dec:</em> {dec:.3f}<br/><em>PA:</em> {pa:.1f}<br/>(J2000)'}})]);\n" ) 
    
    #Make all field boxes
    band2_corners_aladin = []
    band2_corner_line =  'band2_overlay.addFootprints([A.polygon([{0}])]);\n'
    for field in field_corners:
        corner1 = ', '.join( [ '[{0},{1}]'.format(field[i*2+1],field[i*2+2]) for i in range(4) ] )
        band2_corners_aladin.append( band2_corner_line.format( corner1 ) ) 



    #Pick out tiles that are double-completed:
    tile_corners_aladin = []
    for tile_num in np.where(tile_map == 8)[0]:
        coords=np.array(hp.vec2ang(hp.boundaries(32,tile_num,step=1).T,lonlat=True)).T.tolist()
        corner = ','.join([f'[{ra},{dec}]' for ra,dec in iter(coords)])
        tile_corners_aladin.append('tile_overlay.addFootprints([A.polygon([{0}])]);\n'.format(corner))


    page_name = 'Band 1/2 Overlap regions'


    html_header = '''
    <!DOCTYPE html>
    <html>
    <head>
        <!-- Mandatory when setting up Aladin Lite v3 for a smartphones/tablet usage -->
        <meta name="viewport" content="width=device-width, height=device-height, initial-scale=1.0, user-scalable=no">
    </head>
    <body>
        <h1> POSSUM Status Monitor </h1>
        <h2> ''' + page_name + ''' </h2>
        <br><br>
    <!-- Aladin Lite has a dependency on the jQuery library -->
    <script src="https://code.jquery.com/jquery-1.10.1.min.js"></script>

    <!-- Aladin Lite container at requested dimensions -->
    <div id="aladin-lite-div" style="width:700px;height:600px;"></div>

    <!-- Aladin Lite JS code -->
    <script type="text/javascript" src="https://aladin.cds.unistra.fr/AladinLite/api/v3/latest/aladin.js" charset="utf-8"></script>

    <!-- Creation of Aladin Lite instance with initial parameters -->

    <script type="text/javascript">
        let aladin;
        A.init.then(() => {
            aladin = A.aladin('#aladin-lite-div', {fov: 240, projection: "SIN", target:'90 -90', cooFrame: 'equatorial', showReticle:false, showCooGridControl: true, showSimbadPointerControl: true, showCooGrid: false});
            let chipass = aladin.createImageSurvey("CHIPASS", "ov-gso/P/CHIPASS", "http://cade.irap.omp.eu/documents/Ancillary/4Aladin/CHIPASS/", "galactic", 3, {imgFormat: 'png'}); // setting a custom HiPS
            let emu = aladin.createImageSurvey('EMU', 'ASKAP-EMU Full Survey', 'http://www.aoc.nrao.edu/~jmarvil/emuHiPS', 'equatorial', 3, {imgFormat: 'png'}); // setting a custom HiPS
            aladin.setImageSurvey(aladin.createImageSurvey("RACS", "RACS-low-I", "https://casda.csiro.au/hips/RACS/low/I", "equatorial", 8, {imgFormat: 'png'})); // setting a custom HiPS
            let racs = aladin.getBaseImageLayer();
            racs.setColormap('grayscale',{stretch: 'linear', reversed: false});
            aladin.setOverlayImageLayer(chipass,'CHIPASS');
            let chimap = aladin.getOverlayImageLayer('CHIPASS');
            chimap.setOpacity(0.5);
            chimap.setColormap('grayscale',{stretch: 'linear', reversed: false});
            chimap.setBlendingConfig(true);
            racs.setBlendingConfig(true);

            var band1_overlay = A.graphicOverlay({color: 'blue', lineWidth: 1, name:'Band 1 Validated Fields'});
            aladin.addOverlay(band1_overlay);

            var band1_cat = A.catalog({name: 'Band 1 field markers', sourceSize: 10, color: '#eee'});
            aladin.addCatalog(band1_cat);
            band1_cat.hide();

            var band2_overlay = A.graphicOverlay({color: 'yellow', lineWidth: 1, name:'Band 2 Validated Fields'});
            aladin.addOverlay(band2_overlay);

            var band2_cat = A.catalog({name: 'Band 2 field markers', sourceSize: 10, color: '#eee'});
            aladin.addCatalog(band2_cat);
            band2_cat.hide();


            var tile_overlay = A.graphicOverlay({color: 'magenta', lineWidth: 1, name:'Both-band processed tiles'});
            aladin.addOverlay(tile_overlay);
            



    '''

    footer3 = '''
    });

</script>

        <h3> Built by Cameron Van Eck. Based on the EMU Survey Progress page by Josh Marvil </h3> 


</body>
</html>
    '''

    outfile= basename+'aladin_overlap.html'
    with open(outfile,'w') as out1:
        out1.writelines(html_header)
        out1.writelines(band1_corners_aladin)
        out1.writelines(band2_corners_aladin)
        out1.writelines(band1_markers_aladin)
        out1.writelines(band2_markers_aladin)
        out1.writelines(tile_corners_aladin)        
        out1.writelines(footer3)





def auto_update(ps,db_auth_file):
    """Update the sheet with the latest observation statuses, generate a new
    set of figures for all surveys, generate new Aladin Lite pages, and
    upload all results to CANFAR for use on the website.
    
    """
    print('Updating sheet.')
    update_observed_fields(ps,db_auth_file)

    print('Updating plots.')
#    create_plots(ps,'p1','./Pilot_band1_status_')
#    create_plots(ps,'p2','./Pilot_band2_status_')
    create_plots(ps,'1','./Survey_band1_status_')
    create_plots(ps,'2','./Survey_band2_status_')
    create_overlap_plots(ps, './Survey_')
    
    print('Updating Aladin overlays.')
    aladin_webpage(ps,'p1','aladin_pilot_band1.html')
    aladin_webpage(ps,'p2','aladin_pilot_band2.html')
    aladin_webpage(ps,'1','aladin_survey_band1.html')
    aladin_webpage(ps,'2','aladin_survey_band2.html')
    set_update_date_on_main_page('./status_page.html')


def cli():
    """Command line interface. Using flags, select the operation desired."""
    import argparse

    descStr = """
    Access the POSSUM status information. Select the flag appropriate for the
    desired operation. User must supply a valid Google Sheets API key.
    """

    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("Google_token_file", metavar="token.json",
                        help="JSON file containing Google API key.")
    parser.add_argument("db_auth_file", metavar="db_details.json",
                        help="JSON file containing AusSRC database details.")
    parser.add_argument("-s", dest="survey", default='1',type=str,
                        help="Specify survey and band: p1/p2 for pilot, 1/2 for full survey band 1/2. Default: 1")    
    parser.add_argument("-o", dest="update_obs", action="store_true",
                        help="Action: update observation sheet with new SBs.")
    parser.add_argument("-a", dest="aussrc_field", type=int,metavar='SB',
                        help="Action: update observation processed by AusSRC (supply SB#).")
    parser.add_argument("-t", dest="aussrc_tile", type=int,metavar='tile',
                        help="Action: update tile processed by AusSRC (supply tile#).")
    parser.add_argument("-p", dest='make_plots',type=str,metavar='FILENAME',
                        help="Action: create plots (supply filename/path)")
    parser.add_argument("-w", dest='make_webpage',type=str,metavar='aladin.html',
                        help="Action: create Aladin Lite webpage (supply filename/path)")
    parser.add_argument("-O", dest='overlap',type=str,metavar='FILENAME',
                        help="Action: create 2-band overlap plots and Aladin page.")
    
    parser.add_argument("-u", dest='auto_update',action="store_true",
                        help="Action: update for webpage (runs observation update, plots, and Aladin, outputting to current directory)")



    args = parser.parse_args()


    ps = open_sheet(args.Google_token_file)
    verify_sheet(ps)

    if args.update_obs:
        update_observed_fields(ps,args.db_auth_file)

    if args.aussrc_field is not None:
        update_aussrc_field_processed(ps,args.survey,args.aussrc_field)

    if args.aussrc_tile is not None:
        update_aussrc_tile_processed(ps,args.survey,args.aussrc_tile)

    if args.make_plots is not None:
        create_plots(ps, args.survey, args.make_plots)

    if args.make_webpage is not None:
        aladin_webpage(ps, args.survey,args.make_webpage)
        
    if args.overlap is not None:
        create_overlap_plots(ps, args.overlap)

    
    if args.auto_update is True:
        auto_update(ps,args.db_auth_file)






if __name__ == "__main__":
    cli()
