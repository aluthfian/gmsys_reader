import pandas as pd
import re
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from pyproj import Proj, Transformer
from scipy.optimize import minimize

# Reference
# GM-SYS User's Guide v4.9 at https://docplayer.net/19843267-Gm-sys-user-s-guide-version-4-9.html

# Remove number from list
def remove_numbers_from_list(input_list):
    """Remove numbers from a list of strings."""
    return [word for word in input_list if not any(char.isdigit() for char in word)]

# Function to parse BLK file and extract data for a specific layer
def parse_blk_file(filename):
    data = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        for i in range(len(lines)):
            line = lines[i].strip('\n')
            try:
                body_index = int(line.split()[0])
            except:
                if re.match(r'^[A-Za-z]', line):
                    model_units = line
                    density = 1e+3*float(lines[i+1].split()[0])
                    inclination = float(lines[i+1].split()[1])
                    declination = float(lines[i+1].split()[2])
                    magsus = np.round(4*np.pi*float(lines[i+1].split()[3]),4)
                    remmag = 1e+3*float(lines[i+1].split()[4])
                    color_hex = lines[i+2].split()[4].split('/')[0][-6:]  # Extract last 6 characters after '#'
                    color_rgb = mcolors.hex2color('#' + color_hex)  # Convert hex color to RGB tuple
                    data.append([body_index, model_units, density, inclination, declination, magsus, remmag, color_rgb])
                    i += 3
    
    df = pd.DataFrame(data, columns=['Body index', 'Model units', 'Density (kg/m^3)', 'Inclination (deg)', 'Declination (deg)', 'Magsus (SI)', 'RemMag (A/m)', 'Colour (RGB)'])
    df.set_index('Body index', inplace=True)
    return df

# Function to parse SUR file and extract data for all specific layer
def parse_sur_file(filename):
    data = []
    polygon_dict={}
    polygon_dummy={}
    how_much_pts=0
    how_much_skip=0
    index_unit=[]
    with open(filename, 'r') as file:
        lines = file.readlines()
        for i in range(2,len(lines)):
            if (i!=0)&(i<how_much_skip+1):
                continue
            else:
                xy_vertex_now = []
                sfc_index = int(lines[i])
                index_unit_left_block = lines[i+1].split()[0]
                index_unit_right_block = lines[i+1].split()[1]
                index_unit.append(index_unit_left_block)
                index_unit.append(index_unit_right_block)
                how_much_pts = int(lines[i+4].strip())
                how_much_skip = i+4+how_much_pts
                #x,z coordinate
                for k in range(how_much_pts):
                    xy_vertex_now.append([float(coord) for coord in lines[i+5+k].split()])
                xy_vertex_now = 1e+3*np.array(xy_vertex_now)
                xy_vertex_now[:,1] = -1*xy_vertex_now[:,1]
                xy_vertex_flip = np.flipud(xy_vertex_now).copy()
                try:
                    vstack1=[]
                    vstack_dummy=[]
                    topbottom = np.linalg.norm(xy_vertex_now[0] - polygon_dict[index_unit_left_block][-1])
                    toptop = np.linalg.norm(xy_vertex_now[0] - polygon_dict[index_unit_left_block][0])
                    bottomtop = np.linalg.norm(xy_vertex_now[-1] - polygon_dict[index_unit_left_block][0])
                    bottombottom =np.linalg.norm(xy_vertex_now[-1] - polygon_dict[index_unit_left_block][-1])
                    min_value = np.array([topbottom, toptop, bottomtop, bottombottom]).min()
                    if min_value != 0:
                        try:
                            topbottom_dummy = np.linalg.norm(xy_vertex_now[0] - polygon_dummy[index_unit_left_block][-1])
                            toptop_dummy = np.linalg.norm(xy_vertex_now[0] - polygon_dummy[index_unit_left_block][0])
                            bottomtop_dummy = np.linalg.norm(xy_vertex_now[-1] - polygon_dummy[index_unit_left_block][0])
                            bottombottom_dummy =np.linalg.norm(xy_vertex_now[-1] - polygon_dummy[index_unit_left_block][-1])
                            min_value_dummy = np.array([topbottom_dummy, toptop_dummy, bottomtop_dummy, bottombottom_dummy]).min()
                            if min_value_dummy == topbottom_dummy:
                                vstack_dummy = np.vstack([polygon_dummy[index_unit_left_block],xy_vertex_now])
                            elif min_value_dummy == toptop_dummy:
                                vstack_dummy = np.vstack([xy_vertex_flip,polygon_dummy[index_unit_left_block]])
                            elif min_value_dummy == bottomtop_dummy:
                                vstack_dummy = np.vstack([xy_vertex_now,polygon_dummy[index_unit_left_block]])
                            elif min_value_dummy == bottombottom_dummy:
                                vstack_dummy = np.vstack([polygon_dummy[index_unit_left_block],xy_vertex_flip])
                            if np.all(np.in1d(polygon_dummy[index_unit_left_block], vstack1)):
                                pass
                            polygon_dummy[index_unit_left_block]=vstack_dummy
                        except:
                            polygon_dummy[index_unit_left_block]=xy_vertex_now
                    else:
                        if min_value == topbottom:
                            vstack1 = np.vstack([polygon_dict[index_unit_left_block],xy_vertex_now])
                        elif min_value == toptop:
                            vstack1 = np.vstack([xy_vertex_flip,polygon_dict[index_unit_left_block]])
                        elif min_value == bottomtop:
                            vstack1 = np.vstack([xy_vertex_now,polygon_dict[index_unit_left_block]])
                        elif min_value == bottombottom:
                            vstack1 = np.vstack([polygon_dict[index_unit_left_block],xy_vertex_flip])
                        if np.all(np.in1d(polygon_dict[index_unit_left_block], vstack1)):
                            pass
                        polygon_dict[index_unit_left_block]=vstack1
                except:
                    polygon_dict[index_unit_left_block]=xy_vertex_now
                    
                try:
                    vstack2=[]
                    vstack2_dummy=[]
                    topbottom = np.linalg.norm(xy_vertex_flip[0] - polygon_dict[index_unit_right_block][-1])
                    toptop = np.linalg.norm(xy_vertex_flip[0] - polygon_dict[index_unit_right_block][0])
                    bottomtop = np.linalg.norm(xy_vertex_flip[-1] - polygon_dict[index_unit_right_block][0])
                    bottombottom =np.linalg.norm(xy_vertex_flip[-1] - polygon_dict[index_unit_right_block][-1])
                    min_value = np.array([topbottom, toptop, bottomtop, bottombottom]).min()
                    if min_value != 0:
                        try:
                            topbottom_dummy = np.linalg.norm(xy_vertex_flip[0] - polygon_dummy[index_unit_right_block][-1])
                            toptop_dummy = np.linalg.norm(xy_vertex_flip[0] - polygon_dummy[index_unit_right_block][0])
                            bottomtop_dummy = np.linalg.norm(xy_vertex_flip[-1] - polygon_dummy[index_unit_right_block][0])
                            bottombottom_dummy =np.linalg.norm(xy_vertex_flip[-1] - polygon_dummy[index_unit_right_block][-1])
                            min_value_dummy = np.array([topbottom_dummy, toptop_dummy, bottomtop_dummy, bottombottom_dummy]).min()
                            if min_value_dummy == topbottom_dummy:
                                vstack2_dummy = np.vstack([polygon_dummy[index_unit_right_block],xy_vertex_flip])
                            elif min_value_dummy == toptop_dummy:
                                vstack2_dummy = np.vstack([xy_vertex_now,polygon_dummy[index_unit_right_block]])
                            elif min_value_dummy == bottomtop_dummy:
                                vstack2_dummy = np.vstack([xy_vertex_flip,polygon_dummy[index_unit_right_block]])
                            elif min_value_dummy == bottombottom_dummy:
                                vstack2_dummy = np.vstack([polygon_dummy[index_unit_right_block],xy_vertex_now])
                            if np.all(np.in1d(polygon_dummy[index_unit_right_block], vstack1)):
                                pass
                            polygon_dummy[index_unit_right_block]=vstack2_dummy
                        except:
                            polygon_dummy[index_unit_right_block]=xy_vertex_flip
                    else:
                        if min_value == topbottom:
                            vstack2 = np.vstack([polygon_dict[index_unit_right_block],xy_vertex_flip])
                        elif min_value == toptop:
                            vstack2 = np.vstack([xy_vertex_now,polygon_dict[index_unit_right_block]])
                        elif min_value == bottomtop:
                            vstack2 = np.vstack([xy_vertex_flip,polygon_dict[index_unit_right_block]])
                        elif min_value == bottombottom:
                            vstack2 = np.vstack([polygon_dict[index_unit_right_block],xy_vertex_now])
                        if np.all(np.in1d(polygon_dict[index_unit_right_block], vstack2)):
                            pass
                        polygon_dict[index_unit_right_block]=vstack2
                except:
                    polygon_dict[index_unit_right_block]=xy_vertex_flip
   
    index_unit = np.unique(index_unit)
    for units in polygon_dummy.keys():
        vstack_v3 = []
        current_xy = polygon_dict[units]
        dummy_xy = polygon_dummy[units]
        dummy_xy_reverse = np.flipud(dummy_xy)
        topbottom_v3 = np.linalg.norm(current_xy[0] - dummy_xy[-1])
        toptop_v3 = np.linalg.norm(current_xy[0] - dummy_xy[0])
        bottomtop_v3 = np.linalg.norm(current_xy[-1] - dummy_xy[0])
        bottombottom_v3 =np.linalg.norm(current_xy[-1] - dummy_xy[-1])
        min_value_v3 = np.array([topbottom_v3, toptop_v3, bottomtop_v3, bottombottom_v3]).min()
        if min_value_v3 == topbottom_v3:
            vstack_v3 = np.vstack([dummy_xy, current_xy])
        elif min_value_v3 == toptop_v3:
            vstack_v3 = np.vstack([dummy_xy_reverse, current_xy])
        elif min_value_v3 == bottomtop_v3:
            vstack_v3 = np.vstack([current_xy, dummy_xy])
        elif min_value_v3 == bottombottom_v3:
            vstack_v3 = np.vstack([current_xy, dummy_xy_reverse])
        polygon_dict[units] = vstack_v3
    
    
    model_geom = pd.DataFrame(index=[int(unit) for unit in index_unit], columns=['Geometry (Xm,Ym)'])
    model_geom.index.name='Body index'
    for units in index_unit:
        model_geom.at[int(units), 'Geometry (Xm,Ym)']=[vertex.tolist() for vertex in polygon_dict[str(units)]]
    return model_geom

# Function to parse the coordinate system of the model
def parse_ecs_file(filename):
    proj_str = ''
    with open(filename, 'r') as file:
        file_lines = file.readlines()
        for line in range(len(file_lines)):
            if 'Origin' in file_lines[line]:
                text_to_split = file_lines[line].replace('=',' ').replace('\n','')
                line_comp = text_to_split.split()
                model_x0 = float(line_comp[1])
                model_y0 = float(line_comp[2])
            elif 'CosTheta' in file_lines[line]:
                text_to_split = file_lines[line].replace('=',' ').replace('\n','')
                line_comp = text_to_split.split()
                azimuth = 90-np.rad2deg(np.arccos(float(line_comp[1])))
            if 'Projection=' in file_lines[line]:
                text_to_split = file_lines[line].replace('=',';').replace('*','').replace(' / ',';').replace('\n','')
                line_comp = text_to_split.split(';')
                projection_name = line_comp[2].split()
                projection_name = remove_numbers_from_list(projection_name)
                projection_abbreviation = ''.join([word[0].lower() for word in projection_name])
                if 'tm' in projection_abbreviation:
                    projection_abbreviation = 'tmerc'
                proj_str += '+proj='+projection_abbreviation
            if 'Datum' in file_lines[line]:
                #example Datum=NZGD49,6378388,0.0819918899790298,0
                text_to_split = file_lines[line].replace('=',',').replace('\n','')
                line_comp = text_to_split.split(',')
                proj_str += ' +a='+line_comp[2]+' +e='+line_comp[3]
            if 'Method' in file_lines[line]:
                #example Method="New Zealand Map Grid",-41,173,2510000,6023150
                text_to_split = file_lines[line].replace('=',',').replace('\n','')
                line_comp = text_to_split.split(',')
                if len(line_comp)==6:
                    proj_str += ' +lat_0='+line_comp[2]+' +lon_0='+line_comp[3]+' +x_0='+line_comp[4]+' +y_0='+line_comp[5]
                elif len(line_comp)==7:
                    proj_str += ' +lat_0='+line_comp[2]+' +lon_0='+line_comp[3]+' +k='+line_comp[4]+' +x_0='+line_comp[5]+' +y_0='+line_comp[6]
            if 'Unit' in file_lines[line]:
                #example Unit=m,1
                text_to_split = file_lines[line].replace('=',',').replace('\n','')
                line_comp = text_to_split.split(',')
                proj_str += ' +units='+line_comp[1]
    
    crs_of_model = Proj(proj_str, preserve_units=False)
    return model_x0, model_y0, azimuth, crs_of_model

# Function to parse the Well file
def well_file_reader(filename):
    WellData = pd.DataFrame(columns=['Well name','Profile X (m)','Elevation (m)','Units','Depths (m)'])
    with open(filename, 'r') as file:
        file_lines = file.readlines()
        file_lines = [everyline.replace('\n','') for everyline in file_lines]
        idx_now = 0
        for line in range(len(file_lines)):
            try:
                split_line = file_lines[line].split()
                to_float = [float(numbers) for numbers in split_line]
                if len(to_float) ==2:
                    WellData.loc[idx_now,'Well name']=file_lines[line-2]
                    x_dist = round(1e+3*float(file_lines[line-1].split()[0]),2)
                    y_dist = round(-1e+3*to_float[1],2)
                    WellData.loc[idx_now,'Profile X (m)']=round(np.sign(x_dist)*np.sqrt(x_dist**2 + y_dist**2),2)
                    well_z0 = round(-1e+3*to_float[0],2)
                    WellData.loc[idx_now,'Elevation (m)']=well_z0
                    num_of_unit = int(file_lines[line+1])
                    WellData.at[idx_now,'Units']=file_lines[line+3:line+3+(num_of_unit)*2:2]
                    WellData.at[idx_now,'Depths (m)']=[round(well_z0+(1e+3*float(depth)),2) for depth in file_lines[line+2:line+2+(num_of_unit)*2:2]]
                    idx_now += 1
            except:
                pass
    return WellData

# reading Gravity Data file
def grav_file_reader(filename):
    offset_status =''
    offset_val=np.nan
    GravData = pd.DataFrame(columns=['Profile Dist (m)','Offset from Profile (m)','Elevation (m)','Observed (mGal)','Modelled (mGal)'])
    with open(filename, 'r') as file:
        file_lines = file.readlines()
        offset_info = file_lines[0].split()
        if int(offset_info[3])==0:
            if int(offset_info[0])==0:
                absolute_offset = float(offset_info[2])
                offset_status += 'manual'
                offset_val = absolute_offset
            elif int(offset_info[2])==0:
                idx_sta_for_offset = int(offset_info[0])
                offset_status += 'by station'
                offset_val = idx_sta_for_offset
                data_at_sta = np.array([float(value) for value in file_lines[idx_sta_for_offset+1].split()])
                obsGrav_at_sta = data_at_sta[2]
                modGrav_at_sta = data_at_sta[3]
        elif int(offset_info[3])==1:
            autoOffset = True
            offset_status += 'automatic'
            # Initial guess for the offset
            init_offset_val = 0.0
            # Perform optimization
            data_values = pd.read_csv(filename, delimiter='\s+', skiprows=2).to_numpy()
            objective_function= lambda x: np.sum((data_values[:,-2] - (data_values[:,-3] + x))**2)
            optimising_offset = minimize(objective_function, init_offset_val)
            # Get the optimized offset value
            offset_val = optimising_offset.x[0]
        how_many_pts = int(file_lines[1])
        for idx in range(how_many_pts):
            values_inplace = np.array([float(value) for value in file_lines[idx+2].split()])
            values_inplace[:2] = 1e+3*values_inplace[:2] #Data location in the profile distance, and data elevation, positive downward
            values_inplace[-1] = 1e+3*values_inplace[-1] #Data offset at the direction perpendicular to the profile
            # Move the data offset from the original position
            values_inplace = np.insert(values_inplace, 1, values_inplace[-1]) #Insertion before idx=1
            values_inplace = np.delete(values_inplace, -1) #Deletion of redundant data offset at the end of the array
            # Elevation data position moved from (1) to (2)
            values_inplace[2] = -1*values_inplace[2]
            # Gravity data position moved from (-3:-2) to (-2:-1)
            if (offset_status == 'manual')|(offset_status == 'automatic'):
                if offset_val>=0:
                    values_inplace[-1] = round(values_inplace[-1]-offset_val,5)
                else:
                    values_inplace[-1] = round(offset_val+values_inplace[-1],5)
            elif offset_status == 'by station':
                #data@sta - (model@sta - model@pts)
                values_inplace[-1] = round(obsGrav_at_sta - (modGrav_at_sta - values_inplace[-1]),5)
            GravData.loc[idx] = values_inplace
    return GravData, offset_status, offset_val

# Reading Magnetic Data File
def mag_file_reader(filename):
    offset_status =''
    offset_val=np.nan
    MagData = pd.DataFrame(columns=['Profile Dist (m)','Offset from Profile (m)','Elevation (m)','Observed (nT)','Modelled (nT)'])
    with open(filename, 'r') as file:
        file_lines = file.readlines()
        offset_info = file_lines[0].split()
        if int(offset_info[3])==0:
            if int(offset_info[0])==0:
                absolute_offset = float(offset_info[2])
                offset_status += 'manual'
                offset_val = absolute_offset
            elif int(offset_info[2])==0:
                idx_sta_for_offset = int(offset_info[0])
                offset_status += 'by station'
                offset_val = idx_sta_for_offset
                data_at_sta = np.array([float(value) for value in file_lines[idx_sta_for_offset+2].split()])
                obsMag_at_sta = data_at_sta[2]
                modMag_at_sta = data_at_sta[3]
        elif int(offset_info[3])==1:
            autoOffset = True
            offset_status += 'automatic'
            # Initial guess for the offset
            init_offset_val = 0.0
            # Perform optimization
            data_values = pd.read_csv(filename, delimiter='\s+', skiprows=3).to_numpy()
            objective_function= lambda x: np.sum((data_values[:,-2] - (data_values[:,-3] + x))**2)
            optimising_offset = minimize(objective_function, init_offset_val)
            # Get the optimized offset value
            offset_val = optimising_offset.x[0]
        how_many_pts = int(file_lines[1])
        ambient_mag_param = file_lines[2].split()
        inclination = float(ambient_mag_param[0])
        declination = float(ambient_mag_param[1])
        mag_field = float(ambient_mag_param[2])
        ambient_mag_props = {'inc_deg':inclination, 'dec_deg':declination, 'mag_nT':mag_field}
        for idx in range(how_many_pts):
            values_inplace = np.array([float(value) for value in file_lines[idx+3].split()])
            values_inplace[:2] = 1e+3*values_inplace[:2]
            values_inplace[-1] = 1e+3*values_inplace[-1] #Data offset at the direction perpendicular to the profile
            # Move the data offset from the original position
            values_inplace = np.insert(values_inplace, 1, values_inplace[-1]) #Insertion before idx=1
            values_inplace = np.delete(values_inplace, -1) #Deletion of redundant data offset at the end of the array
            # Elevation data position moved from (1) to (2)
            values_inplace[2] = -1*values_inplace[2]
            # Gravity data position moved from (-3:-2) to (-2:-1)
            if (offset_status == 'manual')|(offset_status == 'automatic'):
                if offset_val>=0:
                    values_inplace[-1] = round(values_inplace[-1]-offset_val,5)
                else:
                    values_inplace[-1] = round(offset_val+values_inplace[-1],5)
            elif offset_status == 'by station':
                values_inplace[-1] = round(obsMag_at_sta - (modMag_at_sta - values_inplace[-1]),5)
            MagData.loc[idx] = values_inplace
    return MagData, offset_status, offset_val, ambient_mag_props
