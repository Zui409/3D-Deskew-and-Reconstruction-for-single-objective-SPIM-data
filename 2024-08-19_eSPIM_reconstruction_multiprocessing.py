# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 13:19:54 2024

@author: zz409
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 15:42:08 2024

@author: zz409
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 13:01:18 2023

@author: zz409
"""



from skimage import io
import numpy as np
from scipy import ndimage
from scipy.ndimage import gaussian_filter1d
import tifffile
import multiprocessing
from pathos.multiprocessing import ProcessingPool as Pool
from functools import partial
import os 
import re
import imagej

# Main directory of the data where your positions in (eg. folder where the Pos0, Pos1 belong)
path = r'D:\Folder'

# Where your fiji is. For cropping & projection 
path_fiji = r'D:\Fiji\Fiji.app'              

# What channel do you want to process?
channel = ['638','488','561']               #  488 561 or 638

# Have you done these processing steps?
renamed = False                    # Rename and reorganise the file
cropped = False                    # Fiji cropping + subtract background + reverse
processed = False                  # Deskew
projectioned = False               # Fiji 3D projection
combined = False                   # Combine frames into videos

# ROI rotation angle for multicolour channel matching. 
angle_638 = -5.5
angle_561 = -0.2     
angle_488 = 0

# ROI selection in Fiji (can be retreieved using the Record function in FIJI)
# No binnning
#ROI_638 = "211, 815, 770, 310"                 #They MUST be in string formats
#ROI_561 = "212, 430, 770, 310"
#ROI_488 = "206, 42, 770, 310"

# 2x2 binning
ROI_638 = "0, 0, 212, 155"   
ROI_561 = "0, 62, 242, 155"
ROI_488 = "10, 0, 212, 164"

# Rolling ball background subtraction
sub_background = "50"

# Parameter for deskewing
angle = 30
spacing_factor = 2.29           # no binning 6.5

blur = False       
blur_sigma = 3

# Parameter for 3D projection:
direction = "Y"       # X or Y. Must be CAPITAL letter
angle_increment = 10

# Parameter for combining frames
framenum = [19, 28]

# Normally no need to change
# For Multi/Parallel Processing
param = [angle, spacing_factor, blur, blur_sigma]
num_cpu = multiprocessing.cpu_count()-2          
num_workers = num_cpu

#%% Rename

# Find all directories of the positions in the data_path folder in a list
def findPos(path):
    pos_path_list = [f.path for f in os.scandir(path) if f.is_dir()]
    return pos_path_list

# Find all directories of the slices in a list
def generate_slice(pos_path):
    cycle_path_list = [f.path for f in os.scandir(pos_path) if f.is_dir()]
    
    return cycle_path_list

# Find all directories of the TIFF files in a list
def getListFiles(path):
    filelist = [] 
    for root, dirs, files in os.walk(path):  
        for filespath in files: 
            if filespath[-4:] == '.tif':
                filelist.append(os.path.join(root,filespath)) 
    return filelist

# Rename and save images
def concat_img(filelist,con_name):
    num_imgs = len(filelist)
    for i in range(0, num_imgs):
        im = io.imread(filelist[i])
        tifffile.imsave(con_name, im)

def rename_img(main_path, img_name, i):
    im_newname = main_path + r'\\' + str(i) + '.tif'
    if not os.path.exists(im_newname):
        os.rename(img_name[0], im_newname)
    return im_newname
        
def sort_order(a_list):
    index_list = []
    for element in a_list:
        index = re.findall(r'\d+', element)
        index_list.append(int(index[-1]))
    sorted_list = [x for _,x in sorted(zip(index_list, a_list))] 
    return sorted_list

def rename(path):
    pos_path_list = findPos(path)
    for main_path in pos_path_list:
        slice_list = generate_slice(main_path)
        for i in range(0,len(slice_list)):
            #con_name = main_path + '/' + str(i) +'.tif'
            img_name = getListFiles(slice_list[i])
    return 0

def get_img_list(path):
    pos_path_list = findPos(path)
    pos_name_list = []
    for main_path in pos_path_list:
        img_name = getListFiles(main_path)
        sorted_name_list = sort_order(img_name)
        pos_name_list.append(sorted_name_list)
    return pos_name_list

    
#%% Create analysis folder

# Create sub folder(s) with a given name
def create_sub_folder(path, new_folder_name):
    if isinstance(path, list):
        new_dir_list = []
        for element in path:
            new_dir_name = element + '/' + new_folder_name
            new_dir_list.append(new_dir_name)
            if not os.path.exists(new_dir_name):
                os.makedirs(new_dir_name)
    elif isinstance(path, str):
        new_dir_name = path + '\\' + new_folder_name
        if not os.path.exists(new_dir_name):
            os.makedirs(new_dir_name)
        new_dir_list = [new_dir_name]
    else:
        print("Error: Invalid input path.")
    return new_dir_list

def get_last_name(path):
    if isinstance(path, list):
        name_list = []
        for path_k in path:
            name_k = path_k.split("\\")[-1]
            name_list.append(name_k)
    elif isinstance(path, str):
        name_list = [path.split("\\")[-1]]
    return name_list

# Get a list of direct subfolder of the path folder
def get_sub_folder(path):
    sub_folder_list = [f.path for f in os.scandir(path) if f.is_dir()]
    return sub_folder_list

def create_analysis_folder(path, channel):
    pos_list = get_sub_folder(path)
    pos_name = get_last_name(pos_list)
    ana_pos_list = []
    for i in range(len(pos_list)):
        new_analysis_path = path + "_Analysis" + "/" + pos_name[i]
        ana_pos_list.append(new_analysis_path)
        if not os.path.exists(new_analysis_path):
            os.makedirs(new_analysis_path)
    if '488' in channel:
        folder_488 = create_sub_folder(ana_pos_list, '488')
        folder_488_process = create_sub_folder(ana_pos_list, '488_Processed')
        folder_488_projection = create_sub_folder(ana_pos_list, '488_Projection')
    else: 
        folder_488 = []
        folder_488_process = []
        folder_488_projection = []
    if '561' in channel:
        folder_561 = create_sub_folder(ana_pos_list, '561')
        folder_561_process = create_sub_folder(ana_pos_list, '561_Processed')
        folder_561_projection = create_sub_folder(ana_pos_list, '561_Projection')
    else: 
        folder_561 = []
        folder_561_process = []
        folder_561_projection = []
    if '638' in channel:
        folder_638 = create_sub_folder(ana_pos_list, '638')
        folder_638_process = create_sub_folder(ana_pos_list, '638_Processed')
        folder_638_projection = create_sub_folder(ana_pos_list, '638_Projection')
    else: 
        folder_638 = []
        folder_638_process = []
        folder_638_projection = []
    return ana_pos_list, folder_488, folder_561, folder_638, folder_488_process, folder_488_projection, folder_561_process, folder_561_projection, folder_638_process, folder_638_projection

#%% Crop images in Fiji

def crop_images(path_fiji, pos_name_list, channel, angle_488, angle_561, angle_638, ROI_488, ROI_561, ROI_638, sub_background, folder_488, folder_561, folder_638):
    IJ = imagej.init(path_fiji, mode = 'interactive')
    IJ.ui().showUI()
    for i in range(len(pos_name_list)):
        for j in range(len(pos_name_list[i])):
            image = pos_name_list[i][j]
            image = image.replace("\\", "/")
            #imp =  IJ.io().open(image)     
            macro_new1 = '''
            open("''' +image + '''");
            
            
            '''
            IJ.py.run_macro(macro_new1)
            if '488' in channel:
                process_488 = folder_488[i] + "/" + str(j) + ".tif"
                process_488 = process_488.replace("\\", "/")
                macro_new = '''
                selectWindow("''' +image.split('/')[-1] +'''");
                run("Rotate... ", "angle='''+ str(angle_488) + ''' grid=1 interpolation=Bilinear stack");
                makeRectangle(''' + ROI_488 + ''');
                run("Duplicate...", "duplicate");      
                run("Subtract Background...", "rolling='''+ sub_background + ''' stack");
                run("Reverse");
                saveAs("Tiff", "'''+ process_488 + '''");
                close();
                '''
                IJ.py.run_macro(macro_new)
                
            
                    
            if '638' in channel:      
                process_638 = folder_638[i] + "/" + str(j) + ".tif"
                process_638 = process_638.replace("\\", "/")
                macro_new = '''
                selectWindow("''' +image.split('/')[-1] +'''");
                run("Rotate... ", "angle='''+ str(angle_638) + ''' grid=1 interpolation=Bilinear stack");
                makeRectangle(''' + ROI_638 + ''');
                run("Duplicate...", "duplicate");
                run("Subtract Background...", "rolling=50 stack");
                run("Reverse");
                saveAs("Tiff", "'''+ process_638 + '''");
                close();
                '''
                IJ.py.run_macro(macro_new)
                
            if '561' in channel:     
                process_561 = folder_561[i] + "/" + str(j) + ".tif"
                process_561 = process_561.replace("\\", "/")
                macro_new_561 = '''
                selectWindow("''' +image.split('/')[-1] +'''");
                run("Select All");
                run("Rotate... ", "angle='''+ str(angle_561) + ''' grid=1 interpolation=Bilinear stack");
                run("Scale...", "x=1.017 y=1.017 interpolation=Bilinear average process create");
                selectWindow("''' +image.split('/')[-1][:-4] +'''-1.tif");
                makeRectangle(''' + ROI_561 + ''');
                run("Duplicate...", "duplicate");
                run("Subtract Background...", "rolling='''+ sub_background + ''' stack");
                run("Reverse");
                saveAs("Tiff", "'''+ process_561 + '''");
                close();
                close();
                '''
                IJ.py.run_macro(macro_new_561)
                
            macro_new2 = '''
            close();
            '''
            IJ.py.run_macro(macro_new2)
    #del IJ
    return IJ

def fake_ij(cropped):
    if cropped == True:
        IJ = "No"
    return IJ
            
#%% Read cropped images

def read_cropped_images(ana_pos_list, channel):    
    if '488' in channel:
        SB_488_list = [] 
        for main_path in ana_pos_list:
            folder_488 = main_path + r'/488'
            img_name = getListFiles(folder_488)
            sorted_name_list = sort_order(img_name)
            SB_488_list.append(sorted_name_list)
    else:
        SB_488_list = [] 
        
        
    if '561' in channel:
        SB_561_list = [] 
        
        for main_path in ana_pos_list:
            folder_561 = main_path + r'/561'
            img_name = getListFiles(folder_561)
            sorted_name_list = sort_order(img_name)
            SB_561_list.append(sorted_name_list)
    else:
        SB_561_list = [] 
        
    if '638' in channel:
        SB_638_list = [] 
        for main_path in ana_pos_list:
            folder_638 = main_path + r'/638'
            img_name = getListFiles(folder_638)
            sorted_name_list = sort_order(img_name)
            SB_638_list.append(sorted_name_list)
    else:
        SB_638_list = [] 
    return SB_488_list, SB_561_list, SB_638_list


def sort_order_third(a_list):
    index_list = []
    for element in a_list:
        index = re.findall(r'\d+', element)
        index_list.append(int(index[-4]))
    sorted_list = [x for _,x in sorted(zip(index_list, a_list))] 
    return sorted_list

def read_processed_images(ana_pos_list, channel):    
    if '488' in channel:
        process_488_list = [] 
        process_488_last_list = []
        for main_path in ana_pos_list:
            folder_488 = main_path + r'/488_Processed'
            img_name = getListFiles(folder_488)
            sorted_name_list = sort_order_third(img_name)
            last_img = sorted_name_list[0]
            process_488_last_list.append(last_img)
            process_488_list.append(sorted_name_list)
    else:
        process_488_list = [] 
        process_488_last_list = []
        
    if '561' in channel:
        process_561_list = [] 
        process_561_last_list = []
        for main_path in ana_pos_list:
            folder_561 = main_path + r'/561_Processed'
            img_name = getListFiles(folder_561)
            sorted_name_list = sort_order_third(img_name)
            last_img = sorted_name_list[0]
            process_561_last_list.append(last_img)
            process_561_list.append(sorted_name_list)
    else:
        process_561_list = [] 
        process_561_last_list = []
        
            
    if '638' in channel:
        process_638_list = [] 
        process_638_last_list = []
        for main_path in ana_pos_list:
            folder_638 = main_path + r'/638_Processed'
            img_name = getListFiles(folder_638)
            sorted_name_list = sort_order_third(img_name)
            last_img = sorted_name_list[0]
            process_638_last_list.append(last_img)
            process_638_list.append(sorted_name_list)
    else:
        process_638_list = [] 
        process_638_last_list = []
    return process_488_list, process_561_list, process_638_list, process_488_last_list, process_561_last_list, process_638_last_list

def read_projection_images(ana_pos_list, channel):    
    if '488' in channel:
        projection_488_list = [] 
        for main_path in ana_pos_list:
            folder_488 = main_path + r'/488_Projection'
            img_name = getListFiles(folder_488)
            sorted_name_list = sort_order(img_name)
            projection_488_list.append(sorted_name_list)
    else:
        projection_488_list = [] 
        
    if '561' in channel:
        projection_561_list = [] 
        for main_path in ana_pos_list:
            folder_561 = main_path + r'/561_Projection'
            img_name = getListFiles(folder_561)
            sorted_name_list = sort_order(img_name)
            projection_561_list.append(sorted_name_list)
    else:
        projection_561_list = [] 
        
            
    if '638' in channel:
        projection_638_list = [] 
        for main_path in ana_pos_list:
            folder_638 = main_path + r'/638_Projection'
            img_name = getListFiles(folder_638)
            sorted_name_list = sort_order(img_name)
            projection_638_list.append(sorted_name_list)
    else:
        projection_638_list = [] 
    return projection_488_list, projection_561_list, projection_638_list


def get_intensity_range(channel, process_488_list, process_561_list, process_638_list, process_488_last_list, process_561_last_list, process_638_last_list):
    if '488' in channel:
        im_range_488 = []
        for imname in process_488_last_list:
            img = io.imread(imname)
            im_max = int(np.max(img))
            im_min = np.min(img)
            range_pair = "(" + str(im_min) +", "+ str(im_max) + ")"
            im_range_488.append(range_pair)
        # range_488 = process_488_list
        # for i in range(len(process_488_list)):
        #     for k in range(len(process_488_list[i])):
        #         range_488[i][k] = im_range_488[i]        
    else:
        im_range_488 = []
            
    if '561' in channel:
        im_range_561 = []
        for imname in process_561_last_list:
            img = io.imread(imname)
            im_max = int(np.max(img)*0.75)
            im_min = np.min(img)
            range_pair = "(" + str(im_min) +", "+ str(im_max) + ")"
            im_range_561.append(range_pair)
        # range_561 = process_561_list
        # for i in range(len(process_561_list)):
        #     for k in range(len(process_561_list[i])):
        #         range_561[i][k] = im_range_561[i]
    else:
        im_range_561 = []
            
    if '638' in channel:
        im_range_638 = []
        for imname in process_638_last_list:
            img = io.imread(imname)
            im_max = int(np.max(img))
            im_min = np.min(img)
            range_pair = "(" + str(im_min) +", "+ str(im_max) + ")"
            im_range_638.append(range_pair)
            print(im_max)
            print(im_min)
        # range_638 = process_638_list
        # for i in range(len(process_638_list)):
        #     for k in range(len(process_638_list[i])):
        #         range_638[i][k] = im_range_638[i]
    else:
        im_range_638 = []
    return im_range_488, im_range_561, im_range_638
            
            
#%%
def get_all_sb_image_list(SB_488_list, SB_561_list, SB_638_list):
    img_list = []
    total_488_list = []
    total_561_list = []
    total_638_list = []
    for i in range(len(SB_488_list)):
        total_488_list += SB_488_list[i]
    for j in range(len(SB_561_list)):
        total_561_list += SB_561_list[j]
    for k in range(len(SB_638_list)):
        total_638_list += SB_638_list[k]
    img_list = total_488_list + total_561_list + total_638_list
    return img_list

def get_save_path(img_list):
    name_list = get_last_name(img_list)
    path_name_list = []
    new_save_path_list = []
    for path_k in img_list:
        name_k = path_k.split("\\")[-1]
        length =len(name_k) +1
        name_folder = path_k[:-length]
        new_path_name = name_folder + '_Processed'
        path_name_list.append(new_path_name)
    for i in range(len(name_list)):
        save_path = path_name_list[i] + "\\" +  name_list[i]
        new_save_path_list.append(save_path)
    return new_save_path_list

def get_save_path_image(imname):
    name_k = imname.split("\\")[-1]
    length =len(name_k) +1
    name_folder = imname[:-length]
    new_path_name = name_folder + '_Processed'
    save_path = new_path_name + "\\" +  name_k
    return save_path
    
#%% Multiprocessing Deskew

def process_img(imname, param):
    angle = param[0]
    spacing_factor = param[1]
    blur = param[2]
    blur_sigma = param[3]
    
    #imname = path+'\\'+str(img_index)+".tif"
    name_k = imname.split("\\")[-1]
    length =len(name_k) +1
    name_folder = imname[:-length]
    new_path_name = name_folder + '_Processed'
    save_path = new_path_name + "\\" +  name_k
    
    if os.path.exists(imname):
        img0 = io.imread(imname)
        shear_factor = 1/np.tan(np.radians(angle))
        #img = ndimage.zoom(img0, (spacing_factor*np.sin(np.radians(angle)), np.sin(np.radians(angle)), 1))
        img = ndimage.zoom(img0, (spacing_factor, np.sin(np.radians(angle)), 1))
        
        (frame, x, y) = img.shape
        transform_matrix = [[1, shear_factor, 0],
                        [0, 1, 0],
                        [0, 0, 1]]
        sheared_img = ndimage.affine_transform(img, transform_matrix, offset = (-x*shear_factor, 0, 0), output_shape = (int(frame+x*shear_factor), x, y))
        #rotate_img = ndimage.zoom(sheared_img, (np.sin(np.radians(angle)), np.sin(np.radians(angle)), 1))
    
        if blur == True:
            blur_img = gaussian_filter1d(sheared_img, sigma = blur_sigma, axis = 1)
            blur_imname = save_path[:-4]+'_rotate_blur_angle_'+str(angle)+'_sf_'+str(spacing_factor)+'.tif'
            tifffile.imsave(blur_imname, blur_img)
        else:
            rotate_imname = save_path[:-4]+'_rotate_angle_'+str(angle)+'_sf_'+str(spacing_factor)+'.tif'
            tifffile.imsave(rotate_imname, sheared_img)


#%%  3D Projection
def projection(cropped, IJ, path_fiji, channel, process_488_list, process_561_list, process_638_list, folder_488_projection, folder_561_projection, folder_638_projection, direction, angle_increment, im_range_488, im_range_561, im_range_638):
    if cropped == True:
        IJ = imagej.init(path_fiji, mode = 'interactive')
        IJ.ui().showUI()
        
    if '488' in channel:
        for i in range(len(process_488_list)):
            for j in range(len(process_488_list[i])):
                image = process_488_list[i][j]
                image = image.replace("\\", "/")
                #imp =  IJ.io().open(image)         
                process_488 = folder_488_projection[i] + "/" + direction + "-" + str(j) + ".tif"
                process_488 = process_488.replace("\\", "/")
                macro_new2 = '''
                open("''' +image + '''");
                selectWindow("''' +image.split('/')[-1] +'''");
                setMinAndMax'''+im_range_488[i]+''';
                run("3D Project...", "projection=[Brightest Point] axis=''' + direction + '''-Axis slice=1 initial=0 total=360 rotation=''' + str(angle_increment) +''' lower=1 upper=255 opacity=0 surface=100 interior=50 interpolate");
                saveAs("Tiff", "'''+ process_488 + '''");
                close();
                close();
                '''
                IJ.py.run_macro(macro_new2)
                
    if '561' in channel:
        for i in range(len(process_561_list)):
            for j in range(len(process_561_list[i])):
                image = process_561_list[i][j]
                image = image.replace("\\", "/")
                #imp =  IJ.io().open(image)         
                process_561 = folder_561_projection[i] + "/" + direction + "-" + str(j) + ".tif"
                process_561 = process_561.replace("\\", "/")
                macro_new2 = '''
                open("''' +image + '''");
                selectWindow("''' +image.split('/')[-1] +'''");
                setMinAndMax'''+im_range_561[i]+''';
                run("3D Project...", "projection=[Brightest Point] axis=''' + direction + '''-Axis slice=1 initial=0 total=360 rotation=''' + str(angle_increment) +''' lower=1 upper=255 opacity=0 surface=100 interior=50 interpolate");
                saveAs("Tiff", "'''+ process_561 + '''");
                close();
                close();
                '''
                IJ.py.run_macro(macro_new2)
                
    if '638' in channel:
        for i in range(len(process_638_list)):
            for j in range(len(process_638_list[i])):
                image = process_638_list[i][j]
                image = image.replace("\\", "/")
                #imp =  IJ.io().open(image)         
                process_638 = folder_638_projection[i] + "/" + direction + "-" + str(j) + ".tif"
                process_638 = process_638.replace("\\", "/")
                macro_new2 = '''
                open("''' +image + '''");
                selectWindow("''' +image.split('/')[-1] +'''");
                setMinAndMax'''+im_range_638[i]+''';
                run("3D Project...", "projection=[Brightest Point] axis=''' + direction + '''-Axis slice=1 initial=0 total=360 rotation=''' + str(angle_increment) +''' lower=1 upper=255 opacity=0 surface=100 interior=50 interpolate");
                saveAs("Tiff", "'''+ process_638 + '''");
                close();
                close();
                '''
                IJ.py.run_macro(macro_new2)

    return 0

#%% Combining frames 


def combine_frame(ana_pos_list, projection_488_list, projection_561_list, projection_638_list, framenum, channel):
    new_dir_list = create_sub_folder(ana_pos_list, "Combine_frames")
    if '488' in channel:
        for i in range(len(projection_488_list)):
            for j in range(len(framenum)):     
                num_files = len(projection_488_list[i])
                sample_img = io.imread(projection_488_list[i][0])
                im_shape = (num_files, sample_img.shape[1], sample_img.shape[2])
                stack = np.zeros(im_shape, dtype = np.uint16)
                for l in range(0, num_files):
                    frame = io.imread(projection_488_list[i][l])[framenum[j]-1]
                    stack[l, :, :] = frame
                tifffile.imsave(new_dir_list[i] + r'/488_'+str(framenum[j]) + r'_Combined.tif', stack)
    if '561' in channel:
        for i in range(len(projection_561_list)):
            for j in range(len(framenum)):     
                num_files = len(projection_561_list[i])
                sample_img = io.imread(projection_561_list[i][0])
                im_shape = (num_files, sample_img.shape[1], sample_img.shape[2])
                stack = np.zeros(im_shape, dtype = np.uint16)
                for l in range(0, num_files):
                    frame = io.imread(projection_561_list[i][l])[framenum[j]-1]
                    stack[l, :, :] = frame
                tifffile.imsave(new_dir_list[i] + r'/561_'+str(framenum[j]) + r'_Combined.tif', stack)
    if '638' in channel:
        for i in range(len(projection_638_list)):
            for j in range(len(framenum)):     
                num_files = len(projection_638_list[i])
                sample_img = io.imread(projection_638_list[i][0])
                im_shape = (num_files, sample_img.shape[1], sample_img.shape[2])
                stack = np.zeros(im_shape, dtype = np.uint16)
                for l in range(0, num_files):
                    frame = io.imread(projection_638_list[i][l])[framenum[j]-1]
                    stack[l, :, :] = frame
                tifffile.imsave(new_dir_list[i] + r'/638_'+str(framenum[j]) + r'_Combined.tif', stack)               
                
    return 0


#%%


if __name__ == '__main__':
    print("Rename all files..")
    if renamed == False:
        rename(path)
        pos_name_list = get_img_list(path)
    else:
        pos_name_list = get_img_list(path)
    print("Cropping in Fiji..")
    ana_pos_list, folder_488, folder_561, folder_638, folder_488_process, folder_488_projection, folder_561_process, folder_561_projection, folder_638_process, folder_638_projection = create_analysis_folder(path, channel)
    if cropped == False:
        IJ = crop_images(path_fiji, pos_name_list, channel, angle_488, angle_561, angle_638, ROI_488, ROI_561, ROI_638, sub_background, folder_488, folder_561, folder_638)
    else:
        IJ = fake_ij(cropped)
    print("Deskewing..")
    SB_488_list, SB_561_list, SB_638_list = read_cropped_images(ana_pos_list, channel)
    img_list = get_all_sb_image_list(SB_488_list, SB_561_list, SB_638_list)
    
    if processed == False:
        partial_func = partial(process_img, param = param)
        if __name__ == '__main__':
            pool = Pool(num_workers)
            pool.map(partial_func, img_list)
            pool.close()
            pool.join()
    else:
        pass
    
    print("3D projection")   
    
    process_488_list, process_561_list, process_638_list, process_488_last_list, process_561_last_list, process_638_last_list = read_processed_images(ana_pos_list, channel)
    im_range_488, im_range_561, im_range_638 = get_intensity_range(channel, process_488_list, process_561_list, process_638_list, process_488_last_list, process_561_last_list, process_638_last_list)
    
    if projectioned == False:
        projection(cropped, IJ, path_fiji, channel, process_488_list, process_561_list, process_638_list, folder_488_projection, folder_561_projection, folder_638_projection, direction, angle_increment, im_range_488, im_range_561, im_range_638)
        pass
    print("Combining frame")
    projection_488_list, projection_561_list, projection_638_list = read_projection_images(ana_pos_list, channel)
    if combined == False:
        combine_frame(ana_pos_list, projection_488_list, projection_561_list, projection_638_list, framenum, channel)
    else: 
        pass
    print("Completed")
    
    
        
        
    
    
