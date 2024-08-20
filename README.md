# 3D-Deskew-and-Reconstruction-for-single-objective-SPIM-data

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

3D image reconstruction for simultaneous multicolour SPIM/OPM live data acquired with a tilted light sheet.  

Software
------------
ImageJ

Spyder (any pyhton IDE): Packages can be installed using the reconstruction.yml  

Overview
------------

This repository provides a comprehensive analysis pipeline for single-objective LSFM/SPIM modalities, particularly those based on the oblique plane microscopy (OPM) configuration. In OPM, a single objective lens is used for both illumination and detection, with the sample illuminated by an oblique light-sheet, resulting in a tilted imaging plane. This setup typically requires a remote refocusing detection system composed of two additional objective lenses aligned with the tilted light sheet. To achieve 3D imaging, the light sheet is either scanned across the sample or the sample is moved through the light sheet. Consequently, the resulting image stacks contain tomographic data that is inherently 'tilted' or 'skewed.' Simply stacking these images results in a distorted volumetric image, which usually requires the use of 'deskewing' or 'reconstruction' algorithms.

Key features
------------

This code incorporates a custom deskewing algorithm and offers additional functionality for analysing live or fixed imaging data acquired from multiple emission channels. Key features include:

1. Image Stack Renaming and Organisation: The script automatically renames and reorganises image stacks based on time order. Users specify the folder containing the imaging data, and the code handles the rest.

2. Cropping and Background Subtraction: Using ImageJ, the script crops the images and performs rolling ball background subtraction. Users can easily process data from multiple emission channels in simultaneous multicolor imaging by simply defining the regions of interest (ROI) for each channel.

3. Deskewing: The script applies the deskewing algorithm to the cropped images, with multiprocessing in Python used to accelerate this step.

4. 3D Reconstruction: Utilising the '3D projection' feature in ImageJ, the script reconstructs volumetric images based on the user-defined rotation axis.

5. Time-Series Video Creation: For live imaging data, the script combines reconstructed XY and XZ/YZ projections into a time-series video, enhancing visualization of dynamic processes.

This tool provides a seamless workflow from raw image stacks to fully reconstructed images, making it an essential resource for researchers working with LSFM/SPIM modalities.

