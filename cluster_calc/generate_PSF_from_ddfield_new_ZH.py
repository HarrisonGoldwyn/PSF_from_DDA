from __future__ import print_function
from __future__ import division

import os
import sys
import yaml
import numpy as np
import scipy.io as sio
import scipy.optimize as opt

from textwrap import dedent

# working directory
work_dir = os.getcwd()

## import diffrantion integral solver from Optics folder
optics_folder = os.path.join(work_dir, 'Optics')
sys.path.append(optics_folder)
## diffraction integral solver
import diffraction_int as diffi
## fibonacci algoithm for spherical integral simplification
import fibonacci as fib

## import field functions
field_module_folder = os.path.join(work_dir, 'field_formulas')             
sys.path.append(field_module_folder)
import far_fields as fi

## import physical constants
opened_pc_folder = open('physical_constants.yaml','r')
physical_constants = yaml.load(opened_pc_folder)
hbar = physical_constants['hbar']
c = physical_constants['c']
nm = physical_constants['nm']

## open ddfield specific parameter file
open_param_file = open('parameters.yaml')
param = yaml.load(open_param_file)

##units: nm/dip_space * cm/nm = cm/dip_space 
dipole_spacing = param['ddfield']['dipole_spacing'] * nm  

## Load fields and points from ddfield.E
print('loading ddfield.E from directory {}'.format(work_dir))
loaded_ddfields = np.loadtxt('ddfield.E', skiprows=23)
raw_E = loaded_ddfields[:,-6:]
Ex = raw_E[:,0] + raw_E[:,1]*1j
Ey = raw_E[:,2] + raw_E[:,3]*1j
Ez = raw_E[:,4] + raw_E[:,5]*1j
scattered_E_field_on_sph = np.array([Ex, Ey, Ez]).T

## Load points on sphere 
raw_points_dipole_spacing = loaded_ddfields[:,:3] 
    ## in units of dip_space
##units: dip_space * cm/sip_space = cm
raw_points_in_cm = raw_points_dipole_spacing * dipole_spacing

## calculate focal length used in ddfield points
each_r_sqrd = np.sum(raw_points_in_cm**2., axis=1)
avrg_r_sqrd = each_r_sqrd.sum()/raw_points_in_cm.shape[0]
avrg_r = avrg_r_sqrd**0.5 
    ## r should really be square rooted before summing, but thats probably fine 
    ## as long as they are all the same anyway. 

## rotate coordinates and fields so optical axis is z, assuming x is optical exis
z_oriented_coords = np.zeros_like(raw_points_in_cm)
z_oriented_coords[:,0] = -raw_points_in_cm[:,2] ## = -(-cart_points_on_sph_in_nm[:,0])
z_oriented_coords[:,1] = raw_points_in_cm[:,1] 
z_oriented_coords[:,2] = raw_points_in_cm[:,0] ## = cart_points_on_sph_in_nm[:,2]

z_oriented_fields = np.zeros_like(scattered_E_field_on_sph)
z_oriented_fields[:,0] = -scattered_E_field_on_sph[:,2]
z_oriented_fields[:,1] = scattered_E_field_on_sph[:,1]
z_oriented_fields[:,2] = scattered_E_field_on_sph[:,0]

sphere_points = fib.cart_to_sphere(
    z_oriented_coords[:,0],
    z_oriented_coords[:,1],
    z_oriented_coords[:,2]
    ).T  

thetas_and_phis = np.array([sphere_points[:,1], sphere_points[:,2]]).T

## simulated image 
numerical_aperture = param['optics']['numerical_aperture']
max_theta = np.arcsin(numerical_aperture) 

sensor_size = param['optics']['photo_detector_width'] * nm

pixels = param['optics']['pixels']   # image grid pixels

## focal length determined by points
obj_f = avrg_r

lens_points = raw_points_in_cm.shape[0]

magnification = 1

drive_hbar_omega = 1240/param['system']['drive_wavelength']



## build image sensor
eye = diffi.observation_points( #returns 
    x_min= -sensor_size/2, 
    x_max= sensor_size/2,
    y_min= -sensor_size/2, 
    y_max= sensor_size/2, 
    points= pixels
    )
eye1 = (eye[1]/(nm * magnification)).tolist()
eye2 = (eye[2]/(nm * magnification)).tolist()

## diffract Elliots fields
diffracted_E_field = diffi.perform_integral(
    scattered_E=z_oriented_fields, 
    scattered_sph_coords=thetas_and_phis, 
    obser_pts=eye[0], 
    z=0, 
    obj_f=obj_f, 
    tube_f=magnification*obj_f, 
    k=(drive_hbar_omega/hbar)/c,
    alpha_1_max=max_theta
    )



pixels_on_2 = int(pixels/2) #pixels/2 takes the middle row of pixels
eye2wpixel = eye[2][pixels_on_2,:]
eye2wpixel2n = eye[2][pixels_on_2,:, None]
sio.savemat('field_after_int_B',{'field_after_int_B':diffracted_E_field})
sio.savemat('coordinates_B',{'coordinates_B':eye[2][pixels_on_2,:, None]})
print(eye[2])
mag_E_squared = np.real(
    np.sum(
        np.multiply(diffracted_E_field,
            np.conj(diffracted_E_field)
            ),
        axis=-1
        )
    )
max_E_sq = np.max(mag_E_squared)
min_E_sq = np.min(mag_E_squared)

reshaped_power_flux = mag_E_squared.reshape(eye[1].shape)
rpf = reshaped_power_flux.tolist()
power_slice = reshaped_power_flux[pixels_on_2,:]
ps = power_slice
psn = power_slice[:,None] 



## save image and pixel coordinates 
coordinates = eye[0]
# print('coordinates.shape', coordinates.shape)
np.savetxt(
    'pixel_coordinates_and_image_intensity.txt', 
    np.hstack((coordinates, mag_E_squared[:,None])),
    header='[x, y] pixel coordinates (cm); diffracted and focused |E|^2'
    )

np.savetxt(
    '1dGauss_B',
    np.hstack((eye2wpixel2n, psn)),
    )
