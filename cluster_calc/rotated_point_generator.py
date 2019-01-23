from __future__ import print_function
from __future__ import division

import sys
import os
import numpy as np

work_dir = os.getcwd()
optics_folder = os.path.join(work_dir, 'Optics')
sys.path.append(optics_folder)
## diffraction integral solver
import diffraction_int as diffi
## fibonacci algoithm for spherical integral simplification
import fibonacci as fib

import yaml

## Import field functions
opened_pc_folder = open('physical_constants.yaml','r')
physical_constants = yaml.load(opened_pc_folder)
hbar = physical_constants['hbar']
c = physical_constants['c']
nm = physical_constants['nm']

## open ddfield specific parameter file
open_param_file = open('parameters.yaml')
param = yaml.load(open_param_file)
dipole_spacing = param['ddfield']['dipole_spacing'] * nm

## number of points on total sphere, which is then truncated by the max_angle
lens_points = param['lens_points']['number']  

## relationship between numerical aperture and max_angle for vacuum
numerical_aperture = param['optics']['numerical_aperture']
max_theta = np.arcsin(numerical_aperture) 

## objective lens focal length
obj_f =   param['optics']['focal_length']

sphere_points = fib.fib_alg_k_filter(
    num_points=lens_points, 
    max_ang=max_theta
    )

cart_points_on_sph = fib.sphere_to_cart(
    sphere_points[:,0],
    sphere_points[:,1],
    obj_f*np.ones(np.shape(sphere_points[:,0]))
    )

cart_points_on_sph_in_nm = cart_points_on_sph / dipole_spacing
# print(cart_points_on_sph_in_nm)
## rotate, x-> -z, z->x

rotated_point_array = np.copy(cart_points_on_sph_in_nm)

rotated_point_array[:,0] = cart_points_on_sph_in_nm[:,2]
rotated_point_array[:,2] = -cart_points_on_sph_in_nm[:,0]


doubled_w_ones = np.hstack(
    (
        rotated_point_array, 
        rotated_point_array, 
        np.ones((np.shape(rotated_point_array)[0],1))
    ))

np.savetxt('ddfield.in',doubled_w_ones, fmt='%10d', header='w000r000k000.pol1 =name of file with stored polarization\n5.00e-10 =gamma (interaction cutoff parameter)',comments='')
