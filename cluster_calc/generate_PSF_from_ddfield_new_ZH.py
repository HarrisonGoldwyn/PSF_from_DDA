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
raw_points_in_cm = raw_points_dipole_spacing * dipole_spacing

## calculate focal length used in ddfield points
each_r_sqrd = np.sum(raw_points_in_cm**2., axis=1)
avrg_r_sqrd = each_r_sqrd.sum()/raw_points_in_cm.shape[0]
avrg_r = avrg_r_sqrd**0.5

## rotate coordinates and fields so optical axis is z, assuming x is optical exis
z_oriented_coords = np.zeros_like(raw_points_in_cm)
z_oriented_coords[:,0] = -raw_points_in_cm[:,2]
z_oriented_coords[:,1] = raw_points_in_cm[:,1]
z_oriented_coords[:,2] = raw_points_in_cm[:,0]

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

## test with dipole

di_E_field_on_sph = (
        fi.E_dipole_complex(
            w=drive_hbar_omega/hbar, 
            dipole_magnitude= 1e-14, 
            dipole_phase=1, 
            dipole_ori_unit=[1,0,0], 
            r=z_oriented_coords
            )
        )

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

## diffract test dipole fields
test_diffracted_E_field = diffi.perform_integral(
    scattered_E=di_E_field_on_sph, 
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

## test
test_mag_E_squared = np.real(
    np.sum(
        np.multiply(test_diffracted_E_field,
            np.conj(test_diffracted_E_field)
            ),
        axis=-1
        )
    )
max_E_sq_t = np.max(test_mag_E_squared)
min_E_sq_t = np.min(test_mag_E_squared)

test_reshaped_power_flux = test_mag_E_squared.reshape(eye[1].shape)
trpf = test_reshaped_power_flux.tolist()
test_power_slice = test_reshaped_power_flux[:,pixels_on_2]

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

def generate_script():
    with open('plot.py', "w") as plot_file:
        plot_file.write(write_plot())
    plot_file.close()

def write_plot():
    str = dedent('''\
    import numpy as np
    import matplotlib.pyplot as plt
    from numpy import array
    from mpl_toolkits import axes_grid1
    from mpl_toolkits.mplot3d import Axes3D
    
    #######################################################################
    ## colorbar stuff for plots
    def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
        """Add a vertical color bar to an image plot."""
        divider = axes_grid1.make_axes_locatable(im.axes)
        width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
        pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
        current_ax = plt.gca()
        cax = divider.append_axes("right", size=width, pad=pad)
        fig = plt.sca(current_ax)
        return im.axes.figure.colorbar(im, cax=cax, **kwargs)
    ########################################################################
    
    ##############################
    ## make ff plots
    ##############################
    fig , (test_ax, diff_ax) = plt.subplots(1,2, figsize=(12,6))
    diff_ax.set_xlabel('nm')
    diff_ax.set_aspect('equal', adjustable='box')
    diff_ax.set_title(r'Image $|E|^2$ from ddfield.E')
    ''')
    str += ("ff_max_color = %r\n" % max_E_sq)
    str += ("ff_min_color = %r\n" % min_E_sq)
    str += ("diff_color_plot = diff_ax.pcolor(")
    str += ("%r,\n" % eye1) 
    str += ("%r,\n" % eye2) 
    str += ("%r,\n" % rpf)
    str += ("vmin=ff_min_color,\n")
    str += ("vmax=ff_max_color,\n")
    str += ("cmap=plt.cm.jet,\n")
    str += (")\n")
         
    str += ("test_ax.set_aspect('equal', adjustable='box')\n")
    str += ("test_ax.set_title('PSF, test dipole (x oriented), analytic fields')\n")
    str += ("test_ax.set_ylabel('nm')\n")
    str += ("test_ax.set_xlabel('nm')\n")

    str += ("test_max_color = %r\n" % max_E_sq_t)
    str += ("test_min_color = %r\n" % min_E_sq_t)
    
    str += ("test_color_plot = test_ax.pcolor(")
    str += ("%r, \n" % eye1)
    str += ("%r, \n" % eye2)
    str += ("%r, \n" % trpf)
    str += ("vmin=test_min_color,\n")
    str += ("vmax=test_max_color,\n")
    str += ("cmap=plt.cm.jet,\n")
    str += (")\n")
         
    str += ("for pic in [diff_color_plot, test_color_plot]:\n")
    str += ("    add_colorbar(pic)\n")
    str += ("fig.savefig('heat_map.png', transparent=False, format='png')\n")
    str += ("fig1 = plt.figure()\n")
    str += ("plt.plot(%r, %r)\n" % (eye2wpixel, ps))
    str += ("fig1.savefig('gauss_b.png', transparent=False, format='png')\n")
    return str

if param['plot']['to_plot'] == True:
   # generate script
   generate_script()

