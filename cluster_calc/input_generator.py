#! /usr/bin/env python

import sys
import math

def main():
    """\
    """
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-d', '--ddscat', action='store_const', const='ddscat',
                        dest='mode', help='Generate the ddscat.par file.')
    parser.add_argument('-v', '--var', action='store_const', const='var',
                        dest='mode', help='Generate the var.par file.')
    parser.add_argument('-dN', '--ddscatN', action='store_const', const='ddscatN',
                        dest='mode', help='Generate the ddscat.par_new file.')
    parser.add_argument('-y', '--yaml', action='store_const', const='yaml',
                        dest='mode', help='Generate the parameters.yaml file.')
    parser.add_argument('parIn', help='The parameter.input file.')
    parser.add_argument('sTemp', help='The shape_temp.dat file.')
    parser.add_argument('dTemp', help='The ddscat_temp.par file.')
    parser.add_argument('vTemp', help='The var_temp.par file.')
    parser.add_argument('sNew',  nargs = '?',  help='The shape.dat_new file.')
    parser.add_argument('dFill', nargs = '?', help='The filler_ddscat.par file.')
    args = parser.parse_args()
    
    # Obtain all the parameters
    par = []
    val = []
    with open(args.parIn) as file:
         for line in file:
             if "#" not in line and len(line.split()) is not 0:
                par.append(line.split()[0][:-1])
                if 'mem_allo' in line:
                   val.append(line.split()[1:4])
                else:
                   val.append(line.split()[1])
    
    # Template files to be used
    shape_temp = args.sTemp
    ddscat_temp = args.dTemp
    var_temp = args.vTemp
    shape_new = args.sNew
    filler = args.dFill

    # ddscatpar
    if args.mode == 'ddscat':
       ddscat_out = generate_ddscat(par, val, shape_temp, ddscat_temp) 
       # write to file
       with open('ddscat.par', 'w') as file:
            file.writelines(ddscat_out)

    # varpar
    if args.mode == 'var':
       var_out = generate_var(par, val, var_temp) 
       # write to file
       with open('var.par', 'w') as file:
            file.writelines(var_out)

    # ddscatpar new
    if args.mode == 'ddscatN':
       filler, ddscatN_out1, ddscatN_out2 = generate_ddscatN(par, val, shape_new, filler, ddscat_temp)
       # write to file
       with open('ddscat.par_new', 'w') as file:
            file.writelines(ddscatN_out1)
            file.writelines(filler)
#            file.writelines('\n')
            file.writelines(ddscatN_out2)

    # yaml par
    if args.mode == 'yaml':
       yaml_out = generate_yaml(par, val)
       # write to file
       with open('parameters.yaml', 'w') as file:
            file.writelines(yaml_out)

# Generate "ddscat.par" based on "shape_temp.dat" and "ddscat_temp.par"
def generate_ddscat(par, val, shape_temp, ddscat_temp):
    # dipole number
    file = open(shape_temp)
    for i, line in enumerate(file):
        if i == 1:
           dip_num = int(line.split()[0]) 
           break
    file.close()

    # effective radius
    dip_space = val[par.index("dip_space")]
    effR = (3*dip_num/(4*math.pi))**(1/3.0) * int(dip_space) * 10**(-3)
    effR = "{0:.4f}".format(effR)
    
    # make changes 
    with open(ddscat_temp, 'r') as file:
         data = file.readlines()
    for i in range(len(data)): 
        if "aeff" in data[i]:
           data[i] = data[i].split(' ')
           data[i][0] = data[i][1] = str(effR)
           data[i] = " ".join(data[i])
    
    return data

# Generate "var.par" based on "var_temp.par"
def generate_var(par, val, var_temp):
    # make changes 
    with open(var_temp, 'r') as file:
         data = file.readlines()
    for i in range(len(data)): 
        # num_k
        if "num_k" in data[i]:
           num_k = val[par.index("num_k")]
           data[i] = 'num_k: ' + str(num_k) + '\n'
        # k_out
        if "k_out" in data[i]:
           k_out = val[par.index("k_out")] 
           data[i] = 'k_out: ' + str(k_out) + '\n'
        # k_in
        if "k_in" in data[i]:
           k_in = val[par.index("k_in")]
           data[i] = 'k_in: ' + str(k_in) + '\n'
        # lambda
        if "lambda" in data[i]:
           lambda_ = val[par.index("lambda")]
           data[i] = 'lambda: ' + str(lambda_) + '\n'
        # n_m
        if "n_m" in data[i]:
           n_m = val[par.index("n_m")]
           data[i] = 'n_m: ' + str(n_m) + '\n'
        # I_0
        if "I_0" in data[i]:
           I_0 = val[par.index("I_0")]
           data[i] = 'I_0: ' + str(I_0) + '\n'
        # unit
        if "unit" in data[i]:
           unit = val[par.index("unit")]
           data[i] = 'unit: ' + str("{0:.1f}".format(float(unit))) + '\n'
        # d
        if "d:" in data[i]:
           d = val[par.index("d")]
           data[i] = 'd: ' + str(d) + '\n'
        # x min/max
        if "x_min" in data[i]:
           x_min = val[par.index("x_min")] 
           data[i] = 'x_min: ' + str(x_min) + '\n'
        if "x_max" in data[i]:
           x_max = val[par.index("x_max")] 
           data[i] = 'x_max: ' + str(x_max) + '\n'
        # y min/max
        if "y_min" in data[i]:
           y_min = val[par.index("y_min")] 
           data[i] = 'y_min: ' + str(y_min) + '\n'
        if "y_max" in data[i]:
           y_max = val[par.index("y_max")] 
           data[i] = 'y_max: ' + str(y_max) + '\n'
        # z min/max
        if "z_min" in data[i]:
           z_min = val[par.index("z_min")] 
           data[i] = 'z_min: ' + str(z_min) + '\n'
        if "z_max" in data[i]:
           z_max = val[par.index("z_max")] 
           data[i] = 'z_max: ' + str(z_max) + '\n'
        # x_plane_min/max
        if "x_plane_min" in data[i]:
           x_plane_min = val[par.index("x_plane_min")] 
           data[i] = 'x_plane_min: ' + str(x_plane_min) + '\n'
        if "x_plane_max" in data[i]:
           x_plane_max = val[par.index("x_plane_max")]
           data[i] = 'x_plane_max: ' + str(x_plane_max) + '\n'
    
    return data

# Generate "ddscat.par_new" based on "shape.dat_new" and "filler_ddscat.par"
def generate_ddscatN(par, val, shape_new, filler, ddscat_temp):
    # dipole number
    file = open(shape_new)
    for i, line in enumerate(file):
        if i == 1:
           dip_num = int(line.split()[0]) 
           break
    file.close()

    # effective radius
    dip_space = val[par.index("dip_space")]
    effR = (3*dip_num/(4*math.pi))**(1/3.0) * int(dip_space) * 10**(-3)
    effR = "{0:.4f}".format(effR)
    
    # memory allocation and incident frequency
    mem_allo = val[par.index("mem_allo")] 
    in_freq = val[par.index("in_freq")] 

    # the filler and its length
    with open(filler, 'r') as file:
         filler = file.readlines()
         filler_L = len(filler)

    # make changes 
    with open(ddscat_temp, 'r') as file:
         data = file.readlines()
    for i in range(len(data)):
        if "dimensioning" in data[i]:
           data[i] = data[i].split(' ')
           data[i][0:3] = mem_allo
           data[i] = " ".join(data[i])
        if "NCOMP" in data[i]:
           data[i] = data[i].split(' ')
           data[i][0] = str(filler_L)
           data[i] = " ".join(data[i])
           insertPos = i
        if "wavelengths" in data[i]:
           data[i] = data[i].split(' ')
           data[i][0] = data[i][1] = str(in_freq)
           data[i] = " ".join(data[i])
        if "aeff" in data[i]:
           data[i] = data[i].split(' ')
           data[i][0] = data[i][1] = str(effR)
           data[i] = " ".join(data[i])
        if "tab" in data[i]:
           indexToBeDel = i
    del data[indexToBeDel]

    data_1st = data[:insertPos+1]
    data_2nd = data[insertPos+1:]

    return filler, data_1st, data_2nd

# Generate "parameters.yaml"
def generate_yaml(par, val):
    number     = val[par.index("number")] 
    dip_space  = val[par.index("dip_space")]
    f_len      = val[par.index("f_len")] 
    num_ape    = val[par.index("num_ape")] 
    p_det_wid  = val[par.index("p_det_wid")] 
    pixels     = val[par.index("pixels")] 
    in_freq    = val[par.index("in_freq")] 
    in_freq_nm = int(float(in_freq) * 1000)

    str  = ("lens_points:\n")
    str += ("  number: %r\n\n" % int(number))
    str += ("ddfield:\n")
    str += ("  dipole_spacing: %r  ## in nm\n\n" % int(dip_space)) 
    str += ("optics:\n")
    str += ("  focal_length: %r  ## in cm, 0.1cm = 1 mm = 1000 um\n" % float(f_len))
    str += ("  numerical_aperture: %r ## unitless\n" % float(num_ape))
    str += ("  photo_detector_width: %r  ## in nm\n" % int(p_det_wid))
    str += ("  pixels: !!int %r  ## pixels per dimension of detector\n\n" % int(pixels))
    str += ("system:\n") 
    str += ("  drive_wavelength: %r  ## in nm, as long as 1240/drive_wavelength ~ eV\n\n" % in_freq_nm)
    str += ("plot:\n")
    str += ("  to_plot: True")

    return str

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
