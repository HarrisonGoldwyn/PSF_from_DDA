#!/bin/bash
## Job Name
#SBATCH --job-name=centered
## Allocation Definition
#SBATCH --account=chem
#SBATCH --partition=chem
## Resources
## Nodes
#SBATCH --nodes=1
## Tasks per node (28 is Slurm default)
#SBATCH --ntasks-per-node=28
## Walltime (days-HH:MM:SS)
#SBATCH --time=12:00:00
## Memory per node
#SBATCH --mem=20G
## Specify the working directory for the job
##SBATCH --workdir=/usr/lusers/clairew4/test/thermal_pipeline_cluster
##Output file
#SBTACH --output cfs.out

module load anaconda3_4.3.1
module load matlab_2017a

# DDA_background for thermal
gfortran bricks_ZH.f90 && ./a.out
python input_generator.py -d parameters.input shape_temp.dat ddscat_temp.par var_temp.par 
python input_generator.py -v parameters.input shape_temp.dat ddscat_temp.par var_temp.par 
cp shape_temp.dat shape.dat
#/gscratch/chem/masiello_group/srcPW_tDDAb/ddscat
#/gscratch/chem/masiello_group/srcPW_tDDAb/ddfield

# tDDA
#/gscratch/chem/masiello_group/tDDA_0712/Lattice_Diffusion /gscratch/chem/masiello_group/myGreen.num_300 var.par tdda_input temp.out

# DDA_background for psf
#python check_temp.py temp.out
#matlab -r -nodisplay makemetalpsf
#python input_generator.py -dN parameters.input shape_temp.dat ddscat_temp.par var_temp.par shape.dat_new filler_ddscat.par
#cp ddscat.par_new ddscat.par && cp shape.dat_new shape.dat
python input_generator.py -y parameters.input shape_temp.dat ddscat_temp.par var_temp.par
python rotated_point_generator.py
/gscratch/chem/masiello_group/srcPW_background_psf/ddscat
/gscratch/chem/masiello_group/srcPW_background_psf/ddfield

# PSF generation
python generate_PSF_from_ddfield_new_ZH.py
