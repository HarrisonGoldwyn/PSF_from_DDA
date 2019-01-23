# Part 1.
# top input "bricks.f90"
# top output "a.out", "shape_temp.dat"
rm a.out
rm shape_temp.dat
# middle input "ddscat_temp.par", "shape.dat" (obtained from "shape_temp.dat")
# middle ouput "ddscat.par"
rm ddscat.par
# bottom input "ddscat.par", "shape.dat", "var_temp.par"
# bottom output "var.par"
rm var.par

# Part 2.
# top input "ddscat.par", "shape.dat"
# top output "w000r000k000.pol1" (obtained by setting WPOL=1), "w000r000.avg", "qtable2", "qtable", "mtable", "ddfield.in" (version-specific), "Pol.dat" (specific ddscat version-specific)
rm w000r000k000.pol1
rm w000r000.avg
rm qtable2
rm qtable
rm mtable
rm ddfield.in
rm Pol.dat

# middle input "ddfield.in", "Pot.dat"
# middle output "tdda_input" 
rm tdda_input

# bottom input "myGreen.num_300" (obtained before the simulation), "var.par", "tdda_input"
# bottom output "temp.out"
rm temp.out

# Part 3.
# upper input "temp.out"
# upper output "shape.dat_new", "filler_ddscat.par"
rm shape.dat_new
rm filler_ddscat.par

# lower input "shape.dat_new", "filler_ddscat.par", "ddscat.par"
# lower output "ddscat.par_new", "shape.dat" (obtained by renaming), "ddscat.par" (obtained by renaming)
rm ddscat.par_new
rm shape.dat
# The following is DUPLICATE
#rm ddscat.par

# Part 4.
# top input None
# top output "ddfield.in"

# middle input "shape.dat", "ddscat.par" 
# middle output "w000r000k000.pol1" (obtained by setting WPOL=1), "w000r000.avg", "qtable2", "qtable", "mtable"
# The followings are DUPLICATE
#rm w000r000k000.pol1
#rm w000r000.avg
#rm qtable2
#rm qtable
#rm mtable

# bottom input "ddfield.in"
# bottom output "ddfield.E", "ddfield.B"
rm ddfield.E
rm ddfield.B

# Part 5.
rm parameters.yaml
# input "ddfield.E"
# output "pixel_coordinates_and_image_intensity.txt", "field_after_int_B.mat", "coordinates_B.mat", "1dGauss_B"
rm pixel_coordinates_and_image_intensity.txt
rm field_after_int_B.mat
rm coordinates_B.mat
rm 1dGauss_B
rm plot.py

# extra file
find . -name "*.pyc" -exec rm -rf {} \;

# cluster only
rm core.*
rm slurm*.out

rm *~ 
rm \#*
