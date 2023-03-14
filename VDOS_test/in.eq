processors      * * * grid numa
units		metal
boundary	p p p

variable 	lat equal 5.43
lattice		diamond ${lat}
region		box prism 0 1.0 0 1.0 0 1.0 0.0 0.0 0.0
create_box	1 box
create_atoms	1 box

mass 		1 28.085

pair_style	sw
pair_coeff 	* * Si.sw Si

neighbor 	0.5 bin
neigh_modify 	once no every 1 delay 0 check yes

# Compute initial state
# Define minimization parameters
variable 	etol equal 1.0e-6
variable 	ftol equal 1.0e-10
variable 	maxiter equal 1000
variable 	maxeval equal 10000
variable 	dmax equal 1.0e-2
variable 	vmax equal 1.0e-2

# Setup minimization style
min_style	cg
min_modify	dmax ${dmax} line quadratic

fix      	brelax all box/relax iso 0.0 vmax ${vmax} fixedpoint 0 0 0
minimize 	${etol} ${ftol} ${maxiter} ${maxeval}
write_data 	relax_unitcell.dat

replicate 4 4 4

variable	T equal 300
velocity	all create  ${T} 73 dist gaussian
fix		NVT all nvt temp ${T} ${T} 1

# 500ps for equilibration
timestep 	0.001
run 		500000

write_data 	supercell_4.dat
