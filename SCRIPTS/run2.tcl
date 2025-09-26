wipe;
source UNITS.tcl;

model basic -ndm 3 -ndf 6; 

set dir run2;
file mkdir $dir;

node 1 0. 0. 0.;
node 2 210 0. 0.;

fix 1 1 1 1 1 1 1;
fix 2 0 0 0 1 1 1;

set poisson 0.3;
set G 11200.0;
set J 15.6;
set GJ [expr $G*$J];
set Cw 129000.0;
set E 29000.0;
set A 50.9;
set Iz 8230.0;
set Iy 598.0;

set ColTransfTag 1;
geomTransf Corotational $ColTransfTag 0 0 1;

element elasticBeam 	1 1 2 $A $E $G $J $Iy $Iz $ColTransfTag;

recorder Node -file $dir/DFree.out -closeOnWrite -time -node 2 -dof 1 2 3 4 5 6 disp; 
recorder Node -file $dir/React.out -closeOnWrite -time -node 1 -dof 1 2 3 4 5 6 reaction; 
recorder Element -file $dir/eleForces.out -ele 1 forces

set P [expr -4066.0*$kN] 
pattern Plain 2 "Constant" { 
 load 2 $P 0.0 0.0 0.0 0.0 0.0
} 

integrator LoadControl 0.1

system BandGeneral
test NormDispIncr 1.0e-4 2000 4
numberer Plain 
constraints Plain 
algorithm Linear
analysis Static 
analyze 10 

set patternTag 3;
pattern Plain $patternTag Linear {
	load 2 0.0 1.0 0.0 0.0 0.0 0.0;
}

set lambda 0.2;
integrator DisplacementControl 2 2 $lambda; 

set NSteps [expr int(24./$lambda)]; 
analyze $NSteps;

puts "Done!"
wipe;