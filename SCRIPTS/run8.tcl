wipe;
source UNITS.tcl;

model basic -ndm 3 -ndf 6; 

set dir run8;
file mkdir $dir;

#node 1 0. 0. 0.;
#node 2 210. 0. 0.;

node 1 0. 0. 0.;
node 2 10. 0. 0.;
node 3 20. 0. 0.;
node 4 30. 0. 0.;
node 5 40. 0. 0.;
node 6 50. 0. 0.;
node 7 60. 0. 0.;
node 8 70. 0. 0.;
node 9 80. 0. 0.;
node 10 90. 0. 0.;
node 11 100. 0. 0.;
node 12 110. 0. 0.;
node 13 120. 0. 0.;
node 14 130. 0. 0.;
node 15 140. 0. 0.;
node 16 150. 0. 0.;
node 17 160. 0. 0.;
node 18 170. 0. 0.;
node 19 180. 0. 0.;
node 20 190. 0. 0.;
node 21 200. 0. 0.;
node 22 210. 0. 0.;

fix 1 1 1 1 1 1 1;
#fix 2 0 0 0 1 1 1;
fix 22 0 0 0 1 1 1;

set poisson 0.3;
set G 11200.0;
set J 15.6;
set GJ [expr $G*$J];
set Cw 129000.0;
set E 29000.0;
set A 50.9;
set Iz 8230.0;
set Iy 598.0;

set d 30.4;
set bf 15.0;
set tf 1.07;
set tw 0.655;
set Avy [expr ($d-(2*$tf)) * $tw]
set Avz [expr ($d-(2*$tf)) * $tw]

set ColTransfTag 1;
geomTransf Corotational $ColTransfTag 0 0 1;

#element ElasticTimoshenkoBeam $eleTag $iNode $jNode $E $G $A $Jx $Iy $Iz $Avy $Avz $transfTag <-mass $massDens> <-cMass>

element ElasticTimoshenkoBeam 	1  1  2  $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	2  2  3  $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	3  3  4  $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	4  4  5  $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	5  5  6  $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	6  6  7  $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	7  7  8  $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	8  8  9  $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	9  9  10 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	10 10 11 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	11 11 12 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	12 12 13 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	13 13 14 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	14 14 15 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	15 15 16 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	16 16 17 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	17 17 18 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	18 18 19 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	19 19 20 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	20 20 21 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;
element ElasticTimoshenkoBeam 	21 21 22 $E $G $A $J $Iy $Iz $Avy $Avz $ColTransfTag;

#recorder Node -file $dir/DFree.out -closeOnWrite -time -node 22 -dof 1 2 3 4 5 6 disp; 

recorder Node -file $dir/DFree.out -closeOnWrite -time -node 22 -dof 1 2 3 4 5 6 disp; 
recorder Node -file $dir/React.out -closeOnWrite -time -node 1 -dof 1 2 3 4 5 6 reaction; 
recorder Element -file $dir/eleForces.out -ele 1 forces

set P [expr 4066.0*$kN] 
pattern Plain 2 "Constant" { 
 load 22 $P 0.0 0.0 0.0 0.0 0.0
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
	load 22 0.0 1.0 0.0 0.0 0.0 0.0;
}

set lambda 0.2;
integrator DisplacementControl 22 2 $lambda; 

set NSteps [expr int(24./$lambda)]; 
analyze $NSteps;

puts "Done!"
wipe;