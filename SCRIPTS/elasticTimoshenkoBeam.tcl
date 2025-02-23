# APPENDIX 1: TCL SCRIPT FOR CANTILEVER BEAM WITH END TORQUE (ELASTIC
# BEAM-COLUMN WARPING ELEMENT)

# SET UP
wipe;
# source Wsection.tcl; 
# source LibUnits.tcl;
model basic -ndm 3 -ndf 6; 
set dir Cantilever_endtorque;
file mkdir $dir;
# clear opensess model
# include definition for I-section
# include units
# 3 dimensions, 7 dof per node
# create data directory

# define GEOMETRY 
# nodal coordinates 
node 1 0 0 0;
node 2 12 0 0;
node 3 24 0 0;
node 4 36 0 0;
node 5 48 0 0;
node 6 60 0 0;
node 7 72 0 0;
node 8 84 0 0;
node 9 96 0 0;
node 10 108 0 0;
node 11 120 0 0;
node 12 132 0 0;
node 13 144 0 0;
node 14 156 0 0;
node 15 168 0 0;
node 16 180 0 0;
node 17 192 0 0;
node 18 204 0 0;
node 19 216 0 0;
node 20 228 0 0;
node 21 240 0 0;
# Single point constraints - - Boundary Conditions 
fix 1 1 1 1 1 1 1;

# define material and section 
set poisson 0.3;
set G 11200.0;
set J 5.861;
set GJ [expr $G*$J];
set Cw 9902.0;
set E 29000.0;
set A 27.3;
set Iz 2070.0;
set Iy 92.9;
#W21x93 section
set Avy [expr 0.58*(21.6-(2*0.93))]; #Assume web area
set Avz [expr 2*0.93*8.42]; #Assume flange areas...

# DEFINE GEOMETERIC TRANSFORMATION FOR ELEMENTS ---------------------------------------------------------
set ColTransfTag 1;
#geomTransf Linear $transfTag $vecxzX $vecxzY $vecxzZ <-jntOffset $dXi $dYi $dZi $dXj $dYj $dZj>
geomTransf Corotational $ColTransfTag 0 0 1;

# DEFINE ELEMENTS ---------------------------------------------------------------------
#element elasticBeamColumn $eleTag 	  $iNode $jNode $A $E $G $J $Iy $Iz $transfTag <-mass $massDens> <-cMass>
#element ElasticTimoshenkoBeam3d $tag $iNode $jNode <$E $G $A $J $Iz $Iy $Avy $Avz>or<$sectionTag> $transTag <-mass $m> <-cMass>
element elasticTimoshenkoBeam 	1 1 2    $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag    ;
element elasticTimoshenkoBeam 	2 2 3    $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag    ;
element elasticTimoshenkoBeam 	3 3 4    $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag    ;
element elasticTimoshenkoBeam 	4 4 5    $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag    ;
element elasticTimoshenkoBeam 	5 5 6    $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag    ;
element elasticTimoshenkoBeam 	6 6 7    $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag    ;
element elasticTimoshenkoBeam 	7 7 8    $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag    ;
element elasticTimoshenkoBeam 	8 8 9    $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag    ;
element elasticTimoshenkoBeam 	9 9 10   $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag   ;
element elasticTimoshenkoBeam 	10 10 11 $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag ;
element elasticTimoshenkoBeam 	11 11 12 $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag ;
element elasticTimoshenkoBeam 	12 12 13 $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag ;
element elasticTimoshenkoBeam 	13 13 14 $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag ;
element elasticTimoshenkoBeam 	14 14 15 $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag ;
element elasticTimoshenkoBeam 	15 15 16 $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag ;
element elasticTimoshenkoBeam 	16 16 17 $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag ;
element elasticTimoshenkoBeam 	17 17 18 $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag ;
element elasticTimoshenkoBeam 	18 18 19 $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag ;
element elasticTimoshenkoBeam 	19 19 20 $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag ;
element elasticTimoshenkoBeam 	20 20 21 $E $G $A $J $Iz $Iy $Avy $Avz $ColTransfTag ;

# RECORD AND SAVE OUTPUT (TO BE SET BEFORE ANALYZE COMMAND) -------------------------------------------------------------
#recorder Node -file $dir/DFree.out -time -node 21 -dof 1 2 3 4 5 6 disp; 
#recorder Node -file $dir/DFree.out -node 21 -dof 1 2 3 4 5 6 disp; 
recorder Node -file $dir/DFree.out -closeOnWrite -time -node 21 -dof 1 2 3 4 5 6 disp; 
recorder Node -file $dir/React.out -closeOnWrite -time -node 1 -dof 1 2 3 4 5 6 reaction; 
recorder Element -file $dir/eleForces.out -ele 1 forces

# Records displacement at node 21

# DEFINE LOAD PATTERN (End Torque)----------------------------------------------------------------
set patternTag 1;

pattern Plain $patternTag Linear {
	#load $nodeTag (ndf $LoadValues) #kips
	load 21 0.0 0.0 0.0 1.0 0.0 1.0;
}

# CREATE THE CONSTRAINT HANDLER ------------------------------------------------------
constraints Plain;

# CREATE THE DOF NUMBERER ------------------------------------------------------------
numberer Plain;

# CREATE THE SYSTEM OF EQUATIONS -----------------------------------------------------
system BandGeneral;

# CREATE THE CONVERGENCE TEST --------------------------------------------------------
test NormDispIncr 1.0e-8 10 1; # The norm of the displacement increment with a tolerance of 1e-5 and a max number of iterations of 1000. The "1" or "0" at the end shows/doesn't show all iterations.

# CHOOSE ALGORITHM -------------------------------------------------------------------
# Create the solution algorithm. Choose between Newton, ModifiedNewton and ModifiedNewton -initial 
algorithm NewtonLineSearch; 

# CREATE THE INTEGRATION SCHEME ------------------------------------------------------
set lambda 0.01; # Set the load factor increment. A value of 1 indicates no further divison of load levels into steps. A value of 0.1, for example, would mean subdivision of each load step into 10 further steps.
integrator DisplacementControl 21 4 $lambda; 

# CREATE THE ANALYSIS OBJECT ---------------------------------------------------------
analysis Static; 

# ANALYZE ----------------------------------------------------------------------------
set NSteps [expr int(2./$lambda)]; # Number of steps in which the load, previously defined in pattern, is applied and the structure is analyzed. int() converts floating number into integer.
analyze $NSteps;

puts "Done!"

wipe;