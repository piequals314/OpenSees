# APPENDIX 1: TCL SCRIPT FOR CANTILEVER BEAM WITH END TORQUE (ELASTIC
# BEAM-COLUMN WARPING ELEMENT)

# SET UP
wipe;


model basic -ndm 3 -ndf 6; 
set dir Cantilever_endtorque;
file mkdir $dir;

source LibUnits.tcl;
#source Ex9b.build.WSection3D.tcl; 

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


# DEFINE GEOMETERIC TRANSFORMATION FOR ELEMENTS ---------------------------------------------------------
set ColTransfTag 1;
#geomTransf Linear $transfTag $vecxzX $vecxzY $vecxzZ <-jntOffset $dXi $dYi $dZi $dXj $dYj $dZj>
geomTransf Corotational $ColTransfTag 0 0 1;

set numIP 2;

# DEFINE ELEMENTS ---------------------------------------------------------------------
#element elasticBeamColumn $eleTag 	  $iNode $jNode $A $E $G $J $Iy $Iz $transfTag <-mass $massDens> <-cMass>
#element ElasticTimoshenkoBeam3d $tag $iNode $jNode <$E $G $A $J $Iz $Iy $Avy $Avz>or<$sectionTag> $transTag <-mass $m> <-cMass>
#element nonlinearBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens> <-iter $maxIters $tol>



set Fy [expr 60.0*$ksi]
set Es [expr 29000*$ksi];		# Steel Young's Modulus
set nu 0.3;
set Gs [expr $Es/2./[expr 1+$nu]];  # Torsional stiffness Modulus
set Hiso 0
set Hkin 1000
set matIDhard 1
uniaxialMaterial Hardening  $matIDhard $Es $Fy   $Hiso  $Hkin

# Structural-Steel W-section properties -------------------------------------------------------------------
set SecTag 1
set WSec W27x114

# from Steel Manuals:
# in × lb/ft 	Area (in2) 	d (in) 	bf (in) 	tf (in) 	tw (in) 	Ixx (in4) 	Iyy (in4)
# W27x114  	33.5 		27.29 	10.07 	0.93 	0.57 	4090 	159
set d [expr 27.29*$in];	# nominal depth
set tw [expr 0.57*$in];	# web thickness
set bf [expr 0.57*$in];	# flange width
set tf [expr 10.07*$in];	# flange thickness
set nfdw 16;	# number of fibers along web depth 
set nftw 4;	# number of fibers along web thickness
set nfbf 16;	# number of fibers along flange width (you want this many in a bi-directional loading)
set nftf 4;	# number of fibers along flange thickness
  
  set dw [expr $d - 2 * $tf]
  set y1 [expr -$d/2]
  set y2 [expr -$dw/2]
  set y3 [expr  $dw/2]
  set y4 [expr  $d/2]
  
  set z1 [expr -$bf/2]
  set z2 [expr -$tw/2]
  set z3 [expr  $tw/2]
  set z4 [expr  $bf/2]
  
  #     
puts "Pre section"  
  section fiberSec  $SecTag -GJ 1.0e7 {
     #                     nfIJ  nfJK    yI  zI    yJ  zJ    yK  zK    yL  zL
     patch quadr  $matIDhard $nfbf $nftf   $y1 $z4   $y1 $z1   $y2 $z1   $y2 $z4
     patch quadr  $matIDhard $nftw $nfdw   $y2 $z3   $y2 $z2   $y3 $z2   $y3 $z3
     patch quadr  $matIDhard $nfbf $nftf   $y3 $z4   $y3 $z1   $y4 $z1   $y4 $z4
  }
puts "Post section"  
# assign torsional Stiffness for 3D Model
set SecTagTorsion 99;		# ID tag for torsional section behavior
set SecTag3D 3;			# ID tag for combined behavior for 3D model
uniaxialMaterial Elastic $SecTagTorsion $Ubig;	# define elastic torsional stiffness
section Aggregator $SecTag3D $SecTagTorsion T -section $SecTag;	# combine section properties










puts "Model Built Not"




element nonlinearBeamColumn 	1 1 2    $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	2 2 3    $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	3 3 4    $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	4 4 5    $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	5 5 6    $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	6 6 7    $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	7 7 8    $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	8 8 9    $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	9 9 10   $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	10 10 11 $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	11 11 12 $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	12 12 13 $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	13 13 14 $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	14 14 15 $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	15 15 16 $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	16 16 17 $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	17 17 18 $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	18 18 19 $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	19 19 20 $numIP $SecTag3D $ColTransfTag ;
element nonlinearBeamColumn 	20 20 21 $numIP $SecTag3D $ColTransfTag ;

puts "Element Created"

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
	load 21 0.0 0.0 0.0 0.0 0.0 1.0;
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
integrator DisplacementControl 21 6 $lambda; 

# CREATE THE ANALYSIS OBJECT ---------------------------------------------------------
analysis Static; 

# ANALYZE ----------------------------------------------------------------------------
set NSteps [expr int(2./$lambda)]; # Number of steps in which the load, previously defined in pattern, is applied and the structure is analyzed. int() converts floating number into integer.
analyze $NSteps;

puts "Done!"

wipe;