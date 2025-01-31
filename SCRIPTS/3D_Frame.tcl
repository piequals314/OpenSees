# SE 201B: NONLINEAR STRUCTURAL ANALYSIS (WI 2017)
# 3D LINEAR ELASTIC FRAME EXAMPLE - TUTORIAL 1

#Always start with
wipe; # Clear memory of all past model definitions

# UNITS: kip, in, sec (OpenSees doesn't have units. So be consistent!)
model BasicBuilder -ndm 3 -ndf 6; # Define the model builder, ndm=#dimension, ndf=#dofs

# SETUP DATA DIRECTORY FOR SAVING OUTPUTS --------------------------------------------
set dataDir "Results";	# Set up name of data directory
file mkdir $dataDir; # Create data directory

source UNITS.tcl;

# DEFINE NODES -----------------------------------------------------------------------
set Lx [expr 24.0*$ft];
set Ly [expr 15.0*$ft];
set Lz [expr 18.0*$ft];

#node $nodeTag (ndm $coords) <-mass (ndf $massValues)>
node 1	0.0	0.0	$Lz;
node 2	$Lx	0.0	$Lz;
node 3	$Lx	0.0	0.0;
node 4	0.0	0.0	0.0;
node 5	0.0	$Ly	$Lz;
node 6	$Lx	$Ly	$Lz;
node 7	$Lx	$Ly	0.0;
node 8	0.0	$Ly	0.0;

# APPLY CONSTRAINTS-------------------------------------------------------------------
#fix $nodeTag (ndf $constrValues)
fix 1 1 1 1 1 1 1;
fix 2 1 1 1 1 1 1;
fix 3 1 1 1 1 1 1;
fix 4 1 1 1 1 1 1;

# DEFINE MATERIAL & SECTION PARAMETERS ---------------------------------------------------------
set E		29000.0;#ksi
set mu		0.3;
set G		[expr $E/(2.0*(1.0+$mu))];

set Acol 	10.4; 	#in2
set Izcol 	50.5; 	#in4
set Iycol 	50.5; 	#in4
set Jcol 	85.6; 	#in4

set Abeam 	17.4; 	#in2
set Izbeam 	153.0; 	#in4
set Iybeam 	153.0; 	#in4
set Jbeam 	47.2; 	#in4



# DEFINE GEOMETERIC TRANSFORMATION FOR ELEMENTS ---------------------------------------------------------
set transfTagCol 1;
set transfTagBeam 2;
#geomTransf Linear $transfTag $vecxzX $vecxzY $vecxzZ <-jntOffset $dXi $dYi $dZi $dXj $dYj $dZj>
geomTransf Linear $transfTagCol 0 0 1;
geomTransf Linear $transfTagBeam 0 1 0;

# DEFINE ELEMENT ---------------------------------------------------------------------
#element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag <-mass $massDens> <-cMass>
element elasticBeamColumn 	1 1 5 $Acol $E $G $Jcol $Iycol $Izcol $transfTagCol;
element elasticBeamColumn 	2 2 6 $Acol $E $G $Jcol $Iycol $Izcol $transfTagCol;
element elasticBeamColumn 	3 3 7 $Acol $E $G $Jcol $Iycol $Izcol $transfTagCol;
element elasticBeamColumn 	4 4 8 $Acol $E $G $Jcol $Iycol $Izcol $transfTagCol;
element elasticBeamColumn 	5 5 6 $Abeam $E $G $Jbeam $Iybeam $Izbeam $transfTagBeam;
element elasticBeamColumn 	6 6 7 $Abeam $E $G $Jbeam $Iybeam $Izbeam $transfTagBeam;
element elasticBeamColumn 	7 7 8 $Abeam $E $G $Jbeam $Iybeam $Izbeam $transfTagBeam;
element elasticBeamColumn 	8 8 5 $Abeam $E $G $Jbeam $Iybeam $Izbeam $transfTagBeam;


# DEFINE LOAD PATTERN ----------------------------------------------------------------
set patternTag 1;

pattern Plain $patternTag "Linear" {
	#load $nodeTag (ndf $LoadValues) #kips
	load 5 10.0 -100.0 0.0 0.0 0.0 0.0;
	load 6 0.0 -100.0 -10.0 0.0 0.0 0.0;
	load 7 0.0 -100.0 0.0 0.0 0.0 0.0;
	load 8 0.0 -100.0 0.0 0.0 0.0 0.0;
}

# DEFINE ANALYSIS PARAMETERS ---------------------------------------------------------
# CREATE THE SYSTEM OF EQUATIONS -----------------------------------------------------
system BandGeneral;

# CREATE THE CONSTRAINT HANDLER ------------------------------------------------------
constraints Plain; 

# CREATE THE DOF NUMBERER ------------------------------------------------------------
numberer Plain; 

# CREATE THE CONVERGENCE TEST --------------------------------------------------------
test NormDispIncr 1.0e-5 1000 1; # The norm of the displacement increment with a tolerance of 1e-5 and a max number of iterations of 1000. The "1" or "0" at the end shows/doesn't show all iterations.

# CHOOSE ALGORITHM -------------------------------------------------------------------
# Create the solution algorithm. Choose between Newton, ModifiedNewton and ModifiedNewton -initial 
algorithm Newton; 

# CREATE THE INTEGRATION SCHEME ------------------------------------------------------
set lambda 1.; # Set the load factor increment. A value of 1 indicates no further divison of load levels into steps. A value of 0.1, for example, would mean subdivision of each load step into 10 further steps.
integrator LoadControl $lambda; # The LoadControl scheme

# CREATE THE ANALYSIS OBJECT ---------------------------------------------------------
analysis Static; 


# RECORD AND SAVE OUTPUT (TO BE SET BEFORE ANALYZE COMMAND) -------------------------------------------------------------
recorder Node -file $dataDir/disp.txt -node 6 -dof 1 disp; # Records displacement at node 6 along DOF 1
recorder Node -file $dataDir/reaction.txt -node 1 2 3 4 -dof 1 2 reaction; # Records reactions at node 1,2,3,4 along all DOFs (Specify DOF number after -dof to get reactions along a specific DOF)
recorder Element -file $dataDir/Ele1ForceLocal.txt -ele 1 localForce; # Records Element Forces in Element 1 (Element Local Coordinates)
recorder Element -file $dataDir/Ele1ForceGlobal.txt -ele 1 globalForce; # Records Element Forces in Element 1 (Global Coordinates)

# ANALYZE ----------------------------------------------------------------------------
set NSteps [expr int(1./$lambda)]; # Number of steps in which the load, previously defined in pattern, is applied and the structure is analyzed. int() converts floating number into integer.
analyze $NSteps;

# Don't forget to
wipe;
