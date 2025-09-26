wipe;
source UNITS.tcl;

model basic -ndm 3 -ndf 11; 

set dir run15;
file mkdir $dir;

node 1 0. 0. 0.;
node 2 210. 0. 0.;

#node 1 0. 0. 0.;
#node 2 10. 0. 0.;
#node 3 20. 0. 0.;
#node 4 30. 0. 0.;
#node 5 40. 0. 0.;
#node 6 50. 0. 0.;
#node 7 60. 0. 0.;
#node 8 70. 0. 0.;
#node 9 80. 0. 0.;
#node 10 90. 0. 0.;
#node 11 100. 0. 0.;
#node 12 110. 0. 0.;
#node 13 120. 0. 0.;
#node 14 130. 0. 0.;
#node 15 140. 0. 0.;
#node 16 150. 0. 0.;
#node 17 160. 0. 0.;
#node 18 170. 0. 0.;
#node 19 180. 0. 0.;
#node 20 190. 0. 0.;
#node 21 200. 0. 0.;
#node 22 210. 0. 0.;

fix 1 1 1 1 1 1 1 1 1 1 1 1;
fix 2 0 0 0 1 1 1 1 1 1 1 1;
#fix 22 0 0 0 1 1 1;

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
set Avy [expr ($d-(2*$tf)) * $tw];
set Avz [expr ($bf) * $tf * 2];

set dw [expr $d - 2 * $tf];
set yfb [expr -$d/2 + $tf/2];

set zfb1 [expr -$bf/8*3];
set zfb2 [expr -$bf/8];
set zfb3 [expr  $bf/8];
set zfb4 [expr  $bf/8*3];

set yft [expr $d/2 - $tf/2];

set zft1 [expr -$bf/8*3];
set zft2 [expr -$bf/8];
set zft3 [expr  $bf/8];
set zft4 [expr  $bf/8*3];

set Af [expr $Avz/8]

set yw1 [expr $dw/8];
set yw2 [expr $dw/8*3];
set yw3 [expr -$dw/8];
set yw4 [expr -$dw/8*3];

set zw 0.0;

set Aw [expr $Avy/4]

set Fy [expr 60.0*$ksi];
set bRatio 1.0;
set Hiso 0;
set Hkin 1000;
set matID 7;
set matIDPlate 9;

uniaxialMaterial Elastic $matID $E;
#uniaxialMaterial Steel01 $matID $Fy $E $bRatio
#uniaxialMaterial Hardening $matID $E $Fy $Hiso $Hkin
uniaxialMaterial Elastic $matIDPlate $E;

set SecTag 18;
set nfdw 16;	# number of fibers along web depth 
set nftw 4;	# number of fibers along web thickness
set nfbf 16;	# number of fibers along flange width (you want this many in a bi-directional loading)
set nftf 4;	# number of fibers along flange thickness

#set matIDTorsion 99;		# ID tag for torsional section behavior
#uniaxialMaterial Elastic $matIDTorsion $GJ;	# define elastic torsional stiffness

puts "Pre section"  
   
  section taperedFiberSec $SecTag -hRatio 1.0 -hVal 1.0 -G $G -v $poisson {
     #                     nfIJ  nfJK    yI  zI    yJ  zJ    yK  zK    yL  zL
	 #TaperedFiber yLoc zLoc area matTag plateMatTag tP plFlag

     TaperedFiber $yfb $zfb1 $Af $matID $matIDPlate $tf 1
     TaperedFiber $yfb $zfb2 $Af $matID $matIDPlate $tf 1
	 TaperedFiber $yfb $zfb3 $Af $matID $matIDPlate $tf 1
     TaperedFiber $yfb $zfb4 $Af $matID $matIDPlate $tf 1
	 
	 TaperedFiber $yw1 $zw $Aw $matID $matIDPlate $tw 2
     TaperedFiber $yw2 $zw $Aw $matID $matIDPlate $tw 2
	 TaperedFiber $yw3 $zw $Aw $matID $matIDPlate $tw 2
     TaperedFiber $yw4 $zw $Aw $matID $matIDPlate $tw 2
	 
	 TaperedFiber $yft $zft1 $Af $matID $matIDPlate $tf 3
     TaperedFiber $yft $zft2 $Af $matID $matIDPlate $tf 3
	 TaperedFiber $yft $zft3 $Af $matID $matIDPlate $tf 3
     TaperedFiber $yft $zft4 $Af $matID $matIDPlate $tf 3	
  }

puts "Post section"  

#puts "Pre section"  
#   
#  section Fiber $SecTag -torsion $matIDTorsion {
#     #                     nfIJ  nfJK    yI  zI    yJ  zJ    yK  zK    yL  zL
#     patch quad $matID 4 1 $z1 $y1   $z4 $y1   $z4 $y2   $z1 $y2
#     patch quad $matID 1 4 $z2 $y2   $z3 $y2   $z3 $y3   $z2 $y3
#     patch quad $matID 4 1 $z1 $y3   $z4 $y3   $z4 $y4   $z1 $y4
#  }
#
#puts "Post section"  

#set matIDTorsion 99;		# ID tag for torsional section behavior
#set SecTag3D 3;			# ID tag for combined behavior for 3D model
#uniaxialMaterial Elastic $matIDTorsion $GJ;	# define elastic torsional stiffness
#section Aggregator $SecTag3D $matIDTorsion T -section $SecTag;	# combine section properties


set numIP 2;

set ColTransfTag 1;
#geomTransf Linear $ColTransfTag 0 0 1;
#geomTransf Corotational22 $ColTransfTag 0 -1 0;
geomTransf Corotational22 $ColTransfTag 0 0 1;

element TapereddispBeamColumnS 	1  1  2  $numIP -sections $SecTag $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	2  2  3  $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	3  3  4  $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	4  4  5  $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	5  5  6  $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	6  6  7  $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	7  7  8  $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	8  8  9  $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	9  9  10 $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	10 10 11 $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	11 11 12 $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	12 12 13 $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	13 13 14 $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	14 14 15 $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	15 15 16 $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	16 16 17 $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	17 17 18 $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	18 18 19 $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	19 19 20 $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	20 20 21 $numIP $SecTag $ColTransfTag ;
#element nonlinearBeamColumn 	21 21 22 $numIP $SecTag $ColTransfTag ;

#recorder Node -file $dir/DFree.out -closeOnWrite -time -node 2 -dof 1 2 3 4 5 6 7 8 9 10 11 disp; 
recorder Node -file $dir/DFree.out -closeOnWrite -time -node 2 disp; 

#recorder Node -file $dir/DFree.out -closeOnWrite -time -node 2 -dof 1 2 3 4 5 6 disp; 
#recorder Node -file $dir/React.out -closeOnWrite -time -node 1 -dof 1 2 3 4 5 6 7 8 9 10 11 reaction; 
recorder Node -file $dir/React.out -closeOnWrite -time -node 1 reaction; 
recorder Element -file $dir/eleForces.out -ele 1 forces
#recorder Element -file $dir/ele21Forces.out -ele 21 forces

set P [expr -4066.0*$kN] 
pattern Plain 2 "Constant" { 
 load 2 $P 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
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
	load 2 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
}

set lambda 0.2;
integrator DisplacementControl 2 2 $lambda; 

set NSteps [expr int(24./$lambda)]; 
analyze $NSteps;

puts "Done!"
wipe;