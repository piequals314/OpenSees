 /* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.15 $
// $Date: 2010-06-01 23:44:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/PDeltaCrdTransf3d22.cpp,v $


// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Modified: 04/2005 Andreas Schellenberg (getBasicTrialVel, getBasicTrialAccel)
//
// Purpose: This file contains the implementation for the
// PDeltaCrdTransf3d22 class. PDeltaCrdTransf3d22 is a linear
// transformation for a planar frame between the global
// and basic coordinate systems


#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>

#include <PDeltaCrdTransf3d22.h>


// constructor:
PDeltaCrdTransf3d22::PDeltaCrdTransf3d22(int tag, const Vector &vecInLocXZPlane):
CrdTransf(tag, CRDTR_TAG_PDeltaCrdTransf3d22),
nodeIPtr(0), nodeJPtr(0),
L(0), ul112(0), ul213(0),
nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 3; j++)
            R[i][j] = 0.0;
       
        R[2][0] = vecInLocXZPlane(0);
        R[2][1] = vecInLocXZPlane(1);
        R[2][2] = vecInLocXZPlane(2);
       
        // Does nothing
}


// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
PDeltaCrdTransf3d22::PDeltaCrdTransf3d22():
CrdTransf(0, CRDTR_TAG_PDeltaCrdTransf3d22),
nodeIPtr(0), nodeJPtr(0),
L(0), ul112(0), ul213(0),
nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R[i][j] = 0.0;
}


// destructor:
PDeltaCrdTransf3d22::~PDeltaCrdTransf3d22()
{
    if (nodeIInitialDisp != 0)
        delete [] nodeIInitialDisp;
    if (nodeJInitialDisp != 0)
        delete [] nodeJInitialDisp;
}


int
PDeltaCrdTransf3d22::commitState(void)
{
    return 0;
}


int
PDeltaCrdTransf3d22::revertToLastCommit(void)
{
    return 0;
}


int
PDeltaCrdTransf3d22::revertToStart(void)
{
    return 0;
}


int
PDeltaCrdTransf3d22::initialize(Node *nodeIPointer, Node *nodeJPointer)
{      
    int error;
   
    nodeIPtr = nodeIPointer;
    nodeJPtr = nodeJPointer;
   
    if ((!nodeIPtr) || (!nodeJPtr))
    {
        opserr << "\nPDeltaCrdTransf3d22::initialize";
        opserr << "\ninvalid pointers to the element nodes\n";
        return -1;
    }
   
    // see if there is some initial displacements at nodes
    if (initialDispChecked == false) {
        const Vector &nodeIDisp = nodeIPtr->getDisp();
        const Vector &nodeJDisp = nodeJPtr->getDisp();

        for (int i=0; i<11; i++)
            if (nodeIDisp(i) != 0.0) {
                nodeIInitialDisp = new double [11];
                for (int j=0; j<11; j++)
                     nodeIInitialDisp[j] = nodeIDisp(j);
                i = 11;
            }
           
	for (int j=0; j<11; j++)
	    if (nodeJDisp(j) != 0.0) {
		nodeJInitialDisp = new double [11];
		for (int i=0; i<11; i++)
		    nodeJInitialDisp[i] = nodeJDisp(i);
		j = 11;
	    }
	   
	initialDispChecked = true;
    }
   
    // get element length and orientation
    if ((error = this->computeElemtLengthAndOrient()))
        return error;
   
    static Vector XAxis(3);
    static Vector YAxis(3);
    static Vector ZAxis(3);
   
    // get 3by3 rotation matrix
    if ((error = this->getLocalAxes(XAxis, YAxis, ZAxis)))      
        return error;
   
    return 0;
}


int
PDeltaCrdTransf3d22::update(void)
{
    const Vector &disp1 = nodeIPtr->getTrialDisp();
    const Vector &disp2 = nodeJPtr->getTrialDisp();
   
    static double ug[22];
    for (int i = 0; i < 11; i++) {
        ug[i]    = disp1(i);
        ug[i+11] = disp2(i);
    }
   
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<11; j++)
            ug[j] -= nodeIInitialDisp[j];
    }
   
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<11; j++)
            ug[j+11] -= nodeJInitialDisp[j];
    }
   
    double ul1, ul12, ul2, ul13;
   
    ul1 = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
    ul2 = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];
   
    ul12 = R[1][0]*ug[11] + R[1][1]*ug[12] + R[1][2]*ug[13];
    ul13 = R[2][0]*ug[11] + R[2][1]*ug[12] + R[2][2]*ug[13];
   
    ul112 = ul1-ul12;
    ul213 = ul2-ul13;
   
    return 0;
}


int
PDeltaCrdTransf3d22::computeElemtLengthAndOrient()
{
    // element projection
    static Vector dx(3);
   
    const Vector &ndICoords = nodeIPtr->getCrds();
    const Vector &ndJCoords = nodeJPtr->getCrds();
   
    dx(0) = ndJCoords(0) - ndICoords(0);
    dx(1) = ndJCoords(1) - ndICoords(1);
    dx(2) = ndJCoords(2) - ndICoords(2);
   
    if (nodeIInitialDisp != 0) {
        dx(0) -= nodeIInitialDisp[0];
        dx(1) -= nodeIInitialDisp[1];
        dx(2) -= nodeIInitialDisp[2];
    }
   
    if (nodeJInitialDisp != 0) {
        dx(0) += nodeJInitialDisp[0];
        dx(1) += nodeJInitialDisp[1];
        dx(2) += nodeJInitialDisp[2];
    }
   
    // calculate the element length
    L = dx.Norm();
   
    if (L == 0.0) {
        opserr << "\nPDeltaCrdTransf3d22::computeElemtLengthAndOrien: 0 length\n";
        return -2;  
    }
   
    // calculate the element local x axis components (direction cossines)
    // wrt to the global coordinates
    R[0][0] = dx(0)/L;
    R[0][1] = dx(1)/L;
    R[0][2] = dx(2)/L;
   
    return 0;
}


int
PDeltaCrdTransf3d22::getLocalAxes(Vector &XAxis, Vector &YAxis, Vector &ZAxis)
{
    // Compute y = v cross x
    // Note: v(i) is stored in R[2][i]
    static Vector vAxis(3);
    vAxis(0) = R[2][0]; vAxis(1) = R[2][1];     vAxis(2) = R[2][2];
   
    static Vector xAxis(3);
    xAxis(0) = R[0][0]; xAxis(1) = R[0][1];     xAxis(2) = R[0][2];
    XAxis(0) = xAxis(0);    XAxis(1) = xAxis(1);    XAxis(2) = xAxis(2);
   
    static Vector yAxis(3);
   
    yAxis(0) = vAxis(1)*xAxis(2) - vAxis(2)*xAxis(1);
    yAxis(1) = vAxis(2)*xAxis(0) - vAxis(0)*xAxis(2);
    yAxis(2) = vAxis(0)*xAxis(1) - vAxis(1)*xAxis(0);
   
    double ynorm = yAxis.Norm();
   
    if (ynorm == 0) {
        opserr << "\nPDeltaCrdTransf3d22::getLocalAxes";
        opserr << "\nvector v that defines plane xz is parallel to x axis\n";
        return -3;
    }
   
    yAxis /= ynorm;
   
    YAxis(0) = yAxis(0);    YAxis(1) = yAxis(1);    YAxis(2) = yAxis(2);
   
    // Compute z = x cross y
    static Vector zAxis(3);
   
    zAxis(0) = xAxis(1)*yAxis(2) - xAxis(2)*yAxis(1);
    zAxis(1) = xAxis(2)*yAxis(0) - xAxis(0)*yAxis(2);
    zAxis(2) = xAxis(0)*yAxis(1) - xAxis(1)*yAxis(0);
    ZAxis(0) = zAxis(0);    ZAxis(1) = zAxis(1);    ZAxis(2) = zAxis(2);
   
    // Fill in transformation matrix
    R[1][0] = yAxis(0);
    R[1][1] = yAxis(1);
    R[1][2] = yAxis(2);
   
    R[2][0] = zAxis(0);
    R[2][1] = zAxis(1);
    R[2][2] = zAxis(2);
   
    return 0;
}


double
PDeltaCrdTransf3d22::getInitialLength(void)
{
    return L;
}


double
PDeltaCrdTransf3d22::getDeformedLength(void)
{
    return L;
}


const Vector &
PDeltaCrdTransf3d22::getBasicTrialDisp (void)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getTrialDisp();
    const Vector &disp2 = nodeJPtr->getTrialDisp();
   
    static double ug[22];
    for (int i = 0; i < 11; i++) {
        ug[i]    = disp1(i);
        ug[i+11] = disp2(i);
    }
   
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<11; j++)
            ug[j] -= nodeIInitialDisp[j];
    }
   
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<11; j++)
            ug[j+11] -= nodeJInitialDisp[j];
    }
   
    double oneOverL = 1.0/L;
   
    static Vector ub(17);
   
    static double ul[22];
    ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
    ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
    ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];
      
    ul[3]  = R[0][0]*ug[3] + R[0][1]*ug[4] + R[0][2]*ug[5];
    ul[4]  = R[1][0]*ug[3] + R[1][1]*ug[4] + R[1][2]*ug[5];
    ul[5]  = R[2][0]*ug[3] + R[2][1]*ug[4] + R[2][2]*ug[5];
      
    ul[6]   = ug[6];  // do not transform warping
    ul[7]   = ug[7];  // do not transform top flange rotation
    ul[8]   = ug[8];  // do not transform top flange curvature
    ul[9]   = ug[9];  // do not transform bottom flange rotation
    ul[10]  = ug[10]; // do not transform bottom flange curvature
  
    ul[11]   = R[0][0]*ug[11] + R[0][1]*ug[12] + R[0][2]*ug[13];
    ul[12]   = R[1][0]*ug[11] + R[1][1]*ug[12] + R[1][2]*ug[13];
    ul[13]   = R[2][0]*ug[11] + R[2][1]*ug[12] + R[2][2]*ug[13];
      
    ul[14]  = R[0][0]*ug[14] + R[0][1]*ug[15] + R[0][2]*ug[16];
    ul[15]  = R[1][0]*ug[14] + R[1][1]*ug[15] + R[1][2]*ug[16];
    ul[16]  = R[2][0]*ug[14] + R[2][1]*ug[15] + R[2][2]*ug[16];
      
    ul[17]  = ug[17]; // do not transform warping
    ul[18]  = ug[18]; // do not transform top flange rotation
    ul[19]  = ug[19]; // do not transform top flange curvature
    ul[20]  = ug[20]; // do not transform bottom flange rotation
    ul[21]  = ug[21]; // do not transform bottom flange curvature 
   
    
   
    ub(16) = ul[11] - ul[0];
    double tmp;
    // Theta_z
    tmp = oneOverL*(ul[1]-ul[12]);
    ub(1) = ul[5] + tmp;
    ub(9) = ul[16] + tmp;
    // Theta_y
    tmp = oneOverL*(ul[13]-ul[2]);
    ub(2)  = ul[4] + tmp;
    ub(10) = ul[15] + tmp;
    // Theta_x = Torsion
    ub(0)  = (-ul[14] + ul[3])/2;
    ub(8)  = -ub(0);
    // Bi-moment
    ub(3)  = ul[6];
    ub(11) = ul[17];
    // Top Flange Rotations & Curvatures
    ub(4)  = ul[7];
    ub(5)  = ul[8];
    ub(12) = ul[18];
    ub(13) = ul[19];
    // Bottom Flange Rotations & Curvatures
    ub(6)  = ul[9];
    ub(7)  = ul[10];
    ub(14) = ul[20];
    ub(15) = ul[21]; 
   
    return ub;
}


const Vector &
PDeltaCrdTransf3d22::getBasicIncrDisp (void)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getIncrDisp();
    const Vector &disp2 = nodeJPtr->getIncrDisp();
   
    static double ug[22];
    for (int i = 0; i < 11; i++) {
        ug[i]    = disp1(i);
        ug[i+11] = disp2(i);
    }
   
    double oneOverL = 1.0/L;
   
    static Vector ub(17);
   
    static double ul[22];
    ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
    ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
    ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];
      
    ul[3]  = R[0][0]*ug[3] + R[0][1]*ug[4] + R[0][2]*ug[5];
    ul[4]  = R[1][0]*ug[3] + R[1][1]*ug[4] + R[1][2]*ug[5];
    ul[5]  = R[2][0]*ug[3] + R[2][1]*ug[4] + R[2][2]*ug[5];
      
    ul[6]   = ug[6];  // do not transform warping
    ul[7]   = ug[7];  // do not transform top flange rotation
    ul[8]   = ug[8];  // do not transform top flange curvature
    ul[9]   = ug[9];  // do not transform bottom flange rotation
    ul[10]  = ug[10]; // do not transform bottom flange curvature
  
    ul[11]   = R[0][0]*ug[11] + R[0][1]*ug[12] + R[0][2]*ug[13];
    ul[12]   = R[1][0]*ug[11] + R[1][1]*ug[12] + R[1][2]*ug[13];
    ul[13]   = R[2][0]*ug[11] + R[2][1]*ug[12] + R[2][2]*ug[13];
      
    ul[14]  = R[0][0]*ug[14] + R[0][1]*ug[15] + R[0][2]*ug[16];
    ul[15]  = R[1][0]*ug[14] + R[1][1]*ug[15] + R[1][2]*ug[16];
    ul[16]  = R[2][0]*ug[14] + R[2][1]*ug[15] + R[2][2]*ug[16];
      
    ul[17]  = ug[17]; // do not transform warping
    ul[18]  = ug[18]; // do not transform top flange rotation
    ul[19]  = ug[19]; // do not transform top flange curvature
    ul[20]  = ug[20]; // do not transform bottom flange rotation
    ul[21]  = ug[21]; // do not transform bottom flange curvature 
   

    ub(16) = ul[11] - ul[0];
    double tmp;
    // Theta_z
    tmp = oneOverL*(ul[1]-ul[12]);
    ub(1) = ul[5] + tmp;
    ub(9) = ul[16] + tmp;
    // Theta_y
    tmp = oneOverL*(ul[13]-ul[2]);
    ub(2)  = ul[4] + tmp;
    ub(10) = ul[15] + tmp;
    // Theta_x = Torsion
    ub(0)  = (-ul[14] + ul[3])/2;
    ub(8)  = -ub(0);
    // Bi-moment
    ub(3)  = ul[6];
    ub(11) = ul[17];
    // Top Flange Rotations & Curvatures
    ub(4)  = ul[7];
    ub(5)  = ul[8];
    ub(12) = ul[18];
    ub(13) = ul[19];
    // Bottom Flange Rotations & Curvatures
    ub(6)  = ul[9];
    ub(7)  = ul[10];
    ub(14) = ul[20];
    ub(15) = ul[21]; 
   
    return ub;
}


const Vector &
PDeltaCrdTransf3d22::getBasicIncrDeltaDisp(void)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getIncrDeltaDisp();
    const Vector &disp2 = nodeJPtr->getIncrDeltaDisp();
   
    static double ug[22];
    for (int i = 0; i < 11; i++) {
        ug[i]    = disp1(i);
        ug[i+11] = disp2(i);
    }
   
    double oneOverL = 1.0/L;
   
    static Vector ub(17);
   
    static double ul[22];
    ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
    ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
    ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];
      
    ul[3]  = R[0][0]*ug[3] + R[0][1]*ug[4] + R[0][2]*ug[5];
    ul[4]  = R[1][0]*ug[3] + R[1][1]*ug[4] + R[1][2]*ug[5];
    ul[5]  = R[2][0]*ug[3] + R[2][1]*ug[4] + R[2][2]*ug[5];
      
    ul[6]   = ug[6];  // do not transform warping
    ul[7]   = ug[7];  // do not transform top flange rotation
    ul[8]   = ug[8];  // do not transform top flange curvature
    ul[9]   = ug[9];  // do not transform bottom flange rotation
    ul[10]  = ug[10]; // do not transform bottom flange curvature
  
    ul[11]   = R[0][0]*ug[11] + R[0][1]*ug[12] + R[0][2]*ug[13];
    ul[12]   = R[1][0]*ug[11] + R[1][1]*ug[12] + R[1][2]*ug[13];
    ul[13]   = R[2][0]*ug[11] + R[2][1]*ug[12] + R[2][2]*ug[13];
      
    ul[14]  = R[0][0]*ug[14] + R[0][1]*ug[15] + R[0][2]*ug[16];
    ul[15]  = R[1][0]*ug[14] + R[1][1]*ug[15] + R[1][2]*ug[16];
    ul[16]  = R[2][0]*ug[14] + R[2][1]*ug[15] + R[2][2]*ug[16];
      
    ul[17]  = ug[17]; // do not transform warping
    ul[18]  = ug[18]; // do not transform top flange rotation
    ul[19]  = ug[19]; // do not transform top flange curvature
    ul[20]  = ug[20]; // do not transform bottom flange rotation
    ul[21]  = ug[21]; // do not transform bottom flange curvature 
   

    ub(16) = ul[11] - ul[0];
    double tmp;
    // Theta_z
    tmp = oneOverL*(ul[1]-ul[12]);
    ub(1) = ul[5] + tmp;
    ub(9) = ul[16] + tmp;
    // Theta_y
    tmp = oneOverL*(ul[13]-ul[2]);
    ub(2)  = ul[4] + tmp;
    ub(10) = ul[15] + tmp;
    // Theta_x = Torsion
    ub(0)  = (-ul[14] + ul[3])/2;
    ub(8)  = -ub(0);
    // Bi-moment
    ub(3)  = ul[6];
    ub(11) = ul[17];
    // Top Flange Rotations & Curvatures
    ub(4)  = ul[7];
    ub(5)  = ul[8];
    ub(12) = ul[18];
    ub(13) = ul[19];
    // Bottom Flange Rotations & Curvatures
    ub(6)  = ul[9];
    ub(7)  = ul[10];
    ub(14) = ul[20];
    ub(15) = ul[21]; 
   
    return ub;
}


const Vector &
PDeltaCrdTransf3d22::getBasicTrialVel(void)
{
        // determine global velocities
        const Vector &vel1 = nodeIPtr->getTrialVel();
        const Vector &vel2 = nodeJPtr->getTrialVel();
       
        static double vg[22];
        for (int i = 0; i < 11; i++) {
                vg[i]    = vel1(i);
                vg[i+11] = vel2(i);
        }
       
        double oneOverL = 1.0/L;
       
        static Vector vb(17);
       
        static double vl[22];
       
    vl[0]  = R[0][0]*vg[0] + R[0][1]*vg[1] + R[0][2]*vg[2];
    vl[1]  = R[1][0]*vg[0] + R[1][1]*vg[1] + R[1][2]*vg[2];
    vl[2]  = R[2][0]*vg[0] + R[2][1]*vg[1] + R[2][2]*vg[2];
      
    vl[3]  = R[0][0]*vg[3] + R[0][1]*vg[4] + R[0][2]*vg[5];
    vl[4]  = R[1][0]*vg[3] + R[1][1]*vg[4] + R[1][2]*vg[5];
    vl[5]  = R[2][0]*vg[3] + R[2][1]*vg[4] + R[2][2]*vg[5];
      
    vl[6]  = vg[6];  // do not transform warping
    vl[7]  = vg[7];  // do not transform top flange rotation
    vl[8]  = vg[8];  // do not transform top flange curvature
    vl[9]  = vg[9];  // do not transform bottom flange rotation
    vl[10] = vg[10]; // do not transform bottom flange curvature
  
    vl[11] = R[0][0]*vg[11] + R[0][1]*vg[12] + R[0][2]*vg[13];
    vl[12] = R[1][0]*vg[11] + R[1][1]*vg[12] + R[1][2]*vg[13];
    vl[13] = R[2][0]*vg[11] + R[2][1]*vg[12] + R[2][2]*vg[13];
      
    vl[14] = R[0][0]*vg[14] + R[0][1]*vg[15] + R[0][2]*vg[16];
    vl[15] = R[1][0]*vg[14] + R[1][1]*vg[15] + R[1][2]*vg[16];
    vl[16] = R[2][0]*vg[14] + R[2][1]*vg[15] + R[2][2]*vg[16];
      
    vl[17]  = vg[17]; // do not transform warping
    vl[18]  = vg[18]; // do not transform top flange rotation
    vl[19]  = vg[19]; // do not transform top flange curvature
    vl[20]  = vg[20]; // do not transform bottom flange rotation
    vl[21]  = vg[21]; // do not transform bottom flange curvature 
       

    vb(16) = vl[11] - vl[0];
    double tmp;
    // Theta_z
    tmp = oneOverL*(vl[1]-vl[12]);
    vb(1) = vl[5] + tmp;
    vb(9) = vl[16] + tmp;
    // Theta_y
    tmp = oneOverL*(vl[13]-vl[2]);
    vb(2)  = vl[4] + tmp;
    vb(10) = vl[15] + tmp;
    // Theta_x = Torsion
    vb(0)  = (-vl[14] + vl[3])/2;
    vb(8)  = -vb(0);
    // Bi-moment
    vb(3)  = vl[6];
    vb(11) = vl[17];
    // Top Flange Rotations & Curvatures
    vb(4)  = vl[7];
    vb(5)  = vl[8];
    vb(12) = vl[18];
    vb(13) = vl[19];
    // Bottom Flange Rotations & Curvatures
    vb(6)  = vl[9];
    vb(7)  = vl[10];
    vb(14) = vl[20];
    vb(15) = vl[21]; 
	
        return vb;
}


const Vector &
PDeltaCrdTransf3d22::getBasicTrialAccel(void)
{
    // determine global accelerations
    const Vector &accel1 = nodeIPtr->getTrialAccel();
    const Vector &accel2 = nodeJPtr->getTrialAccel();
   
    static double ag[22];
    for (int i = 0; i < 11; i++) {
	ag[i]    = accel1(i);
	ag[i+11] = accel2(i);
    }
   
    double oneOverL = 1.0/L;
   
    static Vector ab(17);
   
    static double al[22];
       
    al[0]  = R[0][0]*ag[0] + R[0][1]*ag[1] + R[0][2]*ag[2];
    al[1]  = R[1][0]*ag[0] + R[1][1]*ag[1] + R[1][2]*ag[2];
    al[2]  = R[2][0]*ag[0] + R[2][1]*ag[1] + R[2][2]*ag[2];
      
    al[3]  = R[0][0]*ag[3] + R[0][1]*ag[4] + R[0][2]*ag[5];
    al[4]  = R[1][0]*ag[3] + R[1][1]*ag[4] + R[1][2]*ag[5];
    al[5]  = R[2][0]*ag[3] + R[2][1]*ag[4] + R[2][2]*ag[5];
      
    al[6]  = ag[6];  // do not transform warping
    al[7]  = ag[7];  // do not transform top flange rotation
    al[8]  = ag[8];  // do not transform top flange curvature
    al[9]  = ag[9];  // do not transform bottom flange rotation
    al[10] = ag[10]; // do not transform bottom flange curvature
  
    al[11] = R[0][0]*ag[11] + R[0][1]*ag[12] + R[0][2]*ag[13];
    al[12] = R[1][0]*ag[11] + R[1][1]*ag[12] + R[1][2]*ag[13];
    al[13] = R[2][0]*ag[11] + R[2][1]*ag[12] + R[2][2]*ag[13];
     
    al[14] = R[0][0]*ag[14] + R[0][1]*ag[15] + R[0][2]*ag[16];
    al[15] = R[1][0]*ag[14] + R[1][1]*ag[15] + R[1][2]*ag[16];
    al[16] = R[2][0]*ag[14] + R[2][1]*ag[15] + R[2][2]*ag[16];
      
    al[17]  = ag[17]; // do not transform warping
    al[18]  = ag[18]; // do not transform top flange rotation
    al[19]  = ag[19]; // do not transform top flange curvature
    al[20]  = ag[20]; // do not transform bottom flange rotation
    al[21]  = ag[21]; // do not transform bottom flange curvature 
       

    ab(16) = al[11] - al[0];
    double tmp;
    // Theta_z
    tmp = oneOverL*(al[1]-al[12]);
    ab(1) = al[5] + tmp;
    ab(9) = al[16] + tmp;
    // Theta_y
    tmp = oneOverL*(al[13]-al[2]);
    ab(2)  = al[4] + tmp;
    ab(10) = al[15] + tmp;
    // Theta_x = Torsion
    ab(0)  = (-al[14] + al[3])/2;
    ab(8)  = -ab(0);
    // Bi-moment
    ab(3)  = al[6];
    ab(11) = al[17];
    // Top Flange Rotations & Curvatures
    ab(4)  = al[7];
    ab(5)  = al[8];
    ab(12) = al[18];
    ab(13) = al[19];
    // Bottom Flange Rotations & Curvatures
    ab(6)  = al[9];
    ab(7)  = al[10];
    ab(14) = al[20];
    ab(15) = al[21]; 
       
        return ab;
}


const Vector &
PDeltaCrdTransf3d22::getGlobalResistingForce(const Vector &pb, const Vector &p0)
{
    // transform resisting forces from the basic system to local coordinates
    static double pl[22];
      
    double q0  = pb(0);  //T1
    double q1  = pb(1);  //Mz1
    double q2  = pb(2);  //My1
    double q3  = pb(3);  //Bx1
    double q4  = pb(4);  //alpha1
    double q5  = pb(5);  //alpha'1
    double q6  = pb(6);  //beta1
    double q7  = pb(7);  //beta'1
    double q8  = pb(8);  //T2
    double q9  = pb(9);  //Mz2
    double q10 = pb(10); //My2
    double q11 = pb(11); //Bx2
    double q12 = pb(12); //alpha2
    double q13 = pb(13); //alpha'2
    double q14 = pb(14); //beta2
    double q15 = pb(15); //beta'2
    double q16 = pb(16); //P 
   
    double oneOverL = 1.0/L;
   
    pl[0]  = -q16;               //P1
    pl[1]  =  oneOverL*(q1+q9);  //Vy1
    pl[2]  = -oneOverL*(q2+q10); //Vz1
    pl[3]  =  0.5*(q0-q8);       //T1
    pl[4]  =  q2;                //My1
    pl[5]  =  q1;                //Mz1
    pl[6]  =  q3;                //Bx1
    pl[7]  =  q4; //alpha1
    pl[8]  =  q5; //alpha'1
    pl[9]  =  q6; //beta1
    pl[10] =  q7; //beta'1
    pl[11] =  q16;   //P2
    pl[12] = -pl[1]; //Vy2 
    pl[13] = -pl[2]; //Vz2 
    pl[14] = -pl[3]; //T2
    pl[15] =  q10;   //My2
    pl[16] =  q9;    //Mz2
    pl[17] =  q11;   //Bx2
    pl[18] =  q12;   //alpha2
    pl[19] =  q13;   //alpha'2
    pl[20] =  q14;   //beta2
    pl[21] =  q15;   //beta'2 

   
    // Include leaning column effects (P-Delta)
    double NoverL;
    NoverL  = ul112*q16*oneOverL;            
    pl[1]  += NoverL;
    pl[12] -= NoverL;
    NoverL  = ul213*q16*oneOverL;
    pl[2]  += NoverL;
    pl[13] -= NoverL;
   
    // transform resisting forces  from local to global coordinates
    static Vector pg(22);
      
    pg(0)  = R[0][0]*pl[0] + R[1][0]*pl[1] + R[2][0]*pl[2];
    pg(1)  = R[0][1]*pl[0] + R[1][1]*pl[1] + R[2][1]*pl[2];
    pg(2)  = R[0][2]*pl[0] + R[1][2]*pl[1] + R[2][2]*pl[2];
      
    pg(3)  = R[0][0]*pl[3] + R[1][0]*pl[4] + R[2][0]*pl[5];
    pg(4)  = R[0][1]*pl[3] + R[1][1]*pl[4] + R[2][1]*pl[5];
    pg(5)  = R[0][2]*pl[3] + R[1][2]*pl[4] + R[2][2]*pl[5];
      
    pg(6)  = pl[6];
    pg(7)  = pl[7];
    pg(8)  = pl[8];
    pg(9)  = pl[9];
    pg(10) = pl[10];
  
    pg(11) = R[0][0]*pl[11] + R[1][0]*pl[12] + R[2][0]*pl[13];
    pg(12) = R[0][1]*pl[11] + R[1][1]*pl[12] + R[2][1]*pl[13];
    pg(13) = R[0][2]*pl[11] + R[1][2]*pl[12] + R[2][2]*pl[13];
      
    pg(14) = R[0][0]*pl[14] + R[1][0]*pl[15] + R[2][0]*pl[16];
    pg(15) = R[0][1]*pl[14] + R[1][1]*pl[15] + R[2][1]*pl[16];
    pg(16) = R[0][2]*pl[14] + R[1][2]*pl[15] + R[2][2]*pl[16];
  
    pg(17) = pl[17];
    pg(18) = pl[18];
    pg(19) = pl[19];
    pg(20) = pl[20];
    pg(21) = pl[21]; 
   
    return pg;
}


const Matrix &
PDeltaCrdTransf3d22::getGlobalStiffMatrix (const Matrix &KB, const Vector &pb)
{
    static Matrix kg(22,22);    // Global stiffness for return
    static double kb[17][17];   // Basic stiffness
    static double kl[22][22];   // Local stiffness
    static double tmp[22][22];  // Temporary storage
      
    double oneOverL = 1.0/L;
      
    int i,j;
    for (i = 0; i < 17; i++)
        for (j = 0; j < 17; j++)
            kb[i][j] = KB(i,j); 
       
    // Transform basic stiffness to local system
    // First compute kb*T_{bl}
    for (i = 0; i < 17; i++) {                      //basic -> local
    tmp[i][0]  = -kb[i][16];                        //e     -> x1
    tmp[i][1]  =  oneOverL*(kb[i][1]+kb[i][9]);     //rz1/2 -> y1
    tmp[i][2]  = -oneOverL*(kb[i][2]+kb[i][10]);    //ry1/2 -> z1
    tmp[i][3]  =  0.5*(kb[i][0]-kb[i][8]);          //rx1/2 -> T1
    tmp[i][4]  =  kb[i][2];             //ry1   -> ry1
    tmp[i][5]  =  kb[i][1];             //rz1   -> rz1
    tmp[i][6]  =  kb[i][3];             //Bi1   -> Bi1   
    tmp[i][7]  =  kb[i][4];             //al1   -> al1
    tmp[i][8]  =  kb[i][5];             //al'1  -> al'1
    tmp[i][9]  =  kb[i][6];             //be1   -> be1
    tmp[i][10] =  kb[i][7];             //be'1  -> be'1
  
    tmp[i][11] =  kb[i][16];            //e     -> x2
    tmp[i][12] = -tmp[i][1];            //     y2
    tmp[i][13] = -tmp[i][2];            //     z2
    tmp[i][14] = -tmp[i][3];            //     T2
    tmp[i][15] =  kb[i][10];            //ry2   -> ry2
    tmp[i][16] =  kb[i][9];             //rz2   -> rz2
    tmp[i][17] =  kb[i][11];            //Bi2   -> Bi2   
    tmp[i][18] =  kb[i][12];            //al2   -> al2
    tmp[i][19] =  kb[i][13];            //al'2  -> al'2
    tmp[i][20] =  kb[i][14];            //be2   -> be2
    tmp[i][21] =  kb[i][15];            //be'2  -> be'2
    }
  
    // Now compute T'_{bl}*(kb*T_{bl})
    for (i = 0; i < 22; i++) {
    kl[0][i]  = -tmp[16][i];
    kl[1][i]  =  oneOverL*(tmp[1][i]+tmp[9][i]);
    kl[2][i]  = -oneOverL*(tmp[2][i]+tmp[10][i]);
    kl[3][i]  =  0.5*(tmp[0][i]-tmp[8][i]);
    kl[4][i]  =  tmp[2][i];
    kl[5][i]  =  tmp[1][i];
    kl[6][i]  =  tmp[3][i];
    kl[7][i]  =  tmp[4][i];
    kl[8][i]  =  tmp[5][i];
    kl[9][i]  =  tmp[6][i];
    kl[10][i] =  tmp[7][i];
  
    kl[11][i] =  tmp[16][i];
    kl[12][i] = -kl[1][i];
    kl[13][i] = -kl[2][i];
    kl[14][i] = -kl[3][i];
    kl[15][i] =  tmp[10][i];
    kl[16][i] =  tmp[9][i];
    kl[17][i] =  tmp[11][i];
    kl[18][i] =  tmp[12][i];
    kl[19][i] =  tmp[13][i];
    kl[20][i] =  tmp[14][i];
    kl[21][i] =  tmp[15][i];
    } 
       
    // Include geometric stiffness effects in local system
    double NoverL = pb(16)*oneOverL;
    kl[1][1]   += NoverL;
    kl[2][2]   += NoverL;
    kl[12][12] += NoverL;
    kl[13][13] += NoverL;
    kl[1][12]  -= NoverL;
    kl[12][1]  -= NoverL;
    kl[2][13]  -= NoverL;
    kl[13][2]  -= NoverL;
       

       
    // Transform local stiffness to global system
    // First compute kl*T_{lg}
    int m;
    for (m = 0; m < 22; m++) {
	tmp[m][0] = kl[m][0]*R[0][0] + kl[m][1]*R[1][0]  + kl[m][2]*R[2][0];
	tmp[m][1] = kl[m][0]*R[0][1] + kl[m][1]*R[1][1]  + kl[m][2]*R[2][1];
	tmp[m][2] = kl[m][0]*R[0][2] + kl[m][1]*R[1][2]  + kl[m][2]*R[2][2];
       
	tmp[m][3] = kl[m][3]*R[0][0] + kl[m][4]*R[1][0]  + kl[m][5]*R[2][0];
	tmp[m][4] = kl[m][3]*R[0][1] + kl[m][4]*R[1][1]  + kl[m][5]*R[2][1];
	tmp[m][5] = kl[m][3]*R[0][2] + kl[m][4]*R[1][2]  + kl[m][5]*R[2][2];
       
	tmp[m][6]  = kl[m][6];
	tmp[m][7]  = kl[m][7];
	tmp[m][8]  = kl[m][8];
	tmp[m][9]  = kl[m][9];
	tmp[m][10] = kl[m][10];

	tmp[m][11] = kl[m][11]*R[0][0] + kl[m][12]*R[1][0] + kl[m][13]*R[2][0];
	tmp[m][12] = kl[m][11]*R[0][1] + kl[m][12]*R[1][1] + kl[m][13]*R[2][1];
	tmp[m][13] = kl[m][11]*R[0][2] + kl[m][12]*R[1][2] + kl[m][13]*R[2][2];
	  
	tmp[m][14] = kl[m][14]*R[0][0] + kl[m][15]*R[1][0] + kl[m][16]*R[2][0];
	tmp[m][15] = kl[m][14]*R[0][1] + kl[m][15]*R[1][1] + kl[m][16]*R[2][1];
	tmp[m][16] = kl[m][14]*R[0][2] + kl[m][15]*R[1][2] + kl[m][16]*R[2][2]; 

	tmp[m][17] = kl[m][17];
	tmp[m][18] = kl[m][18];
	tmp[m][19] = kl[m][19];
	tmp[m][20] = kl[m][20];
	tmp[m][21] = kl[m][21]; 
    }
       
    // Now compute T'_{lg}*(kl*T_{lg})
    for (m = 0; m < 22; m++) {
	kg(0,m) = R[0][0]*tmp[0][m] + R[1][0]*tmp[1][m]  + R[2][0]*tmp[2][m];
	kg(1,m) = R[0][1]*tmp[0][m] + R[1][1]*tmp[1][m]  + R[2][1]*tmp[2][m];
	kg(2,m) = R[0][2]*tmp[0][m] + R[1][2]*tmp[1][m]  + R[2][2]*tmp[2][m];
       
	kg(3,m) = R[0][0]*tmp[3][m] + R[1][0]*tmp[4][m]  + R[2][0]*tmp[5][m];
	kg(4,m) = R[0][1]*tmp[3][m] + R[1][1]*tmp[4][m]  + R[2][1]*tmp[5][m];
	kg(5,m) = R[0][2]*tmp[3][m] + R[1][2]*tmp[4][m]  + R[2][2]*tmp[5][m];
       
	kg(6,m)  = tmp[6][m];
	kg(7,m)  = tmp[7][m];
	kg(8,m)  = tmp[8][m];
	kg(9,m)  = tmp[9][m];
	kg(10,m) = tmp[10][m];
  
	kg(11,m) = R[0][0]*tmp[11][m] + R[1][0]*tmp[12][m]  + R[2][0]*tmp[13][m];
	kg(12,m) = R[0][1]*tmp[11][m] + R[1][1]*tmp[12][m]  + R[2][1]*tmp[13][m];
	kg(13,m) = R[0][2]*tmp[11][m] + R[1][2]*tmp[12][m]  + R[2][2]*tmp[13][m];
	  
	kg(14,m) = R[0][0]*tmp[14][m] + R[1][0]*tmp[15][m] + R[2][0]*tmp[16][m];
	kg(15,m) = R[0][1]*tmp[14][m] + R[1][1]*tmp[15][m] + R[2][1]*tmp[16][m];
	kg(16,m) = R[0][2]*tmp[14][m] + R[1][2]*tmp[15][m] + R[2][2]*tmp[16][m]; 

	kg(17,m) = tmp[17][m];
	kg(18,m) = tmp[18][m];
	kg(19,m) = tmp[19][m];
	kg(20,m) = tmp[20][m];
	kg(21,m) = tmp[21][m]; 
    }
       
    return kg;
}


const Matrix &
PDeltaCrdTransf3d22::getInitialGlobalStiffMatrix (const Matrix &KB)
{
    static Matrix kg(22,22);    // Global stiffness for return
    static double kb[17][17];   // Basic stiffness
    static double kl[22][22];   // Local stiffness
    static double tmp[22][22];  // Temporary storage
      
    double oneOverL = 1.0/L; 
   
    int i,j;
    for (i = 0; i < 17; i++)
        for (j = 0; j < 17; j++)
            kb[i][j] = KB(i,j);
          
    // Transform basic stiffness to local system
    // First compute kb*T_{bl}
    for (i = 0; i < 17; i++) {
	tmp[i][0]  = -kb[i][16];            	        //e     -> x1
	tmp[i][1]  =  oneOverL*(kb[i][1]+kb[i][9]);     //rz1/2 -> y1
	tmp[i][2]  = -oneOverL*(kb[i][2]+kb[i][10]);    //ry1/2 -> z1
	tmp[i][3]  =  0.5*(kb[i][0]-kb[i][8]);          //rx1/2 -> T1
	tmp[i][4]  =  kb[i][2];             //ry1   -> ry1
	tmp[i][5]  =  kb[i][1];             //rz1   -> rz1
	tmp[i][6]  =  kb[i][3];             //Bi1   -> Bi1   
	tmp[i][7]  =  kb[i][4];             //al1   -> al1
	tmp[i][8]  =  kb[i][5];             //al'1  -> al'1
	tmp[i][9]  =  kb[i][6];             //be1   -> be1
	tmp[i][10] =  kb[i][7];             //be'1  -> be'1

	tmp[i][11] =  kb[i][16];            //e     -> x2
	tmp[i][12] = -tmp[i][1];            //     y2
	tmp[i][13] = -tmp[i][2];            //     z2
	tmp[i][14] = -tmp[i][3];            //     T2
	tmp[i][15] =  kb[i][10];            //ry2   -> ry2
	tmp[i][16] =  kb[i][9];             //rz2   -> rz2
	tmp[i][17] =  kb[i][11];            //Bi2   -> Bi2   
	tmp[i][18] =  kb[i][12];            //al2   -> al2
	tmp[i][19] =  kb[i][13];            //al'2  -> al'2
	tmp[i][20] =  kb[i][14];            //be2   -> be2
	tmp[i][21] =  kb[i][15];            //be'2  -> be'2
    }
      
    // Now compute T'_{bl}*(kb*T_{bl})
    for (i = 0; i < 22; i++) {
	kl[0][i]  = -tmp[16][i];
	kl[1][i]  =  oneOverL*(tmp[1][i]+tmp[9][i]);
	kl[2][i]  = -oneOverL*(tmp[2][i]+tmp[10][i]);
	kl[3][i]  =  0.5*(tmp[0][i]-tmp[8][i]);
	kl[4][i]  =  tmp[2][i];
	kl[5][i]  =  tmp[1][i];
	kl[6][i]  =  tmp[3][i];
	kl[7][i]  =  tmp[4][i];
	kl[8][i]  =  tmp[5][i];
	kl[9][i]  =  tmp[6][i];
	kl[10][i] =  tmp[7][i];

	kl[11][i] =  tmp[16][i];
	kl[12][i] = -kl[1][i];
	kl[13][i] = -kl[2][i];
	kl[14][i] = -kl[3][i];
	kl[15][i] =  tmp[10][i];
	kl[16][i] =  tmp[9][i];
	kl[17][i] =  tmp[11][i];
	kl[18][i] =  tmp[12][i];
	kl[19][i] =  tmp[13][i];
	kl[20][i] =  tmp[14][i];
	kl[21][i] =  tmp[15][i];
    } 
   
 
   
   // Transform local stiffness to global system
    // First compute kl*T_{lg}
    int m;
    for (m = 0; m < 22; m++) {
	tmp[m][0] = kl[m][0]*R[0][0] + kl[m][1]*R[1][0]  + kl[m][2]*R[2][0];
	tmp[m][1] = kl[m][0]*R[0][1] + kl[m][1]*R[1][1]  + kl[m][2]*R[2][1];
	tmp[m][2] = kl[m][0]*R[0][2] + kl[m][1]*R[1][2]  + kl[m][2]*R[2][2];
       
	tmp[m][3] = kl[m][3]*R[0][0] + kl[m][4]*R[1][0]  + kl[m][5]*R[2][0];
	tmp[m][4] = kl[m][3]*R[0][1] + kl[m][4]*R[1][1]  + kl[m][5]*R[2][1];
	tmp[m][5] = kl[m][3]*R[0][2] + kl[m][4]*R[1][2]  + kl[m][5]*R[2][2];
       
	tmp[m][6]  = kl[m][6];
	tmp[m][7]  = kl[m][7];
	tmp[m][8]  = kl[m][8];
	tmp[m][9]  = kl[m][9];
	tmp[m][10] = kl[m][10];
  
	tmp[m][11] = kl[m][11]*R[0][0] + kl[m][12]*R[1][0] + kl[m][13]*R[2][0];
	tmp[m][12] = kl[m][11]*R[0][1] + kl[m][12]*R[1][1] + kl[m][13]*R[2][1];
	tmp[m][13] = kl[m][11]*R[0][2] + kl[m][12]*R[1][2] + kl[m][13]*R[2][2];
	  
	tmp[m][14] = kl[m][14]*R[0][0] + kl[m][15]*R[1][0] + kl[m][16]*R[2][0];
	tmp[m][15] = kl[m][14]*R[0][1] + kl[m][15]*R[1][1] + kl[m][16]*R[2][1];
	tmp[m][16] = kl[m][14]*R[0][2] + kl[m][15]*R[1][2] + kl[m][16]*R[2][2]; 

	tmp[m][17] = kl[m][17];
	tmp[m][18] = kl[m][18];
	tmp[m][19] = kl[m][19];
	tmp[m][20] = kl[m][20];
	tmp[m][21] = kl[m][21]; 
       
    }
   
    // Now compute T'_{lg}*(kl*T_{lg})
    for (m = 0; m < 22; m++) {
	kg(0,m) = R[0][0]*tmp[0][m] + R[1][0]*tmp[1][m]  + R[2][0]*tmp[2][m];
	kg(1,m) = R[0][1]*tmp[0][m] + R[1][1]*tmp[1][m]  + R[2][1]*tmp[2][m];
	kg(2,m) = R[0][2]*tmp[0][m] + R[1][2]*tmp[1][m]  + R[2][2]*tmp[2][m];
       
	kg(3,m) = R[0][0]*tmp[3][m] + R[1][0]*tmp[4][m]  + R[2][0]*tmp[5][m];
	kg(4,m) = R[0][1]*tmp[3][m] + R[1][1]*tmp[4][m]  + R[2][1]*tmp[5][m];
	kg(5,m) = R[0][2]*tmp[3][m] + R[1][2]*tmp[4][m]  + R[2][2]*tmp[5][m];
       
	kg(6,m)  = tmp[6][m];
	kg(7,m)  = tmp[7][m];
	kg(8,m)  = tmp[8][m];
	kg(9,m)  = tmp[9][m];
	kg(10,m) = tmp[10][m];
  
	kg(11,m) = R[0][0]*tmp[11][m] + R[1][0]*tmp[12][m]  + R[2][0]*tmp[13][m];
	kg(12,m) = R[0][1]*tmp[11][m] + R[1][1]*tmp[12][m]  + R[2][1]*tmp[13][m];
	kg(13,m) = R[0][2]*tmp[11][m] + R[1][2]*tmp[12][m]  + R[2][2]*tmp[13][m];
	  
	kg(14,m) = R[0][0]*tmp[14][m] + R[1][0]*tmp[15][m] + R[2][0]*tmp[16][m];
	kg(15,m) = R[0][1]*tmp[14][m] + R[1][1]*tmp[15][m] + R[2][1]*tmp[16][m];
	kg(16,m) = R[0][2]*tmp[14][m] + R[1][2]*tmp[15][m] + R[2][2]*tmp[16][m]; 

	kg(17,m) = tmp[17][m];
	kg(18,m) = tmp[18][m];
	kg(19,m) = tmp[19][m];
	kg(20,m) = tmp[20][m];
	kg(21,m) = tmp[21][m]; 
    }
   
    return kg;
}


CrdTransf *
PDeltaCrdTransf3d22::getCopy3d(void)
{
    // create a new instance of PDeltaCrdTransf3d22
   
    PDeltaCrdTransf3d22 *theCopy;
   
    static Vector xz(3);
    xz(0) = R[2][0];
    xz(1) = R[2][1];
    xz(2) = R[2][2];
   
    theCopy = new PDeltaCrdTransf3d22(this->getTag(), xz);
   
    theCopy->nodeIPtr = nodeIPtr;
    theCopy->nodeJPtr = nodeJPtr;
    theCopy->L = L;
    theCopy->ul112 = ul112;
    theCopy->ul213 = ul213;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            theCopy->R[i][j] = R[i][j];
       
       
        return theCopy;
}


int
PDeltaCrdTransf3d22::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
   
    static Vector data(27);
    data(0) = this->getTag();
    data(1) = L;
   
    if (nodeIInitialDisp != 0) {
        data(2)  = nodeIInitialDisp[0];
        data(3)  = nodeIInitialDisp[1];
        data(4)  = nodeIInitialDisp[2];
        data(5)  = nodeIInitialDisp[3];
        data(6)  = nodeIInitialDisp[4];
        data(7)  = nodeIInitialDisp[5];
	data(8)  = nodeIInitialDisp[6];
	data(9)  = nodeIInitialDisp[7];
	data(10) = nodeIInitialDisp[8];
	data(11) = nodeIInitialDisp[9];
	data(12) = nodeIInitialDisp[10];
    } else {
        data(2)  = 0.0;
        data(3)  = 0.0;
        data(4)  = 0.0;
        data(5)  = 0.0;
        data(6)  = 0.0;
        data(7)  = 0.0;
	data(8)  = 0.0;
	data(9)  = 0.0;
	data(10) = 0.0;
	data(11) = 0.0;
	data(12) = 0.0;
    }
   
    if (nodeJInitialDisp != 0) {
        data(13) = nodeJInitialDisp[0];
        data(14) = nodeJInitialDisp[1];
        data(15) = nodeJInitialDisp[2];
        data(16) = nodeJInitialDisp[3];
        data(17) = nodeJInitialDisp[4];
        data(18) = nodeJInitialDisp[5];
	data(19) = nodeJInitialDisp[6];
	data(20) = nodeJInitialDisp[7];
	data(21) = nodeJInitialDisp[8];
	data(22) = nodeJInitialDisp[9];
	data(23) = nodeJInitialDisp[10];
    } else {
        data(13) = 0.0;
        data(14) = 0.0;
        data(15) = 0.0;
        data(16) = 0.0;
        data(17) = 0.0;
        data(18) = 0.0;
	data(19) = 0.0;
	data(20) = 0.0;
	data(21) = 0.0;
	data(22) = 0.0;
	data(23) = 0.0;
    }
   
    data(24) = R[2][0];
    data(25) = R[2][1];
    data(26) = R[2][2];
   
    res += theChannel.sendVector(this->getDbTag(), cTag, data);  
    if (res < 0) {
        opserr << "PDeltaCrdTransf3d22::sendSelf - failed to send Vector\n";
        return res;
    }
   
    return res;
}


int
PDeltaCrdTransf3d22::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
   
    static Vector data(27);
   
    res += theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0) {
        opserr << "PDeltaCrdTransf3d22::recvSelf - failed to receive Vector\n";
        return res;
    }
   
    this->setTag((int)data(0));
    L = data(1);

    int flag;
    int i, j;

    flag = 0;
    for (i=2; i<=12; i++)
	if (data(i) != 0.0)
	    flag = 1;

    if (flag == 1) {
	if (nodeIInitialDisp == 0)
	    nodeIInitialDisp = new double[11];
	for (i=2, j=0; i<=12; i++, j++)
	    nodeIInitialDisp[j] = data(i);
    }
       
    flag = 0;
    for (i=13; i<=23; i++)
	if (data(i) != 0.0)
	    flag = 1;

    if (flag == 1) {
	if (nodeJInitialDisp == 0)
	    nodeJInitialDisp = new double [11];
	for (i=13, j=0; i<=23; i++, j++)
	    nodeJInitialDisp[j] = data(i);
    }
       
    R[2][0] = data(24);
    R[2][1] = data(25);
    R[2][2] = data(26);
   
    initialDispChecked = true;
                   
    return res;
}


const Vector &
PDeltaCrdTransf3d22::getPointGlobalCoordFromLocal(const Vector &xl)
{
    static Vector xg(3);
   
    xg = nodeIPtr->getCrds();
   
    if (nodeIInitialDisp != 0) {
        xg(0) -= nodeIInitialDisp[0];
        xg(1) -= nodeIInitialDisp[1];
        xg(2) -= nodeIInitialDisp[2];
    }
   
    // xg = xg + Rlj'*xl
    //xg.addMatrixTransposeVector(1.0, Rlj, xl, 1.0);
    xg(0) += R[0][0]*xl(0) + R[1][0]*xl(1) + R[2][0]*xl(2);
    xg(1) += R[0][1]*xl(0) + R[1][1]*xl(1) + R[2][1]*xl(2);
    xg(2) += R[0][2]*xl(0) + R[1][2]*xl(1) + R[2][2]*xl(2);
   
    return xg;  
}


const Vector &
PDeltaCrdTransf3d22::getPointGlobalDisplFromBasic (double xi, const Vector &uxb)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getTrialDisp();
    const Vector &disp2 = nodeJPtr->getTrialDisp();
   
    static double ug[22];
    for (int i = 0; i < 11; i++)
    {
        ug[i]    = disp1(i);
        ug[i+11] = disp2(i);
    }
   
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<11; j++)
            ug[j] -= nodeIInitialDisp[j];
    }
   
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<11; j++)
            ug[j+11] -= nodeJInitialDisp[j];
    }
   
    // transform global end displacements to local coordinates
    //ul.addMatrixVector(0.0, Tlg,  ug, 1.0);       //  ul = Tlg *  ug;
    static double ul[22];
   
    ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
    ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
    ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];
   
    ul[12]  = R[1][0]*ug[11] + R[1][1]*ug[12] + R[1][2]*ug[13];
    ul[13]  = R[2][0]*ug[11] + R[2][1]*ug[12] + R[2][2]*ug[13];
   
      
    // compute displacements at point xi, in local coordinates
    static double uxl[3];
    static Vector uxg(3);
   
    uxl[0] = uxb(0) +        ul[0];
    uxl[1] = uxb(1) + (1-xi)*ul[1] + xi*ul[12];
    uxl[2] = uxb(2) + (1-xi)*ul[2] + xi*ul[13];
   
    // rotate displacements to global coordinates
    // uxg = Rlj'*uxl
    //uxg.addMatrixTransposeVector(0.0, Rlj, uxl, 1.0);
    uxg(0) = R[0][0]*uxl[0] + R[1][0]*uxl[1] + R[2][0]*uxl[2];
    uxg(1) = R[0][1]*uxl[0] + R[1][1]*uxl[1] + R[2][1]*uxl[2];
    uxg(2) = R[0][2]*uxl[0] + R[1][2]*uxl[1] + R[2][2]*uxl[2];
   
    return uxg;  
}


void
PDeltaCrdTransf3d22::Print(OPS_Stream &s, int flag)
{
    s << "\nCrdTransf: " << this->getTag() << " Type: PDeltaCrdTransf3d22" << endln;
}
 

