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

// $Revision: 1.31 $
// $Date: 2009/09/28 22:48:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/TaperedFiberSectionSmoothing3d.cpp,v $

// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of TaperedFiberSectionSmoothing3d.
// Modified by Brighton Laiman from University of California, San Diego (include warping degrees of freedom 
// AND Local Buckling). Refer to Formulation and Implementation of Three-dimensional Beam-Column 
// Analyses with Warping Effects 

//#define DEBUG

#include <stdlib.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <TaperedFiberSectionSmoothing3d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <iostream>
#include <fstream>
using std::string;
using namespace std;


ID TaperedFiberSectionSmoothing3d::code(SEC_TAG_TaperedFiberSectionSmoothing3d);

// constructors:
TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d(int tag, int num,
                               Fiber ** fibers, double hratio, double hvalue, double g, double poissonRatio):SectionForceDeformation(tag,
                                                               SEC_TAG_TaperedFiberSectionSmoothing3d),
numFibers(num), theMaterials(0), thexAxisMaterials(0), theyAxisMaterials(0), matData(0), 
yBar(0.0), zBar(0.0), G(g), hRatio(hratio), hVal(hvalue), v(poissonRatio), e(18), eCommit(18), s(0), ks(0)
{

    if (numFibers != 0) {
        theMaterials 	  = new UniaxialMaterial *[numFibers];
        thexAxisMaterials = new UniaxialMaterial *[numFibers];
        theyAxisMaterials = new UniaxialMaterial *[numFibers];

        if (theMaterials == 0 || thexAxisMaterials == 0 || theyAxisMaterials == 0) {
            opserr <<
                "TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d -- failed to allocate Material pointers\n";
            exit(-1);
        }

        matData = new double[numFibers * 5];

        if (matData == 0) {
            opserr <<
                "TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d -- failed to allocate double array for material data\n";
            exit(-1);
        }

        double Qz = 0.0;
        double Qy = 0.0;
        double A  = 0.0;

        for (int i = 0; i < numFibers; i++) {
            Fiber *theFiber = fibers[i];
            double yLoc, zLoc, Area, tP; 
	    int PlFlag;
            theFiber->getFiberLocation(yLoc, zLoc);
	    tP     = theFiber->gettP();                // MODIFIED FOR TAPERED STUDY
	    PlFlag = theFiber->getPlFlag();	       // added for FLB study
            Area   = theFiber->getArea();

            Qz += yLoc * Area;
            Qy += zLoc * Area;
            A  += Area;

            matData[i * 5]     = yLoc;
            matData[i * 5 + 1] = zLoc;
            matData[i * 5 + 2] = Area;
	    matData[i * 5 + 3] = tP;     // plate thickness for the fiber...this assumes fibers are located at mid-thickness only!
	    matData[i * 5 + 4] = PlFlag; // Flag to distinguish which plate fiber belongs to

            UniaxialMaterial *theMat      = theFiber->getMaterial();
            UniaxialMaterial *thePlateMat = theFiber->getPlateMaterial();
            theMaterials[i] 	 = theMat->getCopy();		// used for uniaxial strain of whole cross-section (original default)
            thexAxisMaterials[i] = thePlateMat->getCopy();	// used for Mx material for flanges & web
            theyAxisMaterials[i] = thePlateMat->getCopy();	// only used for My web material

            if (theMaterials[i] == 0 || thexAxisMaterials[i] == 0 || theyAxisMaterials[i] == 0) {
                opserr <<
                    "TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d -- failed to get copy of a Material\n";
                exit(-1);
            }
        }

        yBar = Qz / A;
        zBar = Qy / A;
    }

    s  = new Vector(sData, 36);
    ks = new Matrix(kData, 12, 12);

    for (int i = 0; i < 36; i++)
	sData[i] = 0.0;

    for (int i = 0; i < 144; i++)
        kData[i] = 0.0;

    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;

    // AddingSensitivity:BEGIN ////////////////////////////////////
    parameterID = 0;
    SHVs = 0;
    // AddingSensitivity:END //////////////////////////////////////
}

// constructor for blank object that recvSelf needs to be invoked upon
TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d():
SectionForceDeformation(0, SEC_TAG_TaperedFiberSectionSmoothing3d),
numFibers(0), theMaterials(0), thexAxisMaterials(0), theyAxisMaterials(0), matData(0),
yBar(0.0), zBar(0.0), G(0.0), hRatio(0.0), hVal(0.0), v(0.0), e(18), eCommit(18), s(0), ks(0)
{
    s = new Vector(sData, 36);
    ks = new Matrix(kData, 12, 12);

    for (int i = 0; i < 36; i++)
        sData[i] = 0.0;

    for (int i = 0; i < 144; i++)
        kData[i] = 0.0;

	// what does this next piece do?
    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;

    // AddingSensitivity:BEGIN ////////////////////////////////////
    parameterID = 0;
    SHVs = 0;
    // AddingSensitivity:END //////////////////////////////////////
}

int TaperedFiberSectionSmoothing3d::addFiber(Fiber & newFiber)
{
    // need to create a larger array
    int newSize = numFibers + 1;

    UniaxialMaterial **newArray = new UniaxialMaterial *[newSize];
    double *newMatData = new double[5 * newSize];

    if (newArray == 0 || newMatData == 0) {
        opserr <<
            "TaperedFiberSectionSmoothing3d::addFiber -- failed to allocate Fiber pointers\n";
        exit(-1);
    }

    // copy the old pointers
    int i;
    for (i = 0; i < numFibers; i++) {
        newArray[i] 	      = theMaterials[i];
        newMatData[5 * i]     = matData[5 * i];
        newMatData[5 * i + 1] = matData[5 * i + 1];
        newMatData[5 * i + 2] = matData[5 * i + 2];
        newMatData[5 * i + 3] = matData[5 * i + 3];
	newMatData[5 * i + 4] = matData[5 * i + 4];
    }
    // set the new pointers
    double yLoc, zLoc, Area, tP; 
    int PlFlag;
    newFiber.getFiberLocation(yLoc, zLoc);
    tP     = newFiber.gettP();
    PlFlag = newFiber.getPlFlag();
    Area   = newFiber.getArea();

    newMatData[numFibers * 5]     = yLoc;
    newMatData[numFibers * 5 + 1] = zLoc;
    newMatData[numFibers * 5 + 2] = Area;
    newMatData[numFibers * 5 + 3] = tP;
    newMatData[numFibers * 5 + 4] = PlFlag;

    UniaxialMaterial *theMat = newFiber.getMaterial();
    newArray[numFibers] = theMat->getCopy();

    if (newArray[numFibers] == 0) {
        opserr <<
            "TaperedFiberSectionSmoothing3d::addFiber -- failed to get copy of a Material\n";
        exit(-1);

        delete[]newArray;
        delete[]newMatData;
        return -1;
    }

    numFibers++;

    if (theMaterials != 0) {
        delete[]theMaterials;
        delete[]matData;
    }

    theMaterials = newArray;
    matData = newMatData;

    double Qz = 0.0;
    double Qy = 0.0;
    double A  = 0.0;

    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
        yLoc   = matData[5 * i];
        zLoc   = matData[5 * i + 1];
        Area   = matData[5 * i + 2];

        A  += Area;
        Qz += yLoc * Area;
        Qy += zLoc * Area;
    }

    yBar = Qz / A;
    zBar = Qy / A;

    return 0;
}

// destructor:
TaperedFiberSectionSmoothing3d::~TaperedFiberSectionSmoothing3d()
{
    if (theMaterials != 0) {
        for (int i = 0; i < numFibers; i++)
            if (theMaterials[i] != 0)
                delete theMaterials[i];

        delete[]theMaterials;
    }

    if (thexAxisMaterials != 0) {
        for (int i = 0; i < numFibers; i++)
            if (thexAxisMaterials[i] != 0)
                delete thexAxisMaterials[i];

        delete[]thexAxisMaterials;
    }

    if (theyAxisMaterials != 0) {
        for (int i = 0; i < numFibers; i++)
            if (theyAxisMaterials[i] != 0)
                delete theyAxisMaterials[i];

        delete[]theyAxisMaterials;
    }

    if (matData != 0)
        delete[]matData;

    if (s != 0)
        delete s;

    if (ks != 0)
        delete ks;
}

int TaperedFiberSectionSmoothing3d::setTrialSectionDeformation(const Vector &deforms)
{
    #ifdef DEBUG
	fprintf(stdout, "TaperedFiberSectionSmoothing3d::setTrialSectionDeformation\n");
	fprintf(stdout, "The number of fibers = %3d\n", numFibers);
	fprintf(stdout, "The tag is           = %3d\n", code);
	fprintf(stdout, "The value of G       = %5.2f\n", G);
	fprintf(stdout, "The value of hVal    = %5.2f\n", hVal);
	fprintf(stdout, "The value of hRatio  = %5.2f\n", hRatio);
	fprintf(stdout, "The value of v       = %5.2f\n", v);
	fprintf(stdout, "\n");
    #endif

    int res = 0;
    e = deforms; // this is the section deformations

    for (int i =  0; i < 144; i++)
	kData[i] = 0.0;

    for (int i =  0; i < 36; i++)
	sData[i] = 0.0;

    // deforms is a vector of generalized strains
    // this is e coming from TaperedDispBeamColumnSmoothing3d
    // for smoothing the order is [u' rx1 rz1 ry1 rx2 rz2 ry2 v'' w'' theta theta' theta'' ...
    // 				   alpha alpha' alpha'' beta beta' beta'']
    double d0  = deforms(0);
    double d1  = deforms(1);
    double d2  = deforms(2);
    double d3  = deforms(3);
    double d4  = deforms(4);
    double d5  = deforms(5);
    double d6  = deforms(6);
    double d7  = deforms(7);
    double d8  = deforms(8);
    double d9  = deforms(9);
    double d10 = deforms(10);
    double d11 = deforms(11);
    double d12 = deforms(12); // alpha
    double d13 = deforms(13); // alpha'
    double d14 = deforms(14); // alpha''
    double d15 = deforms(15); // beta
    double d16 = deforms(16); // beta' 
    double d17 = deforms(17); // beta''

    int loc = 0;

    for (int i = 0; i < numFibers; i++) {
        UniaxialMaterial *theMat      = theMaterials[i];
        UniaxialMaterial *theXMat     = thexAxisMaterials[i];	 // material for flange/web x-axis bending
        UniaxialMaterial *theYMat     = theyAxisMaterials[i];	 // material for web y-axis bending
        double y   = matData[loc++];
        double z   = matData[loc++];
        double A   = matData[loc++];
	double tP  = matData[loc++];     
	int PlFlag = matData[loc++];  // Flag to distinguish which plate fiber belongs to: Top = 1, Bottom = 2, Web = 3
        double yP  = y * hRatio;      // y'
	double hP  = hVal * hRatio;   // h0'

	// calculate psi function
	double psi = 2 * yP * z;
	
        // calculate sectorial area
	double omig = y * z; // NOTE: This has to have the right sign from input

        // determine material strain and set it, include second order terms
	// this is only the axial strain, shear strain to be handled seperately
	// NOTE: These are taking into consideration smoothing effects
	double strain = d0 - y * d7 - z * d8 -omig*d11 - psi * d10 + z * d7 * d9
			-y * d8 * d9 + 0.5 * (y * y + z * z) * d10 * d10 + d3 * (d3/15 - d6/60)
			+ d6 * (d6/15 - d3/60) + d2 * (d2/15 - d5/60) + d5 * (d5/15 - d2/60);

	// shear assumed to behave elastically always
	double gamma = -tP * d10;
	double tau   = G * gamma;   	// NEED TO ADD G, same for whole section and element

	// retrieve axial tangent modulus and stress
        double tangent, stress;
        res += theMat->setTrial(strain, stress, tangent);

	//opserr << "Uniaxial Strain  = " << strain << endln;
	//opserr << "Uniaxial Stress  = " << stress << endln;

        double value   = tangent * A;
        double vas1    = y * value;
        double vas2    = z * value;
	double vas1as2 = vas1 * z;

        // section stiffness matrix k, refer to Alemdar
	kData[0]  += value;
	kData[1]  += -1 * vas1;
	kData[2]  += vas2;
	kData[3]  += (y * y + z * z) * value;
	kData[4]  += -omig * value;
	kData[5]  += -psi * value;
	kData[12] = kData[1];
	kData[13] += y * vas1;
	kData[14] += -1 * vas1as2;
	kData[15] += -1 * vas1 * (y * y + z * z);
	kData[16] += vas1 * omig;
	kData[17] += vas1 * psi;
	kData[24] = kData[2];
	kData[25] = kData[14];
	kData[26] += z * vas2; 
	kData[27] += vas2 * (y * y + z * z);
	kData[28] += -vas2 * omig;
	kData[29] += -vas2 * psi;
	kData[36] = kData[3];
	kData[37] = kData[15];
	kData[38] = kData[27];
	kData[39] += (y * y + z * z) * (y * y + z * z) * value;
	kData[40] += -omig * value * (y * y + z * z);
	kData[41] += -psi * value * (y * y + z * z);
	kData[48] = kData[4];
	kData[49] = kData[16];
	kData[50] = kData[28];
	kData[51] = kData[40];
	kData[52] += omig * omig * value;
	kData[53] += omig * psi * value;
	kData[60] = kData[5];
	kData[61] = kData[17];
	kData[62] = kData[29];
	kData[63] = kData[41];
	kData[64] = kData[53];
	kData[65] += psi * psi * value + (G * A * tP * tP)/3;  // the second term replaces GJ

        // Local Buckling additions
	// Fiber Section Constants
	// eta = constant for each fiber 
	double eta  =  y / hVal;
	double eta2 =  eta * eta;
	double eta3 =  eta2 * eta;
	double C1   =  1/8 - eta/4 - eta2/2 + eta3;
	double C2   = -1/8 - eta/4 + eta2/2 + eta3;
	double C3   = -1/4 - eta + 3 * eta2;
	double C4   = -1/4 + eta + 3 * eta2;
	double dC3  = (-1 + 6 * eta) / hVal;
	double dC4  = (1 + 6 * eta) / hVal;
	double h2   =  hVal * hVal;
	double w    =  A / tP; 			// width of fiber

	double DTx_tan, M_Tx, k_Tx;
	double DBx_tan, M_Bx, k_Bx;
	double DWx_tan, DWy_tan, M_Wx, M_Wy, k_Wx, k_Wy;

	if (PlFlag == 1) { // Top Flange Case
	    // material curvatures for plate theory
	    k_Tx  = -z * d14;
	    res += theXMat->setTrial(k_Tx, M_Tx, DTx_tan);

	    kData[104] += -DTx_tan * z * z * w;
	} else if (PlFlag == 2) { // Bottom Flange Case
	    // material curvatures for plate theory
	    k_Bx  = -z * d17;
	    res += theXMat->setTrial(k_Bx, M_Bx, DBx_tan);

	    kData[143] += -DBx_tan * z * z * w;
	} else if (PlFlag == 3) { // Web Plate Case
	    // material curvatures for plate theory
	    k_Wx  = (2 * hP * C2) * d13 + (hVal * C2) * d14 + (2 * hP * C1) * d16 + (hVal * C1) * d17;
	    k_Wy  = (dC4 / hVal) * d12 + (dC3 / hVal) * d15;
	    res += theXMat->setTrial(k_Wx, M_Wx, DWx_tan);
	    res += theYMat->setTrial(k_Wy, M_Wy, DWy_tan);

	    kData[78]  += -DWy_tan * dC4 * dC4 * w/h2;
	    kData[81]  += -DWy_tan * dC3 * dC4 * w/h2;
	    kData[91]  += -DWx_tan * (2 * hP * C2) * (2 * hP * C2) * w;
	    kData[92]  += -DWx_tan * 2 * hVal * hP * C2 * C2 * w;
	    kData[94]  += -DWx_tan * (2 * hP) * (2 * hP) * C1 * C2 * w;
	    kData[95]  += -DWx_tan * 2 * hVal * hP * C1 * C2 * w;
	    kData[103]  =  kData[92];
	    kData[104] += -DWx_tan * (hVal * C2) * (hVal * C2) * w;
	    kData[106] += -DWx_tan * 2 * hVal * hP * C1 * C2 * w;
	    kData[107] += -DWx_tan * h2 * C1 * C2 * w;
	    kData[114]  =  kData[81];
	    kData[117] += -DWy_tan * dC3 * dC3 * w/h2;
	    kData[126]  =  kData[82];
	    kData[127]  =  kData[94]; 
	    kData[128]  =  kData[106];
	    kData[130] += -DWx_tan * (2 * hP * C1) * (2 * hP * C1) * w;
	    kData[131] += -DWx_tan * 2 * hVal * hP * C1 * C1 * w;
	    kData[139]  =  kData[95];
	    kData[140]  =  kData[107];
	    kData[142]  =  kData[131];
	    kData[143] += -DWx_tan * (hVal * C1) * (hVal * C1) * w; 
	}

        // section force vector D, refer to Alemdar
        double fs0  = stress * A;
	double fs0A = fs0 * A;

        sData[0] +=  fs0;
        sData[1] += -1.0 * fs0 * y;
        sData[2] +=  1.0 * fs0 * z;
        sData[3] +=  fs0 * (y * y + z * z);
        sData[4] += -fs0 * omig;
	sData[5] += -fs0 * psi - (tau * A * tP)/3;  // second term adds pure st. venant torsion
	// Section Forces for Plate Bending AND Gmax Components
	if (PlFlag == 1) { // Top Flange Case
	    sData[8] += -z * M_Tx * w;

	    // Gmax Components
	    sData[12] += -fs0A * z;
	    sData[23] +=  fs0A * z * z;
	    sData[30] +=  fs0A * z * z;
	} else if (PlFlag == 2) { // Bottom Flange Case
	    sData[11] += -z * M_Bx * w;

	    // Gmax Components
	    sData[13] += -fs0A * z ;
	    sData[25] +=  fs0A * z * z;
	    sData[35] +=  fs0A * z * z;
	} else if (PlFlag == 3) { // Web Plate Case
	    sData[6]  += (dC4 * M_Wy / hVal) * w;
	    sData[7]  += 2 * hP * C2 * M_Wx * w;
	    sData[8]  += hVal * C2 * M_Wx * w;
	    sData[9]  += (dC3 * M_Wy / hVal) * w;
	    sData[10] += 2 * hP * C1 * M_Wx * w;
	    sData[11] += hVal * C1 * M_Wx * w;

	    // Gmax Components
	    sData[14] += fs0A * hP * C2;
	    sData[15] += fs0A * hVal * C2;
	    sData[16] += fs0A * hP * C1;
	    sData[17] += fs0A * hVal * C1;
	    sData[18] += fs0A * hP * yP * C2;
	    sData[19] += fs0A * hVal * yP * C2;
	    sData[20] += fs0A * hP * yP * C1;
	    sData[21] += fs0A * hVal * yP * C1;
	    sData[22] += fs0A * y * hP * C2;
	    sData[23] += fs0A * y * hVal * C2;
	    sData[24] += fs0A * y * hP * C1;
	    sData[25] += fs0A * y * hVal * C1;
	    sData[26] += fs0A * hP * hP * C2 * C2;
	    sData[27] += fs0A * hVal * hP * C2 * C2;
	    sData[28] += fs0A * hP * hP * C1 * C2;
	    sData[29] += fs0A * hVal * hP * C1 * C2;
	    sData[30] += fs0A * h2 * C2 * C2;
	    sData[31] += fs0A * hVal * hP * C1 * C2;
	    sData[32] += fs0A * h2 * C1 * C2;
	    sData[33] += fs0A * hP * hP * C1 * C1;
	    sData[34] += fs0A * hVal * hP * C1 * C1;
	    sData[35] += fs0A * h2 * C1 * C1;
	}

    }

    return res;
}

const Matrix &TaperedFiberSectionSmoothing3d::getInitialTangent(void)
{
    
    static double kInitialData[144];
    static Matrix kInitial(kInitialData, 12, 12);

    for (int i = 0; i < 144; i++)
        kInitialData[i] = 0.0;

    int loc = 0;

    for (int i = 0; i < numFibers; i++) {
        UniaxialMaterial *theMat      = theMaterials[i];
        UniaxialMaterial *theXMat     = thexAxisMaterials[i];	 // material for flange/web x-axis bending
        UniaxialMaterial *theYMat     = theyAxisMaterials[i];	 // material for web y-axis bending
        double y   = matData[loc++];
        double z   = matData[loc++];
        double A   = matData[loc++];
	double tP  = matData[loc++];
	int PlFlag = matData[loc++]; // Flag to distinguish which plate fiber belongs to: Top = 1, Bottom = 2, Web = 3
        double yP  = y * hRatio;     // y'
	double hP  = hVal * hRatio;  // h0'

        // calculate sectorial area
        double omig = y * z;

	// calculate psi function
	double psi = 2 * yP * z;

        double tangent = theMat->getInitialTangent();

        double value   = tangent * A;
        double vas1    = y * value;
        double vas2    = z * value;
        double vas1as2 = vas1 * z;

        // section stiffness matrix k, refer to Alemdar
	kInitialData[0]  += value;
	kInitialData[1]  += -1*vas1;
	kInitialData[2]  += vas2;
	kInitialData[3]  += (y*y+z*z)*value;
	kInitialData[4]  += -omig*value;
	kInitialData[5]  += -psi*value;
	kInitialData[12]  = kInitialData[1];
	kInitialData[13]  += y*vas1;
	kInitialData[14]  += -1*vas1as2;
	kInitialData[15]  += -1*vas1*(y*y+z*z);
	kInitialData[16] += vas1*omig;
	kInitialData[17] += vas1*psi;
	kInitialData[24] = kInitialData[2];
	kInitialData[25] = kInitialData[14];
	kInitialData[26] += z*vas2; 
	kInitialData[27] += vas2*(y*y+z*z);
	kInitialData[28] += -vas2*omig;
	kInitialData[29] += -vas2*psi;
	kInitialData[36] = kInitialData[3];
	kInitialData[37] = kInitialData[15];
	kInitialData[38] = kInitialData[27];
	kInitialData[39] += (y*y+z*z)*(y*y+z*z)*value;
	kInitialData[40] += -omig*value*(y*y+z*z);
	kInitialData[41] += -psi*value*(y*y+z*z);
	kInitialData[48] = kInitialData[4];
	kInitialData[49] = kInitialData[16];
	kInitialData[50] = kInitialData[28];
	kInitialData[51] = kInitialData[40];
	kInitialData[52] += omig*omig*value;
	kInitialData[53] += omig*psi*value;
	kInitialData[60] = kInitialData[5];
	kInitialData[61] = kInitialData[17];
	kInitialData[62] = kInitialData[29];
	kInitialData[63] = kInitialData[41];
	kInitialData[64] = kInitialData[53];
	kInitialData[65] += psi * psi * value + (G * A * tP * tP)/3;

        // Local Buckling additions 
	// Fiber Section Constants
	// eta = constant for each fiber 
	double eta  = y / hVal;
	double eta2 = eta * eta;
	double eta3 = eta2 * eta;
	double C1   = 1/8 - eta/4 - eta2/2 + eta3;
	double C2   = -1/8 - eta/4 + eta2/2 + eta3;
	double C3   = -1/4 - eta + 3 * eta2;
	double C4   = -1/4 + eta + 3 * eta2;
	double dC3  = (-1 + 6 * eta) / hVal;
	double dC4  = (1 + 6 * eta) / hVal;
	double h2   = hVal * hVal;
	double w    = A / tP; 			// width of fiber

	if (PlFlag == 1) { // Top Flange Case
	    double DTx_tan = theXMat->getInitialTangent(); 

	    kInitialData[104] += -DTx_tan * z * z * w;
	} else if (PlFlag == 2) { // Bottom Flange Case
	    double DBx_tan = theXMat->getInitialTangent();  

	    kInitialData[143] += -DBx_tan * z * z * w;
	} else if (PlFlag == 3) { // Web Plate Case
	    double DWx_tan = theXMat->getInitialTangent(); 
	    double DWy_tan = theYMat->getInitialTangent(); 

	    kInitialData[78]  += -DWy_tan * dC4 * dC4 * w/h2;
	    kInitialData[81]  += -DWy_tan * dC3 * dC4 * w/h2;
	    kInitialData[91]  += -DWx_tan * (2 * hP * C2) * (2 * hP * C2) * w;
	    kInitialData[92]  += -DWx_tan * 2 * hVal * hP * C2 * C2 * w;
	    kInitialData[94]  += -DWx_tan * (2 * hP) * (2 * hP) * C1 * C2 * w;
	    kInitialData[95]  += -DWx_tan * 2 * hVal * hP * C1 * C2 * w;
	    kInitialData[103]  =  kInitialData[92];
	    kInitialData[104] += -DWx_tan * (hVal * C2) * (hVal * C2) * w;
	    kInitialData[106] += -DWx_tan * 2 * hVal * hP * C1 * C2 * w;
	    kInitialData[107] += -DWx_tan * h2 * C1 * C2 * w;
	    kInitialData[114]  =  kInitialData[81];
	    kInitialData[117] += -DWy_tan * dC3 * dC3 * w/h2;
	    kInitialData[126]  =  kInitialData[82];
	    kInitialData[127]  =  kInitialData[94]; 
	    kInitialData[128]  =  kInitialData[106];
	    kInitialData[130] += -DWx_tan * (2 * hP * C1) * (2 * hP * C1) * w;
	    kInitialData[131] += -DWx_tan * 2 * hVal * hP * C1 * C1 * w;
	    kInitialData[139]  =  kInitialData[95];
	    kInitialData[140]  =  kInitialData[107];
	    kInitialData[142]  =  kInitialData[131];
	    kInitialData[143] += -DWx_tan * (hVal * C1) * (hVal * C1) * w; 
	}
    }

    return kInitial; 
}

const Vector & TaperedFiberSectionSmoothing3d::getSectionDeformation(void)
{
    return e;
}

const Matrix & TaperedFiberSectionSmoothing3d::getSectionTangent(void)
{
    return *ks;
}

const Vector & TaperedFiberSectionSmoothing3d::getStressResultant(void)
{
    return *s;
}

SectionForceDeformation *TaperedFiberSectionSmoothing3d::getCopy(void)
{

    TaperedFiberSectionSmoothing3d *theCopy = new TaperedFiberSectionSmoothing3d();
    theCopy->setTag(this->getTag());

    theCopy->numFibers = numFibers;

    if (numFibers != 0) {
        theCopy->theMaterials        = new UniaxialMaterial *[numFibers];
        theCopy->thexAxisMaterials   = new UniaxialMaterial *[numFibers];
        theCopy->theyAxisMaterials   = new UniaxialMaterial *[numFibers];

        if (theCopy->theMaterials == 0 || theCopy->thexAxisMaterials == 0 || theCopy->theyAxisMaterials == 0) {
            opserr <<
                "TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d -- failed to allocate Material pointers\n";
            exit(-1);
        }

        theCopy->matData = new double[numFibers * 5];

        if (theCopy->matData == 0) {
            opserr <<
                "TaperedFiberSectionSmoothing3d::TaperedFiberSectionSmoothing3d -- failed to allocate double array for material data\n";
            exit(-1);
        }

        for (int i = 0; i < numFibers; i++) {
            theCopy->matData[i * 5]     = matData[i * 5];
            theCopy->matData[i * 5 + 1] = matData[i * 5 + 1];
            theCopy->matData[i * 5 + 2] = matData[i * 5 + 2];
            theCopy->matData[i * 5 + 3] = matData[i * 5 + 3];
            theCopy->matData[i * 5 + 4] = matData[i * 5 + 4];
            theCopy->theMaterials[i] 	  = theMaterials[i]->getCopy();
            theCopy->thexAxisMaterials[i] = thexAxisMaterials[i]->getCopy();
            theCopy->theyAxisMaterials[i] = theyAxisMaterials[i]->getCopy();

            if (theCopy->theMaterials[i] == 0 || theCopy->thexAxisMaterials[i] == 0 || theCopy->theyAxisMaterials[i] == 0) {
                opserr <<
                    "TaperedFiberSectionSmoothing3d::getCopy -- failed to get copy of a Material\n";
                exit(-1);
            }
        }
    }

    theCopy->eCommit = eCommit;
    theCopy->e       = e;
    theCopy->yBar    = yBar;
    theCopy->zBar    = zBar;
    theCopy->G       = G;
    theCopy->hRatio  = hRatio;
    theCopy->hVal    = hVal;
    theCopy->v 	     = v;

    for (int i = 0; i < 144; i++)
        theCopy->kData[i] = kData[i];

    for (int i = 0; i < 36; i++)
	theCopy->sData[i] = sData[i];

    return theCopy;
}

const ID & TaperedFiberSectionSmoothing3d::getType()
{
    return code;
}

int TaperedFiberSectionSmoothing3d::getOrder() const const
{
    return 12; 
}

int TaperedFiberSectionSmoothing3d::commitState(void)
{
    
    int err = 0;

    for (int i = 0; i < numFibers; i++) {
        err += theMaterials[i]->commitState();
        err += thexAxisMaterials[i]->commitState();
        err += theyAxisMaterials[i]->commitState();
    }

    eCommit = e;

    //opserr << "commited strains = " << eCommit << endln;

    return err;
}

int TaperedFiberSectionSmoothing3d::revertToLastCommit(void)
{
    
    int err = 0;

    // Last committed section deformations
    e = eCommit;

    for (int i = 0; i < 144; i++)
	kData[i] = 0.0;

    for (int i = 0; i < 36; i++)
	sData[i] = 0.0;

    int loc = 0;

    for (int i = 0; i < numFibers; i++) {
        UniaxialMaterial *theMat      = theMaterials[i];
        UniaxialMaterial *theXMat     = thexAxisMaterials[i];	 // material for flange/web x-axis bending
        UniaxialMaterial *theYMat     = theyAxisMaterials[i];	 // material for web y-axis bending
        double y   = matData[loc++];
        double z   = matData[loc++];
        double A   = matData[loc++];
        double tP  = matData[loc++];  
	int PlFlag = matData[loc++];    // Flag to distinguish which plate fiber belongs to: Top = 1, Bottom = 2, Web = 3
        double yP  = y * hRatio;	// y'
        double hP  = hVal * hRatio;	// h0'

        // invoke revertToLast on the material
        err += theMat->revertToLastCommit();

        double tangent = theMat->getTangent();
        double stress  = theMat->getStress();

	double gamma = -tP * e[10];
	double tau   = G * gamma;   // NEED TO ADD G, same for whole section and element

        double value   = tangent * A;
        double vas1    = y * value;
        double vas2    = z * value;
        double vas1as2 = vas1 * z;
        double omig    = y * z;
	
	// calculate psi function
	double psi = 2 * yP * z;

	kData[0]  += value;
	kData[1]  += -1 * vas1;
	kData[2]  += vas2;
	kData[3]  += (y * y + z * z) * value;
	kData[4]  += -omig * value;
	kData[5]  += -psi * value;
	kData[12] = kData[1];
	kData[13] += y * vas1;
	kData[14] += -1 * vas1as2;
	kData[15] += -1 * vas1 * (y * y + z * z);
	kData[16] += vas1 * omig;
	kData[17] += vas1 * psi;
	kData[24] = kData[2];
	kData[25] = kData[14];
	kData[26] += z * vas2; 
	kData[27] += vas2 * (y * y + z * z);
	kData[28] += -vas2 * omig;
	kData[29] += -vas2 * psi;
	kData[36] = kData[3];
	kData[37] = kData[15];
	kData[38] = kData[27];
	kData[39] += (y * y + z * z) * (y * y + z * z) * value;
	kData[40] += -omig * value * (y * y + z * z);
	kData[41] += -psi * value * (y * y + z * z);
	kData[48] = kData[4];
	kData[49] = kData[16];
	kData[50] = kData[28];
	kData[51] = kData[40];
	kData[52] += omig * omig * value;
	kData[53] += omig * psi * value;
	kData[60] = kData[5];
	kData[61] = kData[17];
	kData[62] = kData[29];
	kData[63] = kData[41];
	kData[64] = kData[53];
	kData[65] += psi * psi * value + (G * A * tP * tP)/3;  // the second term replaces GJ

        // Local Buckling addition
	// Fiber Section Constants
	// eta = constant for each fiber 
	double eta  = y / hVal;
	double eta2 = eta * eta;
	double eta3 = eta2 * eta;
	double C1   = 1/8 - eta/4 - eta2/2 + eta3;
	double C2   = -1/8 - eta/4 + eta2/2 + eta3;
	double C3   = -1/4 - eta + 3 * eta2;
	double C4   = -1/4 + eta + 3 * eta2;
	double dC3  = (-1 + 6 * eta) / hVal;
	double dC4  = (1 + 6 * eta) / hVal;
	double h2   = hVal * hVal;
	double w    = A / tP; 			// width of fiber

	double DTx_tan, M_Tx;
	double DBx_tan, M_Bx;
	double DWx_tan, DWy_tan, M_Wx, M_Wy;

	if (PlFlag == 1) { // Top Flange Case
	    err += theXMat->revertToLastCommit();
	    DTx_tan = theXMat->getTangent();
	    M_Tx    = theXMat->getStress();

	    kData[104] += -DTx_tan * z * z * w;
	} else if (PlFlag == 2) { // Bottom Flange Case
	    err += theXMat->revertToLastCommit();
	    DBx_tan = theXMat->getTangent();
	    M_Bx    = theXMat->getStress();

	    kData[143] += -DBx_tan * z * z * w;
	} else if (PlFlag == 3) { // Web Plate Case
	    err += theXMat->revertToLastCommit();
	    DWx_tan = theXMat->getTangent();
	    M_Wx    = theXMat->getStress();

	    err += theYMat->revertToLastCommit();
	    DWy_tan = theYMat->getTangent();
	    M_Wy    = theYMat->getStress();

	    kData[78]  += -DWy_tan * dC4 * dC4 * w/h2;
	    kData[81]  += -DWy_tan * dC3 * dC4 * w/h2;
	    kData[91]  += -DWx_tan * (2 * hP * C2) * (2 * hP * C2) * w;
	    kData[92]  += -DWx_tan * 2 * hVal * hP * C2 * C2 * w;
	    kData[94]  += -DWx_tan * (2 * hP) * (2 * hP) * C1 * C2 * w;
	    kData[95]  += -DWx_tan * 2 * hVal * hP * C1 * C2 * w;
	    kData[103]  =  kData[92];
	    kData[104] += -DWx_tan * (hVal * C2) * (hVal * C2) * w;
	    kData[106] += -DWx_tan * 2 * hVal * hP * C1 * C2 * w;
	    kData[107] += -DWx_tan * h2 * C1 * C2 * w;
	    kData[114]  =  kData[81];
	    kData[117] += -DWy_tan * dC3 * dC3 * w/h2;
	    kData[126]  =  kData[82];
	    kData[127]  =  kData[94]; 
	    kData[128]  =  kData[106];
	    kData[130] += -DWx_tan * (2 * hP * C1) * (2 * hP * C1) * w;
	    kData[131] += -DWx_tan * 2 * hVal * hP * C1 * C1 * w;
	    kData[139]  =  kData[95];
	    kData[140]  =  kData[107];
	    kData[142]  =  kData[131];
	    kData[143] += -DWx_tan * (hVal * C1) * (hVal * C1) * w; 
	}

        double fs0  = stress * A;  
	double fs0A = fs0 * A;

        sData[0] +=  fs0;
        sData[1] += -1.0 * fs0 * y;
        sData[2] +=  1.0 * fs0 * z;
        sData[3] +=  fs0 * (y * y + z * z);
        sData[4] += -fs0 * omig;
	sData[5] += -fs0 * psi - (tau * A * tP)/3;
	// Section Forces for Plate Bending AND Gmax Components
	if (PlFlag == 1) { // Top Flange Case
	    sData[8] += -z * M_Tx * w;

	    // Gmax Components
	    sData[12] += -fs0A * z;
	    sData[23] +=  fs0A * z * z;
	    sData[30] +=  fs0A * z * z;
	} else if (PlFlag == 2) { // Bottom Flange Case
	    sData[11] += -z * M_Bx * w;

	    // Gmax Components
	    sData[13] += -fs0A * z ;
	    sData[25] +=  fs0A * z * z;
	    sData[35] +=  fs0A * z * z;
	} else if (PlFlag == 3) { // Web Plate Case
	    sData[6]  += (dC4 * M_Wy / hVal) * w;
	    sData[7]  += 2 * hP * C2 * M_Wx * w;
	    sData[8]  += hVal * C2 * M_Wx * w;
	    sData[9]  += (dC3 * M_Wy / hVal) * w;
	    sData[10] += 2 * hP * C1 * M_Wx * w;
	    sData[11] += hVal * C1 * M_Wx * w;

	    // Gmax Components
	    sData[14] += fs0A * hP * C2;
	    sData[15] += fs0A * hVal * C2;
	    sData[16] += fs0A * hP * C1;
	    sData[17] += fs0A * hVal * C1;
	    sData[18] += fs0A * hP * yP * C2;
	    sData[19] += fs0A * hVal * yP * C2;
	    sData[20] += fs0A * hP * yP * C1;
	    sData[21] += fs0A * hVal * yP * C1;
	    sData[22] += fs0A * y * hP * C2;
	    sData[23] += fs0A * y * hVal * C2;
	    sData[24] += fs0A * y * hP * C1;
	    sData[25] += fs0A * y * hVal * C1;
	    sData[26] += fs0A * hP * hP * C2 * C2;
	    sData[27] += fs0A * hVal * hP * C2 * C2;
	    sData[28] += fs0A * hP * hP * C1 * C2;
	    sData[29] += fs0A * hVal * hP * C1 * C2;
	    sData[30] += fs0A * h2 * C2 * C2;
	    sData[31] += fs0A * hVal * hP * C1 * C2;
	    sData[32] += fs0A * h2 * C1 * C2;
	    sData[33] += fs0A * hP * hP * C1 * C1;
	    sData[34] += fs0A * hVal * hP * C1 * C1;
	    sData[35] += fs0A * h2 * C1 * C1;
	}

    }

    return err;
}

int TaperedFiberSectionSmoothing3d::revertToStart(void)
{
    // revert the fibers to start    
    int err = 0;

    for (int i = 0; i < 144; i++)
	kData[i] = 0.0;

    for (int i = 0; i < 36; i++)
	sData[i] = 0.0;

    int loc = 0;

    for (int i = 0; i < numFibers; i++) {
        UniaxialMaterial *theMat      = theMaterials[i];
        UniaxialMaterial *theXMat     = thexAxisMaterials[i];	 // material for flange/web x-axis bending
        UniaxialMaterial *theYMat     = theyAxisMaterials[i];	 // material for web y-axis bending
        double y      = matData[loc++];
        double z      = matData[loc++];
        double A      = matData[loc++];
        double tP     = matData[loc++]; 
	int PlFlag    = matData[loc++]; // Flag to distinguish which plate fiber belongs to: Top = 1, Bottom = 2, Web = 3
        double yP     = y * hRatio;	// y'
        double hP     = hVal * hRatio;	// h0'

        double omig = y * z;

	// calculate psi function 
	double psi = 2 * yP * z;

        // invoke revertToStart on the material
        err += theMat->revertToStart();

        double tangent = theMat->getTangent();
        double stress  = theMat->getStress();

	double gamma = -tP * e[10];

	double tau = G * gamma;   

        double value   = tangent * A;
        double vas1    = y * value;
        double vas2    = z * value;
        double vas1as2 = vas1 * z;

	kData[0]  += value;
	kData[1]  += -1 * vas1;
	kData[2]  += vas2;
	kData[3]  += (y * y + z * z) * value;
	kData[4]  += -omig * value;
	kData[5]  += -psi * value;
	kData[12] = kData[1];
	kData[13] += y * vas1;
	kData[14] += -1 * vas1as2;
	kData[15] += -1 * vas1 * (y * y + z * z);
	kData[16] += vas1 * omig;
	kData[17] += vas1 * psi;
	kData[24] = kData[2];
	kData[25] = kData[14];
	kData[26] += z * vas2; 
	kData[27] += vas2 * (y * y + z * z);
	kData[28] += -vas2 * omig;
	kData[29] += -vas2 * psi;
	kData[36] = kData[3];
	kData[37] = kData[15];
	kData[38] = kData[27];
	kData[39] += (y * y + z * z) * (y * y + z * z) * value;
	kData[40] += -omig * value * (y * y + z * z);
	kData[41] += -psi * value * (y * y + z * z);
	kData[48] = kData[4];
	kData[49] = kData[16];
	kData[50] = kData[28];
	kData[51] = kData[40];
	kData[52] += omig * omig * value;
	kData[53] += omig * psi * value;
	kData[60] = kData[5];
	kData[61] = kData[17];
	kData[62] = kData[29];
	kData[63] = kData[41];
	kData[64] = kData[53];
	kData[65] += psi * psi * value + (G * A * tP * tP)/3;  // the second term replaces GJ

        // Local Buckling Addition
	// Fiber Section Constants
	// eta = constant for each fiber 
	double eta  = y / hVal;
	double eta2 = eta * eta;
	double eta3 = eta2 * eta;
	double C1   = 1/8 - eta/4 - eta2/2 + eta3;
	double C2   = -1/8 - eta/4 + eta2/2 + eta3;
	double C3   = -1/4 - eta + 3 * eta2;
	double C4   = -1/4 + eta + 3 * eta2;
	double dC3  = (-1 + 6 * eta) / hVal;
	double dC4  = (1 + 6 * eta) / hVal;
	double h2   = hVal * hVal;
	double w    = A / tP; 			// width of fiber

	double DTx_tan, M_Tx;
	double DBx_tan, M_Bx;
	double DWx_tan, DWy_tan, M_Wx, M_Wy;

	if (PlFlag == 1) { // Top Flange Case
	    err += theXMat->revertToStart();
	    DTx_tan = theXMat->getTangent();
	    M_Tx    = theXMat->getStress();

	    kData[104] += -DTx_tan * z * z * w;
	} else if (PlFlag == 2) { // Bottom Flange Case
	    err += theXMat->revertToStart();
	    DBx_tan = theXMat->getTangent();
	    M_Bx    = theXMat->getStress();

	    kData[143] += -DBx_tan * z * z * w;
	} else if (PlFlag == 3) { // Web Plate Case
	    err += theXMat->revertToStart();
	    DWx_tan = theXMat->getTangent();
	    M_Wx    = theXMat->getStress();

	    err += theYMat->revertToStart();
	    DWy_tan = theYMat->getTangent();
	    M_Wy    = theYMat->getStress();

	    kData[78]  += -DWy_tan * dC4 * dC4 * w/h2;
	    kData[81]  += -DWy_tan * dC3 * dC4 * w/h2;
	    kData[91]  += -DWx_tan * (2 * hP * C2) * (2 * hP * C2) * w;
	    kData[92]  += -DWx_tan * 2 * hVal * hP * C2 * C2 * w;
	    kData[94]  += -DWx_tan * (2 * hP) * (2 * hP) * C1 * C2 * w;
	    kData[95]  += -DWx_tan * 2 * hVal * hP * C1 * C2 * w;
	    kData[103]  =  kData[92];
	    kData[104] += -DWx_tan * (hVal * C2) * (hVal * C2) * w;
	    kData[106] += -DWx_tan * 2 * hVal * hP * C1 * C2 * w;
	    kData[107] += -DWx_tan * h2 * C1 * C2 * w;
	    kData[114]  =  kData[81];
	    kData[117] += -DWy_tan * dC3 * dC3 * w/h2;
	    kData[126]  =  kData[82];
	    kData[127]  =  kData[94]; 
	    kData[128]  =  kData[106];
	    kData[130] += -DWx_tan * (2 * hP * C1) * (2 * hP * C1) * w;
	    kData[131] += -DWx_tan * 2 * hVal * hP * C1 * C1 * w;
	    kData[139]  =  kData[95];
	    kData[140]  =  kData[107];
	    kData[142]  =  kData[131];
	    kData[143] += -DWx_tan * (hVal * C1) * (hVal * C1) * w; 
	}

        double fs0 = stress * A;  
	double fs0A = fs0 * A;

        sData[0] +=  fs0;
        sData[1] += -1.0 * fs0 * y;
        sData[2] +=  1.0 * fs0 * z;
        sData[3] +=  fs0 * (y * y + z * z);
        sData[4] += -fs0 * omig;
	sData[5] += -fs0 * psi - (tau * A * tP)/3;
	// Section Forces for Plate Bending AND Gmax Components
	if (PlFlag == 1) { // Top Flange Case
	    sData[8] += -z * M_Tx * w;

	    // Gmax Components
	    sData[12] += -fs0A * z;
	    sData[23] +=  fs0A * z * z;
	    sData[30] +=  fs0A * z * z;
	} else if (PlFlag == 2) { // Bottom Flange Case
	    sData[11] += -z * M_Bx * w;

	    // Gmax Components
	    sData[13] += -fs0A * z ;
	    sData[25] +=  fs0A * z * z;
	    sData[35] +=  fs0A * z * z;
	} else if (PlFlag == 3) { // Web Plate Case
	    sData[6]  += (dC4 * M_Wy / hVal) * w;
	    sData[7]  += 2 * hP * C2 * M_Wx * w;
	    sData[8]  += hVal * C2 * M_Wx * w;
	    sData[9]  += (dC3 * M_Wy / hVal) * w;
	    sData[10] += 2 * hP * C1 * M_Wx * w;
	    sData[11] += hVal * C1 * M_Wx * w;

	    // Gmax Components
	    sData[14] += fs0A * hP * C2;
	    sData[15] += fs0A * hVal * C2;
	    sData[16] += fs0A * hP * C1;
	    sData[17] += fs0A * hVal * C1;
	    sData[18] += fs0A * hP * yP * C2;
	    sData[19] += fs0A * hVal * yP * C2;
	    sData[20] += fs0A * hP * yP * C1;
	    sData[21] += fs0A * hVal * yP * C1;
	    sData[22] += fs0A * y * hP * C2;
	    sData[23] += fs0A * y * hVal * C2;
	    sData[24] += fs0A * y * hP * C1;
	    sData[25] += fs0A * y * hVal * C1;
	    sData[26] += fs0A * hP * hP * C2 * C2;
	    sData[27] += fs0A * hVal * hP * C2 * C2;
	    sData[28] += fs0A * hP * hP * C1 * C2;
	    sData[29] += fs0A * hVal * hP * C1 * C2;
	    sData[30] += fs0A * h2 * C2 * C2;
	    sData[31] += fs0A * hVal * hP * C1 * C2;
	    sData[32] += fs0A * h2 * C1 * C2;
	    sData[33] += fs0A * hP * hP * C1 * C1;
	    sData[34] += fs0A * hVal * hP * C1 * C1;
	    sData[35] += fs0A * h2 * C1 * C1;
	}
    }

    return err;
}

int TaperedFiberSectionSmoothing3d::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    // create an id to send objects tag and numFibers, 
    // size 3 so no conflict with matData below if just 1 fiber
    static ID data(3);
    data(0) = this->getTag();
    data(1) = numFibers;
    int dbTag = this->getDbTag();
    res += theChannel.sendID(dbTag, commitTag, data);
    if (res < 0) {
        opserr <<
            "FiberSection3d::sendSelf - failed to send ID data\n";
        return res;
    }

    if (numFibers != 0) {

        // create an id containingg classTag and dbTag for each material & send it
        ID materialData(2 * numFibers);
        for (int i = 0; i < numFibers; i++) {
            UniaxialMaterial *theMat = theMaterials[i];
            materialData(2 * i) = theMat->getClassTag();
            int matDbTag = theMat->getDbTag();
            if (matDbTag == 0) {
                matDbTag = theChannel.getDbTag();
                if (matDbTag != 0)
                    theMat->setDbTag(matDbTag);
            }
            materialData(2 * i + 1) = matDbTag;
        }

        res += theChannel.sendID(dbTag, commitTag, materialData);
        if (res < 0) {
            opserr <<
                "FiberSection3d::sendSelf - failed to send material data\n";
            return res;
        }

        // send the fiber data, i.e. area and loc
        Vector fiberData(matData, 5 * numFibers);
        res += theChannel.sendVector(dbTag, commitTag, fiberData);
        if (res < 0) {
            opserr <<
                "FiberSection2d::sendSelf - failed to send material data\n";
            return res;
        }

        // now invoke send(0 on all the materials
        for (int j = 0; j < numFibers; j++)
            theMaterials[j]->sendSelf(commitTag, theChannel);
    }

    return res;
}

int TaperedFiberSectionSmoothing3d::recvSelf(int commitTag, Channel & theChannel,
                             FEM_ObjectBroker & theBroker)
{
    int res = 0;

    static ID data(3);

    int dbTag = this->getDbTag();
    res += theChannel.recvID(dbTag, commitTag, data);

    if (res < 0) {
        opserr <<
            "FiberSection2d::sendSelf - failed to recv ID data\n";
        return res;
    }

    this->setTag(data(0));

    // recv data about materials objects, classTag and dbTag
    if (data(1) != 0) {
        ID materialData(2 * data(1));
        res += theChannel.recvID(dbTag, commitTag, materialData);
        if (res < 0) {
            opserr <<
                "FiberSection2d::sendSelf - failed to send material data\n";
            return res;
        }

        // if current arrays not of correct size, release old and resize
        if (theMaterials == 0 || numFibers != data(1)) {
            // delete old stuff if outa date
            if (theMaterials != 0) {
                for (int i = 0; i < numFibers; i++)
                    delete theMaterials[i];
                delete[]theMaterials;
                if (matData != 0)
                    delete[]matData;
                matData = 0;
                theMaterials = 0;
            }

            // create memory to hold material pointers and fiber data
            numFibers = data(1);
            if (numFibers != 0) {

                theMaterials = new UniaxialMaterial *[numFibers];

                if (theMaterials == 0) {
                    opserr <<
                        "FiberSection3d::recvSelf -- failed to allocate Material pointers\n";
                    exit(-1);
                }

                for (int j = 0; j < numFibers; j++)
                    theMaterials[j] = 0;

                matData = new double[numFibers * 5];

                if (matData == 0) {
                    opserr <<
                        "FiberSection2d::recvSelf  -- failed to allocate double array for material data\n";
                    exit(-1);
                }
            }
        }

        Vector fiberData(matData, 5 * numFibers);
        res += theChannel.recvVector(dbTag, commitTag, fiberData);
        if (res < 0) {
            opserr <<
                "FiberSection2d::sendSelf - failed to send material data\n";
            return res;
        }

        int i;
        for (i = 0; i < numFibers; i++) {
            int classTag = materialData(2 * i);
            int dbTag = materialData(2 * i + 1);

            // if material pointed to is blank or not of corrcet type, 
            // release old and create a new one
            if (theMaterials[i] == 0)
                theMaterials[i] =
                    theBroker.getNewUniaxialMaterial(classTag);
            else if (theMaterials[i]->getClassTag() != classTag) {
                delete theMaterials[i];
                theMaterials[i] =
                    theBroker.getNewUniaxialMaterial(classTag);
            }

            if (theMaterials[i] == 0) {
                opserr <<
                    "FiberSection2d::recvSelf -- failed to allocate double array for material data\n";
                exit(-1);
            }

            theMaterials[i]->setDbTag(dbTag);
            res +=
                theMaterials[i]->recvSelf(commitTag, theChannel,
                                          theBroker);
        }

        double Qz = 0.0;
        double Qy = 0.0;
        double A = 0.0;
        double yLoc, zLoc, Area;

        // Recompute centroid
        for (i = 0; i < numFibers; i++) {
            yLoc   = matData[5 * i];
            zLoc   = matData[5 * i + 1];
            Area   = matData[5 * i + 2];

            A  += Area;
            Qz += yLoc * Area;
            Qy += zLoc * Area;
        }

        yBar = Qz / A;
        zBar = Qy / A;
    }

    return res;
}

void TaperedFiberSectionSmoothing3d::Print(OPS_Stream & s, int flag)
{
    if (flag == 2) {
        for (int i = 0; i < numFibers; i++) {
            s << -matData[5 * i] << " " << matData[5 * i +
                                                   1] << " " <<
                matData[5 * i + 2] << " ";
            s << theMaterials[i]->
                getStress() << " " << theMaterials[i]->
                getStrain() << endln;
        }
    }
    else {
        s << "\nTaperedFiberSectionSmoothing3d, tag: " << this->getTag() << endln;
        s << "\tSection code: " << code;
        s << "\tNumber of Fibers: " << numFibers << endln;
        s << "\tCentroid: (" << yBar << ", " << zBar << ')' << endln;

        if (flag == 1) {
            for (int i = 0; i < numFibers; i++) {
                s << "\nLocation (y, z) = (" << matData[6 *
                                                         i] << ", " <<
                    matData[5 * i + 1] << ")";
                s << "\nArea = " << matData[5 * i + 2] << endln;
                theMaterials[i]->Print(s, flag);
            }
        }
    }
}

Response *TaperedFiberSectionSmoothing3d::setResponse(const char **argv, int argc,
                                      OPS_Stream & output)
{

    const ID & type = this->getType();
    int typeSize = this->getOrder();

    Response *theResponse = 0;

    output.tag("SectionOutput");
    output.attr("secType", this->getClassType());
    output.attr("secTag", this->getTag());

    // deformations
    if (strcmp(argv[0], "deformations") == 0
        || strcmp(argv[0], "deformation") == 0) {
        for (int i = 0; i < typeSize; i++) {
            int code = type(i);
            switch (code) {
            case SECTION_RESPONSE_MZ:
                output.tag("ResponseType", "kappaZ");
                break;
            case SECTION_RESPONSE_P:
                output.tag("ResponseType", "eps");
                break;
            case SECTION_RESPONSE_VY:
                output.tag("ResponseType", "gammaY");
                break;
            case SECTION_RESPONSE_MY:
                output.tag("ResponseType", "kappaY");
                break;
            case SECTION_RESPONSE_VZ:
                output.tag("ResponseType", "gammaZ");
                break;
            case SECTION_RESPONSE_T:
                output.tag("ResponseType", "theta");
                break;
            default:
                output.tag("ResponseType", "Unknown");
            }
        }
        theResponse =
            new MaterialResponse(this, 1,
                                 this->getSectionDeformation());

        // forces
    }
    else if (strcmp(argv[0], "forces") == 0
             || strcmp(argv[0], "force") == 0) {
        for (int i = 0; i < typeSize; i++) {
            int code = type(i);
            switch (code) {
            case SECTION_RESPONSE_MZ:
                output.tag("ResponseType", "Mz");
                break;
            case SECTION_RESPONSE_P:
                output.tag("ResponseType", "P");
                break;
            case SECTION_RESPONSE_VY:
                output.tag("ResponseType", "Vy");
                break;
            case SECTION_RESPONSE_MY:
                output.tag("ResponseType", "My");
                break;
            case SECTION_RESPONSE_VZ:
                output.tag("ResponseType", "Vz");
                break;
            case SECTION_RESPONSE_T:
                output.tag("ResponseType", "T");
                break;
            default:
                output.tag("ResponseType", "Unknown");
            }
        }
        theResponse =
            new MaterialResponse(this, 2, this->getStressResultant());

        // force and deformation
    }
    else if (strcmp(argv[0], "forceAndDeformation") == 0) {
        for (int i = 0; i < typeSize; i++) {
            int code = type(i);
            switch (code) {
            case SECTION_RESPONSE_MZ:
                output.tag("ResponseType", "kappaZ");
                break;
            case SECTION_RESPONSE_P:
                output.tag("ResponseType", "eps");
                break;
            case SECTION_RESPONSE_VY:
                output.tag("ResponseType", "gammaY");
                break;
            case SECTION_RESPONSE_MY:
                output.tag("ResponseType", "kappaY");
                break;
            case SECTION_RESPONSE_VZ:
                output.tag("ResponseType", "gammaZ");
                break;
            case SECTION_RESPONSE_T:
                output.tag("ResponseType", "theta");
                break;
            default:
                output.tag("ResponseType", "Unknown");
            }
        }
        for (int j = 0; j < typeSize; j++) {
            int code = type(j);
            switch (code) {
            case SECTION_RESPONSE_MZ:
                output.tag("ResponseType", "Mz");
                break;
            case SECTION_RESPONSE_P:
                output.tag("ResponseType", "P");
                break;
            case SECTION_RESPONSE_VY:
                output.tag("ResponseType", "Vy");
                break;
            case SECTION_RESPONSE_MY:
                output.tag("ResponseType", "My");
                break;
            case SECTION_RESPONSE_VZ:
                output.tag("ResponseType", "Vz");
                break;
            case SECTION_RESPONSE_T:
                output.tag("ResponseType", "T");
                break;
            default:
                output.tag("ResponseType", "Unknown");
            }
        }

        theResponse =
            new MaterialResponse(this, 4,
                                 Vector(2 * this->getOrder()));

    }

    else {
        if (argc > 2 || strcmp(argv[0], "TaperedFiber") == 0) {

            int key = numFibers;
            int passarg = 2;


            if (argc <= 3) {    // fiber number was input directly

                key = atoi(argv[1]);

            }
            else if (argc > 4) {        // find fiber closest to coord. with mat tag
                int matTag = atoi(argv[3]);
                double yCoord = atof(argv[1]);
                double zCoord = atof(argv[2]);
                double closestDist;
                double ySearch, zSearch, dy, dz;
                double distance;
                int j;

                // Find first fiber with specified material tag
                for (j = 0; j < numFibers; j++) {
                    if (matTag == theMaterials[j]->getTag()) {
                        ySearch = matData[5 * j];
                        zSearch = matData[5 * j + 1];
                        dy = ySearch - yCoord;
                        dz = zSearch - zCoord;
                        closestDist = sqrt(dy * dy + dz * dz);
                        key = j;
                        break;
                    }
                }

                // Search the remaining fibers
                for (; j < numFibers; j++) {
                    if (matTag == theMaterials[j]->getTag()) {
                        ySearch = matData[5 * j];
                        zSearch = matData[5 * j + 1];
                        dy = ySearch - yCoord;
                        dz = zSearch - zCoord;
                        distance = sqrt(dy * dy + dz * dz);
                        if (distance < closestDist) {
                            closestDist = distance;
                            key = j;
                        }
                    }
                }
                passarg = 4;
            }

            else {              // fiber near-to coordinate specified
                double yCoord = atof(argv[1]);
                double zCoord = atof(argv[2]);
                double closestDist;
                double ySearch, zSearch, dy, dz;
                double distance;
                ySearch = matData[0];
                zSearch = matData[1];
                dy = ySearch - yCoord;
                dz = zSearch - zCoord;
                closestDist = sqrt(dy * dy + dz * dz);
                key = 0;
                for (int j = 1; j < numFibers; j++) {
                    ySearch = matData[5 * j];
                    zSearch = matData[5 * j + 1];
                    dy = ySearch - yCoord;
                    dz = zSearch - zCoord;
                    distance = sqrt(dy * dy + dz * dz);
                    if (distance < closestDist) {
                        closestDist = distance;
                        key = j;
                    }
                }
                passarg = 3;
            }

            if (key < numFibers && key >= 0) {
                output.tag("FiberOutput");
                output.attr("yLoc", -matData[5 * key]);
                output.attr("zLoc", matData[5 * key + 1]);
                output.attr("area", matData[5 * key + 2]);

                theResponse =
                    theMaterials[key]->setResponse(&argv[passarg],
                                                   argc - passarg,
                                                   output);

                output.endTag();
            }
        }
    }

    output.endTag();
    return theResponse;
}

int TaperedFiberSectionSmoothing3d::getResponse(int responseID,
                                Information & sectInfo)
{
    // Just call the base class method ... don't need to define
    // this function, but keeping it here just for clarity
    return SectionForceDeformation::getResponse(responseID, sectInfo);
}

int TaperedFiberSectionSmoothing3d::setParameter(const char **argv, int argc,
                                 Parameter & param)
{
    if (argc < 3)
        return -1;


    int result = 0;

    // A material parameter
    if (strstr(argv[0], "material") != 0) {

        // Get the tag of the material
        int paramMatTag = atoi(argv[1]);

        // Loop over fibers to find the right material(s)
        int ok = 0;
        for (int i = 0; i < numFibers; i++)
            if (paramMatTag == theMaterials[i]->getTag()) {
                ok = theMaterials[i]->setParameter(&argv[2], argc - 2,
                                                   param);
                if (ok != -1)
                    result = ok;
            }

        return result;
    }

    int ok = 0;

    // loop over every material
    for (int i = 0; i < numFibers; i++) {
        ok = theMaterials[i]->setParameter(argv, argc, param);
        if (ok != -1)
            result = ok;
    }

    return result;
}

const Vector &
    TaperedFiberSectionSmoothing3d::getSectionDeformationSensitivity(int gradIndex)
{
    static Vector dummy(3);
    dummy.Zero();
    if (SHVs != 0) {
        dummy(0) = (*SHVs) (0, gradIndex);
        dummy(1) = (*SHVs) (1, gradIndex);
        dummy(2) = (*SHVs) (2, gradIndex);
    }
    return dummy;
}


const Vector &
    TaperedFiberSectionSmoothing3d::getStressResultantSensitivity(int gradIndex,
                                                  bool conditional)
{

    static Vector ds(3);

    ds.Zero();

    double stressGradient;
    int loc = 0;


    for (int i = 0; i < numFibers; i++) {
        UniaxialMaterial *theMat = theMaterials[i];
        double y = matData[loc++] - yBar;
        double z = matData[loc++] - zBar;
        double A = matData[loc++];
        stressGradient =
            theMaterials[i]->getStressSensitivity(gradIndex,
                                                  conditional);
        stressGradient *= A;
        ds(0) += stressGradient;
        ds(1) += stressGradient * y;
        ds(2) += stressGradient * z;

    }

    return ds;
}

const Matrix &
    TaperedFiberSectionSmoothing3d::getSectionTangentSensitivity(int gradIndex)
{
    static Matrix something(2, 2);

    something.Zero();

    return something;
}

int TaperedFiberSectionSmoothing3d::commitSensitivity(const Vector & defSens,
                                      int gradIndex, int numGrads)
{

    // here add SHVs to store the strain sensitivity.

    if (SHVs == 0) {
        SHVs = new Matrix(3, numGrads);
    }

    (*SHVs) (0, gradIndex) = defSens(0);
    (*SHVs) (1, gradIndex) = defSens(1);
    (*SHVs) (2, gradIndex) = defSens(2);

    int loc = 0;

    double d0 = defSens(0);
    double d1 = defSens(1);
    double d2 = defSens(2);

    for (int i = 0; i < numFibers; i++) {
        UniaxialMaterial *theMat = theMaterials[i];
        double y = matData[loc++] - yBar;
        double z = matData[loc++] - zBar;
        loc++;                  // skip A data.

        double strainSens = d0 + y * d1 + z * d2;



        theMat->commitSensitivity(strainSens, gradIndex, numGrads);
    }

    return 0;
}

// AddingSensitivity:END ///////////////////////////////////
