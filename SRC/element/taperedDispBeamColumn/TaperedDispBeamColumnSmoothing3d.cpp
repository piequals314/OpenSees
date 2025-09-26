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

// $Revision: 1.25 $
// $Date: 2008/11/04 21:32:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/TaperedDispBeamColumnSmoothing3d.cpp,v $

// Description: This file contains the class definition for TaperedDispBeamColumnSmoothing3d which 
// include warping degrees of freedom and local buckling.
// Modified by Brighton Laiman from University of California, San Diego (include warping degrees of freedom 
// AND Local Buckling). Refer to Formulation and Implementation of Three-dimensional Beam-Column 
// Analyses with Warping Effects 

#define DEBUG
//#define STIFFMATRIX
//#define N*MATRIX

#include <TaperedDispBeamColumnSmoothing3d.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <BeamIntegration.h>
#include <Parameter.h>
#include <math.h>
#include <elementAPI.h>
#include <string>

using std::string;
using namespace std;

Matrix TaperedDispBeamColumnSmoothing3d::K(22, 22);
Vector TaperedDispBeamColumnSmoothing3d::P(22);
double TaperedDispBeamColumnSmoothing3d::workArea[800];

void* OPS_TaperedDispBeamColumnSmoothing3d()
{
    int dampingTag = 0;
    Damping* theDamping = 0;
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "insufficient arguments:eleTag,iNode,jNode,transfTag,integrationTag <-mass mass> <-cmass>\n";
        return 0;
    }

    // inputs: 
    int iData[5];
    int numData = 5;
    if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
        opserr << "WARNING: invalid integer inputs\n";
        return 0;
    }

    // options
    double mass = 0.0;
    int cmass = 0;
    numData = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* type = OPS_GetString();
        if (strcmp(type, "-cMass") == 0) {
            cmass = 1;
        }
        else if (strcmp(type, "-mass") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numData, &mass) < 0) {
                    opserr << "WARNING: invalid mass\n";
                    return 0;
                }
            }
 /*       }
        else if (strcmp(type, "-damp") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
                theDamping = OPS_getDamping(dampingTag);
                if (theDamping == 0) {
                    opserr << "damping not found\n";
                    return 0;
                }
            }*/
        }
    }

    // check transf
    CrdTransf* theTransf = OPS_getCrdTransf(iData[3]);
    if (theTransf == 0) {
        opserr << "coord transfomration not found\n";
        return 0;
    }

    // check beam integrataion
    BeamIntegrationRule* theRule = OPS_getBeamIntegrationRule(iData[4]);
    if (theRule == 0) {
        opserr << "beam integration not found\n";
        return 0;
    }
    BeamIntegration* bi = theRule->getBeamIntegration();
    if (bi == 0) {
        opserr << "beam integration is null\n";
        return 0;
    }

    // check sections
    const ID& secTags = theRule->getSectionTags();
    SectionForceDeformation** sections = new SectionForceDeformation * [secTags.Size()];
    for (int i = 0; i < secTags.Size(); i++) {
        sections[i] = OPS_getSectionForceDeformation(secTags(i));
        if (sections[i] == 0) {
            opserr << "section " << secTags(i) << "not found\n";
            delete[] sections;
            return 0;
        }
    }

    Element* theEle = new TaperedDispBeamColumnSmoothing3d(iData[0], iData[1], iData[2], secTags.Size(), sections,
        *bi, *theTransf);
    delete[] sections;
    return theEle;
}

void* OPS_TaperedDispBeamColumnSmoothing3d(const ID& info) {
    // data needed
    int iData[5];
    double mass = 0.0;
    int cmass = 0;
    int numData;
    int dampingTag = 0;
    Damping* theDamping = 0;
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if (ndm != 3 || ndf != 6) {
        opserr << "ndm must be 3 and ndf must be 6\n";
        return 0;
    }

    // 1. regular elements
    if (info.Size() == 0) {
        numData = 3;
        if (OPS_GetNumRemainingInputArgs() < numData) {
            opserr << "insufficient "
                "arguments:eleTag,iNode,jNode\n";
            return 0;
        }
        if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
            opserr << "WARNING invalid int inputs\n";
            return 0;
        }
    }

    // 2. regular elements or save data
    if (info.Size() == 0 || info(0) == 1) {
        numData = 2;
        if (OPS_GetNumRemainingInputArgs() < numData) {
            opserr << "insufficient "
                "arguments:transfTag,integrationTag\n";
            return 0;
        }
        if (OPS_GetIntInput(&numData, &iData[3]) < 0) {
            opserr << "WARNING invalid int inputs\n";
            return 0;
        }

        numData = 1;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char* type = OPS_GetString();
            if (strcmp(type, "-cMass") == 0) {
                cmass = 1;
            }
            else if (strcmp(type, "-mass") == 0) {
                if (OPS_GetNumRemainingInputArgs() > 0) {
                    if (OPS_GetDoubleInput(&numData, &mass) < 0) {
                        opserr << "WARNING: invalid mass\n";
                        return 0;
                    }
                }
             }
         /*   else if (strcmp(type, "-damp") == 0) {

                if (OPS_GetNumRemainingInputArgs() > 0) {
                    if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
                    theDamping = OPS_getDamping(dampingTag);
                    if (theDamping == 0) {
                        opserr << "damping not found\n";
                        return 0;
                    }
                }
            }*/
        }
    }

    // 3: save data
    static std::map<int, Vector> meshdata;
    if (info.Size() > 0 && info(0) == 1) {
        if (info.Size() < 2) {
            opserr << "WARNING: need info -- inmesh, meshtag\n";
            return 0;
        }

        // save the data for a mesh
        Vector& mdata = meshdata[info(1)];
        mdata.resize(4);
        mdata(0) = iData[3];
        mdata(1) = iData[4];
        mdata(2) = mass;
        mdata(3) = cmass;
        return &meshdata;
    }

    // 4: load data
    if (info.Size() > 0 && info(0) == 2) {
        if (info.Size() < 5) {
            opserr << "WARNING: need info -- inmesh, meshtag, "
                "eleTag, nd1, nd2\n";
            return 0;
        }

        // get the data for a mesh
        Vector& mdata = meshdata[info(1)];
        if (mdata.Size() < 4) return 0;

        iData[0] = info(2);
        iData[1] = info(3);
        iData[2] = info(4);
        iData[3] = mdata(0);
        iData[4] = mdata(1);
        mass = mdata(2);
        cmass = mdata(3);
    }

    // 5: create element
    CrdTransf* theTransf = OPS_getCrdTransf(iData[3]);
    if (theTransf == 0) {
        opserr << "coord transfomration not found\n";
        return 0;
    }

    // check beam integrataion
    BeamIntegrationRule* theRule =
        OPS_getBeamIntegrationRule(iData[4]);
    if (theRule == 0) {
        opserr << "beam integration not found\n";
        return 0;
    }
    BeamIntegration* bi = theRule->getBeamIntegration();
    if (bi == 0) {
        opserr << "beam integration is null\n";
        return 0;
    }

    // check sections
    const ID& secTags = theRule->getSectionTags();
    SectionForceDeformation** sections =
        new SectionForceDeformation * [secTags.Size()];
    for (int i = 0; i < secTags.Size(); i++) {
        sections[i] = OPS_getSectionForceDeformation(secTags(i));
        if (sections[i] == 0) {
            opserr << "section " << secTags(i) << "not found\n";
            delete[] sections;
            return 0;
        }
    }

    Element* theEle = new TaperedDispBeamColumnSmoothing3d(
        iData[0], iData[1], iData[2], secTags.Size(), sections, *bi,
        *theTransf);
    delete[] sections;
    return theEle;
}



TaperedDispBeamColumnSmoothing3d::TaperedDispBeamColumnSmoothing3d(int tag, int nd1,
                                                 int nd2, int numSec,
                                                 SectionForceDeformation
                                                 ** s,
                                                 BeamIntegration & bi,
                                                 CrdTransf &
                                                 coordTransf,
                                                 double r)
:Element(tag, ELE_TAG_TaperedDispBeamColumnSmoothing3d), numSections(numSec),
theSections(0), crdTransf(0), beamInt(0), connectedExternalNodes(2),
Q(22), q(17), rho(r), parameterID(0)
{
    // Allocate arrays of pointers to SectionForceDeformations
    theSections = new SectionForceDeformation *[numSections];

    if (theSections == 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::TaperedDispBeamColumnSmoothing3d - failed to allocate section model pointer\n";
        exit(-1);
    }

    for (int i = 0; i < numSections; i++) {

        // Get copies of the material model for each integration point
        theSections[i] = s[i]->getCopy();

        // Check allocation
        if (theSections[i] == 0) {
            opserr <<
                "TaperedDispBeamColumnSmoothing3d::TaperedDispBeamColumnSmoothing3d -- failed to get a copy of section model\n";
            exit(-1);
        }
    }

    beamInt = bi.getCopy();

    if (beamInt == 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::TaperedDispBeamColumnSmoothing3d - failed to copy beam integration\n";
        exit(-1);
    }

    crdTransf = coordTransf.getCopy3d();

    if (crdTransf == 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::TaperedDispBeamColumnSmoothing3d - failed to copy coordinate transformation\n";
        exit(-1);
    }

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;


    theNodes[0] = 0;
    theNodes[1] = 0;

    for (int i = 0; i < 17; i++)
	q0[i] = 0.0;

    for (int i = 0; i < 5; i++)
	p0[i] = 0.0;

}

TaperedDispBeamColumnSmoothing3d::TaperedDispBeamColumnSmoothing3d()
:  
Element(0, ELE_TAG_TaperedDispBeamColumnSmoothing3d),
numSections(0), theSections(0), crdTransf(0), beamInt(0),
connectedExternalNodes(2), Q(22), q(17), rho(0.0), parameterID(0)
{

    for (int i = 0; i < 17; i++)
	q0[i] = 0.0;

    for (int i = 0; i < 5; i++)
	p0[i] = 0.0;

    theNodes[0] = 0;
    theNodes[1] = 0;
}

TaperedDispBeamColumnSmoothing3d::~TaperedDispBeamColumnSmoothing3d()
{
    for (int i = 0; i < numSections; i++) {
        if (theSections[i])
            delete theSections[i];
    }

    // Delete the array of pointers to SectionForceDeformation pointer arrays
    if (theSections)
        delete[]theSections;

    if (crdTransf)
        delete crdTransf;

    if (beamInt != 0)
        delete beamInt;
}

int TaperedDispBeamColumnSmoothing3d::getNumExternalNodes() const 
{

    return 2;
}

const ID & TaperedDispBeamColumnSmoothing3d::getExternalNodes()
{
    
    return connectedExternalNodes;
}

Node **TaperedDispBeamColumnSmoothing3d::getNodePtrs()
{

    return theNodes;
}

int TaperedDispBeamColumnSmoothing3d::getNumDOF()
{

    return 22;
}

void TaperedDispBeamColumnSmoothing3d::setDomain(Domain * theDomain)
{
    // Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        theNodes[0] = 0;
        theNodes[1] = 0;
        return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);

    if (theNodes[0] == 0 || theNodes[1] == 0) {
        opserr <<
            "FATAL ERROR TaperedDispBeamColumnSmoothing3d (tag: %d), node not found in domain",
            this->getTag();

        return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

    if (crdTransf->initialize(theNodes[0], theNodes[1])) {
	opserr << "\nInitialization of TaperedDispBeamColumnSmoothing3d failed\n";
	return;
    }

    double L = crdTransf->getInitialLength();

    if (L == 0.0) {
	opserr << "\nINITIAL LENGTH OF ELEMENT IS INVALID!\n";
	return;
    }

    this->DomainComponent::setDomain(theDomain);

    this->update();
}

int TaperedDispBeamColumnSmoothing3d::commitState()
{

    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::commitState () - failed in base class";
    }

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
        retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    return retVal;
}

int TaperedDispBeamColumnSmoothing3d::revertToLastCommit()
{
    #ifdef DEBUG
    fprintf(stdout, "TaperedDispBeamColumnSmoothing3d::revertToLastCommit\n");
    #endif

    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
        retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    return retVal;
}

int TaperedDispBeamColumnSmoothing3d::revertToStart()
{
    #ifdef DEBUG
    fprintf(stdout, "TaperedDispBeamColumnSmoothing3d::revertToStart\n");
    #endif

    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
        retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();

    return retVal;
}

int TaperedDispBeamColumnSmoothing3d::update(void)
{
    int err = 0;

    // Update the transformation
    crdTransf->update();

    // Get basic deformations
    // this is collecting qn (Chang) | n (Syd), which is the natural
    // displacements (no rigid body)
    const Vector & v = crdTransf->getBasicTrialDisp();

    double L = crdTransf->getInitialLength();
    double oneOverL       = 1.0 / L;
    double oneOverLsquare = oneOverL / L;

    //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

        // this is about to take the qn vector and create v
        // this is the same e as used in fibersection3d
        Vector e(workArea, 18);  

	double xi6  = 6.0 * xi[i];
	double xi12 = 12.0 * xi[i];
	double xi1  = xi[i];
	double x3   = xi[i] * xi[i] * xi[i];
	double x2   = xi[i] * xi[i];    

	// N means shape function, dN means 1st derivative, ddN 2nd 
	// N* - * is the order of the shape function 3-cubic 1-linear
	// N*x -  x is the dof of the shape function
	// xi is bringing in x/L already
	double N31   = 1 - 3 * x2 + 2 * x3;
	double N32   = xi[i] * L * (1 - xi[i]) * (1 - xi[i]);
	double N33   = 3 * x2 - 2 * x3;
	double N34   = -xi[i] * L * (xi[i] - x2);
	double dN31  = 6 * oneOverL * (x2 - xi[i]); 
	double dN32  = 1 + 3 * x2 - 4 * xi[i];
	double dN33  = 6 * oneOverL * (xi[i] - x2);
	double dN34  = 3 * x2 - 2 * xi[i];
	double ddN31 = oneOverLsquare * (xi12 - 6);
	double ddN32 = oneOverL * (xi6 - 4);
	double ddN33 = oneOverLsquare * (6 - xi12);
	double ddN34 = oneOverL * (xi6 - 2);
	double dN12  = oneOverL;

        // here total strain (not incremental) is used
        // e is v in Sydney's code
        // NOTE: This probably means that FiberSection3d needs reordering
	// The v referred to is the nodal displacements in the correct
	// frame. It has this order [rx1 rz1 ry1 ddrx1 alpha1 dalpha1 beta1 dbeta1 rx2 rz2 ry2 ddrx2 alpha2 dalpha2 beta2 dbeta2 e]
	// This is the new e:
	
	e(0)  = dN12 * v(16); 				      		     // u prime
	e(1)  = v(0); 						       	     // rx1
	e(2)  = v(1); 					               	     // rz1
	e(3)  = v(2);  						       	     // ry1
	e(4)  = v(8); 						       	     // rx2
	e(5)  = v(9); 						       	     // rz2
	e(6)  = v(10); 						    	     // ry2
	e(7)  =  ddN32 * v(1) + ddN34 * v(9); 			             // v dprime
	e(8)  = -ddN32 * v(2) - ddN34 * v(10); 			             // w dprime
	e(9)  =   N31 * v(0) +   N32 * v(3) +   N33 * v(8)  +   N34 * v(11); // Theta
	e(10) =  dN31 * v(0) +  dN32 * v(3) +  dN33 * v(8)  +  dN34 * v(11); // Theta prime
	e(11) = ddN31 * v(0) + ddN32 * v(3) + ddN33 * v(8)  + ddN34 * v(11); // Theta dprime
	e(12) =   N31 * v(4) +   N32 * v(5) +   N33 * v(12) +   N34 * v(13); // alpha
	e(13) =  dN31 * v(4) +  dN32 * v(5) +  dN33 * v(12) +  dN34 * v(13); // alpha'
	e(14) = ddN31 * v(4) + ddN32 * v(5) + ddN33 * v(12) + ddN34 * v(13); // alpha''
	e(15) =   N31 * v(6) +   N32 * v(7) +   N33 * v(14) +   N34 * v(15); // beta
	e(16) =  dN31 * v(6) +  dN32 * v(7) +  dN33 * v(14) +  dN34 * v(15); // beta'
	e(17) = ddN31 * v(6) + ddN32 * v(7) + ddN33 * v(14) + ddN34 * v(15); // beta''

	#ifdef DEBUG
	/*
	fprintf(stdout, "Printing the total strain of the BASIC system for Section %2d:\n\n", i+1);
	for (int a = 0; a < 18; ++a)
	    fprintf(stdout, "e(%2d): %2.1f\n", a+1, e(a));
	fprintf(stdout, "\n");
	*/
	#endif
	
        // Set the section deformations
	err += theSections[i]->setTrialSectionDeformation(e);
    }

    if (err != 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::update() - failed setTrialSectionDeformations()\n";
        return err;
    }

    return 0;
}

const Matrix & TaperedDispBeamColumnSmoothing3d::getTangentStiff()
{
    // kb is the element stiffness matrix, natural
    // N1 is Q in Chang and Ndeld1 in Sydney
    // N2sec is B in Chang and Ndeld2 in Sydney
    // N3 is made up here for calc purposes and is N1_transpose * ks * N1
    // N4sec is to transform the Gmax matrix into the natural system DOF
    // kbPart1 is the elastic stiffness matrix part of kb
    // kbPart2 is the geometric stiffness matrix
    static Matrix kb(17,17);
    static Matrix N1(12,18);
    static Matrix N2(18,17);
    static Matrix N3(18,18);
    static Matrix N4(14,17);
    static Matrix kbPart1(17,17);
    static Matrix Gmax(14,14); 
    static Matrix kbPart2(17,17);

    const Vector & v = crdTransf->getBasicTrialDisp();

    // Zero for integral
    kb.Zero();
    q.Zero();

    double L = crdTransf->getInitialLength();
    double oneOverL       = 1.0 / L;
    double oneOverLsquare = oneOverL / L;

    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);
    double wt[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wt);

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

        N1.Zero();
        N2.Zero();
        N3.Zero();
        N4.Zero();
        kbPart1.Zero();
        Gmax.Zero();
        kbPart2.Zero();

	double xi6  = 6.0 * xi[i];
	double xi12 = 12.0 * xi[i];
	double xi1  = xi[i];
	double x3   = xi[i] * xi[i] * xi[i];
	double x2   = xi[i] * xi[i];    

	// N means shape function, dN means 1st derivative, ddN 2nd 
	// N* - * is the order of the shape function 3-cubic 1-linear
	// N*x -  x is the dof of the shape function
	// xi is bringing in x/L already
	double N31   = 1 - 3 * x2 + 2 * x3;
	double N32   = xi[i] * L * (1 - xi[i]) * (1 - xi[i]);
	double N33   = 3 * x2 - 2 * x3;
	double N34   = -xi[i] * L * (xi[i] - x2);
	double dN31  = 6 * oneOverL * (x2 - xi[i]); 
	double dN32  = 1 + 3 * x2 - 4 * xi[i];
	double dN33  = 6 * oneOverL * (xi[i] - x2);
	double dN34  = 3 * x2 - 2 * xi[i];
	double ddN31 = oneOverLsquare * (xi12 - 6);
	double ddN32 = oneOverL * (xi6 - 4);
	double ddN33 = oneOverLsquare * (6 - xi12);
	double ddN34 = oneOverL * (xi6 - 2);
	double dN12  = oneOverL;

	N1(0,0)   = 1.0;
        N1(0,2)   = (4 * v( 1) - v( 9)) / 30;
        N1(0,3)   = (4 * v( 2) - v(10)) / 30;
        N1(0,5)   = (4 * v( 9) - v( 1)) / 30;
        N1(0,6)   = (4 * v(10) - v( 2)) / 30;
        N1(1,7)   = 1.0;
        N1(1,8)   =    N31 * v(0) +   N32 * v(3) + N33 * v(8) + N34 * v(11);  // Theta
        N1(1,9)   = -ddN32 * v(2) - ddN34 * v(10); 			      // w dprime
        N1(2,7)   = N1(1,8);
        N1(2,8)   = -1.0;
        N1(2,9)   = ddN32 * v(1) + ddN34 * v(9); 		       	      // v dprime
        N1(3,10)  =  dN31 * v(0) +  dN32 * v(3) + dN33 * v(8) + dN34 * v(11); // Theta prime
        N1(4,11)  = 1.0;
        N1(5,10)  = 1.0; 
	N1(6,12)  = 1.0; 
	N1(7,13)  = 1.0; 
	N1(8,14)  = 1.0; 
	N1(9,15)  = 1.0; 
	N1(10,16) = 1.0; 
	N1(11,17) = 1.0; 

	N2(0, 16)  =   dN12;
	N2(1, 0)   =      1;
	N2(2, 1)   =      1;
	N2(3, 2)   =      1;
	N2(4, 8)   =      1;
	N2(5, 9)   =      1;
	N2(6, 10)  =      1;
	N2(7, 1)   =  ddN32;  // rz1->v dprime
	N2(7, 9)   =  ddN34;  // rz2->v dprime
	N2(8, 2)   = -ddN32;  // ry1->w dprime
	N2(8, 10)  = -ddN34;  // ry2->w dprime
	N2(9, 0)   =    N31;  // rx1->theta
	N2(9, 3)   =    N32;  // rx1 prime->theta
	N2(9, 8)   =    N33;  // rx2->theta
	N2(9, 11)  =    N34;  // rx2 prime->theta
	N2(10, 0)  =   dN31;  // rx1->theta prime
	N2(10, 3)  =   dN32;  // rx1 prime->theta prime
	N2(10, 8)  =   dN33;  // rx2->theta prime
	N2(10, 11) =   dN34;  // rx2 prime->theta prime
	N2(11, 0)  =  ddN31;  // rx1->theta dprime
	N2(11, 3)  =  ddN32;  // rx1 prime->theta dprime
	N2(11, 8)  =  ddN33;  // rx2->theta dprime
	N2(11, 11) =  ddN34;  // rx2 prime->theta dprime 
	N2(12, 4)  =    N31;
	N2(12, 5)  =    N32;
	N2(12, 12) =    N33;
	N2(12, 13) =    N34;
	N2(13, 4)  =   dN31;
	N2(13, 5)  =   dN32;
	N2(13, 12) =   dN33;
	N2(13, 13) =   dN34;
	N2(14, 4)  =  ddN31;
	N2(14, 5)  =  ddN32;
	N2(14, 12) =  ddN33;
	N2(14, 13) =  ddN34;
	N2(15, 6)  =    N31;
	N2(15, 7)  =    N32;
	N2(15, 14) =    N33;
	N2(15, 15) =    N34;
	N2(16, 6)  =   dN31;
	N2(16, 7)  =   dN32;
	N2(16, 14) =   dN33;
	N2(16, 15) =   dN34;
	N2(17, 6)  =  ddN31;
	N2(17, 7)  =  ddN32;
	N2(17, 14) =  ddN33;
	N2(17, 15) =  ddN34;

	N4(0,1)   =   dN32;
	N4(0,9)   =   dN34;
	N4(1,1)   =      1;
	N4(2,2)   =      1;
	N4(3,2)   =  -dN32;
	N4(3,10)  =  -dN34;
	N4(4,9)   =      1;
	N4(4,10)  =      1;
	N4(6,1)   =  ddN32; // rz1->v dprime
	N4(6,9)   =  ddN34; // rz2->v dprime
	N4(7,2)   = -ddN32; // ry1->w dprime
	N4(7,10)  = -ddN34; // ry2->w dprime
	N4(8,0)   =    N31; // rx1->theta
	N4(8,3)   =    N32; // rx1 prime->theta
	N4(8,8)   =    N33; // rx2->theta
	N4(8,11)  =    N34; // rx2 prime->theta
	N4(9,0)   =   dN31; // rx1->theta prime
	N4(9,3)   =   dN32; // rx1 prime->theta prime
	N4(9,8)   =   dN33; // rx2->theta prime
	N4(9,11)  =   dN34; // rx2 prime->theta prime
	N4(10,4)  =    N31;
	N4(10,5)  =    N32;
	N4(10,12) =    N33;
	N4(10,13) =    N34;
	N4(11,4)  =   dN31;
	N4(11,5)  =   dN32;
	N4(11,12) =   dN33;
	N4(11,13) =   dN34;
	N4(12,6)  =    N31;
	N4(12,7)  =    N32;
	N4(12,14) =    N33;
	N4(12,15) =    N34;
	N4(13,6)  =   dN31;
	N4(13,7)  =   dN32;
	N4(13,14) =   dN33;
	N4(13,15) =   dN34;

	#ifdef N*MATRIX
	fprintf(stdout, "TaperedDispBeamColumnSmoothing3d::getTangentStiff:\n");
	fprintf(stdout, "Printing N1 for Section %2d:\n\n", i);
	for (int a = 0; a < 12; ++a) {
	    for (int b = 0; b < 18; ++b)
		fprintf(stdout, "%4.1f ", N1(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "Printing N2 for Section %2d:\n\n", i);
	for (int a = 0; a < 18; ++a) {
	    for (int b = 0; b < 17; ++b)
		fprintf(stdout, "%4.1f ", N2(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "Printing N4 for Section %2d:\n\n", i);
	for (int a = 0; a < 14; ++a) {
	    for (int b = 0; b < 17; ++b)
		fprintf(stdout, "%4.1f ", N4(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");
	#endif

        // Get the section tangent stiffness and stress resultant
        const Matrix & ks = theSections[i]->getSectionTangent();

        //calculate kb, refer to Alemdar
        N3.addMatrixTripleProduct(0.0, N1, ks, 1.0);

        kbPart1.addMatrixTripleProduct(0.0, N2, N3, 1.0);

        const Vector &s = theSections[i]->getStressResultant();

	// s needs to be in the following order
	// [P Mz My W B Tsv]
	 Gmax(1,1)   = (4/30) * s(0);
	 Gmax(2,2)   = Gmax(1,1);
	 Gmax(4,4)   = Gmax(1,1);
	 Gmax(5,5)   = Gmax(1,1);
	 Gmax(1,4)   = -s(0) / 30;
	 Gmax(4,1)   = Gmax(1,4);
	 Gmax(2,5)   = Gmax(1,4);
	 Gmax(5,2)   = Gmax(1,4);
	 Gmax(6,8)   = s(2);
	 Gmax(8,6)   = Gmax(6,8);
	 Gmax(7,8)   = s(1);
	 Gmax(8,7)   = Gmax(7,8);
	 Gmax(10,10) = s(3); 
	 // Couple Portion of Gmax
	 Gmax(11,0) = s(12);
	 Gmax(0,11) = Gmax(11,0);
	 Gmax(13,0) = s(13);
	 Gmax(0,13) = Gmax(13,0);
	 Gmax(10,3) = s(14);
	 Gmax(3,10) = Gmax(10,3);
	 Gmax(11,3) = s(15);
	 Gmax(3,11) = Gmax(11,3);
	 Gmax(12,3) = s(16);
	 Gmax(3,12) = Gmax(12,3);
	 Gmax(13,3) = s(17);
	 Gmax(3,13) = Gmax(13,3);
	 Gmax(10,8) = s(18);
	 Gmax(8,10) = Gmax(10,8);
	 Gmax(11,8) = s(19);
	 Gmax(8,11) = Gmax(11,8);
	 Gmax(12,8) = s(20);
	 Gmax(8,12) = Gmax(12,8);
	 Gmax(13,8) = s(21);
	 Gmax(8,13) = Gmax(13,8);
	 Gmax(10,9) = s(22);
	 Gmax(9,10) = Gmax(10,9);
	 Gmax(11,9) = s(23);
	 Gmax(9,11) = Gmax(11,9);
	 Gmax(12,9) = s(24);
	 Gmax(9,12) = Gmax(12,9);
	 Gmax(13,9) = s(25);
	 Gmax(9,13) = Gmax(13,9);
	 Gmax(10,10) = s(26);
	 Gmax(11,10) = s(27);
	 Gmax(10,11) = Gmax(11,10);
	 Gmax(12,10) = s(28);
	 Gmax(10,12) = Gmax(12,10);
	 Gmax(13,9)  = s(29);
	 Gmax(9,13)  = Gmax(13,9);
	 Gmax(11,11) = s(30);
	 Gmax(12,11) = s(31);
	 Gmax(11,12) = Gmax(12,11);
	 Gmax(13,11) = s(32);
	 Gmax(11,13) = Gmax(13,11);
	 Gmax(12,12) = s(33);
	 Gmax(13,12) = s(34);
	 Gmax(12,13) = Gmax(13,12);
	 Gmax(13,13) = s(35);

	 kbPart2.addMatrixTripleProduct(0.0, N4, Gmax, 1.0); 

        // Perform numerical integration
        double wti = wt[i];

        for (int j = 0; j < 17; j++) {
            for (int k = 0; k < 17; k++) {
                kb(j, k) += (kbPart1(j, k) + kbPart2(j, k)) * L * wti;
            }
        }

	// Establish the section forces before transforming to nodal forces
	static Vector sForces(12);
	sForces.Zero();
	for (int j = 0; j < 12; j++)
	    sForces(j) = s(j);

	#ifdef DEBUG
	fprintf(stdout, "TaperedDispBeamColumnSmoothing3d::getTangentStiff:\n");
	fprintf(stdout, "Section Forces for Section: %2d\n", i+1);
	for(int c = 0; c < 12; ++c) 
	    fprintf(stdout, "The value of s[%d] = %10.2e\n", c, sForces(c));
	fprintf(stdout, "\n");
	#endif

	#ifdef STIFFMATRIX
	fprintf(stdout, "Printing Section Stiffness Matrix for Section %2d:\n", i+1);
	for (int a = 0; a < 12; ++a) {
	    for (int b = 0; b < 12; ++b)
		fprintf(stdout, "%10.2e ", ks(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "Printing Material Stiffness Matrix for Section %2d:\n", i+1);
	for (int a = 0; a < 17; ++a) {
	    for (int b = 0; b < 17; ++b)
		fprintf(stdout, "%10.2e ", kbPart1(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "Printing Geometric Stiffness Matrix for Section %2d:\n", i+1);
	for (int a = 0; a < 17; ++a) {
	    for (int b = 0; b < 17; ++b)
		fprintf(stdout, "%10.2e ", kbPart2(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "Printing BASIC Stiffness Matrix with %2d sections:\n", i+1);
	for (int a = 0; a < 17; ++a) {
	    for (int b = 0; b < 17; ++b)
		fprintf(stdout, "%10.2e ", kb(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");
	#endif

        static Vector qProduct1(18);
        static Vector qProduct2(17);
        qProduct1.Zero();
        qProduct2.Zero();
        qProduct1.addMatrixTransposeVector(0.0, N1, sForces, 1.0);
        qProduct2.addMatrixTransposeVector(0.0, N2, qProduct1, 1.0);

        for (int j = 0; j < 17; j++) {
            q(j) += qProduct2(j) * L * wti;
        }
    }

    // Transform to global stiffness
    K = crdTransf->getGlobalStiffMatrix(kb, q);

    //opserr << "TaperedDispBeamColumnSmoothing3d (tag: %d): Element Stiffness Matrix", this->getTag() << endln;
    /*
    fprintf(stdout, "TaperedDispBeamColumnSmoothing3d (tag: %d): Element Stiffness Matrix\n", this->getTag());
    for (int a = 0; a < 22; ++a) {
	for (int b = 0; b < 22; ++b)
	    fprintf(stdout, "%10.2e ", K(a,b));
	fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
    */

    return K;
}

const Matrix &TaperedDispBeamColumnSmoothing3d::getInitialBasicStiff()
{
    #ifdef DEBUG
    fprintf(stdout, "TaperedDispBeamColumnSmoothing3d::getInitialBasicStiff\n");
    #endif

    static Matrix kb(17,17);
    static Matrix N1(12,18);
    static Matrix N2(18,17);
    static Matrix N3(18,18);
    static Matrix N4(14,17);
    static Matrix kbPart1(17,17);
    static Matrix Gmax(14,14); 
    static Matrix kbPart2(17,17);

    // Zero for integral
    kb.Zero();
    const Vector & v = crdTransf->getBasicTrialDisp();

    double L = crdTransf->getInitialLength();
    double oneOverL       = 1.0 / L;
    double oneOverLsquare = oneOverL / L;

    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);
    double wt[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wt);

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

        N1.Zero();
        N2.Zero();
        N3.Zero();
        N4.Zero();
        kbPart1.Zero();
        Gmax.Zero();
        kbPart2.Zero();

	double xi6  = 6.0 * xi[i];
	double xi12 = 12.0 * xi[i];
	double xi1  = xi[i];
	double x3   = xi[i] * xi[i] * xi[i];
	double x2   = xi[i] * xi[i];    

	// N means shape function, dN means 1st derivative, ddN 2nd 
	// N* - * is the order of the shape function 3-cubic 1-linear
	// N*x -  x is the dof of the shape function
	// xi is bringing in x/L already
	double N31   = 1 - 3 * x2 + 2 * x3;
	double N32   = xi[i] * L * (1 - xi[i]) * (1 - xi[i]);
	double N33   = 3 * x2 - 2 * x3;
	double N34   = -xi[i] * L * (xi[i] - x2);
	double dN31  = 6 * oneOverL * (x2 - xi[i]); 
	double dN32  = 1 + 3 * x2 - 4 * xi[i];
	double dN33  = 6 * oneOverL * (xi[i] - x2);
	double dN34  = 3 * x2 - 2 * xi[i];
	double ddN31 = oneOverLsquare * (xi12 - 6);
	double ddN32 = oneOverL * (xi6 - 4);
	double ddN33 = oneOverLsquare * (6 - xi12);
	double ddN34 = oneOverL * (xi6 - 2);
	double dN12  = oneOverL;

	N1(0,0)   = 1.0;
        N1(0,2)   = (4 * v( 1) - v( 9)) / 30;
        N1(0,3)   = (4 * v( 2) - v(10)) / 30;
        N1(0,5)   = (4 * v( 9) - v( 1)) / 30;
        N1(0,6)   = (4 * v(10) - v( 2)) / 30;
        N1(1,7)   = 1.0;
        N1(1,8)   =    N31 * v(0) +   N32 * v(3) + N33 * v(8) + N34 * v(11);  // Theta
        N1(1,9)   = -ddN32 * v(2) - ddN34 * v(10); 			      // w dprime
        N1(2,7)   = N1(1,8);
        N1(2,8)   = -1.0;
        N1(2,9)   = ddN32 * v(1) + ddN34 * v(9); 		       	      // v dprime
        N1(3,10)  =  dN31 * v(0) +  dN32 * v(3) + dN33 * v(8) + dN34 * v(11); // Theta prime
        N1(4,11)  = 1.0;
        N1(5,10)  = 1.0; 
	N1(6,12)  = 1.0; 
	N1(7,13)  = 1.0; 
	N1(8,14)  = 1.0; 
	N1(9,15)  = 1.0; 
	N1(10,16) = 1.0; 
	N1(11,17) = 1.0; 

	N2(0, 16)  =   dN12;
	N2(1, 0)   =      1;
	N2(2, 1)   =      1;
	N2(3, 2)   =      1;
	N2(4, 8)   =      1;
	N2(5, 9)   =      1;
	N2(6, 10)  =      1;
	N2(7, 1)   =  ddN32;  // rz1->v dprime
	N2(7, 9)   =  ddN34;  // rz2->v dprime
	N2(8, 2)   = -ddN32;  // ry1->w dprime
	N2(8, 10)  = -ddN34;  // ry2->w dprime
	N2(9, 0)   =    N31;  // rx1->theta
	N2(9, 3)   =    N32;  // rx1 prime->theta
	N2(9, 8)   =    N33;  // rx2->theta
	N2(9, 11)  =    N34;  // rx2 prime->theta
	N2(10, 0)  =   dN31;  // rx1->theta prime
	N2(10, 3)  =   dN32;  // rx1 prime->theta prime
	N2(10, 8)  =   dN33;  // rx2->theta prime
	N2(10, 11) =   dN34;  // rx2 prime->theta prime
	N2(11, 0)  =  ddN31;  // rx1->theta dprime
	N2(11, 3)  =  ddN32;  // rx1 prime->theta dprime
	N2(11, 8)  =  ddN33;  // rx2->theta dprime
	N2(11, 11) =  ddN34;  // rx2 prime->theta dprime 
	N2(12, 4)  =    N31;
	N2(12, 5)  =    N32;
	N2(12, 12) =    N33;
	N2(12, 13) =    N34;
	N2(13, 4)  =   dN31;
	N2(13, 5)  =   dN32;
	N2(13, 12) =   dN33;
	N2(13, 13) =   dN34;
	N2(14, 4)  =  ddN31;
	N2(14, 5)  =  ddN32;
	N2(14, 12) =  ddN33;
	N2(14, 13) =  ddN34;
	N2(15, 6)  =    N31;
	N2(15, 7)  =    N32;
	N2(15, 14) =    N33;
	N2(15, 15) =    N34;
	N2(16, 6)  =   dN31;
	N2(16, 7)  =   dN32;
	N2(16, 14) =   dN33;
	N2(16, 15) =   dN34;
	N2(17, 6)  =  ddN31;
	N2(17, 7)  =  ddN32;
	N2(17, 14) =  ddN33;
	N2(17, 15) =  ddN34;

	N4(0,1)   =   dN32;
	N4(0,9)   =   dN34;
	N4(1,1)   =      1;
	N4(2,2)   =      1;
	N4(3,2)   =  -dN32;
	N4(3,10)  =  -dN34;
	N4(4,9)   =      1;
	N4(4,10)  =      1;
	N4(6,1)   =  ddN32; // rz1->v dprime
	N4(6,9)   =  ddN34; // rz2->v dprime
	N4(7,2)   = -ddN32; // ry1->w dprime
	N4(7,10)  = -ddN34; // ry2->w dprime
	N4(8,0)   =    N31; // rx1->theta
	N4(8,3)   =    N32; // rx1 prime->theta
	N4(8,8)   =    N33; // rx2->theta
	N4(8,11)  =    N34; // rx2 prime->theta
	N4(9,0)   =   dN31; // rx1->theta prime
	N4(9,3)   =   dN32; // rx1 prime->theta prime
	N4(9,8)   =   dN33; // rx2->theta prime
	N4(9,11)  =   dN34; // rx2 prime->theta prime
	N4(10,4)  =    N31;
	N4(10,5)  =    N32;
	N4(10,12) =    N33;
	N4(10,13) =    N34;
	N4(11,4)  =   dN31;
	N4(11,5)  =   dN32;
	N4(11,12) =   dN33;
	N4(11,13) =   dN34;
	N4(12,6)  =    N31;
	N4(12,7)  =    N32;
	N4(12,14) =    N33;
	N4(12,15) =    N34;
	N4(13,6)  =   dN31;
	N4(13,7)  =   dN32;
	N4(13,14) =   dN33;
	N4(13,15) =   dN34;

	#ifdef N*MATRIX
	fprintf(stdout, "TaperedDispBeamColumnSmoothing3d::getInitialBasicStiff:\n");
	fprintf(stdout, "Printing N1 for Section %2d:\n\n", i);
	for (int a = 0; a < 12; ++a) {
	    for (int b = 0; b < 18; ++b)
		fprintf(stdout, "%4.1f ", N1(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "Printing N2 for Section %2d:\n\n", i);
	for (int a = 0; a < 18; ++a) {
	    for (int b = 0; b < 17; ++b)
		fprintf(stdout, "%4.1f ", N2(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "Printing N4 for Section %2d:\n\n", i);
	for (int a = 0; a < 14; ++a) {
	    for (int b = 0; b < 17; ++b)
		fprintf(stdout, "%4.1f ", N4(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");
	#endif

        // Get the section tangent stiffness and stress resultant
        const Matrix & ks = theSections[i]->getInitialTangent();

        N3.addMatrixTripleProduct(0.0, N1, ks, 1.0);
        kbPart1.addMatrixTripleProduct(0.0, N2, N3, 1.0);

        const Vector & s = theSections[i]->getStressResultant();

	// s needs to be in the following order:
	// [P Mz My W B Tsv]
	 Gmax(1,1)   = (4/30)*s(0);
	 Gmax(2,2)   = Gmax(1,1);
	 Gmax(4,4)   = Gmax(1,1);
	 Gmax(5,5)   = Gmax(1,1);
	 Gmax(1,4)   = -s(0)/30;
	 Gmax(4,1)   = Gmax(1,4);
	 Gmax(2,5)   = Gmax(1,4);
	 Gmax(5,2)   = Gmax(1,4);
	 Gmax(6,8)   = s(2);
	 Gmax(8,6)   = Gmax(6,8);
	 Gmax(7,8)   = s(1);
	 Gmax(8,7)   = Gmax(7,8);
	 Gmax(10,10) = s(3); 
	 // Couple Portion of Gmax
	 Gmax(11,0) = s(12);
	 Gmax(0,11) = Gmax(11,0);
	 Gmax(13,0) = s(13);
	 Gmax(0,13) = Gmax(13,0);
	 Gmax(10,3) = s(14);
	 Gmax(3,10) = Gmax(10,3);
	 Gmax(11,3) = s(15);
	 Gmax(3,11) = Gmax(11,3);
	 Gmax(12,3) = s(16);
	 Gmax(3,12) = Gmax(12,3);
	 Gmax(13,3) = s(17);
	 Gmax(3,13) = Gmax(13,3);
	 Gmax(10,8) = s(18);
	 Gmax(8,10) = Gmax(10,8);
	 Gmax(11,8) = s(19);
	 Gmax(8,11) = Gmax(11,8);
	 Gmax(12,8) = s(20);
	 Gmax(8,12) = Gmax(12,8);
	 Gmax(13,8) = s(21);
	 Gmax(8,13) = Gmax(13,8);
	 Gmax(10,9) = s(22);
	 Gmax(9,10) = Gmax(10,9);
	 Gmax(11,9) = s(23);
	 Gmax(9,11) = Gmax(11,9);
	 Gmax(12,9) = s(24);
	 Gmax(9,12) = Gmax(12,9);
	 Gmax(13,9) = s(25);
	 Gmax(9,13) = Gmax(13,9);
	 Gmax(10,10) = s(26);
	 Gmax(11,10) = s(27);
	 Gmax(10,11) = Gmax(11,10);
	 Gmax(12,10) = s(28);
	 Gmax(10,12) = Gmax(12,10);
	 Gmax(13,9)  = s(29);
	 Gmax(9,13)  = Gmax(13,9);
	 Gmax(11,11) = s(30);
	 Gmax(12,11) = s(31);
	 Gmax(11,12) = Gmax(12,11);
	 Gmax(13,11) = s(32);
	 Gmax(11,13) = Gmax(13,11);
	 Gmax(12,12) = s(33);
	 Gmax(13,12) = s(34);
	 Gmax(12,13) = Gmax(13,12);
	 Gmax(13,13) = s(35);

        kbPart2.addMatrixTripleProduct(0.0, N4, Gmax, 1.0); 

        // Perform numerical integration
        double wti = wt[i];

        for (int j = 0; j < 17; j++) {
            for (int k = 0; k < 17; k++) {
                kb(j, k) += kbPart1(j, k) * L * wti + kbPart2(j, k) * L * wti;
            }
        }

	#ifdef STIFFMATRIX
	fprintf(stdout, "TaperedDispBeamColumnSmoothing3d::getInitialBasicStiff:\n");
	fprintf(stdout, "Printing Section Stiffness Matrix for Section %2d:\n", i);
	for (int a = 0; a < 12; ++a) {
	    for (int b = 0; b < 12; ++b)
		fprintf(stdout, "%10.2e ", ks(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "Printing Material Stiffness Matrix for Section %2d:\n", i);
	for (int a = 0; a < 17; ++a) {
	    for (int b = 0; b < 17; ++b)
		fprintf(stdout, "%10.2e ", kbPart1(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "Printing Geometric Stiffness Matrix for Section %2d:\n", i);
	for (int a = 0; a < 17; ++a) {
	    for (int b = 0; b < 17; ++b)
		fprintf(stdout, "%10.2e ", kbPart2(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "Printing BASIC Stiffness Matrix with %2d sections:\n", i+1);
	for (int a = 0; a < 17; ++a) {
	    for (int b = 0; b < 17; ++b)
		fprintf(stdout, "%10.2e ", kb(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");
	#endif

    }

    return kb;
}

const Matrix &TaperedDispBeamColumnSmoothing3d::getInitialStiff()
{
    const Matrix &kb = this->getInitialBasicStiff();

    // Transform to global stiffness
    K = crdTransf->getInitialGlobalStiffMatrix(kb);

    #ifdef STIFFMATRIX
    fprintf(stdout, "TaperedDispBeamColumnSmoothing3d::getInitialStiff:\n");
    fprintf(stdout, "Printing K (global) Matrixd:\n\n");
    for (int a = 0; a < 22; ++a) {
	for (int b = 0; b < 22; ++b)
	    fprintf(stdout, "%10.2e ", K(a, b));
	fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
    #endif

    return K;
}

const Matrix &TaperedDispBeamColumnSmoothing3d::getMass()
{

    K.Zero();

    if (rho == 0.0)
        return K;

    double L = crdTransf->getInitialLength();
    double m = 0.5 * rho * L; // this just applies half half the mass 
                              // each node

    K(0, 0) = K(1, 1) = K(2, 2) = K(11, 11) = K(12, 12) = K(13, 13) = m;

    return K;
}

void TaperedDispBeamColumnSmoothing3d::zeroLoad(void)
{
    Q.Zero();

    for (int i = 0; i < 17; i++)
	q0[i] = 0.0;

    for (int i = 0; i < 5; i++)
	p0[i] = 0.0;

    return;
}

int TaperedDispBeamColumnSmoothing3d::addLoad(ElementalLoad * theLoad,
                                     double loadFactor)
{

    int type;
    const Vector &data = theLoad->getData(type, loadFactor);
    double L = crdTransf->getInitialLength();

    if (type == LOAD_TAG_Beam3dUniformLoad) {
        double wy = data(0) * loadFactor;       // Transverse
        double wz = data(1) * loadFactor;       // Transverse
        double wx = data(2) * loadFactor;       // Axial (+ve from node I to J)

        double Vy = 0.5 * wy * L;
        double Mz = Vy * L / 6.0;       // wy*L*L/12
        double Vz = 0.5 * wz * L;
        double My = Vz * L / 6.0;       // wz*L*L/12
        double P = wx * L;

        // Reactions in basic system
        p0[0] -= P;
        p0[1] -= Vy;
        p0[2] -= Vy;
        p0[3] -= Vz;
        p0[4] -= Vz;

        // Fixed end forces in basic system
        q0[0] -= 0.5 * P;
        q0[1] -= Mz;
        q0[2] += Mz;
        q0[3] += My;
        q0[4] -= My;
    }
    else if (type == LOAD_TAG_Beam3dPointLoad) {
        double Py = data(0) * loadFactor;
        double Pz = data(1) * loadFactor;
        double N = data(2) * loadFactor;
        double aOverL = data(3);

        if (aOverL < 0.0 || aOverL > 1.0)
            return 0;

        double a = aOverL * L;
        double b = L - a;

        // Reactions in basic system
        p0[0] -= N;
        double V1, V2;
        V1 = Py * (1.0 - aOverL);
        V2 = Py * aOverL;
        p0[1] -= V1;
        p0[2] -= V2;
        V1 = Pz * (1.0 - aOverL);
        V2 = Pz * aOverL;
        p0[3] -= V1;
        p0[4] -= V2;

        double L2 = 1.0 / (L * L);
        double a2 = a * a;
        double b2 = b * b;

        // Fixed end forces in basic system
        q0[0] -= N * aOverL;
        double M1, M2;
        M1 = -a * b2 * Py * L2;
        M2 = a2 * b * Py * L2;
        q0[1] += M1;
        q0[2] += M2;
        M1 = -a * b2 * Pz * L2;
        M2 = a2 * b * Pz * L2;
        q0[3] -= M1;
        q0[4] -= M2;
    }
    else {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::addLoad() -- load type unknown for element with tag: "
            << this->getTag() << endln;
        return -1;
    }

    return 0;
}

int TaperedDispBeamColumnSmoothing3d::
addInertiaLoadToUnbalance(const Vector &accel)
{
    // Check for a quick return
    if (rho == 0.0)
        return 0;

    // Get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);

    if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
        return -1;
    }

    double L = crdTransf->getInitialLength();
    double m = 0.5 * rho * L;

    // Want to add ( - fact * M R * accel ) to unbalance
    // Take advantage of lumped mass matrix
    Q(0)  -= m * Raccel1(0);
    Q(1)  -= m * Raccel1(1);
    Q(2)  -= m * Raccel1(2);
    Q(11) -= m * Raccel2(0);
    Q(12) -= m * Raccel2(1);
    Q(13) -= m * Raccel2(2);

    return 0;
}

const Vector &TaperedDispBeamColumnSmoothing3d::getResistingForce()
{

    double L = crdTransf->getInitialLength();
    double oneOverL       = 1.0 / L;
    double oneOverLsquare = oneOverL / L;

    static Matrix N1(12,18);
    static Matrix N2(18,17);

    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);
    double wt[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wt);
    const Vector & v = crdTransf->getBasicTrialDisp();

    // Zero for integration
    q.Zero();

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

        N1.Zero();
        N2.Zero();

	double xi6  = 6.0 * xi[i];
	double xi12 = 12.0 * xi[i];
	double xi1  = xi[i];
	double x3   = xi[i] * xi[i] * xi[i];
	double x2   = xi[i] * xi[i];    

	// N means shape function, dN means 1st derivative, ddN 2nd 
	// N* - * is the order of the shape function 3-cubic 1-linear
	// N*x -  x is the dof of the shape function
	// xi is bringing in x/L already
	double N31   = 1 - 3 * x2 + 2 * x3;
	double N32   = xi[i] * L * (1 - xi[i]) * (1 - xi[i]);
	double N33   = 3 * x2 - 2 * x3;
	double N34   = -xi[i] * L * (xi[i] - x2);
	double dN31  = 6 * oneOverL * (x2 - xi[i]); 
	double dN32  = 1 + 3 * x2 - 4 * xi[i];
	double dN33  = 6 * oneOverL * (xi[i] - x2);
	double dN34  = 3 * x2 - 2 * xi[i];
	double ddN31 = oneOverLsquare * (xi12 - 6);
	double ddN32 = oneOverL * (xi6 - 4);
	double ddN33 = oneOverLsquare * (6 - xi12);
	double ddN34 = oneOverL * (xi6 - 2);
	double dN12  = oneOverL;

	N1(0,0)   = 1.0;
        N1(0,2)   = (4 * v( 1) - v( 9)) / 30;
        N1(0,3)   = (4 * v( 2) - v(10)) / 30;
        N1(0,5)   = (4 * v( 9) - v( 1)) / 30;
        N1(0,6)   = (4 * v(10) - v( 2)) / 30;
        N1(1,7)   = 1.0;
        N1(1,8)   =    N31 * v(0) +   N32 * v(3) + N33 * v(8) + N34 * v(11);  // Theta
        N1(1,9)   = -ddN32 * v(2) - ddN34 * v(10); 			      // w dprime
        N1(2,7)   = N1(1,8);
        N1(2,8)   = -1.0;
        N1(2,9)   = ddN32 * v(1) + ddN34 * v(9); 		       	      // v dprime
        N1(3,10)  =  dN31 * v(0) +  dN32 * v(3) + dN33 * v(8) + dN34 * v(11); // Theta prime
        N1(4,11)  = 1.0;
        N1(5,10)  = 1.0; 
	N1(6,12)  = 1.0; 
	N1(7,13)  = 1.0; 
	N1(8,14)  = 1.0; 
	N1(9,15)  = 1.0; 
	N1(10,16) = 1.0; 
	N1(11,17) = 1.0; 

	N2(0, 16)  =   dN12;
	N2(1, 0)   =      1;
	N2(2, 1)   =      1;
	N2(3, 2)   =      1;
	N2(4, 8)   =      1;
	N2(5, 9)   =      1;
	N2(6, 10)  =      1;
	N2(7, 1)   =  ddN32;  // rz1->v dprime
	N2(7, 9)   =  ddN34;  // rz2->v dprime
	N2(8, 2)   = -ddN32;  // ry1->w dprime
	N2(8, 10)  = -ddN34;  // ry2->w dprime
	N2(9, 0)   =    N31;  // rx1->theta
	N2(9, 3)   =    N32;  // rx1 prime->theta
	N2(9, 8)   =    N33;  // rx2->theta
	N2(9, 11)  =    N34;  // rx2 prime->theta
	N2(10, 0)  =   dN31;  // rx1->theta prime
	N2(10, 3)  =   dN32;  // rx1 prime->theta prime
	N2(10, 8)  =   dN33;  // rx2->theta prime
	N2(10, 11) =   dN34;  // rx2 prime->theta prime
	N2(11, 0)  =  ddN31;  // rx1->theta dprime
	N2(11, 3)  =  ddN32;  // rx1 prime->theta dprime
	N2(11, 8)  =  ddN33;  // rx2->theta dprime
	N2(11, 11) =  ddN34;  // rx2 prime->theta dprime 
	N2(12, 4)  =    N31;
	N2(12, 5)  =    N32;
	N2(12, 12) =    N33;
	N2(12, 13) =    N34;
	N2(13, 4)  =   dN31;
	N2(13, 5)  =   dN32;
	N2(13, 12) =   dN33;
	N2(13, 13) =   dN34;
	N2(14, 4)  =  ddN31;
	N2(14, 5)  =  ddN32;
	N2(14, 12) =  ddN33;
	N2(14, 13) =  ddN34;
	N2(15, 6)  =    N31;
	N2(15, 7)  =    N32;
	N2(15, 14) =    N33;
	N2(15, 15) =    N34;
	N2(16, 6)  =   dN31;
	N2(16, 7)  =   dN32;
	N2(16, 14) =   dN33;
	N2(16, 15) =   dN34;
	N2(17, 6)  =  ddN31;
	N2(17, 7)  =  ddN32;
	N2(17, 14) =  ddN33;
	N2(17, 15) =  ddN34;

        // Get section stress resultant
        const Vector &s = theSections[i]->getStressResultant();
        //opserr<<"s"<<s << endln;

        double wti = wt[i];

	// Establish the section forces before transforming to nodal forces
	static Vector sForces(12);
	sForces.Zero();

	for (int j = 0; j < 12; j++)
	    sForces(j) = s(j);

	//opserr << "section forces: " << sForces << endln;

        static Vector qProduct1(18);
        static Vector qProduct2(17);
        qProduct1.Zero();
        qProduct2.Zero();
        qProduct1.addMatrixTransposeVector(0.0, N1, sForces, 1.0);
        qProduct2.addMatrixTransposeVector(0.0, N2, qProduct1, 1.0);

        for (int j = 0; j < 17; j++) 
            q(j) += qProduct2(j) * L * wti;

	#ifdef N*MATRIX
	fprintf(stdout, "TaperedDispBeamColumnSmoothing3d::getResistingForce:\n");
	fprintf(stdout, "Printing N1 for Section %2d:\n\n", i);
	for (int a = 0; a < 12; ++a) {
	    for (int b = 0; b < 18; ++b)
		fprintf(stdout, "%4.1f ", N1(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "Printing N2 for Section: %2d\n\n", i);
	for (int a = 0; a < 18; ++a) {
	    for (int b = 0; b < 17; ++b)
		fprintf(stdout, "%4.1f ", N2(a,b));
	    fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");
	#endif

	#ifdef DEBUG
	/*
	fprintf(stdout, "Element %2d: Printing sForce for section %2d\n", this->getTag(), i+1); 
	for (int a = 0; a < 12; ++a) 
	    fprintf(stdout, "sForce[%2d] = %10.2e\n", a+1, sForces(a));
	fprintf(stdout, "\n");

	fprintf(stdout, "Printing N1'*sForce for section: %2d\n", i);
	for (int a = 0; a < 18; ++a) 
	    fprintf(stdout, "qProduct1[%2d] = %.4f\n", a+1, qProduct1(a));
	fprintf(stdout, "\n");
	
	fprintf(stdout, "Printing Basic Forces for section: %2d\n", i+1);
	for (int a = 0; a < 17; ++a) 
	    fprintf(stdout, "qProduct2[%2d] = %.2f\n", a+1, qProduct2(a));
	fprintf(stdout, "\n");
	*/
	#endif
    }

    #ifdef DEBUG
    //opserr << "Element %2d Basic Element Forces: ", this->getTag() << q << endln;
    opserr << "Element " << this->getTag() << " Basic Element Forces: " << q << endln;

    #endif

    // Transform forces
    Vector p0Vec(p0, 5); // This is for uniform loads; to be added in future work
    P = crdTransf->getGlobalResistingForce(q, p0Vec);

    // Subtract other external nodal loads ... P_res = P_int - P_ext
    P.addVector(1.0, Q, -1.0); 

    /*
    fprintf(stdout, "TaperedDispBeamColumnSmoothing3d (tag: %d): Element Resisting Force\n", this->getTag());
    for (int a = 0; a < 22; ++a) 
	fprintf(stdout, "%10.2f\n", P(a));
    fprintf(stdout, "\n");
    */

    return P;
}

const Vector & TaperedDispBeamColumnSmoothing3d::getResistingForceIncInertia()
{

    this->getResistingForce();

    if (rho != 0.0) {
        const Vector & accel1 = theNodes[0]->getTrialAccel();
        const Vector & accel2 = theNodes[1]->getTrialAccel();

        // Compute the current resisting force
        this->getResistingForce();

        double L = crdTransf->getInitialLength();
        double m = 0.5 * rho * L;

        P(0)  += m * accel1(0);
        P(1)  += m * accel1(1);
        P(2)  += m * accel1(2);
        P(11) += m * accel2(0);
        P(12) += m * accel2(1);
        P(13) += m * accel2(2);

        // add the damping forces if rayleigh damping
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0
            || betaKc != 0.0)
            P += this->getRayleighDampingForces();

    }
    else {

        // add the damping forces if rayleigh damping
        if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
            P += this->getRayleighDampingForces();
    }

    
    return P;
}

int TaperedDispBeamColumnSmoothing3d::sendSelf(int commitTag,
                                      Channel & theChannel)
{

    // place the integer data into an ID

    int dbTag = this->getDbTag();
    int i, j;
    int loc = 0;

    static ID idData(7);        // one bigger than needed so no clash later
    idData(0) = this->getTag();
    idData(1) = connectedExternalNodes(0);
    idData(2) = connectedExternalNodes(1);
    idData(3) = numSections;
    idData(4) = crdTransf->getClassTag();
    int crdTransfDbTag = crdTransf->getDbTag();
    if (crdTransfDbTag == 0) {
        crdTransfDbTag = theChannel.getDbTag();
        if (crdTransfDbTag != 0)
            crdTransf->setDbTag(crdTransfDbTag);
    }
    idData(5) = crdTransfDbTag;

    if (alphaM != 0 || betaK != 0 || betaK0 != 0 || betaKc != 0)
        idData(6) = 1;
    else
        idData(6) = 0;


    if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::sendSelf() - failed to send ID data\n";
        return -1;
    }

    if (idData(6) == 1) {
        // send damping coefficients
        static Vector dData(4);
        dData(0) = alphaM;
        dData(1) = betaK;
        dData(2) = betaK0;
        dData(3) = betaKc;
        if (theChannel.sendVector(dbTag, commitTag, dData) < 0) {
            opserr <<
                "TaperedDispBeamColumnSmoothing3d::sendSelf() - failed to send double data\n";
            return -1;
        }
    }

    // send the coordinate transformation
    if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::sendSelf() - failed to send crdTranf\n";
        return -1;
    }

    //
    // send an ID for the sections containing each sections dbTag and classTag
    // if section ha no dbTag get one and assign it
    //

    ID idSections(2 * numSections);
    loc = 0;
    for (i = 0; i < numSections; i++) {
        int sectClassTag = theSections[i]->getClassTag();
        int sectDbTag = theSections[i]->getDbTag();
        if (sectDbTag == 0) {
            sectDbTag = theChannel.getDbTag();
            theSections[i]->setDbTag(sectDbTag);
        }

        idSections(loc) = sectClassTag;
        idSections(loc + 1) = sectDbTag;
        loc += 2;
    }

    if (theChannel.sendID(dbTag, commitTag, idSections) < 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::sendSelf() - failed to send ID data\n";
        return -1;
    }

    //
    // send the sections
    //

    for (j = 0; j < numSections; j++) {
        if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
            opserr << "TaperedDispBeamColumnSmoothing3d::sendSelf() - section "
                << j << "failed to send itself\n";
            return -1;
        }
    }

    return 0;
}

int TaperedDispBeamColumnSmoothing3d::recvSelf(int commitTag,
                                      Channel & theChannel,
                                      FEM_ObjectBroker & theBroker)
{

    //
    // receive the integer data containing tag, numSections and coord transformation info
    //
    int dbTag = this->getDbTag();
    int i;

    static ID idData(7);        // one bigger than needed so no clash with section ID

    if (theChannel.recvID(dbTag, commitTag, idData) < 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::recvSelf() - failed to recv ID data\n";
        return -1;
    }

    this->setTag(idData(0));
    connectedExternalNodes(0) = idData(1);
    connectedExternalNodes(1) = idData(2);

    int crdTransfClassTag = idData(4);
    int crdTransfDbTag = idData(5);

    if (idData(6) == 1) {
        // recv damping coefficients
        static Vector dData(4);
        if (theChannel.recvVector(dbTag, commitTag, dData) < 0) {
            opserr <<
                "TaperedDispBeamColumnSmoothing3d::sendSelf() - failed to recv double data\n";
            return -1;
        }
        alphaM = dData(0);
        betaK = dData(1);
        betaK0 = dData(2);
        betaKc = dData(3);
    }

    // create a new crdTransf object if one needed
    if (crdTransf == 0
        || crdTransf->getClassTag() != crdTransfClassTag) {
        if (crdTransf != 0)
            delete crdTransf;

        crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

        if (crdTransf == 0) {
            opserr << "TaperedDispBeamColumnSmoothing3d::recvSelf() - " <<
                "failed to obtain a CrdTrans object with classTag" <<
                crdTransfClassTag << endln;
            return -2;
        }
    }

    crdTransf->setDbTag(crdTransfDbTag);

    // invoke recvSelf on the crdTransf object
    if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::sendSelf() - failed to recv crdTranf\n";
        return -3;
    }

    //
    // recv an ID for the sections containing each sections dbTag and classTag
    //

    ID idSections(2 * idData(3));
    int loc = 0;

    if (theChannel.recvID(dbTag, commitTag, idSections) < 0) {
        opserr <<
            "TaperedDispBeamColumnSmoothing3d::recvSelf() - failed to recv ID data\n";
        return -1;
    }

    //
    // now receive the sections
    //

    if (numSections != idData(3)) {

        //
        // we do not have correct number of sections, must delete the old and create
        // new ones before can recvSelf on the sections
        //

        // delete the old
        if (numSections != 0) {
            for (int i = 0; i < numSections; i++)
                delete theSections[i];
            delete[]theSections;
        }

        // create a new array to hold pointers
        theSections = new SectionForceDeformation *[idData(3)];
        if (theSections == 0) {
            opserr <<
                "TaperedDispBeamColumnSmoothing3d::recvSelf() - out of memory creating sections array of size"
                << idData(3) << endln;
            exit(-1);
        }

        // create a section and recvSelf on it
        numSections = idData(3);
        loc = 0;

        for (i = 0; i < numSections; i++) {
            int sectClassTag = idSections(loc);
            int sectDbTag = idSections(loc + 1);
            loc += 2;
            theSections[i] = theBroker.getNewSection(sectClassTag);
            if (theSections[i] == 0) {
                opserr <<
                    "TaperedDispBeamColumnSmoothing3d::recvSelf() - Broker could not create Section of class type"
                    << sectClassTag << endln;
                exit(-1);
            }
            theSections[i]->setDbTag(sectDbTag);
            if (theSections[i]->
                recvSelf(commitTag, theChannel, theBroker) < 0) {
                opserr <<
                    "TaperedDispBeamColumnSmoothing3d::recvSelf() - section "
                    << i << "failed to recv itself\n";
                return -1;
            }
        }

    }
    else {

        // 
        // for each existing section, check it is of correct type
        // (if not delete old & create a new one) then recvSelf on it
        //

        loc = 0;
        for (i = 0; i < numSections; i++) {
            int sectClassTag = idSections(loc);
            int sectDbTag = idSections(loc + 1);
            loc += 2;

            // check of correct type
            if (theSections[i]->getClassTag() != sectClassTag) {
                // delete the old section[i] and create a new one
                delete theSections[i];
                theSections[i] =
                    theBroker.getNewSection(sectClassTag);
                if (theSections[i] == 0) {
                    opserr <<
                        "TaperedDispBeamColumnSmoothing3d::recvSelf() - Broker could not create Section of class type"
                        << sectClassTag << endln;
                    exit(-1);
                }
            }

            // recvSelf on it
            theSections[i]->setDbTag(sectDbTag);
            if (theSections[i]->
                recvSelf(commitTag, theChannel, theBroker) < 0) {
                opserr <<
                    "TaperedDispBeamColumnSmoothing3d::recvSelf() - section "
                    << i << "failed to recv itself\n";
                return -1;
            }
        }
    }

    return 0;
}

void TaperedDispBeamColumnSmoothing3d::Print(OPS_Stream &s, int flag)
{
    s << "\nTaperedDispBeamColumnSmoothing3d, element id:  " << this->
        getTag() << endln;
    s << "\tConnected external nodes:  " << connectedExternalNodes;
    s << "\tmass density:  " << rho << endln;

    double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
    double Bi1, Bi2;
    double Alpha1, dAlpha1, Beta1, dBeta1, Alpha2, dAlpha2, Beta2, dBeta2;
    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    N   = q(16);
    Mz1 = q(1);
    Mz2 = q(9);
    Vy  = (Mz1 + Mz2) * oneOverL;
    My1 = q(2);
    My2 = q(10);
    Vz  = -(My1 + My2) * oneOverL;
    T   = q(8);
    // Additional DOF printout values
    Bi1     = q(3);
    Bi2     = q(11);
    Alpha1  = q(4);
    dAlpha1 = q(5);
    Beta1   = q(6);
    dBeta1  = q(7);
    Alpha2  = q(12);
    dAlpha2 = q(13);
    Beta2   = q(14);
    dBeta2  = q(15);


    s << "\tEnd 1 Forces (P Mz Vy My Vz T; Bi TF TF' BF BF'): "
        << -N + p0[0] << ' ' << Mz1 << ' ' << Vy +
        p0[1] << ' ' << My1 << ' ' << Vz +
        p0[3] << ' ' << -T << "\n \t" 
	<< Bi1 << ' ' << Alpha1 << ' ' << dAlpha1 
	<< ' ' << Beta1 << ' ' << dBeta1 << endln;
    s << "\tEnd 2 Forces (P Mz Vy My Vz T; Bi TF TF' BF BF'): " 
	<< N << ' ' << Mz2 << ' ' << -Vy + 
	p0[2] << ' ' << My2 << ' ' << -Vz +
        p0[4] << ' ' << T << "\n \t" 
	<< Bi2 << ' ' << Alpha2 << ' ' << dAlpha2 
	<< ' ' << Beta2 << ' ' << dBeta2 << endln;

}


int TaperedDispBeamColumnSmoothing3d::displaySelf(Renderer & theViewer,
                                         int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    const Vector & end1Crd = theNodes[0]->getCrds();
    const Vector & end2Crd = theNodes[1]->getCrds();

    static Vector v1(3);
    static Vector v2(3);

    if (displayMode >= 0) {
        const Vector & end1Disp = theNodes[0]->getDisp();
        const Vector & end2Disp = theNodes[1]->getDisp();

        for (int i = 0; i < 3; i++) {
            v1(i) = end1Crd(i) + end1Disp(i) * fact;
            v2(i) = end2Crd(i) + end2Disp(i) * fact;
        }
    }
    else {
        int mode = displayMode * -1;
        const Matrix & eigen1 = theNodes[0]->getEigenvectors();
        const Matrix & eigen2 = theNodes[1]->getEigenvectors();
        if (eigen1.noCols() >= mode) {
            for (int i = 0; i < 3; i++) {
                v1(i) = end1Crd(i) + eigen1(i, mode - 1) * fact;
                v2(i) = end2Crd(i) + eigen2(i, mode - 1) * fact;
            }

        }
        else {
            for (int i = 0; i < 3; i++) {
                v1(i) = end1Crd(i);
                v2(i) = end2Crd(i);
            }
        }
    }
    return theViewer.drawLine(v1, v2, 1.0, 1.0);
}

Response *TaperedDispBeamColumnSmoothing3d::setResponse(const char **argv,
                                               int argc,
                                               OPS_Stream & output)
{

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType", "TaperedDispBeamColumnSmoothing3d");
    output.attr("eleTag", this->getTag());
    output.attr("node1", connectedExternalNodes[0]);
    output.attr("node2", connectedExternalNodes[1]);

    //
    // we compare argv[0] for known response types 
    //

    // global force - 
    if (strcmp(argv[0], "forces") == 0
        || strcmp(argv[0], "force") == 0
        || strcmp(argv[0], "globalForce") == 0
        || strcmp(argv[0], "globalForces") == 0) {

        output.tag("ResponseType", "Px_1");
        output.tag("ResponseType", "Py_1");
        output.tag("ResponseType", "Pz_1");
        output.tag("ResponseType", "Mx_1");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "Mz_1");
        output.tag("ResponseType", "Bx_1");
        output.tag("ResponseType", "TF_1");
        output.tag("ResponseType", "dTF_1");
        output.tag("ResponseType", "BF_1");
        output.tag("ResponseType", "dBF_1");

        output.tag("ResponseType", "Px_2");
        output.tag("ResponseType", "Py_2");
        output.tag("ResponseType", "Pz_2");
        output.tag("ResponseType", "Mx_2");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "Mz_2");
	output.tag("ResponseType", "Bx_2");
        output.tag("ResponseType", "TF_2");
        output.tag("ResponseType", "dTF_2");
        output.tag("ResponseType", "BF_2");
        output.tag("ResponseType", "dBF_2");


        theResponse = new ElementResponse(this, 1, P);

    }
    // local force:
    else if (strcmp(argv[0], "localForce") == 0
             || strcmp(argv[0], "localForces") == 0) {

        output.tag("ResponseType", "N_1");
        output.tag("ResponseType", "Vy_1");
        output.tag("ResponseType", "Vz_1");
        output.tag("ResponseType", "T_1");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "Mz_1");
        output.tag("ResponseType", "Bx_1");
        output.tag("ResponseType", "TF_1");
        output.tag("ResponseType", "dTF_1");
        output.tag("ResponseType", "BF_1");
        output.tag("ResponseType", "dBF_1");

        output.tag("ResponseType", "N_2");
        output.tag("ResponseType", "Vy_2");
        output.tag("ResponseType", "Vz_2");
        output.tag("ResponseType", "T_2");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "Mz_2");
        output.tag("ResponseType", "Bx_2");
        output.tag("ResponseType", "TF_2");
        output.tag("ResponseType", "dTF_2");
        output.tag("ResponseType", "BF_2");
        output.tag("ResponseType", "dBF_2");

        theResponse = new ElementResponse(this, 2, P);

    }
    // chord rotation:
    else if (strcmp(argv[0], "chordRotation") == 0
             || strcmp(argv[0], "chordDeformation") == 0
             || strcmp(argv[0], "basicDeformation") == 0) {

        output.tag("ResponseType", "eps");
        output.tag("ResponseType", "thetaZ_1");
        output.tag("ResponseType", "thetaZ_2");
        output.tag("ResponseType", "thetaY_1");
        output.tag("ResponseType", "thetaY_2");
        output.tag("ResponseType", "thetaX");

        theResponse = new ElementResponse(this, 3, Vector(6));

    }
    // plastic rotation:
    else if (strcmp(argv[0], "plasticRotation") == 0
             || strcmp(argv[0], "plasticDeformation") == 0) {

        output.tag("ResponseType", "epsP");
        output.tag("ResponseType", "thetaZP_1");
        output.tag("ResponseType", "thetaZP_2");
        output.tag("ResponseType", "thetaYP_1");
        output.tag("ResponseType", "thetaYP_2");
        output.tag("ResponseType", "thetaXP");

        theResponse = new ElementResponse(this, 4, Vector(6));

    }
    // section response:
    else if (strstr(argv[0], "sectionX") != 0) {
        if (argc > 2) {
            float sectionLoc = atof(argv[1]);

            double xi[maxNumSections];
            double L = crdTransf->getInitialLength();
            beamInt->getSectionLocations(numSections, L, xi);

            sectionLoc /= L;

            float minDistance = fabs(xi[0] - sectionLoc);
            int sectionNum = 0;
            for (int i = 1; i < numSections; i++) {
                if (fabs(xi[i] - sectionLoc) < minDistance) {
                    minDistance = fabs(xi[i] - sectionLoc);
                    sectionNum = i;
                }
            }

            output.tag("GaussPointOutput");
            output.attr("number", sectionNum + 1);
            output.attr("eta", xi[sectionNum] * L);

            theResponse =
                theSections[sectionNum]->setResponse(&argv[2],
                                                     argc - 2,
                                                     output);
        }
    }
    else if (strcmp(argv[0], "section") == 0) {
        if (argc > 2) {

            int sectionNum = atoi(argv[1]);
            if (sectionNum > 0 && sectionNum <= numSections) {

                double xi[maxNumSections];
                double L = crdTransf->getInitialLength();
                beamInt->getSectionLocations(numSections, L, xi);

                output.tag("GaussPointOutput");
                output.attr("number", sectionNum);
                output.attr("eta", xi[sectionNum - 1] * L);

                theResponse =
                    theSections[sectionNum - 1]->setResponse(&argv[2],
                                                             argc - 2,
                                                             output);

                output.endTag();
            }
        }
    }

    output.endTag();
    return theResponse;
}

int TaperedDispBeamColumnSmoothing3d::getResponse(int responseID,
                                         Information &eleInfo)
{

    // Beam Response Forces
    double N, V, M1, M2, T1, T2, Bi1, Bi2;
    // Plate Response Forces
    double TFi, dTFi, TFj, dTFj, BFi, dBFi, BFj, dBFj;

    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    if (responseID == 1)
        return eleInfo.setVector(this->getResistingForce());

    else if (responseID == 2) {
	// NOTE: 22 total local forces
	// Axial
	N     = q(16);
        P(0)  = -N;
        P(11) = N;

        // Torsion
	T1    = q(0);
	T2    = q(8);
        P(3)  = T1;
        P(14) = T2;

        // Moments about z and shears along y
        M1    = q(1);
        M2    = q(9);
        P(5)  = M1;
        P(16) = M2;

        V     = (M1 + M2) * oneOverL;
        P(1)  = V;
        P(12) = -V;

	/*
	opserr << "Mz I: " << M1 << endln;
	opserr << "Mz J: " << M2 << endln;
	opserr << "Vy: " << V << endln;
	*/

        // Moments about y and shears along z
        M1    = q(2);
        M2    = q(10);
        P(4)  = M1;
        P(15) = M2;

        V     = -(M1 + M2) * oneOverL;
        P(2)  = V;
        P(13) = -V;

	/*
	opserr << "My I: " << M1 << endln;
	opserr << "My J: " << M2 << endln;
	opserr << "Vz: " << V << endln;
	*/

	// Bi-Moments about x
	Bi1   = q(3);
	Bi2   = q(11);
	P(6)  = Bi1;
	P(17) = Bi2;

	// Top Flange Forces
	TFi   = q(4);
	dTFi  = q(5);
	TFj   = q(12);
	dTFj  = q(13);
	P(7)  = TFi;
	P(8)  = dTFi;
	P(18) = TFj;
	P(19) = dTFj;
	
	// Bottom Flange Forces
	BFi   = q(6);
	dBFi  = q(7);
	BFj   = q(14);
	dBFj  = q(15);
	P(9)  = BFi;
	P(10) = dBFi;
	P(20) = BFj;
	P(21) = dBFj;

	/*
	opserr << "Axial Load: " << N << endln;
	opserr << "Torsion I: " << T1 << endln;
	opserr << "Torsion J: " << T2 << endln;
	opserr << "BiMoment I: " << Bi1 << endln;
	opserr << "BiMoment J: " << Bi2 << endln;
	opserr << "Alpha  I: " << TFi << endln;
	opserr << "Alpha' I: " << dTFi << endln;
	opserr << "Alpha  J: " << TFj << endln;
	opserr << "Alpha' J: " << dTFj << endln;
	opserr << "Beta  I: " << BFi << endln;
	opserr << "Beta' I: " << dBFi << endln;
	opserr << "Beta  J: " << BFj << endln;
	opserr << "Beta' J: " << dBFj << endln;
	*/

        return eleInfo.setVector(P);
    }

    // Chord rotation
    else if (responseID == 3)
        return eleInfo.setVector(crdTransf->getBasicTrialDisp());

    // Plastic rotation
    else if (responseID == 4) {
        static Vector vp(6);
        static Vector ve(6);
        const Matrix & kb = this->getInitialBasicStiff();
        kb.Solve(q, ve);
        vp = crdTransf->getBasicTrialDisp();
        vp -= ve;
        return eleInfo.setVector(vp);
    }

    else
        return -1;
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int TaperedDispBeamColumnSmoothing3d::setParameter(const char **argv, int argc,
                                          Parameter &param)
{
    if (argc < 1)
        return -1;

    // If the parameter belongs to the element itself
    if (strcmp(argv[0], "rho") == 0)
        return param.addObject(1, this);

    if (strstr(argv[0], "sectionX") != 0) {
        if (argc < 3)
            return -1;

        float sectionLoc = atof(argv[1]);

        double xi[maxNumSections];
        double L = crdTransf->getInitialLength();
        beamInt->getSectionLocations(numSections, L, xi);

        sectionLoc /= L;

        float minDistance = fabs(xi[0] - sectionLoc);
        int sectionNum = 0;
        for (int i = 1; i < numSections; i++) {
            if (fabs(xi[i] - sectionLoc) < minDistance) {
                minDistance = fabs(xi[i] - sectionLoc);
                sectionNum = i;
            }
        }
        return theSections[sectionNum]->setParameter(&argv[2],
                                                     argc - 2, param);
    }
    // If the parameter belongs to a section or lower
    if (strstr(argv[0], "section") != 0) {

        if (argc < 3)
            return -1;

        // Get section and material tag numbers from user input
        int paramSectionTag = atoi(argv[1]);

        // Find the right section and call its setParameter method
        int ok = 0;
        for (int i = 0; i < numSections; i++)
            if (paramSectionTag == theSections[i]->getTag())
                ok +=
                    theSections[i]->setParameter(&argv[2], argc - 2,
                                                 param);

        return ok;
    }

    else if (strstr(argv[0], "integration") != 0) {

        if (argc < 2)
            return -1;

        return beamInt->setParameter(&argv[1], argc - 1, param);
    }

    // Default, send to every object
    int ok = 0;
    for (int i = 0; i < numSections; i++)
        ok += theSections[i]->setParameter(argv, argc, param);
    ok += beamInt->setParameter(argv, argc, param);
    return ok;
}

int TaperedDispBeamColumnSmoothing3d::updateParameter(int parameterID,
                                             Information & info)
{
    if (parameterID == 1) {
        rho = info.theDouble;
        return 0;
    }
    else
        return -1;
}


int TaperedDispBeamColumnSmoothing3d::activateParameter(int passedParameterID)
{
    parameterID = passedParameterID;

    return 0;
}

const Matrix &
    TaperedDispBeamColumnSmoothing3d::getKiSensitivity(int gradNumber)
{
    K.Zero();
    return K;
}

const Matrix &
    TaperedDispBeamColumnSmoothing3d::getMassSensitivity(int gradNumber)
{
    K.Zero();
    return K;
}



const Vector &
    TaperedDispBeamColumnSmoothing3d::
getResistingForceSensitivity(int gradNumber)
{
    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
    //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);
    double wt[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wt);

    // Zero for integration
    static Vector dqdh(6);
    dqdh.Zero();

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

        int order = theSections[i]->getOrder();
        const ID & code = theSections[i]->getType();

        //double xi6 = 6.0*pts(i,0);
        double xi6 = 6.0 * xi[i];
        //double wti = wts(i);
        double wti = wt[i];

        // Get section stress resultant gradient
        const Vector & dsdh =
            theSections[i]->getStressResultantSensitivity(gradNumber,
                                                          true);

        // Perform numerical integration on internal force gradient
        double sensi;
        for (int j = 0; j < order; j++) {
            sensi = dsdh(j) * wti;
            switch (code(j)) {
            case SECTION_RESPONSE_P:
                dqdh(0) += sensi;
                break;
            case SECTION_RESPONSE_MZ:
                dqdh(1) += (xi6 - 4.0) * sensi;
                dqdh(2) += (xi6 - 2.0) * sensi;
                break;
            case SECTION_RESPONSE_MY:
                dqdh(3) += (xi6 - 4.0) * sensi;
                dqdh(4) += (xi6 - 2.0) * sensi;
                break;
            case SECTION_RESPONSE_T:
                dqdh(5) += sensi;
                break;
            default:
                break;
            }
        }
    }

    // Transform forces
    static Vector dp0dh(6);     // No distributed loads

    P.Zero();

    //////////////////////////////////////////////////////////////

    if (crdTransf->isShapeSensitivity()) {

        // Perform numerical integration to obtain basic stiffness matrix
        // Some extra declarations
        static Matrix kbmine(6, 6);
        kbmine.Zero();
        q.Zero();

        double tmp;

        int j, k;

        for (int i = 0; i < numSections; i++) {

            int order = theSections[i]->getOrder();
            const ID & code = theSections[i]->getType();

            //double xi6 = 6.0*pts(i,0);
            double xi6 = 6.0 * xi[i];
            //double wti = wts(i);
            double wti = wt[i];

            const Vector & s = theSections[i]->getStressResultant();
            const Matrix & ks = theSections[i]->getSectionTangent();

            Matrix ka(workArea, order, 6);
            ka.Zero();

            double si;
            for (j = 0; j < order; j++) {
                si = s(j) * wti;
                switch (code(j)) {
                case SECTION_RESPONSE_P:
                    q(0) += si;
                    for (k = 0; k < order; k++) {
                        ka(k, 0) += ks(k, j) * wti;
                    }
                    break;
                case SECTION_RESPONSE_MZ:
                    q(1) += (xi6 - 4.0) * si;
                    q(2) += (xi6 - 2.0) * si;
                    for (k = 0; k < order; k++) {
                        tmp = ks(k, j) * wti;
                        ka(k, 1) += (xi6 - 4.0) * tmp;
                        ka(k, 2) += (xi6 - 2.0) * tmp;
                    }
                    break;
                case SECTION_RESPONSE_MY:
                    q(3) += (xi6 - 4.0) * si;
                    q(4) += (xi6 - 2.0) * si;
                    for (k = 0; k < order; k++) {
                        tmp = ks(k, j) * wti;
                        ka(k, 3) += (xi6 - 4.0) * tmp;
                        ka(k, 4) += (xi6 - 2.0) * tmp;
                    }
                    break;
                case SECTION_RESPONSE_T:
                    q(5) += si;
                    for (k = 0; k < order; k++) {
                        ka(k, 5) += ks(k, j) * wti;
                    }
                    break;
                default:
                    break;
                }
            }
            for (j = 0; j < order; j++) {
                switch (code(j)) {
                case SECTION_RESPONSE_P:
                    for (k = 0; k < 6; k++) {
                        kbmine(0, k) += ka(j, k);
                    }
                    break;
                case SECTION_RESPONSE_MZ:
                    for (k = 0; k < 6; k++) {
                        tmp = ka(j, k);
                        kbmine(1, k) += (xi6 - 4.0) * tmp;
                        kbmine(2, k) += (xi6 - 2.0) * tmp;
                    }
                    break;
                case SECTION_RESPONSE_MY:
                    for (k = 0; k < 6; k++) {
                        tmp = ka(j, k);
                        kbmine(3, k) += (xi6 - 4.0) * tmp;
                        kbmine(4, k) += (xi6 - 2.0) * tmp;
                    }
                    break;
                case SECTION_RESPONSE_T:
                    for (k = 0; k < 6; k++) {
                        kbmine(5, k) += ka(j, k);
                    }
                    break;
                default:
                    break;
                }
            }
        }

        const Vector & A_u = crdTransf->getBasicTrialDisp();
        double dLdh = crdTransf->getdLdh();
        double d1overLdh = -dLdh / (L * L);
        // a^T k_s dadh v
        dqdh.addMatrixVector(1.0, kbmine, A_u, d1overLdh);

        // k dAdh u
        const Vector & dAdh_u =
            crdTransf->getBasicTrialDispShapeSensitivity();
        dqdh.addMatrixVector(1.0, kbmine, dAdh_u, oneOverL);

        // dAdh^T q
        P += crdTransf->getGlobalResistingForceShapeSensitivity(q,
                                                                dp0dh,
                                                                gradNumber);
    }

    // A^T (dqdh + k dAdh u)
    P += crdTransf->getGlobalResistingForce(dqdh, dp0dh);

    return P;
}



// NEW METHOD
int TaperedDispBeamColumnSmoothing3d::commitSensitivity(int gradNumber,
                                               int numGrads)
{
    // Get basic deformation and sensitivities
    const Vector & v = crdTransf->getBasicTrialDisp();

    static Vector dvdh(6);
    dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);

    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;
    //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
    double xi[maxNumSections];
    beamInt->getSectionLocations(numSections, L, xi);

    // Some extra declarations
    double d1oLdh = crdTransf->getd1overLdh();

    // Loop over the integration points
    for (int i = 0; i < numSections; i++) {

        int order = theSections[i]->getOrder();
        const ID & code = theSections[i]->getType();

        Vector e(workArea, order);

        //double xi6 = 6.0*pts(i,0);
        double xi6 = 6.0 * xi[i];

        for (int j = 0; j < order; j++) {
            switch (code(j)) {
            case SECTION_RESPONSE_P:
                e(j) = oneOverL * dvdh(0)
                    + d1oLdh * v(0);
                break;
            case SECTION_RESPONSE_MZ:
                e(j) =
                    oneOverL * ((xi6 - 4.0) * dvdh(1) +
                                (xi6 - 2.0) * dvdh(2))
                    + d1oLdh * ((xi6 - 4.0) * v(1) +
                                (xi6 - 2.0) * v(2));
                break;
            case SECTION_RESPONSE_MY:
                e(j) =
                    oneOverL * ((xi6 - 4.0) * dvdh(3) +
                                (xi6 - 2.0) * dvdh(4))
                    + d1oLdh * ((xi6 - 4.0) * v(3) +
                                (xi6 - 2.0) * v(4));
                break;
            case SECTION_RESPONSE_T:
                e(j) = oneOverL * dvdh(5)
                    + d1oLdh * v(5);
                break;
            default:
                e(j) = 0.0;
                break;
            }
        }

        // Set the section deformations
        theSections[i]->commitSensitivity(e, gradNumber, numGrads);
    }

    return 0;
}
// AddingSensitivity:END /////////////////////////////////////////////
