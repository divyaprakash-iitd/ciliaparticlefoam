/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
// My inclusions
#include "sphereToCell.H"
#include <set>
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "diracdelta.H"
// Declare the external subroutine
extern "C" {
    void generatecilia(int *noelpts, double *cdl, double *h, double *LX, double *LY, double *LZ);
    // void sayhello();
    void getpositions(double *pposx, double *pposy, double *pposz, int *noelpts);
    // subroutine getpositions(XC,YC,ZC) bind(C)
    void calculateforcesandmoments(double *pfx, double *pfy, double *pfz, double *pmx, double *pmy, double *pmz, int *noelpts);
    // subroutine calculateforces(FXC,FYC,FZC) bind(C)
    void updatepositions(double *pvx, double *pvy, double *pvz, double *mx, double *my, double *mz, double *dt, int *noelpts);
    // subroutine updatepositions(U,V,W,dt) bind(C)
    // void arraycheck(double *pxyz, int* n);
    // void arraycheck(int *pxyz, int* n);
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    // Box dimensions
    double LX,LY,LZ;
    // Info << mesh.bounds().max()[0] << endl;
    LX = mesh.bounds().max()[0];
    LY = mesh.bounds().max()[1];
    LZ = mesh.bounds().max()[2];
    int noelpts;
    double meshwidth = mesh.C()[1][0] - mesh.C()[0][0];
    double cdl;
    generatecilia(&noelpts,&cdl,&meshwidth,&LX,&LY,&LZ);
    Info << "No. of cilia points: " << noelpts << endl;
    Info << "CDL: " << cdl << endl << "Mesh Width: " << meshwidth << endl;
    // // Initial position of the point source
    scalar pzmid = (mesh.C()[1][0] - mesh.C()[0][0])/2.0;
    // Info << "ZZ: " << pzmid << endl;
    vector rr(0.0,0.0,pzmid); // Define the point


    // // // Create vectors for storing the particle's Positions, Forces and Velocity
    std::vector<double> pposx(noelpts,0), pposy(noelpts,0), pposz(noelpts,0);   //Position
    std::vector<double> pfx(noelpts,0), pfy(noelpts,0), pfz(noelpts,0);         // Force
    std::vector<double> pmx(noelpts,0), pmy(noelpts,0), pmz(noelpts,0);         // Moment
    std::vector<double> pvx(noelpts,0), pvy(noelpts,0), pvz(noelpts,0);         // Velocity
    std::vector<double> mx(noelpts,0), my(noelpts,0), mz(noelpts,0);         // Angular Velocity

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Foam::sphereToCell sph(mesh,rr,0.005,0);

        // Info << "SPH: " << sph.typeName() << endl;

	    //# define omega 0.05
	    //const dimensionedVector mySource("mySource", dimensionSet(0,1,-2,0,0,0,0), 1000*Foam::sin(runTime.value()*omega)*vector(0,1,0));

        /////////////////////////////////////////////////////////////////////////////////////////        
        // Get positions of the cilia
        getpositions(pposx.data(),pposy.data(),pposz.data(),&noelpts); 
        // Calculate the forces and moments in the cilia
        calculateforcesandmoments(pfx.data(),pfy.data(),pfz.data(), pmx.data(),pmy.data(),pmz.data(), &noelpts);
        // Info << pfx << endl << pfy << endl << pfz << endl << pmx << endl << pmy << endl << pmz << endl;
	// Calculate the moments in the cilia : Cilia

        // Create a list of lists to store the neighbours for 
        // every node of the particle
        Foam::DynamicList<Foam::labelList> neighborsList;
        for (int inoelpts = 0; inoelpts < noelpts; ++inoelpts) {
            // Get the node's location
            rr[0] = pposx[inoelpts]; rr[1] = pposy[inoelpts]; rr[2] = pposz[inoelpts];
            #include "getNeighbours.H"
            neighborsList.append(appended);
            appended.clear();
        }


        // Initialize the force source term to zero
        F = F*0;
        mden = mden*0;
        // Go through each of the nodes and then their neighbours
        for (int inoelpts = 0; inoelpts < noelpts; ++inoelpts) {
            // Get the node's location
            rr[0] = pposx[inoelpts]; rr[1] = pposy[inoelpts]; rr[2] = pposz[inoelpts];
            // Get the node's force (Placeholder variables)
            vector pf(0.0,0.0,0.0);
            pf[0] = pfx[inoelpts]; pf[1] = pfy[inoelpts]; pf[2] = pfz[inoelpts];
            vector pm(0.0,0.0,0.0);
            pm[0] = pmx[inoelpts]; pm[1] = pmy[inoelpts]; pm[2] = pmz[inoelpts];
            // Get the neighbors list for that node alone
            labelList singleNodeNeighbors = neighborsList[inoelpts];

            // Iterate through each of the neighbouring cells of the selected node 
            forAll(singleNodeNeighbors,idx) {
                int icell = singleNodeNeighbors[idx];   
                #include "interpolateForces.H"
                #include "interpolateMden.H"
            }
        }

        F = F + 0.5*fvc::curl(mden);
        // F = F*0;
        // Info << "cdl: " << cdl << endl;
	// Calculate forces arising from moments: Cilia
	// Create and interpolate mden from cilia nodes and then calculate curl of mden to obtain these forces
        ///////////////////////////////////////////////////////////////////////////////////////////

        // Momentum predictor
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U) 
          - F // IBM Source Term 
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Calculate vorticity from velocity
        volVectorField W = 0.5*fvc::curl(U);
        // Go through each of the nodes and then their neighbours
        for (int inoelpts = 0; inoelpts < noelpts; ++inoelpts) {
            // Get the node's location
            rr[0] = pposx[inoelpts];
            rr[1] = pposy[inoelpts];
            rr[2] = pposz[inoelpts];
            // Get the neighbors list for that node alone
            labelList singleNodeNeighbors = neighborsList[inoelpts];

            // Define a variable to store the point source velocity
            vector pu(0,0,0);
            // Define a variable to store the point source moment
            vector pmm(0,0,0);
            // Iterate through each of the neighbouring cells of the selected node 
            forAll(singleNodeNeighbors,idx) {
                int icell = singleNodeNeighbors[idx];   
                #include "interpolateVelocity.H"
                #include "interpolateMoment.H"
            }
            // Store the interpolated velocity at each node
            pvx[inoelpts] = pu[0]; 
            pvy[inoelpts] = pu[1];
            pvz[inoelpts] = pu[2];
            // Store the interpolated moment at each node
            pmx[inoelpts] = pmm[0]; 
            pmy[inoelpts] = pmm[1];
            pmz[inoelpts] = pmm[2];
        }


	// Calculate angular velocity of fluid and interpolate it to the cilia nodes: Cilia

        // Clear the list of lists of neighbours
        neighborsList.clear();
        double simdt = runTime.deltaTValue();

        // Update the position of all the nodes
        updatepositions(pvx.data(),pvy.data(),pvz.data(),pmx.data(),pmy.data(),pmz.data(),&simdt,&noelpts);
        // Update the orientation of all the cilia nodes: Cilia
        ///////////////////////////////////////////////////////////////////////////////////////////

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
