/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    chamberFoam

Group
    grpCompressibleSolvers

Description
    Extension of rhoCentralFoam tailored for chamber dynamics studies.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
// #include "psiThermo.H"
// #include "rhoThermo.H"
//#include "myThermo.H"
#include "myRhoThermo.H"
//#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "lookupTables2D.H"
#include "myRadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Density-based compressible flow solver based on"
        " central-upwind schemes of Kurganov and Tadmor with"
        " support for mesh-motion and topology changes."
    );


    #define NO_CONTROL
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"

    //turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    const dimensionedScalar v_zero(dimVolume/dimTime, Zero);

    bool hydrodynamics = runTime.controlDict().getOrDefault<bool>
        (
            "hydrodynamics",
            true
        );

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    /*IOdictionary thermoProperties
    (
        IOobject
        (
            "thermoPhysicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dictionary& dict(thermoProperties);
    scalarLookupTable2D TTable_(dict.subDict("thermodynamics"), "rho", "e", "T");*/

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        Info << "p[1] " << p[1] << endl;
        Info << "T[1] " << T[1] << endl;
        Info << "rho[1] " << rho[1] << endl;
        Info << "e[1] " << e[1] << endl;


        #include "readTimeControls.H"

        if (!LTS)
        {
            #include "setDeltaT.H"

            ++runTime;

            // Do any mesh changes
            mesh.update();
        }

        

        // --- Directed interpolation of primitive fields onto faces

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        //volScalarField rPsi("rPsi", 1.0/psi);
        //surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        //surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

        surfaceScalarField p_pos(interpolate(p, neg));
        surfaceScalarField p_neg(interpolate(p, neg));

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        // Note: extracted out the orientation so becomes unoriented
        phiv_pos.setOriented(false);
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());
        phiv_neg.setOriented(false);

        // Make fluxes relative to mesh-motion
        if (mesh.moving())
        {
            surfaceScalarField meshPhi(mesh.phi());
            meshPhi.setOriented(false);
            phiv_pos -= meshPhi;
            phiv_neg -= meshPhi;
        }

        //- Sound velocity
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, pos, T.name())*mesh.magSf()
        );

        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, neg, T.name())*mesh.magSf()
        );

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );

        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos("a_pos", ap/(ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf("aSf", am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "centralCourantNo.H"

        if (LTS)
        {
            #include "setRDeltaT.H"

            ++runTime;
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        surfaceVectorField phiU(aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg);
        // Note: reassembled orientation from the pos and neg parts so becomes
        // oriented
        phiU.setOriented(true);

        surfaceVectorField phiUp(phiU + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf());

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg
        );

        // Make flux for pressure-work absolute
        if (mesh.moving())
        {
            surfaceScalarField meshPhi(mesh.phi());
            meshPhi.setOriented(false);
            phiEp += meshPhi*(a_pos*p_pos + a_neg*p_neg);
        }

        //volScalarField muEff("muEff", turbulence->muEff());
        //volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        // Info << " --------------------------- Fluxes --------------------------------- " << endl;

        // Info << "phi \n " << phi << endl;
        // Info << "phiU \n " << phiU << endl;
        // Info << "phiUp \n " << phiUp << endl; 
        // Info << "phiEp \n " << phiEp << endl; 


        // Info << "-------------------------- Conserved quantities --------------------- " << endl;

        // Info << "rho \n " << rho << endl;
        // Info << "rhoU \n " << rhoU << endl;
        // Info << "rhoE \n " << rhoE << endl;

        // Info << " -------------------------------------------------------------------- " << endl;

        if (hydrodynamics)
        {
            

            // --- Solve density

            Info << "Solving for density \n" << endl;
            solve(fvm::ddt(rho) + fvc::div(phi));
            rho.correctBoundaryConditions();


            // --- Solve momentum
            Info << "Solving for momentum \n " << endl;
            solve(fvm::ddt(rhoU) + fvc::div(phiUp));

            U.ref() =
                rhoU()
               /rho();
            U.correctBoundaryConditions();
            rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

            // Turbulence models depend on the thermophysicalModels library.
            // Turned off for now.

            // if (!inviscid)
            // {
            //     solve
            //     (
            //         fvm::ddt(rho, U) - fvc::ddt(rho, U)
            //       - fvm::laplacian(muEff, U)
            //       - fvc::div(tauMC)
            //     );
            //     rhoU = rho*U;
            // }

            // --- Solve energy
            // surfaceScalarField sigmaDotU
            // (
            //     "sigmaDotU",
            //     (
            //         fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
            //       + fvc::dotInterpolate(mesh.Sf(), tauMC)
            //     )
            //   & (a_pos*U_pos + a_neg*U_neg)
            // );



            // For some reason the radiation source cannot be added here because there would be an incompatibility
            // This is the error that I get (for future reference):

            // Incompatible fields for operation
            // [rhoE] - [thermo:e]

            Info << "Solving for energy " << endl;
            solve
            (
                fvm::ddt(rhoE)
              + fvc::div(phiEp)
              //- fvc::div(sigmaDotU)
            );

            e = rhoE/rho - 0.5*magSqr(U);

            // Info << "--------------------Solution------------ " << endl;
            // Info << "rho \n " << rho << endl;
            // Info << "rhoU \n " << rhoU << endl;
            // Info << "rhoE \n " << rhoE << endl;
            // Info << "-----------------------------------------" << endl;

            //e.correctBoundaryConditions();

            Info << "Correcting thermo quantities " << endl;

            // Calculate new temperature, pressure, thermophysical properties
            thermo.correctFromRhoE();
            rhoE.boundaryFieldRef() ==
                rho.boundaryField()*
                (
                    e.boundaryField() + 0.5*magSqr(U.boundaryField())
                );
        }

        if (radiation->radiation())
        {
            Info << "Radiation is on \n" << endl;

            Info << "Radiation correct \n " << endl;
            radiation->correct();

            Info << "Correcting internal energy to account for radiation \n " << endl;
            solve
            (
                fvm::ddt(rho, e) 
                //- fvc::ddt(rho, e) I don't know what fvc is supposed to do here
                - radiation->Sh(thermo, e)
              //- fvm::laplacian(turbulence->alphaEff(), e)
            );

            thermo.correctTemperature();
        }

        if (hydrodynamics)
        {
            thermo.correctFromRhoE();
            rhoE = rho*(e + 0.5*magSqr(U));    
        }
       
        /*if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
              - fvm::laplacian(turbulence->alphaEff(), e)
            );
            thermo.correct();
            rhoE = rho*(e + 0.5*magSqr(U));
        }*/

        //Info << "thermo quantities after correction \n " << endl;

        // Info << "p " << p << endl;
        // Info << "rho " << rho << endl;
        // Info << "T " << T << endl;
        // Info << "c " << c << endl;


        // Pressure correction of internal field and boundary is performed in correctFromRhoE()

        /*p.ref() =
            rho()
           /psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();*/

        

        //turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);


      
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
