#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;


/*autoPtr<rhoThermo> pThermo
(
    rhoThermo::New(mesh)
);

rhoThermo& thermo = pThermo();*/

/*autoPtr<myThermo> pThermo
(
    myThermo::New(mesh)
);

myThermo& thermo = pThermo();*/

autoPtr<myRhoThermo> pThermo
(
    myRhoThermo::New(mesh)
);

myRhoThermo& thermo = pThermo();


volScalarField& e = thermo.eRef();

volScalarField& rho = thermo.getRho();

Info<< "Reading field U\n" << endl;


volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// volVectorField internalEnergy
// (
//     IOobject
//     (
//         "internalEnergy",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     e
// );

// Info << "field named e created " << endl;
// Info << internalEnergy << endl;

// volScalarField rho
// (
//     IOobject
//     (
//         "rho",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     thermo.rho()
// );

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

surfaceScalarField phi("phi", fvc::flux(rhoU));


#include "createRadiationModel.H"


Info<< "Turbulence model turned off for now\n" << endl;
// autoPtr<compressible::turbulenceModel> turbulence
// (
//     compressible::turbulenceModel::New
//     (
//         rho,
//         U,
//         phi,
//         thermo
//     )
// );
