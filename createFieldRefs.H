volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const volScalarField& c = thermo.speedOfSound();
const volScalarField& mu = thermo.mu();


Info << "p[1] = " << p[1] << endl;
Info << "T[1] = " << T[1] << endl;
Info << "c[1] = " << c[1] << endl;


bool inviscid(true);
if (max(mu.primitiveField()) > 0.0)
{
    inviscid = false;
}
