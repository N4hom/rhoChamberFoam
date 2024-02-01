/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "lookupTables2D.H"

// * * * * * * * * * * * * * * Scalar Functions * * * * * * * * * * * * * * //

template<>
Foam::labelList Foam::lookupTable2D<Foam::scalar>::boundi
(
    const scalar& f
) const
{
    label& i = ij_.x();
    const label j = ij_.y();
    if (f < data_[0][j])
    {
        return labelList(1, 0);
    }

    DynamicList<label> Is(data_.size());
    for (i = 0; i < data_.size(); i++)
    {
        if (f > data_[i][j] && f < data_[i][j+1])
        {
            Is.append(i);
        }
    }
    if (!Is.size())
    {
        return labelList(1, data_.size() - 2);
    }
    if (Is.size() == 1)
    {
        return std::move(Is);
    }
    DynamicList<label> newIs(Is.size());
    for (i = 0; i < Is.size() - 1; i++)
    {
        if (Is[i]+1 == Is[i+1])
        {
            newIs.append(Is[i]);
        }
    }
    return std::move(newIs);
}


template<>
Foam::labelList Foam::lookupTable2D<Foam::scalar>::boundj
(
    const scalar& f
) const
{
    const label i = ij_.x();
    label& j = ij_.y();
    if (data_[i][0] > f)
    {
        return labelList(1, 0);
    }

    DynamicList<label> Js(data_[i].size());
    for (j = 0; j < data_[i].size(); j++)
    {
        if (f > data_[i][j] && f < data_[i+1][j])
        {
            Js.append(j);
        }
    }
    if (!Js.size())  // if Js.size() == 0
    {
        return labelList(1, data_[i].size() - 2);
    }
    if (Js.size() == 1)
    {
        return std::move(Js);
    }
    DynamicList<label> newJs(Js.size());
    for (j = 0; j < Js.size() - 1; j++)
    {
        if (Js[j]+1 == Js[j+1])
        {
            newJs.append(Js[j]);
        }
    }
    return std::move(newJs);
}


template<>
Foam::scalar Foam::lookupTable2D<Foam::scalar>::reverseLookupX
(
    const scalar& fin,
    const scalar y
) const
{
    ////Info << "fin " <<  fin << endl;
    ////Info << "y " <<  y << endl;
    scalar f(mod_()(fin));
    ij_.y() = yIndexing_->findIndex(modY_()(y));
    labelList Is(boundi(f));

    const label j = ij_.y();
    label& i = ij_.x();
    scalar fy = interpolationWeight1D::linearWeight
    (
        modY_()(y),
        yModValues_[j],
        yModValues_[j+1]
    );
    scalar fx = 1.0;

    if (Is.size() == 1)
    {
        i = Is[0];
        const scalar mm(data_[i][j]);
        const scalar pm(data_[i+1][j]);
        const scalar mp(data_[i][j+1]);
        const scalar pp(data_[i+1][j+1]);

        ////Info << (f + fy*(mm - mp) - mm)/(fy*(mm - mp - pm + pp) - mm + pm) << endl;
        fx =
            (f + fy*(mm - mp) - mm)
           /(fy*(mm - mp - pm + pp) - mm + pm);
        ////Info << "i = " << i << endl;
        ////Info << "fx = " << fx << endl;
        ////Info << "xModValues_ = " << xModValues_ << endl;
        ////Info << "modX_->inv(getValue(i, fx, xModValues_)) " << modX_->inv(getValue(i, fx, xModValues_)) << endl;

        return modX_->inv(getValue(i, fx, xModValues_));
    }
    

    //- If multiple indicies meet criteria, check for closest
    Field<scalar> xTrys(Is.size());
    Field<scalar> errors(Is.size(), GREAT);
    forAll(Is, I)
    {
        i = Is[I];
        const scalar mm(data_[i][j]);
        const scalar pm(data_[i+1][j]);
        const scalar mp(data_[i][j+1]);
        const scalar pp(data_[i+1][j+1]);
        fx =
            (f + fy*(mm - mp) - mm)
           /(fy*(mm - mp - pm + pp) - mm + pm);
        xTrys[I] = modX_->inv(getValue(i, fx, xModValues_));
        errors[I] = mag(fin - lookup(xTrys[I], y));
    }
    i = findMin(errors);


    ////Info << "xTrys[j] " << xTrys[i] << endl;
    return xTrys[i];
}


template<>
Foam::scalar Foam::lookupTable2D<Foam::scalar>::reverseLookupY
(
    const scalar& fin,
    const scalar x
) const
{
    ////Info << "reverseLookupY " << endl;

    ////Info << "x = " << x << endl;
    ////Info << "fin = " << fin << endl;
    scalar f(mod_()(fin));
    ////Info << "fmod = " << f << endl;

    ////Info << "f " << f << endl;
    ////Info << "fin " << fin << endl;

    ij_.x() = xIndexing_->findIndex(modX_()(x));
    ////Info << "x index associated with x = " << x << ": " << ij_.x() << endl;
    labelList Js(boundj(f));

    ////Info << "labelList Js that bounds fmod " << Js << endl;

    const label i = ij_.x();
    ////Info << "label i " << i << endl;

    ////Info << "i " << i << endl;
    label& j = ij_.y();

    ////Info << "j " << j << endl;
    scalar fx = interpolationWeight1D::linearWeight
    (
        modX_()(x),
        xModValues_[i],
        xModValues_[i+1]
    );
    ////Info << "modX_()(x) = " << modX_()(x) << endl;
    ////Info << "xModValues_ = " << xModValues_ << endl;
    ////Info << "xModValues_[i]= " << xModValues_[i] << endl;
    
    ////Info << "yModValues_ " << yModValues_ << endl;
    
    ////Info << "fx " << fx << endl;


    ////Info << "fx " << fx << endl;
    scalar fy = 1.0;

    if (Js.size() == 1)
    {
      //  ////Info << "Jd.size() == 1 " << endl;
        j = Js[0];
        const scalar mm(data_[i][j]);
        const scalar pm(data_[i+1][j]);
        const scalar mp(data_[i][j+1]);
        const scalar pp(data_[i+1][j+1]);
        
        // fy is used to interpolate the y values of the table (e.g. T if f = f(p,T))
        ////Info << "fy " << fy << endl;
        fy =
            (f + fx*(mm  - pm) - mm)
           /(fx*(mm - pm - mp + pp) - mm + mp);
        ////Info << "fy " << fy << endl;

        // getValue gives the y value at which the function f(x,y) will be evaluated
        ////Info << getValue(j, fy, yModValues_) << endl;

        // the function inv essentially performs the inverse operation defined by the modifier (e.g. exp(getValue) if the modifier is ln) 
        ////Info << modY_->inv(getValue(j, fy, yModValues_)) << endl;

        return modY_->inv(getValue(j, fy, yModValues_));
    }

    //- If multiple indicies meet criteria, check for closest
    Field<scalar> yTrys(Js.size());
    ////Info <<"yTrys " << yTrys  << endl;
    Field<scalar> errors(Js.size(), GREAT);
    forAll(Js, J)
    {
      //  ////Info << "Js index " << J << endl;
        j = Js[J];
        const scalar mm(data_[i][j]);
      //  ////Info << "mm " << mm << endl;
        const scalar pm(data_[i+1][j]);
      //  ////Info << "pm " << pm << endl;
        const scalar mp(data_[i][j+1]);
      //  ////Info << "mp " << mp << endl;
        const scalar pp(data_[i+1][j+1]);
      //  ////Info << "pp " << pp << endl;

        fy =
            (f + fx*(mm  - pm) - mm)
           /(fx*(mm - pm - mp + pp) - mm + mp);

      //  ////Info << "fy " << (f + fx*(mm  - pm) - mm)
      //     /(fx*(mm - pm - mp + pp) - mm + mp) << endl;
        yTrys[J] = modY_->inv(getValue(j, fy, yModValues_));
      //  ////Info << "yTrys[J] " << yTrys[J] << endl;
        errors[J] = mag(fin - lookup(x, yTrys[J]));
    }
    j = findMin(errors);

    ////Info << "yTrys[j] " << yTrys << endl;

    return yTrys[j];
}


// ************************************************************************* //
