/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
     \\/     M anipulation  |
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

Class
    Foam::univariateMixingI

Description
    Solve a univariate mixing equation with the extended quadrature
    method of moments.

SourceFiles
    univariateMixing.C

\*---------------------------------------------------------------------------*/



// * * * * * * * * * * * * Memeber Functon Declaration * * * * * * * * * * * //

//- Constant Access to Concentration A
const Foam::volScalarField& 
Foam::univariateMixing::Ca() const
{
    return Ca_;
}

//- Constant Access to Concentration B
const Foam::volScalarField& 
Foam::univariateMixing::Cb() const
{
    return Cb_;
}

//- Constant Access to Concentration R
const Foam::volScalarField& 
Foam::univariateMixing::Cr() const
{
    return Cr_;
}

//- Constant Access to Concentration S
const Foam::volScalarField& 
Foam::univariateMixing::Cs() const
{
    return Cs_;
}

//- Constant Access to nodes
const Foam::PtrList<Foam::basicVolScalarNode>& 
Foam::univariateMixing::nodes() const
{
    return nodes_();
}
// ************************************************************************* //
