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

\*---------------------------------------------------------------------------*/

#include "reactionKernel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionSubModels
{
    defineTypeNameAndDebug(reactionKernel, 0);

    defineRunTimeSelectionTable(reactionKernel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionSubModels::reactionKernel::reactionKernel
(
    const dictionary& dict
)
:
    dict_(dict),
    Ca0_
    (
        dict.lookupOrDefault
        (
            "Ca0", 
            dimensionedScalar("zero", dimless, 0.0)  
        )
    ),
    Cb0_
    (
        dict.lookupOrDefault
        (
            "Cb0",
            dimensionedScalar("zero", dimless, 0.0)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionSubModels::reactionKernel::~reactionKernel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::reactionSubModels::reactionKernel::Ca0() const
{
    return Ca0_;
}

Foam::dimensionedScalar Foam::reactionSubModels::reactionKernel::Cb0() const
{
    return Cb0_;
}
// ************************************************************************* //
