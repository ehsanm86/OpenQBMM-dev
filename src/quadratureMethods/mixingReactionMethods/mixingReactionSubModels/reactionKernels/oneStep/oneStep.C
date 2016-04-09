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

#include "oneStep.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionSubModels
{
namespace reactionKernels
{
    defineTypeNameAndDebug(oneStep, 0);

    addToRunTimeSelectionTable
    (
        reactionKernel,
        oneStep,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionSubModels::reactionKernels::oneStep
::oneStep
(
    const dictionary& dict
)
:
    reactionKernel(dict),
    k1_(dict.lookup("k1"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionSubModels::reactionKernels::oneStep
::~oneStep()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::reactionSubModels::reactionKernels::oneStep::R1
(
    const volScalarField& abscissa,
    const volScalarField& Y1,
    const volScalarField& Y2
) const
{   
    dimensionedScalar mixFraction0 = Ca0_/(Ca0_ + Cb0_);

    return mixFraction0*k1_*Cb0_*((1.0 - abscissa)/(1.0 - mixFraction0) - Y1)
        *(abscissa/mixFraction0 - Y1);
}

Foam::tmp<Foam::volScalarField>
Foam::reactionSubModels::reactionKernels::oneStep::R2
(
    const volScalarField& abscissa,
    const volScalarField& Y1,
    const volScalarField& Y2
) const
{   
    tmp<volScalarField> r2
    (
        new volScalarField
        (
            IOobject
            (
                "r2",
                Y1.mesh().time().timeName(),
                Y1.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            Y1.mesh(),
            dimensionedScalar("r2", inv(dimTime), 0.0)
        )
    );
    
    return r2;
}

//- Update concentrations
void Foam::reactionSubModels::reactionKernels::oneStep::updateConcentrations
(
    volScalarField& Ca,
    volScalarField& Cb,
    volScalarField& Cr,
    volScalarField& Cs,
    const momentFieldSet
    <
        basicVolVectorMoment,
        basicVolVectorNode
    >& moments
)
{
    const basicVolVectorMoment& mixF(moments(1,0,0));
    const basicVolVectorMoment& Y1(moments(0,1,0));
    
    dimensionedScalar mixF0 = Ca0_/(Ca0_ + Cb0_);
    
    Ca == Ca0_*(1.0 - mixF - (1.0 - mixF0)*Y1);
    Cb == Cb0_*(mixF - mixF0*Y1);
    Cr == Cb0_*mixF0*Y1;
    
    return;
}
// ************************************************************************* //
