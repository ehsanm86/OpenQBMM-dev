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

#include "competitiveConsecutive.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionSubModels
{
namespace reactionKernels
{
    defineTypeNameAndDebug(competitiveConsecutive, 0);

    addToRunTimeSelectionTable
    (
        reactionKernel,
        competitiveConsecutive,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionSubModels::reactionKernels::competitiveConsecutive
::competitiveConsecutive
(
    const dictionary& dict
)
:
    reactionKernel(dict),
    k1_(dict.lookup("k1")),
    k2_(dict.lookup("k2"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionSubModels::reactionKernels::competitiveConsecutive
::~competitiveConsecutive()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::reactionSubModels::reactionKernels::competitiveConsecutive::R1
(
    const volScalarField& abscissa,
    const volScalarField& Y1,
    const volScalarField& Y2
) const
{   
    dimensionedScalar mixFraction0 = Ca0_/(Ca0_ + Cb0_);

    tmp<volScalarField> r1 = 
        mixFraction0*k1_*Cb0_*((1.0 - abscissa)/(1.0 - mixFraction0) - Y1)
        *(abscissa/mixFraction0 - Y1 - Y2);
        
//     Info << "abscissa=" << abscissa << endl;
//     Info << "Y1=" << Y1 << endl;
//     Info << "Y2=" << Y2 << endl;
//     Info << "k1_" << k1_ << endl;
//     Info << "abscissa/mixFraction0" << abscissa/mixFraction0 << endl;

    return r1;
}

Foam::tmp<Foam::volScalarField>
Foam::reactionSubModels::reactionKernels::competitiveConsecutive::R2
(
    const volScalarField& abscissa,
    const volScalarField& Y1,
    const volScalarField& Y2
) const
{
    dimensionedScalar mixFraction0 = Ca0_/(Ca0_ + Cb0_);

    tmp<volScalarField> r2 = mixFraction0*k2_*Cb0_*(Y1 - Y2)
        *(abscissa/mixFraction0 - Y1 - Y2);

    return r2;
}

void Foam::reactionSubModels::reactionKernels
::competitiveConsecutive::updateConcentrations
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
    const basicVolVectorMoment& Y2(moments(0,0,1));
    
    dimensionedScalar mixF0 = Ca0_/(Ca0_ + Cb0_);
    
//     Ca == Ca0_*(moments(0,0,0) - mixF - (1.0 - mixF0)*Y1)/moments(0,0,0);
//     Cb == Cb0_*(mixF - mixF0*(Y1 + Y2))/moments(0,0,0);
//     Cr == Cb0_*mixF0*(Y1 - Y2)/moments(0,0,0);
//     Cs == Cb0_*mixF0*Y2/moments(0,0,0);
    
    Ca == Ca0_*(1.0 - mixF - (1.0 - mixF0)*Y1);
    Cb == Cb0_*(mixF - mixF0*(Y1 + Y2));
    Cr == Cb0_*mixF0*(Y1 - Y2);
    Cs == Cb0_*mixF0*Y2;
    
    return;
}
// ************************************************************************* //
