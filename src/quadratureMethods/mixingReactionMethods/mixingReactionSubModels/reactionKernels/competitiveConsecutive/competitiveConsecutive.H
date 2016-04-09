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
    Foam::mixingReactionModels::reactionKernels::competitiveConsecutive

Description
    competitiveConsecutive reaction kernel.
    
    Reference
    \verbatim
        
    \endverbatim    

SourceFiles
    competitiveConsecutive.C

\*---------------------------------------------------------------------------*/

#ifndef competitiveConsecutive_H
#define competitiveConsecutive_H

#include "reactionKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionSubModels
{
namespace reactionKernels
{

/*---------------------------------------------------------------------------*\
                    Class competitiveConsecutive Declaration
\*---------------------------------------------------------------------------*/

class competitiveConsecutive
:
    public reactionKernel
{
private:
    //- Private Data
    
        //- Reaction rate constant 1
        dimensionedScalar k1_;
        
        //- Reaction rate constant 2
        dimensionedScalar k2_;
        
public:

        //- Runtime type information
        TypeName("competitiveConsecutive");


    // Constructors

        //- Construct from components
        competitiveConsecutive(const dictionary& dict);


    //- Destructor
    virtual ~competitiveConsecutive();


    // Member Functions

        //- reaction kernel 1
        virtual tmp<volScalarField> R1
        (
            const volScalarField& abscissa,
            const volScalarField& Y1,
            const volScalarField& Y2
        ) const;
        
        //- reaction kernel 2
        virtual tmp<volScalarField> R2
        (
            const volScalarField& abscissa,
            const volScalarField& Y1,
            const volScalarField& Y2
        ) const;
        
        //- Update concentrations
        virtual void updateConcentrations
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
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionKernels
} // End namespace reactionSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
