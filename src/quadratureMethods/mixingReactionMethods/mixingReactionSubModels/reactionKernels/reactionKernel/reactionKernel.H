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
    Foam::mixingReactionSubModels::reactionKernel

Description
    Abstract class for reaction kernels.
    
SourceFiles
    reactionKernel.C

\*---------------------------------------------------------------------------*/

#ifndef reactionKernel_H
#define reactionKernel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "fvCFD.H"
#include "quadratureNodes.H"
#include "moments.H"
#include "momentFieldSet.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionSubModels
{

/*---------------------------------------------------------------------------*\
                    Class reactionKernel Declaration
\*---------------------------------------------------------------------------*/

class reactionKernel
{   
    // Private member functions

        //- Disallow default bitwise copy construct
        reactionKernel(const reactionKernel&);

        //- Disallow default bitwise assignment
        void operator=(const reactionKernel&);


protected:

    // Protected data

        const dictionary& dict_;
        
        //- Initial concentration 1
        dimensionedScalar Ca0_;
        
        //- Initial concentration 2
        dimensionedScalar Cb0_;


public:

    //- Runtime type information
    TypeName("reactionKernel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        reactionKernel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        reactionKernel(const dictionary& dict);


    // Selectors

        static autoPtr<reactionKernel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~reactionKernel();


    // Member Functions

        //- Constant access to Ca0_
        virtual dimensionedScalar Ca0() const;
        
        //- constant access to Cb0_
        virtual dimensionedScalar Cb0() const;
        
        //- reaction kernel 1
        virtual tmp<volScalarField> R1
        (
            const volScalarField& abscissa,
            const volScalarField& Y1,
            const volScalarField& Y2
        ) const = 0;
        
        //- reaction kernel 2
        virtual tmp<volScalarField> R2
        (
            const volScalarField& abscissa,
            const volScalarField& Y1,
            const volScalarField& Y2
        ) const = 0;
        
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
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
