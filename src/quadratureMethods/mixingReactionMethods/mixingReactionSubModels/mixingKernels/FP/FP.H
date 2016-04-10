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
    Foam::mixingModels::mixingKernels::FP

Description
    Fokker-Planck(FP) kernel. It includes the effect of differential diffusion.
    
    \verbatim
        "Computational methods for turbulent reacting flows"
        Rodney O. Fox
        Cambridge university press, 2003
    \endverbatim

SourceFiles
    FP.C

\*---------------------------------------------------------------------------*/

#ifndef FP_H
#define FP_H

#include "mixingKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
namespace mixingKernels
{

/*---------------------------------------------------------------------------*\
                          Class FP Declaration
\*---------------------------------------------------------------------------*/

class FP
:
    public mixingKernel
{
  //- Private Data
        //- Cphi constant 
        dimensionedScalar Cphi_;
public:

        //- Runtime type information
        TypeName("FP");


    // Constructors

        //- Construct from components
        FP(const dictionary& dict);


    //- Destructor
    virtual ~FP();


    // Member Functions

        //- Mixing kernel
        virtual tmp<fvScalarMatrix> Km
        (
            const labelList& order,
            const momentFieldSet<basicVolUnivariateMoment,basicVolScalarNode>& moments
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixingKernels
} // End namespace mixingSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
