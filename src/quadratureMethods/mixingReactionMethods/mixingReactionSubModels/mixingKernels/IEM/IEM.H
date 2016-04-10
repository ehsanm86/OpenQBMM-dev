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
    Foam::mixingModels::mixingKernels::IEM

Description
    Interaction by exchange with the mean(IEM) kernel. It assumes a linear
    relaxation of the scalar towards its mean vlaue.
    
    \verbatim
        "Computational methods for turbulent reacting flows"
        Rodney O. Fox
        Cambridge university press, 2003
    \endverbatim

SourceFiles
    IEM.C

\*---------------------------------------------------------------------------*/

#ifndef IEM_H
#define IEM_H

#include "mixingKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
namespace mixingKernels
{

/*---------------------------------------------------------------------------*\
                          Class IEM Declaration
\*---------------------------------------------------------------------------*/

class IEM
:
    public mixingKernel
{
    //- Private Data
        //- Cphi constant 
        dimensionedScalar Cphi_;
        
public: 
        //- Runtime type information
        TypeName("IEM");


    // Constructors

        //- Construct from components
        IEM(const dictionary& dict);


    //- Destructor
    virtual ~IEM();


    // Member Functions

        //- Mixing kernel
        virtual tmp<fvScalarMatrix> Km
        (
            const basicVolUnivariateMoment& moment,
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
