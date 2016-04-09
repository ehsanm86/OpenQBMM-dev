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

#include "IEM.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
namespace mixingKernels
{
    defineTypeNameAndDebug(IEM, 0);

    addToRunTimeSelectionTable
    (
        mixingKernel,
        IEM,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernels::IEM::IEM
(
    const dictionary& dict
)
:
    mixingKernel(dict),
    Cphi_(dict.lookup("Cphi"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernels::IEM::~IEM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::mixingSubModels::mixingKernels::IEM::Km
(
    const labelList& order,
    const momentFieldSet<basicVolVectorMoment,basicVolVectorNode>& moments
)
{
    label a = order[0];
    label b = order[1];
    label c = order[2];
    
    //incompressibleSolver
//     const incompressible::turbulenceModel& fluidTurb =
//         moments(0,0,0).mesh().lookupObject<incompressible::turbulenceModel>
//         (
//             turbulenceModel::propertiesName
//         );
    // compressibleSolver
    const compressible::turbulenceModel& fluidTurb =
        moments(0,0,0).mesh().lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );
    
    tmp<fvScalarMatrix> IEMK
    (
        new fvScalarMatrix
        (
            moments(a,b,c),
            moments(a,b,c).dimensions()*dimVol/dimTime
        )
    );

    if (a == 0)
    {
        return IEMK;
    }
    else
    {
        IEMK.ref() += Cphi_*scalar(a)*fluidTurb.epsilon()/fluidTurb.k()
            *(moments(a-1,b,c)*moments(1,0,0))
            - fvm::SuSp(Cphi_*scalar(a)*fluidTurb.epsilon()
            /fluidTurb.k(), moments(a,b,c));
        
        return IEMK;
    }
}

// ************************************************************************* //
