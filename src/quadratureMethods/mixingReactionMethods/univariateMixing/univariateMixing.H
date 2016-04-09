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
    Foam::univariateMixing

Description
    Solve a univariate mixing equation with the extended quadrature
    method of moments.

SourceFiles
    univariateMixing.C

\*---------------------------------------------------------------------------*/

#ifndef univariateMixing_H
#define univariateMixing_H

#include "mixingKernel.H"
#include "reactionKernel.H"
#include "diffusionModel.H"
#include "univariateMomentSet.H"
#include "fvCFD.H"
#include "IOdictionary.H"
#include "PtrList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "quadratureNodes.H"
#include "moments.H"
#include "momentFieldSet.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class univariateMixing Declaration
\*---------------------------------------------------------------------------*/

class univariateMixing
:
    public IOdictionary
{
    // Private data
  
        //- Name of the quadrature approximation
        const word name_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- AutoPtr to the PtrList of quadrature nodes
        autoPtr<PtrList<basicVolVectorNode> > nodes_;

        //- List of moments
        momentFieldSet<basicVolVectorMoment, basicVolVectorNode> moments_;
        
        //- Number of primary quadrature nodes
        label nPrimaryNodes_;
        
        //- List of interpolated nodes (neighbour)
        autoPtr<PtrList<basicSurfaceVectorNode> > nodesNei_;
        
        //- List of interpolated nodes (owner)
        autoPtr<PtrList<basicSurfaceVectorNode> > nodesOwn_;
        
        //- Dimensionality of the distribution function
        const label nDimensions_;

        //- Number of moments
        const label nMoments_;
        
        //- List of interpolated moments (neighbour)
        momentFieldSet<basicSurfaceVectorMoment, basicSurfaceVectorNode> momentsNei_;

        //- List of interpolated moments (owner)
        momentFieldSet<basicSurfaceVectorMoment, basicSurfaceVectorNode> momentsOwn_;
        
        //- Support of the distribution function
        word support_;
        
        //- Switch indicating the status of the mixing model
        Switch mixing_;
        
        // Switch indicating the status of the reaction model
        Switch reaction_;
        
        //- Concentration A
        volScalarField Ca_;
        
        //- Concentration B
        volScalarField Cb_;
        
        //- Concentration R
        volScalarField Cr_;
        
        //- Concentration S
        volScalarField Cs_;
        
        //- mixing kernel model
        autoPtr
        <
            Foam::mixingSubModels::mixingKernel
        > mixingKernel_;
        
        //- reaction kernel model
        autoPtr
        <
            Foam::reactionSubModels::reactionKernel
        > reactionKernel_;
        
        //- Diffusion model
        autoPtr
        <
            Foam::populationBalanceSubModels::diffusionModel
        > diffusionModel_;
        
    // Private member functions

        //- Interpolate nodes
        void interpolateNodes();
        
        //- Update quadrature nodes on fixed-value boundaries
        void updateBoundaryQuadrature
        (
            const momentFieldSet
            <
                basicVolVectorMoment,
                basicVolVectorNode
            >& moments,
            PtrList<basicVolVectorNode>& nodes
        );
        
        //- Recalculate the quadrature nodes from the moments
        void updateQuadrature
        (
            const momentFieldSet
            <
                basicVolVectorMoment,
                basicVolVectorNode
            >& moments,
            PtrList<basicVolVectorNode>& nodes
        );
        
        //- Update fields required to compute kinetic fluxes
        void updateAdvection
        (
            surfaceScalarField& phiOwn, 
            surfaceScalarField& phiNei
        );
        
        //- Calculate explicit advection term 
        tmp<volScalarField> advectMoment
        (
            const basicVolVectorMoment& moment,
            const surfaceScalarField& phiOwn,
            const surfaceScalarField& phiNei
        );
        
        //- Calculate mixing source terms
        tmp<fvScalarMatrix> mixingSource
        (
            const basicVolVectorMoment& moment,
            const momentFieldSet<basicVolVectorMoment,basicVolVectorNode>& moments
        );
        
        //- Calculate reaction source terms
        tmp<volScalarField> reactionSource
        (
            const basicVolVectorMoment& moment,
            const PtrList<basicVolVectorNode>& nodes
        );

protected:
        
    // Protected data

        //- Fluid phase cell-centered velocity
        const volVectorField& U_;
        
        //- Fluid-phase face-normal velocity
        const surfaceScalarField& phi_;
        
public:

    //- Runtime type information
    TypeName("univariateMixing");


    // Constructors

        //- Construct from components
        univariateMixing
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh,
            const volVectorField& U_,
            const surfaceScalarField& phi_,
            const word support
        );


        //- Select null constructed
        static autoPtr<univariateMixing> New
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh,
            const volVectorField& U,
            const surfaceScalarField& phi_,
            const word support
        );
        
    //- Destructor
    ~univariateMixing();


    // Member Functions

        //- Constant Access to Concentration A
        inline const volScalarField& Ca() const;
        
        //- Constant Access to Concentration B
        inline const volScalarField& Cb() const;
        
        //- Constant Access to Concentration R
        inline const volScalarField& Cr() const;
        
        //- Constant Access to Concentration S
        inline const volScalarField& Cs() const;
        
        //- Constant Access to nodes
        inline const PtrList<basicVolVectorNode>& nodes() const;
        
        //- Solve univariate mixing equation
        void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateMixingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
