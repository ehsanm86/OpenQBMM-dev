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

#include "univariateMixing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(univariateMixing, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateMixing::univariateMixing
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word support
)
:
    IOdictionary
    (
        Foam::IOobject
        (
            "quadratureProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    mesh_(mesh),
    nodes_(),
    moments_(name_, *this, mesh_, nodes_),
    nPrimaryNodes_(0),
    nodesNei_(),
    nodesOwn_(),
    nDimensions_(1),                 
    nMoments_(moments_.size()),
    momentsNei_
    (
        name_, nMoments_, nodesNei_, nDimensions_, moments_.momentMap()
    ),
    momentsOwn_
    (
        name_, nMoments_, nodesOwn_, nDimensions_, moments_.momentMap()
    ),
    support_(support),
    mixing_(dict.lookup("mixing")),
    Ca_
    (
        Foam::IOobject
        (
            "Ca",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("Ca", dimless, 0.0)
    ),
    Cb_
    (
        Foam::IOobject
        (
            "Cb",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("Cb", dimless, 0.0)
    ),
    Cr_
    (
        Foam::IOobject
        (
            "Cr",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("Cr", dimless, 0.0)
    ),
    Cs_
    (
        Foam::IOobject
        (
            "Cs",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("Cs", dimless, 0.0)
    ),
    mixingKernel_
    (
        Foam::mixingSubModels::mixingKernel::New
        (
            dict.subDict("mixingKernel")
        )
    ),
    diffusionModel_
    (
        Foam::populationBalanceSubModels::diffusionModel::New
        (
            dict.subDict("diffusionModel")
        )
    ),
    U_(U),
    phi_(phi)
{
    // Allocating nodes
    nodes_ = autoPtr<PtrList<basicVolScalarNode> >
    (
        new PtrList<basicVolScalarNode>
        (
            lookup("nodes"), 
            Foam::basicVolScalarNode::iNew
            (
                name_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions(),
                moments_[0].boundaryField().types()
            )
        )
    );
   
    nPrimaryNodes_ = nodes_().size();
    
    nodesNei_ = autoPtr<PtrList<basicSurfaceScalarNode> >
    (
        new PtrList<basicSurfaceScalarNode>(nPrimaryNodes_)
    );
    
    nodesOwn_ = autoPtr<PtrList<basicSurfaceScalarNode> >
    (
        new PtrList<basicSurfaceScalarNode>(nPrimaryNodes_)
    );
    
    PtrList<basicVolScalarNode>& nodes = nodes_();
    PtrList<basicSurfaceScalarNode>& nodesNei = nodesNei_();
    PtrList<basicSurfaceScalarNode>& nodesOwn = nodesOwn_();

    // Populating interpolated nodes.
    forAll(nodes, pNodeI)
    {
        basicVolScalarNode& node(nodes[pNodeI]);
              
        nodesNei.set
        (
            pNodeI,
            new basicSurfaceScalarNode
            (
                node.name() + "Nei",
                name_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions()
            )
        );

        nodesOwn.set
        (
            pNodeI,
            new basicSurfaceScalarNode
            (
                node.name() + "Own",
                name_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions()
            )
        );
    }
    
    // Setting face values of moments
    forAll(momentsNei_, mI)
    {
        momentsNei_.set
        (
            mI,
            new Foam::basicSurfaceUnivariateMoment
            (
                name_,
                moments_[mI].cmptOrders(),
                nodesNei_,
                fvc::interpolate(moments_[mI])
            )
        );
        
        momentsOwn_.set
        (
            mI,
            new Foam::basicSurfaceUnivariateMoment
            (
                name_,
                moments_[mI].cmptOrders(),
                nodesOwn_,
                fvc::interpolate(moments_[mI])
            )
        );
    }
    
    updateQuadrature();
    interpolateNodes();
    updateBoundaryQuadrature();
    momentsNei_.update();
    momentsOwn_.update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateMixing::~univariateMixing()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::univariateMixing::interpolateNodes()
{
    surfaceScalarField nei
    (
        IOobject
        (
            "nei",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("nei", dimless, -1.0)
    );

    surfaceScalarField own
    (
        IOobject
        (
            "own",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("own", dimless, 1.0)
    );
    
    const PtrList<basicVolScalarNode>& nodes = nodes_();
    PtrList<basicSurfaceScalarNode>& nodesOwn = nodesOwn_();
    PtrList<basicSurfaceScalarNode>& nodesNei = nodesNei_();
        
    forAll(nodes, pNodeI)
    {
        const basicVolScalarNode& node(nodes[pNodeI]);
        basicSurfaceScalarNode& nodeOwn(nodesOwn[pNodeI]);
        basicSurfaceScalarNode& nodeNei(nodesNei[pNodeI]);
               
        nodeOwn.primaryWeight() = 
            fvc::interpolate
            (
                node.primaryWeight(),
                own,
                "reconstruct(primaryWeight)"
            );
        
        nodeOwn.primaryAbscissa() = 
            fvc::interpolate
            (
                node.primaryAbscissa(), 
                own, 
                "reconstruct(primaryAbscissa)"
            );

        nodeNei.primaryWeight() = 
            fvc::interpolate
            (
                node.primaryWeight(),
                nei,
                "reconstruct(primaryWeight)"
            );
        
        nodeNei.primaryAbscissa() = 
            fvc::interpolate
            (
                node.primaryAbscissa(), 
                nei, 
                "reconstruct(primaryAbscissa)"
            );
    }
}

void Foam::univariateMixing::updateBoundaryQuadrature()
{
    // Recover reference to boundaryField of zero-order moment.
    // All moments will share the same BC types at a given boundary.
    
    volScalarField::GeometricBoundaryField& bf = moments_().boundaryField();
    
    forAll(bf, patchI)
    {
        fvPatchScalarField& m0Patch = bf[patchI];
        
        if (m0Patch.fixesValue())
        {
            forAll(m0Patch, faceI)
            {
                univariateMomentSet momentsToInvert(2.0*nPrimaryNodes_, 0, support_);

                // Copying moments from a face
                momentsToInvert[0] = moments_[0].boundaryField()[patchI][faceI];
                    
                momentsToInvert[1] = moments_[1].boundaryField()[patchI][faceI];
                    
                momentsToInvert[2] = moments_[2].boundaryField()[patchI][faceI];
                    
                momentsToInvert[3] = moments_[3].boundaryField()[patchI][faceI];

                // Inverting moments
                momentsToInvert.invert();
                
                // Recovering primary primaryWeight and primaryAbscissae from moment inverter
                const scalarDiagonalMatrix& pWeights
                (
                    momentsToInvert.weights()
                );
                
                const scalarDiagonalMatrix& pAbscissae
                (
                    momentsToInvert.abscissae()
                );
                
                // Copying quadrature data to boundary face
                for (label pNodeI = 0; pNodeI < nPrimaryNodes_; pNodeI++)
                {
                    basicVolScalarNode& node = nodes_()[pNodeI];
                    
                    node.primaryWeight().boundaryField()[patchI][faceI] =
                        pWeights[pNodeI];
            
                    node.primaryAbscissa().boundaryField()[patchI][faceI] =
                        pAbscissae[pNodeI];
                }
            }
        }
    } 
}

void Foam::univariateMixing::updateQuadrature()
{
    const volScalarField& m0(moments_[0]);
    
    PtrList<basicVolScalarNode>& nodes(nodes_());
    
    forAll(m0, cellI)
    {
        univariateMomentSet momentsToInvert(2.0*nPrimaryNodes_, 0, support_);

        // Copying moment set from a cell to univariateMomentSet
        momentsToInvert[0] = moments_[0][cellI];
        momentsToInvert[1] = moments_[1][cellI];
        momentsToInvert[2] = moments_[2][cellI];
        momentsToInvert[3] = moments_[3][cellI];
        
        // Inverting moments
        momentsToInvert.invert();
        
        // Recovering primary primaryWeight and primaryAbscissae from moment inverter
        const scalarDiagonalMatrix& pWeights
        (
            momentsToInvert.weights()
        );
        
        const scalarDiagonalMatrix& pAbscissae
        (
            momentsToInvert.abscissae()
        );
        
        // Copying to fields
        for (label pNodeI = 0; pNodeI < nPrimaryNodes_; pNodeI++)
        {
            basicVolScalarNode& node(nodes[pNodeI]);

            // Copy primary node
            node.primaryWeight()[cellI] = pWeights[pNodeI];
            node.primaryAbscissa()[cellI] = pAbscissae[pNodeI];
        }
    }

    // Updating boundary conditions
    forAll(nodes, pNodeI)
    {
        basicVolScalarNode& pNode(nodes[pNodeI]);
        pNode.primaryWeight().correctBoundaryConditions();
        pNode.primaryAbscissa().correctBoundaryConditions();
    }
    
    updateBoundaryQuadrature();
    moments_.update();
}

void Foam::univariateMixing::updateAdvection
(
    surfaceScalarField& phiOwn,
    surfaceScalarField& phiNei
)
{
    surfaceScalarField nei
    (
        IOobject
        (
            "nei",
            U_.mesh().time().timeName(),
            U_.mesh()
        ),
        U_.mesh(),
        dimensionedScalar("nei", dimless, -1.0)
    );

    surfaceScalarField own
    (
        IOobject
        (
            "own",
            U_.mesh().time().timeName(),
            U_.mesh()
        ),
        U_.mesh(),
        dimensionedScalar("own", dimless, 1.0)
    );
    
    phiOwn = fvc::interpolate(U_, own, "reconstruct(U)") & U_.mesh().Sf();
    phiNei = fvc::interpolate(U_, nei, "reconstruct(U)") & U_.mesh().Sf();
  
    // Update interpolated nodes
    interpolateNodes();
    
    // Updated reconstructed moments
    momentsNei_.update();
    momentsOwn_.update();
}
 
Foam::tmp<Foam::volScalarField>
Foam::univariateMixing::advectMoment
(
    const basicVolUnivariateMoment& moment,
    const surfaceScalarField& phiOwn,
    const surfaceScalarField& phiNei
)
{
    dimensionedScalar zeroPhi("zero", phiNei.dimensions(), 0.0);
    
    tmp<volScalarField> divMoment
    (
        new volScalarField
        (
            IOobject
            (
                "divMoment",
                U_.mesh().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    
    label order = moment.order();
    
    surfaceScalarField mFlux
    (   
        momentsNei_[order]*min(phiNei, zeroPhi) 
      + momentsOwn_[order]*max(phiOwn, zeroPhi)
    );
 
    fvc::surfaceIntegrate(divMoment.ref(), mFlux);
    divMoment.ref().dimensions().reset(moment.dimensions()/dimTime);
    
    return divMoment;
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::univariateMixing::mixingSource
(
    const basicVolUnivariateMoment& moment,
    const momentFieldSet<basicVolUnivariateMoment,basicVolScalarNode>& moments_
)
{
    return mixingKernel_->Km(moment, moments_);
}

void Foam::univariateMixing::solve()
{
    surfaceScalarField phiOwn("phiOwn", fvc::interpolate(U_) & U_.mesh().Sf());
    surfaceScalarField phiNei("phiNei", fvc::interpolate(U_) & U_.mesh().Sf());
    updateAdvection(phiOwn, phiNei);
    
    // List of moment transport equations
    PtrList<fvScalarMatrix> momentEqns(nMoments_);
    
    forAll(moments_, mI)
    {
        basicVolUnivariateMoment& m = moments_[mI];
        
        momentEqns.set
        (
            mI,
            new fvScalarMatrix
            (
                fvm::ddt(m)
              + advectMoment(m, phiOwn, phiNei)
              - diffusionModel_->momentDiff(m)
              ==
                mixingSource(m,moments_)
            )
        );
    }
    
    forAll (momentEqns, mEqnI)
    {
        momentEqns[mEqnI].relax();
        momentEqns[mEqnI].solve();
    }
    
    updateQuadrature();
    
    if (U_.mesh().time().outputTime() )
    {
        const basicVolScalarNode& node0 = nodes_()[0];
        const basicVolScalarNode& node1 = nodes_()[1];
        node0.primaryAbscissa().write();
        node0.primaryWeight().write();
        node1.primaryAbscissa().write();
        node1.primaryWeight().write();
    }
}
// ************************************************************************* //
