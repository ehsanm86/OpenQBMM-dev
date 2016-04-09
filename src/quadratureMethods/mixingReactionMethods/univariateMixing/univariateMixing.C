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
    nDimensions_(3),                 
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
    reaction_(dict.lookup("reaction")),
    
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
    
    reactionKernel_
    (
        Foam::reactionSubModels::reactionKernel::New
        (
            dict.subDict("reactionKernel")
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
    nodes_ = autoPtr<PtrList<basicVolVectorNode> >
    (
        new PtrList<basicVolVectorNode>
        (
            lookup("nodes"), 
            Foam::basicVolVectorNode::iNew
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
    
    nodesNei_ = autoPtr<PtrList<basicSurfaceVectorNode> >
    (
        new PtrList<basicSurfaceVectorNode>(nPrimaryNodes_)
    );
    
    nodesOwn_ = autoPtr<PtrList<basicSurfaceVectorNode> >
    (
        new PtrList<basicSurfaceVectorNode>(nPrimaryNodes_)
    );
    
    PtrList<basicVolVectorNode>& nodes = nodes_();
    PtrList<basicSurfaceVectorNode>& nodesNei = nodesNei_();
    PtrList<basicSurfaceVectorNode>& nodesOwn = nodesOwn_();

    // Populating interpolated nodes.
    forAll(nodes, pNodeI)
    {
        basicVolVectorNode& node(nodes[pNodeI]);
              
        nodesNei.set
        (
            pNodeI,
            new basicSurfaceVectorNode
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
            new basicSurfaceVectorNode
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
            new Foam::basicSurfaceVectorMoment
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
            new Foam::basicSurfaceVectorMoment
            (
                name_,
                moments_[mI].cmptOrders(),
                nodesOwn_,
                fvc::interpolate(moments_[mI])
            )
        );
    }
    
    updateQuadrature(moments_, nodes_());
    interpolateNodes();
    updateBoundaryQuadrature(moments_, nodes_());
    
    momentsNei_.update();
    momentsOwn_.update();
    reactionKernel_->updateConcentrations(Ca_, Cb_, Cr_, Cs_, moments_);
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
    
    const PtrList<basicVolVectorNode>& nodes = nodes_();
    PtrList<basicSurfaceVectorNode>& nodesNei = nodesNei_();
    PtrList<basicSurfaceVectorNode>& nodesOwn = nodesOwn_();
        
    forAll(nodes, pNodeI)
    {
        const basicVolVectorNode& node(nodes[pNodeI]);
        basicSurfaceVectorNode& nodeOwn(nodesOwn[pNodeI]);
        basicSurfaceVectorNode& nodeNei(nodesNei[pNodeI]);
               
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

void Foam::univariateMixing::updateBoundaryQuadrature
(
    const momentFieldSet
    <
        basicVolVectorMoment,
        basicVolVectorNode
    >& moments,
    PtrList<basicVolVectorNode>& nodes
)
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
                momentsToInvert[0] = 
                    moments(0,0,0).boundaryField()[patchI][faceI];
                    
                momentsToInvert[1] = 
                    moments(1,0,0).boundaryField()[patchI][faceI];
                    
                momentsToInvert[2] = 
                    moments(2,0,0).boundaryField()[patchI][faceI];
                    
                momentsToInvert[3] = 
                    moments(3,0,0).boundaryField()[patchI][faceI];

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
                
                if (momentsToInvert.isDegenerate())
                {
                    basicVolVectorNode& node = nodes[0];
                    
                    node.primaryWeight().boundaryField()[patchI][faceI] =
                            pWeights[0];
                
                    node.primaryAbscissa().boundaryField()[patchI][faceI].component(0) =
                        pAbscissae[0];
                        
                    node.primaryAbscissa().boundaryField()[patchI][faceI].component(1) =
                        scalar(0.0);
                        
                    node.primaryAbscissa().boundaryField()[patchI][faceI].component(2) =
                        scalar(0.0);
                }
                else
                {
                    // Copying quadrature data to boundary face
                    for (label pNodeI = 0; pNodeI < nPrimaryNodes_; pNodeI++)
                    {
                        basicVolVectorNode& node = nodes[pNodeI];
                        
                        node.primaryWeight().boundaryField()[patchI][faceI] =
                            pWeights[pNodeI];
                
                        node.primaryAbscissa().boundaryField()[patchI][faceI].component(0) =
                            pAbscissae[pNodeI];
                    }

                    // Calculate conditional abscissa for reaction
                    if (reaction_)
                    {
                        basicVolVectorNode& node0 = nodes[0];
                        
                        Foam::vector& abscissa0 =
                            node0.primaryAbscissa().boundaryField()[patchI][faceI];
                            
                        basicVolVectorNode& node1 = nodes[1];
                        
                        Foam::vector& abscissa1 =
                            node1.primaryAbscissa().boundaryField()[patchI][faceI];
                        
                        abscissa0.component(1) = (pAbscissae[1]
                          *moments(0,1,0).boundaryField()[patchI][faceI] 
                          - moments(1,1,0).boundaryField()[patchI][faceI])
                          /max((pWeights[0]*(pAbscissae[1] - pAbscissae[0])),1e-12);
                            
                        abscissa1.component(1) = (-pAbscissae[0]
                          *moments(0,1,0).boundaryField()[patchI][faceI] 
                          + moments(1,1,0).boundaryField()[patchI][faceI])
                          /max((pWeights[1]*(pAbscissae[1] - pAbscissae[0])),1e-12);
                            
                        abscissa0.component(2) = (pAbscissae[1]
                          *moments(0,0,1).boundaryField()[patchI][faceI] 
                          - moments(1,0,1).boundaryField()[patchI][faceI])
                          /max((pWeights[0]*(pAbscissae[1] - pAbscissae[0])),1e-12);
                            
                        abscissa1.component(2) = (-pAbscissae[0]
                          *moments(0,0,1).boundaryField()[patchI][faceI] 
                          + moments(1,0,1).boundaryField()[patchI][faceI])
                          /max((pWeights[1]*(pAbscissae[1] - pAbscissae[0])),1e-12);
                    }
                }
            }
        }
    } 
}

void Foam::univariateMixing::updateQuadrature
(
    const momentFieldSet
    <
        basicVolVectorMoment,
        basicVolVectorNode
    >& moments,
    PtrList<basicVolVectorNode>& nodes
)
{
    const volScalarField& m0(moments(0,0,0));
    forAll(m0, cellI)
    {
        univariateMomentSet momentsToInvert(2.0*nPrimaryNodes_, 0, support_);

        // Copying moment set from a cell to univariateMomentSet
        momentsToInvert[0] = moments(0,0,0)[cellI];
        momentsToInvert[1] = moments(1,0,0)[cellI];
        momentsToInvert[2] = moments(2,0,0)[cellI];
        momentsToInvert[3] = moments(3,0,0)[cellI];
        
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
        
//         Info << "weights=" << pWeights << endl;
//         Info << "abscissae=" << pAbscissae << endl;
        
        if (momentsToInvert.isDegenerate())
        {
            basicVolVectorNode& node = nodes[0];
            node.primaryWeight()[cellI] = pWeights[0];
            node.primaryAbscissa()[cellI].component(0) = pAbscissae[0];
            node.primaryAbscissa()[cellI].component(1) = scalar(0);
            node.primaryAbscissa()[cellI].component(2) = scalar(0);
        }
        else
        {
            // Copying to fields
            for (label pNodeI = 0; pNodeI < nPrimaryNodes_; pNodeI++)
            {
                basicVolVectorNode& node = nodes[pNodeI];

                // Copy primary node
                node.primaryWeight()[cellI] = pWeights[pNodeI];
                node.primaryAbscissa()[cellI].component(0) = pAbscissae[pNodeI];
            }
        
            // Calculate conditional abscissa for reaction
            if (reaction_)
            {
                basicVolVectorNode& node0 = nodes[0];
                Foam::vector& abscissa0 = node0.primaryAbscissa()[cellI];
                
                basicVolVectorNode& node1 = nodes[1];
                Foam::vector& abscissa1 = node1.primaryAbscissa()[cellI];
                
                abscissa0.component(1) =
                    (pAbscissae[1]*moments(0,1,0)[cellI]
                  - moments(1,1,0)[cellI])
                  /max((pWeights[0]*(pAbscissae[1] - pAbscissae[0])),1e-12);
                    
                abscissa1.component(1) =
                    (-pAbscissae[0]*moments(0,1,0)[cellI] 
                  + moments(1,1,0)[cellI])
                  /max((pWeights[1]*(pAbscissae[1] - pAbscissae[0])),1e-12);
                    
                abscissa0.component(2) =
                    (pAbscissae[1]*moments(0,0,1)[cellI] 
                  - moments(1,0,1)[cellI])
                  /max((pWeights[0]*(pAbscissae[1] - pAbscissae[0])),1e-12);
                    
                abscissa1.component(2) =
                    (-pAbscissae[0]*moments(0,0,1)[cellI] 
                  + moments(1,0,1)[cellI])
                  /max((pWeights[1]*(pAbscissae[1] - pAbscissae[0])),1e-12);
            }
        }
    }

    // Updating boundary conditions
    forAll(nodes, pNodeI)
    {
        basicVolVectorNode& pNode(nodes[pNodeI]);
        pNode.primaryWeight().correctBoundaryConditions();
        pNode.primaryAbscissa().correctBoundaryConditions();
    }
    
    updateBoundaryQuadrature(moments, nodes);
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
    const basicVolVectorMoment& moment,
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
    
    const labelList& order = moment.cmptOrders();
    
    label a = order[0];
    label b = order[1];
    label c = order[2];
    
    surfaceScalarField mFlux
    (   
        momentsNei_(a,b,c)*min(phiNei, zeroPhi) 
      + momentsOwn_(a,b,c)*max(phiOwn, zeroPhi)
    );
 
    fvc::surfaceIntegrate(divMoment.ref(), mFlux);
    divMoment.ref().dimensions().reset(moment.dimensions()/dimTime);
    
    return divMoment;
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::univariateMixing::mixingSource
(
    const basicVolVectorMoment& moment,
    const momentFieldSet<basicVolVectorMoment,basicVolVectorNode>& moments
)
{
    return mixingKernel_->Km(moment.cmptOrders(), moments);
}

Foam::tmp<Foam::volScalarField>
Foam::univariateMixing::reactionSource
(
    const basicVolVectorMoment& moment,
    const PtrList<basicVolVectorNode>& nodes
)
{
    tmp<volScalarField> rSource
    (
        new volScalarField
        (
            IOobject
            (
                "rSource",
                U_.mesh().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U_.mesh(),
            dimensionedScalar("zero", moment.dimensions()/dimTime, 0.0)
        )
    );

    if (!reaction_)
    {
        rSource.ref().dimensions().reset(moment.dimensions()/dimTime);
        return rSource;
    }
    
    volScalarField& reactionSource = rSource.ref();
    
    const labelList& order = moment.cmptOrders();
    label a = order[0];
    label b = order[1];
    label c = order[2];
    
    // Set pointer lists so quadrature can be run in loop
    const basicVolVectorNode& node0 = nodes[0];
    const basicVolVectorNode& node1 = nodes[1];
    
    PtrList<volScalarField> primaryAbscissa(2);
    primaryAbscissa.set(0, node0.primaryAbscissa().component(0));
    primaryAbscissa.set(1, node1.primaryAbscissa().component(0));
    
    PtrList<volScalarField> Y1(2);
    Y1.set(0, node0.primaryAbscissa().component(1));
    Y1.set(1, node1.primaryAbscissa().component(1));
    
    PtrList<volScalarField> Y2(2);
    Y2.set(0, node0.primaryAbscissa().component(2));
    Y2.set(1, node1.primaryAbscissa().component(2));
    
    PtrList<volScalarField> primaryWeight(2);
    primaryWeight.set(0, node0.primaryWeight());
    primaryWeight.set(1, node1.primaryWeight());
    
    // Compute the reaction kernel
    for (label nodeI = 0; nodeI < nPrimaryNodes_; nodeI++)
    {
        reactionSource == reactionSource
            + primaryWeight[nodeI]*reactionKernel_->R1
            (
                primaryAbscissa[nodeI],
                Y1[nodeI],
                Y2[nodeI]
            )*scalar(b)*pow(primaryAbscissa[nodeI], scalar(a));

        reactionSource == reactionSource
            + primaryWeight[nodeI]*reactionKernel_->R2
            (
                primaryAbscissa[nodeI],
                Y1[nodeI],
                Y2[nodeI]
            )*scalar(c)*pow(primaryAbscissa[nodeI], scalar(a));
    }
    
//     forAll(reactionSource, cellI)
//     {
//         if (Cb_[cellI] < 1.0e-5)
//         {
//             reactionSource[cellI] = scalar(0);
//         }
//         
//         if (Ca_[cellI] < 1.0e-5)
//         {
//             reactionSource[cellI] = scalar(0);
//         }
//         
//         if (Foam::mag(primaryAbscissa[1][cellI]-primaryAbscissa[0][cellI]) < 1.0e-4)
//         {
//             reactionSource[cellI] = scalar(0);
//         }
//     }

//         forAll(reactionSource, cellI)
//     {
//         if (Ca_[cellI] < 0.0)
//         {
//             reactionSource[cellI] = scalar(0);
//             
//             FatalErrorIn
//         (
//             "Foam::mixingReaction\n"
//         )
//             << "momentsOrder=" << moment.cmptOrders()
//             << "moment=" << moment[cellI]
//             << "primaryAbscissae=" << primaryAbscissa[0][cellI]
//             << "primaryAbscissae=" << primaryAbscissa[1][cellI]
//             << "primaryWeight=" << primaryWeight[0][cellI]
//             << "primaryWeight=" << primaryWeight[1][cellI]
//             << "reactionSource=" << reactionSource[cellI]
//             << "cellI=" << cellI
//             << abort(FatalError);
//         }
//     }
    rSource.ref().dimensions().reset(moment.dimensions()/dimTime);
    
    return rSource;
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
        basicVolVectorMoment& m = moments_[mI];
        
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
              + reactionSource(m, nodes_)
            )
        );
    }
    
    forAll (momentEqns, mEqnI)
    {
        momentEqns[mEqnI].relax();
        momentEqns[mEqnI].solve();
    }
    
    updateQuadrature(moments_, nodes_());
    
    reactionKernel_->updateConcentrations(Ca_, Cb_, Cr_, Cs_, moments_);
    
    if (U_.mesh().time().outputTime() )
    {
        const basicVolVectorNode& node0 = nodes_()[0];
        const basicVolVectorNode& node1 = nodes_()[1];
        node0.primaryAbscissa().write();
        node0.primaryWeight().write();
        node1.primaryAbscissa().write();
        node1.primaryWeight().write();
    }
}
// ************************************************************************* //
