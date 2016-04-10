/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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

#include "univariateQuadratureApproximation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateQuadratureApproximation::univariateQuadratureApproximation
(
    const word& name,
    const fvMesh& mesh,
    const word support
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("quadratureProperties", name),
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
    support_(support)
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

    if (nMoments_ != 2*nPrimaryNodes_)
    {
        FatalErrorInFunction
            << "Number of moments from dictionary different from number" << endl
            << "of moments calculated from primary quadrature nodes."
            << abort(FatalError);
    }

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

Foam::univariateQuadratureApproximation::~univariateQuadratureApproximation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::univariateQuadratureApproximation::interpolateNodes()
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
    PtrList<basicSurfaceScalarNode>& nodesNei = nodesNei_();
    PtrList<basicSurfaceScalarNode>& nodesOwn = nodesOwn_();

    forAll(nodes, pNodeI)
    {
        const basicVolScalarNode& node(nodes[pNodeI]);
        basicSurfaceScalarNode& nodeOwn(nodesOwn[pNodeI]);
        basicSurfaceScalarNode& nodeNei(nodesNei[pNodeI]);

        nodeOwn.primaryWeight() =
            fvc::interpolate(node.primaryWeight(), own, "reconstruct(weight)");

        nodeOwn.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                own,
                "reconstruct(abscissa)"
            );

        nodeNei.primaryWeight() =
            fvc::interpolate(node.primaryWeight(), nei, "reconstruct(weight)");

        nodeNei.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                nei,
                "reconstruct(abscissa)"
            );
    }
}

void Foam::univariateQuadratureApproximation::updateBoundaryQuadrature()
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
                univariateMomentSet momentsToInvert(nMoments_, 0, support_);

                // Copying moments from a face
                forAll(momentsToInvert, mI)
                {
                    momentsToInvert[mI]
                        = moments_[mI].boundaryField()[patchI][faceI];
                }

                // Inverting them
                momentsToInvert.invert();

                // Copying quadrature data to boundary face
                for (label pNodeI = 0; pNodeI < nPrimaryNodes_; pNodeI++)
                {
                    basicVolScalarNode& node = nodes_()[pNodeI];

                    node.primaryWeight().boundaryField()[patchI][faceI]
                        = momentsToInvert.weights()[pNodeI];

                    node.primaryAbscissa().boundaryField()[patchI][faceI]
                        = momentsToInvert.abscissae()[pNodeI];
                }

            }
        }
    }
}

void Foam::univariateQuadratureApproximation::updateQuadrature()
{
    const volScalarField& m0(moments_[0]);

    PtrList<basicVolScalarNode>& nodes(nodes_());

    forAll(m0, cellI)
    {
        univariateMomentSet momentsToInvert(nMoments_, 0.0, support_);

        // Copying moment set from a cell to univariateMomentSet
        forAll(momentsToInvert, mI)
        {
            momentsToInvert[mI] = moments_[mI][cellI];
        }

        // Inverting moments and updating secondary quadrature
        momentsToInvert.invert();

        // Recovering primary weights and abscissae from moment inverter
        const scalarDiagonalMatrix& pWeights(momentsToInvert.weights());

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
    updateMoments();
}


void Foam::univariateQuadratureApproximation::updateMoments()
{
    moments_.update();
}


// ************************************************************************* //
