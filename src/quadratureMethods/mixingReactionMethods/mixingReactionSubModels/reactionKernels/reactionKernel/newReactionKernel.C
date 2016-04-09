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

#include "reactionKernel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::reactionSubModels::reactionKernel>
Foam::reactionSubModels::reactionKernel::New
(
    const dictionary& dict
)
{
    word reactionKernelType(dict.lookup("reactionKernel"));

    Info<< "Selecting reactionKernel "
        << reactionKernelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(reactionKernelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "reactionKernel::New(const dictionary&) : " << endl
            << "    unknown reactionKernelType type "
            << reactionKernelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid reactionKernelType types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc() << abort(FatalError);
    }

    return autoPtr<reactionKernel>(cstrIter()(dict));
}


// ************************************************************************* //
