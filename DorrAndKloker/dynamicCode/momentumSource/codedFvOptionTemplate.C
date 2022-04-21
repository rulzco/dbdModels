/*---------------------------------------------------------------------------*  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "codedFvOptionTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"

//{{{ begin codeInclude
#line 29 "/home/rzco/OpenFOAM/rzco-8/run/kloker/kloker/system/fvOptions/momentumSource"

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = d74c4e5e5c1f24eb9d4eb164ed020ba975892247
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void momentumSource_d74c4e5e5c1f24eb9d4eb164ed020ba975892247(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//makeRemovablePatchTypeField
//(
//    fvPatch,
//    momentumSourceFvOptionvectorSource
//);
defineTypeNameAndDebug(momentumSourceFvOptionvectorSource, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    momentumSourceFvOptionvectorSource,
    dictionary
);


const char* const momentumSourceFvOptionvectorSource::SHA1sum =
    "d74c4e5e5c1f24eb9d4eb164ed020ba975892247";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

momentumSourceFvOptionvectorSource::
momentumSourceFvOptionvectorSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh)
{
    if (false)
    {
        Info<<"construct momentumSource sha1: d74c4e5e5c1f24eb9d4eb164ed020ba975892247"
            " from components\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

momentumSourceFvOptionvectorSource::
~momentumSourceFvOptionvectorSource()
{
    if (false)
    {
        Info<<"destroy momentumSource sha1: d74c4e5e5c1f24eb9d4eb164ed020ba975892247\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void momentumSourceFvOptionvectorSource::correct
(
    GeometricField<vector, fvPatchField, volMesh>& fld
)
{
    if (false)
    {
        Info<<"momentumSourceFvOptionvectorSource::correct()\n";
    }

//{{{ begin code
    #line 34 "/home/rzco/OpenFOAM/rzco-8/run/kloker/kloker/system/fvOptions/momentumSource"
//Pout<< "**codeCorrect**" << endl;
//}}} end code
}


void momentumSourceFvOptionvectorSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<<"momentumSourceFvOptionvectorSource::addSup()\n";
    }

//{{{ begin code
    #line 39 "/home/rzco/OpenFOAM/rzco-8/run/kloker/kloker/system/fvOptions/momentumSource"
//const Time& time = mesh().time();
        const scalarField& V = mesh_.V();
        const vectorField& cells = mesh_.C();
        vectorField& uSource = eqn.source();
        //const scalarField& x = mesh_.C().component(0);
        //const scalarField& y = mesh_.C().component(1);
        
		//x_pa = 0.0
		float a0 = 55;
		float a1 = 8;
		float a2 = 10; 
		float b0 = 34;
		float b1 = 2.7;
		float b2 = 0.7;
		float cx = 80;
		
		const scalar xmin = 0;
		const scalar xmax = 0.1;
		const scalar ymax = 0.1;
		float rho0 = 1.225;
		float u0 = 34.922;
		float L = 0.04;
		float magF = pow(u0,2)*rho0/L;
        
        //cellSet selectedCells (mesh_, cellSetName_);
        //labelList cells = selectedCells.toc ();
                    
        DimensionedField <scalar, volMesh> bForce
        (
            IOobject
            (
                "bForce",
                mesh_.time().timeName(),
                mesh_,
                IOobject :: NO_READ,
                IOobject :: AUTO_WRITE
            ),
            mesh_,
            dimensioned <scalar>
            (
                 "0",
                dimAcceleration,
                pTraits <scalar> :: zero
            )
        );
        
        
        if(mesh().time().value()>0)
        {
            forAll(cells, i)
            {
                const scalar x = cells[i][0];
                const scalar y = cells[i][1];
                if( x > xmin && x < xmax && y < ymax)
                {

					bForce[i] = cx *(a0*a1*x + pow(a0,2)*a2*pow(x,2)) * exp(-a0*x)*
					(b1*y + b2*pow(y,2)) * exp(-b0*pow(y,0.4))*magF;

                    uSource[i][0] -= bForce[i] *V[i];
                    //uSource[i][1] -= Ey * ec * rhoC * omega * deltaT * V[i];
                }
                else{bForce[i] = 0;}
            }
            
        }
                
        
		if(mesh_.time().writeTime())
        {
			bForce.write();
			Info << "magF = " << magF << endl;
		}
//}}} end code
}


void momentumSourceFvOptionvectorSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<<"momentumSourceFvOptionvectorSource::addSup()\n";
    }

//{{{ begin code
    #line 39 "/home/rzco/OpenFOAM/rzco-8/run/kloker/kloker/system/fvOptions/momentumSource"
//const Time& time = mesh().time();
        const scalarField& V = mesh_.V();
        const vectorField& cells = mesh_.C();
        vectorField& uSource = eqn.source();
        //const scalarField& x = mesh_.C().component(0);
        //const scalarField& y = mesh_.C().component(1);
        
		//x_pa = 0.0
		float a0 = 55;
		float a1 = 8;
		float a2 = 10; 
		float b0 = 34;
		float b1 = 2.7;
		float b2 = 0.7;
		float cx = 80;
		
		const scalar xmin = 0;
		const scalar xmax = 0.1;
		const scalar ymax = 0.1;
		float rho0 = 1.225;
		float u0 = 34.922;
		float L = 0.04;
		float magF = pow(u0,2)*rho0/L;
        
        //cellSet selectedCells (mesh_, cellSetName_);
        //labelList cells = selectedCells.toc ();
                    
        DimensionedField <scalar, volMesh> bForce
        (
            IOobject
            (
                "bForce",
                mesh_.time().timeName(),
                mesh_,
                IOobject :: NO_READ,
                IOobject :: AUTO_WRITE
            ),
            mesh_,
            dimensioned <scalar>
            (
                 "0",
                dimAcceleration,
                pTraits <scalar> :: zero
            )
        );
        
        
        if(mesh().time().value()>0)
        {
            forAll(cells, i)
            {
                const scalar x = cells[i][0];
                const scalar y = cells[i][1];
                if( x > xmin && x < xmax && y < ymax)
                {

					bForce[i] = cx *(a0*a1*x + pow(a0,2)*a2*pow(x,2)) * exp(-a0*x)*
					(b1*y + b2*pow(y,2)) * exp(-b0*pow(y,0.4))*magF;

                    uSource[i][0] -= bForce[i] *V[i];
                    //uSource[i][1] -= Ey * ec * rhoC * omega * deltaT * V[i];
                }
                else{bForce[i] = 0;}
            }
            
        }
                
        
		if(mesh_.time().writeTime())
        {
			bForce.write();
			Info << "magF = " << magF << endl;
		}
//}}} end code
}


void momentumSourceFvOptionvectorSource::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<<"momentumSourceFvOptionvectorSource::constrain()\n";
    }

//{{{ begin code
    #line 118 "/home/rzco/OpenFOAM/rzco-8/run/kloker/kloker/system/fvOptions/momentumSource"
//Pout<< "**codeSetValue**" << endl;
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv
// ************************************************************************* //

