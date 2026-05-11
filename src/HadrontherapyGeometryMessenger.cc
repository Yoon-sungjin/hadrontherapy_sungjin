//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

//
#include "HadrontherapyGeometryMessenger.hh"
#include "HadrontherapyGeometryController.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "HadrontherapyMatrix.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyGeometryMessenger::HadrontherapyGeometryMessenger(HadrontherapyGeometryController* controller)
:hadrontherapyGeometryController(controller)

{
    changeTheGeometryDir = new G4UIdirectory("/geometrySetup/");
    changeTheGeometryDir -> SetGuidance("Geometry setup");
    
    changeTheGeometryCmd = new G4UIcmdWithAString("/geometrySetup/selectGeometry",this);
    changeTheGeometryCmd -> SetGuidance("Select the geometry you wish to use");
    changeTheGeometryCmd -> SetParameterName("Geometry",false);
    changeTheGeometryCmd -> AvailableForStates(G4State_PreInit);
    
    /*  changeTheDetectorCmd = new G4UIcmdWithAString("/geometrySetup/selectDetector",this);
    changeTheDetectorCmd -> SetGuidance("Select the detector you wish to use");
    changeTheDetectorCmd -> SetParameterName("Detector",false);
    changeTheDetectorCmd -> AvailableForStates(G4State_PreInit);*/

    // 생성자(Constructor) 내부에 커맨드 정의 추가
    sphereRadiusCmd = new G4UIcmdWithADoubleAndUnit("/changeDetector/sphereRadius", this);
    sphereRadiusCmd->SetGuidance("Set the radius of the sensitive sphere for Lineal Energy calculation");
    sphereRadiusCmd->SetParameterName("radius", false);
    sphereRadiusCmd->SetDefaultUnit("um");
    sphereRadiusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    envelopeRadiusCmd = new G4UIcmdWithADoubleAndUnit("/changeDetector/envelopeRadius", this);
    envelopeRadiusCmd->SetGuidance("Set the radius of the envelope for variance reduction");
    envelopeRadiusCmd->SetParameterName("envRadius", false);
    envelopeRadiusCmd->SetDefaultUnit("mm");
    envelopeRadiusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyGeometryMessenger::~HadrontherapyGeometryMessenger()
{
    delete changeTheGeometryDir;
    delete changeTheGeometryCmd;
    delete sphereRadiusCmd;
    delete envelopeRadiusCmd;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyGeometryMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

    if( command == changeTheGeometryCmd )
    {
        hadrontherapyGeometryController -> SetGeometry (newValue);
    }

    if (command == sphereRadiusCmd)
    {
        HadrontherapyMatrix::sphereRadius = sphereRadiusCmd->GetNewDoubleValue(newValue);
    }
    if (command == envelopeRadiusCmd)
    {
        HadrontherapyMatrix::envelopeRadius = envelopeRadiusCmd->GetNewDoubleValue(newValue);
    }
}
