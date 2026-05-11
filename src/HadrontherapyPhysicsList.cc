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
//
//    ******      SUGGESTED PHYSICS FOR ACCURATE SIMULATIONS    *********
//    ******            IN MEDICAL PHYSICS APPLICATIONS         *********
//
// 'HADRONTHERAPY_1' and 'HADRONTHERAPY_2' are both suggested;
// It can be activated inside any macro file using the command:
// /Physics/addPhysics HADRONTHERAPY_1 (HADRONTHERAPY_2)

#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyPhysicsListMessenger.hh"
#include "HadrontherapyStepMax.hh"
#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessManager.hh"
#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4AutoDelete.hh"
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"

// for penelope model
#include "G4EmConfigurator.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4PenelopeBremsstrahlungModel.hh"
#include "G4PenelopeAnnihilationModel.hh"
#include "G4PenelopePhotoElectricModel.hh"
#include "G4PenelopeComptonModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4PenelopeRayleighModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "HadrontherapyMatrix.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyPhysicsList::HadrontherapyPhysicsList() : G4VModularPhysicsList()
{
    G4LossTableManager::Instance();
    defaultCutValue = 1.*mm;
    cutForGamma     = defaultCutValue;
    cutForElectron  = defaultCutValue;
    cutForPositron  = defaultCutValue;
    
    pMessenger = new HadrontherapyPhysicsListMessenger(this);
    SetVerboseLevel(1);
    decay_List = new G4DecayPhysics();
    // Elecromagnetic physics
    //
    emPhysicsList = new G4EmStandardPhysics_option4();
    emName = "standard_opt4";
    
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyPhysicsList::~HadrontherapyPhysicsList()
{
    delete pMessenger;
    delete emPhysicsList;
    delete decay_List;
    //delete radioactiveDecay_List;
    hadronPhys.clear();
    for(size_t i=0; i<hadronPhys.size(); i++)
    {
        delete hadronPhys[i];
    }
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::ConstructParticle()
{
    decay_List -> ConstructParticle();
    
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::ConstructProcess()
{
    // Transportation
    //
    AddTransportation();
    
    decay_List -> ConstructProcess();
    emPhysicsList -> ConstructProcess();
    
    
    //em_config.AddModels();
    
    // Hadronic physics
    //
    // Penelope Hybrid mode
    if (emName == "EmPenelope") {
        G4EmConfigurator* emConfig = G4LossTableManager::Instance()->EmConfigurator();
        // 1. 엄마 구역 (DetectorLog)에 적용
        G4String regionName1 = "DetectorLog"; 
        G4Region* region1 = G4RegionStore::GetInstance()->GetRegion(regionName1, false);
        
        // 🌟 [핵심 추가] 자식 구역 (EnvelopeLog)에도 똑같이 적용 🌟
        G4String regionName2 = "EnvelopeLog"; 
        G4Region* region2 = G4RegionStore::GetInstance()->GetRegion(regionName2, false);

        if (region1 || region2) {
            G4double highLimit = 1.0 * GeV; 

            // --- Electron Models ---
            // Ionisation: Penelope model + Universal Fluctuation
            G4PenelopeIonisationModel* eIoni = new G4PenelopeIonisationModel();
            eIoni->SetHighEnergyLimit(highLimit);
            if(region1) emConfig->SetExtraEmModel("e-", "eIoni", eIoni, regionName1, 0.0, highLimit, new G4UniversalFluctuation());
            if(region2) emConfig->SetExtraEmModel("e-", "eIoni", eIoni, regionName2, 0.0, highLimit, new G4UniversalFluctuation());

            // Bremsstrahlung
            G4PenelopeBremsstrahlungModel* eBrem = new G4PenelopeBremsstrahlungModel();
            eBrem->SetHighEnergyLimit(highLimit);
            if(region1) emConfig->SetExtraEmModel("e-", "eBrem", eBrem, regionName1, 0.0, highLimit);
            if(region2) emConfig->SetExtraEmModel("e-", "eBrem", eBrem, regionName2, 0.0, highLimit);
            
            // --- Positron Models ---
            // Ionisation
            G4PenelopeIonisationModel* pIoni = new G4PenelopeIonisationModel();
            pIoni->SetHighEnergyLimit(highLimit);
            if(region1) emConfig->SetExtraEmModel("e+", "eIoni", pIoni, regionName1, 0.0, highLimit, new G4UniversalFluctuation());
            if(region2) emConfig->SetExtraEmModel("e+", "eIoni", pIoni, regionName2, 0.0, highLimit, new G4UniversalFluctuation());
            // Bremsstrahlung
            G4PenelopeBremsstrahlungModel* pBrem = new G4PenelopeBremsstrahlungModel();
            pBrem->SetHighEnergyLimit(highLimit);
            if(region1) emConfig->SetExtraEmModel("e+", "eBrem", pBrem, regionName1, 0.0, highLimit);
            if(region2) emConfig->SetExtraEmModel("e+", "eBrem", pBrem, regionName2, 0.0, highLimit);
            // Annihilation
            G4PenelopeAnnihilationModel* pAnn = new G4PenelopeAnnihilationModel();
            pAnn->SetHighEnergyLimit(highLimit);
            if(region1) emConfig->SetExtraEmModel("e+", "annihil", pAnn, regionName1, 0.0, highLimit);
            if(region2) emConfig->SetExtraEmModel("e+", "annihil", pAnn, regionName2, 0.0, highLimit);
            // --- Gamma Models ---
            // PhotoElectric
            G4PenelopePhotoElectricModel* gPhot = new G4PenelopePhotoElectricModel();
            gPhot->SetHighEnergyLimit(highLimit);
            if(region1) emConfig->SetExtraEmModel("gamma", "phot", gPhot, regionName1, 0.0, highLimit);
            if(region2) emConfig->SetExtraEmModel("gamma", "phot", gPhot, regionName2, 0.0, highLimit);
            // Compton
            G4PenelopeComptonModel* gComp = new G4PenelopeComptonModel();
            gComp->SetHighEnergyLimit(highLimit);
            if(region1) emConfig->SetExtraEmModel("gamma", "compt", gComp, regionName1, 0.0, highLimit);
            if(region2) emConfig->SetExtraEmModel("gamma", "compt", gComp, regionName2, 0.0, highLimit);
            // Gamma Conversion
            G4PenelopeGammaConversionModel* gConv = new G4PenelopeGammaConversionModel();
            gConv->SetHighEnergyLimit(highLimit);
            if(region1) emConfig->SetExtraEmModel("gamma", "conv", gConv, regionName1, 0.0, highLimit);
            if(region2) emConfig->SetExtraEmModel("gamma", "conv", gConv, regionName2, 0.0, highLimit);
            // Rayleigh Scattering
            G4PenelopeRayleighModel* gRayl = new G4PenelopeRayleighModel();
            gRayl->SetHighEnergyLimit(highLimit);
            if(region1) emConfig->SetExtraEmModel("gamma", "Rayl", gRayl, regionName1, 0.0, highLimit);
            if(region2) emConfig->SetExtraEmModel("gamma", "Rayl", gRayl, regionName2, 0.0, highLimit);
            G4cout << ">>> [HadrontherapyPhysicsList] PENELOPE EM Physics Activated for Regions: DetectorLog and EnvelopeLog <<<" << G4endl;
            
        }
        else {
            G4cout << ">>> [HadrontherapyPhysicsList] WARNING: Regions NOT FOUND. Using Standard Opt4. <<<" << G4endl;
        }
    }


    for(size_t i=0; i < hadronPhys.size(); i++)
    {
        hadronPhys[i] -> ConstructProcess();
    }
    
    // step limitation (as a full process)
    //
    AddStepMax();
    
    //Parallel world sensitivity
    //
    G4ParallelWorldPhysics* pWorld = new G4ParallelWorldPhysics("DetectorROGeometry");
    pWorld->ConstructProcess();
    
    return;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::AddPhysicsList(const G4String& name)
{
    if (verboseLevel>1) {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
    }
    if (name == emName) return;
    
    ///////////////////////////////////
    //   ELECTROMAGNETIC MODELS
    ///////////////////////////////////
    if (name == "standard_opt4") {
        emName = name;
        delete emPhysicsList;
        hadronPhys.clear();
        emPhysicsList = new G4EmStandardPhysics_option4();
        G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
        G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option4" << G4endl;
        
        ////////////////////////////////////////
        //   ELECTROMAGNETIC + HADRONIC MODELS
        ////////////////////////////////////////
        
    }  else if (name == "HADRONTHERAPY_1") {
        
        AddPhysicsList("standard_opt4");
        hadronPhys.push_back( new G4DecayPhysics());
        hadronPhys.push_back( new G4RadioactiveDecayPhysics());
        hadronPhys.push_back( new G4IonBinaryCascadePhysics());
        hadronPhys.push_back( new G4EmExtraPhysics());
        hadronPhys.push_back( new G4HadronElasticPhysicsHP());
        hadronPhys.push_back( new G4StoppingPhysics());
        hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC_HP());
        hadronPhys.push_back( new G4NeutronTrackingCut());
        
        G4cout << "HADRONTHERAPY_1 PHYSICS LIST has been activated" << G4endl;
    }
    
    else if (name == "HADRONTHERAPY_2") {
        // HP models are switched off
        AddPhysicsList("standard_opt4");
        hadronPhys.push_back( new G4DecayPhysics());
        hadronPhys.push_back( new G4RadioactiveDecayPhysics());
        hadronPhys.push_back( new G4IonBinaryCascadePhysics());
        hadronPhys.push_back( new G4EmExtraPhysics());
        hadronPhys.push_back( new G4HadronElasticPhysics());
        hadronPhys.push_back( new G4StoppingPhysics());
        hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC());
        hadronPhys.push_back( new G4NeutronTrackingCut());
        
        G4cout << "HADRONTHERAPY_2 PHYSICS LIST has been activated" << G4endl;    }

    else if (name == "EmPenelope") {
        emName = name;
        delete emPhysicsList;
        emPhysicsList = new G4EmStandardPhysics_option4();
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();

        G4cout << ">>> Requesting PENELOPE EM Physics (Hybrid Mode) <<<" << G4endl;
        G4cout << ">>> Base Physics: G4EmStandardPhysics_option4 / Region Overlay: DetectorLog <<<" << G4endl;
    }

    else {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
        << " is not defined"
        << G4endl;
    }
    
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyPhysicsList::AddStepMax()
{
    // Step limitation seen as a process
    // This process must exist in all threads.
    //
    HadrontherapyStepMax* stepMaxProcess  = new HadrontherapyStepMax();
    
    auto particleIterator = GetParticleIterator();
    particleIterator->reset();
    while ((*particleIterator)()){
        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        
        if (stepMaxProcess->IsApplicable(*particle) && pmanager)
        {
            pmanager ->AddDiscreteProcess(stepMaxProcess);
        }
    }
}

void HadrontherapyPhysicsList::SetCuts()
{
    // 1. 기본 Cut 설정 (매크로의 /run/setCut 명령어로 입력된 값이 defaultCutValue로 들어옴)
    if (verboseLevel > 0) {
        G4cout << "HadrontherapyPhysicsList::SetCuts: default cut length : " 
               << G4BestUnit(defaultCutValue,"Length") << G4endl;
    }

    SetCutValue(cutForGamma, "gamma");
    SetCutValue(cutForElectron, "e-");
    SetCutValue(cutForPositron, "e+");

    // // 2. 모드에 따른 분기 처리 (Lineal 모드일 때만)
    // G4String currentMode = "DoseLET"; // 기본값
    // if (HadrontherapyMatrix::GetInstance()) {
    //     currentMode = HadrontherapyMatrix::GetInstance()->GetCalculationMode();
    // }

    // if (currentMode == "Lineal") 
    // {
    //     // Lineal 모드: 매크로의 /run/setCutForRegion EnvelopeLog 0.6 um 명령어가 
    //     // Geant4 코어에 의해 이 시점 이후에 적용됩니다. 
    //     G4cout << "===> [Lineal Mode] Physics Culling Strategy Active." << G4endl;
    //     G4cout << "===> Awaiting macro /run/setCutForRegion EnvelopeLog commands..." << G4endl;
    // }
    // else 
    // {
    //     // DoseLET 모드: Envelope 영역에 미세 Cut이 적용되면 안 되므로 방어 코드를 작성합니다.
    //     G4Region* envRegion = G4RegionStore::GetInstance()->GetRegion("EnvelopeLog", false);
    //     if (envRegion) 
    //     {
    //         // DoseLET 모드에서는 Envelope 영역의 Cut을 
    //         // 물 팬텀 전체의 기본 Cut(defaultCutValue)과 강제로 동일하게 맞춥니다.
    //         G4ProductionCuts* stdCuts = new G4ProductionCuts();
    //         stdCuts->SetProductionCut(defaultCutValue, "gamma");
    //         stdCuts->SetProductionCut(defaultCutValue, "e-");
    //         stdCuts->SetProductionCut(defaultCutValue, "e+");
            
    //         envRegion->SetProductionCuts(stdCuts);
            
    //         G4cout << "===> [DoseLET Mode] Envelope region detected but fine cuts are disabled." << G4endl;
    //         G4cout << "===> Envelope cut forced to default: " << G4BestUnit(defaultCutValue, "Length") << G4endl;
    //     }
    // }

    if (verboseLevel > 0) DumpCutValuesTable();
}
