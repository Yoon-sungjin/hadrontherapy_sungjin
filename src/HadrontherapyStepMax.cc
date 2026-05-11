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

#include "HadrontherapyStepMax.hh"
#include "HadrontherapyStepMaxMessenger.hh"
#include "HadrontherapyMatrix.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyStepMax::HadrontherapyStepMax(const G4String& processName)
 : G4VDiscreteProcess(processName), 
   MaxChargedStep(DBL_MAX), 
   EnvelopeMaxStep(DBL_MAX), 
   SphereMaxStep(DBL_MAX)
{
  pMess = new HadrontherapyStepMaxMessenger(this);
}
 
/////////////////////////////////////////////////////////////////////////////
HadrontherapyStepMax::~HadrontherapyStepMax() { delete pMess; }

/////////////////////////////////////////////////////////////////////////////
G4bool HadrontherapyStepMax::IsApplicable(const G4ParticleDefinition& particle) 
{ 
  return (particle.GetPDGCharge() != 0.);
}

/////////////////////////////////////////////////////////////////////////////    
void HadrontherapyStepMax::SetMaxStep(G4double step) {MaxChargedStep = step;}
void HadrontherapyStepMax::SetEnvelopeMaxStep(G4double step) {EnvelopeMaxStep = step;}
void HadrontherapyStepMax::SetSphereMaxStep(G4double step) {SphereMaxStep = step;}

/////////////////////////////////////////////////////////////////////////////
G4double HadrontherapyStepMax::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                                  G4double,
                                                  G4ForceCondition* condition )
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  G4double ProposedStep = DBL_MAX;
  
  // 체적 정보가 없으면 패스
  if(aTrack.GetVolume() == 0) return ProposedStep;

  G4String volName = aTrack.GetVolume()->GetName();

  // 1. 현재 모드 확인
  G4String currentMode = "DoseLET"; 
  if (HadrontherapyMatrix::GetInstance()) {
      currentMode = HadrontherapyMatrix::GetInstance()->GetCalculationMode();
  }

  // 2. 모드에 따른 분기 처리
  if (currentMode == "Lineal") 
  {
      // [Lineal 모드] 미시 영역(Envelope, Sphere)은 독자적인 Step 제어
      if (volName == "SensitiveSpherePhys") {
          if (SphereMaxStep < DBL_MAX && SphereMaxStep > 0.) {
              ProposedStep = SphereMaxStep; // 1순위: Sphere 전용 값이 있으면 사용
          } else if (EnvelopeMaxStep < DBL_MAX && EnvelopeMaxStep > 0.) {
              ProposedStep = EnvelopeMaxStep; // 2순위: 값이 없으면 바깥쪽 Envelope 값을 빌려 씀
          }
      }
      else if (volName == "EnvelopeSpherePhys" && EnvelopeMaxStep < DBL_MAX && EnvelopeMaxStep > 0.) {
          ProposedStep = EnvelopeMaxStep;
      }
      // 그 외 거시 영역은 기존 매크로의 MaxChargedStep 적용
      else if ((MaxChargedStep > 0.) && (volName == "DetectorPhys" || volName == "InternalChamber" || volName == "CollimatorHole" || volName == "PhysFourthTQuad" || volName == "PhysThirdTQuad" || volName == "PhysSecondTQuad" || volName == "PhysFirstTQuad" || volName =="physQuadChamber" || volName =="PVirtualMag" || volName =="PhysicCup")) {
          ProposedStep = MaxChargedStep;
      }
  }
  else 
  {
      // [DoseLET 모드] 기존 로직 그대로 유지 (Envelope, Sphere 전용 Step 무시)
      if((MaxChargedStep > 0.) && 
         (volName == "DetectorPhys" || volName == "InternalChamber" || volName == "CollimatorHole" || volName == "PhysFourthTQuad" || volName == "PhysThirdTQuad" || volName == "PhysSecondTQuad" || volName == "PhysFirstTQuad" || volName =="physQuadChamber" || volName =="PVirtualMag" || volName =="PhysicCup"))
      {
          ProposedStep = MaxChargedStep;
      }
  }

  return ProposedStep;
}

/////////////////////////////////////////////////////////////////////////////
G4VParticleChange* HadrontherapyStepMax::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

