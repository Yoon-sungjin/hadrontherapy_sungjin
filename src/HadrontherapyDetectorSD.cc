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

#include "HadrontherapyDetectorSD.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4TrackVector.hh"
#include "G4TransportationManager.hh"
#include "G4UserEventAction.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VTouchable.hh"
#include "G4ios.hh"

#include "HadrontherapyDetectorHit.hh"
#include "HadrontherapyLet.hh"
#include "HadrontherapyMatrix.hh"
#include "HadrontherapyRBE.hh"
#include "HadrontherapyRunAction.hh"
#include "HadrontherapySteppingAction.hh"
#include <G4AccumulableManager.hh>

#include "HadrontherapyDetectorConstruction.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorSD::HadrontherapyDetectorSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname = "HadrontherapyDetectorHitsCollection");
  HitsCollection = NULL;
  sensitiveDetectorName = name;
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorSD::~HadrontherapyDetectorSD() {}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorSD::Initialize(G4HCofThisEvent*)
{
  HitsCollection =
    new HadrontherapyDetectorHitsCollection(sensitiveDetectorName, collectionName[0]);
}

G4bool HadrontherapyDetectorSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4String volName = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  // Matrix에서 현재 매크로로 설정된 모드를 실시간으로 가져옵니다.
  // G4String currentMode = HadrontherapyMatrix::GetInstance()->GetCalculationMode();

  G4String currentMode = HadrontherapyMatrix::GetInstance()->GetCalculationMode();

  G4double edep = aStep->GetTotalEnergyDeposit();
  // if (edep > 0.) {
  //     G4cout << "[DEBUG-HIT] Vol: " << volName 
  //            << " | Mode: " << currentMode 
  //            << " | Edep: " << edep/keV << " keV" << G4endl;
  // }

  // 1. Lineal 모드 켰을 때
  if (currentMode == "Lineal" && volName == "SensitiveSpherePhys")
  {
    // G4cout << "  ---> [SUCCESS] ProcessHitsLineal success!" << G4endl;
    // 구 안에 있을 때만 에너지 기록! (스텝 2만 통과)
    return ProcessHitsLineal(aStep);
  }
  // 2. DoseLET 모드 켰을 때
  else if (currentMode == "DoseLET" && volName == "RODetectorZDivisionPhys")
  {
    // 구 밖에 있고 복셀 안에 있을 때만 에너지 기록! (스텝 1, 3만 통과)
    return ProcessHitsDoseLET(aStep);
  }
  // 그 외의 경우(모드가 안 맞거나 다른 물체인 경우)는 기록하지 않고 무시!
  return false;
}

/////////////////////////////////////////////////////////////////////////////
G4bool HadrontherapyDetectorSD::ProcessHitsLineal(G4Step* aStep)
{
  // 민감 체적(Sensitive Volume) 확인
  // if (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "SensitiveSpherePhys")
  //   return false;

  // 기본 물리량 추출
  G4Track* theTrack = aStep->GetTrack();
  G4ParticleDefinition* particleDef = theTrack->GetDefinition();
  G4int pdg = particleDef->GetPDGEncoding();
  G4int trackID = theTrack->GetTrackID();
  G4int Z = particleDef->GetAtomicNumber();
  G4int A = particleDef->GetAtomicMass();
  G4double energyDeposit = aStep->GetTotalEnergyDeposit();
  G4double DX = aStep->GetStepLength();
  G4double kineticEnergy = theTrack->GetKineticEnergy();

  HadrontherapyMatrix* matrix = HadrontherapyMatrix::GetInstance();
  if (!matrix) return false;

  // ************************************************************
  // Replica 계층 구조를 타고 올라가서 X, Y, Z 인덱스 추출 
  // ************************************************************
  const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();

  // 현재 입자가 SensitiveSphere 안에 있을 때의 볼륨 깊이(Depth)
  // Depth 0: SensitiveSpherePhys (현재 위치)
  // Depth 1: EnvelopeSpherePhys (엄마)
  // Depth 2: VoxelZPhys (할머니: Z축 복셀)
  // Depth 3: RepYPhys (증조할머니: Y축 복셀)
  // Depth 4: RepXPhys (고조할머니: X축 복셀)
  
  G4int k = touchable->GetReplicaNumber(2); // Z 축 인덱스
  G4int j = touchable->GetReplicaNumber(3); // Y 축 인덱스
  G4int i = touchable->GetReplicaNumber(4); // X 축 인덱스

  // 4. Matrix (Dose & Fluence) 데이터 처리
  if (matrix)
  {
    G4int* hitTrack = matrix->GetHitTrack(i, j, k);

    // Fluence 업데이트 (이벤트당 트랙 ID가 처음 등장할 때만)
    if (*hitTrack != trackID)
    {
      *hitTrack = trackID;
      if (Z >= 1) matrix->Fill(trackID, particleDef, i, j, k, 0, true);
    }

    // Dose 업데이트 (직접 에너지 침전 + 생성된 전자 에너지 합산)
    if (energyDeposit > 0.)
    {
      // [트랙 1] 이온별(Ion-specific) 선량 데이터 (전자, 감마 등은 제외)
      if (Z >= 1)
      {
        matrix->Fill(trackID, particleDef, i, j, k, energyDeposit);
      }

      // [트랙 2] Total Dose 및 시각화용 Hit 생성
      // Z >= 1 조건이 없으므로, 2차 전자가 남긴 energyDeposit도 무사히 Hit으로 생성되어 저장됨!
      HadrontherapyDetectorHit* detectorHit = new HadrontherapyDetectorHit();
      detectorHit->SetEdepAndPosition(i, j, k, energyDeposit);
      HitsCollection->insert(detectorHit);
    }
  }

  // 5. RBE (Relative Biological Effectiveness) 처리
  auto rbe = HadrontherapyRBE::GetInstance();
  if (rbe->IsCalculationEnabled() && A > 0)
  {
    if (!fRBEAccumulable)
    {
      fRBEAccumulable = dynamic_cast<HadrontherapyRBEAccumulable*>(
        G4AccumulableManager::Instance()->GetAccumulable("RBE"));
    }

    if (fRBEAccumulable)
    {
      fRBEAccumulable->Accumulate(kineticEnergy / A, energyDeposit, DX, Z, i, j, k);
    }
  }

  return true;
}

G4bool HadrontherapyDetectorSD::ProcessHitsDoseLET(G4Step* aStep)
{
  // 기본 물리량 추출
  G4Track* theTrack = aStep->GetTrack();
  G4ParticleDefinition* particleDef = theTrack->GetDefinition();
  G4int pdg = particleDef->GetPDGEncoding();
  G4int trackID = theTrack->GetTrackID();
  G4int Z = particleDef->GetAtomicNumber();
  G4int A = particleDef->GetAtomicMass();
  G4double energyDeposit = aStep->GetTotalEnergyDeposit();
  G4double DX = aStep->GetStepLength();
  G4double kineticEnergy = theTrack->GetKineticEnergy();

  // 복셀 인덱스 추출 (Dose/LET Geometry 기준)
  // [중요] Lineal과 달리 복셀 순서가 다릅니다. (X=2, Y=1, Z=0)
  const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int k = touchable->GetReplicaNumber(0);  // Z
  G4int j = touchable->GetReplicaNumber(1);  // Y
  G4int i = touchable->GetReplicaNumber(2);  // X

  // 싱글톤 인스턴스 가져오기 (한 번만 선언!)
  HadrontherapyMatrix* matrix = HadrontherapyMatrix::GetInstance();
  HadrontherapyLet* let = HadrontherapyLet::GetInstance();

  // 3. LET (Linear Energy Transfer) 데이터 처리
  if (let && Z >= 1 && energyDeposit > 0. && DX > 0.)
  {
    // 중성자, 감마, 전자를 제외한 이온만 처리
    if (pdg != 22 && pdg != 11 && !(Z == 0 && A == 1))
    {
      G4double eKinPre = aStep->GetPreStepPoint()->GetKineticEnergy();
      G4double eKinPost = aStep->GetPostStepPoint()->GetKineticEnergy();
      G4double eKinMean = (eKinPre + eKinPost) * 0.5;
      const G4Material* materialStep = aStep->GetPreStepPoint()->GetMaterial();

      G4double energySecondaryElectrons = 0.;
      const std::vector<const G4Track*>* secondary = aStep->GetSecondaryInCurrentStep();

      if (secondary && !secondary->empty())
      {
        for (size_t numsec = 0; numsec < secondary->size(); numsec++)
        {
          if ((*secondary)[numsec]->GetDefinition()->GetPDGEncoding() == 11)  // Electron (e-)
          {
            energySecondaryElectrons += (*secondary)[numsec]->GetKineticEnergy();
          }
        }
      }

      let->FillEnergySpectrum(trackID, particleDef, eKinMean, materialStep, energyDeposit,
                              energySecondaryElectrons, DX, i, j, k);
    }
  }

  // 4. Matrix (Dose & Fluence) 데이터 처리
  if (matrix)
  {
    G4int* hitTrack = matrix->GetHitTrack(i, j, k);

    // Fluence 업데이트
    if (*hitTrack != trackID)
    {
      *hitTrack = trackID;
      if (Z >= 1) matrix->Fill(trackID, particleDef, i, j, k, 0, true);
    }

    // Dose 업데이트
    if (energyDeposit > 0.)
    {
      // [트랙 1] 이온별(Ion-specific) 선량 데이터 (전자, 감마 등은 제외)
      if (Z >= 1)
      {
        matrix->Fill(trackID, particleDef, i, j, k, energyDeposit);
      }

      // [트랙 2] Total Dose 및 시각화용 Hit 생성
      // Z >= 1 조건이 없으므로, 2차 전자가 남긴 energyDeposit도 무사히 Hit으로 생성되어 저장됨!
      HadrontherapyDetectorHit* detectorHit = new HadrontherapyDetectorHit();
      detectorHit->SetEdepAndPosition(i, j, k, energyDeposit);
      HitsCollection->insert(detectorHit);
    }
  }

  // 5. RBE (Relative Biological Effectiveness) 처리
  auto rbe = HadrontherapyRBE::GetInstance();
  if (rbe->IsCalculationEnabled() && A > 0)
  {
    if (!fRBEAccumulable)
    {
      fRBEAccumulable = dynamic_cast<HadrontherapyRBEAccumulable*>(
        G4AccumulableManager::Instance()->GetAccumulable("RBE"));
    }

    if (fRBEAccumulable)
    {
      fRBEAccumulable->Accumulate(kineticEnergy / A, energyDeposit, DX, Z, i, j, k);
    }
  }

  return true;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if (HCID < 0)
  {
    HCID = GetCollectionID(0);
  }

  HCE->AddHitsCollection(HCID, HitsCollection);
}
