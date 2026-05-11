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

#include "HadrontherapyDetectorROGeometry.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RotationMatrix.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyDummySD.hh"
#include "HadrontherapyMatrix.hh"
#include "G4Orb.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "HadrontherapyDetectorConstruction.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorROGeometry::HadrontherapyDetectorROGeometry(G4String aString)
  : G4VUserParallelWorld(aString),
    RODetector(0),
    RODetectorXDivision(0),
    RODetectorYDivision(0),
    RODetectorZDivision(0),
    worldLogical(0),
    RODetectorLV(0),
    RODetectorXDivisionLV(0),
    RODetectorYDivisionLV(0),
    RODetectorZDivisionLV(0),
    sensitiveLogicalVolume(0)
{
  isBuilt = false;
  isInitialized = false;
}

/////////////////////////////////////////////////////////////////////////////

void HadrontherapyDetectorROGeometry::Initialize(G4ThreeVector pos, G4double detectorDimX,
                                                 G4double detectorDimY, G4double detectorDimZ,
                                                 G4int numberOfVoxelsX, G4int numberOfVoxelsY,
                                                 G4int numberOfVoxelsZ)
{
  detectorToWorldPosition = pos;
  detectorSizeX = detectorDimX;
  detectorSizeY = detectorDimY;
  detectorSizeZ = detectorDimZ;
  numberOfVoxelsAlongX = numberOfVoxelsX;
  numberOfVoxelsAlongY = numberOfVoxelsY;
  numberOfVoxelsAlongZ = numberOfVoxelsZ;

  isInitialized = true;
}

void HadrontherapyDetectorROGeometry::UpdateROGeometry()
{
  // Nothing happens if the RO geometry is not built. But parameters are properly set.
  if (!isBuilt)
  {
    // G4Exception("HadrontherapyDetectorROGeometry::UpdateROGeometry","had001",
    //		  JustWarning,"Cannot update geometry before it is built");
    return;
  }

  // 1) Update the dimensions of the G4Boxes
  G4double halfDetectorSizeX = detectorSizeX;
  G4double halfDetectorSizeY = detectorSizeY;
  G4double halfDetectorSizeZ = detectorSizeZ;

  RODetector->SetXHalfLength(halfDetectorSizeX);
  RODetector->SetYHalfLength(halfDetectorSizeY);
  RODetector->SetZHalfLength(halfDetectorSizeZ);

  G4double halfXVoxelSizeX = halfDetectorSizeX / numberOfVoxelsAlongX;
  G4double halfXVoxelSizeY = halfDetectorSizeY;
  G4double halfXVoxelSizeZ = halfDetectorSizeZ;
  G4double voxelXThickness = 2 * halfXVoxelSizeX;

  RODetectorXDivision->SetXHalfLength(halfXVoxelSizeX);
  RODetectorXDivision->SetYHalfLength(halfXVoxelSizeY);
  RODetectorXDivision->SetZHalfLength(halfXVoxelSizeZ);

  G4double halfYVoxelSizeX = halfXVoxelSizeX;
  G4double halfYVoxelSizeY = halfDetectorSizeY / numberOfVoxelsAlongY;
  G4double halfYVoxelSizeZ = halfDetectorSizeZ;
  G4double voxelYThickness = 2 * halfYVoxelSizeY;

  RODetectorYDivision->SetXHalfLength(halfYVoxelSizeX);
  RODetectorYDivision->SetYHalfLength(halfYVoxelSizeY);
  RODetectorYDivision->SetZHalfLength(halfYVoxelSizeZ);

  G4double halfZVoxelSizeX = halfXVoxelSizeX;
  G4double halfZVoxelSizeY = halfYVoxelSizeY;
  G4double halfZVoxelSizeZ = halfDetectorSizeZ / numberOfVoxelsAlongZ;
  G4double voxelZThickness = 2 * halfZVoxelSizeZ;

  RODetectorZDivision->SetXHalfLength(halfZVoxelSizeX);
  RODetectorZDivision->SetYHalfLength(halfZVoxelSizeY);
  RODetectorZDivision->SetZHalfLength(halfZVoxelSizeZ);

  // Delete and re-build the relevant physical volumes
  G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();

  // Delete...
  G4VPhysicalVolume* myVol = store->GetVolume("RODetectorPhys");
  store->DeRegister(myVol);
  //..and rebuild
  G4VPhysicalVolume* RODetectorPhys = new G4PVPlacement(0, detectorToWorldPosition, RODetectorLV,
                                                        "RODetectorPhys", worldLogical, false, 0);

  myVol = store->GetVolume("RODetectorXDivisionPhys");
  store->DeRegister(myVol);
  G4VPhysicalVolume* RODetectorXDivisionPhys =
    new G4PVReplica("RODetectorXDivisionPhys", RODetectorXDivisionLV, RODetectorPhys, kXAxis,
                    numberOfVoxelsAlongX, voxelXThickness);
  myVol = store->GetVolume("RODetectorYDivisionPhys");
  store->DeRegister(myVol);
  G4VPhysicalVolume* RODetectorYDivisionPhys =
    new G4PVReplica("RODetectorYDivisionPhys", RODetectorYDivisionLV, RODetectorXDivisionPhys,
                    kYAxis, numberOfVoxelsAlongY, voxelYThickness);

  myVol = store->GetVolume("RODetectorZDivisionPhys");
  store->DeRegister(myVol);
  new G4PVReplica("RODetectorZDivisionPhys", RODetectorZDivisionLV, RODetectorYDivisionPhys, kZAxis,
                  numberOfVoxelsAlongZ, voxelZThickness);

  // G4Orb* solidEnvelope = dynamic_cast<G4Orb*>(G4SolidStore::GetInstance()->GetSolid("EnvelopeSphere", false));
  // if (solidEnvelope) {
  //     solidEnvelope->SetRadius(HadrontherapyMatrix::envelopeRadius);
  // }

  // G4Orb* solidSphere = dynamic_cast<G4Orb*>(G4SolidStore::GetInstance()->GetSolid("SensitiveSphere", false));
  // if (solidSphere) {
  //     solidSphere->SetRadius(HadrontherapyMatrix::sphereRadius);
  // }

  return;
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorROGeometry::~HadrontherapyDetectorROGeometry()
{
  ;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorROGeometry::Construct()
{
  // A dummy material is used to fill the volumes of the readout geometry.
  // (It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)

  //
  // World
  //
  G4VPhysicalVolume* ghostWorld = GetWorld();
  worldLogical = ghostWorld->GetLogicalVolume();

  if (!isInitialized)
  {
    G4Exception("HadrontherapyDetectorROGeometry::Construct", "had001", FatalException,
                "Parameters of the RO geometry are not initialized");
    return;
  }

  G4double halfDetectorSizeX = detectorSizeX;
  G4double halfDetectorSizeY = detectorSizeY;
  G4double halfDetectorSizeZ = detectorSizeZ;

  // World volume of ROGeometry ... SERVE SOLO PER LA ROG

  // Detector ROGeometry
  RODetector = new G4Box("RODetector", halfDetectorSizeX, halfDetectorSizeY, halfDetectorSizeZ);

  RODetectorLV = new G4LogicalVolume(RODetector, 0, "RODetectorLV", 0, 0, 0);

  G4VPhysicalVolume* RODetectorPhys = new G4PVPlacement(0, detectorToWorldPosition, RODetectorLV,
                                                        "RODetectorPhys", worldLogical, false, 0);

  // Division along X axis: the detector is divided in slices along the X axis

  G4double halfXVoxelSizeX = halfDetectorSizeX / numberOfVoxelsAlongX;
  G4double halfXVoxelSizeY = halfDetectorSizeY;
  G4double halfXVoxelSizeZ = halfDetectorSizeZ;
  G4double voxelXThickness = 2 * halfXVoxelSizeX;

  RODetectorXDivision =
    new G4Box("RODetectorXDivision", halfXVoxelSizeX, halfXVoxelSizeY, halfXVoxelSizeZ);

  RODetectorXDivisionLV =
    new G4LogicalVolume(RODetectorXDivision, 0, "RODetectorXDivisionLV", 0, 0, 0);

  G4VPhysicalVolume* RODetectorXDivisionPhys =
    new G4PVReplica("RODetectorXDivisionPhys", RODetectorXDivisionLV, RODetectorPhys, kXAxis,
                    numberOfVoxelsAlongX, voxelXThickness);

  // Division along Y axis: the slices along the X axis are divided along the Y axis

  G4double halfYVoxelSizeX = halfXVoxelSizeX;
  G4double halfYVoxelSizeY = halfDetectorSizeY / numberOfVoxelsAlongY;
  G4double halfYVoxelSizeZ = halfDetectorSizeZ;
  G4double voxelYThickness = 2 * halfYVoxelSizeY;

  RODetectorYDivision =
    new G4Box("RODetectorYDivision", halfYVoxelSizeX, halfYVoxelSizeY, halfYVoxelSizeZ);

  RODetectorYDivisionLV =
    new G4LogicalVolume(RODetectorYDivision, 0, "RODetectorYDivisionLV", 0, 0, 0);

  G4VPhysicalVolume* RODetectorYDivisionPhys =
    new G4PVReplica("RODetectorYDivisionPhys", RODetectorYDivisionLV, RODetectorXDivisionPhys,
                    kYAxis, numberOfVoxelsAlongY, voxelYThickness);

  // Division along Z axis: the slices along the Y axis are divided along the Z axis

  G4double halfZVoxelSizeX = halfXVoxelSizeX;
  G4double halfZVoxelSizeY = halfYVoxelSizeY;
  G4double halfZVoxelSizeZ = halfDetectorSizeZ / numberOfVoxelsAlongZ;
  G4double voxelZThickness = 2 * halfZVoxelSizeZ;

  RODetectorZDivision =
    new G4Box("RODetectorZDivision", halfZVoxelSizeX, halfZVoxelSizeY, halfZVoxelSizeZ);

  RODetectorZDivisionLV =
    new G4LogicalVolume(RODetectorZDivision, 0, "RODetectorZDivisionLV", 0, 0, 0);

  new G4PVReplica("RODetectorZDivisionPhys", RODetectorZDivisionLV, RODetectorYDivisionPhys, kZAxis,
                  numberOfVoxelsAlongZ, voxelZThickness);

  // // ############################################################
  // // construct sphere and envelope for lineal energy
  // // ############################################################

  // G4double sphereRadius = HadrontherapyMatrix::sphereRadius;
  // // 새로 추가: Envelope 반경 (Matrix 클래스에서 받아옴)
  // G4double envelopeRadius = HadrontherapyMatrix::envelopeRadius; 

  // // 1. Envelope 체적 생성 (가장 바깥쪽 버퍼)
  // G4Orb* solidEnvelope = new G4Orb("EnvelopeSphere", envelopeRadius);
  // G4LogicalVolume* envelopeSphereLV = new G4LogicalVolume(solidEnvelope, 
  //                                                         0, // Readout geometry 더미 물질
  //                                                         "EnvelopeSphereLV");

  // // Envelope을 Voxel(RODetectorZDivisionLV)의 중앙에 배치
  // new G4PVPlacement(0, G4ThreeVector(0, 0, 0), envelopeSphereLV, "EnvelopeSpherePhys", 
  //                   RODetectorZDivisionLV, false, 0);

  // // 2. Sensitive Sphere 체적 생성 (실제 센서)
  // G4Orb* solidSphere = new G4Orb("SensitiveSphere", sphereRadius);
  // sensitiveSphereLV = new G4LogicalVolume(solidSphere, 
  //                                         0, // Readout geometry 더미 물질
  //                                         "SensitiveSphereLV");

  // // Sensitive Sphere를 Voxel이 아닌 'Envelope(envelopeSphereLV)' 내부에 배치 (계층화)
  // new G4PVPlacement(0, G4ThreeVector(0, 0, 0), sensitiveSphereLV, "SensitiveSpherePhys", 
  //                   envelopeSphereLV, false, 0);

  // // 3. 나중에 Cut 값을 따로 주기 위해 Envelope을 G4Region으로 등록
  // G4Region* envelopeRegion = G4RegionStore::GetInstance()->GetRegion("EnvelopeLog", false);
  // if (!envelopeRegion) {
  //     envelopeRegion = new G4Region("EnvelopeLog");
  // }
  // envelopeSphereLV->SetRegion(envelopeRegion);
  // envelopeRegion->AddRootLogicalVolume(envelopeSphereLV);

  isBuilt = true;
}

void HadrontherapyDetectorROGeometry::ConstructSD()
{
  G4String sensitiveDetectorName = "RODetector";
  HadrontherapyDetectorSD* detectorSD = new HadrontherapyDetectorSD(sensitiveDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(detectorSD);

  // 🌟 [핵심 수정] 🌟
  // 매크로 실행 전이라 모드를 알 수 없으므로, 두 체적에 무조건 센서를 모두 붙입니다!
  // 어차피 빔을 쏠 때 ProcessHits 함수가 현재 모드를 보고 알아서 차단해 줄 것입니다.
  RODetectorZDivisionLV->SetSensitiveDetector(detectorSD);

  G4cout << "===> ROGeometry: Sensitive Detector attached to Voxels." << G4endl;

  // 🌟 [수정 2] 질량 세계의 Sphere를 가져와서 SD 부착 (Lineal 모드용)
  G4LogicalVolume* sensSphereLV = G4LogicalVolumeStore::GetInstance()->GetVolume("SensitiveSphereLV", false);
  if(sensSphereLV) {
      sensSphereLV->SetSensitiveDetector(detectorSD);
      G4cout << "===> ROGeometry: SD safely attached to Mass World Sphere via Store!" << G4endl;
  }
}