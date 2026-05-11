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

#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyDetectorROGeometry.hh"
#include "HadrontherapyDetectorMessenger.hh"
#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyMatrix.hh"
#include "HadrontherapyLet.hh"
#include "PassiveProtonBeamLine.hh"
#include "BESTPassiveProtonBeamLine.hh"
#include "HadrontherapyMatrix.hh"

#include "HadrontherapyRBE.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

#include "G4Orb.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PVReplica.hh"


HadrontherapyDetectorConstruction* HadrontherapyDetectorConstruction::instance = 0;
/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorConstruction::HadrontherapyDetectorConstruction(G4VPhysicalVolume* physicalTreatmentRoom)
: motherPhys(physicalTreatmentRoom), // pointer to WORLD volume
detectorSD(0), detectorROGeometry(0), matrix(0),
phantom(0), detector(0),
phantomLogicalVolume(0), detectorLogicalVolume(0),
phantomPhysicalVolume(0), detectorPhysicalVolume(0),
aRegion(0)
{
    
    
    /* NOTE! that the HadrontherapyDetectorConstruction class
     * does NOT inherit from G4VUserDetectorConstruction G4 class
     * So the Construct() mandatory virtual method is inside another geometric class
     * like the passiveProtonBeamLIne, ...
     */
    instance = this;
    // Messenger to change parameters of the phantom/detector geometry
    detectorMessenger = new HadrontherapyDetectorMessenger(this);
    
    // Default detector voxels size
    // slabs 
    sizeOfVoxelAlongX = 1 *cm;
    sizeOfVoxelAlongY = 1 *cm;
    sizeOfVoxelAlongZ = 1 *cm;
    
    // Define here the material of the water phantom and of the detector
    SetPhantomMaterial("G4_WATER");
    //SetPhantomMaterial("G4_AIR");
    
    // Construct geometry (messenger commands)
    // SetDetectorSize(4.*cm, 4.*cm, 4.*cm);
    SetDetectorSize(4. *cm, 4. *cm, 4. *cm);
    SetPhantomSize(40. *cm, 40. *cm, 40. *cm);
    
    SetPhantomPosition(G4ThreeVector(20. *cm, 0. *cm, 0. *cm));
    SetDetectorToPhantomPosition(G4ThreeVector(0. *cm, 0. *cm, 0. *cm));
    SetDetectorPosition();
    //GetDetectorToWorldPosition();
    
    // Write virtual parameters to the real ones and check for consistency
    UpdateGeometry();
    
    
    
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorConstruction::~HadrontherapyDetectorConstruction()
{
    delete detectorROGeometry;
    delete matrix;
    delete detectorMessenger;
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorConstruction* HadrontherapyDetectorConstruction::GetInstance()
{
    return instance;
}

/////////////////////////////////////////////////////////////////////////////
// ConstructPhantom() is the method that construct a water box (called phantom
// (or water phantom) in the usual Medical physicists slang).
// A water phantom can be considered a good approximation of a an human body.
void HadrontherapyDetectorConstruction::ConstructPhantom()
{
    // Definition of the solid volume of the Phantom
    phantom = new G4Box("Phantom",
                        phantomSizeX/2,
                        phantomSizeY/2,
                        phantomSizeZ/2);
    
    // Definition of the logical volume of the Phantom
    phantomLogicalVolume = new G4LogicalVolume(phantom,
                                               phantomMaterial,
                                               "phantomLog", 0, 0, 0);
    
    // Definition of the physics volume of the Phantom
    phantomPhysicalVolume = new G4PVPlacement(0,
                                              phantomPosition,
                                              "phantomPhys",
                                              phantomLogicalVolume,
                                              motherPhys,
                                              false,
                                              0);
    
    // Visualisation attributes of the phantom
    red = new G4VisAttributes(G4Colour(255/255., 0/255. ,0/255.));
    red -> SetVisibility(true);
    red -> SetForceSolid(true);
    red -> SetForceWireframe(true);
    phantomLogicalVolume -> SetVisAttributes(red);
}

/////////////////////////////////////////////////////////////////////////////
// ConstructDetector() is the method the reconstruct a detector region
// inside the water phantom. It is a volume, located inside the water phantom.
//
//           **************************
//           *    water phantom       *
//           *                        *
//           *                        *
//           *---------------         *
//  Beam     *              -         *
//  ----->   * detector     -         *
//           *              -         *
//           *---------------         *
//           *                        *
//           *                        *
//           *                        *
//           **************************
//
// The detector can be dived in slices or voxelized
// and inside it different quantities (dose distribution, fluence distribution, LET, etc)
// can be stored.
void HadrontherapyDetectorConstruction::ConstructDetector()

{
    // Definition of the solid volume of the Detector
    detector = new G4Box("Detector",
                         
                         detectorSizeX/2,
                         
                         detectorSizeY/2,
                         
                         detectorSizeZ/2);
    
    // Definition of the logic volume of the Phantom
    detectorLogicalVolume = new G4LogicalVolume(detector,
                                                detectorMaterial,
                                                "DetectorLog",
                                                0,0,0);
    // Definition of the physical volume of the Phantom
    detectorPhysicalVolume = new G4PVPlacement(0,
                                               detectorPosition, // Setted by displacement
                                               "DetectorPhys",
                                               detectorLogicalVolume,
                                               phantomPhysicalVolume,
                                               false,0);
    
    // Visualisation attributes of the detector
    skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
    skyBlue -> SetVisibility(true);
    skyBlue -> SetForceSolid(true);
    //skyBlue -> SetForceWireframe(true);
    detectorLogicalVolume -> SetVisAttributes(skyBlue);

    // ************************************************************
    // G4PVReplica를 이용한 3단 격자 분할 및 구슬 배치
    // ************************************************************
    // Voxel 크기의 절반값 계산
    G4double halfXVoxelSize = sizeOfVoxelAlongX / 2.0;
    G4double halfYVoxelSize = sizeOfVoxelAlongY / 2.0;
    G4double halfZVoxelSize = sizeOfVoxelAlongZ / 2.0;

    // (1) X축 분할
    G4Box* repX = new G4Box("RepX", halfXVoxelSize, detectorSizeY / 2.0, detectorSizeZ / 2.0);
    G4LogicalVolume* repXLV = new G4LogicalVolume(repX, detectorMaterial, "RepXLV");
    new G4PVReplica("RepXPhys", repXLV, detectorLogicalVolume, kXAxis, numberOfVoxelsAlongX, sizeOfVoxelAlongX);

    // (2) Y축 분할 (X 분할된 것을 다시 Y로 분할)
    G4Box* repY = new G4Box("RepY", halfXVoxelSize, halfYVoxelSize, detectorSizeZ / 2.0);
    G4LogicalVolume* repYLV = new G4LogicalVolume(repY, detectorMaterial, "RepYLV");
    new G4PVReplica("RepYPhys", repYLV, repXLV, kYAxis, numberOfVoxelsAlongY, sizeOfVoxelAlongY);

    // (3) Z축 분할 (최종 단위 복셀 VoxelZ 완성)
    G4Box* voxelZ = new G4Box("VoxelZ", halfXVoxelSize, halfYVoxelSize, halfZVoxelSize);
    G4LogicalVolume* voxelZLV = new G4LogicalVolume(voxelZ, detectorMaterial, "VoxelZLV");
    new G4PVReplica("VoxelZPhys", voxelZLV, repYLV, kZAxis, numberOfVoxelsAlongZ, sizeOfVoxelAlongZ);

    // ------------------------------------------------------------
    // 단위 복셀(VoxelZ) 중심에 Envelope & Sphere 배치 (for문 없이 단 1번 배치)
    // ------------------------------------------------------------
    G4double sphereRadius = HadrontherapyMatrix::sphereRadius;
    G4double envelopeRadius = HadrontherapyMatrix::envelopeRadius; 

    // Envelope 체적 생성 (가장 바깥쪽 버퍼)
    G4Orb* solidEnvelope = new G4Orb("EnvelopeSphere", envelopeRadius);
    envelopeSphereLV = new G4LogicalVolume(solidEnvelope, detectorMaterial, "EnvelopeSphereLV");

    // Sensitive Sphere 체적 생성 (실제 센서)
    G4Orb* solidSphere = new G4Orb("SensitiveSphere", sphereRadius);
    sensitiveSphereLV = new G4LogicalVolume(solidSphere, detectorMaterial, "SensitiveSphereLV");

    // Voxel 중앙에 Envelope 배치
    new G4PVPlacement(0, G4ThreeVector(0,0,0), envelopeSphereLV, "EnvelopeSpherePhys", voxelZLV, false, 0);

    // Envelope 중앙에 Sensitive Sphere 배치
    new G4PVPlacement(0, G4ThreeVector(0,0,0), sensitiveSphereLV, "SensitiveSpherePhys", envelopeSphereLV, false, 0);
    // ************************************************************


    // **************
    // Cut per Region
    // **************
    if (!aRegion)
    {
        aRegion = new G4Region("DetectorLog");
        detectorLogicalVolume -> SetRegion(aRegion);
        aRegion->AddRootLogicalVolume( detectorLogicalVolume );
    }

    // 🌟 [추가] Envelope에 미세 Cut을 적용하기 위한 별도 Region 등록
    G4Region* envelopeRegion = G4RegionStore::GetInstance()->GetRegion("EnvelopeLog", false);
    if (!envelopeRegion) {
        envelopeRegion = new G4Region("EnvelopeLog");
    }
    envelopeSphereLV->SetRegion(envelopeRegion);
    envelopeRegion->AddRootLogicalVolume(envelopeSphereLV);
}

///////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::InitializeDetectorROGeometry(
                                                                     HadrontherapyDetectorROGeometry* RO,
                                                                     G4ThreeVector
                                                                     detectorToWorldPosition)
{
    RO->Initialize(detectorToWorldPosition,
                   detectorSizeX/2,
                   detectorSizeY/2,
                   detectorSizeZ/2,
                   numberOfVoxelsAlongX,
                   numberOfVoxelsAlongY,
                   numberOfVoxelsAlongZ);
}
void HadrontherapyDetectorConstruction::VirtualLayer(G4bool Varbool)
{
   
    //Virtual  plane
    VirtualLayerPosition = G4ThreeVector(0*cm,0*cm,0*cm);
    NewSource= Varbool;
    if(NewSource == true)
    {
       // std::cout<<"trr"<<std::endl;
        G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
        
        solidVirtualLayer = new G4Box("VirtualLayer",
                                      1.*um,
                                      20.*cm,
                                      40.*cm);
        
        logicVirtualLayer = new G4LogicalVolume(
                                                solidVirtualLayer,
                                                airNist,
                                                "VirtualLayer");
        
        physVirtualLayer= new G4PVPlacement(0,VirtualLayerPosition,
                                            "VirtualLayer",
                                            logicVirtualLayer,
                                            motherPhys,
                                            false,
                                            0);
        
        logicVirtualLayer -> SetVisAttributes(skyBlue);
    }
    
    
    
    
}


///////////////////////////////////////////////////////////////////////
void  HadrontherapyDetectorConstruction::ParametersCheck()
{
    // Check phantom/detector sizes & relative position
    if (!IsInside(detectorSizeX,
                  detectorSizeY,
                  detectorSizeZ,
                  phantomSizeX,
                  phantomSizeY,
                  phantomSizeZ,
                  detectorToPhantomPosition
                  ))
        G4Exception("HadrontherapyDetectorConstruction::ParametersCheck()", "Hadrontherapy0001", FatalException, "Error: Detector is not fully inside Phantom!");
    
    // Check Detector sizes respect to the voxel ones
    
    if ( detectorSizeX < sizeOfVoxelAlongX) {
        G4Exception("HadrontherapyDetectorConstruction::ParametersCheck()", "Hadrontherapy0002", FatalException, "Error:  Detector X size must be bigger or equal than that of Voxel X!");
    }
    if ( detectorSizeY < sizeOfVoxelAlongY) {
        G4Exception(" HadrontherapyDetectorConstruction::ParametersCheck()", "Hadrontherapy0003", FatalException, "Error:  Detector Y size must be bigger or equal than that of Voxel Y!");
    }
    if ( detectorSizeZ < sizeOfVoxelAlongZ) {
        G4Exception(" HadrontherapyDetectorConstruction::ParametersCheck()", "Hadrontherapy0004", FatalException, "Error:  Detector Z size must be bigger or equal than that of Voxel Z!");
    }
}

///////////////////////////////////////////////////////////////////////
G4bool HadrontherapyDetectorConstruction::SetPhantomMaterial(G4String material)
{
    
    if (G4Material* pMat = G4NistManager::Instance()->FindOrBuildMaterial(material, false) )
    {
        phantomMaterial  = pMat;
        detectorMaterial = pMat;
        if (detectorLogicalVolume && phantomLogicalVolume)
        {
            detectorLogicalVolume -> SetMaterial(pMat);
            phantomLogicalVolume ->  SetMaterial(pMat);
            
            G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
            G4RunManager::GetRunManager() -> GeometryHasBeenModified();
            G4cout << "The material of Phantom/Detector has been changed to " << material << G4endl;
        }
    }
    else
    {
        G4cout << "WARNING: material \"" << material << "\" doesn't exist in NIST elements/materials"
        " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl;
        G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl;
        return false;
    }
    
    return true;
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) phantomSizeX = sizeX;
    if (sizeY > 0.) phantomSizeY = sizeY;
    if (sizeZ > 0.) phantomSizeZ = sizeZ;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetDetectorSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) {detectorSizeX = sizeX;}
    if (sizeY > 0.) {detectorSizeY = sizeY;}
    if (sizeZ > 0.) {detectorSizeZ = sizeZ;}
    SetVoxelSize(sizeOfVoxelAlongX, sizeOfVoxelAlongY, sizeOfVoxelAlongZ);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetVoxelSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) {sizeOfVoxelAlongX = sizeX;}
    if (sizeY > 0.) {sizeOfVoxelAlongY = sizeY;}
    if (sizeZ > 0.) {sizeOfVoxelAlongZ = sizeZ;}
}

///////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetPhantomPosition(G4ThreeVector pos)
{
    phantomPosition = pos;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetDetectorToPhantomPosition(G4ThreeVector displ)
{
    detectorToPhantomPosition = displ;
}

void HadrontherapyDetectorConstruction::SetVirtualLayerPosition(G4ThreeVector position)
{
    
    VirtualLayerPosition = position;
    physVirtualLayer->SetTranslation(VirtualLayerPosition);
    
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::UpdateGeometry()
{
    /*
     * Check parameters consistency
     */
    ParametersCheck();
    
    G4GeometryManager::GetInstance() -> OpenGeometry();

    // Round to nearest integer number of voxel
    
    numberOfVoxelsAlongX = G4lrint(detectorSizeX / sizeOfVoxelAlongX);
    sizeOfVoxelAlongX = ( detectorSizeX / numberOfVoxelsAlongX );
    numberOfVoxelsAlongY = G4lrint(detectorSizeY / sizeOfVoxelAlongY);
    sizeOfVoxelAlongY = ( detectorSizeY / numberOfVoxelsAlongY );
    numberOfVoxelsAlongZ = G4lrint(detectorSizeZ / sizeOfVoxelAlongZ);
    sizeOfVoxelAlongZ = ( detectorSizeZ / numberOfVoxelsAlongZ );

    if (phantom)
    {
        phantom -> SetXHalfLength(phantomSizeX/2);
        phantom -> SetYHalfLength(phantomSizeY/2);
        phantom -> SetZHalfLength(phantomSizeZ/2);
        
        phantomPhysicalVolume -> SetTranslation(phantomPosition);
    }
    else   ConstructPhantom();
    
    
    // Get the center of the detector
    SetDetectorPosition();
    if (detector)
    {
        
        detector -> SetXHalfLength(detectorSizeX/2);
        detector -> SetYHalfLength(detectorSizeY/2);
        detector -> SetZHalfLength(detectorSizeZ/2);
        
        detectorPhysicalVolume -> SetTranslation(detectorPosition);
    }
    else    ConstructDetector();
    
    //std::cout<<NewSource<<std::endl;
    /*if(NewSource)
     {
     std::cout<<"via"<<std::endl;
     }*/
    
    
    // std::cout<<"i"<<std::endl;
    // std::cout<<VirtualLayerPosition<<std::endl;
    // physVirtualLayer->SetTranslation(VirtualLayerPosition);
    
    

    // ************************************************************
    // RO 격자에 맞춰 질량 세계에 구슬 배열(Array) 자동 살포 
    // ************************************************************
    G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
    G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();

    // 1. 기존에 만들어둔 분할 볼륨 및 구슬의 물리적 배치 정보(Phys) 삭제
    std::vector<G4String> physNames = {"RepXPhys", "RepYPhys", "VoxelZPhys", "EnvelopeSpherePhys", "SensitiveSpherePhys"};
    for (const auto& name : physNames) {
        G4VPhysicalVolume* pVol = store->GetVolume(name, false);
        if (pVol) {
            store->DeRegister(pVol);
            delete pVol;
        }
    }

    // 2. 만들어두었던 논리 볼륨(LV) 포인터 가져오기
    G4LogicalVolume* repXLV = lvStore->GetVolume("RepXLV", false);
    G4LogicalVolume* repYLV = lvStore->GetVolume("RepYLV", false);
    G4LogicalVolume* voxelZLV = lvStore->GetVolume("VoxelZLV", false);
    
    // 전역 포인터(클래스 멤버 변수)로 가지고 있는 envelopeSphereLV, sensitiveSphereLV 사용
    if (repXLV && repYLV && voxelZLV && envelopeSphereLV && sensitiveSphereLV) {
        
        G4double halfXVoxelSize = sizeOfVoxelAlongX / 2.0;
        G4double halfYVoxelSize = sizeOfVoxelAlongY / 2.0;
        G4double halfZVoxelSize = sizeOfVoxelAlongZ / 2.0;

        // 3. 기하학 Solid(모양)의 크기를 매크로에서 입력한 새로운 값으로 갱신
        dynamic_cast<G4Box*>(repXLV->GetSolid())->SetXHalfLength(halfXVoxelSize);
        dynamic_cast<G4Box*>(repXLV->GetSolid())->SetYHalfLength(detectorSizeY / 2.0);
        dynamic_cast<G4Box*>(repXLV->GetSolid())->SetZHalfLength(detectorSizeZ / 2.0);

        dynamic_cast<G4Box*>(repYLV->GetSolid())->SetXHalfLength(halfXVoxelSize);
        dynamic_cast<G4Box*>(repYLV->GetSolid())->SetYHalfLength(halfYVoxelSize);
        dynamic_cast<G4Box*>(repYLV->GetSolid())->SetZHalfLength(detectorSizeZ / 2.0);

        dynamic_cast<G4Box*>(voxelZLV->GetSolid())->SetXHalfLength(halfXVoxelSize);
        dynamic_cast<G4Box*>(voxelZLV->GetSolid())->SetYHalfLength(halfYVoxelSize);
        dynamic_cast<G4Box*>(voxelZLV->GetSolid())->SetZHalfLength(halfZVoxelSize);

        // 구슬 반경도 새로운 매크로 값으로 갱신
        dynamic_cast<G4Orb*>(envelopeSphereLV->GetSolid())->SetRadius(HadrontherapyMatrix::envelopeRadius);
        dynamic_cast<G4Orb*>(sensitiveSphereLV->GetSolid())->SetRadius(HadrontherapyMatrix::sphereRadius);

        // 4. Replica 수학적 분할 재구축 (1초 컷)
        new G4PVReplica("RepXPhys", repXLV, detectorLogicalVolume, kXAxis, numberOfVoxelsAlongX, sizeOfVoxelAlongX);
        new G4PVReplica("RepYPhys", repYLV, repXLV, kYAxis, numberOfVoxelsAlongY, sizeOfVoxelAlongY);
        new G4PVReplica("VoxelZPhys", voxelZLV, repYLV, kZAxis, numberOfVoxelsAlongZ, sizeOfVoxelAlongZ);

        // 5. 단위 복셀(VoxelZ) 내부에 Envelope과 Sphere 1번 배치 (자동으로 모든 복셀에 복제됨)
        new G4PVPlacement(0, G4ThreeVector(0,0,0), envelopeSphereLV, "EnvelopeSpherePhys", voxelZLV, false, 0);
        new G4PVPlacement(0, G4ThreeVector(0,0,0), sensitiveSphereLV, "SensitiveSpherePhys", envelopeSphereLV, false, 0);
    }
    // ***********************************************************


    PassiveProtonBeamLine *ppbl= (PassiveProtonBeamLine*)
    
    G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    
    HadrontherapyDetectorROGeometry* RO = (HadrontherapyDetectorROGeometry*) ppbl->GetParallelWorld(0);
    
    //Set parameters, either for the Construct() or for the UpdateROGeometry()
    RO->Initialize(GetDetectorToWorldPosition(),
                   detectorSizeX/2,
                   detectorSizeY/2,
                   detectorSizeZ/2,
                   numberOfVoxelsAlongX,
                   numberOfVoxelsAlongY,
                   numberOfVoxelsAlongZ);
    
    //This method below has an effect only if the RO geometry is already built.
    RO->UpdateROGeometry();
    
    
    
    volumeOfVoxel = sizeOfVoxelAlongX * sizeOfVoxelAlongY * sizeOfVoxelAlongZ;
    massOfVoxel = detectorMaterial -> GetDensity() * volumeOfVoxel;
    //  This will clear the existing matrix (together with all data inside it)!
    matrix = HadrontherapyMatrix::GetInstance(numberOfVoxelsAlongX,
                                              numberOfVoxelsAlongY,
                                              numberOfVoxelsAlongZ,
                                              massOfVoxel);
    
    
    // Initialize RBE
    HadrontherapyRBE::CreateInstance(numberOfVoxelsAlongX, numberOfVoxelsAlongY, numberOfVoxelsAlongZ, massOfVoxel);
    
    // Comment out the line below if let calculation is not needed!
    // Initialize LET with energy of primaries and clear data inside
    if ( (let = HadrontherapyLet::GetInstance(this)) )
    {
        HadrontherapyLet::GetInstance() -> Initialize();
    }
    
    
    // Initialize analysis
    // Inform the kernel about the new geometry
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
    
    PrintParameters();
    
    // CheckOverlaps();
}

/////////////////////////////////////////////////////////////////////////////
//Check of the geometry
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::CheckOverlaps()
{
    G4PhysicalVolumeStore* thePVStore = G4PhysicalVolumeStore::GetInstance();
    G4cout << thePVStore->size() << " physical volumes are defined" << G4endl;
    G4bool overlapFlag = false;
    G4int res=1000;
    G4double tol=0.; //tolerance
    for (size_t i=0;i<thePVStore->size();i++)
    {
        //overlapFlag = (*thePVStore)[i]->CheckOverlaps(res,tol,fCheckOverlapsVerbosity) | overlapFlag;
        overlapFlag = (*thePVStore)[i]->CheckOverlaps(res,tol,true) | overlapFlag;    }
    if (overlapFlag)
        G4cout << "Check: there are overlapping volumes" << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::PrintParameters()
{
    
    G4cout << "The (X,Y,Z) dimensions of the phantom are : (" <<
    G4BestUnit( phantom -> GetXHalfLength()*2., "Length") << ',' <<
    G4BestUnit( phantom -> GetYHalfLength()*2., "Length") << ',' <<
    G4BestUnit( phantom -> GetZHalfLength()*2., "Length") << ')' << G4endl;
    
    G4cout << "The (X,Y,Z) dimensions of the detector are : (" <<
    G4BestUnit( detector -> GetXHalfLength()*2., "Length") << ',' <<
    G4BestUnit( detector -> GetYHalfLength()*2., "Length") << ',' <<
    G4BestUnit( detector -> GetZHalfLength()*2., "Length") << ')' << G4endl;
    
    G4cout << "Displacement between Phantom and World is: ";
    G4cout << "DX= "<< G4BestUnit(phantomPosition.getX(),"Length") <<
    "DY= "<< G4BestUnit(phantomPosition.getY(),"Length") <<
    "DZ= "<< G4BestUnit(phantomPosition.getZ(),"Length") << G4endl;
    
    G4cout << "The (X,Y,Z) sizes of the Voxels are: (" <<
    G4BestUnit(sizeOfVoxelAlongX, "Length")  << ',' <<
    G4BestUnit(sizeOfVoxelAlongY, "Length")  << ',' <<
    G4BestUnit(sizeOfVoxelAlongZ, "Length") << ')' << G4endl;
    
    G4cout << "The number of Voxels along (X,Y,Z) is: (" <<
    numberOfVoxelsAlongX  << ',' <<
    numberOfVoxelsAlongY  <<','  <<
    numberOfVoxelsAlongZ  << ')' << G4endl;
}

// void HadrontherapyDetectorConstruction::ConstructSDandField()
// {
//     // 이름이 겹치지 않게 질량 세계 전용 이름을 줍니다.
//     G4String sensitiveDetectorName = "MassDetectorSD"; 
//     HadrontherapyDetectorSD* massDetectorSD = new HadrontherapyDetectorSD(sensitiveDetectorName);
//     G4SDManager::GetSDMpointer()->AddNewDetector(massDetectorSD);

//     // 🌟 G4VUserDetectorConstruction이 제공하는 가장 안전한 SD 부착 함수입니다! 🌟
//     // 워커 스레드 환경에서도 이름만으로 정확히 찾아가서 찰싹 달라붙습니다.
//     SetSensitiveDetector("SensitiveSphereLV", massDetectorSD);

//     G4cout << "===> Mass World: SD perfectly attached to SensitiveSphere!" << G4endl;
// }
