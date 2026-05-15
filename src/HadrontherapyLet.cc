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

#include "HadrontherapyLet.hh"

#include "G4AnalysisManager.hh"
#include "G4AutoLock.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyInteractionParameters.hh"
#include "HadrontherapyMatrix.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"

#include <cmath>
#include <map>
#include <set>

namespace
{
G4Mutex letMutex = G4MUTEX_INITIALIZER;

struct LocalLetData
{
    G4double dn{0.}, dd{0.}, tn{0.}, td{0.};
};

// [수정] 로컬 맵 최적화: 4개의 맵을 구조체 1개로 통합
thread_local std::unordered_map<G4int, LocalLetData> locTotalLet;
thread_local std::map<G4int, std::unordered_map<G4int, LocalLetData>> locIonLet;
thread_local std::map<G4int, G4ParticleDefinition*> locParticleDefLet;
thread_local std::map<G4int, G4bool> locIsPrimaryLet;
}  // namespace

HadrontherapyLet* HadrontherapyLet::instance = NULL;
G4bool HadrontherapyLet::CalculationDone = false;

HadrontherapyLet* HadrontherapyLet::GetInstance(HadrontherapyDetectorConstruction* pDet)
{
  if (instance) delete instance;
  instance = new HadrontherapyLet(pDet);
  return instance;
}

HadrontherapyLet* HadrontherapyLet::GetInstance()
{
  return instance;
}

HadrontherapyLet::HadrontherapyLet(HadrontherapyDetectorConstruction* pDet)
  : matrix(0)  // Default output filename
{
  matrix = HadrontherapyMatrix::GetInstance();
  if (!matrix)
    G4Exception("HadrontherapyLet::HadrontherapyLet", "Hadrontherapy0005", FatalException,
                "HadrontherapyMatrix not found. Firstly create an instance of it.");

  nVoxels = matrix->GetNvoxel();
  numberOfVoxelAlongX = matrix->GetNumberOfVoxelAlongX();
  numberOfVoxelAlongY = matrix->GetNumberOfVoxelAlongY();
  numberOfVoxelAlongZ = matrix->GetNumberOfVoxelAlongZ();

  detectorMat = pDet->GetDetectorLogicalVolume()->GetMaterial();
  density = detectorMat->GetDensity();

  // [수정] 거대 포인터 배열들 모두 삭제됨
}

HadrontherapyLet::~HadrontherapyLet()
{
  Clear();
}

// Fill energy spectrum for every voxel (local energy spectrum)
void HadrontherapyLet::Initialize()
{
  totalLetStats.clear();
  Clear();
}
/**
 * Clear all stored data
 */
void HadrontherapyLet::Clear()
{
  ionLetStore.clear();
}
void HadrontherapyLet::FillEnergySpectrum(G4int trackID, G4ParticleDefinition* particleDef,
                                          G4double ekinMean, const G4Material* mat, G4double DE,
                                          G4double DEEletrons, G4double DX, G4int i, G4int j,
                                          G4int k)
{
  if (DE <= 0. || DX <= 0.) return;
  if (!CalculationDone) return;

  G4int Z = particleDef->GetAtomicNumber();
  if (Z < 1) return;

  G4int PDGencoding = particleDef->GetPDGEncoding();
  PDGencoding -= PDGencoding % 10;

  G4int voxel = matrix->Index(i, j, k);

  G4EmCalculator emCal;
  G4double Lsn = emCal.ComputeElectronicDEDX(ekinMean, particleDef, mat);

  G4double dN = (DE + DEEletrons) * Lsn;
  G4double dD = DE + DEEletrons;
  G4double tN = DX * Lsn;
  G4double tD = DX;

  // 로컬 버퍼 누적
  locTotalLet[voxel].dn += dN;
  locTotalLet[voxel].dd += dD;
  locTotalLet[voxel].tn += tN;
  locTotalLet[voxel].td += tD;

  locIonLet[PDGencoding][voxel].dn += dN;
  locIonLet[PDGencoding][voxel].dd += dD;
  locIonLet[PDGencoding][voxel].tn += tN;
  locIonLet[PDGencoding][voxel].td += tD;

  locParticleDefLet[PDGencoding] = particleDef;
  locIsPrimaryLet[PDGencoding] = (trackID == 1);
}

void HadrontherapyLet::StoreLETEventData()
{
  {
    G4AutoLock l(&letMutex);

    // 1. Total LET Flush
    for (const auto& pair : locTotalLet)
    {
      G4int v = pair.first;
      const LocalLetData& d = pair.second;
      if (d.dd > 0. || d.dn > 0. || d.td > 0. || d.tn > 0.)
      {
        totalLetStats[v].fill(d.dn, d.dd, d.tn, d.td);
      }
    }

    // 2. Ion LET Flush
    std::set<G4int> activePDGs;
    for (const auto& pair : locIonLet)
      activePDGs.insert(pair.first);

    for (G4int pdg : activePDGs)
    {
      ionLet* targetIon = nullptr;
      for (size_t idx = 0; idx < ionLetStore.size(); idx++)
      {
        if (ionLetStore[idx].PDGencoding == pdg)
        {
          targetIon = &ionLetStore[idx];
          break;
        }
      }

      if (!targetIon)
      {
        G4ParticleDefinition* pDef = locParticleDefLet[pdg];
        G4String fullName = pDef->GetParticleName();
        G4String name = fullName.substr(0, fullName.find("["));

        ionLet newIon;
        newIon.isPrimary = locIsPrimaryLet[pdg];
        newIon.PDGencoding = pdg;
        newIon.fullName = fullName;
        newIon.name = name;
        newIon.Z = pDef->GetAtomicNumber();
        newIon.A = pDef->GetAtomicMass();

        ionLetStore.push_back(newIon);
        targetIon = &ionLetStore.back();
      }

      for (const auto& pair : locIonLet[pdg])
      {
        G4int v = pair.first;
        const LocalLetData& d = pair.second;
        if (d.dd > 0. || d.dn > 0. || d.td > 0. || d.tn > 0.)
        {
          targetIon->sparseLet[v].fill(d.dn, d.dd, d.tn, d.td);
        }
      }
    }
  }  // <--- Lock 해제

  locTotalLet.clear();
  locIonLet.clear();
  locParticleDefLet.clear();
  locIsPrimaryLet.clear();
}

void HadrontherapyLet::LetOutput()
{
  G4int nEvts = 1;
  if (matrix) nEvts = matrix->GetTotalEvents();
  G4double N = (G4double)nEvts;

  // 1. Total LET 최종 계산
  for (auto& pair : totalLetStats)
  {
    LetStat& stat = pair.second;

    if (stat.dD > 0.)
    {
      if (N > 1.0)
      {
        G4double R = stat.dN / stat.dD;
        G4double s2_Y = (stat.dN2 - (stat.dN * stat.dN) / N) / (N - 1.0);
        G4double s2_X = (stat.dD2 - (stat.dD * stat.dD) / N) / (N - 1.0);
        G4double s_XY = (stat.dCov - (stat.dN * stat.dD) / N) / (N - 1.0);
        G4double varR = (N / (stat.dD * stat.dD)) * (s2_Y + R * R * s2_X - 2.0 * R * s_XY);
        if (varR > 0.) stat.letD_SE = std::sqrt(varR);
      }
      stat.letD = stat.dN / stat.dD;
    }

    if (stat.tD > 0.)
    {
      if (N > 1.0)
      {
        G4double R = stat.tN / stat.tD;
        G4double s2_Y = (stat.tN2 - (stat.tN * stat.tN) / N) / (N - 1.0);
        G4double s2_X = (stat.tD2 - (stat.tD * stat.tD) / N) / (N - 1.0);
        G4double s_XY = (stat.tCov - (stat.tN * stat.tD) / N) / (N - 1.0);
        G4double varR = (N / (stat.tD * stat.tD)) * (s2_Y + R * R * s2_X - 2.0 * R * s_XY);
        if (varR > 0.) stat.letT_SE = std::sqrt(varR);
      }
      stat.letT = stat.tN / stat.tD;
    }
  }

  std::sort(ionLetStore.begin(), ionLetStore.end());

  // 2. Ion LET 최종 계산
  for (size_t ion = 0; ion < ionLetStore.size(); ion++)
  {
    for (auto& pair : ionLetStore[ion].sparseLet)
    {
      LetStat& stat = pair.second;

      if (stat.dD > 0.)
      {
        if (N > 1.0)
        {
          G4double R = stat.dN / stat.dD;
          G4double s2_Y = (stat.dN2 - (stat.dN * stat.dN) / N) / (N - 1.0);
          G4double s2_X = (stat.dD2 - (stat.dD * stat.dD) / N) / (N - 1.0);
          G4double s_XY = (stat.dCov - (stat.dN * stat.dD) / N) / (N - 1.0);
          G4double varR = (N / (stat.dD * stat.dD)) * (s2_Y + R * R * s2_X - 2.0 * R * s_XY);
          if (varR > 0.) stat.letD_SE = std::sqrt(varR);
        }
        stat.letD = stat.dN / stat.dD;
      }

      if (stat.tD > 0.)
      {
        if (N > 1.0)
        {
          G4double R = stat.tN / stat.tD;
          G4double s2_Y = (stat.tN2 - (stat.tN * stat.tN) / N) / (N - 1.0);
          G4double s2_X = (stat.tD2 - (stat.tD * stat.tD) / N) / (N - 1.0);
          G4double s_XY = (stat.tCov - (stat.tN * stat.tD) / N) / (N - 1.0);
          G4double varR = (N / (stat.tD * stat.tD)) * (s2_Y + R * R * s2_X - 2.0 * R * s_XY);
          if (varR > 0.) stat.letT_SE = std::sqrt(varR);
        }
        stat.letT = stat.tN / stat.tD;
      }
    }
  }
}

void HadrontherapyLet::StoreLetBinary(G4String filename)
{
  // if (ionLetStore.size() == 0) return;
  
  // [수정 1] HadrontherapyMatrix로부터 출력 포맷 정보 가져오기
  G4String outputFormat = "binary";
  if (matrix) {
      outputFormat = matrix->GetOutputFormat(); 
  }

  // [수정 2] 포맷에 따른 파일 확장자 및 열기 모드 변경
  if (outputFormat == "csv") {
    ofs.open(filename + ".csv", std::ios::out);
  } else {
    ofs.open(filename + ".out", std::ios::out | std::ios::binary);
  }

  if (ofs.is_open())
  {
    G4int numIons = ionLetStore.size();

    // [수정 3] 헤더 작성 분기 처리
    if (outputFormat == "csv") {
      ofs << "i,j,k,tDoseLet,tDoseLetSE,tTrackLet,tTrackLetSE";
      for (size_t l = 0; l < numIons; l++) {
        G4String a = (ionLetStore[l].isPrimary) ? "_1" : "";
        G4String ionName = ionLetStore[l].name + a;
        ofs << "," << ionName << "_DoseLet," << ionName << "_DoseLetSE"
            << "," << ionName << "_TrackLet," << ionName << "_TrackLetSE";
      }
      ofs << "\n";
    } else {
      ofs.write(reinterpret_cast<const char*>(&numberOfVoxelAlongX), sizeof(G4int));
      ofs.write(reinterpret_cast<const char*>(&numberOfVoxelAlongY), sizeof(G4int));
      ofs.write(reinterpret_cast<const char*>(&numberOfVoxelAlongZ), sizeof(G4int));
      ofs.write(reinterpret_cast<const char*>(&numIons), sizeof(G4int));

      for (size_t l = 0; l < numIons; l++)
      {
        G4String a = (ionLetStore[l].isPrimary) ? "_1" : "";
        G4String ionName = ionLetStore[l].name + a;
        G4int nameLen = ionName.length();
        ofs.write(reinterpret_cast<const char*>(&nameLen), sizeof(G4int));
        ofs.write(ionName.c_str(), nameLen);
      }
    }

    // 2. 본문 저장
    G4double conv = 1.0 / (keV / um);

    for (G4int i = 0; i < numberOfVoxelAlongX; i++)
      for (G4int j = 0; j < numberOfVoxelAlongY; j++)
        for (G4int k = 0; k < numberOfVoxelAlongZ; k++)
        {
          G4int v = matrix->Index(i, j, k);

          // Total LET 값이 기록된 복셀만 저장
          if (totalLetStats.find(v) != totalLetStats.end())
          {
            const LetStat& tStat = totalLetStats[v];

            G4double tDoseLet = tStat.letD * conv;
            G4double tDoseLetSE = tStat.letD_SE * conv;
            G4double tTrackLet = tStat.letT * conv;
            G4double tTrackLetSE = tStat.letT_SE * conv;

            // [수정 4] 좌표와 Total LET 데이터 쓰기 분기
            if (outputFormat == "csv") {
              ofs << i << "," << j << "," << k << "," 
                  << tDoseLet << "," << tDoseLetSE << "," 
                  << tTrackLet << "," << tTrackLetSE;
            } else {
              ofs.write(reinterpret_cast<const char*>(&i), sizeof(G4int));
              ofs.write(reinterpret_cast<const char*>(&j), sizeof(G4int));
              ofs.write(reinterpret_cast<const char*>(&k), sizeof(G4int));
              ofs.write(reinterpret_cast<const char*>(&tDoseLet), sizeof(G4double));
              ofs.write(reinterpret_cast<const char*>(&tDoseLetSE), sizeof(G4double));
              ofs.write(reinterpret_cast<const char*>(&tTrackLet), sizeof(G4double));
              ofs.write(reinterpret_cast<const char*>(&tTrackLetSE), sizeof(G4double));
            }

            for (size_t l = 0; l < numIons; l++)
            {
              G4double iDoseLet = 0., iDoseLetSE = 0., iTrackLet = 0., iTrackLetSE = 0.;

              if (ionLetStore[l].sparseLet.find(v) != ionLetStore[l].sparseLet.end())
              {
                const LetStat& iStat = ionLetStore[l].sparseLet[v];
                iDoseLet = iStat.letD * conv;
                iDoseLetSE = iStat.letD_SE * conv;
                iTrackLet = iStat.letT * conv;
                iTrackLetSE = iStat.letT_SE * conv;
              }

              // [수정 5] Ion 별 LET 데이터 쓰기 분기
              if (outputFormat == "csv") {
                ofs << "," << iDoseLet << "," << iDoseLetSE << "," 
                    << iTrackLet << "," << iTrackLetSE;
              } else {
                ofs.write(reinterpret_cast<const char*>(&iDoseLet), sizeof(G4double));
                ofs.write(reinterpret_cast<const char*>(&iDoseLetSE), sizeof(G4double));
                ofs.write(reinterpret_cast<const char*>(&iTrackLet), sizeof(G4double));
                ofs.write(reinterpret_cast<const char*>(&iTrackLetSE), sizeof(G4double));
              }
            }
            // CSV일 경우 줄바꿈 추가
            if (outputFormat == "csv") {
                ofs << "\n";
            }
          }
        }
    ofs.close();
  }
}