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

#include "HadrontherapyMatrix.hh"

#include "G4AutoLock.hh"
#include "G4ParticleGun.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "globals.hh"

#include "HadrontherapyAnalysisFileMessenger.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapySteppingAction.hh"
#include <time.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

namespace
{
G4Mutex matrixMutex = G4MUTEX_INITIALIZER;  // This is for avoiding multi thread core dumped

// independent individual storage  (no need for Lock)
// PDGencoding -> (VoxelIndex -> value)
thread_local std::map<G4int, std::map<G4int, G4double>> localIonDose;
thread_local std::map<G4int, std::map<G4int, unsigned int>> localIonFluence;

thread_local std::map<G4int, G4ParticleDefinition*> localParticleDef;
thread_local std::map<G4int, G4bool> localIsPrimary;
thread_local std::unordered_map<G4int, G4int> localHitTrack;

}  // namespace

HadrontherapyAnalysis* HadrontherapyAnalysis::instance = 0;
/////////////////////////////////////////////////////////////////////////////

HadrontherapyAnalysis::HadrontherapyAnalysis()
{
  fMess = new HadrontherapyAnalysisFileMessenger(this);
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysis::~HadrontherapyAnalysis()
{
  delete fMess;
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysis* HadrontherapyAnalysis::GetInstance()
{
  if (instance == 0) instance = new HadrontherapyAnalysis;
  return instance;
}

HadrontherapyMatrix* HadrontherapyMatrix::instance = NULL;
G4bool HadrontherapyMatrix::secondary = false;
G4bool HadrontherapyMatrix::storeRawStats = false;
G4double HadrontherapyMatrix::sphereRadius = 5.0 * um;
G4double HadrontherapyMatrix::envelopeRadius = 1.0 * mm;

// Only return a pointer to matrix
HadrontherapyMatrix* HadrontherapyMatrix::GetInstance()
{
  return instance;
}

/////////////////////////////////////////////////////////////////////////////
// This STATIC method delete (!) the old matrix and rewrite a new object returning a pointer to it
// TODO A check on the parameters is required!
HadrontherapyMatrix* HadrontherapyMatrix::GetInstance(G4int voxelX, G4int voxelY, G4int voxelZ,
                                                      G4double mass)
{
  if (instance) delete instance;
  instance = new HadrontherapyMatrix(voxelX, voxelY, voxelZ, mass);
  instance->Initialize();
  return instance;
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyMatrix::HadrontherapyMatrix(G4int voxelX, G4int voxelY, G4int voxelZ, G4double mass)
  : doseUnit(gray), totalEvents(0)
{
  // Number of the voxels of the phantom
  // For Y = Z = 1 the phantom is divided in slices (and not in voxels)
  // orthogonal to the beam axis
  numberOfVoxelAlongX = voxelX;
  numberOfVoxelAlongY = voxelY;
  numberOfVoxelAlongZ = voxelZ;
  massOfVoxel = mass;

  // Create the dose matrix
  matrix = new G4double[numberOfVoxelAlongX * numberOfVoxelAlongY * numberOfVoxelAlongZ];
  matrixSquare = new G4double[numberOfVoxelAlongX * numberOfVoxelAlongY * numberOfVoxelAlongZ];
  matrixCube = new G4double[numberOfVoxelAlongX * numberOfVoxelAlongY * numberOfVoxelAlongZ];
  matrixQuad = new G4double[numberOfVoxelAlongX * numberOfVoxelAlongY * numberOfVoxelAlongZ];
  matrixSumSq = new G4double[numberOfVoxelAlongX * numberOfVoxelAlongY * numberOfVoxelAlongZ];

  for (int i = 0; i < numberOfVoxelAlongX * numberOfVoxelAlongY * numberOfVoxelAlongZ; i++)
  {
    matrix[i] = 0;
    matrixSquare[i] = 0;
    matrixCube[i] = 0;
    matrixQuad[i] = 0;
    matrixSumSq[i] = 0;
  }

  if (matrix)
  {
    G4cout << "HadrontherapyMatrix: Memory space to store physical dose into "
           << numberOfVoxelAlongX * numberOfVoxelAlongY * numberOfVoxelAlongZ
           << " voxels has been allocated " << G4endl;
  }

  else
    G4Exception("HadrontherapyMatrix::HadrontherapyMatrix()", "Hadrontherapy0005", FatalException,
                "Can't allocate memory to store physical dose!");

  // Hit voxel (TrackID) marker
  // This array mark the status of voxel, if a hit occur, with the trackID of the particle
  // Must be initialized

  hitTrack = new G4int[numberOfVoxelAlongX * numberOfVoxelAlongY * numberOfVoxelAlongZ];
  ClearHitTrack();
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyMatrix::~HadrontherapyMatrix()
{
  delete[] matrix;
  delete[] matrixSquare;
  delete[] matrixCube;
  delete[] matrixQuad;
  delete[] matrixSumSq;

  Clear();
}

void HadrontherapyMatrix::Clear()
{
  ionlinStore.clear();
  ionStore.clear();
}

/////////////////////////////////////////////////////////////////////////////
// Initialise the elements of the matrix to zero

void HadrontherapyMatrix::Initialize()
{
  // Clear ions store
  Clear();
  // Clear dose
  for (int i = 0; i < numberOfVoxelAlongX * numberOfVoxelAlongY * numberOfVoxelAlongZ; i++)
  {
    matrix[i] = 0;
    matrixSquare[i] = 0;  //
    matrixCube[i] = 0;
    matrixQuad[i] = 0;
    matrixSumSq[i] = 0;
  }
}

/////////////////////////////////////////////////////////////////////////////
// Print generated nuclides list

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyMatrix::PrintNuclides()
{
  for (size_t i = 0; i < ionlinStore.size(); i++)
  {
    G4cout << ionlinStore[i].name << G4endl;
  }
}

/////////////////////////////////////////////////////////////////////////////
// Clear Hit voxel (TrackID) markers

void HadrontherapyMatrix::ClearHitTrack()
{
  // [수정] 락(Lock)을 걸 필요도 없고, 루프를 돌 필요도 없습니다!
  localHitTrack.clear();
}

G4int* HadrontherapyMatrix::GetHitTrack(G4int i, G4int j, G4int k)
{
  // [수정] 맵은 해당 키가 없으면 0으로 자동 생성 후 그 주소를 반환합니다.
  return &localHitTrack[Index(i, j, k)];
}
/////////////////////////////////////////////////////////////////////////////
// Dose methods...
// Fill DOSE/fluence matrix for secondary particles:
// If fluence parameter is true (default value is FALSE) then fluence at voxel (i, j, k) is
// increased. The energyDeposit parameter fill the dose matrix for voxel (i,j,k)
/////////////////////////////////////////////////////////////////////////////
G4bool HadrontherapyMatrix::Fill(G4int trackID, G4ParticleDefinition* particleDef, G4int i, G4int j,
                                 G4int k, G4double energyDeposit, G4bool fluence)
{
  /*G4AutoLock l(&matrixMutex);
  G4int tid = G4Threading::G4GetThreadId();*/

  // for server Segfault error
  // if i, j, k , 0 then trash
  if (i < 0 || i >= numberOfVoxelAlongX || j < 0 || j >= numberOfVoxelAlongY || k < 0
      || k >= numberOfVoxelAlongZ)
  {
    return false;
  }

  if ((energyDeposit <= 0. && !fluence) || !secondary) return false;

  // Get Particle Data Group particle ID
  G4int n = Index(i, j, k);
  G4int PDGencoding = particleDef->GetPDGEncoding();
  PDGencoding -= PDGencoding % 10;

  //
  if (energyDeposit > 0.)
  {
    localIonDose[PDGencoding][n] += energyDeposit;
  }
  if (fluence)
  {
    localIonFluence[PDGencoding][n]++;
  }

  //
  localParticleDef[PDGencoding] = particleDef;
  localIsPrimary[PDGencoding] = (trackID == 1);

  return true;
}

void HadrontherapyMatrix::StoreEventDataLineal(const std::map<G4int, G4double>& eventEnergyMap)
{
  {
    G4AutoLock l(&matrixMutex);
    totalEvents++;  // Event count

    // G4int tid = G4Threading::G4GetThreadId();
    //  Total Lineal Energy
    for (auto const& [voxelIndex, ep] : eventEnergyMap)
    {
      if (ep > 0.)
      {
        matrix[voxelIndex] += ep;
        matrixSquare[voxelIndex] += (ep * ep);
        matrixCube[voxelIndex] += (ep * ep * ep);
        matrixQuad[voxelIndex] += (ep * ep * ep * ep);
      }
    }
    // globalDoseBuffer.clear();

    // Secondary Ion
    std::set<G4int> activePDGs;
    for (const auto& pair : localIonDose)
      activePDGs.insert(pair.first);
    for (const auto& pair : localIonFluence)
      activePDGs.insert(pair.first);

    for (G4int pdg : activePDGs)
    {
      //
      ion_lineal* targetIon = nullptr;
      for (size_t idx = 0; idx < ionlinStore.size(); idx++)
      {
        if (ionlinStore[idx].PDGencoding == pdg)
        {
          targetIon = &ionlinStore[idx];
          break;
        }
      }

      //
      if (!targetIon)
      {
        G4ParticleDefinition* pDef = localParticleDef[pdg];
        G4String fullName = pDef->GetParticleName();
        G4String name = fullName.substr(0, fullName.find("["));

        ion_lineal newIon;
        newIon.isPrimary = localIsPrimary[pdg];
        newIon.PDGencoding = pdg;
        newIon.name = name;
        newIon.len = name.length();
        newIon.Z = pDef->GetAtomicNumber();
        newIon.A = pDef->GetAtomicMass();
        // [수정: new G4double[] 할당을 완전히 삭제하여 수백 MB 메모리 폭발 방지]

        ionlinStore.push_back(newIon);
        targetIon = &ionlinStore.back();
      }

      for (auto const& [voxelIndex, ep] : localIonDose[pdg])
      {
        unsigned int count = 0;
        if (localIonFluence[pdg].count(voxelIndex) > 0)
        {
          count = localIonFluence[pdg][voxelIndex];
        }
        if (ep > 0.) targetIon->sparseData[voxelIndex].fill(ep, count);
      }
    }
  }

  //
  localIonDose.clear();
  localIonFluence.clear();
  localParticleDef.clear();
  localIsPrimary.clear();

  localHitTrack.clear();
}

void HadrontherapyMatrix::StoreEventDataDose(const std::map<G4int, G4double>& eventEnergyMap)
{
  {
    G4AutoLock l(&matrixMutex);
    totalEvents++;

    // 1. Total Dose
    for (auto const& [voxelIndex, ep] : eventEnergyMap)
    {
      if (ep > 0.)
      {
        matrix[voxelIndex] += ep;
        matrixSumSq[voxelIndex] += (ep * ep);
      }
    }

    // 2. Secondary Ion
    std::set<G4int> activePDGs;
    for (const auto& pair : localIonDose)
      activePDGs.insert(pair.first);
    for (const auto& pair : localIonFluence)
      activePDGs.insert(pair.first);

    for (G4int pdg : activePDGs)
    {
      ion* targetIon = nullptr;
      for (size_t idx = 0; idx < ionStore.size(); idx++)
      {
        if (ionStore[idx].PDGencoding == pdg)
        {
          targetIon = &ionStore[idx];
          break;
        }
      }

      if (!targetIon)
      {
        G4ParticleDefinition* pDef = localParticleDef[pdg];
        G4String fullName = pDef->GetParticleName();
        G4String name = fullName.substr(0, fullName.find("["));

        ion newIon;
        newIon.isPrimary = localIsPrimary[pdg];
        newIon.PDGencoding = pdg;
        newIon.name = name;
        newIon.len = name.length();
        newIon.Z = pDef->GetAtomicNumber();
        newIon.A = pDef->GetAtomicMass();

        ionStore.push_back(newIon);
        targetIon = &ionStore.back();
      }

      // [수정: Buffer Flush - G4StatDouble의 fill 함수로 위임]
      for (auto const& [voxelIndex, ep] : localIonDose[pdg])
      {
        if (ep > 0.) targetIon->sparseStats[voxelIndex].fill(ep);
      }
      for (auto const& [voxelIndex, count] : localIonFluence[pdg])
      {
        targetIon->sparseFluence[voxelIndex] += count;
      }
    }
  }

  localIonDose.clear();
  localIonFluence.clear();
  localParticleDef.clear();
  localIsPrimary.clear();

  localHitTrack.clear();
}

void HadrontherapyMatrix::StoreDoseFluenceLinealBinary(G4String file)
{
  filename = file;
  std::sort(ionlinStore.begin(), ionlinStore.end());

  // [수정] 1. 포맷에 따른 파일 열기 모드 변경
  if (outputFormat == "csv")
  {
    ofs.open(filename + ".csv", std::ios::out);
  }
  else
  {
    ofs.open(filename + ".out", std::ios::out | std::ios::binary);
  }

  G4double r_um = HadrontherapyMatrix::sphereRadius / um;
  G4double meanChordLength_um = (4.0 / 3.0) * r_um;

  // // [디버깅 메세지 추가]
  // // 계산이 시작되기 직전에 터미널에 현재 사용 중인 반경 값을 출력합니다.
  // G4cout << "\n==================================================" << G4endl;
  // G4cout << "[DEBUG] Lineal Energy Calculation Output Started" << G4endl;
  // G4cout << "[DEBUG] Applied Sphere Radius : " << r_um << " um" << G4endl;
  // G4cout << "[DEBUG] Mean Chord Length   : " << meanChordLength_um << " um" << G4endl;
  // G4cout << "==================================================\n" << G4endl;

  if (ofs.is_open())
  {
    G4int numIons = ionlinStore.size();
    
    if (storeRawStats == false)
    {
      // [수정] 2. 헤더 작성 분기 처리
      if (outputFormat == "csv")
      {
        ofs << "i,j,k,total_yD,total_yD_SE,total_dose,total_dose_SE";
        for (size_t l = 0; l < ionlinStore.size(); l++)
        {
          G4String a = (ionlinStore[l].isPrimary) ? "_1" : "";
          G4String ionName = ionlinStore[l].name + a;
          ofs << "," << ionName << "_yD," << ionName << "_yD_SE," << ionName << "_dose," << ionName
              << "_dose_SE," << ionName << "_f";
        }
        ofs << "\n";
      }
      else
      {
        ofs.write(reinterpret_cast<const char*>(&numberOfVoxelAlongX), sizeof(G4int));
        ofs.write(reinterpret_cast<const char*>(&numberOfVoxelAlongY), sizeof(G4int));
        ofs.write(reinterpret_cast<const char*>(&numberOfVoxelAlongZ), sizeof(G4int));
        ofs.write(reinterpret_cast<const char*>(&numIons), sizeof(G4int));

        for (size_t l = 0; l < ionlinStore.size(); l++)
        {
          G4String a = (ionlinStore[l].isPrimary) ? "_1" : "";
          G4String ionName = ionlinStore[l].name + a;
          G4int nameLen = ionName.length();
          ofs.write(reinterpret_cast<const char*>(&nameLen), sizeof(G4int));
          ofs.write(ionName.c_str(), nameLen);
        }
      }

      G4double N = (G4double)totalEvents;
      G4double C = 1000.0 / meanChordLength_um;

      for (G4int i = 0; i < numberOfVoxelAlongX; i++)
        for (G4int j = 0; j < numberOfVoxelAlongY; j++)
          for (G4int k = 0; k < numberOfVoxelAlongZ; k++)
          {
            G4int n = Index(i, j, k);

            if (matrix[n] > 0)
            {
              // --- [1] Total yD 및 SE 계산 ---
              G4double total_yD = (matrixSquare[n] / matrix[n]) * C;
              G4double total_yD_SE = 0.;

              if (N > 1.0)
              {
                G4double sumY = matrixSquare[n];
                G4double sumX = matrix[n];
                G4double sumY2 = matrixQuad[n];
                G4double sumX2 = matrixSquare[n];
                G4double sumXY = matrixCube[n];

                G4double R = sumY / sumX;
                G4double s2_Y = (sumY2 - (sumY * sumY) / N) / (N - 1.0);
                G4double s2_X = (sumX2 - (sumX * sumX) / N) / (N - 1.0);
                G4double s_XY = (sumXY - (sumY * sumX) / N) / (N - 1.0);

                G4double varR = (N / (sumX * sumX)) * (s2_Y + R * R * s2_X - 2.0 * R * s_XY);
                if (varR > 0.) total_yD_SE = std::sqrt(varR) * C;
              }

              // --- [2] Total Dose 및 SE 계산 ---
              G4double total_dose = (matrix[n] / massOfVoxel / N) / doseUnit;
              G4double total_dose_SE = 0.;
              if (N > 1.0)
              {
                G4double sumX = matrix[n];
                G4double sumX2 = matrixSquare[n];
                G4double variance = (sumX2 - (sumX * sumX) / N) / (N - 1.0);
                if (variance > 0.)
                {
                  G4double errTotalEnergy = std::sqrt(variance / N);
                  total_dose_SE = (errTotalEnergy / massOfVoxel) / doseUnit;
                }
              }

              // [수정] 3. 좌표와 Total 데이터 쓰기 분기
              if (outputFormat == "csv")
              {
                ofs << i << "," << j << "," << k << "," << total_yD << "," << total_yD_SE << ","
                    << total_dose << "," << total_dose_SE;
              }
              else
              {
                ofs.write(reinterpret_cast<const char*>(&i), sizeof(G4int));
                ofs.write(reinterpret_cast<const char*>(&j), sizeof(G4int));
                ofs.write(reinterpret_cast<const char*>(&k), sizeof(G4int));
                ofs.write(reinterpret_cast<const char*>(&total_yD), sizeof(G4double));
                ofs.write(reinterpret_cast<const char*>(&total_yD_SE), sizeof(G4double));
                ofs.write(reinterpret_cast<const char*>(&total_dose), sizeof(G4double));
                ofs.write(reinterpret_cast<const char*>(&total_dose_SE), sizeof(G4double));
              }

              // 각 Ion 별 데이터 계산 및 쓰기
              for (size_t l = 0; l < ionlinStore.size(); l++)
              {
                G4double ion_yD = 0.;
                G4double ion_yD_SE = 0.;
                G4double ion_dose = 0.;
                G4double ion_dose_SE = 0.;
                unsigned int ion_f = 0;

                auto it = ionlinStore[l].sparseData.find(n);
                if (it != ionlinStore[l].sparseData.end())
                {
                  const LinealStat& stat = it->second;
                  ion_f = stat.f;

                  if (stat.e > 0)
                  {
                    // --- Ion yD 계산 ---
                    ion_yD = (stat.e2 / stat.e) * C;
                    if (N > 1.0)
                    {
                      G4double sumY = stat.e2;
                      G4double sumX = stat.e;
                      G4double sumY2 = stat.e4;
                      G4double sumX2 = stat.e2;
                      G4double sumXY = stat.e3;

                      G4double R = sumY / sumX;
                      G4double s2_Y = (sumY2 - (sumY * sumY) / N) / (N - 1.0);
                      G4double s2_X = (sumX2 - (sumX * sumX) / N) / (N - 1.0);
                      G4double s_XY = (sumXY - (sumY * sumX) / N) / (N - 1.0);

                      G4double varR = (N / (sumX * sumX)) * (s2_Y + R * R * s2_X - 2.0 * R * s_XY);
                      if (varR > 0.) ion_yD_SE = std::sqrt(varR) * C;
                    }

                    // --- Ion Dose 계산 ---
                    ion_dose = (stat.e / massOfVoxel / N) / doseUnit;
                    if (N > 1.0)
                    {
                      G4double sumX = stat.e;
                      G4double sumX2 = stat.e2;
                      G4double variance = (sumX2 - (sumX * sumX) / N) / (N - 1.0);
                      if (variance > 0.)
                      {
                        G4double err_ion = std::sqrt(variance / N);
                        ion_dose_SE = (err_ion / massOfVoxel) / doseUnit;
                      }
                    }
                  }
                }

                // [수정] 4. Ion 별 데이터 쓰기 분기
                if (outputFormat == "csv")
                {
                  ofs << "," << ion_yD << "," << ion_yD_SE << "," << ion_dose << "," << ion_dose_SE
                      << "," << ion_f;
                }
                else
                {
                  ofs.write(reinterpret_cast<const char*>(&ion_yD), sizeof(G4double));
                  ofs.write(reinterpret_cast<const char*>(&ion_yD_SE), sizeof(G4double));
                  ofs.write(reinterpret_cast<const char*>(&ion_dose), sizeof(G4double));
                  ofs.write(reinterpret_cast<const char*>(&ion_dose_SE), sizeof(G4double));
                  ofs.write(reinterpret_cast<const char*>(&ion_f), sizeof(unsigned int));
                }
              }
              if (outputFormat == "csv")
              {
                ofs << "\n";
              }
            }
          }
    }
    else
    {
      ofs << "," <<"_f";
    }
  }
  ofs.close();
}

void HadrontherapyMatrix::StoreDoseFluenceBinary(G4String file)
{
  filename = file;
  std::sort(ionStore.begin(), ionStore.end());

  // [수정] 1. 포맷에 따른 파일 열기 모드 변경
  if (outputFormat == "csv")
  {
    ofs.open(filename + ".csv", std::ios::out);
  }
  else
  {
    ofs.open(filename + ".bin", std::ios::out | std::ios::binary);
  }

  if (ofs.is_open())
  {
    G4int numIons = secondary ? ionStore.size() : 0;

    // [수정] 2. 헤더 작성 분기 처리
    if (outputFormat == "csv")
    {
      ofs << "i,j,k,Dose,Dose_SE";
      if (secondary)
      {
        for (size_t l = 0; l < ionStore.size(); l++)
        {
          G4String a = (ionStore[l].isPrimary) ? "_1" : "";
          G4String ionName = ionStore[l].name + a;
          ofs << "," << ionName << "_Dose," << ionName << "_Dose_SE," << ionName << "_f";
        }
      }
      ofs << "\n";
    }
    else
    {
      ofs.write(reinterpret_cast<const char*>(&numberOfVoxelAlongX), sizeof(G4int));
      ofs.write(reinterpret_cast<const char*>(&numberOfVoxelAlongY), sizeof(G4int));
      ofs.write(reinterpret_cast<const char*>(&numberOfVoxelAlongZ), sizeof(G4int));
      ofs.write(reinterpret_cast<const char*>(&numIons), sizeof(G4int));

      if (secondary)
      {
        for (size_t l = 0; l < ionStore.size(); l++)
        {
          G4String a = (ionStore[l].isPrimary) ? "_1" : "";
          G4String ionName = ionStore[l].name + a;
          G4int nameLen = ionName.length();
          ofs.write(reinterpret_cast<const char*>(&nameLen), sizeof(G4int));
          ofs.write(ionName.c_str(), nameLen);
        }
      }
    }

    G4double N = (G4double)totalEvents;

    for (G4int i = 0; i < numberOfVoxelAlongX; i++)
      for (G4int j = 0; j < numberOfVoxelAlongY; j++)
        for (G4int k = 0; k < numberOfVoxelAlongZ; k++)
        {
          G4int n = Index(i, j, k);

          if (matrix[n] > 0)
          {
            // Total Dose 및 SE 계산
            G4double dose = (matrix[n] / massOfVoxel / N) / doseUnit;
            G4double dose_SE = 0.;

            if (N > 1.0)
            {
              G4double sumX = matrix[n];
              G4double sumX2 = matrixSumSq[n];
              G4double variance = (sumX2 - (sumX * sumX) / N) / (N - 1.0);
              if (variance > 0.)
              {
                G4double errTotalEnergy = std::sqrt(variance / N);
                dose_SE = (errTotalEnergy / massOfVoxel) / doseUnit;
              }
            }

            // [수정] 3. 좌표와 Total Dose 데이터 쓰기 분기
            if (outputFormat == "csv")
            {
              ofs << i << "," << j << "," << k << "," << dose << "," << dose_SE;
            }
            else
            {
              ofs.write(reinterpret_cast<const char*>(&i), sizeof(G4int));
              ofs.write(reinterpret_cast<const char*>(&j), sizeof(G4int));
              ofs.write(reinterpret_cast<const char*>(&k), sizeof(G4int));
              ofs.write(reinterpret_cast<const char*>(&dose), sizeof(G4double));
              ofs.write(reinterpret_cast<const char*>(&dose_SE), sizeof(G4double));
            }

            if (secondary)
            {
              for (size_t l = 0; l < ionStore.size(); l++)
              {
                G4double ionDose = 0.;
                G4double ionDose_SE = 0.;
                unsigned int ion_f = 0;

                auto statIt = ionStore[l].sparseStats.find(n);
                if (statIt != ionStore[l].sparseStats.end())
                {
                  G4StatDouble& stat = statIt->second;

                  ionDose = (stat.sum_wx() / N) / massOfVoxel / doseUnit;

                  if (stat.n() > 1)
                  {
                    G4double err_ion = (stat.rms() * std::sqrt(stat.n())) / N;
                    ionDose_SE = (err_ion / massOfVoxel) / doseUnit;
                  }
                }

                auto fluIt = ionStore[l].sparseFluence.find(n);
                if (fluIt != ionStore[l].sparseFluence.end())
                {
                  ion_f = fluIt->second;
                }

                // [수정] 4. Ion 별 데이터 쓰기 분기
                if (outputFormat == "csv")
                {
                  ofs << "," << ionDose << "," << ionDose_SE << "," << ion_f;
                }
                else
                {
                  ofs.write(reinterpret_cast<const char*>(&ionDose), sizeof(G4double));
                  ofs.write(reinterpret_cast<const char*>(&ionDose_SE), sizeof(G4double));
                  ofs.write(reinterpret_cast<const char*>(&ion_f), sizeof(unsigned int));
                }
              }
            }
            if (outputFormat == "csv")
            {
              ofs << "\n";
            }
          }
        }
    ofs.close();
  }
}