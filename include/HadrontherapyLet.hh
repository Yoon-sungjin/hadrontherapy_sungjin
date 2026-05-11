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

#ifndef HadrontherapyLet_h
#  define HadrontherapyLet_h 1
#endif

#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include "HadrontherapyMatrix.hh"

#include <fstream>
#include <string>
#include <vector>
struct LetStat
{
    G4double dN{0.}, dD{0.}, tN{0.}, tD{0.};  // 누적 (Numerator, Denominator)
    G4double dN2{0.}, dD2{0.}, tN2{0.}, tD2{0.};  // 분산 계산용 제곱 합
    G4double dCov{0.}, tCov{0.};  // 공분산

    // 계산된 최종 LET 값과 오차
    G4double letD{0.}, letD_SE{0.};
    G4double letT{0.}, letT_SE{0.};

    void fill(G4double nDose, G4double dDose, G4double nTrack, G4double dTrack)
    {
      dN += nDose;
      dD += dDose;
      tN += nTrack;
      tD += dTrack;

      dN2 += nDose * nDose;
      dD2 += dDose * dDose;
      tN2 += nTrack * nTrack;
      tD2 += dTrack * dTrack;

      dCov += nDose * dDose;
      tCov += nTrack * dTrack;
    }
};

struct ionLet
{
    G4bool isPrimary;
    G4int PDGencoding;
    G4String fullName;
    G4String name;
    G4int Z;
    G4int A;

    // [수정] 수십 개의 거대한 포인터 배열을 하나의 Sparse Map으로 압축
    std::unordered_map<G4int, LetStat> sparseLet;

    G4bool operator<(const ionLet& a) const
    {
      return (this->Z == a.Z) ? this->A < a.A : this->Z < a.Z;
    }
};

class G4Material;
class HadrontherapyMatrix;
class HadrontherapyPrimaryGeneratorAction;
class HadrontherapyInteractionParameters;
class HadrontherapyDetectorConstruction;

class HadrontherapyLet
{
  private:
    HadrontherapyLet(HadrontherapyDetectorConstruction*);

  public:
    ~HadrontherapyLet();
    static HadrontherapyLet* GetInstance(HadrontherapyDetectorConstruction*);
    static HadrontherapyLet* GetInstance();
    static G4bool CalculationDone;
    void Initialize();
    void Clear();

    void Fill(G4int i, G4int j, G4int k, G4double DE, G4double DX);
    void FillEnergySpectrum(G4int trackID, G4ParticleDefinition* particleDef, G4double ekinMean,
                            const G4Material* mat, G4double DE, G4double DEEletrons, G4double DX,
                            G4int i, G4int j, G4int k);
    void LetOutput();
    void StoreLetBinary(G4String filename = "LET.out");
    void StoreLETEventData();

  private:
    static HadrontherapyLet* instance;
    HadrontherapyPrimaryGeneratorAction* pPGA;
    G4Material* detectorMat;
    G4double density;
    G4String filename;
    std::ofstream ofs;
    HadrontherapyMatrix* matrix;

    G4int nVoxels, numberOfVoxelAlongX, numberOfVoxelAlongY, numberOfVoxelAlongZ;

    // [수정] 글로벌 배열들을 하나의 Sparse Map으로 대체
    std::unordered_map<G4int, LetStat> totalLetStats;

    std::vector<ionLet> ionLetStore;
};