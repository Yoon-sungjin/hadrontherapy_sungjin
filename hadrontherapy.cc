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
// ----------------------------------------------------------------------------
//                         GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
//
//                                     MAIN AUTHORS
//                                  ====================
//         G.A.P. Cirrone(a)*, L. Pandola(a), G. Petringa(a), S.Fattori(a), A.Sciuto(a)
//              *Corresponding author, email to pablo.cirrone@lns.infn.it
//
//
//                 ==========>   PAST CONTRIBUTORS  <==========
//
//                    R.Calcagno(a), G.Danielsen (b), F.Di Rosa(a),
//                    S.Guatelli(c), A.Heikkinen(b), P.Kaitaniemi(b),
//                    A.Lechner(d), S.E.Mazzaglia(a), Z. Mei(h), G.Milluzzo(a),
//                    M.G.Pia(e), F.Romano(a), G.Russo(a,g),
//                    M.Russo(a), A.Tramontana (a), A.Varisano(a)
//
//              (a) Laboratori Nazionali del Sud of INFN, Catania, Italy
//              (b) Helsinki Institute of Physics, Helsinki, Finland
//              (c) University of Wallongong, Australia
//              (d) CERN, Geneve, Switzwerland
//              (e) INFN Section of Genova, Genova, Italy
//              (f) Physics and Astronomy Department, Univ. of Catania, Catania, Italy
//              (g) CNR-IBFM, Italy
//              (h) Institute of Applied Electromagnetic Engineering(IAEE)
//                  Huazhong University of Science and Technology(HUST), Wuhan, China
//
//
//                                          WEB
//                                      ===========
//       https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy
//
// ----------------------------------------------------------------------------

#include "G4ParallelWorldPhysics.hh"
#include "G4PhysListFactory.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4ScoringManager.hh"
#include "G4Timer.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "HadrontherapyActionInitialization.hh"
#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyEventAction.hh"
#include "HadrontherapyGeometryController.hh"
#include "HadrontherapyGeometryMessenger.hh"
#include "HadrontherapyInteractionParameters.hh"
#include "HadrontherapyLet.hh"
#include "HadrontherapyMatrix.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyRunAction.hh"
#include "HadrontherapySteppingAction.hh"
// #include "QBBC.hh"
#include <time.h>

//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4String fileName = "";
  G4String filePrefix = "";  // 파일 이름 앞에 붙을 접두사

  if (argc >= 2) fileName = argv[1];  // 첫 번째 인자는 매크로 파일
  if (argc >= 3) filePrefix = argv[2];  // 두 번

  G4UIExecutive* ui = nullptr;


  // Instantiate the G4Timer object, to monitor the CPU time spent for
  // the entire execution
  G4Timer* theTimer = new G4Timer();
  // Start the benchmark
  theTimer->Start();

  // Set the Random engine
  // The following guarantees random generation also for different runs
  // in multithread
  CLHEP::RanluxEngine defaultEngine(1234567, 4);
  G4Random::setTheEngine(&defaultEngine);
  G4int seed = (G4int)time(NULL);
  G4Random::setTheSeed(seed);

  // Construct the default run manager
  //
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // Define the number of threads for the simulation runs
  G4int nThreads = G4Threading::G4GetNumberOfCores();
  runManager->SetNumberOfThreads(nThreads-1);

  // Geometry controller is responsible for instantiating the geometries.
  //
  HadrontherapyGeometryController* geometryController = new HadrontherapyGeometryController();

  // Connect the geometry controller to the G4 user interface
  //
  HadrontherapyGeometryMessenger* geometryMessenger =
    new HadrontherapyGeometryMessenger(geometryController);

  G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManager();
  scoringManager->SetVerboseLevel(1);

  // Initialize the default Hadrontherapy geometry
  geometryController->SetGeometry("default");

  // Initialize the physics
  G4PhysListFactory factory;
  G4VModularPhysicsList* physicsList = 0;
  physicsList = new HadrontherapyPhysicsList();
  //   G4VModularPhysicsList* physicsList = new QBBC();
  //   runManager->SetUserInitialization(physicsList);

  if (physicsList)
  {
    G4cout << "Going to register G4ParallelWorldPhysics" << G4endl;
    physicsList->RegisterPhysics(new G4ParallelWorldPhysics("DetectorROGeometry"));
  }

  // Initialisations of physics
  runManager->SetUserInitialization(physicsList);

  // Initialisation of the Actions
  runManager->SetUserInitialization(new HadrontherapyActionInitialization);

  // Initialize command based scoring
  G4ScoringManager::GetScoringManager();

  // Interaction data: stopping powers
  //
  HadrontherapyInteractionParameters* pInteraction = new HadrontherapyInteractionParameters(true);

  // Initialize analysis
  //
  HadrontherapyAnalysis* analysis = HadrontherapyAnalysis::GetInstance();

  // Initialise the Visualisation
  //
  auto visManager = new G4VisExecutive(argc, argv);
  visManager->Initialize();

  //** Get the pointer to the User Interface manager
  //
  auto UImanager = G4UImanager::GetUIpointer();

  if (!(argc == 1))
  {
    // batch mode
    G4String command = "/control/execute ../macro/";
    UImanager->ApplyCommand(command + fileName);

    // 2. 핵심: 매크로 실행 후, 시각화 창(Viewer)이 켜졌는지 검사합니다.
    bool manualOpenUI = false;
    if (manualOpenUI)
    {
      G4cout << "=== Visualization viewer detected! Entering interactive session... ===" << G4endl;
      
      // 매크로에서 /vis/open을 통해 창을 켰다면, 
      // 프로그램이 종료되지 않도록 늦게라도 UI 세션을 시작해 대기 상태로 만듭니다.
      ui = new G4UIExecutive(argc, argv);
      ui->SessionStart();
      // 사용자가 창을 닫거나 exit를 치면 비로소 아래로 넘어갑니다.
      
    }
  }

  else
  {
    // UI mode
    UImanager->ApplyCommand("/control/execute ../macro/carbon_beamline.mac");
    ui = new G4UIExecutive(argc, argv);
    ui->SessionStart();
    // delete ui;
  }

  

  // Stop the benchmark here
  theTimer->Stop();

  G4cout << "The simulation took: " << theTimer->GetRealElapsed() << " s to run (real time)"
         << G4endl;

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  const G4UserEventAction* userEventAction = runManager->GetUserEventAction();

  // 2. 우리가 만든 HadrontherapyEventAction 타입으로 형변환(Casting) 합니다.
  // auto* evAction = (HadrontherapyEventAction*)userEventAction;

  if (HadrontherapyMatrix* pMatrix = HadrontherapyMatrix::GetInstance())
  {
    // evAction 대신 pMatrix에서 모드를 확인합니다.
    if (pMatrix->GetCalculationMode() == "Lineal")
    {
      G4cout << "Saving " << filePrefix << "LinealEnergy.out" << G4endl;
      pMatrix->StoreDoseFluenceLinealBinary(filePrefix + "LinealEnergy");
    }
    else
    {
      G4cout << "Saving " << filePrefix << "Dose.out" << G4endl;
      pMatrix->StoreDoseFluenceBinary(filePrefix + "Dose");

      HadrontherapyLet* let = HadrontherapyLet::GetInstance();
      if (let && let->CalculationDone)
      {
        G4cout << "Saving " << filePrefix << "LET.out" << G4endl;
        let->LetOutput();
        let->StoreLetBinary(filePrefix + "LET");
      }
    }
  }

  if (ui){
    delete ui;
  }
  delete visManager;

  delete geometryMessenger;
  delete geometryController;
  delete pInteraction;
  delete runManager;
  delete analysis;
  return 0;
}
