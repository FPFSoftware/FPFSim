#include "DetectorConstruction.hh"
#include "DetectorConstructionMaterial.hh"
#include "DetectorConstructionMessenger.hh"
#include "LArBoxSD.hh"
#include "FASER2TrackerSD.hh"

#include "geometry/FASER2DetectorConstruction.hh"
#include "geometry/GeometricalParameters.hh"
#include "geometry/FASERnu2DetectorConstruction.hh"
#include "geometry/FORMOSADetectorConstruction.hh"
#include "geometry/FLArETPCDetectorConstruction.hh"
#include "geometry/FLArEHadCatcherMuonFinderConstruction.hh"
#include "geometry/BabyMINDDetectorConstruction.hh"

#include "fields/BabyMINDMagneticField.hh"

#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <G4VisAttributes.hh>
#include <G4Box.hh>
#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4SDManager.hh>
#include <G4RunManager.hh>
#include <G4GeometryManager.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4SolidStore.hh>
#include <G4UserLimits.hh>
#include <G4GlobalMagFieldMessenger.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include <G4GDMLParser.hh>

using namespace std;

G4ThreadLocal G4UniformMagField* DetectorConstruction::magField = 0;
G4ThreadLocal G4FieldManager* DetectorConstruction::fieldMgr = 0;
G4ThreadLocal G4UniformMagField* DetectorConstruction::magFieldFASER2 = 0;
G4ThreadLocal G4FieldManager* DetectorConstruction::fieldMgrFASER2 = 0;
G4ThreadLocal BabyMINDMagneticField* DetectorConstruction::babyMINDField = 0;
G4ThreadLocal G4FieldManager* DetectorConstruction::babyMINDFieldMgr = 0;

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(), 
    m_addFLArE(true), m_addFORMOSA(true), m_addFASERnu2(true), m_addFASER2(true), m_useBabyMIND(false)
{
  DefineMaterial();
  messenger = new DetectorConstructionMessenger(this);
  m_saveGdml = false;
  m_fileGdml = "FPF_geo.gdml";
  fCheckOverlap = false;
}

DetectorConstruction::~DetectorConstruction() 
{
  delete messenger;
}

void DetectorConstruction::DefineMaterial() {
  //-----------------------------
  // construction of materials
  //-----------------------------
  
  LArBoxMaterials = DetectorConstructionMaterial::GetInstance();
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  auto worldBox = new G4Box("worldBox", 30*m/2, 30*m/2, 150*m/2);
  worldLV = new G4LogicalVolume(worldBox, LArBoxMaterials->Material("Rock"), "worldLV");
  auto worldPV = new G4PVPlacement(nullptr,
                                   {},
                                   worldLV,
                                   "worldPV",
                                   nullptr,
                                   false, 
                                   0);

  // FPF long paper: https://dx.doi.org/10.1088/1361-6471/ac865e
  // section 2.1.1, figure 5, and figure 7
  G4double hallSizeX  = 9.4 * m;
  G4double hallSizeY  = 7.6 * m;
  G4double hallSizeZ  = 64.6 * m;
  
  // this offset accounts for: 
  // - distance between the entrance wall of the hall and the first detector, so the first detector 
  //   starts at the center of the global coordinate
  // - position of the cavern center w.r.t. the line of sight, since it's not in the exact middle
  G4ThreeVector hallOffset( GeometricalParameters::Get()->GetHallOffsetX(), 
                            GeometricalParameters::Get()->GetHallOffsetY(), 
                            hallSizeZ/2 - GeometricalParameters::Get()->GetHallHeadDistance()); 
                                                           
  auto hallBox = new G4Box("hallBox", hallSizeX/2, hallSizeY/2, hallSizeZ/2);
  hallLV = new G4LogicalVolume(hallBox, LArBoxMaterials->Material("Air"), "hallLV");
  auto hallPV = new G4PVPlacement(nullptr, hallOffset, hallLV, "hallPV", worldLV, false, 0, fCheckOverlap);

  //-----------------------------------
  // FLArE TPC volume

  if (m_addFLArE) {
    GeometricalParameters::Get()->SetAddFLArE(m_addFLArE);    
    FLArETPCDetectorConstruction *FLArETPCAssembler = new FLArETPCDetectorConstruction();
    G4double lArSizeZ               = GeometricalParameters::Get()->GetTPCSizeZ();
    G4double TPCInsulationThickness = GeometricalParameters::Get()->GetTPCInsulationThickness();

    G4LogicalVolume* FLArETPCAssembly = FLArETPCAssembler->GetFLArETPCAssembly();
    TPCModuleLogical = FLArETPCAssembler->GetFLArETPCVolume();

    // positioning
    G4double lengthFLArE = 2*TPCInsulationThickness + lArSizeZ; 
    G4ThreeVector FLArEPos = GeometricalParameters::Get()->GetFLArEPosition();
    FLArEPos -= hallOffset;
    new G4PVPlacement(nullptr, FLArEPos, FLArETPCAssembly, "FLArETPCPhysical", hallLV, false, 0, fCheckOverlap);

    G4cout << "Length of FLArE     : " << lengthFLArE << G4endl;
    G4cout << "Center of FLArE TPC : " << FLArEPos+hallOffset << G4endl; // w.r.t the global coordinate
  
    //-----------------------------------
    // FLArE HadCal/MuonCatcher or BabyMIND

    if( m_useBabyMIND ){   /// use BabyMIND
  
      BabyMINDDetectorConstruction *BabyMINDAssembler = new BabyMINDDetectorConstruction();
      G4LogicalVolume *BabyMINDAssembly = BabyMINDAssembler->GetBabyMINDAssembly();
  
      BabyMINDMagnetPlateLogical = BabyMINDAssembler->GetMagnetPlate();
      BabyMINDVerticalBar = BabyMINDAssembler->GetVerticalBar();
      BabyMINDHorizontalBar = BabyMINDAssembler->GetHorizontalBar();
    
      G4double babyMINDLengthZ  = GeometricalParameters::Get()->GetBabyMINDTotalSizeZ();
      G4ThreeVector babyMINDPos = GeometricalParameters::Get()->GetFLArEPosition() +
                                  G4ThreeVector(0.,0.,lArSizeZ/2.+TPCInsulationThickness) +
                                  G4ThreeVector(0.,0.,babyMINDLengthZ/2.);
      babyMINDPos -= hallOffset;
      new G4PVPlacement(nullptr, babyMINDPos, BabyMINDAssembly, "BabyMINDPhysical", hallLV, false, 0, fCheckOverlap);
      
      G4cout << "Length of BabyMIND : " << babyMINDLengthZ << G4endl;
      G4cout << "Center of BabyMIND : " << babyMINDPos+hallOffset << G4endl; // w.r.t the global coordinate
    
    }
    else{  //legacy HadCal/MuonCatcher
    
      FLArEHadCatcherMuonFinderConstruction *HadCatMuonFindAssembler = new FLArEHadCatcherMuonFinderConstruction();
      G4double HadCatcherLength       = GeometricalParameters::Get()->GetHadCalLength();
      G4double MuonFinderLength       = GeometricalParameters::Get()->GetMuonCatcherLength();
    
      G4LogicalVolume* HadCatMuonFindAssembly = HadCatMuonFindAssembler->GetHadCatcherMuonFinderAssembly();
      HadCalXCellLogical          = HadCatMuonFindAssembler->GetHadCalXVolume();
      HadCalYCellLogical          = HadCatMuonFindAssembler->GetHadCalYVolume();
      HadAbsorLayersLogical       = HadCatMuonFindAssembler->GetHadCalAbsorbVolume();
      MuonFinderXCellLogical      = HadCatMuonFindAssembler->GetMuonCatcherXVolume();
      MuonFinderYCellLogical      = HadCatMuonFindAssembler->GetMuonCatcherYVolume();
      MuonFinderAbsorLayersLogical= HadCatMuonFindAssembler->GetMuonCatcherAbsorbVolume();
      
      G4double HadCatMuonFindLengthZ  = HadCatcherLength + MuonFinderLength;
      G4ThreeVector HadCatMuonFindPos = GeometricalParameters::Get()->GetFLArEPosition() +
	                                      G4ThreeVector(0.,0.,lArSizeZ/2.+TPCInsulationThickness) +
				                                G4ThreeVector(0.,0.,HadCatMuonFindLengthZ/2.);
    
      HadCatMuonFindPos -= hallOffset;
      new G4PVPlacement(nullptr, HadCatMuonFindPos, HadCatMuonFindAssembly, "HadCatMuonFindPhysical", hallLV, false, 0, fCheckOverlap);
      
      G4cout << "Length of HadCatherMuonFinder : " << HadCatMuonFindLengthZ << G4endl;
      G4cout << "Center of HadCatherMuonFinder : " << HadCatMuonFindPos+hallOffset << G4endl; // w.r.t the global coordinate
    }
  }
  
  //-----------------------------------
  // FORMOSA

  if (m_addFORMOSA) {
    GeometricalParameters::Get()->SetAddFORMOSA(m_addFORMOSA);
    FORMOSADetectorConstruction *FORMOSAAssembler = new FORMOSADetectorConstruction();
    G4LogicalVolume* FORMOSAAssembly = FORMOSAAssembler->GetFORMOSAAssembly();
    FORMOSAScintillatorBarLogical = FORMOSAAssembler->GetScintillatorBar();

    // positioning
    G4double lengthFORMOSA = GeometricalParameters::Get()->GetFORMOSATotalSizeZ();
    G4ThreeVector FORMOSAPos = GeometricalParameters::Get()->GetFORMOSAPosition();
    FORMOSAPos -= hallOffset;
    new G4PVPlacement(nullptr, FORMOSAPos, FORMOSAAssembly, "FORMOSAPhysical", hallLV, false, 0, fCheckOverlap);

    G4cout<<"Length of FORMOSA : "<<lengthFORMOSA<<G4endl;
    G4cout<<"Center of FORMOSA : "<<FORMOSAPos+hallOffset<<G4endl; // w.r.t the global coordinate
  }
                         
  //-----------------------------------
  // FASERnu2 Emulsion Detector

  if (m_addFASERnu2) {
    GeometricalParameters::Get()->SetAddFASERnu2(m_addFASERnu2);
    FASERnu2DetectorConstruction *FASERnu2Assembler = new FASERnu2DetectorConstruction();
    FASERnu2EmulsionLogical = FASERnu2Assembler->GetEmulsionFilm();
    FASERnu2VetoInterfaceLogical = FASERnu2Assembler->GetVetoInterfaceDetector();
    G4LogicalVolume* FASERnu2Assembly = FASERnu2Assembler->GetFASERnu2Assembly();
    
    // positioning
    G4double lengthFASERnu2 = GeometricalParameters::Get()->GetFASERnu2TotalSizeZ();
    G4ThreeVector FASERnu2Pos = GeometricalParameters::Get()->GetFASERnu2Position();
    FASERnu2Pos -= hallOffset;
    new G4PVPlacement(nullptr, FASERnu2Pos, FASERnu2Assembly, "FASERnu2Physical", hallLV, false, 0, fCheckOverlap);

    G4cout<<"Length of FASERnu2 : "<<lengthFASERnu2<<G4endl;
    G4cout<<"Center of FASERnu2 : "<<FASERnu2Pos+hallOffset<<G4endl; // w.r.t the global coordinate
  }

  //-----------------------------------
  // FASER2 Magnet + Tracking stations

  if (m_addFASER2) {
    GeometricalParameters::Get()->SetAddFASER2(m_addFASER2);
    FASER2DetectorConstruction *magnetAssembler = new FASER2DetectorConstruction();
    FASER2MagnetLogical = magnetAssembler->GetMagneticVolume(); //need to assign B field
    FASER2TrackingLogical = magnetAssembler->GetTrackingStations();
    FASER2VetoLogical = magnetAssembler->GetVetoStations();
    FASER2MuonLogical = magnetAssembler->GetMuonStations();
    G4LogicalVolume* FASER2Assembly = magnetAssembler->GetFASER2Assembly();

    // positioning
    G4double lengthFASER2Assembly = GeometricalParameters::Get()->GetFASER2TotalSizeZ();
    G4ThreeVector FASER2Pos = GeometricalParameters::Get()->GetFASER2Position();
    GeometricalParameters::Get()->SetMagnetZPosition(FASER2Pos.z()); // save for momentum analysis
    FASER2Pos -= hallOffset;
    new G4PVPlacement(nullptr, FASER2Pos, FASER2Assembly, "FASER2Physical", hallLV, false, 0, fCheckOverlap);

    //detectorGapLength = 1.2*m;
    G4cout<<"Length of FASER2 Spectrometer : "<<lengthFASER2Assembly<<G4endl;
    G4cout<<"Center of FASER2 Spectrometer : "<<FASER2Pos+hallOffset<<G4endl; // w.r.t the global coordinate
  }
  
  //-------------------------------------------------------------------

  // visualization
  G4VisAttributes* worldVis = new G4VisAttributes(G4Colour(167./255, 168./255, 189./255));
  worldVis->SetVisibility(false);
  worldLV->SetVisAttributes(worldVis);

  G4VisAttributes* hallVis = new G4VisAttributes(G4Colour(167./255, 168./255, 189./255));
  hallVis->SetVisibility(true);
  hallLV->SetVisAttributes(hallVis);

  if (m_saveGdml) {
    G4GDMLParser fParser;
    G4cout << "Exporting geometry to " << m_fileGdml << G4endl;
    fParser.Write(m_fileGdml, worldPV, false);
  }

  return worldPV;
}

void DetectorConstruction::ConstructSDandField() {
  
  G4SDManager* sdManager = G4SDManager::GetSDMpointer();
  sdManager->SetVerboseLevel(2);
  int SDIdx = 0;

  if (m_addFLArE) {

    LArBoxSD* TPCModuleSD = new LArBoxSD("lArBoxSD");
    TPCModuleLogical->SetSensitiveDetector(TPCModuleSD);
    sdManager->AddNewDetector(TPCModuleSD);
    GeometricalParameters::Get()->AddSD2List(SDIdx, "lArBoxSD/lar_box");
    SDIdx++;
    
    if (m_useBabyMIND) {
    
      LArBoxSD* BabyMINDHorBarSD = new LArBoxSD("BabyMINDHorBarSD");
      BabyMINDHorizontalBar->SetSensitiveDetector(BabyMINDHorBarSD);
      sdManager->AddNewDetector(BabyMINDHorBarSD);
      GeometricalParameters::Get()->AddSD2List(SDIdx, "BabyMINDHorBarSD/lar_box");
      SDIdx++;
      
      LArBoxSD* BabyMINDVerBarSD = new LArBoxSD("BabyMINDVerBarSD");
      BabyMINDVerticalBar->SetSensitiveDetector(BabyMINDVerBarSD);
      sdManager->AddNewDetector(BabyMINDVerBarSD);
      GeometricalParameters::Get()->AddSD2List(SDIdx, "BabyMINDVerBarSD/lar_box");
      SDIdx++;

      // magnetic field for BabyMIND
      babyMINDField = new BabyMINDMagneticField();
      babyMINDFieldMgr = new G4FieldManager();
      babyMINDFieldMgr->SetDetectorField(babyMINDField);
      babyMINDFieldMgr->CreateChordFinder(babyMINDField);
      BabyMINDMagnetPlateLogical->SetFieldManager(babyMINDFieldMgr, true);
    
    } else {
    
      LArBoxSD* HadCalXSD = new LArBoxSD("HadCalXSD");
      HadCalXCellLogical->SetSensitiveDetector(HadCalXSD);
      sdManager->AddNewDetector(HadCalXSD);
      GeometricalParameters::Get()->AddSD2List(SDIdx, "HadCalXSD/lar_box");
      SDIdx++;

      LArBoxSD* HadCalYSD = new LArBoxSD("HadCalYSD");
      HadCalYCellLogical->SetSensitiveDetector(HadCalYSD);
      sdManager->AddNewDetector(HadCalYSD);
      GeometricalParameters::Get()->AddSD2List(SDIdx, "HadCalYSD/lar_box");
      SDIdx++;

      LArBoxSD* MuonFinderXSD = new LArBoxSD("MuonFinderXSD");
      MuonFinderXCellLogical->SetSensitiveDetector(MuonFinderXSD);
      sdManager->AddNewDetector(MuonFinderXSD);
      GeometricalParameters::Get()->AddSD2List(SDIdx, "MuonFinderXSD/lar_box");
      SDIdx++;

      LArBoxSD* MuonFinderYSD = new LArBoxSD("MuonFinderYSD");
      MuonFinderYCellLogical->SetSensitiveDetector(MuonFinderYSD);
      sdManager->AddNewDetector(MuonFinderYSD);
      GeometricalParameters::Get()->AddSD2List(SDIdx, "MuonFinderYSD/lar_box");
      SDIdx++;

      LArBoxSD* HadAbsorbSD = new LArBoxSD("HadAbsorbSD");
      HadAbsorLayersLogical->SetSensitiveDetector(HadAbsorbSD);
      sdManager->AddNewDetector(HadAbsorbSD);
      GeometricalParameters::Get()->AddSD2List(SDIdx, "HadAbsorbSD/lar_box");
      SDIdx++;

      LArBoxSD* MuonFinderAbsorbSD = new LArBoxSD("MuonFinderAbsorbSD");
      MuonFinderAbsorLayersLogical->SetSensitiveDetector(MuonFinderAbsorbSD);
      sdManager->AddNewDetector(MuonFinderAbsorbSD);
      GeometricalParameters::Get()->AddSD2List(SDIdx, "MuonFinderAbsorbSD/lar_box");
      SDIdx++;
    
      // magnetic field for HadCatcher + MuonFinder
      G4ThreeVector fieldValue = G4ThreeVector(0,fFieldValue, 0);
      magField = new G4UniformMagField(fieldValue);
      fieldMgr = new G4FieldManager();
      fieldMgr->SetDetectorField(magField);
      fieldMgr->CreateChordFinder(magField);
      HadCalXCellLogical->SetFieldManager(fieldMgr, true);
      HadCalYCellLogical->SetFieldManager(fieldMgr, true);
      HadAbsorLayersLogical->SetFieldManager(fieldMgr, true);
      MuonFinderXCellLogical->SetFieldManager(fieldMgr, true);
      MuonFinderYCellLogical->SetFieldManager(fieldMgr, true);
      MuonFinderAbsorLayersLogical->SetFieldManager(fieldMgr, true);
    }
  }

  if (m_addFORMOSA) {
    LArBoxSD* ScintillatorBarSD = new LArBoxSD("FORMOSAScinBarSD");
    FORMOSAScintillatorBarLogical->SetSensitiveDetector(ScintillatorBarSD);
    sdManager->AddNewDetector(ScintillatorBarSD);
    GeometricalParameters::Get()->AddSD2List(SDIdx, "FORMOSAScinBarSD/lar_box");
    SDIdx++;
  }

  if (m_addFASERnu2) {
    LArBoxSD* EmulsionFilmSD = new LArBoxSD("FASERnu2EmulsionSD");
    FASERnu2EmulsionLogical->SetSensitiveDetector(EmulsionFilmSD);
    sdManager->AddNewDetector(EmulsionFilmSD);
    GeometricalParameters::Get()->AddSD2List(SDIdx, "FASERnu2EmulsionSD/lar_box");
    SDIdx++;
    
    LArBoxSD* VetoInterfaceSD = new LArBoxSD("FASERnu2VetoInterfaceSD");
    FASERnu2VetoInterfaceLogical->SetSensitiveDetector(VetoInterfaceSD);
    sdManager->AddNewDetector(VetoInterfaceSD);
    GeometricalParameters::Get()->AddSD2List(SDIdx, "FASERnu2VetoInterfaceSD/lar_box");
    SDIdx++;
  }

  if (m_addFASER2) {
    
    FASER2TrackerSD* TrkScinSD = new FASER2TrackerSD("FASER2TrackerSD");
    FASER2TrackingLogical->SetSensitiveDetector(TrkScinSD);
    sdManager->AddNewDetector(TrkScinSD);
    GeometricalParameters::Get()->AddSD2List(SDIdx, "FASER2TrackerSD/tracker_box");
    SDIdx++;

    // LArBoxSD* FASER2VetoSD = new LArBoxSD("FASER2VetoSD"); //TODO: Should implement a new SD class for FASER2Veto
    // FASER2VetoLogical->SetSensitiveDetector(FASER2VetoSD);
    // sdManager->AddNewDetector(FASER2VetoSD);
    // GeometricalParameters::Get()->AddSD2List(SDIdx, "FASER2VetoSD/lar_box");
    // SDIdx++;

    // LArBoxSD* FASER2MuonDetSD = new LArBoxSD("FASER2MuonSD"); //TODO: Should probably be a FASER2TrackerSD
    // FASER2MuonLogical->SetSensitiveDetector(FASER2MuonDetSD);
    // sdManager->AddNewDetector(FASER2MuonDetSD);
    // GeometricalParameters::Get()->AddSD2List(SDIdx, "FASER2MuonSD/lar_box");
    // SDIdx++;

    // FASER2 magnetic field
    G4ThreeVector fieldValueFASER2 = GeometricalParameters::Get()->GetFASER2MagnetField();
    magFieldFASER2 = new G4UniformMagField(fieldValueFASER2);
    fieldMgrFASER2 = new G4FieldManager();
    fieldMgrFASER2->SetDetectorField(magFieldFASER2);
    fieldMgrFASER2->CreateChordFinder(magFieldFASER2);
    FASER2MagnetLogical->SetFieldManager(fieldMgrFASER2, true);
  }
}

//void DetectorConstruction::UpdateGeometry() {
  //// clean up previous geometry
  //G4GeometryManager::GetInstance()->OpenGeometry();

  //G4PhysicalVolumeStore::GetInstance()->Clean();
  //G4LogicalVolumeStore::GetInstance()->Clean();
  //G4SolidStore::GetInstance()->Clean();
  //G4LogicalSkinSurface::CleanSurfaceTable();
  //G4LogicalBorderSurface::CleanSurfaceTable();
  //G4SurfaceProperty::CleanSurfacePropertyTable();

  //// define new one
  //G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
//  G4RunManager::GetRunManager()->GeometryHasBeenModified();
//}
