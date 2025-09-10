////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeReco
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyzeReco_module.cc
//
// Generated at Tue Sep  2 03:17:05 2025 by Carlos Martin Morales using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

// All the includes that I will need

// Basic ones:
#include "art/Framework/Core/FileBlock.h"             
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Additional framework includes
#include "art_root_io/TFileService.h"

// ROOT includes
#include <TTree.h>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"
#include <larcorealg/Geometry/Exceptions.h> 
#include <vector>
#include <limits>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>                                                                                                                                                                       
#include <string>

// raw type variables include
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"      // ??

// POT includes
#include "larcoreobj/SummaryData/POTSummary.h"

// Generator includes (MCT)
#include "nusimdata/SimulationBase/MCTruth.h"       // MCTruth
#include "nusimdata/SimulationBase/MCNeutrino.h"    // MCNeutrino
#include "lardataobj/MCBase/MCTrack.h"              // MAYBE I need it for truth matching

// g4 includes (MCP)
#include "nusimdata/SimulationBase/MCParticle.h"    // MCParticle

// PFP includes
#include "lardataobj/RecoBase/PFParticle.h"              // PFP
#include "lardataobj/RecoBase/PFParticleMetadata.h"      // nu_score, is_clear_cosmic, ...
#include "lardataobj/RecoBase/Slice.h"                   // Slice
#include "lardataobj/RecoBase/Track.h"                   // Track
#include "lardataobj/RecoBase/Shower.h"                  // Shower
#include "lardataobj/RecoBase/Vertex.h"                  // Vertex
#include "lardataobj/RecoBase/Hit.h"                     // Hit
#include "lardataobj/AnalysisBase/ParticleID.h"          // I am not sure why it is included :(

// Calorimetry includes
#include "lardataobj/AnalysisBase/Calorimetry.h"    // Calorimetry
#include "lardataobj/RecoBase/SpacePoint.h"         // Calorimetry with hits

// Service includes (clock, detector properties,...)
#include "lardata/DetectorInfoServices/DetectorClocksService.h"         // Time
#include "lardataobj/AnalysisBase/T0.h"                                 // Time (t_0)
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"             // Detector services
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Geometry includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

// Associations and pointers includes
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// TruthMatching includes
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"

// I define the minimum possible int
constexpr int def_int = std::numeric_limits<int>::min();
// I define the maximum possible float
constexpr float def_float   = -std::numeric_limits<float>::max();


namespace sbnd {
  class AnalyzeReco;
}


class sbnd::AnalyzeReco : public art::EDAnalyzer {
public:
  explicit AnalyzeReco(fhicl::ParameterSet const& p);            // This is the constructor
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeReco(AnalyzeReco const&) = delete;                 // Deletes the method to copy this class to another variable
  AnalyzeReco(AnalyzeReco&&) = delete;                      // Deletes the method to move this class
  AnalyzeReco& operator=(AnalyzeReco const&) = delete;      // Deletes the copy by assigment operator
  AnalyzeReco& operator=(AnalyzeReco&&) = delete;           // Deletes the move by assigment operator

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.

  // Begin/end functions
  void beginJob() override;
  void endJob() override;
  void beginSubRun(const art::SubRun &sr);
  void endSubRun(const art::SubRun &sr);

  // Methods used for the analysis

  // Files method
  void respondToOpenInputFile(art::FileBlock const&); 

  // Method to get the total events
  int getTotalGenEvents(const art::Event &e);

  //Reestart variables
  void clearMaps();
  void setupMaps(const art::Event &e, const art::ValidHandle<std::vector<recob::Hit>> &hit_handle,
                                       const art::ValidHandle<std::vector<recob::PFParticle>> &PFP_Handle);

  void clearVars();
  void clearSubRunVars();
  void setDefaultGenVars();
  void setDefaultG4Vars();
  void setDefaultRecoVars();

  //Inside analyse functions
  // Analysis of truth (generator)
  void analyseTruth_gen(const art::Ptr<simb::MCTruth> MCT);

  void analyseTruth_g4(const std::vector<art::Ptr<simb::MCParticle>> MCP_vec);

  // Analysis of PFPs (reco, pandora)
  void analysePFPs(const art::Event &e,
          const std::vector<art::Ptr<simb::MCTruth>> MCT_vec,
          const std::vector<art::Ptr < simb::MCParticle>> MCP_vec,
          const art::Ptr<recob::PFParticle> &prim,                                       
          const std::vector<art::Ptr<recob::PFParticle>> &PFP_vec,
					const art::ValidHandle<std::vector<recob::PFParticle>> &PFP_handle, 
					const art::ValidHandle<std::vector<recob::Track>> &track_handle,
          const art::ValidHandle<std::vector<recob::Shower>> &shower_handle
					);

  // Analysis of tracks (reco, pandora)
  void analyseTrack(const art::Event &e, const art::Ptr<recob::Track> &track,
					  const art::ValidHandle<std::vector<recob::Track>> &track_handle);

  // Analysis of showers (reco, pandora)
  void analyseShower(const art::Event &e, const art::Ptr<recob::Shower> &shower,
					  const art::ValidHandle<std::vector<recob::Shower>> &shower_handle);

  // Analysis of calorimetry (reco, pandora)
  void getCalo(const art::Ptr<anab::Calorimetry> &calo);

  // For truth matching
  void getTrackTruthMatch(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> hits);

  void getSliceTruthMatch(const art::Event &e);

  art::Ptr<recob::PFParticle> getPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &PFP_vec);

  // Methods to calculate the completeness and the purity
  float completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);
  float purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);

  // Method to get the chi^2
  void getChi2PID(const art::Ptr<anab::ParticleID> &chi2pid);


private:
  // Vraibles that I need for truthmatching
  art::ServiceHandle<cheat::BackTrackerService> backTracker;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;

  // Declare member data here.

  // I will include the variables of my analysis and the TTrees here
  // Tree for the subruns (save POT of each subrun)
  TTree *fSubRunTree;

  //SubRun TTree variables
  double fPOT;                // POT of each subrun
  double fSpill;              // Number of neutrino spills of each subrun
  unsigned int fSubrun_ID;    // Subrun ID
  int fNum_gen_evts;          // Number of events generated in each subrun

  // output TTree
  TTree* fTree;

  // output TTree variables
  // General variables
  std::string fFilename;      // Filename that has been analyzed
  unsigned int fEvent_ID;     // Event ID
  unsigned int fRun_ID;       // Run ID

  //MCTruth variables
  //Consult https://internal.dunescience.org/doxygen/classsimb_1_1MCNeutrino.html#a8e213917487f1f5cd52ca5280d33a777 for information on the numbers meaning
  int fNu_PDG;                          // PDGcode of the neutrino
  double fNu_E0;                        // Energy of the neutrino
  double fNu_weight;                    // Weight of the interaction (1 always for GENIE; =!1 for GiBUU)
  int fGen_index;                       // Generator index, it should be 1 for genie, 2 for corsika and 0 for unknown
  int fNu_interaction_mode;    // interaction mode of the neutrino
  int fNu_interaction_type;    // Interaction type of the neutrino (expansion of interaction mode)
  int fNu_CC_NC;                        // CC(0) or NC(1) interaction
  int fNu_target;                       // PDG of the nuclear target
  int fNu_HitNuc;                       // PDG of the interacting nucleons (proton or neutron)
  int fNu_HitQuark;                     // PDG of the hit quark (for DIS interactions)
  double fNu_W;                         // Hadronic invariant mass (in GeV)
  double fNu_X;                         // Bjorken variable
  double fNu_Y;                         // Inelasticity
  double fNu_Q2;                        // Momentumm transfer Q^2 in GeV^2

  //MCTruth Particle Variables
  std::vector<int> fGen_part_trackID;       // ID of the track of this particle
  std::vector<int> fGen_part_statusCode;    // statusCode of the particle. Stable final state particles have statusCode=1
  std::vector<int> fGen_part_mother;        // ID of the mother of the particle
  std::vector<int> fGen_part_PDGcode;       // PDGcode of the particle
  std::vector<double> fGen_part_mass;       // Mass of the particle
  std::vector<double> fGen_part_E0;         // Initial energy of the particle
  std::vector<double> fGen_part_StartPos_x, fGen_part_StartPos_y, fGen_part_StartPos_z;   // Start position of the particle
  std::vector<double> fGen_part_P0_x, fGen_part_P0_y, fGen_part_P0_z;                    // Final position of the particle

  //Geant4 Particles variables
  std::vector<int> fG4_part_trackID;    // ID given to the particle
  std::vector<int> fG4_part_mother;     // Mother of the particle
  std::vector<int> fG4_part_PDGcode;    // PDGcode of the particle
  std::vector<double> fG4_part_mass;    // mass of the particle
  std::vector<double> fG4_part_E0;      // Initial energy of the particle
  std::vector<double> fG4_part_Ef;      // Final energyh of the particle
    
  std::vector<double> fG4_part_StartPos_x, fG4_part_StartPos_y, fG4_part_StartPos_z;    // Initial position of the particle
  std::vector<double> fG4_part_EndPos_x, fG4_part_EndPos_y, fG4_part_EndPos_z;          // Final position of the particle

  std::vector<double> fG4_part_P0_x, fG4_part_P0_y, fG4_part_P0_z;    // Initial and final momentum of the particle
  std::vector<double> fG4_part_Pf_x, fG4_part_Pf_y, fG4_part_Pf_z;    // Final momentum of the particle
  
  std::vector<std::string> fG4_part_process;      // Process that creates the particle
  std::vector<std::string> fG4_part_endProcess;   // Process that makes the particle stop propagating
  
  //Slices related parameters / pandora variables
  unsigned int fReco_event_nSlices;              // Number of slices in the event

  unsigned int fReco_slice_ID;                   // ID of the slice
  unsigned int fReco_slice_nPFParticles;         // Number of particles (PFPs) in the slice
  unsigned int fReco_slice_prim_ID;              // ID of the primnary particle (neutrino, if not cosmic)
  unsigned int fReco_slice_prim_PDG;             // pdg (12 or 14) of the primary particle 
  unsigned int fReco_slice_prim_nDaughters;      // number of daughters of the primary particle. In case of a neutrino interaction these are the primary particles

  unsigned int fReco_slice_nTracks;              // Number of tracks in the slice
  unsigned int fReco_slice_nShowers;             // Number of showers in the slice
  unsigned int fReco_slice_nPrimTracks;          // Number of primary tracks in the slice
  unsigned int fReco_slice_nPrimShowers;         // Number of primary showers in the slice

  double fReco_slice_vertexPos_x, fReco_slice_vertexPos_y, fReco_slice_vertexPos_z;    // Reconstructed vertex of the slice

  std::vector<int> fReco_part_ID;                   // ID of the particle of the slice
  std::vector<int> fReco_part_mother;               // Mother of the particle of the slice
  std::vector<int> fReco_part_PDGCode;              // PDGcode of the 'neutrino' particle of the slice. It only says if it is electron-like or muon-like trace
  std::vector<bool> fReco_part_isPrimaryChildren;   // true for primary particles (daughters of the neutrino)

  std::vector<double> fReco_part_vertexPos_x, fReco_part_vertexPos_y, fReco_part_vertexPos_z;        // Vertex of the particle of the slice

  // Track variables
  std::vector<double> fReco_part_track_score;       // Track score of the particle of the slice
  double fReco_nu_score;                            // Neutrino score of the particle of the slice
  bool fReco_is_clear_cosmic;                       // Says if a slice corresponds to a clear cosmic

  std::vector<double> fReco_part_track_length;      // Track length of the particle of the slice
  std::vector<double> fReco_part_track_theta;       // AAA Angle of the particle with AAA
  std::vector<double> fReco_part_track_phi;         // AAA Angle of the particle with AAA

  std::vector<double> fReco_part_track_startPos_x, fReco_part_track_startPos_y, fReco_part_track_startPos_z;    // Start position of the particle of the slice
  std::vector<double> fReco_part_track_EndPos_x, fReco_part_track_EndPos_y, fReco_part_track_EndPos_z;          // End position of the particle of the slice
  std::vector<double> fReco_part_startDir_x, fReco_part_startDir_y, fReco_part_startDir_z;                      // (initial) Direction of the particle of the slice

  // Calorimetry variables --> they will be empty for the neutrino because it does not leave a trace
  std::vector<double> fReco_track_kineticEnergy;       // Kinetic energ of the particle of the slice
  std::vector<double> fReco_track_visibleEnergy;       // Should also be the kinetic energy of the particle of the slice

  std::vector<float> fReco_track_completeness;                 // Completeness of the track of the particle of the slice
  std::vector<float> fReco_track_purity;                       // Purity of the track of the particle of the slice

  std::vector<std::vector<float>> fReco_track_dEdx;               // dE/dx for each step of the deposition of energy of the particle of the slice
  std::vector<std::vector<float>> fReco_track_residualRange;      // AAA length (in cm) of each step of the deposition of energy of the particle of the slice
  std::vector<std::vector<float>> fReco_track_pitch;              // AAA length (in cm) of each step of the deposition of energy of the particle of the slice
  std::vector<float> fReco_track_range;                           // AAA total length of the track of a particle

  std::vector<float> fReco_part_chi2_muon;            // chi^2 of the deposition of energy of the particle compared with a muon
  std::vector<float> fReco_part_chi2_pion;            // chi^2 of the deposition of energy of the particle compared with a pion
  std::vector<float> fReco_part_chi2_kaon;            // chi^2 of the deposition of energy of the particle compared with a kaon
  std::vector<float> fReco_part_chi2_proton;          // chi^2 of the deposition of energy of the particle compared with a proton

  // Shower variables --> they will be empty for the neutrino because it does not leave a trace
  std::vector<double> fReco_shower_dir_x, fReco_shower_dir_y, fReco_shower_dir_z;                     // Direction of the shower of the slice

  std::vector<double> fReco_shower_start_x, fReco_shower_start_y, fReco_shower_start_z;               // Start position of the shower of the slice

  std::vector<double> fReco_shower_openAngle;     // Opening angle of the shower of the slice
  std::vector<double> fReco_shower_length;        // Length of the shower of the slice
  std::vector<double> fReco_shower_energy;        // Energy deposited by the shower of the slice

  // Truth matching variables
  std::vector<int> fReco_track_trueTrackID;           // True ID IN g4 of the particle of the slice. ALWAYS USE IN ABSOLUTE VALUE
  std::vector<int> fReco_track_truePDG;               // True PDGcode of the particle of the slice
  std::vector<double> fReco_true_E0;              // True initial Energy of the particle of the slice
  std::vector<double> fReco_true_Ef;              // True final Energy of the particle of the slice
  std::vector<double> fReco_true_p0_x, fReco_true_p0_y, fReco_true_p0_z;   // True momentum of the particle of the slice
  std::vector<double> fReco_true_pf_x, fReco_true_pf_y, fReco_true_pf_z;   // True momentum of the particle of the slice
  std::vector<std::string> fReco_true_endProcess;   // True Process that makes the particle stop propagating
  std::vector<double> fReco_true_StartPos_x, fReco_true_StartPos_y, fReco_true_StartPos_z;   // True start position of the particle
  std::vector<bool> fReco_hasTrack;                                                         // Variable that says if the particle has a track (photons, neutrinos do not leave track)
  std::vector<bool> fReco_isNeutrino;                                                      // Variable that says is this PFP track is the neutrino

  // Maps --> YET TO KNOW WHAT THEY ALL DO
  std::map<int, int> fHitsMap;                                                    // track ID , nHits
  std::map<const art::Ptr<simb::MCTruth>, int> fMCTruthHitsMap;                   // MCTruth object, nHits of the truth particle
  std::map<int, art::Ptr<recob::PFParticle>> fPFPMap;                             // PFP ID, PFP
  std::map<int, std::set<art::Ptr<recob::PFParticle>>> fRecoPFPMap;               // NOT USED

  // Declaration of fhicl parameters
  const std::string fMCTruthLabel;
  const std::string fMCParticleLabel;
  const std::string fPOTModuleLabel;
  const std::string fSliceLabel;
  const std::string fSliceModuleLabel;
  const std::string fPFParticleLabel;
  const std::string fPFParticleMetadataLabel;
  const std::string fVertexLabel;
  const std::string fTrackLabel;
  const std::string fShowerLabel;
  const std::string fCalorimetryLabel;
  const std::string fChi2ModuleLabel;
  const std::string fHitLabel;

  const bool fVerbose;

  const bool fSave_truth_gen_particles = false;
  const bool fSave_truth_g4 = false;
  const bool fSave_reco = false;
  const bool fSave_TruthMatching = false;

  const double fTrackScoreLimit;

};


sbnd::AnalyzeReco::AnalyzeReco(fhicl::ParameterSet const& p)                     // This is the constructor of the class
  : EDAnalyzer{p},

  //From here, I get all the labels
    fMCTruthLabel(p.get<std::string>("MCTruthLabel")) ,
    fMCParticleLabel(p.get<std::string>("MCParticleLabel")),
    fPOTModuleLabel(p.get<std::string>("POTLabel")),
    fSliceLabel(p.get<std::string>("SliceLabel")) ,
    fSliceModuleLabel(p.get<std::string>("SliceModuleLabel")) ,
    fPFParticleLabel(p.get<std::string>("PFParticleLabel")) ,
    fPFParticleMetadataLabel(p.get<std::string>("PFParticleMetadataLabel")) ,
    fVertexLabel(p.get<std::string>("VertexLabel")) ,
    fTrackLabel(p.get<std::string>("TrackLabel")) ,
    fShowerLabel(p.get<std::string>("ShowerLabel")) ,
    fCalorimetryLabel(p.get<std::string>("CalorimetryLabel")) ,
    fChi2ModuleLabel(p.get<std::string>("Chi2ModuleLabel")) ,
    fHitLabel(p.get<std::string>("HitLabel")) ,

    fVerbose(p.get<bool>("Verbose")),

    fSave_truth_gen_particles(p.get<bool>("save_truth_gen_particles")),
    fSave_truth_g4(p.get<bool>("save_truth_g4")),
    fSave_reco(p.get<bool>("save_reco")),
    fSave_TruthMatching(p.get<bool>("save_TruthMatching")),

    fTrackScoreLimit(p.get<double>("TrackScoreLimit"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::AnalyzeReco::analyze(art::Event const& e)                             // This is the analyze method/function. it is called for every event
{
  // Implementation of required member function here.

  // I get the event ID, subrun ID and  total number of events
  fEvent_ID = e.id().event();
  fRun_ID = e.id().run();
  fSubrun_ID = e.id().subRun();
  fNum_gen_evts = getTotalGenEvents(e);     // In this case I am using a method that I have defined for this

  if (fVerbose) {
    std::cout << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl;
    std::cout << "Event number: " << fEvent_ID << std::endl;
    std::cout << "Run number: " << fRun_ID << std::endl;
    std::cout << "Subrun number: " << fSubrun_ID;
    std::cout << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl;
  }


  // I get access to the handles of the different data product that I will be using
  // Truth handles
  // generator
  art::ValidHandle<std::vector<simb::MCTruth>> MCT_handle = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);                // generator
  // g4
  art::ValidHandle<std::vector<simb::MCParticle>> MCP_handle = e.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleLabel);       // largeant

  // Pandora (reco) handles
  // Slices
  art::ValidHandle<std::vector<recob::Slice>> slice_handle = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);                  // pandora
  // PFParticles
  art::ValidHandle<std::vector<recob::PFParticle>> PFP_handle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);   // pandora
  // Tracks
  art::ValidHandle<std::vector<recob::Track>> track_handle = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);                  // pandoraTrack
  // Showers
  art::ValidHandle<std::vector<recob::Shower>> shower_handle = e.getValidHandle<std::vector<recob::Shower>>(fShowerLabel);                // pandoraShowerSBN 
  // Vertex
  art::ValidHandle<std::vector<recob::Vertex>> vertex_handle = e.getValidHandle<std::vector<recob::Vertex>>(fVertexLabel);                // pandora 
  // Hits
  art::ValidHandle<std::vector<recob::Hit>> hit_handle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);                          // gaushit

  // Now I get the associations using FindManyP
  // Association between slices and PFParticles
  art::FindManyP<recob::PFParticle> SlicePfpAssoc(slice_handle, e, fPFParticleLabel);
  // Association between PFP and vertex
  art::FindOneP<recob::Vertex> PfpVertexAssoc(PFP_handle, e, fVertexLabel);
  // Association between shower and hit
  art::FindManyP<recob::Hit> ShowerHitAssoc(shower_handle, e, fShowerLabel);
  // Association between slice and hit
  art::FindManyP<recob::Hit> SliceHitAssoc(slice_handle, e, fSliceModuleLabel);

  // Get the vectors of gen, g4, slices
  // Get gen vector
  std::vector<art::Ptr < simb::MCTruth>> MCT_vec;
  art::fill_ptr_vector(MCT_vec, MCT_handle);
  

  // Get g4 vector
  std::vector<art::Ptr < simb::MCParticle>> MCP_vec;
  art::fill_ptr_vector(MCP_vec, MCP_handle);
  
  //Get slice vector
  std::vector<art::Ptr < recob::Slice>> slice_vec;
  art::fill_ptr_vector(slice_vec, slice_handle);
  // Other associated vectors will be found inside the functions

  /*
  THINGS THAT I HAVE NOT DONE YET
  - Create the methods to clear all variables and to empty them
  - Write the config.fcl with all the labels and variables
  - Make sure the analyzer works if some module is missing in the analyzed .root file (for that, make sure to have the correct ifs when accessing each module)
  - 
  - Write the analysis methods for reco
  */


  // MAKE A LOOP OVER AL SLICES
  /*
  What I will do:
  - loop over all slices
  - Save the reco info of each slice: tracks, showers, hits, ...
  - I also have to save the info of truth information: I CAN'T DO IT FOR EACH SLICE or I would repeat information way too many times
  - I will save truth information ONLY when the events seems to be a neutrino event:
    analyseTruth_gen(MCT_vec) and analyseTruth_g4(MCP_vec) methods will only appear when the events seem to be a neutrino event
    Therefore, I need a 'is_clear_cosmic' variable
      --> I will analyze the generator info (analyseTruth_gen(MCT_vec)) of the:
          1. Neutrino events
          2. Matches the best (on truth matching) a certain slice (because this analyzer will work for spill, each MCT_vec might have more than 1 entry, it might have 1 entry for each neutrino event in the same spill)
      --> I will analyze the g4 info only when the slice seems to be a neutrino event. I will use 'fReco_is_clear_cosmic' to judge that.
  */


  /*
  WORKFLOW

  1. I save basic information (run_ID, event_ID, ...)
  2. I loop over all the slices
    2.1. I save the g4 information
    2.2. I save basic slice information (slice ID, number of PFPs in the slice, ...)
    2.3. I save basic information about the primary particle (ID, pdg, number of daughters, ...)
  3. I save the reco information of all slices (PFPs, track, shower, hits, calo)
    3.1. I save track info
      3.1.1. I save the calorimetry
      3.1.2. I save the TruthMatching of the tracks
        3.1.2.1. I save the generator info in for the best MCTruth fit to the reco slice that is found
    3.2. I save shower info
  4. Fill the TTree of this slice

  FINISH
  */

 // I loop over all the slices of each event:
 for (art::Ptr<recob::Slice> &slice : slice_vec){
  // I clear the variables first
  clearVars();
  clearMaps();

  // Get the slice ID
  fReco_slice_ID = slice->ID();

  // I analyze the g4 stage
  analyseTruth_g4(MCP_vec);

  // Get the PFPs in the slice using the assciation betwen slice and PFPs
  std::vector<art::Ptr<recob::PFParticle>> PFP_vec(SlicePfpAssoc.at(slice.key()));
  fReco_slice_nPFParticles = PFP_vec.size();

  // I make sure I only continue if the slice is not empty
  if(PFP_vec.size()!=0){

    // Define 'prim'--> the primary particle of the slice (the neutrino, if there are no cosmics)
    art::Ptr<recob::PFParticle> prim = getPrimaryPFP(PFP_vec);
    // I make sure we only continue if a primary particle is found
    if(!prim.isNull()){

      // I save the primary particle information
      fReco_slice_prim_ID = prim->Self();
      fReco_slice_prim_PDG = prim->PdgCode();
      fReco_slice_prim_nDaughters = prim->NumDaughters();

      // Save the primary vertex position
      art::Ptr<recob::Vertex> vertex = PfpVertexAssoc.at(prim.key());
      if(vertex.isNonnull()){
        fReco_slice_vertexPos_x = vertex->position().X();
        fReco_slice_vertexPos_y = vertex->position().Y();
        fReco_slice_vertexPos_z = vertex->position().Z();
      }

      // Setup the maps
      setupMaps(e, hit_handle, PFP_handle);

      // Analyze the data: PFPs, tracks, showers, calorimetry, truth matching and generator
      analysePFPs(e,MCT_vec, MCP_vec, prim, PFP_vec, PFP_handle, track_handle, shower_handle);

    }// Checking if there is a primary particle
  }// Checking the size of the vector of PFPs is not 0

  // At the end, I have to fill the trees. I will have one entry of the tree for each slice
  fTree->Fill();

 }// Finish of the loop over slices
}

//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Auxiliary functions //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// Reset the values of all variables and empty vectors
void sbnd::AnalyzeReco::clearVars()
{
  //Set default neutrinos values
    fNu_PDG = -1;
    fNu_E0 = -1;
    fNu_weight = -1;
    fNu_interaction_mode = -1;
    fNu_interaction_type = -1;
    fNu_CC_NC = -1;
    fNu_target = -1;
    fNu_HitNuc = -1;
    fNu_HitQuark = -1;
    fNu_W = -1;
    fNu_X = -1;
    fNu_Y = -1;
    fNu_Q2 = -1;
    fGen_index = -1;

    //Reset the vectors with MCTruth particles information
    fGen_part_trackID.clear();
    fGen_part_statusCode.clear();
    fGen_part_mother.clear();
    fGen_part_PDGcode.clear();
    fGen_part_mass.clear();
    fGen_part_E0.clear();
    fGen_part_StartPos_x.clear();
    fGen_part_StartPos_y.clear();
    fGen_part_StartPos_z.clear();
    fGen_part_P0_x.clear();
    fGen_part_P0_y.clear();
    fGen_part_P0_z.clear();


    //Reset the vectors with geant4 particles information
    fG4_part_trackID.clear();
    fG4_part_mother.clear();
    fG4_part_PDGcode.clear();
    fG4_part_mass.clear();
    fG4_part_E0.clear();
    fG4_part_Ef.clear();

    fG4_part_StartPos_x.clear();
    fG4_part_StartPos_y.clear();
    fG4_part_StartPos_z.clear();
    fG4_part_EndPos_x.clear();
    fG4_part_EndPos_y.clear();
    fG4_part_EndPos_z.clear();

    fG4_part_P0_x.clear();
    fG4_part_P0_y.clear();
    fG4_part_P0_z.clear();
    fG4_part_Pf_x.clear();
    fG4_part_Pf_y.clear();
    fG4_part_Pf_z.clear();

    fG4_part_process.clear();
    fG4_part_endProcess.clear();

    //Reset reco particle vectors
    fReco_slice_ID = -1;
    fReco_slice_nPFParticles = 0;
    fReco_slice_prim_ID = 0;
    fReco_slice_prim_PDG = 0;
    fReco_slice_prim_nDaughters = 0;

    fReco_slice_nTracks = 0;
    fReco_slice_nShowers = 0;
    fReco_slice_nPrimTracks = 0;
    fReco_slice_nPrimShowers = 0;

    fReco_slice_vertexPos_x = -999;
    fReco_slice_vertexPos_y = -999;
    fReco_slice_vertexPos_z = -999;

    fReco_part_ID.clear();
    fReco_part_mother.clear();
    fReco_part_PDGCode.clear();
    fReco_part_isPrimaryChildren.clear();

    fReco_part_vertexPos_x.clear();
    fReco_part_vertexPos_y.clear();
    fReco_part_vertexPos_z.clear();

    fReco_part_track_score.clear();
    fReco_nu_score = -1;
    fReco_is_clear_cosmic=false;
    fReco_part_track_length.clear();
    fReco_part_track_theta.clear();
    fReco_part_track_phi.clear();

    fReco_part_track_startPos_x.clear();
    fReco_part_track_startPos_y.clear();
    fReco_part_track_startPos_z.clear();
    fReco_part_track_EndPos_x.clear();
    fReco_part_track_EndPos_y.clear();
    fReco_part_track_EndPos_z.clear();

    fReco_part_startDir_x.clear();
    fReco_part_startDir_y.clear();
    fReco_part_startDir_z.clear();

    fReco_track_kineticEnergy.clear();
    fReco_track_visibleEnergy.clear();
    fReco_track_dEdx.clear();
    fReco_track_residualRange.clear();
    fReco_track_pitch.clear();
    fReco_track_range.clear();

    fReco_part_chi2_muon.clear();
    fReco_part_chi2_pion.clear();
    fReco_part_chi2_kaon.clear();
    fReco_part_chi2_proton.clear();

    fReco_track_trueTrackID.clear();
    fReco_track_truePDG.clear();
    fReco_true_E0.clear();
    fReco_true_Ef.clear();
    fReco_true_p0_x.clear();
    fReco_true_p0_y.clear();
    fReco_true_p0_z.clear();
    fReco_true_pf_x.clear();
    fReco_true_pf_y.clear();
    fReco_true_pf_z.clear();
    fReco_true_endProcess.clear();
    fReco_true_StartPos_x.clear();
    fReco_true_StartPos_y.clear();
    fReco_true_StartPos_z.clear();
    fReco_hasTrack.clear();
    fReco_isNeutrino.clear();   

    fReco_track_completeness.clear();
    fReco_track_purity.clear();

    fReco_shower_dir_x.clear();
    fReco_shower_dir_y.clear();
    fReco_shower_dir_z.clear();

    fReco_shower_start_x.clear();
    fReco_shower_start_y.clear();
    fReco_shower_start_z.clear();

    fReco_shower_openAngle.clear();
    fReco_shower_length.clear();
    fReco_shower_energy.clear();
}

// Reset subrun variables
void sbnd::AnalyzeReco::clearSubRunVars()
{
  fPOT = 0.; 
  fSpill = 0; 
  fSubrun_ID = 0;
  fNum_gen_evts = 0;
}

// Set a default value for gen vector variables
void sbnd::AnalyzeReco::setDefaultGenVars()
{
  fGen_part_trackID.push_back(-1);
  fGen_part_statusCode.push_back(-1);
  fGen_part_mother.push_back(-1);
  fGen_part_PDGcode.push_back(-1);
  fGen_part_mass.push_back(-999);
  fGen_part_E0.push_back(-999);
  fGen_part_StartPos_x.push_back(-999);
  fGen_part_StartPos_y.push_back(-999);
  fGen_part_StartPos_z.push_back(-999);
  fGen_part_P0_x.push_back(-999);
  fGen_part_P0_y.push_back(-999);
  fGen_part_P0_z.push_back(-999);
}

// Set a default value for g4 vector variables
void sbnd::AnalyzeReco::setDefaultG4Vars()
{
  fG4_part_trackID.push_back(-1);
  fG4_part_mother.push_back(-1);
  fG4_part_PDGcode.push_back(-1);
  fG4_part_mass.push_back(-999);
  fG4_part_E0.push_back(-999);
  fG4_part_Ef.push_back(-999);

  fG4_part_StartPos_x.push_back(-999);
  fG4_part_StartPos_y.push_back(-999);
  fG4_part_StartPos_z.push_back(-999);
  fG4_part_EndPos_x.push_back(-999);
  fG4_part_EndPos_y.push_back(-999);
  fG4_part_EndPos_z.push_back(-999);

  fG4_part_P0_x.push_back(-999);
  fG4_part_P0_y.push_back(-999);
  fG4_part_P0_z.push_back(-999);
  fG4_part_Pf_x.push_back(-999);
  fG4_part_Pf_y.push_back(-999);
  fG4_part_Pf_z.push_back(-999);

  fG4_part_process.push_back("");
  fG4_part_endProcess.push_back("");
}

// Set a default value for reco vector variables
void sbnd::AnalyzeReco::setDefaultRecoVars()
{
  fReco_part_ID.push_back(-1);
  fReco_part_mother.push_back(-1);
  fReco_part_PDGCode.push_back(-1);
  fReco_part_isPrimaryChildren.push_back(false);

  fReco_part_track_score.push_back(-1);
  fReco_part_vertexPos_x.push_back(-999);
  fReco_part_vertexPos_y.push_back(-999);
  fReco_part_vertexPos_z.push_back(-999);

  fReco_part_track_length.push_back(-999);
  fReco_part_track_theta.push_back(-999);
  fReco_part_track_phi.push_back(-999);


  fReco_part_startDir_x.push_back(-999);
  fReco_part_startDir_y.push_back(-999);
  fReco_part_startDir_z.push_back(-999);

  fReco_part_track_startPos_x.push_back(-999);
  fReco_part_track_startPos_y.push_back(-999);
  fReco_part_track_startPos_z.push_back(-999);
  fReco_part_track_EndPos_x.push_back(-999);
  fReco_part_track_EndPos_y.push_back(-999);
  fReco_part_track_EndPos_z.push_back(-999);

  fReco_track_kineticEnergy.push_back(-999);
  fReco_track_visibleEnergy.push_back(-999);
  fReco_track_dEdx.push_back({});
  fReco_track_residualRange.push_back({});
  fReco_track_pitch.push_back({});
  fReco_track_range.push_back(-999);

  fReco_track_trueTrackID.push_back(-1);
  fReco_track_truePDG.push_back(-1);
  fReco_true_E0.push_back(-999);
  fReco_true_Ef.push_back(-999);
  fReco_true_p0_x.push_back(-999);
  fReco_true_p0_y.push_back(-999);
  fReco_true_p0_z.push_back(-999);
  fReco_true_pf_x.push_back(-999);
  fReco_true_pf_y.push_back(-999);
  fReco_true_pf_z.push_back(-999);
  fReco_true_endProcess.push_back("");
  fReco_true_StartPos_x.push_back(-999);
  fReco_true_StartPos_y.push_back(-999);
  fReco_true_StartPos_z.push_back(-999);
  fReco_hasTrack.push_back(false);     // true for default, I only change it to false if it the track is found
  fReco_isNeutrino.push_back(false);     // true for default, I only change it to false if it the track is found

  fReco_track_completeness.push_back(-1);
  fReco_track_purity.push_back(-1);

  fReco_part_chi2_muon.push_back(-1);
  fReco_part_chi2_pion.push_back(-1);
  fReco_part_chi2_kaon.push_back(-1);
  fReco_part_chi2_proton.push_back(-1);

  fReco_shower_dir_x.push_back(-999);
  fReco_shower_dir_y.push_back(-999);
  fReco_shower_dir_z.push_back(-999);

  fReco_shower_start_x.push_back(-999);
  fReco_shower_start_y.push_back(-999);
  fReco_shower_start_z.push_back(-999);

  fReco_shower_openAngle.push_back(-999);
  fReco_shower_length.push_back(-999);

  fReco_shower_energy.push_back(-999);
}

void sbnd::AnalyzeReco::clearMaps()
{
  fHitsMap.clear();
  fMCTruthHitsMap.clear();
  fPFPMap.clear();
  fRecoPFPMap.clear();
}


void sbnd::AnalyzeReco::setupMaps(const art::Event &e, const art::ValidHandle<std::vector<recob::Hit>> &hit_handle,
                                       const art::ValidHandle<std::vector<recob::PFParticle>> &PFP_handle)
                                       {
  // I get the neccesary info for truth matching
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  // I create the vector with all the hits 
  std::vector<art::Ptr<recob::Hit>> hit_vec;
  art::fill_ptr_vector(hit_vec, hit_handle);
  // I loop over all the hits
  for(auto const& hit : hit_vec){
    // I find the trackID of the truth particle that made the hit
    const int trackID = TruthMatchUtils::TrueParticleID(clockData,hit,true);
    // I fill the map: trackID of the particle that created the hit --> number of hits
    // Counts how many hits are associated with each true particle.
    fHitsMap[trackID]++;

    // I define an MCTruth object
    art::Ptr<simb::MCTruth> MCT;

    // I check that the trackID that was found is correct. In case it is correct, I fill the MCT object
    if (trackID == def_int) {
      MCT = art::Ptr<simb::MCTruth>(); // empty pointer
    }
    else {
      MCT = particleInv->TrackIdToMCTruth_P(trackID);   // MCT object
    }

    // I fill the map fMCTruthHitsMap of the MCT object
    // Counts how many hits are associated with each MCTruth interaction
    fMCTruthHitsMap[MCT]++;
  }

  // Create the vector of PFPs
  std::vector<art::Ptr<recob::PFParticle>> PFP_vec;
  art::fill_ptr_vector(PFP_vec, PFP_handle);

  // I loop over all PFParticles
  for(auto const& PFP : PFP_vec)
    // It gives to each ID of a PFParticle the corresponding PFParticle. This makes it easy to look up the PFParticle pointer given just its ID
    fPFPMap[PFP->Self()] = PFP;

}

// Method to find the total number of events inside the file
/*
Receives: the event 'e'
Returns: the total number of events
*/
int sbnd::AnalyzeReco::getTotalGenEvents(const art::Event &e)
{
  int nGenEvt = 0; // This will count the event number
  for (const art::ProcessConfiguration &process: e.processHistory()) {                 // Loop over all the history of all events
    std::optional<fhicl::ParameterSet> genConfig = e.getProcessParameterSet(process.processName());
    if (genConfig && genConfig->has_key("source") && genConfig->has_key("source.maxEvents") && genConfig->has_key("source.module_type") ) {
      int maxEvents = genConfig->get<int>("source.maxEvents");                        // Get the maximum number of events
      std::string moduleType = genConfig->get<std::string>("source.module_type");     // When the 'moduleType' is an empty event
      if (moduleType == "EmptyEvent") {
        nGenEvt += maxEvents;               // sum the events
      }
    }
  }

  return nGenEvt;
}

// Method to calculate the completeness of a trace. If there are no hits --> purity = def_float (happens for photons)
float sbnd::AnalyzeReco::completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID){

  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::map<int, int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

  float comp = (fHitsMap[trackID] == 0) ? def_float : objectHitsMap[trackID]/static_cast<float>(fHitsMap[trackID]);

  if(fVerbose){
    std::cout << "completeness: " << comp << std::endl;
  }
  return comp;
}

// Method to calculate the purity of a trace. If there are no hits --> purity = 0 (happens for photons)
float sbnd::AnalyzeReco::purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID){
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::map<int, int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

  float pur = (objectHits.size() == 0) ? def_float : objectHitsMap[trackID]/static_cast<float>(objectHits.size());

  if(fVerbose){
    std::cout << "purity: " << pur << std::endl;
  }

  return pur;
}

// Method to calculate the chi2 with different particles
void sbnd::AnalyzeReco::getChi2PID(const art::Ptr<anab::ParticleID> &chi2pid){

  const std::vector<anab::sParticleIDAlgScores> AlgScoresVec = chi2pid->ParticleIDAlgScores();
  std::vector<std::pair<int, float>> chi2s;

  for(size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
  {
    const anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);

    if(AlgScore.fAlgName == "Chi2")
    {
      switch(AlgScore.fAssumedPdg)
      {
        case 13:
          fReco_part_chi2_muon.back() = AlgScore.fValue;
          break;
        case 211:
          fReco_part_chi2_pion.back() = AlgScore.fValue;
          break;
        case 321:
          fReco_part_chi2_kaon.back() = AlgScore.fValue;
          break;
        case 2212:
          fReco_part_chi2_proton.back() = AlgScore.fValue;
          break;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Analyze functions //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void sbnd::AnalyzeReco::analyseTruth_gen(const art::Ptr<simb::MCTruth> MCT)  
{
  // Get generator information
  // Get neutrino information

  fNu_PDG = MCT->GetNeutrino().Nu().PdgCode();
  fNu_E0 = MCT->GetNeutrino().Nu().E(0);
  fNu_weight = MCT->GetNeutrino().Nu().Weight();

  fNu_interaction_mode = MCT->GetNeutrino().Mode();
  fNu_interaction_type = MCT->GetNeutrino().InteractionType();
  fNu_CC_NC = MCT->GetNeutrino().CCNC();

  fNu_HitNuc = MCT->GetNeutrino().HitNuc();
  fNu_target = MCT->GetNeutrino().Target();
  fNu_HitQuark = MCT->GetNeutrino().HitQuark();
  fNu_W = MCT->GetNeutrino().W();
  fNu_X = MCT->GetNeutrino().X();
  fNu_Y = MCT->GetNeutrino().Y();
  fNu_Q2 = MCT->GetNeutrino().QSqr();

  fGen_index = MCT->Origin();

  if (fVerbose) {
    std::cout << "---- event " << fEvent_ID << " ----" << std::endl;
    std::cout << "Neutrino Information: " << std::endl;
    std::cout << "Neutrino properties: " << std::endl;
    std::cout << "PDG: " << fNu_PDG << " Energy: " << fNu_E0 << " weight: " << fNu_weight << std::endl;
    std::cout << "Interaction mode: " << fNu_interaction_mode << " Interaction type: " << fNu_interaction_type << " CCorNC: " << fNu_CC_NC << std::endl;
    std::cout << "Target: " << fNu_target << " HitNuc: " << fNu_HitNuc << " HitQuark: " << fNu_HitQuark << std::endl;
    std::cout << " W: " << fNu_W << " X: " << fNu_X << " Y: " << fNu_Y << " Qsqr: " << fNu_Q2 << std::endl;
    std::cout << std::endl;
  }

  // Get particles information: I loop over all the particles
  for(int i=0; i<MCT->NParticles(); i++){
    setDefaultGenVars();

    simb::MCParticle MCTpart = MCT->GetParticle(i);   // I store in the variable MCTpart all the information of one particle

    fGen_part_trackID.back() = MCTpart.TrackId();
    fGen_part_statusCode.back() = MCTpart.StatusCode();
    fGen_part_mother.back() = MCTpart.Mother();
    fGen_part_PDGcode.back() = MCTpart.PdgCode();
    fGen_part_mass.back() = MCTpart.Mass();
    fGen_part_E0.back() = MCTpart.E(0);
    fGen_part_StartPos_x.back() = MCTpart.Vx(0);
    fGen_part_StartPos_y.back() = MCTpart.Vy(0);
    fGen_part_StartPos_z.back() = MCTpart.Vz(0);
    fGen_part_P0_x.back() = MCTpart.Px(0);
    fGen_part_P0_y.back() = MCTpart.Py(0);
    fGen_part_P0_z.back() = MCTpart.Pz(0);

    if (fVerbose) {
      std::cout << "ID: " << MCTpart.TrackId() << " Status code: "<< MCTpart.StatusCode() <<  " Mother: " << MCTpart.Mother()<< std::endl;
      std::cout << "   PDGCode: " << MCTpart.PdgCode() << " Mass: " << MCTpart.Mass() << " Energy: " << MCTpart.E(0) << std::endl;
      std::cout << "   Start Pos (x,y,z) : (" << MCTpart.Vx(0) << ", " << MCTpart.Vy(0) << ", " << MCTpart.Vz(0)  << ")" << std::endl;
      std::cout << "   Start Momentum (x,y,z) : (" << MCTpart.Px(0) << ", " << MCTpart.Py(0) << ", " << MCTpart.Pz(0)  << ")" <<std::endl;
    }
  }
}

void sbnd::AnalyzeReco::analyseTruth_g4(const std::vector<art::Ptr < simb::MCParticle>> MCP_vec)
{
  // Get g4 information 
  for (const art::Ptr <simb::MCParticle> &MCP: MCP_vec) {    // I loop over all the information in the vector. This vector runs over all particles
    setDefaultG4Vars();

    // I get the number of trajectory points (the number of depositions of energy)
    const size_t numberTrajectoryPoints = MCP->NumberTrajectoryPoints();
    const int last = numberTrajectoryPoints - 1;                            // I will use this to access the characteristics of a particle at the end of its propagation

    if(numberTrajectoryPoints>0){
      // I store all the variables
      fG4_part_trackID.back() = MCP->TrackId();
      fG4_part_mother.back() = MCP->Mother();
      fG4_part_PDGcode.back() = MCP->PdgCode();
      fG4_part_mass.back() = MCP->Mass();
      fG4_part_E0.back() = MCP->E(0);
      fG4_part_Ef.back() = MCP->E(last);

      fG4_part_StartPos_x.back() = MCP->Vx(0);
      fG4_part_StartPos_y.back() = MCP->Vy(0);
      fG4_part_StartPos_z.back() = MCP->Vz(0);
      fG4_part_EndPos_x.back() = MCP->Vx(last);
      fG4_part_EndPos_y.back() = MCP->Vy(last);
      fG4_part_EndPos_z.back() = MCP->Vz(last);

      fG4_part_P0_x.back() = MCP->Px(0);
      fG4_part_P0_y.back() = MCP->Py(0);
      fG4_part_P0_z.back() = MCP->Pz(0);
      fG4_part_Pf_x.back() = MCP->Px(last);
      fG4_part_Pf_y.back() = MCP->Py(last);
      fG4_part_Pf_z.back() = MCP->Pz(last);

      fG4_part_process.back() = MCP->Process();
      fG4_part_endProcess.back() = MCP->EndProcess();
    }
    
    if (fVerbose) {
      std::cout << "ID: " << MCP->TrackId()<< " Mother: " << MCP->Mother()<<std::endl;
      std::cout << "   PDGCode: " << MCP->PdgCode()<< " Mass: " << MCP->Mass()
                << " Initial Energy: " << MCP->E(0) << " Final Energy: " << MCP->E(last) << std::endl;
      std::cout << "   Start Pos (x,y,z,t) : (" <<  MCP->Vx(0) << ", " << MCP->Vy(0) << ", " << MCP->Vz(0) << ", " << MCP->T(0) << ")" <<std::endl;
      std::cout << "   End Pos (x,y,z,t) : (" << MCP->Vx(last)<< ", " << MCP->Vy(last)<< ", "<< MCP->Vz(last)<< ", " << MCP->T(last)<< ")" <<std::endl;
      std::cout << "   Start Momentum (x,y,z) : (" << MCP->Px(0) << ", " << MCP->Py(0) << ", "<< MCP->Pz(0) << ")" <<std::endl;
      std::cout << "   End Momentum (x,y,z) : (" << MCP->Px(last)<< ", " << MCP->Py(last)<< ", "<< MCP->Pz(last)<< ")" <<std::endl;
      std::cout << "   Process: " << MCP->Process()<<std::endl;
      std::cout << "   End_process: " << MCP->EndProcess()<< std::endl;
      std::cout << std::endl;
    }
  }
}


/*
analysePFPs
RECEIVES:
- The art::Event 'e'
- The MCTruth vector to run it during the TruthMatching
- The MCParticle vector to run it after finding the is_clear_cosmic
- The PFParticle 'prim' is the neutrino (or, in case of cosmics, the primary cosmic). 
  Pandora selects its pdg (although it may be wrong). All other particles are its daughters
- The vector with ALL PFParticles PFP_vec
- The handle for the class recob::Track
- The handle for the class recob::Shower
*/
void sbnd::AnalyzeReco::analysePFPs(const art::Event &e, 
                    const std::vector<art::Ptr<simb::MCTruth>> MCT_vec,
                    const std::vector<art::Ptr < simb::MCParticle>> MCP_vec,
                    const art::Ptr<recob::PFParticle> &prim,
                    const std::vector<art::Ptr<recob::PFParticle>> &PFP_vec,
                    const art::ValidHandle<std::vector<recob::PFParticle>> &PFP_handle,
                    const art::ValidHandle<std::vector<recob::Track>> &track_handle,
                    const art::ValidHandle<std::vector<recob::Shower>> &shower_handle)
{
  // I get the assoociations that I need:
  // FindOneP --> If we know one object will ONLY have one associated object of that class (one PFP will only have 1 track)
  // FindManyP --> If the object may have more than one associated object (one track may have many hits, one slice may have many tracks (PFPs), ...)
  art::FindOneP<recob::Track> PfpTrackAssoc(PFP_handle, e, fTrackLabel);          // Assoc between PFP and tracks
  art::FindOneP<recob::Shower> PfpShowerAssoc(PFP_handle, e, fShowerLabel);      // Assoc between PFP and showers
  art::FindManyP<larpandoraobj::PFParticleMetadata> PfpPfpMetaAssoc(PFP_handle, e, fPFParticleMetadataLabel);     // Assoc between PFP annd its metadata
  art::FindOneP<recob::Vertex> PfpVertexAssoc(PFP_handle, e, fVertexLabel);       // Assoc between PFP and vertex

  // I initialize counters of:
  int nTracks=0, nShowers=0;              // Number of tracks and showers
  int nPrimTracks=0, nPrimShowers=0;      // Numer of PRIMARY tracks ans showers

  // I loop over the particles of the slice ( I use the PFP_vec that I will find in the analyze() method and pass it to this method)
  for(const art::Ptr<recob::PFParticle> &PFP : PFP_vec){
    // I set the default variables of the reco vectors
    setDefaultRecoVars();

    //  I check that the ID of the neutrino == ID of the parent of the particle --> PRIMARY PARTICLE
    bool primaryChild = ( prim->Self() == PFP->Parent() );
    // I store the particle information now that I know if it is primary or not
    fReco_part_ID.back() = PFP->Self();                     // ID of the particle
    fReco_part_mother.back() = PFP->Parent();               // Parent of the particle
    fReco_part_isPrimaryChildren.back() = primaryChild;     // If it is primary or not
    fReco_part_PDGCode.back() = PFP->PdgCode();             // PdgCode of the particle --> ONLY SAYS if the trace is electron-like (11) or muon-like (13)

    // Get track_score
    // I get the association between PFPs and PFPMetadata. From here I get the track_score, nu_score, ...
    // Each PFP has a vector of PfpMetadata
    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > PfpMeta_vec = PfpPfpMetaAssoc.at(PFP.key());

    if(!PfpMeta_vec.empty()){   // In case the vector of metadata is NOT empty (otherwise all variables stay with value 1)
      // I extract a map from the PfpMeta_vec which has the important information
      std::map<std::string, float> Pfp_PropMap = PfpMeta_vec.at(0)->GetPropertiesMap();

      // -----TrackScore-----
      if( !(( Pfp_PropMap.find("TrackScore") == Pfp_PropMap.end() ) || ( Pfp_PropMap["TrackScore"]<0 )) ){
        fReco_part_track_score.back() = Pfp_PropMap["TrackScore"];
      }
    }

    // I find again the PfpPfpMeta assoc ONLY if the particle is primary
    if(PFP->IsPrimary()){
      const art::Ptr<recob::PFParticle> PFP_prim = PFP;

      std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > PfpMeta_vec = PfpPfpMetaAssoc.at(PFP_prim.key());
      std::map<std::string, float> Pfp_PropMap = PfpMeta_vec.at(0)->GetPropertiesMap();

      if(!PfpMeta_vec.empty()){
        if( ( std::abs(PFP_prim->PdgCode())==12 ) || ( std::abs(PFP_prim->PdgCode())==14 ) ){
          // -----NuScore----- --> I only find it for the primary PFP
          if(!(( Pfp_PropMap.find("NuScore") == Pfp_PropMap.end() ) || ( Pfp_PropMap["NuScore"]<0 ))){
            fReco_nu_score = Pfp_PropMap["NuScore"];
          }

          // -----IsClearCosmic----- --> I only find it for the primary PFP
          if(  Pfp_PropMap.find("IsClearCosmic") != Pfp_PropMap.end()  ){
            fReco_is_clear_cosmic = Pfp_PropMap["IsClearCosmic"];
          }
        }
      }
    }

    // From the trackScore I can count the number of showers or tracks
    if( ( fReco_part_PDGCode.back()!=12 ) && ( fReco_part_PDGCode.back()!=14 ) ){   // I do not count the neutrino that pandora gives to the slice
      if(fReco_part_track_score.back() > fTrackScoreLimit){
        nTracks++;
        if(primaryChild){
          nPrimTracks++;
        }
      }else{
        nShowers++;
        if(primaryChild){
          nPrimShowers++;
        }
      }
    }
    
    // Get the PFP vertex associated to this ONE PFParticle
    art::Ptr<recob::Vertex> vertex = PfpVertexAssoc.at(PFP.key());

    if(vertex.isNonnull()){
      fReco_part_vertexPos_x.back() = vertex->position().X();
      fReco_part_vertexPos_y.back() = vertex->position().Y();
      fReco_part_vertexPos_z.back() = vertex->position().Z();
    }

    // Get the track information associated to this ONE PFParticle
    const art::Ptr<recob::Track> track = PfpTrackAssoc.at(PFP.key());

    // I use another method to analyse the track
    if(track.isNonnull()){
      fReco_hasTrack.back() = true;
      fReco_isNeutrino.back() = false;
      analyseTrack(e, track, track_handle);
    }else{
      fReco_hasTrack.back() = false;
      if( ( fReco_part_PDGCode.back()==12 ) && ( fReco_part_PDGCode.back()==14 ) ){
        fReco_isNeutrino.back() = true;
      }else{
        fReco_isNeutrino.back() = false;
      }

      if(fVerbose){
        std::cout << std::endl << std::endl << std::endl << "                 PARTICLE WITH NO TRACK" << std::endl << std::endl << std::endl;
      }
      // AAA I CAN ADD SOMETHING HERE SO THAT I DO NOT HAVE A 'DEFAULT' ENTRY IN THE VECTOR FOR THE PFP=NEUTRINO
    }

    // Get the shower information associated to this ONE PFParticle
    const art::Ptr<recob::Shower> shower = PfpShowerAssoc.at(PFP.key());

    // I use another method to analyse the shower
    if(shower.isNonnull()){
      analyseShower(e, shower, shower_handle);
    }
  } // Finish loop over PFParticles

  // I find the MCT object that matches the best this slice
  getSliceTruthMatch(e);

  // I can now save the counting of tracks and showers
  fReco_slice_nTracks = nTracks;
  fReco_slice_nShowers = nShowers;
  fReco_slice_nPrimTracks = nPrimTracks;
  fReco_slice_nPrimShowers = nPrimShowers;

  // For debugging, I print some of the variables:
  if(fVerbose){
    std::cout << std::endl << " ---------- SLICE INFO ----------" << std::endl << std::endl; 
    std::cout << "ID of the slice: " << fReco_slice_ID << std::endl;
    std::cout << "  Number of tracks: " << fReco_slice_nTracks << ". Number of primary tracks " << fReco_slice_nPrimTracks << std::endl;
    std::cout << "  Number of showers: " << fReco_slice_nShowers << ". Number of primary showers " << fReco_slice_nPrimShowers << std::endl;
    std::cout << "  Nu score: " << fReco_nu_score << ". Is it clear cosmic? " << fReco_is_clear_cosmic << std::endl<<std::endl;

    // I loop over the particles 
    std::cout << " ---------- PARTICLE INFO ----------" << std::endl << std::endl; 
    for(long unsigned int i_pfp=0; i_pfp<fReco_part_ID.size(); i_pfp++){
      std::cout << "Id of the particle: " << fReco_part_ID.at(i_pfp) << std::endl;
      std::cout << "  It is primary? " << fReco_part_isPrimaryChildren.at(i_pfp) << std::endl;
      std::cout << "  Id of the mother of the particle: " << fReco_part_mother.at(i_pfp) << std::endl;
      std::cout << "  Pandora PdgCode: " << fReco_part_PDGCode.at(i_pfp) << std::endl;
      std::cout << "  Track score: " << fReco_part_track_score.at(i_pfp) << std::endl;
      std::cout << "  Track length of the particle: " << fReco_part_track_length.at(i_pfp) << std::endl;
      std::cout << "  Energy of the particle: " << fReco_track_kineticEnergy.at(i_pfp) << " . Visible energy: " << fReco_track_visibleEnergy.at(i_pfp) << std::endl;
      std::cout << "  reco direction (x, y, z): (" << fReco_part_startDir_x.at(i_pfp) << ", " << fReco_part_startDir_y.at(i_pfp) << ", " << fReco_part_startDir_z.at(i_pfp) << ")" << std::endl;
      std::cout << "  reco direction (theta, phi): (" << fReco_part_track_theta.at(i_pfp) << ", " << fReco_part_track_phi.at(i_pfp) << ")" << std::endl<<std::endl;

    }

    std::cout << " ---------- TRUTH MATCHING INFO ----------" << std::endl << std::endl; 
    for(long unsigned int i_pfp=0; i_pfp<fReco_track_trueTrackID.size(); i_pfp++){
      std::cout << "g4 Id of the particle: " << fReco_track_trueTrackID.at(i_pfp) << std::endl;
      std::cout << "  Does it have a track " << fReco_hasTrack.at(i_pfp) << std::endl;
      std::cout << "  True PdgCode: " << fReco_track_truePDG.at(i_pfp) << std::endl;
      std::cout << "  True energy: " << fReco_true_E0.at(i_pfp) << std::endl;
      std::cout << "  True end process: " << fReco_true_endProcess.at(i_pfp) << std::endl;
      std::cout << "  True start position: (" << fReco_true_StartPos_x.at(i_pfp) << ", " << fReco_true_StartPos_y.at(i_pfp) << ", " << fReco_true_StartPos_z.at(i_pfp) << ")" << std::endl<<std::endl;
    }
  }
}

/*
analyseTrack
RECEIVES:
- The art::Event 'e'
- The recob::Track variable 'track'
- The track handle to to associations to the track (calorimetry, for example)
*/
void sbnd::AnalyzeReco::analyseTrack(const art::Event &e, 
                                    const art::Ptr<recob::Track> &track, 
                                    const art::ValidHandle<std::vector<recob::Track>> &track_handle)
                                    {

  // I extract variables directly from the track
  fReco_part_track_length.back() = track->Length();
  fReco_part_track_theta.back() = track->Theta();
  fReco_part_track_phi.back() = track->Phi();

  fReco_part_startDir_x.back() = track->StartDirection().X();
  fReco_part_startDir_y.back() = track->StartDirection().Y();
  fReco_part_startDir_z.back() = track->StartDirection().Z();

  fReco_part_track_startPos_x.back() = track->Start().X();
  fReco_part_track_startPos_y.back() = track->Start().Y();
  fReco_part_track_startPos_z.back() = track->Start().Z();

  fReco_part_track_EndPos_x.back() = track->End().X();
  fReco_part_track_EndPos_y.back() = track->End().Y();
  fReco_part_track_EndPos_z.back() = track->End().Z();
  
  // I get the Assoc that I need:
  art::FindManyP<anab::Calorimetry> TrackCaloAssoc(track_handle, e, fCalorimetryLabel);     // calorimetry of each track
  art::FindManyP<recob::Hit> TrackHitAssoc(track_handle, e, fTrackLabel);                     // Hits of each track
  art::FindManyP<anab::ParticleID> TrackChi2Assoc(track_handle, e, fChi2ModuleLabel);       // Values of chi2 of each track

  // Get the vector of calorimetry associated to the track
  const std::vector<art::Ptr<anab::Calorimetry>> calos = TrackCaloAssoc.at(track.key());
  // Find the best plane
  // Look for the plane with the biggest number of hits: (0, 1, 2), (U, V, Y)
  /*
  If the number of planes is not 3(what it should be) --> set maxHits to -1
  Compare the number of hits of each plane and save the number of hits of the plane ith the biggest number of hits
  Find the plane with the biggest number of hits
  */
  long unsigned int maxHits = 0;
  int bestPlane = -1;
  if(calos.size()==3){    // There are 3 planes where calorimetry is measured
    maxHits = std::max({calos[0]->dEdx().size(), calos[1]->dEdx().size(), calos[2]->dEdx().size()});
    if(calos[0]->dEdx().size() == maxHits){
      bestPlane=0;
    }else if(calos[1]->dEdx().size() == maxHits){
      bestPlane=1;
    }else if(calos[2]->dEdx().size() == maxHits){
      bestPlane=2;
    }
  // Get the calorimetry from the best plane
    getCalo(calos[bestPlane]);
  }

  // Get true particle of track: TRUTH MATCHING
  const std::vector<art::Ptr<recob::Hit>> trackHits = TrackHitAssoc.at(track.key());
  if (!trackHits.empty()){
    getTrackTruthMatch(e, trackHits);
  }

  // Get the chi2 values
  const std::vector<art::Ptr<anab::ParticleID>> chi2s = TrackChi2Assoc.at(track.key());
  if(chi2s.size() == 3){
    getChi2PID(chi2s[bestPlane]); 
  }

}


/*
analyseShower
RECEIVES:
- The art::Event 'e'
- The recob::Track variable 'shower'
- The track handle to to associations to the track (calorimetry, for example). IT IS NOT USED
*/
void sbnd::AnalyzeReco::analyseShower(const art::Event &e, 
                                    const art::Ptr<recob::Shower> &shower, 
                                    const art::ValidHandle<std::vector<recob::Shower>> &shower_handle){

  // I store some of the variables:
  fReco_shower_dir_x.back() = shower->Direction().X();
  fReco_shower_dir_y.back() = shower->Direction().Y();
  fReco_shower_dir_z.back() = shower->Direction().Z();

  fReco_shower_start_x.back() = shower->ShowerStart().X();
  fReco_shower_start_y.back() = shower->ShowerStart().Y();
  fReco_shower_start_z.back() = shower->ShowerStart().Z();

  // I find the opening angle and length of the shower
  if (shower->has_open_angle()){
    fReco_shower_openAngle.back() = shower->OpenAngle();
  }

  if (shower->has_length()){
    fReco_shower_length.back() = shower->Length();
  }

  // I use what pandora says is the best plane to get the simple calorimetry:
  int bestPlane = shower->best_plane();
  //Get shower calorimetry at bestPlane
  fReco_shower_energy.back() = shower->Energy()[bestPlane];
}


/*
getCalo
RECEIVES:
- The CALORIMETRY pointer
Saves calorimetry information from the associated track
*/
void sbnd::AnalyzeReco::getCalo(const art::Ptr<anab::Calorimetry> &calo)
{
  // I extract the info that I want
  fReco_track_kineticEnergy.back() = calo->KineticEnergy();
  fReco_track_dEdx.back() = calo->dEdx();
  fReco_track_residualRange.back() = calo->ResidualRange();
  fReco_track_range.back() = calo->Range();

  // I also get the 'visible' energy, This should be the same as the kinetic energy
  double visible_energy = 0;

  for(size_t i = 0; i < calo->dEdx().size(); i++) {
    if (calo->dEdx()[i] <= 1000.){
      visible_energy     += calo->dEdx()[i] * calo->TrkPitchVec()[i];
    }
  }
  fReco_track_visibleEnergy.back() = visible_energy;
  fReco_track_pitch.back() = calo->TrkPitchVec();

}


/*
getTracktruthMatch
RECEIVES:
- The art::Event 'e'
- The vector of hits, that will be used for the truth matching
FINDS:
the truth information of the tracks (based on their hits) of pandora
*/

void sbnd::AnalyzeReco::getTrackTruthMatch(const art::Event &e, 
                                          const std::vector<art::Ptr<recob::Hit>> hits)
{

  // I get the detector timing info (the time of its clock)
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  // Finds the particle that is responsible for the track of this PFP:
  /*
  - Looks at all the hits of this track
  - Sees which truth particle deposited energy is these same points as the reconstrcted particle
  - Returns the g4 ID
  */
  const int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, false);
  fReco_track_trueTrackID.back() = trackID;

  // Calculate the completeness and the purity of the trace comparng the reco hits and the truth hits
  // I only find it for particles that leave a non-0 number of hits
  if(fHitsMap[trackID] != 0){
    fReco_track_completeness.back() = completeness(e, hits, trackID);
  }
  
  fReco_track_purity.back() = purity(e, hits, trackID);
  
  // I save the important g4 information from the track
  for (size_t i_p = 0; i_p < fG4_part_trackID.size(); i_p++){
    if (fG4_part_trackID.at(i_p) == std::abs(trackID)){
      fReco_track_truePDG.back() = fG4_part_PDGcode.at(i_p);

      // Energy
      fReco_true_E0.back() = fG4_part_E0.at(i_p);
      fReco_true_Ef.back() = fG4_part_Ef.at(i_p);

      // Momentum
      fReco_true_p0_x.back() = fG4_part_P0_x.at(i_p);
      fReco_true_p0_y.back() = fG4_part_P0_y.at(i_p);
      fReco_true_p0_z.back() = fG4_part_P0_z.at(i_p);
      fReco_true_pf_x.back() = fG4_part_Pf_x.at(i_p);
      fReco_true_pf_y.back() = fG4_part_Pf_y.at(i_p);
      fReco_true_pf_z.back() = fG4_part_Pf_z.at(i_p);

      // End process
      fReco_true_endProcess.back() = fG4_part_endProcess.at(i_p);

      // Start position
      fReco_true_StartPos_x.back() = fG4_part_StartPos_x.at(i_p);
      fReco_true_StartPos_y.back() = fG4_part_StartPos_y.at(i_p);
      fReco_true_StartPos_z.back() = fG4_part_StartPos_z.at(i_p);
    }
  }
}

/*
getTracktruthMatch
RECEIVES:
- The art::Event 'e'
(- the map fHitsMap: g4-trackID --> number of hits of this track)
(- the map fMCTruthHitsMap: MCTruth object --> number of hitsb of this whole object)
FINDS:
the MCT that corresponds the best to this slice and analyzes it
*/

void sbnd::AnalyzeReco::getSliceTruthMatch(const art::Event &e)
{

  // I save the generator information only for the MCTruth object that matches the slice
  // Do the truth matching of the slice (to know if it is cosmic or other things)

  // I set minimum values for variables I will change later  
  int max_hits = def_int;
  int best_track_ID = def_int;

  // I look for the 'best_track_ID'--> The track id of the track with the most hits in this slice
  for(auto const& [track_ID, nhits] : fHitsMap) {
    if(nhits > max_hits){
      best_track_ID = track_ID;
      max_hits=nhits;
    }
  }

  if(fVerbose){
    std::cout << "Best track ID: " << best_track_ID << std::endl;
  }

  // Now I find the best MCTruth object: the one with the most hits in this slice
  max_hits=def_int;
  art::Ptr<simb::MCTruth> best_MCT = art::Ptr<simb::MCTruth>();
  for(auto const& [MCT, nhits] : fMCTruthHitsMap){
    if(nhits > max_hits) {
      max_hits = nhits;
      best_MCT  = MCT;
    }
  }

  // With this information, I can now save the generator information of ONLY the ONE MCTruth object that matches the best this slice
  if(!best_MCT.isNull()){
    analyseTruth_gen(best_MCT);
  }
}

/*
getPrimaryPFP
RECEIVES:
- The vector of PFParticles
RETURNS:
- ONE PFParticle if it is primary (the neutrino/primary cosmic ray)
- Empty object if there is not any primary PFP found
*/
art::Ptr<recob::PFParticle> sbnd::AnalyzeReco::getPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &PFP_vec)
{
  // Finds and return the primary particle
  for(auto PFP : PFP_vec){
    if(PFP->IsPrimary()){
      return PFP;
    }
  }
    
  // If none found, returns an empty (nullptr) object --> NO primary particle was found
  return art::Ptr<recob::PFParticle>();
}

//Change the name of the file saved on the tree when a new file is opened
void sbnd::AnalyzeReco::respondToOpenInputFile(art::FileBlock const& fb) { 
  fFilename = fb.fileName();
}

void sbnd::AnalyzeReco::beginSubRun(const art::SubRun &sr)
{
  clearSubRunVars();

  // Get POT
  art::Handle<sumdata::POTSummary> pot_handle;
  sr.getByLabel(fPOTModuleLabel, pot_handle);
  if(!pot_handle.isValid()){
    if(fVerbose){
      std::cout << "POT product " << fPOTModuleLabel << " not found..." << std::endl;
    }
    fPOT = 0;
    fSpill = 1;
  }else {
    fPOT = pot_handle->totpot;
    fSpill = pot_handle->totspills;
  } 
}


void sbnd::AnalyzeReco::endSubRun(const art::SubRun &sr)
{
  fSubRunTree->Fill();
}


void sbnd::AnalyzeReco::beginJob()                                               // Implementation of the method beginJob
{
  // Get the TFileService to create the output TTree
  art::ServiceHandle<art::TFileService> tfs;

  fSubRunTree = tfs->make<TTree>("subrun_tree", "SubRun TTree");
  // Branches of the subrun TTree
  fSubRunTree->Branch("POT", &fPOT);
  fSubRunTree->Branch("spill", &fSpill);
  fSubRunTree->Branch("num_gen_evts", &fNum_gen_evts);
  fSubRunTree->Branch("subrun_ID", &fSubrun_ID);

  fTree = tfs->make<TTree>("tree", "Output TTree");
  // Branches of the output TTree
  fTree->Branch("filename", &fFilename);
  fTree->Branch("event_ID", &fEvent_ID);
  fTree->Branch("run_ID", &fRun_ID);
  fTree->Branch("gen_index", &fGen_index);

  // MCTruth variables
  fTree->Branch("nu_PDG", &fNu_PDG);
  fTree->Branch("nu_E0", &fNu_E0);
  fTree->Branch("nu_weight", &fNu_weight); 
  fTree->Branch("nu_interaction_mode", &fNu_interaction_mode);
  fTree->Branch("nu_interaction_type", &fNu_interaction_type);
  fTree->Branch("nu_CC_NC", &fNu_CC_NC);
  fTree->Branch("nu_target", &fNu_target);
  fTree->Branch("nu_HitNuc", &fNu_HitNuc);
  fTree->Branch("nu_HitQuark", &fNu_HitQuark);
  fTree->Branch("nu_W", &fNu_W);
  fTree->Branch("nu_X", &fNu_X);
  fTree->Branch("nu_Y", &fNu_Y);
  fTree->Branch("nu_Q2", &fNu_Q2);  

  //MCTruth Particles variables
  if (fSave_truth_gen_particles) {
    fTree->Branch("gen_part_trackID", &fGen_part_trackID);
    fTree->Branch("gen_part_statusCode", &fGen_part_statusCode);
    fTree->Branch("gen_part_mother", &fGen_part_mother);
    fTree->Branch("gen_part_PDGcode", &fGen_part_PDGcode);
    fTree->Branch("gen_part_mass", &fGen_part_mass);
    fTree->Branch("gen_part_E0", &fGen_part_E0);
    fTree->Branch("gen_part_StartPos_x", &fGen_part_StartPos_x);
    fTree->Branch("gen_part_StartPos_y", &fGen_part_StartPos_y);
    fTree->Branch("gen_part_StartPos_z", &fGen_part_StartPos_z);
    fTree->Branch("gen_part_P0_x", &fGen_part_P0_x);
    fTree->Branch("gen_part_P0_y", &fGen_part_P0_y);
    fTree->Branch("gen_part_P0_z", &fGen_part_P0_z);
  }

  //Geant4 Particles variables
  if(fSave_truth_g4){
    fTree->Branch("g4_part_trackID", &fG4_part_trackID);
    fTree->Branch("g4_part_mother", &fG4_part_mother);
    fTree->Branch("g4_part_PDGcode", &fG4_part_PDGcode);
    fTree->Branch("g4_part_mass", &fG4_part_mass);
    fTree->Branch("g4_part_E0", &fG4_part_E0);
    fTree->Branch("g4_part_Ef", &fG4_part_Ef);
      
    fTree->Branch("g4_part_startPos_x", &fG4_part_StartPos_x);
    fTree->Branch("g4_part_startPos_y", &fG4_part_StartPos_y);
    fTree->Branch("g4_part_startPos_z", &fG4_part_StartPos_z);
    fTree->Branch("g4_part_EndPos_x", &fG4_part_EndPos_x);
    fTree->Branch("g4_part_EndPos_y", &fG4_part_EndPos_y);
    fTree->Branch("g4_part_EndPos_z", &fG4_part_EndPos_z);

    fTree->Branch("g4_part_P0_x", &fG4_part_P0_x);
    fTree->Branch("g4_part_P0_y", &fG4_part_P0_y);
    fTree->Branch("g4_part_P0_z", &fG4_part_P0_z);
    fTree->Branch("g4_part_Pf_x", &fG4_part_Pf_x);
    fTree->Branch("g4_part_Pf_y", &fG4_part_Pf_y);
    fTree->Branch("g4_part_Pf_z", &fG4_part_Pf_z);
    
    fTree->Branch("g4_part_process", &fG4_part_process);
    fTree->Branch("g4_part_endProcess", &fG4_part_endProcess);
  }

  //Pandora Particles variables

  if(fSave_reco){
    fTree->Branch("reco_event_nSlices", &fReco_event_nSlices);

    fTree->Branch("reco_slice_ID", &fReco_slice_ID);
    fTree->Branch("reco_slice_nPFParticles", &fReco_slice_nPFParticles);
    fTree->Branch("reco_slice_prim_ID", &fReco_slice_prim_ID);
    fTree->Branch("reco_slice_prim_PDG", &fReco_slice_prim_PDG);
    fTree->Branch("reco_slice_prim_nDaughters", &fReco_slice_prim_nDaughters);

    fTree->Branch("reco_slice_nTracks", &fReco_slice_nTracks);
    fTree->Branch("reco_slice_nShowers", &fReco_slice_nShowers);
    fTree->Branch("reco_slice_nPrimTracks", &fReco_slice_nPrimTracks);
    fTree->Branch("reco_slice_nPrimShowers", &fReco_slice_nPrimShowers);

    fTree->Branch("reco_slice_vertexPos_x", &fReco_slice_vertexPos_x);
    fTree->Branch("reco_slice_vertexPos_y", &fReco_slice_vertexPos_y);
    fTree->Branch("reco_slice_vertexPos_z", &fReco_slice_vertexPos_z);

    fTree->Branch("reco_hasTrack", &fReco_hasTrack);
    fTree->Branch("reco_isNeutrino", &fReco_isNeutrino);
    fTree->Branch("reco_part_ID", &fReco_part_ID);
    fTree->Branch("reco_part_mother", &fReco_part_mother);
    fTree->Branch("reco_part_PDGCode", &fReco_part_PDGCode);
    fTree->Branch("reco_part_isPrimaryChildren", &fReco_part_isPrimaryChildren);

    fTree->Branch("reco_part_vertexPosX", &fReco_part_vertexPos_x);
    fTree->Branch("reco_part_vertexPosY", &fReco_part_vertexPos_y);
    fTree->Branch("reco_part_vertexPosZ", &fReco_part_vertexPos_z);
    
    fTree->Branch("reco_part_track_score", &fReco_part_track_score);
    fTree->Branch("reco_nu_score", &fReco_nu_score);
    fTree->Branch("reco_is_clear_cosmic", &fReco_is_clear_cosmic);
    fTree->Branch("reco_part_track_length", &fReco_part_track_length);
    fTree->Branch("reco_part_track_theta", &fReco_part_track_theta);
    fTree->Branch("reco_part_track_phi", &fReco_part_track_phi);

    fTree->Branch("reco_part_track_startPosX", &fReco_part_track_startPos_x);
    fTree->Branch("reco_part_track_startPosY", &fReco_part_track_startPos_y);
    fTree->Branch("reco_part_track_startPosZ", &fReco_part_track_startPos_z);

    fTree->Branch("reco_part_track_EndPosX", &fReco_part_track_EndPos_x);
    fTree->Branch("reco_part_track_EndPosY", &fReco_part_track_EndPos_y);
    fTree->Branch("reco_part_track_EndPosZ", &fReco_part_track_EndPos_z);

    fTree->Branch("reco_part_startDirX", &fReco_part_startDir_x);
    fTree->Branch("reco_part_startDirY", &fReco_part_startDir_y);
    fTree->Branch("reco_part_startDirZ", &fReco_part_startDir_z);

    fTree->Branch("reco_track_kineticEnergy", &fReco_track_kineticEnergy);
    fTree->Branch("reco_track_visibleEnergy", &fReco_track_visibleEnergy);
    fTree->Branch("reco_track_dEdx", &fReco_track_dEdx);
    fTree->Branch("reco_track_residualRange", &fReco_track_residualRange);
    fTree->Branch("reco_track_range", &fReco_track_range);

    fTree->Branch("reco_track_completeness", &fReco_track_completeness);
    fTree->Branch("reco_track_purity", &fReco_track_purity); 

    fTree->Branch("reco_part_chi2_muon", &fReco_part_chi2_muon);
    fTree->Branch("reco_part_chi2_pion", &fReco_part_chi2_pion);
    fTree->Branch("reco_part_chi2_kaon", &fReco_part_chi2_kaon);
    fTree->Branch("reco_part_chi2_proton", &fReco_part_chi2_proton);

    fTree->Branch("reco_shower_dirX", &fReco_shower_dir_x);
    fTree->Branch("reco_shower_dirY", &fReco_shower_dir_y);
    fTree->Branch("reco_shower_dirZ", &fReco_shower_dir_z);

    fTree->Branch("reco_shower_startX", &fReco_shower_start_x);
    fTree->Branch("reco_shower_startY", &fReco_shower_start_y);
    fTree->Branch("reco_shower_startZ", &fReco_shower_start_z);

    fTree->Branch("reco_shower_openAngle", &fReco_shower_openAngle);
    fTree->Branch("reco_shower_length", &fReco_shower_length);

    fTree->Branch("reco_shower_energy", &fReco_shower_energy);
  }

  // Truth Matching variables
  if(fSave_TruthMatching){
    fTree->Branch("reco_track_trueTrackID", &fReco_track_trueTrackID);
    fTree->Branch("reco_track_truePDG", &fReco_track_truePDG);
    fTree->Branch("reco_true_E0", &fReco_true_E0);
    fTree->Branch("reco_true_Ef", &fReco_true_Ef);
    fTree->Branch("reco_true_p0_x", &fReco_true_p0_x);
    fTree->Branch("reco_true_p0_y", &fReco_true_p0_y);
    fTree->Branch("reco_true_p0_z", &fReco_true_p0_z);
    fTree->Branch("reco_true_pf_x", &fReco_true_pf_x);
    fTree->Branch("reco_true_pf_y", &fReco_true_pf_y);
    fTree->Branch("reco_true_pf_z", &fReco_true_pf_z);
    fTree->Branch("reco_true_endProcess", &fReco_true_endProcess);
    fTree->Branch("reco_true_startPos_x", &fReco_true_StartPos_x);
    fTree->Branch("reco_true_startPos_y", &fReco_true_StartPos_y);
    fTree->Branch("reco_true_startPos_z", &fReco_true_StartPos_z);
  }

}

void sbnd::AnalyzeReco::endJob()                                                 // Implementation of the method beginJob
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::AnalyzeReco)