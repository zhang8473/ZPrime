#define HepMCGenParticle
#ifdef HepMCGenParticle
//#define TrackingParticles
#endif


// system include files 
#include <memory>
#include <vector>
#include <map>
#include <cmath>
#include <cstring>
#include <algorithm>

// Framework Core
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Muon
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//Tracks
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

//Tracking Particles, Match to simulated tracks
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimMuon/MCTruth/interface/MuonAssociatorByHits.h"

//SimTracks
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
//SimVertex
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

//Generation level particles
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//pdgTable
//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

//HLT
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

/*
//CSC rechit
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

//Geometry
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
*/

//ROOT
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

#define maxFilterObjects 30
#define MaxHLTObjDeviation2 1e-3 // the maximum deviation between pt and eta of reco objs and HLT objs to be accepted as a HLT objs. this value is the veto cone of muon squared
using namespace std;

class DiMuonSelector : public edm::EDFilter {
   public:
      explicit DiMuonSelector(const edm::ParameterSet&);
      ~DiMuonSelector();
   private:
      virtual void beginJob();
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual reco::Muon::ArbitrationType MuonArbitrationTypeFromString( const std::string &);
      inline void ClearVecs_RECO();
#ifdef TrackingParticles
      virtual void HepMCParentTree(HepMC::GenParticle *);
      virtual void SimTrackDaughtersTree(SimTrack *);
      inline bool IstheSameDChain(vector<int> &,vector<int> &);
      inline Int_t FindSimTrackRecordingPosition( Int_t ParToSimPos ) {
	Int_t count=0;
	vector<Bool_t>::const_iterator IsParInHep_iter = IsParInHep->begin();
	for (; IsParInHep_iter != IsParInHep->end(); IsParInHep_iter++ ) {
	  if (!*IsParInHep_iter) {
	    if (count==ParToSimPos) break;
	    else count++;
	  }
	}
	return IsParInHep_iter-IsParInHep->begin();
      }
      inline Int_t FindHepMCRecordingPosition( Int_t ParToHepPos ) {
	Int_t count=0;
	vector<Bool_t>::const_iterator IsParInHep_iter = IsParInHep->begin();
	for (; IsParInHep_iter != IsParInHep->end(); IsParInHep_iter++ )
	  if (*IsParInHep_iter) {
	    if (count==ParToHepPos) break;
	    else count++;
	  }
	return IsParInHep_iter-IsParInHep->begin();
      }
#endif
      virtual void endJob() ;

//-------------member data------------
   Bool_t FirstEntry;
   //TTree&TFile
   TFile *file;
   string FileName;
   TTree *Muons_Tree,*Summarization_Tree,*ErrorMsg_Tree;
   struct General_EvtInfo
   {
     ULong64_t RUN,EVENT,LS,ORBIT,BX;
     ULong64_t First_Run,First_Event;//First_Run # and First_Event # in this HLT scope
#ifdef HepMCGenParticle
     Byte_t num_ZpMuMuInMC,num_ZMuMuInMC;
#endif
     Bool_t isRealData;
   } Info;
   //Muon
   Float_t MuonPtCut,DiMuonInvarMassCut;
   vector<Float_t> *pt,*eta,*phi,*Vertex_x,*Vertex_y,*Vertex_z,*isoR03sumPt,*isoR03emEt,*isoR03hadEt,*isoR03hoEt,*isoR03nJets,*isoR03nTracks,*isoR05sumPt,*isoR05emEt,*isoR05hadEt,*isoR05hoEt,*isoR05nJets,*isoR05nTracks,*isoemVetoEt,*isohadVetoEt,*isohoVetoEt,*normalizedChi2,*TrkRelChi2,*CaloE_emMax,*CaloE_emS9,*CaloE_emS25,*CaloE_hadMax,*CaloE_hadS9,*Calo_emPos_R,*Calo_emPos_eta,*Calo_emPos_phi,*Calo_hadPos_R,*Calo_hadPos_eta,*Calo_hadPos_phi,*dEdx,*dEdxError,*TrkKink,*GlbKink;
   vector<Bool_t> *chargeMinus,*isGlobalMu,*isTrackerMu;
   vector<Int_t> *dEdx_numberOfSaturatedMeasurements,*dEdx_numberOfMeasurements;
   //Dimuon Mass saved position: num_mu*mu1-mu1*(mu1+1)/2+(mu2-mu1)-1 (mu2>mu1) (mu1 starts at 0)
   vector<Float_t> *DiMuonInvariantMass,*CosThetaStar;
   //Tracks
   edm::InputTag tracksTag;
   string dEdxTag;
   vector<UInt_t> *InnerTrack_nValidTrackerHits,*InnerTrack_nValidPixelHits,*InnerTrack_nLostTrackerHits,*InnerTrack_nLostPixelHits,*InnerTrack_ndof,*GlobalTrack_ndof;
   vector<Float_t> *InnerTrack_chi2,*GlobalTrack_chi2,;
   //Muon Selectors
   unsigned int num_Cuts;
   reco::Muon::ArbitrationType MuonArbitrationType;
   vector<muon::SelectionType> OfficialMuonSelectors;
   vector<Bool_t> * SelectorPassed[30],*MySelector;
   //Muon Chamber Match  0-3 DT, 4-7 CSC
   vector<Float_t> *TrackDistToChamberEdge,*TrackDistToChamberEdgeErr,*DXTrackToSegment,*DYTrackToSegment,*DXErrTrackToSegment,*DYErrTrackToSegment,*DRTrackToSegment,*DRErrTrackToSegment;
   vector<Bool_t> *IsSegmentBelongsToTrackByDR,*IsSegmentBelongsToTrackByCleaning;
   vector<UInt_t> *NumberOfHitsInSegment,*StationMask,*RequiredStationMask;
   double maxChamberDist,maxChamberDistPull;
   //PrimaryVertex
   string PrimaryVerticesTag;
   vector<Float_t> *vx,*vxError,*vy,*vyError,*vz,*vzError;
   //HLT
   vector<string> *HLTNamesSet;
   vector<Bool_t> *HLTacceptance;
   edm::InputTag TriggerResultsTag;
   //HLTObjects
   unsigned int num_HLTsSaveObjs;
   vector<Float_t> *HLTObj_pt[maxFilterObjects],*HLTObj_eta[maxFilterObjects],*HLTObj_phi[maxFilterObjects];
   vector<Bool_t> *isHLTObj[maxFilterObjects];
   vector<edm::InputTag> HLTFilterNames;
   edm::InputTag TriggerEventTag;
#ifdef TrackingParticles // combination of simulated tracks
   UInt_t minTrackHits;
   vector<Float_t> *TrkParticles_pt,*TrkParticles_eta,*TrkParticles_phi;
   vector<Int_t> *TrkParticles_pdgId,*TrkParticles_charge;
   vector<Double_t> *SharedHitsRatio,*MCMatchChi2;
   vector<Int_t> *DChains,*theSameWithMuon;
   //temparory variables
   bool ChainRecord;
   vector< vector<Int_t> > SimChains;
   vector< vector<Int_t> > HepMCChains;
   vector<Int_t> DChain;
   vector<SimTrack *> ParToSim,MaskOut;
   map< SimTrack *, vector<SimTrack *> > Daughters;
   vector<SimVertex> SVC;
   //Others
   ULong64_t NumMisMatch;
#endif

#ifdef HepMCGenParticle
   Float_t HepMCFourVec[4];
   vector<Float_t> *Gen_pt,*Gen_eta,*Gen_phi,*Gen_vx,*Gen_vy,*Gen_vz,*Gen_vt;
   vector<Int_t> *Gen_pdgId;
   vector<Bool_t> *IsParInHep;
   vector<HepMC::GenParticle *> ParToHep;
   string HepMCTag;
   inline void ClearVecs_HepMC();
#endif
   //Summarization
   struct Summary
   {
     ULong64_t Total_Events,Total_TrackerMuons,Total_GlobalMuon,Total_GlobalnotTrackerMuon;
     ULong64_t First_Run,First_Event;
     Double_t CrossSection;
   }Summarization;
   //Errors and Warnings reports
   struct ErrorMsg
   {
     Byte_t ErrorCode;
     ULong64_t Run,Event;
   }Error;
   inline void ClearSummarization() {Summarization.Total_Events=0;Summarization.Total_TrackerMuons=0;Summarization.Total_GlobalMuon=0;Summarization.Total_GlobalnotTrackerMuon=0;Summarization.First_Run=Info.First_Run;Summarization.First_Event=Info.First_Event;};
};

// ------------ method called once each job just before starting event loop  ------------
void DiMuonSelector::beginJob(){
}

// ------------ method called once each job just after ending the event loop  ------------
void DiMuonSelector::endJob() {
}
