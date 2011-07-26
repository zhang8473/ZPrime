#define IsMC
#ifdef IsMC
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

//Simulation Infomation
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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
//Tracking Tools
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

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

//HLT
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

//ROOT
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

#define maxFilterObjects 30
#define MaxHLTObjDeviation2 0.1 // the maximum deviation between pt and eta of reco objs and HLT objs to be accepted as a HLT objs. this value is the veto cone of muon squared
using namespace std;

class MuESelector : public edm::EDFilter {
   public:
      explicit MuESelector(const edm::ParameterSet&);
      ~MuESelector();
   private:
      virtual void beginJob();
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual reco::Muon::ArbitrationType MuonArbitrationTypeFromString( const std::string &);
      inline  void ClearVecs_RECO();
#ifdef TrackingParticles
      virtual void GetDecayChains(TrackingParticleRef tpr,vector<int> *DChains, vector <TheMuonType> &type, HepMC::GenEvent &HepGenEvent);
      virtual void HepMCParentTree(HepMC::GenParticle *);
      virtual void SimTrackDaughtersTree(SimTrack *);
      inline  bool IstheSameDChain(vector<int> &,vector<int> &);
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
   //TTree&TFile
   TFile *file;
   string FileName;
   TTree *Muons_Tree,*HLTNames_Tree,*ErrorMsg_Tree;
   //Event Information
   struct General_EvtInfo
   {
     ULong64_t RUN,EVENT,LS,ORBIT,BX;
     ULong64_t First_Run,First_Event;//First_Run # and First_Event # in this HLT scope
   } Info;
   Bool_t isRealData;
   Double_t GenEventWeight;
   //PU information
   Int_t numberOfPUVertices;
   //Muon
   Float_t MuonPtCut,DiMuonInvarMassCut;
   vector<Float_t> *pt,*eta,*phi,*Vertex_x,*Vertex_y,*Vertex_z,*isoR03sumPt,*isoR03emEt,*isoR03hadEt,*isoR03hoEt,*isoR03nJets,*isoR03nTracks,*isoR05sumPt,*isoR05emEt,*isoR05hadEt,*isoR05hoEt,*isoR05nJets,*isoR05nTracks,*isoemVetoEt,*isohadVetoEt,*isohoVetoEt,*normalizedChi2,*TrkRelChi2,*CaloE_emMax,*CaloE_emS9,*CaloE_emS25,*CaloE_hadMax,*CaloE_hadS9,*Calo_emPos_R,*Calo_emPos_eta,*Calo_emPos_phi,*Calo_hadPos_R,*Calo_hadPos_eta,*Calo_hadPos_phi,*dEdx,*dEdxError,*TrkKink,*GlbKink;
   vector<Bool_t> *chargeMinus,*isGlobalMu,*isTrackerMu;
   vector<Int_t> *dEdx_numberOfSaturatedMeasurements,*dEdx_numberOfMeasurements,*numberOfMatchedSegments,*numberOfMatchedStations;
   //Dimuon Mass saved position: num_mu*mu1-mu1*(mu1+1)/2+(mu2-mu1)-1 (mu2>mu1) (mu1 starts at 0)
   vector<Float_t> *DiMuonInvariantMass,*CosThetaStar,*AngleBetweenDiMuon;
   //Tracks
   edm::InputTag tracksTag;
   string dEdxTag;
   vector<UInt_t> *InnerTrack_nValidTrackerHits,*InnerTrack_nValidPixelHits,*InnerTrack_nLostTrackerHits,*InnerTrack_nLostPixelHits,*InnerTrack_ndof,*GlobalTrack_ndof,*numberOfValidMuonHits;
   vector<Float_t> *InnerTrack_chi2,*GlobalTrack_chi2,*DXYwtBS,*DXYwtPV,*DXYErrwtPV;
   //Muon Selectors
   unsigned int num_Cuts;
   reco::Muon::ArbitrationType MuonArbitrationType;
   vector<muon::SelectionType> OfficialMuonSelectors;
   vector<Bool_t> * SelectorPassed[30],*MySelector;
   //Muon Chamber and Segment Match
   vector<Float_t> *TrackDistToChamberEdge,*TrackDistToChamberEdgeErr,*XTrack,*YTrack,*XErrTrack,*YErrTrack,*XSegment,*YSegment,*XErrSegment,*YErrSegment;
   vector<Bool_t> *IsCSCChamber,*IsSegmentOwnedExclusively,*IsSegmentBelongsToTrackByDR,*IsSegmentBelongsToTrackByCleaning;
   vector<UInt_t> *StationMask,*RequiredStationMask;
   vector<Byte_t> *SectorChamber,*MuonIndex,*NumberOfHitsInSegment;
   vector<Char_t> *StationRing;
   double maxChamberDist,maxChamberDistPull;
   //PrimaryVertex
   string PrimaryVerticesTag;
   vector<Float_t> *vx,*vxError,*vy,*vyError,*vz,*vzError;
   //HLT
   vector<string> *HLTNamesSet;
   vector<Bool_t> *HLTacceptance;
   edm::InputTag TriggerResultsTag;

   //HLTObjects
   UInt_t num_HLTsSaveObjs;
   UInt_t HLTFilterNamesAcceptenceIndex[maxFilterObjects];
   vector<Float_t> *HLTObj_pt[maxFilterObjects],*HLTObj_eta[maxFilterObjects],*HLTObj_phi[maxFilterObjects];
   vector<Bool_t> *isHLTObj[maxFilterObjects];
   vector<string> HLTFilter_HLTNames;
   vector<edm::InputTag> HLTFilterNames;
   edm::InputTag TriggerEventTag;
#ifdef IsMC
   vector<HepMC::GenParticle *> ParToHep;
#endif

#ifdef TrackingParticles // combination of simulated tracks
   UInt_t minTrackHits;
   vector<Float_t> *TrkParticles_pt,*TrkParticles_eta,*TrkParticles_phi;
   vector<Int_t> *TrkParticles_pdgId,*TrkParticles_charge;
   vector<Double_t> *SharedHitsRatio,*MCMatchChi2;
   vector<Int_t> *DChains,*theSameWithMuon;
   vector<Bool_t> *IsParInHep;
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

   Float_t HepMCFourVec[4];
   vector<Float_t> *Gen_pt,*Gen_eta,*Gen_phi,*Gen_vx,*Gen_vy,*Gen_vz,*Gen_vt;
   vector<Int_t> *Gen_pdgId;
   string HepMCTag;
   inline void ClearVecs_HepMC();

   //Errors and Warnings reports
   struct ErrorMsg
   {
     Byte_t ErrorCode;
     ULong64_t Run,Event;
   }Error;
   inline Bool_t IsTheSameSegment(UInt_t Seg1,UInt_t Seg2) {
     if ( (*IsCSCChamber)[Seg1]!=(*IsCSCChamber)[Seg2] ) return false;
     if ( (*XSegment)[Seg1]!=(*XSegment)[Seg2] || (*YSegment)[Seg1]!=(*YSegment)[Seg2] ) return false;
     if ( (*StationRing)[Seg1] != (*StationRing)[Seg2] || (*SectorChamber)[Seg1] != (*SectorChamber)[Seg2] ) return false;
     return true;
   }
};

// ------------ method called once each job just before starting event loop  ------------
void MuESelector::beginJob(){
}

// ------------ method called once each job just after ending the event loop  ------------
void MuESelector::endJob() {
}
