// -*- C++ -*-
//  Total_Events=0;HLTNameSaved = false;

// Package:    Analyzer
// Class:      Analyzer
// 
/**\class Analyzer Analyzer.cc Analyzer/Analyzer/src/Analyzer.cc

 Description: <one line class summary>
 It will save the information about primary vertex, HLT, muons with the original cut (GlobalMuon&&pt>20Gev/c). The muon objects will be saved in descending order of pt.
 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jinzhong Zhang
//         Created:  Wed Sep  9 18:30:00 CEST 2009
// $Id: Analyzer.cc,v 1.7 2009/12/01 10:28:53 zhangjin Exp $
//
//


// C++ include files
#include <memory>
#include <vector>
#include <cmath>
#include <cstring>

//user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Muon
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//HLT
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

//Generation level particles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//ROOT
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include "TFile.h"

//#include "TLorentzVector.h"

//constants decleration
#define maxMuon 50
#define maxHLTnum 200
#define maxFilterObjects 30
#define maxKindsofPV 6
#define MaxHLTObjDeviation 0.003 // the maximum deviation between pt and eta of reco objs and HLT objs to be accepted as a HLT objs
//
// class decleration
//
using namespace edm;
using namespace std;

class Analyzer : public edm::EDAnalyzer {
   public:
      explicit Analyzer(const edm::ParameterSet&);
      ~Analyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  //TTree&TFile
  TFile *file;
  string FileName;
  TTree *fMuon_Tree,*fSummarization_Tree;
  ULong64_t Total_Events;
  //General Event Information
  struct General_EvtInfo
  {
    ULong64_t RunNum,EventNum,LumiBlock;
  } Info;
  //HLT
  unsigned int HLTSize;
  Bool_t HLTacceptance[maxHLTnum],HLTNameSaved;
  //HLTObjects
  unsigned int number_Filters;
  vector<Double_t> *HLTObj_pt[maxFilterObjects],*HLTObj_eta[maxFilterObjects],*HLTObj_phi[maxFilterObjects];
  vector<InputTag> fHLTFilterNames;
  InputTag fTriggerResultsTag,fTriggerEventTag;
  //PrimaryVertex
  unsigned int number_kindsofPV;
  vector<string> fPrimaryVerticesTag;
  vector<Double_t> *vx[maxKindsofPV],*vxError[maxKindsofPV],*vy[maxKindsofPV],*vyError[maxKindsofPV],*vz[maxKindsofPV],*vzError[maxKindsofPV];
  //Generation Level Muon
  vector<Float_t> *Gen_pt,*Gen_eta,*Gen_phi,*AntiGen_pt,*AntiGen_eta,*AntiGen_phi;
  //Muon
  struct Stuff_Muon
  {
    UInt_t number;
    Float_t pt[maxMuon],eta[maxMuon],phi[maxMuon],Vertex[maxMuon][3],isoR03sumPt[maxMuon],isoR03emEt[maxMuon],isoR03hadEt[maxMuon],isoR03hoEt[maxMuon],isoR03nJets[maxMuon],isoR03nTracks[maxMuon],isoR05sumPt[maxMuon],isoR05emEt[maxMuon],isoR05hadEt[maxMuon],isoR05hoEt[maxMuon],isoR05nJets[maxMuon],isoR05nTracks[maxMuon],isoemVetoEt[maxMuon],isohadVetoEt[maxMuon],isohoVetoEt[maxMuon],normalizedChi2[maxMuon],STARecoChi2[maxMuon],TrkRecoChi2[maxMuon];
    Bool_t isHLTObj[maxHLTnum][maxMuon];
  } muon,antimuon,*muonpointer;
};
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Analyzer::Analyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   FileName = iConfig.getParameter<string>("FileName");
   file=new TFile(FileName.c_str(),"RECREATE");
   fMuon_Tree  = new TTree ("Muon","Muons") ;
   fSummarization_Tree  = new TTree ("TotalEvents","Total_Events") ;
   vector< string > default_PVTags;
   default_PVTags.push_back("offlinePrimaryVerticesWithBS");
   default_PVTags.push_back("offlinePrimaryVertices");
   fPrimaryVerticesTag = iConfig.getUntrackedParameter< vector<string> >("PrimaryVertices",default_PVTags);//"offlinePrimaryVerticesWithBS": Primary vertex reconstructed using the tracks taken from the generalTracks collection, and imposing the offline beam spot as a constraint in the fit of the vertex position. Another possible tag is "offlinePrimaryVertices", which is Primary vertex reconstructed using the tracks taken from the generalTracks collection
   number_kindsofPV=fPrimaryVerticesTag.size();
   if (number_kindsofPV>maxKindsofPV)
     {
	printf("You need to increase the constant maxKindsofPV.\n");
	exit(0);
     }
   const InputTag default_TriggerResultsTag("TriggerResults::HLT");
   fTriggerResultsTag = iConfig.getUntrackedParameter<InputTag>("TriggerResultsTag",default_TriggerResultsTag);
   const InputTag default_TriggerEventTag("hltTriggerSummaryAOD","","HLT");
   fTriggerEventTag = iConfig.getUntrackedParameter<InputTag>("triggerEventTag",default_TriggerEventTag);
   fHLTFilterNames = iConfig.getParameter< vector<InputTag> >("hltFilterNames");
   number_Filters = fHLTFilterNames.size();
   if (number_Filters>maxFilterObjects)
     {
	printf("You need to increase the constant maxFilterObjects.\n");
	exit(0);
     }
}

Analyzer::~Analyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  //  fMuon_Tree->MakeClass("MuonClass");
}

// member functions
//

// ------------ method called to for each event  ------------
void
Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	
   unsigned int i,j,k,l,Num_HLTObjs[maxHLTnum];
   Total_Events++;
   //Event Information
   Info.RunNum = iEvent.run();
   Info.EventNum = iEvent.id().event();
   Info.LumiBlock = iEvent.luminosityBlock();
   //HLT Information
   Handle<TriggerResults> hltTriggerResults;
   iEvent.getByLabel(fTriggerResultsTag,hltTriggerResults);
   if (hltTriggerResults.isValid()) 
     {
       if (!HLTNameSaved)
	 {
	   string HLTTriggerNames;
	   TriggerNames HLTNames;
	   HLTNames.init(*hltTriggerResults);
	   HLTSize = hltTriggerResults->size();
	   if (HLTSize>maxHLTnum)
	     {
	       printf("You need to increase the constant maxHLTnum.\n");
	       exit(0);
	     }
	   HLTTriggerNames=HLTNames.triggerName(0)+"/O:";
	   for (i = 1 ; i < HLTSize-1; i++)
	     HLTTriggerNames+=HLTNames.triggerName(i)+":";
	   HLTTriggerNames+=HLTNames.triggerName(i);
	   fMuon_Tree->Branch("HLTAcceptance",HLTacceptance,HLTTriggerNames.c_str());
	   HLTNameSaved=true;
	 }
      for (i = 0 ; i < HLTSize ; i++) 
	HLTacceptance[i]=hltTriggerResults->accept(i);
    }
   else printf("Invalid TriggerResultsTag:\"%s\".\n",fTriggerResultsTag.label().c_str());
       
   //Muon
   Handle<reco::MuonCollection> Muon;
   iEvent.getByLabel("muons", Muon);
   if (!Muon.isValid()) return;
   unsigned char MuonQueue[2][maxMuon],temp;// MuonQueue[0][i] for Muons,MuonQueue[1][i] for AntiMuons,
   reco::MuonCollection const & muons = *Muon;
   reco::Muon iterMuon;
   i=0;muon.number=0;antimuon.number=0;
   for ( reco::MuonCollection::const_iterator iter = muons.begin(); iter != muons.end() ; ++iter)
     { 
       if (iter->isGlobalMuon())
	 {
	   if (iter->charge()==-1) 
	     {
	       MuonQueue[0][muon.number]=i;
	       muon.number++;
	     }
	   else 
	     {
	       MuonQueue[1][antimuon.number]=i;
	       antimuon.number++;
	     }
	   if (muon.number>=maxMuon||antimuon.number>=maxMuon)
	     {
	       printf("You need to increase the constant maxMuon.\n");
	       break;
	     }
	 }
       i++;
     }
   if (!muon.number||!antimuon.number) return;//only keep the events which contains at least one muon and one antimuon
   //HLT Objects
   Handle<trigger::TriggerEvent> trgEvent;
   iEvent.getByLabel(fTriggerEventTag,trgEvent);
   if(trgEvent.isValid())
     {
       const trigger::TriggerObjectCollection& TOC = trgEvent->getObjects();
       int Total_Filters = trgEvent->sizeFilters();
       for( i =0; i < number_Filters;i++)
	 {
	   trigger::size_type pos_Filter = trgEvent->filterIndex(fHLTFilterNames[i]);
	   if (pos_Filter < Total_Filters)
	     {
	       const trigger::Keys& KEYS(trgEvent->filterKeys(pos_Filter));
	       Num_HLTObjs[i] = int(KEYS.size());
	       for(j = 0; j < Num_HLTObjs[i]; j++)
		 {
		   const trigger::TriggerObject& TO = TOC[KEYS[j]];
		   //printf("pt%f\teta%f\tphi%f\n",TO.pt(),TO.eta(),TO.phi());
		   HLTObj_pt[i]->push_back(TO.pt());
		   HLTObj_eta[i]->push_back(TO.eta());
		   HLTObj_phi[i]->push_back(TO.phi());
		 }
	     }
	   else printf("FilterName \"%s\" is not valid.\n",fHLTFilterNames[i].label().c_str());
	 }
     }
   else printf("TriggerEventTag \"%s\" is not valid.",fTriggerEventTag.label().c_str());
   //Save Muon (need HLT object information)
   muonpointer=&muon;
   for (k=0; k<2; k++)
     {
       for (i=0; i<muonpointer->number; i++)
	 for (j=i+1; j<muonpointer->number; j++) 
	   if (muons[MuonQueue[k][i]].pt()<muons[MuonQueue[k][j]].pt())
	     {
	       temp=MuonQueue[k][i];
	       MuonQueue[k][i]=MuonQueue[k][j];
	       MuonQueue[k][j]=temp;
	     }
       for (i=0; i<muonpointer->number; i++)
	 {
	   iterMuon=muons[MuonQueue[k][i]];
	   muonpointer->pt[i]=iterMuon.pt();
	   muonpointer->eta[i]=iterMuon.eta();
	   muonpointer->phi[i]=iterMuon.phi();
	   muonpointer->Vertex[i][0]=iterMuon.vertex().X();
	   muonpointer->Vertex[i][1]=iterMuon.vertex().Y();
	   muonpointer->Vertex[i][2]=iterMuon.vertex().Z();
	   muonpointer->isoR03sumPt[i]=iterMuon.isolationR03().sumPt;
	   muonpointer->isoR03emEt[i]=iterMuon.isolationR03().emEt;
	   muonpointer->isoR03hadEt[i]=iterMuon.isolationR03().hadEt;
	   muonpointer->isoR03hoEt[i]=iterMuon.isolationR03().hoEt;
	   muonpointer->isoR03nJets[i]=iterMuon.isolationR03().nJets;
	   muonpointer->isoR03nTracks[i]=iterMuon.isolationR03().nTracks;
	   muonpointer->isoR05sumPt[i]=iterMuon.isolationR05().sumPt;
	   muonpointer->isoR05emEt[i]=iterMuon.isolationR05().emEt;
	   muonpointer->isoR05hadEt[i]=iterMuon.isolationR05().hadEt;
	   muonpointer->isoR05hoEt[i]=iterMuon.isolationR05().hoEt;
	   muonpointer->isoR05nJets[i]=iterMuon.isolationR05().nJets;
	   muonpointer->isoR05nTracks[i]=iterMuon.isolationR05().nTracks;
	   muonpointer->isoemVetoEt[i]=iterMuon.isolationR03().emVetoEt;
	   muonpointer->isohadVetoEt[i]=iterMuon.isolationR03().hadVetoEt;
	   muonpointer->isohoVetoEt[i]=iterMuon.isolationR03().hoVetoEt;
	   muonpointer->normalizedChi2[i]=iterMuon.globalTrack()->normalizedChi2();
	   muonpointer->STARecoChi2[i]=iterMuon.combinedQuality().staRelChi2;
	   muonpointer->TrkRecoChi2[i]=iterMuon.combinedQuality().trkRelChi2;
	   for (j=0; j<number_Filters; j++)
	     {
	       muonpointer->isHLTObj[j][i]=false;
	       for (l=0; l<Num_HLTObjs[j]; l++)
		 if ((abs(muonpointer->eta[i]-(*HLTObj_eta[j])[l])<MaxHLTObjDeviation)&&(abs(muonpointer->phi[i]-(*HLTObj_phi[j])[l])<MaxHLTObjDeviation)) muonpointer->isHLTObj[j][i]=true;
	     }
	 }
        muonpointer=&antimuon;
     }
   //Primary Vertex
   for (i=0; i < number_kindsofPV; i++)
     {
       Handle<reco::VertexCollection> recVtxs;
       iEvent.getByLabel(fPrimaryVerticesTag[i].c_str(),recVtxs);
       if (recVtxs.isValid())
	 for(reco::VertexCollection::const_iterator v=recVtxs->begin(); v!=recVtxs->end(); ++v)
	   {
	     vx[i]->push_back(v->x());	   vxError[i]->push_back(v->xError());
	     vy[i]->push_back(v->y());	   vyError[i]->push_back(v->yError());
	     vz[i]->push_back(v->z());	   vzError[i]->push_back(v->zError());
	   }
       else printf("%s information is not valid.\n",fPrimaryVerticesTag[i].c_str());
     }
   //Generation level Muons
   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByLabel("genParticles",genParticles);
   if (genParticles.isValid())
     { 
       for(reco::GenParticleCollection::const_iterator genPar=genParticles->begin(); genPar!=genParticles->end(); ++genPar)
	 {
	   if (genPar->pdgId()==13)
	    {
	      Gen_pt->push_back(genPar->pt());
	      Gen_eta->push_back(genPar->eta());
	      Gen_phi->push_back(genPar->phi());
	    }
	  if (genPar->pdgId()==-13)
	    {
	      AntiGen_pt->push_back(genPar->pt());
	      AntiGen_eta->push_back(genPar->eta());
	      AntiGen_phi->push_back(genPar->phi());
	    }
	 }
     }
     else printf("GenerationLevelParticles information is not valid.\n");
   fMuon_Tree->Fill();
   //clear all vectors
   for( i =0; i < number_Filters;i++)
     {
       HLTObj_pt[i]->clear();
       HLTObj_eta[i]->clear();
       HLTObj_phi[i]->clear();
     }
   for (i=0; i < number_kindsofPV; i++)
     {
       vx[i]->clear();   vxError[i]->clear();
       vy[i]->clear();   vyError[i]->clear();
       vz[i]->clear();   vzError[i]->clear();
     }
   Gen_pt->clear();
   AntiGen_pt->clear();
   Gen_eta->clear();
   AntiGen_eta->clear();
   Gen_phi->clear();
   AntiGen_phi->clear();
}


// ------------ method called once each job just before starting event loop  ------------
void 
Analyzer::beginJob()
{ 
  //at the beginning of event, initialize tree variables
  fMuon_Tree->SetCircular(500000);//the max events in a single root file
  Total_Events=0;HLTNameSaved = false;
  unsigned int i;
  //General Run Information
  fMuon_Tree->Branch("Event_Info",&Info.RunNum,"RunNum/l:EventNum:LumiBlock");
  //Primary Vertices
  for( i =0; i < number_kindsofPV;i++)
    {
      vx[i] = new vector<Double_t>();  vxError[i] = new vector<Double_t>();
      vy[i] = new vector<Double_t>();  vyError[i] = new vector<Double_t>();
      vz[i] = new vector<Double_t>();  vzError[i] = new vector<Double_t>();
      fMuon_Tree->Branch((fPrimaryVerticesTag[i]+"_X").c_str(),&vx[i]);
      fMuon_Tree->Branch((fPrimaryVerticesTag[i]+"_eX").c_str(),&vxError[i]);
      fMuon_Tree->Branch((fPrimaryVerticesTag[i]+"_Y").c_str(),&vy[i]);
      fMuon_Tree->Branch((fPrimaryVerticesTag[i]+"_eY").c_str(),&vyError[i]);
      fMuon_Tree->Branch((fPrimaryVerticesTag[i]+"_Z").c_str(),&vz[i]);
      fMuon_Tree->Branch((fPrimaryVerticesTag[i]+"_eZ").c_str(),&vzError[i]);
    }
  //Generation Level Muons
  Gen_pt = new vector<Float_t>();
  Gen_eta = new vector<Float_t>();
  Gen_phi = new vector<Float_t>();
  fMuon_Tree->Branch("GenMuon_pt",&Gen_pt);
  fMuon_Tree->Branch("GenMuon_eta",&Gen_eta);
  fMuon_Tree->Branch("GenMuon_phi",&Gen_phi);
  AntiGen_pt = new vector<Float_t>();
  AntiGen_eta = new vector<Float_t>();
  AntiGen_phi = new vector<Float_t>();
  fMuon_Tree->Branch("GenAntiMuon_pt",&AntiGen_pt);
  fMuon_Tree->Branch("GenAntiMuon_eta",&AntiGen_eta);
  fMuon_Tree->Branch("GenAntiMuon_phi",&AntiGen_phi);
  //Muon
  fMuon_Tree->Branch("Muon_size",&muon.number,"Mnum/i");
  fMuon_Tree->Branch("Muon_pt",muon.pt,"Mpt[Mnum]/F");
  fMuon_Tree->Branch("Muon_eta",muon.eta,"Meta[Mnum]/F");
  fMuon_Tree->Branch("Muon_phi",muon.phi,"Mphi[Mnum]/F");
  fMuon_Tree->Branch("Muon_Vertex",muon.Vertex,"Mphi[Mnum][3]/F");
  fMuon_Tree->Branch("Muon_isoR03sumPt",muon.isoR03sumPt,"MisoR03sumPt[Mnum]/F");
  fMuon_Tree->Branch("Muon_isoR03emEt",muon.isoR03emEt,"MisoR03emEt[Mnum]/F");
  fMuon_Tree->Branch("Muon_isoR03hadEt",muon.isoR03hadEt,"MisoR03hadEt[Mnum]/F");
  fMuon_Tree->Branch("Muon_isoR03hoEt",muon.isoR03hoEt,"MisoR03hoEt[Mnum]/F");
  fMuon_Tree->Branch("Muon_isoR03nJets",muon.isoR03nJets,"MisoR03nJets[Mnum]/F");
  fMuon_Tree->Branch("Muon_isoR03nTracks",muon.isoR03nTracks,"MisoR03nTracks[Mnum]/F");
  fMuon_Tree->Branch("Muon_isoR05sumPt",muon.isoR05sumPt,"MisoR05sumPt[Mnum]/F");
  fMuon_Tree->Branch("Muon_isoR05emEt",muon.isoR05emEt,"MisoR05emEt[Mnum]/F");
  fMuon_Tree->Branch("Muon_isoR05hadEt",muon.isoR05hadEt,"MisoR05hadEt[Mnum]/F");
  fMuon_Tree->Branch("Muon_isoR05hoEt",muon.isoR05hoEt,"MisoR05hoEt[Mnum]/F");
  fMuon_Tree->Branch("Muon_isoR05nJets",muon.isoR05nJets,"MisoR05nJets[Mnum]/F");
  fMuon_Tree->Branch("Muon_isoR05nTracks",muon.isoR05nTracks,"MisoR05nTracks[Mnum]/F");
  fMuon_Tree->Branch("Muon_isoemVetoEt",muon.isoemVetoEt,"MisoemVetoEt[Mnum]/F");
  fMuon_Tree->Branch("Muon_isohadVetoEt",muon.isohadVetoEt,"MisohadVetoEt[Mnum]/F");
  fMuon_Tree->Branch("Muon_isohoVetoEt",muon.isohoVetoEt,"MisohoVetoEt[Mnum]/F");
  fMuon_Tree->Branch("Muon_normalizedChi2",muon.normalizedChi2,"MnormalizedChi2[Mnum]/F");
  fMuon_Tree->Branch("Muon_STARecoChi2",muon.STARecoChi2,"MSTARecoChi2[Mnum]/F");
  fMuon_Tree->Branch("Muon_TrkRecoChi2",muon.TrkRecoChi2,"MTrkRecoChi2[Mnum]/F");
  fMuon_Tree->Branch("AntiMuon_size",&antimuon.number,"Anum/i");
  fMuon_Tree->Branch("AntiMuon_pt",antimuon.pt,"Apt[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_eta",antimuon.eta,"Aeta[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_phi",antimuon.phi,"Aphi[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_Vertex",antimuon.Vertex,"Aphi[Anum][3]/F");
  fMuon_Tree->Branch("AntiMuon_isoR03sumPt",antimuon.isoR03sumPt,"AisoR03sumPt[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isoR03emEt",antimuon.isoR03emEt,"AisoR03emEt[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isoR03hadEt",antimuon.isoR03hadEt,"AisoR03hadEt[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isoR03hoEt",antimuon.isoR03hoEt,"AisoR03hoEt[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isoR03nJets",antimuon.isoR03nJets,"AisoR03nJets[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isoR03nTracks",antimuon.isoR03nTracks,"AisoR03nTracks[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isoR05sumPt",antimuon.isoR05sumPt,"AisoR05sumPt[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isoR05emEt",antimuon.isoR05emEt,"AisoR05emEt[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isoR05hadEt",antimuon.isoR05hadEt,"AisoR05hadEt[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isoR05hoEt",antimuon.isoR05hoEt,"AisoR05hoEt[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isoR05nJets",antimuon.isoR05nJets,"AisoR05nJets[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isoR05nTracks",antimuon.isoR05nTracks,"AisoR05nTracks[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isoemVetoEt",antimuon.isoemVetoEt,"AisoemVetoEt[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isohadVetoEt",antimuon.isohadVetoEt,"AisohadVetoEt[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_isohoVetoEt",antimuon.isohoVetoEt,"AisohoVetoEt[Anum]/F");
  fMuon_Tree->Branch("AntiMuon_normalizedChi2",muon.normalizedChi2,"AnormalizedChi2[Mnum]/F");
  fMuon_Tree->Branch("AntiMuon_STARecoChi2",antimuon.STARecoChi2,"ASTARecoChi2[Mnum]/F");
  fMuon_Tree->Branch("AntiMuon_TrkRecoChi2",antimuon.TrkRecoChi2,"ATrkRecoChi2[Mnum]/F");
  //HLT Objects
  for( i =0; i < number_Filters;i++)
    {
      HLTObj_pt[i]=new vector<Double_t>();
      HLTObj_eta[i]=new vector<Double_t>();
      HLTObj_phi[i]=new vector<Double_t>();
      fMuon_Tree->Branch((fHLTFilterNames[i].label()+"_pt").c_str(),&HLTObj_pt[i]);
      fMuon_Tree->Branch((fHLTFilterNames[i].label()+"_eta").c_str(),&HLTObj_eta[i]);
      fMuon_Tree->Branch((fHLTFilterNames[i].label()+"_phi").c_str(),&HLTObj_phi[i]);
      fMuon_Tree->Branch(("Muon_Is"+fHLTFilterNames[i].label()+"Obj").c_str(),&muon.isHLTObj[i],("Mis"+fHLTFilterNames[i].label()+"Obj[Mnum]/O").c_str());
      fMuon_Tree->Branch(("AntiMuon_Is"+fHLTFilterNames[i].label()+"Obj").c_str(),&antimuon.isHLTObj[i],("Ais"+fHLTFilterNames[i].label()+"Obj[Mnum]/O").c_str());
    }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analyzer::endJob()
{
   //at the end of the job record the total events went through for calculating the weight number which is lum*cross*FilterEff/Total_Events_Number
   fSummarization_Tree->Branch("Total_Events",&Total_Events,"TotalEvents/l");
   fSummarization_Tree->Fill();
   file->Write();
   file->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
//LogInfo("MuonSelector") << "TEST " <<iter->pt()*cosh(iter->eta())<<"  "<<sqrt(iter->px()*iter->px()+iter->py()*iter->py()+iter->pz()*iter->pz());
