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
// $Id: Analyzer.cc,v 1.6 2009/10/22 15:27:38 zhangjin Exp $
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
   unsigned int Total_Events;
   //HLT
   unsigned int HLTSize;
   bool HLTacceptance[maxHLTnum],HLTNameSaved;
   //HLTObjects
   unsigned int number_Filters;
   vector<double> *pt[maxFilterObjects],*eta[maxFilterObjects],*phi[maxFilterObjects];
   vector<InputTag> fHLTFilterNames;
   InputTag fTriggerResultsTag,fTriggerEventTag;
   //PrimaryVertex
   unsigned int number_kindsofPV;
   vector<string> fPrimaryVerticesTag;
   vector<double> *vx[maxKindsofPV],*vxError[maxKindsofPV],*vy[maxKindsofPV],*vyError[maxKindsofPV],*vz[maxKindsofPV],*vzError[maxKindsofPV];
   //Generation Level Muon
   vector<float> *Gen_pt,*Gen_eta,*Gen_phi,*AntiGen_pt,*AntiGen_eta,*AntiGen_phi;
   //Muon
   struct Stuff_Muon
   {
	unsigned int number;
	float pt[maxMuon],eta[maxMuon],phi[maxMuon],Vertex[maxMuon][3],isoR03sumPt[maxMuon],isoR03emEt[maxMuon],isoR03hadEt[maxMuon],isoR03hoEt[maxMuon],isoR03nJets[maxMuon],isoR03nTracks[maxMuon],isoR05sumPt[maxMuon],isoR05emEt[maxMuon],isoR05hadEt[maxMuon],isoR05hoEt[maxMuon],isoR05nJets[maxMuon],isoR05nTracks[maxMuon],isoemVetoEt[maxMuon],isohadVetoEt[maxMuon],isohoVetoEt[maxMuon];
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
  unsigned int i,j,k,number_Keys=0;
   Total_Events++;
   //int eventNum = iEvent.id().event();
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
       if (iter->isGlobalMuon()&&iter->pt()>20)
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
	 }
        muonpointer=&antimuon;
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
	       number_Keys = int(KEYS.size());
	       for(j = 0; j < number_Keys; j++)
		 {
		   const trigger::TriggerObject& TO = TOC[KEYS[j]];
		   //printf("pt%f\teta%f\tphi%f\n",TO.pt(),TO.eta(),TO.phi());
		   pt[i]->push_back(TO.pt());
		   eta[i]->push_back(TO.eta());
		   phi[i]->push_back(TO.phi());
		 }
	     }
	   else printf("FilterName \"%s\" is not valid.\n",fHLTFilterNames[i].label().c_str());
	 }
     }
   else printf("TriggerEventTag \"%s\" is not valid.",fTriggerEventTag.label().c_str());
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
	      //printf("Got Muon, pt: %f;eta: %f;phi: %f;",genPar->pt(),genPar->eta(),genPar->phi());
	    }
	  if (genPar->pdgId()==-13)
	    {
	      AntiGen_pt->push_back(genPar->pt());
	      AntiGen_eta->push_back(genPar->eta());
	      AntiGen_phi->push_back(genPar->phi());
	      //printf("Got AntiMuon, pt: %f;eta: %f;phi: %f;\n",genPar->pt(),genPar->eta(),genPar->phi());
	    }
	 }
     }
     else printf("GenerationLevelParticles information is not valid.\n");
   fMuon_Tree->Fill();
   //clear all vectors
   for( i =0; i < number_Filters;i++)
     {
       pt[i]->clear();
       eta[i]->clear();
       phi[i]->clear();
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
  //HLT Objects
  unsigned int i;  
  string HLTTriggerNames;
  for( i =0; i < number_Filters;i++)
    {
      pt[i]=new vector<double>();
      eta[i]=new vector<double>();
      phi[i]=new vector<double>();
      HLTTriggerNames=fHLTFilterNames[i].label()+"_pt";
      fMuon_Tree->Branch(HLTTriggerNames.c_str(),&pt[i]);
      HLTTriggerNames=fHLTFilterNames[i].label()+"_eta";
      fMuon_Tree->Branch(HLTTriggerNames.c_str(),&eta[i]);
      HLTTriggerNames=fHLTFilterNames[i].label()+"_phi";
      fMuon_Tree->Branch(HLTTriggerNames.c_str(),&phi[i]);
      }
  //Primary Vertices
  for( i =0; i < number_kindsofPV;i++)
    {
      vx[i] = new vector<double>();  vxError[i] = new vector<double>();
      vy[i] = new vector<double>();  vyError[i] = new vector<double>();
      vz[i] = new vector<double>();  vzError[i] = new vector<double>();
      fMuon_Tree->Branch((fPrimaryVerticesTag[i]+"_X").c_str(),&vx[i]);
      fMuon_Tree->Branch((fPrimaryVerticesTag[i]+"_eX").c_str(),&vxError[i]);
      fMuon_Tree->Branch((fPrimaryVerticesTag[i]+"_Y").c_str(),&vy[i]);
      fMuon_Tree->Branch((fPrimaryVerticesTag[i]+"_eY").c_str(),&vyError[i]);
      fMuon_Tree->Branch((fPrimaryVerticesTag[i]+"_Z").c_str(),&vz[i]);
      fMuon_Tree->Branch((fPrimaryVerticesTag[i]+"_eZ").c_str(),&vzError[i]);
    }
  //Generation Level Muons
  Gen_pt = new vector<float>();
  Gen_eta = new vector<float>();
  Gen_phi = new vector<float>();
  fMuon_Tree->Branch("GenMuon_pt",&Gen_pt);
  fMuon_Tree->Branch("GenMuon_eta",&Gen_eta);
  fMuon_Tree->Branch("GenMuon_phi",&Gen_phi);
  AntiGen_pt = new vector<float>();
  AntiGen_eta = new vector<float>();
  AntiGen_phi = new vector<float>();
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
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analyzer::endJob()
{
   //at the end of the job record the total events went through for calculating the weight number which is lum*cross*FilterEff/Total_Events_Number
   fSummarization_Tree->Branch("Total_Events",&Total_Events,"TotalEvents/i");
   fSummarization_Tree->Fill();
   file->Write();
   file->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
//LogInfo("MuonSelector") << "TEST " <<iter->pt()*cosh(iter->eta())<<"  "<<sqrt(iter->px()*iter->px()+iter->py()*iter->py()+iter->pz()*iter->pz());
