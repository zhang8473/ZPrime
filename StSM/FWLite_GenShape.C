#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include <iostream>
#include <vector>
#include <cmath>
#endif

using namespace std;
#define Muon_Mass 0.105658367

void FWLite_GenShape(char *filename)
{
  TFile  * file = TFile::Open(filename);
  char result_filename[300];
  strncpy(result_filename,filename,strstr(filename,".root")-filename);
  strcpy(result_filename+strlen(filename)-5,"_GenShape.root\0");
  TFile * result = new TFile(result_filename,"RECREATE");
  Float_t M_Zp=0;
  char *mass_start=strstr(filename,"_M")+2,mass_str[40];
  Byte_t mass_len=strstr(mass_start,"_")-mass_start;
  strncpy(mass_str,mass_start,mass_len);
  mass_str[mass_len]='\0';
  printf("mass=%sGeV/c2\n",mass_str);
  if ( EOF==sscanf(mass_str, "%f", &M_Zp) ) return;
  TTree * dimuon_mass = new TTree("T","Dimuons Mass");
  Float_t dimu_mass=0.;
  dimuon_mass->Branch("dimuon_mass",&dimu_mass,"dimu_mass/F");
  TH1D * hist_Dmu = new TH1D("Statistic","Dimuons Mass", 2000, M_Zp-10., M_Zp+10. );
  fwlite::Event ev(file);
  //loop through each event
  for( ev.toBegin();! ev.atEnd();++ev)
    {
      fwlite::Handle<edm::HepMCProduct> HepMCH;
      HepMCH.getByLabel(ev,"generator");
      if (!HepMCH.isValid() ) continue;
      HepMC::GenEvent HepGenEvent(*(HepMCH->GetEvent()));
      //loop through each Muon
      vector<TLorentzVector *> ZprimeMuons;
      for (HepMC::GenEvent::particle_iterator GenParticle_iter = HepGenEvent.particles_begin(); GenParticle_iter != HepGenEvent.particles_end();GenParticle_iter++) 
	if (abs((*GenParticle_iter)->pdg_id())==13) { 
	HepMC::GenVertex *thisVtx=(*GenParticle_iter)->production_vertex();
	if (thisVtx) 
	  for (HepMC::GenVertex::particles_in_const_iterator pgenD = thisVtx->particles_in_const_begin(); pgenD != thisVtx->particles_in_const_end(); ++pgenD)
	    if ((*pgenD)->pdg_id()==32||(*pgenD)->pdg_id()==23) //if the parent particle is Z or Z'
	      if ( abs((*GenParticle_iter)->momentum().eta())<2.4 ) {//if it is inside the acceptance
		ZprimeMuons.push_back( new TLorentzVector() );
		ZprimeMuons.back()->SetXYZM((*GenParticle_iter)->momentum().px(),(*GenParticle_iter)->momentum().py(),(*GenParticle_iter)->momentum().pz(),Muon_Mass);
	      }
	}
      if (ZprimeMuons.size()>=2) {
	if (ZprimeMuons.size()>2) cerr<<"Warning: Something is wrong, more than two zprime muons are found."<<endl;
	else {
	  dimu_mass=( (*ZprimeMuons[0])+(*ZprimeMuons[1]) ).M();
	  hist_Dmu->Fill( dimu_mass );
	  dimuon_mass->Fill();
	}
      }
    }//end event loop
   hist_Dmu->Write();
   dimuon_mass->Write();
   result->Close();
}
