#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include "TApplication.h"
//#include "/home/zhangjin/root/montecarlo/pythia6/inc/TPythia6.h"
#include "TPythia6.h"
#include "TFile.h"
#include "TLegend.h"
//#include "TError.h"
//#include "TTree.h"
//#include "TClonesArray.h"
//#include "TH1.h"
//#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "setTDRStyle.h"


//#include "Riostream.h"
#include <cstdlib>

using namespace ROOT::Math;
using namespace std;

#define SW2 0.2312 //square of Weinberg Angle (theoretical)
//#define MW 80.398  //mass of W (could be determined from Mz and SW2)
#define MZ 91.1876 //mass of Z
#define Mt 171.2//mass of t
#define GF 1.16637E-5
#define gM2 1.414213562*GF*MZ*MZ
#define sqrts 14000
// This function just load the needed libraries if we're executing from
// an interactive session.
void loadLibraries()
{
#ifdef __CINT__
  // Load the Event Generator abstraction library, Pythia 6
  // library, and the Pythia 6 interface library.
  printf("Loading libs....................");
  gSystem->Load("libEG");
  gSystem->Load("$ROOTSYS/../pythia6/libPythia6"); //change to your setup
  gSystem->Load("libEGPythia6");
#endif
}

void pythiaUESettingsBlock(TPythia6* &pythia) {
  //parameters from http://cmslxr.fnal.gov/lxr/source/Configuration/Generator/python/PythiaUESettings_cfi.py, meanings see the same page
  pythia->SetMSTJ(11,3);
  pythia->SetMSTJ(22,2);
  pythia->SetPARJ(71,10);
  pythia->SetMSTP(2,1);
  pythia->SetMSTP(33,0);
  // Exteral PDF CTEQ6L1, I do not have. It gives little difference in the cross section
  //  pythia->SetMSTP(51,10042);
  //  pythia->SetMSTP(52,2);
  pythia->SetMSTP(81,1);
  pythia->SetMSTP(82,4);
  pythia->SetMSTU(21,1);
  pythia->SetPARP(82,1.8387);
  pythia->SetPARP(89,1960.);
  pythia->SetPARP(83,0.5);
  pythia->SetPARP(84,0.4);
  pythia->SetPARP(90,0.16);
  pythia->SetPARP(67,2.5);
  pythia->SetPARP(85,1.0);
  pythia->SetPARP(86,1.0);
  pythia->SetPARP(62,1.25);
  pythia->SetPARP(64,0.2);
  pythia->SetMSTP(91,1);
  pythia->SetPARP(91,2.1);
  pythia->SetPARP(93,15.0);
}

// nEvents is how many events we want.
void makeEventSample(Int_t nEvents)
{
  // Load needed libraries
  loadLibraries();

  setTDRStyle();
  gStyle->SetFillColor(10);//box fill color to white

  double CW2=1-SW2;
  double CW=sqrt(CW2);
  double SW=sqrt(SW2);
  double TW2=SW2/CW2;
  double TW=SW/CW;
  double MZ2=MZ*MZ;
  double MW2=MZ2*CW2;//determined from Mz and SW2

  double MZp[]={200,250,300,350,400,450,500,550,600,650,700,750,800,1000,1200,1500,2000,3000,3500};
  const Byte_t n_MZp=sizeof(MZp)/sizeof(double);
  double Epsilon[]={0.01,0.02,0.03,0.04,0.05,0.06};
  const Byte_t n_Epsilon=sizeof(Epsilon)/sizeof(double);

  double sigma[n_Epsilon][n_MZp];
  TMultiGraph *cross_section_plot = new TMultiGraph("sigmaplot","Steuckelburg Z' cross sections;M_{Z'};#sigma(pp->StuZ'->#mu^{+}#mu^{-}) (fb)");
  const EColor ColorList[]={kRed,kGreen,kBlue,kCyan,EColor(kOrange+7),kMagenta,kGray};
  unsigned char numColor=sizeof(ColorList)/sizeof(EColor);
  TLegend *Legends = new TLegend(0.75,0.92-0.05*n_Epsilon,0.95,0.92);

  for (UInt_t i_Epsilon=0; i_Epsilon<n_Epsilon; i_Epsilon++) {
    double Epsilon2=Epsilon[i_Epsilon]*Epsilon[i_Epsilon];
    for (UInt_t i_MZp=0; i_MZp<n_MZp; i_MZp++) {
      double MZp2=MZp[i_MZp]*MZp[i_MZp];//Mplus2
      double MOne2=MZp2*(MZp2-MZ2)/(MZp2-MZ2+Epsilon2*(MZp2-MW2));
      double MZZp2=MZ2+MOne2*(1+Epsilon2);
      double Mminus2=0.5*(MZZp2-sqrt(MZZp2*MZZp2-4*MOne2*(MZ2+MW2*Epsilon2)));
      double tmp1=MOne2*Epsilon[i_Epsilon]*(MW2-MZp2)/MW2/TW/(MOne2-MZp2);
      double tmp2=(MW2-MZp2)/MW2/TW;
      double tmp3=MOne2*Epsilon[i_Epsilon]*(MW2-Mminus2)/MW2/TW/(MOne2-Mminus2);
      double tmp4=(MW2-Mminus2)/MW2/TW;
      SMatrix<double,3> R;
      R(2,0) = 1/sqrt(tmp1*tmp1+tmp2*tmp2+1);
      R(1,0) = R(2,0)*(MW2-MZp2)/MW2/TW;
      R(2,1) = 1/sqrt(tmp3*tmp3+tmp4*tmp4+1);
      R(1,1) = R(2,1)*(MW2-Mminus2)/MW2/TW;
      SVector <double,4> T3l(0.5,-0.5,0.5,-0.5);//T3L of neutrinos,leptons,u,d
      SVector <double,4> Q(0,-1,2/3.,-1/3.);//charges of neutrinos,leptons,u,d
      SVector <double,4> Vf,Af;
      Vf=(CW*R(2,0)-SW*R(1,0))*T3l+2*Q*SW*R(1,0);
      Af=(CW*R(2,0)-SW*R(1,0))*T3l;
      clog <<"M1: "<<sqrt(MOne2)<<"; Epsilon: "<<Epsilon<<endl;
      clog <<"SW: "<<SW<<"; CW: "<<CW<<"; TW: "<<TW<<endl;
      clog <<"MZp(M+): "<<MZp[i_MZp]<<"; M-: "<<sqrt(Mminus2)<<endl;
      clog<<"R matrix:"<<endl<<R<<endl;
      Vf*=2;Af*=2;//For PYTHIA you have to mutiple them by 2, see http://www.hep.uiuc.edu/home/catutza/nota12.ps
      // for WW, see http://www-spires.slac.stanford.edu/spires/find/hep/www?j=ZEPYA%2CC45%2C109
      // and PYTHIA 6.4 Manual p188
      double theStrengthOfWWcoupling=MZp[i_MZp]>=2*sqrt(MW2)?(MZp2/MW2)*sqrt((1+(1/TW2)*(1+Epsilon2))*TW)*R(2,0):0;
      clog << "af (nu,l,u,d): "<<Af<<endl<<"vf (nu,l,u,d): "<<Vf<<endl;
      clog << "WW ratio:" <<theStrengthOfWWcoupling<<endl;

      // Create an instance of the Pythia event generator ...
      TPythia6* pythia = new TPythia6;

      // ... and initialise it to run p+p at sqrt(700) GeV in Center of Mass System
      pythiaUESettingsBlock(pythia);//CMS pythiaUESettingsBlock 
      pythia->SetMSEL(0); //turn off processes
      pythia->SetMSTU(12,12345); //no title page is written
      pythia->SetMSUB(141,1); //Turn on Drell-Yan SubProcess f+fbar-->gamma/Z0/Zp0, same as MSEL=21
      pythia->SetMSTP(44,3); //gamma/Z/Zp: 3=> Zp only
      pythia->SetPMAS(32,1,MZp[i_MZp]);//M_Zprime for Zprime
      //First generation
      for (int i=3;i>=0;i--) {
	pythia->SetPARU(128-i*2-1,Vf(i));
	pythia->SetPARU(128-i*2,Af(i));
      }
      //Second generation
      for (int i=3;i>=0;i--) {
	pythia->SetPARJ(187-i*2-1,Vf(i));
	pythia->SetPARJ(187-i*2,Af(i));
      }  
      //Third generation
      for (int i=3;i>=0;i--)
	if (i!=2) {
	  pythia->SetPARJ(195-i*2-1,Vf(i));
	  pythia->SetPARJ(195-i*2,Af(i));
	}
	else {
	  pythia->SetPARJ(190,(MZp[i_MZp]>=2*Mt)?Vf(2)*sqrt(sqrt(1-4*Mt*Mt/MZp2)*(1+2*Mt*Mt/MZp2)):0);
	  pythia->SetPARJ(191,(MZp[i_MZp]>=2*Mt)?Af(2)*sqrt(sqrt(1-4*Mt*Mt/MZp2)*(1-4*Mt*Mt/MZp2)):0);
	}
      //WW
      pythia->SetPARU(129,theStrengthOfWWcoupling);
      //Channels
      pythia->SetMDME(289,1,0);//d dbar
      pythia->SetMDME(290,1,0);//u ubar
      pythia->SetMDME(291,1,0);//s sbar
      pythia->SetMDME(292,1,0);//c cbar
      pythia->SetMDME(293,1,0);//b bbar
      pythia->SetMDME(294,1,0);//t tbar
      pythia->SetMDME(295,1,-1);//4th gen Q Qbar
      pythia->SetMDME(296,1,-1);//4th gen Q Qbar
      pythia->SetMDME(297,1,0);//ee
      pythia->SetMDME(298,1,0);//neutrino e e
      pythia->SetMDME(299,1,1);//mu mu
      pythia->SetMDME(300,1,0);//neutrino mu mu
      pythia->SetMDME(301,1,0);//tau tau
      pythia->SetMDME(302,1,0);//neutrino tau tau
      pythia->SetMDME(303,1,-1);//4th generation lepton l- l+'
      pythia->SetMDME(304,1,-1);// 4th generation neutrino nu nubar
      pythia->SetMDME(305,1,-1);// W+ W-'
      pythia->SetMDME(306,1,-1);// H  charged higgs'
      pythia->SetMDME(307,1,-1);// Z0 gamma'
      pythia->SetMDME(308,1,-1);// Z0 h0
      pythia->SetMDME(309,1,-1);// sm higgs
      pythia->SetMDME(310,1,-1);// weird neutral higgs H0A0
      // Now we make some events
      pythia->Initialize("cms", "p", "p", sqrts);
      clog << setfill('0');
      for (Int_t i = 1; i <= nEvents; i++) {
	// Show how far we got
	if ( i % (nEvents/20)  == 0) clog<<"\r"<<"Calculating "<<i<<"/"<<nEvents<<"("<<i/Float_t(nEvents)*100<<"%)";
	pythia->GenerateEvent(); // Make one event.
      }
      printf("\n");
      //Get Cross Section
      pythia->Pystat(1);
      Pyint5_t *pint5=pythia->GetPyint5();
      sigma[i_Epsilon][i_MZp]=pint5->XSEC[2][141]*1E12;//convert mb to pb
      printf("epsilon=%.2f,MZp=%.1f: sigma=%.4ffb\n",Epsilon[i_Epsilon],MZp[i_MZp],sigma[i_Epsilon][i_MZp]);
    }
    TGraph *one_epsilon = new TGraph(n_MZp,MZp,sigma[i_Epsilon]);
    one_epsilon->SetLineColor(ColorList[i_Epsilon%numColor]);
    cross_section_plot->Add(one_epsilon);
    char legendtext[30];
    sprintf(legendtext,"Z'_{StSM} #epsilon=%.2f",Epsilon[i_Epsilon]);
    Legends->AddEntry(one_epsilon,legendtext,"L");
  }
  char namestring[200];
  sprintf(namestring,"CrossSection_%dGeV",sqrts);
  TCanvas *canvas=new TCanvas(namestring,namestring,500,500);
  sprintf(namestring,"%s.pdf",namestring);
  cross_section_plot->Draw("AC*");
  gPad->SetLogy();
  Legends->Draw();
  canvas->SaveAs(namestring);
      //printf("cross section:%E\n",sigma);

      /*
      // Open an output file
      TFile* file = TFile::Open(FILENAME, "RECREATE");
      if (!file || !file->IsOpen()) {
      Error("makeEventSample", "Couldn't open file %s", FILENAME);
      return 1;
      }

      // Make a tree in that file ...
      TTree* tree = new TTree(TREENAME, "Pythia 6 tree");

      // ... and register a the cache of pythia on a branch (It's a
      // TClonesArray of TMCParticle objects. )
      TClonesArray* particles = (TClonesArray*)pythia->GetListOfParticles();
      tree->Branch(BRANCHNAME, &particles);
      */
}

