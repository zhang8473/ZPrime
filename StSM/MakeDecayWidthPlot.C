#include "Math/SVector.h"
#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TLegend.h"
#include "setTDRStyle.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

using namespace ROOT::Math;
using namespace std;

const EColor ColorList[]={kRed,kGreen,kBlue,EColor(kOrange+7),kCyan,kMagenta,kGray};
const UInt_t n_Epsilon=5;
#define SW2 0.2312 //square of Weinberg Angle (theoretical)
//#define MW 80.398 //mass of W (could be determined from Mz and SW2)
#define CW2 (1-SW2)
#define MZ 91.1876 //mass of Z
#define CW sqrt(CW2)
#define SW sqrt(SW2)
#define TW2 (SW2/CW2)
#define TW (SW/CW)
#define MZ2 (MZ*MZ)
#define MW2 (MZ2*CW2) //determined from Mz and SW2
#define Mt 171.2//mass of t
#define GF 1.16637E-5
#define gM2 (1.414213562*GF*MZ*MZ)
#define alphas 0.1184

inline Double_t DecayWidth(double MZp2,double Epsilon) {
  double MOne2=MZp2*(MZp2-MZ2)/(MZp2-MZ2+Epsilon*Epsilon*(MZp2-MW2));
  double tmp1=MOne2*Epsilon*(MW2-MZp2)/MW2/TW/(MOne2-MZp2);
  double tmp2=(MW2-MZp2)/MW2/TW;
  double R31 = 1/sqrt(tmp1*tmp1+tmp2*tmp2+1);
  double R21 = R31*(MW2-MZp2)/MW2/TW;
  SVector <double,4> T3l(0.5,-0.5,0.5,-0.5);//T3L of neutrinos,leptons,u,d
  SVector <double,4> Q(0,-1,2/3.,-1/3.);//charges of neutrinos,leptons,u,d
  SVector <double,4> Vf,Af;
  Vf=(CW*R31-SW*R21)*T3l+2*Q*SW*R21;
  Af=(CW*R31-SW*R21)*T3l;
  //clog <<"===== Mz': "<<sqrt(MZp2)<<"; M1: "<<sqrt(MOne2)<<"; Epsilon: "<<Epsilon<<" ====="<<endl;
  // clog <<"SW: "<<SW<<"; CW: "<<CW<<"; TW: "<<TW<<endl;
  // clog <<"MZp(M+): "<<MZp<<"; M-: "<<sqrt(Mminus2)<<endl;
  // clog<<"R matrix:"<<endl<<R<<endl;
  // clog << "af (nu,l,u,d): "<<Af<<endl<<"vf (nu,l,u,d): "<<Vf<<endl;
  //the following decay width neglect the factor gM2/(12*M_PI)*sqrt(MZp2)
  double DecayWidth_neutrinos=Vf(0)*Vf(0)+Af(0)*Af(0);
  double DecayWidth_ll=Vf(1)*Vf(1)+Af(1)*Af(1);//ee mumu tautau
  double DecayWidth_uu=3*(Vf(2)*Vf(2)+Af(2)*Af(2))*(1+alphas/M_PI);//uu cc
  double DecayWidth_dd=3*(Vf(3)*Vf(3)+Af(3)*Af(3))*(1+alphas/M_PI);//dd ss bb
  double DecayWidth_tt=MZp2>=4*Mt*Mt?3*sqrt(1-4*Mt*Mt/MZp2)*(Vf(2)*Vf(2)*(1+2*Mt*Mt/MZp2)+Af(2)*Af(2)*(1-4*Mt*Mt/MZp2))*(1+alphas/M_PI):0;// tt
  double DecayWidth_WW=MZp2>=4*MW2?CW2*R31*R31/4*MZp2*MZp2/MW2/MW2*pow(1-4*MW2/MZp2,1.5)*(1+20*MW2/MZp2+12*MW2*MW2/MZp2/MZp2):0;// WW
  return gM2/(12*M_PI)*sqrt(MZp2)*(3*DecayWidth_neutrinos+3*DecayWidth_ll+2*DecayWidth_uu+3*DecayWidth_dd+DecayWidth_tt+DecayWidth_WW);
  // clog << "Br(Z'->vv) = " << DecayWidth_neutrinos/DecayWidth_All*100 << "%"<<endl;
  // clog << "Br(Z'->ll) = " << DecayWidth_ll/DecayWidth_All*100 << "%"<<endl;
  // clog << "Br(Z'->uu,cc) = " << DecayWidth_uu/DecayWidth_All*100 << "%"<<endl;
  // clog << "Br(Z'->dd,ss,bb) = " << DecayWidth_dd/DecayWidth_All*100 << "%"<<endl;
  // clog << "Br(Z'->tt) = " << DecayWidth_tt/DecayWidth_All*100 << "%"<<endl;
  //clog << "Br(Z'->WW) = " << gM2/(12*M_PI)*sqrt(MZp2)*DecayWidth_WW <<"/"<< gM2/(12*M_PI)*sqrt(MZp2)*DecayWidth_All<<"="<<DecayWidth_WW/DecayWidth_All*100 << "%"<<endl;
}

void Make(){
  setTDRStyle();
  gStyle->SetFillColor(10);//box fill color to white
  gStyle->SetOptFit(0);
  char filename[200];
  unsigned char numColor=sizeof(ColorList)/sizeof(EColor);
  TLegend *Legends = new TLegend(0.75,0.92-0.05*n_Epsilon,0.95,0.92);
  TCanvas *canvas = new TCanvas("dwplot","dwplot",500,500);
  canvas->DrawFrame(180,0,1500,0.1,"Steuckelburg Z' decay width;M_{Z'} (GeV/c^{2});#Gamma(pp->StuZ'->l^{+}l^{-}) (GeV/c^{2})");
  for (UInt_t i_Epsilon=0; i_Epsilon<n_Epsilon; i_Epsilon++) {
    sprintf(filename,"zprime_stu_ep0%.0f0_cteq6l1_decaywidth.txt",2.+i_Epsilon);
    ifstream in;
    in.open(filename);
    Float_t mass[200],masserr[200],dw[200],dwerr[200],dw_calculated[200];
    Int_t npoints = 0;
    printf("<<<< %s >>>>\n",filename);
    while (1) {
      in >> mass[npoints] >> dw[npoints] >> dwerr[npoints];
      masserr[npoints]=0;
      dw_calculated[npoints]=DecayWidth(mass[npoints]*mass[npoints],(2.+i_Epsilon)*0.01);
      if (!in.good()) {
	printf("mass=%8f, dw=%8f, dwerr=%8f\n",mass[npoints-1],dw[npoints-1],dwerr[npoints-1]);
	break;
      }
      else npoints++;
    }
    TGraphErrors *gr = new TGraphErrors(npoints,mass,dw,masserr,dwerr);
    gr->SetName(filename);
    gr->SetLineColor(ColorList[i_Epsilon%numColor]);
    gr->SetMarkerColor(ColorList[i_Epsilon%numColor]);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.5);
    gr->Draw("same");
    //TF1 *linear = new TF1("linear","[0]+[1]*x", 180, 1500);
    //gr->Fit(linear,"R+");
    TGraph *gr_cal = new TGraph(npoints,mass,dw_calculated);
    gr_cal->SetLineColor(kBlack);
    gr_cal->Draw("L");
    char legendtext[30];
    sprintf(legendtext,"Z'_{StSM} #epsilon=0.0%.0f",2.+i_Epsilon);
    Legends->AddEntry(gr,legendtext,"L");
  }
  //  TCanvas *canvas=new TCanvas("","",500,500);
  Legends->Draw();
}
