#include "Math/SVector.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TPolyLine.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TArc.h"
#include "TText.h"
#include "TColor.h"
#include <cstdlib>
#include "setTDRStyle.h"
#include <fstream>

using namespace ROOT::Math;
using namespace std;

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

#define Cu_ssm 2.43E-3
#define Cd_ssm 3.13E-3
#define Cu_psi 7.9E-4
#define Cd_psi 7.9E-4
//#define RotateCooridnate
//#define DrawObsLimit
#ifdef DrawObsLimit
#define DrawAllowedRegion
#endif
//#define LogScale
#ifdef RotateCooridnate 
  #define Angle 0.295
  #define XMinimum -2E-7
  #define XMaximum 3E-7
  #define YMinimum 1E-10
  #ifdef DrawObsLimit
    #define YMaximum 4E-5
  #else 
    #define YMaximum 3E-5
  #endif
#else
  #define Angle 0
  #define ShowRotateMark 0.295
  #ifdef LogScale
    #define XMinimum 5E-7
    #define XMaximum 2E-5
  #else
    #ifdef ShowRotateMark
       #define XMinimum -3E-7
    #else
       #define XMinimum 0
    #endif
    #define XMaximum 1E-5
  #endif
    #define YMinimum 1E-6
    #define YMaximum 3E-5
#endif

#define NPoints 1000
#define EpsilonIncrement 1E-5

const Double_t MassLimit[]={400,500,600,700,800,900,1000};
const UInt_t NContours=sizeof(MassLimit)/sizeof(Double_t);
const EColor ColorList[]={EColor(kMagenta-8),EColor(kMagenta-7),EColor(kMagenta-3),kMagenta,EColor(kMagenta+3)};
const Byte_t numColor=sizeof(ColorList)/sizeof(EColor);

inline void RotateCoor(Double_t & x,Double_t & y, Double_t thisAngle=Angle) {
  Double_t xprime=cos(thisAngle)*x-sin(thisAngle)*y,yprime=sin(thisAngle)*x+cos(thisAngle)*y;
  x=xprime;y=yprime;
}

#define CdMinimum ( cos(Angle)*XMinimum+sin(Angle)*YMinimum )
#define CuMinimum ( -sin(Angle)*XMinimum+cos(Angle)*YMinimum )

inline void DrawTitle() {
  TLatex *t = new TLatex(XMinimum,YMaximum*1.01,"CMS");
  t->SetTextSize(0.03);
  t->SetTextAlign(11);
  t->Draw();
  t = new TLatex(XMaximum,YMaximum*1.01,"#int Ldt=4.98/5.28fb^{-1}, e^{+}e^{-}/mu^{+}#mu^{-}");
  t->SetTextSize(0.03);
  t->SetTextAlign(31);
  t->Draw();
}

//first is Cd, second is Cu
inline pair<Double_t,Double_t> CuCdCalculation(double MZp2,double Epsilon) {
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
  double DecayWidth_All=3*DecayWidth_neutrinos+3*DecayWidth_ll+2*DecayWidth_uu+3*DecayWidth_dd+DecayWidth_tt+DecayWidth_WW;
  // clog << "Br(Z'->vv) = " << DecayWidth_neutrinos/DecayWidth_All*100 << "%"<<endl;
  // clog << "Br(Z'->ll) = " << DecayWidth_ll/DecayWidth_All*100 << "%"<<endl;
  // clog << "Br(Z'->uu,cc) = " << DecayWidth_uu/DecayWidth_All*100 << "%"<<endl;
  // clog << "Br(Z'->dd,ss,bb) = " << DecayWidth_dd/DecayWidth_All*100 << "%"<<endl;
  // clog << "Br(Z'->tt) = " << DecayWidth_tt/DecayWidth_All*100 << "%"<<endl;
  //clog << "Br(Z'->WW) = " << gM2/(12*M_PI)*sqrt(MZp2)*DecayWidth_WW <<"/"<< gM2/(12*M_PI)*sqrt(MZp2)*DecayWidth_All<<"="<<DecayWidth_WW/DecayWidth_All*100 << "%"<<endl;
  double fator4Cucd=2*gM2*DecayWidth_ll/DecayWidth_All;
  return pair<double,double>(fator4Cucd*(Vf(3)*Vf(3)+Af(3)*Af(3)),fator4Cucd*(Vf(2)*Vf(2)+Af(2)*Af(2)));
}

//I ignore the factor of the Wu and Wd calculated here: Wu(d)_real=Wu(d)_here*48s/M_PI*1E-3(the unit in the files are mbar)
Bool_t CuCdContourLineCal(Float_t mass, Double_t & a, Double_t & b) {
  if (Cu_ssm*Cd_psi==Cd_ssm*Cu_psi) {
    cerr<<"Not able to calculate Wu and Wd based on these two models (Cu_1*Cd_2==Cu_2*Cd_1)"<<endl;
    return 0;
  }
  ifstream file_ssm,file_psi,file_limit;
  //file_limit.open("dilepton_explimit_median_16feb2012v1.txt");
  file_limit.open("ll_mass_cteq_low_5ifb_meidan.txt");
  file_ssm.open("zprime_ssm_cteq6l1_masswindow.txt");
  file_psi.open("zprime_psi_cteq6l1_masswindow.txt");
  //file_ssm.open("zprime_ssm_cteq6l1_all.txt");
  //file_psi.open("zprime_psi_cteq6l1_all.txt");
  Double_t xsec_ssm,xsec_psi,xsec_limit;
  Float_t tmp_mass=-10.;
  while ( !file_ssm.eof() && tmp_mass!=mass )
    file_ssm >> tmp_mass >> xsec_ssm;
  tmp_mass=-10.;
  while ( !file_psi.eof() && tmp_mass!=mass )
    file_psi >> tmp_mass >> xsec_psi;
  tmp_mass=-10.;
  while ( !file_limit.eof() && tmp_mass!=mass )
    file_limit>>tmp_mass>>xsec_limit;
  if ( file_ssm.eof() || file_psi.eof() || file_limit.eof() ) return false;
  Double_t kfactor= 1.23312 - 0.154901*(mass-1000.0)/1000.0 +0.0516781*(mass-1000.0)/1000.0*(mass-1000.0)/1000.0,
           Wu=kfactor*(xsec_ssm/Cd_ssm-xsec_psi/Cd_psi)/(Cu_ssm/Cd_ssm-Cu_psi/Cd_psi),
           Wd=kfactor*(xsec_ssm/Cu_ssm-xsec_psi/Cu_psi)/(Cd_ssm/Cu_ssm-Cd_psi/Cu_psi);
  xsec_limit*=1E-6;//970.E-9;//Z xsec is 970E-9 mbar
  a=xsec_limit/Wu;
  b=Wd/Wu;
  pair<double,double> CuCdValues=CuCdCalculation(mass*mass,0.06);
  double Cd=CuCdValues.first,Cu=CuCdValues.second;
  printf("mass=%8f, kfactor=%.2f, Wu=%E, Wd=%E, Cu=%E-%E*Cd, xsec_limit=%E, NNLO xsec_stuzp=%E\n",mass,kfactor,Wu,Wd,a,b,xsec_limit,Cu*Wu+Cd*Wd);
  return true;
}

TMultiGraph *MLimitContours() {
  TMultiGraph *CuCd = new TMultiGraph("CuCdPlot","StSM Z Prime CuCd Plot;C_{d};C_{u}");
  CuCd->SetMinimum(YMinimum);
  Double_t Cu[NContours][NPoints],Cd[NContours][NPoints];
  Double_t a,b;
  for (UInt_t n=0;n<NContours;n++) 
    if ( CuCdContourLineCal(MassLimit[n],a,b) ) {
      Double_t Increment=((a-CuMinimum)/b-CdMinimum)/Double_t(NPoints-10);
      cerr<<Increment<<endl;
      for (Int_t i=0;i<NPoints;i++) {
	Cd[n][i]=CdMinimum+Increment*i;
	Cu[n][i]=a-b*Cd[n][i];
	RotateCoor(Cd[n][i],Cu[n][i]);
      }
      TGraph *temp=new TGraph(NPoints,Cd[n],Cu[n]);
      //    temp->SetLineWidth(1);
      CuCd->Add(temp);
      char textatpoint[30];
      sprintf(textatpoint," %.0f GeV",MassLimit[n]);
      TText *t = new TText(XMinimum,temp->Eval(XMinimum),textatpoint);
      t->SetTextSize(0.03);
      t->Draw();
    }
  return CuCd;
}

void DrawText(Float_t X,Float_t Y,Float_t MZp,Float_t Epsilon, EColor color=kBlack, Float_t fontsize=0.04, Byte_t align=13) {
  char limittext[30];
  if (MZp>0&&Epsilon>0) sprintf(limittext," (%.0f GeV/c^{2}, %.3f)",MZp,Epsilon);
  else if (Epsilon>0) sprintf(limittext,"#epsilon = %.3f",Epsilon);
  else if (MZp>0) sprintf(limittext,"%.0f",MZp);
  TLatex *t = new TLatex(X,Y,limittext);
  t->SetTextAlign(align);//Top Left
  t->SetTextColor(color);
  t->SetTextSize(fontsize);
  t->Draw();
}

//fixed epsilon lines
inline void DrawFixedEpsilonLines(Double_t *Masspoints, Double_t *Theo_EpsilonLimit, Byte_t NMasspoints) {
  Byte_t iep=0;
  for (Double_t Epsilon=0.02;Epsilon<0.061;Epsilon+=0.01) {
    Double_t X[NMasspoints],Y[NMasspoints];
    Short_t Npoints=0;
    for ( Byte_t n=0; n<NMasspoints; n++ ) 
      if (Epsilon<Theo_EpsilonLimit[n]){
	pair<double,double> CuCdValues=CuCdCalculation(Masspoints[n]*Masspoints[n],Epsilon);
	if (CuCdValues.first>CdMinimum&&CuCdValues.second>CuMinimum) {
	  X[Npoints]=CuCdValues.first;
	  Y[Npoints]=CuCdValues.second;
	  RotateCoor(X[Npoints],Y[Npoints]);
	  Npoints++;
	}
      }
    TPolyLine *epline = new TPolyLine(Npoints,X,Y);
    epline->SetLineColor(ColorList[iep%numColor]);
    epline->SetLineWidth(2);
    epline->Draw("same");
    DrawText(X[0],Y[0],-1,Epsilon,ColorList[iep%numColor]);
    iep++;
  }
}
void OneMassLimit(Float_t Mass)
{
  Double_t a,b; 
  if ( !CuCdContourLineCal(Mass,a,b) ) {
    cerr<<"*** Not enough data to calculate the limit of this mass point. ***"<<endl;
    return;
  }
  setTDRStyle();
  gStyle->SetFillColor(10);//box fill color to white
  TCanvas *c1=new TCanvas("CuCd","CuCd",500,500);
  TLegend *Legend = new TLegend(0.4,0.75,0.6,0.9);
  char tmptxt[100];
#ifdef RotateCooridnate
  sprintf(tmptxt,"StSM Z Prime CuCd Plot;cos(%.3f)C_{d}-sin(%.3f)C_{u};sin(%.3f)C_{d}+cos(%.3f)C_{u}",Angle,Angle,Angle,Angle);
#else
  sprintf(tmptxt,"StSM Z Prime CuCd Plot;C_{u};C_{d}");
#endif
  TH1F * frame=c1->DrawFrame(XMinimum,YMinimum,XMaximum,YMaximum,tmptxt);
  frame->GetXaxis()->SetMoreLogLabels(true);
  frame->GetXaxis()->SetLabelSize(0.03);
  frame->GetXaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetMoreLogLabels(true);
  frame->GetYaxis()->SetLabelSize(0.03);
  frame->GetYaxis()->SetTitleSize(0.04);
  Double_t Cu[Int_t(0.06/EpsilonIncrement)],Cd[Int_t(0.06/EpsilonIncrement)];
  Double_t Increment=((a-CuMinimum)/b-CdMinimum)/Double_t(NPoints-10),TextYPos=YMinimum;
  for (Int_t i=0;i<NPoints;i++) {
    Cd[i]=CdMinimum+Increment*i;
    Cu[i]=a-b*Cd[i];
    RotateCoor(Cd[i],Cu[i]);
    if (Cd[i]<XMinimum) TextYPos=Cu[i];
  }
  TGraph *LimitContour=new TGraph(NPoints,Cd,Cu);
  //    temp->SetLineWidth(1);
  LimitContour->Draw("C");
  char textatpoint[30];
  sprintf(textatpoint," %.0f GeV",Mass);
  TLatex *t = new TLatex(XMinimum,TextYPos,textatpoint);
  t->SetTextSize(0.03);
  t->Draw();
  Double_t Exp_EpsilonLimit=0.062,Theo_EpsilonLimit=0.062,MZp2=Mass*Mass;
  for (Double_t Epsilon=0.001;Epsilon<0.062;Epsilon+=EpsilonIncrement) {
    Double_t MOne2=MZp2*(MZp2-MZ2)/(MZp2-MZ2+Epsilon*Epsilon*(MZp2-MW2));
    if ( Epsilon>0.061*sqrt(1-MZ2/MOne2) ) break;// the theo upperbound on epsilon
    Theo_EpsilonLimit=Epsilon;
    pair<double,double> CuCdValues=CuCdCalculation(MZp2,Epsilon);
    if ( CuCdValues.second<=a-b*CuCdValues.first ) Exp_EpsilonLimit=Epsilon;
    Epsilon+=EpsilonIncrement;
  }
  
  Short_t NPolyGonVertex=0,NLabels=0;
  Double_t LabelX[10],LabelY[10];
  for (Double_t Epsilon=Theo_EpsilonLimit;Epsilon>Exp_EpsilonLimit;Epsilon-=EpsilonIncrement) {
    pair<double,double> CuCdValues=CuCdCalculation(Mass*Mass,Epsilon);
    if (CuCdValues.first>CdMinimum&&CuCdValues.second>CuMinimum) {
      Cd[NPolyGonVertex]=CuCdValues.first;
      Cu[NPolyGonVertex]=CuCdValues.second;
      RotateCoor(Cd[NPolyGonVertex],Cu[NPolyGonVertex]);
      if ( (Epsilon>0.03-EpsilonIncrement/2&&Epsilon<0.03+EpsilonIncrement/2)||
	   (Epsilon>0.04-EpsilonIncrement/2&&Epsilon<0.04+EpsilonIncrement/2)||
	   (Epsilon>0.05-EpsilonIncrement/2&&Epsilon<0.05+EpsilonIncrement/2)||
	   (Epsilon>0.06-EpsilonIncrement/2&&Epsilon<0.06+EpsilonIncrement/2)
	   ) {
	LabelX[NLabels]=Cd[NPolyGonVertex];
	LabelY[NLabels]=Cu[NPolyGonVertex];
	DrawText(LabelX[NLabels],LabelY[NLabels],-1,Epsilon);
	NLabels++;
      }
      NPolyGonVertex++;
    }
  }
  if (NPolyGonVertex>0) {
    TGraph *ExcludedCuCd=new TGraph(NPolyGonVertex,Cd,Cu);
    ExcludedCuCd->SetLineColor(kGray);
    ExcludedCuCd->SetLineWidth(1);
    ExcludedCuCd->Draw("same,C");
    TGraph *EpsilonLabel=new TGraph(NLabels,LabelX,LabelY);
    EpsilonLabel->SetMarkerColor(kGray);
    EpsilonLabel->SetLineColor(kGray);
    EpsilonLabel->SetMarkerStyle(2);
    EpsilonLabel->SetLineWidth(1);
    EpsilonLabel->Draw("same,P");
    Legend->AddEntry(EpsilonLabel,"Excluded Region","LP");
  }

  NPolyGonVertex=0;NLabels=0;
  for (Double_t Epsilon=0.001;Epsilon<=Exp_EpsilonLimit;Epsilon+=EpsilonIncrement) {
    pair<double,double> CuCdValues=CuCdCalculation(Mass*Mass,Epsilon);
    if (CuCdValues.first>CdMinimum&&CuCdValues.second>CuMinimum) {
      Cd[NPolyGonVertex]=CuCdValues.first;
      Cu[NPolyGonVertex]=CuCdValues.second;
      RotateCoor(Cd[NPolyGonVertex],Cu[NPolyGonVertex]);
      if ( (Epsilon>0.03-EpsilonIncrement/2&&Epsilon<0.03+EpsilonIncrement/2)||
	   (Epsilon>0.04-EpsilonIncrement/2&&Epsilon<0.04+EpsilonIncrement/2)||
	   (Epsilon>0.05-EpsilonIncrement/2&&Epsilon<0.05+EpsilonIncrement/2)||
	   (Epsilon>0.06-EpsilonIncrement/2&&Epsilon<0.06+EpsilonIncrement/2)
	   ) {
	LabelX[NLabels]=Cd[NPolyGonVertex];
	LabelY[NLabels]=Cu[NPolyGonVertex];
	DrawText(LabelX[NLabels],LabelY[NLabels],-1,Epsilon);
	NLabels++;
      }
      NPolyGonVertex++;
    }
  }
  if (NPolyGonVertex>0) {
    TGraph *AllowedCuCd=new TGraph(NPolyGonVertex,Cd,Cu);
    AllowedCuCd->SetLineColor(kGreen);
    AllowedCuCd->SetLineWidth(1);
    AllowedCuCd->Draw("same,C");
    TGraph *EpsilonLabel=new TGraph(NLabels,LabelX,LabelY);
    EpsilonLabel->SetMarkerColor(kGreen);
    EpsilonLabel->SetLineColor(kGreen);
    EpsilonLabel->SetMarkerStyle(2);
    EpsilonLabel->SetLineWidth(1);
    EpsilonLabel->Draw("same,P");
    Legend->AddEntry(EpsilonLabel,"Allowed Region","LP");
  }

  Legend->Draw();
  DrawTitle();
  sprintf(tmptxt,"Stueckelberg Model, M_{Z'} = %.0f GeV/c^{2}",Mass);
  t = new TLatex( (XMinimum+XMaximum)/2.,YMaximum*0.98,tmptxt);
  t->SetTextSize(0.02);
  t->SetTextAlign(23);
  t->Draw();
}

void FullLimit()
{
  setTDRStyle();
  gStyle->SetFillColor(10);//box fill color to white
  TCanvas *c1=new TCanvas("CuCd","CuCd",500,500);
  char tmptxt[100];
#ifdef RotateCooridnate
  sprintf(tmptxt,"StSM Z Prime CuCd Plot;cos(%.3f)C_{d}-sin(%.3f)C_{u};sin(%.3f)C_{d}+cos(%.3f)C_{u}",Angle,Angle,Angle,Angle);
#else
  sprintf(tmptxt,"StSM Z Prime CuCd Plot;C_{u};C_{d}");
#endif
  TH1F * frame=c1->DrawFrame(XMinimum,YMinimum,XMaximum,YMaximum,tmptxt);
  frame->GetXaxis()->SetMoreLogLabels(true);
  frame->GetXaxis()->SetLabelSize(0.03);
  frame->GetXaxis()->SetTitleSize(0.04);
  frame->GetXaxis()->SetTitleOffset(1.02);
  frame->GetYaxis()->SetMoreLogLabels(true);
  frame->GetYaxis()->SetLabelSize(0.03);
  frame->GetYaxis()->SetTitleSize(0.04);
#ifdef DrawObsLimit
  TMultiGraph *CuCd = MLimitContours();
  CuCd->Draw("C");
#endif
  //calculate epsilon limit for each mass  
  const UInt_t NMasspoints=(MassLimit[NContours-1]-MassLimit[0])/5+1;
  Double_t Masspoints[NMasspoints],Exp_EpsilonLimit[NMasspoints],Theo_EpsilonLimit[NMasspoints];
  for ( UInt_t n=0; n<NMasspoints; n++ ) {
    Masspoints[n]=MassLimit[0]+5*n;
    Double_t MZp2=Masspoints[n]*Masspoints[n];//Mplus2
    Double_t Epsilon=0.001;
    Theo_EpsilonLimit[n]=0.061;
    Exp_EpsilonLimit[n]=-1.;//<0 means not available
    Double_t a,b;
    if ( CuCdContourLineCal(Masspoints[n],a,b) ) Exp_EpsilonLimit[n]=1.;//it has exp limit
    do {
      Double_t MOne2=MZp2*(MZp2-MZ2)/(MZp2-MZ2+Epsilon*Epsilon*(MZp2-MW2));
      if ( Epsilon>0.061*sqrt(1-MZ2/MOne2) ) break;// the theo upperbound on epsilon
      Theo_EpsilonLimit[n]=Epsilon;
#ifdef DrawAllowedRegion
      if ( Exp_EpsilonLimit[n]>0. ) {
	pair<double,double> CuCdValues=CuCdCalculation(MZp2,Epsilon);
	Double_t Cd=CuCdValues.first,Cu=CuCdValues.second;
	if ( Cu<=a-b*Cd ) Exp_EpsilonLimit[n]=Epsilon;
      }
#else
      Exp_EpsilonLimit[n]=Epsilon;
#endif
      Epsilon+=EpsilonIncrement;
    }while(true);
    if ( Exp_EpsilonLimit[n]>0. ) cerr<<"exp ep lim="<<Exp_EpsilonLimit[n]<<endl;
  }

  //experimental allowed StuZp CuCd or the whole EWK allowed Stu CuCd
  Double_t polyX[Int_t(0.13/EpsilonIncrement)+NMasspoints],polyY[Int_t(0.13/EpsilonIncrement)+NMasspoints];

  Short_t NPolyGonVertex=0;
  //the lowest mass, right boarder
  if ( Exp_EpsilonLimit[0]>0. ) 
    for (Double_t Epsilon=0.001;Epsilon<=Exp_EpsilonLimit[0];Epsilon+=EpsilonIncrement) {
      pair<double,double> CuCdValues=CuCdCalculation(Masspoints[0]*Masspoints[0],Epsilon);
      if (CuCdValues.first>CdMinimum&&CuCdValues.second>CuMinimum) {
	polyX[NPolyGonVertex]=CuCdValues.first;
	polyY[NPolyGonVertex]=CuCdValues.second;
	RotateCoor(polyX[NPolyGonVertex],polyY[NPolyGonVertex]);
  	NPolyGonVertex++;
      }
    }
  // DrawText(polyX[NPolyGonVertex-1],polyY[NPolyGonVertex-1],Masspoints[0],EpsilonLimit[0]);
  // boarder for each mass
  for ( UInt_t n=0; n<NMasspoints; n++ ) 
    if ( Exp_EpsilonLimit[n]>0. ) {
      pair<double,double> CuCdValues=CuCdCalculation(Masspoints[n]*Masspoints[n],Exp_EpsilonLimit[n]);
      polyX[NPolyGonVertex]=CuCdValues.first;
      polyY[NPolyGonVertex]=CuCdValues.second;
      RotateCoor(polyX[NPolyGonVertex],polyY[NPolyGonVertex]);
      NPolyGonVertex++;
    }
  //the highest mass, left boarder
  if ( Exp_EpsilonLimit[NMasspoints-1]>0. )
    for (Double_t Epsilon=Exp_EpsilonLimit[NMasspoints-1];Epsilon>=0.001;Epsilon-=EpsilonIncrement) {
      pair<double,double> CuCdValues=CuCdCalculation(Masspoints[NMasspoints-1]*Masspoints[NMasspoints-1],Epsilon);
      if (CuCdValues.first>CdMinimum&&CuCdValues.second>CuMinimum) {
	polyX[NPolyGonVertex]=CuCdValues.first;
	polyY[NPolyGonVertex]=CuCdValues.second;
	RotateCoor(polyX[NPolyGonVertex],polyY[NPolyGonVertex]);
	NPolyGonVertex++;
      }
    }

  TPolyLine *pline_inclu = new TPolyLine(NPolyGonVertex,polyX,polyY);
  pline_inclu->SetFillColor(kGreen);
  pline_inclu->SetFillStyle(3005);
  pline_inclu->SetLineColor(kGreen);
  pline_inclu->SetLineWidth(1);
#ifdef DrawAllowedRegion
  pline_inclu->Draw("f,same");
#endif
  pline_inclu->Draw("same");

#ifdef ShowRotateMark
  NPolyGonVertex=0;
  //the lowest mass, right boarder
  for (Double_t Epsilon=0.001;Epsilon<=Theo_EpsilonLimit[0];Epsilon+=EpsilonIncrement) {
    pair<double,double> CuCdValues=CuCdCalculation(Masspoints[0]*Masspoints[0],Epsilon);
    if (CuCdValues.first>CdMinimum&&CuCdValues.second>CuMinimum) {
      polyX[NPolyGonVertex]=CuCdValues.first;
      polyY[NPolyGonVertex]=CuCdValues.second;
      RotateCoor(polyX[NPolyGonVertex],polyY[NPolyGonVertex],ShowRotateMark);
      NPolyGonVertex++;
    }
  }
  
  // DrawText(polyX[NPolyGonVertex-1],polyY[NPolyGonVertex-1],Masspoints[0],EpsilonLimit[0]);
  // boarder for each mass
  for ( UInt_t n=0; n<NMasspoints; n++ )  {
    pair<double,double> CuCdValues=CuCdCalculation(Masspoints[n]*Masspoints[n],Theo_EpsilonLimit[n]);
    polyX[NPolyGonVertex]=CuCdValues.first;
    polyY[NPolyGonVertex]=CuCdValues.second;
    RotateCoor(polyX[NPolyGonVertex],polyY[NPolyGonVertex],ShowRotateMark);
    NPolyGonVertex++;
  }
  //the highest mass, left boarder
  for (Double_t Epsilon=Theo_EpsilonLimit[NMasspoints-1];Epsilon>=0.001;Epsilon-=EpsilonIncrement) {
    pair<double,double> CuCdValues=CuCdCalculation(Masspoints[NMasspoints-1]*Masspoints[NMasspoints-1],Epsilon);
    if (CuCdValues.first>CdMinimum&&CuCdValues.second>CuMinimum) {
      polyX[NPolyGonVertex]=CuCdValues.first;
      polyY[NPolyGonVertex]=CuCdValues.second;
      RotateCoor(polyX[NPolyGonVertex],polyY[NPolyGonVertex],ShowRotateMark);
      NPolyGonVertex++;
    }
  }
  TPolyLine *pline_rotated = new TPolyLine(NPolyGonVertex,polyX,polyY);
  pline_rotated->SetLineColor(kGray+3);
  pline_rotated->SetLineStyle(2);
  pline_rotated->Draw("same");
  TArc *Indicator=new TArc(0,0,8E-6,90.1-ShowRotateMark/M_PI*180,90.1);
  Indicator->SetFillStyle(0);
  Indicator->SetLineColor(kGray+3);
  Indicator->SetLineStyle(2);
  Indicator->SetNoEdges(true);
  Indicator->Draw("same");
  sprintf(tmptxt,"#splitline{%.3frad}{(%.1f#circ)}",ShowRotateMark,ShowRotateMark/M_PI*180);
  TArrow *arrow=new TArrow(0,8E-6,1E-8,8E-6,0.02,"<|");
  arrow->SetAngle(20);
  arrow->SetLineColor(kGray+3);
  arrow->Draw();
  TLatex *t_marker = new TLatex( 1.5E-6,9E-6,tmptxt);
  t_marker->SetTextSize(0.035);
  t_marker->SetTextAlign(21);
  t_marker->Draw("same");
#endif

#ifdef DrawAllowedRegion
  //experimental excluded StuZp CuCd
  NPolyGonVertex=0;
  //the lowest mass, right boarder
  if ( Exp_EpsilonLimit[0]>0. ) 
    for (Double_t Epsilon=Exp_EpsilonLimit[0];Epsilon<=Theo_EpsilonLimit[0];Epsilon+=EpsilonIncrement) {
      pair<double,double> CuCdValues=CuCdCalculation(Masspoints[0]*Masspoints[0],Epsilon);
      if (CuCdValues.first>CdMinimum&&CuCdValues.second>CuMinimum) {
	polyX[NPolyGonVertex]=CuCdValues.first;
	polyY[NPolyGonVertex]=CuCdValues.second;
	RotateCoor(polyX[NPolyGonVertex],polyY[NPolyGonVertex]);
	NPolyGonVertex++;
      }
    }
  // DrawText(polyX[NPolyGonVertex-1],polyY[NPolyGonVertex-1],Masspoints[0],EpsilonLimit[0]);
  // upper boarder for each mass
  for ( UInt_t n=0; n<NMasspoints; n++ ) {
    pair<double,double> CuCdValues=CuCdCalculation(Masspoints[n]*Masspoints[n],Theo_EpsilonLimit[n]);
    polyX[NPolyGonVertex]=CuCdValues.first;
    polyY[NPolyGonVertex]=CuCdValues.second;
    RotateCoor(polyX[NPolyGonVertex],polyY[NPolyGonVertex]);
    NPolyGonVertex++;
  }
  //the highest mass, left boarder
  if ( Exp_EpsilonLimit[NMasspoints-1]>0. )
    for (Double_t Epsilon=Theo_EpsilonLimit[NMasspoints-1];Epsilon>=Exp_EpsilonLimit[NMasspoints-1];Epsilon-=EpsilonIncrement) {
      pair<double,double> CuCdValues=CuCdCalculation(Masspoints[NMasspoints-1]*Masspoints[NMasspoints-1],Epsilon);
      if (CuCdValues.first>CdMinimum&&CuCdValues.second>CuMinimum) {
	polyX[NPolyGonVertex]=CuCdValues.first;
	polyY[NPolyGonVertex]=CuCdValues.second;
	RotateCoor(polyX[NPolyGonVertex],polyY[NPolyGonVertex]);
	NPolyGonVertex++;
      }
    }
  // lower boarder for each mass
  for ( Int_t n=NMasspoints-1; n>=0; n-- ) 
    if ( Exp_EpsilonLimit[n]>0. ) {
      pair<double,double> CuCdValues=CuCdCalculation(Masspoints[n]*Masspoints[n],Exp_EpsilonLimit[n]);
      polyX[NPolyGonVertex]=CuCdValues.first;
      polyY[NPolyGonVertex]=CuCdValues.second;
      RotateCoor(polyX[NPolyGonVertex],polyY[NPolyGonVertex]);
      NPolyGonVertex++;
    }

  TPolyLine *pline_exclu = new TPolyLine(NPolyGonVertex,polyX,polyY);
  pline_exclu->SetFillColor(kGray);
  pline_exclu->SetFillStyle(3004);
  pline_exclu->SetLineColor(kGray+3);
  pline_exclu->SetLineWidth(1);
  pline_exclu->Draw("f,same");
  pline_exclu->Draw("same");
#endif

  //fixed mass lines
#ifdef RotateCooridnate
  DrawFixedEpsilonLines(Masspoints,Theo_EpsilonLimit,NMasspoints);
  Bool_t aligntop=true;
#endif
  for (UInt_t n=0;n<NContours;n++) {
    NPolyGonVertex=0;
    Int_t MassPointIndex=NMasspoints-1;
    for (  ;MassPointIndex>=0; MassPointIndex-- )
      if (Masspoints[MassPointIndex]==MassLimit[n]) break;
    if ( MassPointIndex>=0 )
      for (Double_t Epsilon=Theo_EpsilonLimit[MassPointIndex];Epsilon>=0.001;Epsilon-=EpsilonIncrement) {
	pair<double,double> CuCdValues=CuCdCalculation(MassLimit[n]*MassLimit[n],Epsilon);
	if (CuCdValues.first>CdMinimum&&CuCdValues.second>CuMinimum) {
	  polyX[NPolyGonVertex]=CuCdValues.first;
	  polyY[NPolyGonVertex]=CuCdValues.second;
	  RotateCoor(polyX[NPolyGonVertex],polyY[NPolyGonVertex]);
	  NPolyGonVertex++;
	}
      }
    TPolyLine *massline = new TPolyLine(NPolyGonVertex,polyX,polyY);
#ifdef RotateCooridnate
    massline->SetLineColor(kGray+1);
    massline->SetLineWidth(1);
#else
    massline->SetLineColor(kGray+1);
    massline->SetLineWidth(0.5);
#endif
    massline->Draw("same");
#ifdef RotateCooridnate
    if (aligntop) DrawText(polyX[0],polyY[0],MassLimit[n],-1,EColor(kGray+1),0.03,11);
    else DrawText(polyX[0],polyY[0],MassLimit[n],-1,EColor(kGray+1),0.03,13);
    aligntop=!aligntop;
#endif
  }

#ifndef RotateCooridnate
  DrawFixedEpsilonLines(Masspoints,Theo_EpsilonLimit,NMasspoints);
#endif

#ifdef DrawAllowedRegion
  TLegend *Legend = new TLegend(0.4,0.75,0.6,0.9);
  Legend->AddEntry(pline_exclu,"Excluded Region","F");
  Legend->AddEntry(pline_inclu,"Allowed Region","F");
  Legend->Draw();
#endif

#ifdef DrawObsLimit
  DrawTitle();
#endif

  sprintf(tmptxt,"#splitline{Stueckelberg Model}{%.0f GeV/c^{2} #geq M_{Z'} #geq %.0f GeV/c^{2}}",MassLimit[NContours-1],MassLimit[0]);
#ifndef LogScale
  #define TITLEXPOS (XMinimum+XMaximum)/2.
#else
  #define TITLEXPOS sqrt(XMinimum*XMaximum)
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
#endif

  TLatex *t = new TLatex( TITLEXPOS,YMaximum*0.98,tmptxt);
  t->SetTextSize(0.035);
  t->SetTextAlign(23);
  t->Draw();
}
