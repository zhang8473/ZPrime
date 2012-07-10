#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <utility>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

using namespace ROOT::Math;
using namespace std;

//Standard Model Constants
#define MuMu//choices are MuMu,ee,eMu,WW,All
#define SW2 0.2312 //square of Weinberg Angle (theoretical)
//#define MW 80.398  //mass of W (could be determined from Mz and SW2)
#define MZ 91.1876 //mass of Z
#define Mt 171.2//mass of t
 
void Cal(double MZp, double Epsilon)//Stuckelberg Model Constants, To avoid zero denominator, do not set MOne=MZ or Epsilon=0
{
  double CW2=1-SW2;
  double CW=sqrt(CW2);
  double SW=sqrt(SW2);
  double TW2=SW2/CW2;
  double TW=SW/CW;
  double MZ2=MZ*MZ;
  double MW2=MZ2*CW2;//determined from Mz and SW2
  double Epsilon2=Epsilon*Epsilon;
  double MZp2=MZp*MZp;//Mplus2
  double MOne2=MZp2*(MZp2-MZ2)/(MZp2-MZ2+Epsilon2*(MZp2-MW2));
  double MZZp2=MZ2+MOne2*(1+Epsilon2);
  double Mminus2=0.5*(MZZp2-sqrt(MZZp2*MZZp2-4*MOne2*(MZ2+MW2*Epsilon2)));
  double tmp1=MOne2*Epsilon*(MW2-MZp2)/MW2/TW/(MOne2-MZp2);
  double tmp2=(MW2-MZp2)/MW2/TW;
  double tmp3=MOne2*Epsilon*(MW2-Mminus2)/MW2/TW/(MOne2-Mminus2);
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
  clog <<"MZp(M+): "<<MZp<<"; M-: "<<sqrt(Mminus2)<<endl;
  clog<<"R matrix:"<<endl<<R<<endl;
  Vf*=2;Af*=2;//For PYTHIA you have to mutiple them by 2, see http://www.hep.uiuc.edu/home/catutza/nota12.ps  
  double theStrengthOfWWcoupling=MZp>=2*sqrt(MW2)?(MZp2/MW2)*sqrt((1+(1/TW2)*(1+Epsilon2))*TW)*R(2,0):0;
  clog << "af (nu,l,u,d): "<<Af<<endl<<"vf (nu,l,u,d): "<<Vf<<endl;
  clog << "WW ratio:" <<theStrengthOfWWcoupling<<endl;

  //Write a PYTHIA 6.4 CARD
  map<int, pair<string,int> > DecayChannels;
#if defined(MuMu) || defined(ee) || defined(WW)
  DecayChannels[289]=make_pair("d dbar",0);
  DecayChannels[290]=make_pair("u ubar",0);
  DecayChannels[291]=make_pair("s sbar",0);
  DecayChannels[292]=make_pair("c cbar",0);
  DecayChannels[293]=make_pair("b bbar",0);
  DecayChannels[294]=make_pair("t tbar",0);
  DecayChannels[295]=make_pair("4th gen Q Qbar",-1);
  DecayChannels[296]=make_pair("4th gen Q Qbar",-1);
  DecayChannels[298]=make_pair("neutrino e e",0);
  DecayChannels[300]=make_pair("neutrino mu mu",0);
  DecayChannels[301]=make_pair("tau tau",0);
  DecayChannels[302]=make_pair("neutrino tau tau",0);
  DecayChannels[303]=make_pair("4th generation lepton l- l+",-1);
  DecayChannels[304]=make_pair("4th generation neutrino nu nubar",-1);
  DecayChannels[306]=make_pair("H  charged higgs",-1);
  DecayChannels[307]=make_pair("Z0 gamma",-1);
  DecayChannels[308]=make_pair("Z0 h0",-1);
  DecayChannels[309]=make_pair("sm higgs",-1);
  DecayChannels[310]=make_pair("weird neutral higgs H0A0",-1);
#endif

#ifdef MuMu
  DecayChannels[297]=make_pair("e e",0);
  DecayChannels[299]=make_pair("mu mu",1);
  DecayChannels[305]=make_pair("W+ W-",-1);
  const char *name="MuMu";
#endif

#ifdef ee
  DecayChannels[297]=make_pair("e e",1);
  DecayChannels[299]=make_pair("mu mu",0);
  DecayChannels[305]=make_pair("W+ W-",-1);
  const char *name="EE";
#endif

#ifdef eMu
  DecayChannels[297]=make_pair("e e",1);
  DecayChannels[299]=make_pair("mu mu",1);
  DecayChannels[305]=make_pair("W+ W-",-1);
  const char *name="eMu";
#endif

#ifdef WW
  if (MZp<2*MZ*CW) {printf("Off-Shell for WW Channel!\n"); return;}
  DecayChannels[297]=make_pair("e e",0);
  DecayChannels[299]=make_pair("mu mu",0);
  DecayChannels[305]=make_pair("W+ W-",1);
  const char *name="WW";
#endif

#ifdef All
  if (MZp<2*Mt) DecayChannels[294]=make_pair("t tbar",0);
  if (MZp>=2*MZ*CW) DecayChannels[305]=make_pair("W+ W-",1);
  const char *name="All";
#endif

  char ss[255];
  sprintf(ss,"StuZprimeTo%s_M%0.0f_Epsilon%0.3f_7TeV_pythia6_cff.py",name,MZp,Epsilon);
  string CardName(ss);
  CardName.erase(CardName.find(".")-1,2);
  ofstream PYCard(CardName.c_str());
  PYCard<<"#Stueckelburg Z prime pythia6 card"<<endl
        <<"import FWCore.ParameterSet.Config as cms"<<endl
	<<"source = cms.Source(\"EmptySource\")"<<endl
	<<"from Configuration.Generator.PythiaUEZ2Settings_cfi import *"<<endl
	<<"generator = cms.EDFilter(\"Pythia6GeneratorFilter\","<<endl
	<<"pythiaHepMCVerbosity = cms.untracked.bool(False),"<<endl
	<<"maxEventsToPrint = cms.untracked.int32(0),"<<endl
	<<"pythiaPylistVerbosity = cms.untracked.int32(1),"<<endl
	<<"filterEfficiency = cms.untracked.double(1.),"<<endl
	<<"comEnergy = cms.double(7000.0),"<<endl
	<<"crossSection = cms.untracked.double(0),"<<endl
	<<"PythiaParameters = cms.PSet("<<endl
        <<"pythiaUESettingsBlock,"<<endl
	<<"processParameters = cms.vstring("<<endl;
  PYCard <<"'MSEL=0                  !turn off processes',"<<endl;
  PYCard <<"'MSTU(12)=12345          !no title page is written',"<<endl;
  PYCard <<"'MSUB(141)=1             !Turn on Drell-Yan SubProcess f+fbar-->gamma/Z0/Zp0, same as MSEL=21',"<<endl;
  PYCard <<"'PMAS(32,1)="<<MZp<<"          !M_Zprime for Zprime',"<<endl
	 <<"'CKIN(1)="<<MZp*0.6<<"         !lower invariant mass cutoff (GeV)',"<<endl
	 <<"'CKIN(2)="<<MZp*1.4<<"         !upper invariant mass cutoff (GeV)',"<<endl;
  string SMParticles[3][4]={{"Nu_e","e","u","d"},{"Nu_mu","mu","c","s"},{"Nu_tau","tau","t","b"}};
  //First generation
  for (int i=3;i>=0;i--)
    {
      PYCard <<"'PARU("<<128-i*2-1<<")="<<Vf(i)<<"   ! g_v for "<<SMParticles[0][i]<<" to Zp',"<<endl;
      PYCard <<"'PARU("<<128-i*2<<")="<<Af(i)<<"   ! g_a for "<<SMParticles[0][i]<<" to Zp',"<<endl;
    }
#if defined(WW) || defined(All)
  PYCard <<"'PARU(129)="<<theStrengthOfWWcoupling<<"    ! WW decay width ratio ',"<<endl;
#endif
  //Second generation
  for (int i=3;i>=0;i--)
    {
      PYCard <<"'PARJ("<<187-i*2-1<<")="<<Vf(i)<<"   ! g_v for "<<SMParticles[1][i]<<" to Zp',"<<endl;
      PYCard <<"'PARJ("<<187-i*2<<")="<<Af(i)<<"   ! g_a for "<<SMParticles[1][i]<<" to Zp',"<<endl;
    }
  //Third generation
  for (int i=3;i>=0;i--)
    if (i!=2) {
      PYCard <<"'PARJ("<<195-i*2-1<<")="<<Vf(i)<<"   ! g_v for "<<SMParticles[2][i]<<" to Zp',"<<endl;
      PYCard <<"'PARJ("<<195-i*2<<")="<<Af(i)<<"   ! g_a for "<<SMParticles[2][i]<<" to Zp',"<<endl;
    }
    else {
      PYCard <<"'PARJ(190)="<< ((MZp>=2*Mt)?Vf(2)*sqrt(sqrt(1-4*Mt*Mt/MZp2)*(1+2*Mt*Mt/MZp2)):0) <<"   ! g_v for t to Zp',"<<endl;
      PYCard <<"'PARJ(191)="<< ((MZp>=2*Mt)?Af(2)*sqrt(sqrt(1-4*Mt*Mt/MZp2)*(1-4*Mt*Mt/MZp2)):0) <<"   ! g_a for t to Zp',"<<endl;
    }
  PYCard <<"'MDCY(15,1)=0            ! sets tau stable - this is important (for tauola)!',"<<endl;
  //Decay channels
  for (map<int, pair<string,int> >::iterator iter=DecayChannels.begin();iter != DecayChannels.end();iter++)
    PYCard <<"'MDME("<<iter->first<<",1)="<<iter->second.second<<"    ! "<<iter->second.first<<"',"<<endl;
  PYCard << "'MSTP(44)=3              !gamma/Z/Zp: 3=> Zp only'),"<<endl;
  PYCard<<"parameterSets = cms.vstring('pythiaUESettings','processParameters')"<<endl
	<<"))"<<endl
	<<"ProductionFilterSequence = cms.Sequence(generator)"<<endl;
  cout<<CardName<<" is written."<<endl;
}
