#define Analyzer_Nt16_cxx
//Self
#include "Analyzer_Nt16.h"
#include "TotalEvents.h"
//ROOT
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>

//C++
#include <fstream>
#include <cmath>

#define Muon_Mass 0.105658367
//Mu_Id Rule: + for Muon;- for AntiMuon; 0 means ignore the cut (pass).
using namespace std;
bool Selector(vector<char *>& Names,vector<char *>& FileNames, vector<double> &sigma);
bool Analyzer_Nt16::IsHLTObj(int Mu_Id, vector<double>& HLTObj_pt, vector<double>& HLTObj_eta,  vector<double>& HLTObj_phi)
{
  double pt,eta,phi;
  if (!Mu_Id) return true;
  if (Mu_Id>0)
    {
      pt=Muon_pt[Mu_Id-1];
      eta=Muon_eta[Mu_Id-1];
      phi=Muon_phi[Mu_Id-1];
    }
  else
    {
      pt=AntiMuon_pt[Mu_Id-1];
      eta=AntiMuon_eta[-Mu_Id-1];
      phi=AntiMuon_phi[-Mu_Id-1];
    }
  unsigned int Num_HLTObj=HLTObj_pt.size();
  //For determing the parameters
  /*
  clog<<"pt: "<<pt<<";eta: "<<eta<<";phi: "<<phi<<endl<<"-------------------------------"<<endl;
  for (unsigned int i=0; i<Num_HLTObj; i++)
    clog<<"pt: "<<HLTObj_pt[i]<<";eta: "<<HLTObj_eta[i]<<";phi: "<<HLTObj_phi[i]<<endl;
  clog<<"***************************************"<<endl;
  */
  for (unsigned int i=0; i<Num_HLTObj; i++)
    if ((abs(eta-HLTObj_eta[i])<0.003)&&(abs(phi-HLTObj_phi[i])<0.003)) return true;
  return false;
}
bool Analyzer_Nt16::IsoR03Cut(int Mu_Id, double Sumpt)
{
  if (!Mu_Id) return true;
  if (Mu_Id>0)
    if (abs(Muon_isoR03sumPt[Mu_Id-1])<Sumpt) return true;
    else return false;
  else
    if (abs(AntiMuon_isoR03sumPt[-Mu_Id-1])<Sumpt) return true;
    else return false;
}
bool Analyzer_Nt16::ptCut(int Mu_Id, double pt)
{
  if (!Mu_Id) return true;
  if (Mu_Id>0)
    if (abs(Muon_pt[Mu_Id-1])>pt) return true;
    else return false;
  else
    if (abs(AntiMuon_pt[-Mu_Id-1])>pt) return true;
    else return false;
}
bool Analyzer_Nt16::EtaCut(int Mu_Id, double eta)
{
  if (!Mu_Id) return true;
  if (Mu_Id>0)
    if (Muon_eta[Mu_Id-1]>-eta&&Muon_eta[Mu_Id-1]<eta) return true;
    else return false;
  else
    if (AntiMuon_eta[-Mu_Id-1]>-eta&&AntiMuon_eta[-Mu_Id-1]<eta) return true;
    else return false;
}
bool Analyzer_Nt16::Cut(int Mu_Id)
{
  return EtaCut(Mu_Id,2.1)&&ptCut(Mu_Id,50)&&IsoR03Cut(Mu_Id,0.05)&&(IsHLTObj(Mu_Id,*hltSingleMu15L3PreFiltered15_pt,*hltSingleMu15L3PreFiltered15_eta,*hltSingleMu15L3PreFiltered15_phi)||IsHLTObj(Mu_Id,*hltSingleMuIsoL3IsoFiltered9_pt,*hltSingleMuIsoL3IsoFiltered9_eta,*hltSingleMuIsoL3IsoFiltered9_phi));
}
double Analyzer_Nt16::d0(int Mu_Id,unsigned int PV_Id )
{
   double dX,dY;
   if (!Mu_Id) return 0;
   if (Mu_Id>0)
     {
       dX=((double) Muon_Vertex[Mu_Id-1][0])-(*offlinePrimaryVertices_X)[PV_Id];
       dY=((double) Muon_Vertex[Mu_Id-1][1])-(*offlinePrimaryVertices_Y)[PV_Id];
     }
   else
     {
       dX=((double) AntiMuon_Vertex[-Mu_Id-1][0])-(*offlinePrimaryVertices_X)[PV_Id];
       dY=((double) AntiMuon_Vertex[-Mu_Id-1][1])-(*offlinePrimaryVertices_Y)[PV_Id];
     }
   return sqrt(dX*dX+dY*dY);
}
double Analyzer_Nt16::InvarDimuonMass(int Mu_Id1, int Mu_Id2)
{
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > mu1,mu2;
   if (!Mu_Id1||!Mu_Id2) return 0;
   if (Mu_Id1>0)
      mu1=ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >(Muon_pt[Mu_Id1-1],Muon_eta[Mu_Id1-1],Muon_phi[Mu_Id1-1],Muon_Mass);
   else
      mu1=ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >(AntiMuon_pt[-Mu_Id1-1],AntiMuon_eta[-Mu_Id1-1],AntiMuon_phi[-Mu_Id1-1],Muon_Mass);
   if (Mu_Id2>0)
      mu2=ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >(Muon_pt[Mu_Id2-1],Muon_eta[Mu_Id2-1],Muon_phi[Mu_Id2-1],Muon_Mass);
   else
      mu2=ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >(AntiMuon_pt[-Mu_Id2-1],AntiMuon_eta[-Mu_Id2-1],AntiMuon_phi[-Mu_Id2-1],Muon_Mass);
   return (mu1+mu2).mass();
}
void Analyzer_Nt16::Loop()
{
   unsigned int i,Total_Muons=0,Total_AntiMuons=0,Total_Muons_Passed=0,Total_AntiMuons_Passed=0;
   Hist1D[0] = new TH1D("DM","Dimuons Mass",100, 200, 800 );
   Hist1D[1] = new TH1D("VD","Vertex D0^2",100,0,0.03);
   Hist1D[2] = new TH1D("Iso","MuonIsoR03SumPt",100, -2.001, 2.001 );
   vector <int> Passed_Mu,Passed_AntiMu;
   vector <int>::const_iterator Mu1,Mu2;
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
     {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       for (i=0;i<Muon_size;i++)
	 if (Cut(i+1)) Passed_Mu.push_back(i);
       for (i=0;i<AntiMuon_size;i++)
	 if (Cut(-i-1)) Passed_AntiMu.push_back(i);
       Total_Muons+=Muon_size;
       Total_AntiMuons+=AntiMuon_size;
       Total_Muons_Passed+=Passed_Mu.size();
       Total_AntiMuons_Passed+=Passed_AntiMu.size();

       for (Mu1=Passed_Mu.begin();Mu1!=Passed_Mu.end();Mu1++)
	 for (Mu2=Passed_AntiMu.begin();Mu2!=Passed_AntiMu.end();Mu2++)
	   Hist1D[0]->Fill(InvarDimuonMass(*Mu1+1,-*Mu2-1));
       /*
       for (i=0;i<Muon_size;i++)
	 for (j=0;j<AntiMuon_size;j++)
	     Hist1D[0]->Fill(InvarDimuonMass(i+1,-j-1));
       */
        Hist1D[1]->Fill(d0(1,0));
       for (Mu1=Passed_Mu.begin();Mu1!=Passed_Mu.end();Mu1++)
	 Hist1D[2]->Fill(Muon_isoR03sumPt[*Mu1]);
       for (Mu2=Passed_AntiMu.begin();Mu2!=Passed_AntiMu.end();Mu2++)
	 Hist1D[2]->Fill(AntiMuon_isoR03sumPt[*Mu2]);
       Passed_Mu.clear();
       Passed_AntiMu.clear();
   }
   clog<<"Total Number of Muons: "<<Total_Muons<<endl<<"Total Number of AntiMuons: "<<Total_AntiMuons<<endl;
   clog<<"Passed Muons: "<<Total_Muons_Passed<<endl<<"Passed AntiMuons: "<<Total_AntiMuons_Passed<<endl;
   clog<<"Efficiency: "<<(double) (Total_Muons_Passed+Total_AntiMuons_Passed)/(double) (Total_Muons+Total_AntiMuons)<<endl;
}

void MyMain()
{
  unsigned char ColorList[]={2,3,4,5,6,7,8,9,10,41,42,43,44,45,46};
  unsigned char numColor=sizeof(ColorList)/sizeof(unsigned char);
  vector <char *> Names, FileNames;
  vector <double> sigma;
  THStack *THS[MaxHist];
  TLegend *Legends = new TLegend(0.75,0.6,0.89,0.89);
  TTree *tree;
  Long64_t OrgTotEvents;
  double weight;
  Analyzer_Nt16 *Histo;
  TH1D *tempHist;
  THS[0] = new THStack("h1","Invariant Mass of Dimuon (Weighted)");
  THS[1] = new THStack("h2","Invariant Mass of Dimuon (Unweighted)");
  THS[2] = new THStack("h3","d0");
  THS[3] = new THStack("h4","IsoR03SumPt");
  if (!Selector(Names,FileNames,sigma)) exit(0);
  for (unsigned int i=0; i<sigma.size(); i++)
    {
      TFile f(FileNames[i]);
      if (!f.IsOpen())
	{
	  cerr<<"File "<<FileNames[i]<<" cannot be opened."<<endl;
	  continue;
	}
      f.cd("MuonAnalyzer");
      tree = (TTree*)gDirectory->Get("Summerization");
      OrgTotEvents=TotalEvents(tree);
      weight= sigma[i]/ (double) OrgTotEvents;
      clog<<Names[i]<<" : Weight=sigma/Original Total Events="<<sigma[i]<<"/"<<OrgTotEvents<<"="<<weight<<endl;
      tree = (TTree*)gDirectory->Get("Muon");
      Histo = new Analyzer_Nt16(tree);
      Histo->Loop();
      Histo->Hist1D[0]->SetFillColor(ColorList[i%numColor]);
      Histo->Hist1D[1]->SetFillColor(ColorList[i%numColor]);
      Histo->Hist1D[2]->SetFillColor(ColorList[i%numColor]);
      Legends->AddEntry(Histo->Hist1D[0],Names[i],"F");
      tempHist = new TH1D("DM","Dimuons Mass",100, 200, 800 );
      *tempHist = weight*(*Histo->Hist1D[0]);
      THS[0]->Add(tempHist);
      THS[1]->Add(Histo->Hist1D[0]);
      THS[2]->Add(Histo->Hist1D[1]);
      THS[3]->Add(Histo->Hist1D[2]);
      //Ghosts, without them, the hisograms will be empty.
      //It is to redefine the names of the histos to new handles.
      TH1D *Ghost1 = new TH1D("DM","Ghost",10, 0, 1 );
      TH1D *Ghost2 = new TH1D("VD","Ghost",10, 0, 1 );
      TH1D *Ghost3 = new TH1D("Iso","Ghost",10, 0, 1 );
    }
   TCanvas *MyC = new TCanvas("MyC","Test canvas",1);
   MyC->Divide(2,2);
   MyC->cd(1);
   THS[0]->Draw();
   Legends->Draw();
   MyC->cd(2);
   THS[1]->Draw();
   Legends->Draw();
   MyC->cd(3);
   THS[2]->Draw();
   Legends->Draw();
   MyC->cd(4);
   Legends->Draw();
   THS[3]->Draw();
}

bool Selector(vector<char *>& Names,vector<char *>& FileNames, vector<double> &sigma)
{
  unsigned int Default_Choices[]={2,3,4,5,6,7,8,9,10,11}, num_Ntuples=0, temp, i;
  fstream f("/home/zhangjin/DATALIST.csv");
  if (!f) 
    {
      cerr<<"File cannot be opened"<<endl;
      return false;
    }
  f.ignore(2048,'\n');//skip the title line
  while ( (! f.eof()) && (f.get()=='"') )
    {
      Names.push_back(new char [512]);
      f.getline(Names.back(),512,'"');//read names
      f.ignore(1);
      f.ignore(512,',');//skip DBS position
      f.ignore(512,',');//skip castor position
      f.ignore(512,',');//skip Python position
      if (f.get()!=',')
	{
	   FileNames.push_back(new char [512]);
	   f.getline(FileNames.back(),512,'"');//read Ntuples file names
	   f.ignore(1);
	   sigma.push_back(0);
	   f>>sigma.back();//read cross section
	   num_Ntuples++;
	   clog<<num_Ntuples<<" : "<<FileNames.back()<<", "<<sigma.back()<<endl;
	}
      else
	Names.pop_back();
      f.ignore(1024,'\n');//skip the all the others
    }
  f.close();
  clog<<0<<" Default Choices:";
  for (i=0; i<sizeof(Default_Choices)/sizeof(unsigned int);i++)
    clog<<Default_Choices[i]<<",";
  clog<<endl<<"If pause, push <enter>";
  cin.clear(ios::goodbit);
  cin.get();
  clog<<endl<<"Choose some: ";
  while (cin>>temp)
    {
      if (!temp) break;
      if (temp<=num_Ntuples&&temp>0)
	{
	  Names.push_back(Names[temp-1]);
	  FileNames.push_back(FileNames[temp-1]);
	  sigma.push_back(sigma[temp-1]);
	}
      else
	cerr<<"Invalid choice."<<endl;
      clog<<"Choose some: ";
    }
  if (!temp)
    for (i=0; i<sizeof(Default_Choices)/sizeof(unsigned int);i++)
      {
	  Names.push_back(Names[Default_Choices[i]-1]);
	  FileNames.push_back(FileNames[Default_Choices[i]-1]);
	  sigma.push_back(sigma[Default_Choices[i]-1]);
      }
  Names.erase(Names.begin(),Names.begin()+num_Ntuples);
  FileNames.erase(FileNames.begin(),FileNames.begin()+num_Ntuples);
  sigma.erase(sigma.begin(),sigma.begin()+num_Ntuples);
  return true;
}
