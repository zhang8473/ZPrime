//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov  4 16:43:42 2009 by ROOT version 5.25/02
// from TTree Muon/Muon
// found on file: MyNtuple_DYmumu_Mcut200_7TeV_MC31XV8.root
//////////////////////////////////////////////////////////

#ifndef Analyzer_Nt16_h
#define Analyzer_Nt16_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

#define MaxMuon 20
#define MaxHist 9

class Analyzer_Nt16 {
public :
   
   Long64_t Loop();
   TH1D *Hist1D[MaxHist];
   TH1I *Hist1I[MaxHist];
   
private :
   //Analyzer
   bool IsHLTObj(int Mu_Id, vector<double>& HLTObj_pt, vector<double>& HLTObj_eta,  vector<double>& HLTObj_phi);
   bool IsoR03Cut(int Mu_Id, double Sumpt);
   bool ptCut(int Mu_Id, double pt);
   bool EtaCut(int Mu_Id, double eta);
   bool Cut(int Mu_Id);
   double d0(int Mu_Id,unsigned int PV_Id );
   double InvarDimuonMass(int Mu_Id1, int Mu_Id2);

   //TTree
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<double>  *hltSingleMu15L3PreFiltered15_pt;
   vector<double>  *hltSingleMu15L3PreFiltered15_eta;
   vector<double>  *hltSingleMu15L3PreFiltered15_phi;
   vector<double>  *hltSingleMuIsoL3IsoFiltered9_pt;
   vector<double>  *hltSingleMuIsoL3IsoFiltered9_eta;
   vector<double>  *hltSingleMuIsoL3IsoFiltered9_phi;
   vector<double>  *offlinePrimaryVertices_X;
   vector<double>  *offlinePrimaryVertices_eX;
   vector<double>  *offlinePrimaryVertices_Y;
   vector<double>  *offlinePrimaryVertices_eY;
   vector<double>  *offlinePrimaryVertices_Z;
   vector<double>  *offlinePrimaryVertices_eZ;
   vector<double>  *offlinePrimaryVerticesWithBS_X;
   vector<double>  *offlinePrimaryVerticesWithBS_eX;
   vector<double>  *offlinePrimaryVerticesWithBS_Y;
   vector<double>  *offlinePrimaryVerticesWithBS_eY;
   vector<double>  *offlinePrimaryVerticesWithBS_Z;
   vector<double>  *offlinePrimaryVerticesWithBS_eZ;
   vector<float>   *GenMuon_pt;
   vector<float>   *GenMuon_eta;
   vector<float>   *GenMuon_phi;
   vector<float>   *GenAntiMuon_pt;
   vector<float>   *GenAntiMuon_eta;
   vector<float>   *GenAntiMuon_phi;
   UInt_t          Muon_size;
   Float_t         Muon_pt[MaxMuon];   //[Mnum]
   Float_t         Muon_eta[MaxMuon];   //[Mnum]
   Float_t         Muon_phi[MaxMuon];   //[Mnum]
   Float_t         Muon_Vertex[MaxMuon][3];   //[Mnum]
   Float_t         Muon_isoR03sumPt[MaxMuon];   //[Mnum]
   Float_t         Muon_isoR03emEt[MaxMuon];   //[Mnum]
   Float_t         Muon_isoR03hadEt[MaxMuon];   //[Mnum]
   Float_t         Muon_isoR03hoEt[MaxMuon];   //[Mnum]
   Float_t         Muon_isoR03nJets[MaxMuon];   //[Mnum]
   Float_t         Muon_isoR03nTracks[MaxMuon];   //[Mnum]
   Float_t         Muon_isoR05sumPt[MaxMuon];   //[Mnum]
   Float_t         Muon_isoR05emEt[MaxMuon];   //[Mnum]
   Float_t         Muon_isoR05hadEt[MaxMuon];   //[Mnum]
   Float_t         Muon_isoR05hoEt[MaxMuon];   //[Mnum]
   Float_t         Muon_isoR05nJets[MaxMuon];   //[Mnum]
   Float_t         Muon_isoR05nTracks[MaxMuon];   //[Mnum]
   Float_t         Muon_isoemVetoEt[MaxMuon];   //[Mnum]
   Float_t         Muon_isohadVetoEt[MaxMuon];   //[Mnum]
   Float_t         Muon_isohoVetoEt[MaxMuon];   //[Mnum]
   UInt_t          AntiMuon_size;
   Float_t         AntiMuon_pt[MaxMuon];   //[Anum]
   Float_t         AntiMuon_eta[MaxMuon];   //[Anum]
   Float_t         AntiMuon_phi[MaxMuon];   //[Anum]
   Float_t         AntiMuon_Vertex[MaxMuon][3];   //[Anum]
   Float_t         AntiMuon_isoR03sumPt[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isoR03emEt[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isoR03hadEt[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isoR03hoEt[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isoR03nJets[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isoR03nTracks[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isoR05sumPt[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isoR05emEt[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isoR05hadEt[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isoR05hoEt[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isoR05nJets[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isoR05nTracks[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isoemVetoEt[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isohadVetoEt[MaxMuon];   //[Anum]
   Float_t         AntiMuon_isohoVetoEt[MaxMuon];   //[Anum]
   Bool_t          HLTAcceptance_generation_step;
   Bool_t          HLTAcceptance_simulation_step;
   Bool_t          HLTAcceptance_digitisation_step;
   Bool_t          HLTAcceptance_L1simulation_step;
   Bool_t          HLTAcceptance_digi2raw_step;
   Bool_t          HLTAcceptance_HLTriggerFirstPath;
   Bool_t          HLTAcceptance_HLT_L1Jet15;
   Bool_t          HLTAcceptance_HLT_Jet30;
   Bool_t          HLTAcceptance_HLT_Jet50;
   Bool_t          HLTAcceptance_HLT_Jet80;
   Bool_t          HLTAcceptance_HLT_Jet110;
   Bool_t          HLTAcceptance_HLT_Jet140;
   Bool_t          HLTAcceptance_HLT_Jet180;
   Bool_t          HLTAcceptance_HLT_FwdJet40;
   Bool_t          HLTAcceptance_HLT_DiJetAve15U_1E31;
   Bool_t          HLTAcceptance_HLT_DiJetAve30U_1E31;
   Bool_t          HLTAcceptance_HLT_DiJetAve50U;
   Bool_t          HLTAcceptance_HLT_DiJetAve70U;
   Bool_t          HLTAcceptance_HLT_DiJetAve130U;
   Bool_t          HLTAcceptance_HLT_QuadJet30;
   Bool_t          HLTAcceptance_HLT_SumET120;
   Bool_t          HLTAcceptance_HLT_L1MET20;
   Bool_t          HLTAcceptance_HLT_MET35;
   Bool_t          HLTAcceptance_HLT_MET60;
   Bool_t          HLTAcceptance_HLT_MET100;
   Bool_t          HLTAcceptance_HLT_HT200;
   Bool_t          HLTAcceptance_HLT_HT300_MHT100;
   Bool_t          HLTAcceptance_HLT_L1MuOpen;
   Bool_t          HLTAcceptance_HLT_L1Mu;
   Bool_t          HLTAcceptance_HLT_L1Mu20HQ;
   Bool_t          HLTAcceptance_HLT_L1Mu30;
   Bool_t          HLTAcceptance_HLT_L2Mu11;
   Bool_t          HLTAcceptance_HLT_IsoMu9;
   Bool_t          HLTAcceptance_HLT_Mu5;
   Bool_t          HLTAcceptance_HLT_Mu9;
   Bool_t          HLTAcceptance_HLT_Mu11;
   Bool_t          HLTAcceptance_HLT_Mu15;
   Bool_t          HLTAcceptance_HLT_L1DoubleMuOpen;
   Bool_t          HLTAcceptance_HLT_DoubleMu0;
   Bool_t          HLTAcceptance_HLT_DoubleMu3;
   Bool_t          HLTAcceptance_HLT_L1SingleEG5;
   Bool_t          HLTAcceptance_HLT_Ele10_SW_L1R;
   Bool_t          HLTAcceptance_HLT_Ele15_SW_L1R;
   Bool_t          HLTAcceptance_HLT_Ele15_SW_EleId_L1R;
   Bool_t          HLTAcceptance_HLT_Ele15_SW_LooseTrackIso_L1R;
   Bool_t          HLTAcceptance_HLT_Ele15_SiStrip_L1R;
   Bool_t          HLTAcceptance_HLT_Ele15_SC15_SW_LooseTrackIso_L1R;
   Bool_t          HLTAcceptance_HLT_Ele15_SC15_SW_EleId_L1R;
   Bool_t          HLTAcceptance_HLT_Ele20_SW_L1R;
   Bool_t          HLTAcceptance_HLT_Ele20_SiStrip_L1R;
   Bool_t          HLTAcceptance_HLT_Ele20_SC15_SW_L1R;
   Bool_t          HLTAcceptance_HLT_Ele25_SW_L1R;
   Bool_t          HLTAcceptance_HLT_Ele25_SW_EleId_LooseTrackIso_L1R;
   Bool_t          HLTAcceptance_HLT_DoubleEle5_SW_Jpsi_L1R;
   Bool_t          HLTAcceptance_HLT_DoubleEle5_SW_Upsilon_L1R;
   Bool_t          HLTAcceptance_HLT_DoubleEle10_SW_L1R;
   Bool_t          HLTAcceptance_HLT_Photon10_L1R;
   Bool_t          HLTAcceptance_HLT_Photon10_LooseEcalIso_TrackIso_L1R;
   Bool_t          HLTAcceptance_HLT_Photon15_L1R;
   Bool_t          HLTAcceptance_HLT_Photon20_LooseEcalIso_TrackIso_L1R;
   Bool_t          HLTAcceptance_HLT_Photon25_L1R;
   Bool_t          HLTAcceptance_HLT_Photon25_LooseEcalIso_TrackIso_L1R;
   Bool_t          HLTAcceptance_HLT_Photon30_L1R_1E31;
   Bool_t          HLTAcceptance_HLT_Photon70_L1R;
   Bool_t          HLTAcceptance_HLT_DoublePhoton10_L1R;
   Bool_t          HLTAcceptance_HLT_DoublePhoton15_L1R;
   Bool_t          HLTAcceptance_HLT_DoublePhoton15_VeryLooseEcalIso_L1R;
   Bool_t          HLTAcceptance_HLT_SingleIsoTau30_Trk5;
   Bool_t          HLTAcceptance_HLT_DoubleLooseIsoTau15_Trk5;
   Bool_t          HLTAcceptance_HLT_BTagIP_Jet80;
   Bool_t          HLTAcceptance_HLT_BTagIP_Jet120;
   Bool_t          HLTAcceptance_HLT_BTagMu_Jet20;
   Bool_t          HLTAcceptance_HLT_StoppedHSCP_1E31;
   Bool_t          HLTAcceptance_HLT_L1Mu14_L1SingleEG10;
   Bool_t          HLTAcceptance_HLT_L1Mu14_L1SingleJet15;
   Bool_t          HLTAcceptance_HLT_L1Mu14_L1ETM40;
   Bool_t          HLTAcceptance_HLT_L2Mu5_Photon9_L1R;
   Bool_t          HLTAcceptance_HLT_L2Mu9_DiJet30;
   Bool_t          HLTAcceptance_HLT_L2Mu8_HT50;
   Bool_t          HLTAcceptance_HLT_Ele10_SW_L1R_TripleJet30;
   Bool_t          HLTAcceptance_HLT_Ele10_LW_L1R_HT200;
   Bool_t          HLTAcceptance_HLT_ZeroBias;
   Bool_t          HLTAcceptance_HLT_MinBiasHcal;
   Bool_t          HLTAcceptance_HLT_MinBiasEcal;
   Bool_t          HLTAcceptance_HLT_MinBiasPixel;
   Bool_t          HLTAcceptance_HLT_MinBiasPixel_Trk5;
   Bool_t          HLTAcceptance_HLT_CSCBeamHalo;
   Bool_t          HLTAcceptance_HLT_CSCBeamHaloOverlapRing1;
   Bool_t          HLTAcceptance_HLT_CSCBeamHaloOverlapRing2;
   Bool_t          HLTAcceptance_HLT_CSCBeamHaloRing2or3;
   Bool_t          HLTAcceptance_HLT_BackwardBSC;
   Bool_t          HLTAcceptance_HLT_ForwardBSC;
   Bool_t          HLTAcceptance_HLT_TrackerCosmics;
   Bool_t          HLTAcceptance_HLT_IsoTrack_1E31;
   Bool_t          HLTAcceptance_AlCa_HcalPhiSym;
   Bool_t          HLTAcceptance_AlCa_EcalPhiSym;
   Bool_t          HLTAcceptance_AlCa_EcalPi0_1E31;
   Bool_t          HLTAcceptance_AlCa_EcalEta_1E31;
   Bool_t          HLTAcceptance_AlCa_RPCMuonNoHits;
   Bool_t          HLTAcceptance_AlCa_RPCMuonNormalisation;
   Bool_t          HLTAcceptance_HLTriggerFinalPath;
   Bool_t          HLTAcceptance_endjob_step;

   // List of branches
   TBranch        *b_hltSingleMu15L3PreFiltered15_pt;   //!
   TBranch        *b_hltSingleMu15L3PreFiltered15_eta;   //!
   TBranch        *b_hltSingleMu15L3PreFiltered15_phi;   //!
   TBranch        *b_hltSingleMuIsoL3IsoFiltered9_pt;   //!
   TBranch        *b_hltSingleMuIsoL3IsoFiltered9_eta;   //!
   TBranch        *b_hltSingleMuIsoL3IsoFiltered9_phi;   //!
   TBranch        *b_offlinePrimaryVertices_X;   //!
   TBranch        *b_offlinePrimaryVertices_eX;   //!
   TBranch        *b_offlinePrimaryVertices_Y;   //!
   TBranch        *b_offlinePrimaryVertices_eY;   //!
   TBranch        *b_offlinePrimaryVertices_Z;   //!
   TBranch        *b_offlinePrimaryVertices_eZ;   //!
   TBranch        *b_offlinePrimaryVerticesWithBS_X;   //!
   TBranch        *b_offlinePrimaryVerticesWithBS_eX;   //!
   TBranch        *b_offlinePrimaryVerticesWithBS_Y;   //!
   TBranch        *b_offlinePrimaryVerticesWithBS_eY;   //!
   TBranch        *b_offlinePrimaryVerticesWithBS_Z;   //!
   TBranch        *b_offlinePrimaryVerticesWithBS_eZ;   //!
   TBranch        *b_GenMuon_pt;   //!
   TBranch        *b_GenMuon_eta;   //!
   TBranch        *b_GenMuon_phi;   //!
   TBranch        *b_GenAntiMuon_pt;   //!
   TBranch        *b_GenAntiMuon_eta;   //!
   TBranch        *b_GenAntiMuon_phi;   //!
   TBranch        *b_Mnum;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_Vertex;   //!
   TBranch        *b_Muon_isoR03sumPt;   //!
   TBranch        *b_Muon_isoR03emEt;   //!
   TBranch        *b_Muon_isoR03hadEt;   //!
   TBranch        *b_Muon_isoR03hoEt;   //!
   TBranch        *b_Muon_isoR03nJets;   //!
   TBranch        *b_Muon_isoR03nTracks;   //!
   TBranch        *b_Muon_isoR05sumPt;   //!
   TBranch        *b_Muon_isoR05emEt;   //!
   TBranch        *b_Muon_isoR05hadEt;   //!
   TBranch        *b_Muon_isoR05hoEt;   //!
   TBranch        *b_Muon_isoR05nJets;   //!
   TBranch        *b_Muon_isoR05nTracks;   //!
   TBranch        *b_Muon_isoemVetoEt;   //!
   TBranch        *b_Muon_isohadVetoEt;   //!
   TBranch        *b_Muon_isohoVetoEt;   //!
   TBranch        *b_Anum;   //!
   TBranch        *b_AntiMuon_pt;   //!
   TBranch        *b_AntiMuon_eta;   //!
   TBranch        *b_AntiMuon_phi;   //!
   TBranch        *b_AntiMuon_Vertex;   //!
   TBranch        *b_AntiMuon_isoR03sumPt;   //!
   TBranch        *b_AntiMuon_isoR03emEt;   //!
   TBranch        *b_AntiMuon_isoR03hadEt;   //!
   TBranch        *b_AntiMuon_isoR03hoEt;   //!
   TBranch        *b_AntiMuon_isoR03nJets;   //!
   TBranch        *b_AntiMuon_isoR03nTracks;   //!
   TBranch        *b_AntiMuon_isoR05sumPt;   //!
   TBranch        *b_AntiMuon_isoR05emEt;   //!
   TBranch        *b_AntiMuon_isoR05hadEt;   //!
   TBranch        *b_AntiMuon_isoR05hoEt;   //!
   TBranch        *b_AntiMuon_isoR05nJets;   //!
   TBranch        *b_AntiMuon_isoR05nTracks;   //!
   TBranch        *b_AntiMuon_isoemVetoEt;   //!
   TBranch        *b_AntiMuon_isohadVetoEt;   //!
   TBranch        *b_AntiMuon_isohoVetoEt;   //!
   TBranch        *b_HLTAcceptance;   //!

   Analyzer_Nt16(TTree *tree=0);
   virtual ~Analyzer_Nt16();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analyzer_Nt16_cxx
Analyzer_Nt16::Analyzer_Nt16(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MyNtuple_DYmumu_Mcut200_7TeV_MC31XV8.root");
      if (!f) {
         f = new TFile("MyNtuple_DYmumu_Mcut200_7TeV_MC31XV8.root");
         f->cd("MyNtuple_DYmumu_Mcut200_7TeV_MC31XV8.root:/MuonAnalyzer");
      }
      tree = (TTree*)gDirectory->Get("Muon");

   }
   Init(tree);
}

Analyzer_Nt16::~Analyzer_Nt16()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analyzer_Nt16::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analyzer_Nt16::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analyzer_Nt16::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   hltSingleMu15L3PreFiltered15_pt = 0;
   hltSingleMu15L3PreFiltered15_eta = 0;
   hltSingleMu15L3PreFiltered15_phi = 0;
   hltSingleMuIsoL3IsoFiltered9_pt = 0;
   hltSingleMuIsoL3IsoFiltered9_eta = 0;
   hltSingleMuIsoL3IsoFiltered9_phi = 0;
   offlinePrimaryVertices_X = 0;
   offlinePrimaryVertices_eX = 0;
   offlinePrimaryVertices_Y = 0;
   offlinePrimaryVertices_eY = 0;
   offlinePrimaryVertices_Z = 0;
   offlinePrimaryVertices_eZ = 0;
   offlinePrimaryVerticesWithBS_X = 0;
   offlinePrimaryVerticesWithBS_eX = 0;
   offlinePrimaryVerticesWithBS_Y = 0;
   offlinePrimaryVerticesWithBS_eY = 0;
   offlinePrimaryVerticesWithBS_Z = 0;
   offlinePrimaryVerticesWithBS_eZ = 0;
   GenMuon_pt = 0;
   GenMuon_eta = 0;
   GenMuon_phi = 0;
   GenAntiMuon_pt = 0;
   GenAntiMuon_eta = 0;
   GenAntiMuon_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("hltSingleMu15L3PreFiltered15_pt", &hltSingleMu15L3PreFiltered15_pt, &b_hltSingleMu15L3PreFiltered15_pt);
   fChain->SetBranchAddress("hltSingleMu15L3PreFiltered15_eta", &hltSingleMu15L3PreFiltered15_eta, &b_hltSingleMu15L3PreFiltered15_eta);
   fChain->SetBranchAddress("hltSingleMu15L3PreFiltered15_phi", &hltSingleMu15L3PreFiltered15_phi, &b_hltSingleMu15L3PreFiltered15_phi);
   fChain->SetBranchAddress("hltSingleMuIsoL3IsoFiltered9_pt", &hltSingleMuIsoL3IsoFiltered9_pt, &b_hltSingleMuIsoL3IsoFiltered9_pt);
   fChain->SetBranchAddress("hltSingleMuIsoL3IsoFiltered9_eta", &hltSingleMuIsoL3IsoFiltered9_eta, &b_hltSingleMuIsoL3IsoFiltered9_eta);
   fChain->SetBranchAddress("hltSingleMuIsoL3IsoFiltered9_phi", &hltSingleMuIsoL3IsoFiltered9_phi, &b_hltSingleMuIsoL3IsoFiltered9_phi);
   fChain->SetBranchAddress("offlinePrimaryVertices_X", &offlinePrimaryVertices_X, &b_offlinePrimaryVertices_X);
   fChain->SetBranchAddress("offlinePrimaryVertices_eX", &offlinePrimaryVertices_eX, &b_offlinePrimaryVertices_eX);
   fChain->SetBranchAddress("offlinePrimaryVertices_Y", &offlinePrimaryVertices_Y, &b_offlinePrimaryVertices_Y);
   fChain->SetBranchAddress("offlinePrimaryVertices_eY", &offlinePrimaryVertices_eY, &b_offlinePrimaryVertices_eY);
   fChain->SetBranchAddress("offlinePrimaryVertices_Z", &offlinePrimaryVertices_Z, &b_offlinePrimaryVertices_Z);
   fChain->SetBranchAddress("offlinePrimaryVertices_eZ", &offlinePrimaryVertices_eZ, &b_offlinePrimaryVertices_eZ);
   fChain->SetBranchAddress("offlinePrimaryVerticesWithBS_X", &offlinePrimaryVerticesWithBS_X, &b_offlinePrimaryVerticesWithBS_X);
   fChain->SetBranchAddress("offlinePrimaryVerticesWithBS_eX", &offlinePrimaryVerticesWithBS_eX, &b_offlinePrimaryVerticesWithBS_eX);
   fChain->SetBranchAddress("offlinePrimaryVerticesWithBS_Y", &offlinePrimaryVerticesWithBS_Y, &b_offlinePrimaryVerticesWithBS_Y);
   fChain->SetBranchAddress("offlinePrimaryVerticesWithBS_eY", &offlinePrimaryVerticesWithBS_eY, &b_offlinePrimaryVerticesWithBS_eY);
   fChain->SetBranchAddress("offlinePrimaryVerticesWithBS_Z", &offlinePrimaryVerticesWithBS_Z, &b_offlinePrimaryVerticesWithBS_Z);
   fChain->SetBranchAddress("offlinePrimaryVerticesWithBS_eZ", &offlinePrimaryVerticesWithBS_eZ, &b_offlinePrimaryVerticesWithBS_eZ);
   fChain->SetBranchAddress("GenMuon_pt", &GenMuon_pt, &b_GenMuon_pt);
   fChain->SetBranchAddress("GenMuon_eta", &GenMuon_eta, &b_GenMuon_eta);
   fChain->SetBranchAddress("GenMuon_phi", &GenMuon_phi, &b_GenMuon_phi);
   fChain->SetBranchAddress("GenAntiMuon_pt", &GenAntiMuon_pt, &b_GenAntiMuon_pt);
   fChain->SetBranchAddress("GenAntiMuon_eta", &GenAntiMuon_eta, &b_GenAntiMuon_eta);
   fChain->SetBranchAddress("GenAntiMuon_phi", &GenAntiMuon_phi, &b_GenAntiMuon_phi);
   fChain->SetBranchAddress("Muon_size", &Muon_size, &b_Mnum);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_Vertex", Muon_Vertex, &b_Muon_Vertex);
   fChain->SetBranchAddress("Muon_isoR03sumPt", Muon_isoR03sumPt, &b_Muon_isoR03sumPt);
   fChain->SetBranchAddress("Muon_isoR03emEt", Muon_isoR03emEt, &b_Muon_isoR03emEt);
   fChain->SetBranchAddress("Muon_isoR03hadEt", Muon_isoR03hadEt, &b_Muon_isoR03hadEt);
   fChain->SetBranchAddress("Muon_isoR03hoEt", Muon_isoR03hoEt, &b_Muon_isoR03hoEt);
   fChain->SetBranchAddress("Muon_isoR03nJets", Muon_isoR03nJets, &b_Muon_isoR03nJets);
   fChain->SetBranchAddress("Muon_isoR03nTracks", Muon_isoR03nTracks, &b_Muon_isoR03nTracks);
   fChain->SetBranchAddress("Muon_isoR05sumPt", Muon_isoR05sumPt, &b_Muon_isoR05sumPt);
   fChain->SetBranchAddress("Muon_isoR05emEt", Muon_isoR05emEt, &b_Muon_isoR05emEt);
   fChain->SetBranchAddress("Muon_isoR05hadEt", Muon_isoR05hadEt, &b_Muon_isoR05hadEt);
   fChain->SetBranchAddress("Muon_isoR05hoEt", Muon_isoR05hoEt, &b_Muon_isoR05hoEt);
   fChain->SetBranchAddress("Muon_isoR05nJets", Muon_isoR05nJets, &b_Muon_isoR05nJets);
   fChain->SetBranchAddress("Muon_isoR05nTracks", Muon_isoR05nTracks, &b_Muon_isoR05nTracks);
   fChain->SetBranchAddress("Muon_isoemVetoEt", Muon_isoemVetoEt, &b_Muon_isoemVetoEt);
   fChain->SetBranchAddress("Muon_isohadVetoEt", Muon_isohadVetoEt, &b_Muon_isohadVetoEt);
   fChain->SetBranchAddress("Muon_isohoVetoEt", Muon_isohoVetoEt, &b_Muon_isohoVetoEt);
   fChain->SetBranchAddress("AntiMuon_size", &AntiMuon_size, &b_Anum);
   fChain->SetBranchAddress("AntiMuon_pt", AntiMuon_pt, &b_AntiMuon_pt);
   fChain->SetBranchAddress("AntiMuon_eta", AntiMuon_eta, &b_AntiMuon_eta);
   fChain->SetBranchAddress("AntiMuon_phi", AntiMuon_phi, &b_AntiMuon_phi);
   fChain->SetBranchAddress("AntiMuon_Vertex", AntiMuon_Vertex, &b_AntiMuon_Vertex);
   fChain->SetBranchAddress("AntiMuon_isoR03sumPt", AntiMuon_isoR03sumPt, &b_AntiMuon_isoR03sumPt);
   fChain->SetBranchAddress("AntiMuon_isoR03emEt", AntiMuon_isoR03emEt, &b_AntiMuon_isoR03emEt);
   fChain->SetBranchAddress("AntiMuon_isoR03hadEt", AntiMuon_isoR03hadEt, &b_AntiMuon_isoR03hadEt);
   fChain->SetBranchAddress("AntiMuon_isoR03hoEt", AntiMuon_isoR03hoEt, &b_AntiMuon_isoR03hoEt);
   fChain->SetBranchAddress("AntiMuon_isoR03nJets", AntiMuon_isoR03nJets, &b_AntiMuon_isoR03nJets);
   fChain->SetBranchAddress("AntiMuon_isoR03nTracks", AntiMuon_isoR03nTracks, &b_AntiMuon_isoR03nTracks);
   fChain->SetBranchAddress("AntiMuon_isoR05sumPt", AntiMuon_isoR05sumPt, &b_AntiMuon_isoR05sumPt);
   fChain->SetBranchAddress("AntiMuon_isoR05emEt", AntiMuon_isoR05emEt, &b_AntiMuon_isoR05emEt);
   fChain->SetBranchAddress("AntiMuon_isoR05hadEt", AntiMuon_isoR05hadEt, &b_AntiMuon_isoR05hadEt);
   fChain->SetBranchAddress("AntiMuon_isoR05hoEt", AntiMuon_isoR05hoEt, &b_AntiMuon_isoR05hoEt);
   fChain->SetBranchAddress("AntiMuon_isoR05nJets", AntiMuon_isoR05nJets, &b_AntiMuon_isoR05nJets);
   fChain->SetBranchAddress("AntiMuon_isoR05nTracks", AntiMuon_isoR05nTracks, &b_AntiMuon_isoR05nTracks);
   fChain->SetBranchAddress("AntiMuon_isoemVetoEt", AntiMuon_isoemVetoEt, &b_AntiMuon_isoemVetoEt);
   fChain->SetBranchAddress("AntiMuon_isohadVetoEt", AntiMuon_isohadVetoEt, &b_AntiMuon_isohadVetoEt);
   fChain->SetBranchAddress("AntiMuon_isohoVetoEt", AntiMuon_isohoVetoEt, &b_AntiMuon_isohoVetoEt);
   fChain->SetBranchAddress("HLTAcceptance", &HLTAcceptance_generation_step, &b_HLTAcceptance);
   Notify();
}

Bool_t Analyzer_Nt16::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analyzer_Nt16::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
#endif // #ifdef Analyzer_Nt16_cxx
