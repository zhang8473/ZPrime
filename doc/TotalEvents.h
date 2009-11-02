#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

Long64_t TotalEvents(TTree *tree)
{
  TTree   *fChain;//!pointer to the analyzed TTree or TChain
  Int_t   fCurrent; //!current Tree number in a TChain
  UInt_t  Total_Events; // Declaration of leaf types
  TBranch *b_TotalEvents; // List of branches
  if (tree == 0) 
     {
       cerr<<"Failed to open total events tree"<<endl;
	return 0;
     }
  //Init tree
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("Total_Events", &Total_Events, &b_TotalEvents);
  //Loop tree
  if (fChain == 0) return 0;
     Long64_t nentries = fChain->GetEntriesFast();
     Long64_t nbytes = 0, nb = 0, Events_Total = 0;
     for (Long64_t jentry=0; jentry<nentries;jentry++) 
       {
	 //Load Tree
	 if (!fChain) break;
	 else
	   {
	     Long64_t ientry = fChain->LoadTree(jentry);
	     if (ientry < 0) break;
	     if (fChain->InheritsFrom(TChain::Class()))
	       {
		 TChain *chain = (TChain*)fChain;
		 if (chain->GetTreeNumber() != fCurrent) 
		   fCurrent = chain->GetTreeNumber();
	       }
	     nb = fChain->GetEntry(jentry);nbytes += nb;
	     Events_Total+=Total_Events;
	   }
       }
  return Events_Total;
}
