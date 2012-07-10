#include "TFile.h"
#include "TTree.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include <fstream>

using namespace RooFit;

void WidthFit(char *filename){
  TFile * f = TFile::Open(filename);

  Float_t M_Zp=0;
  char *mass_start=strstr(filename,"_M")+2,mass_str[40];
  Byte_t mass_len=strstr(mass_start,"_")-mass_start;
  strncpy(mass_str,mass_start,mass_len);
  mass_str[mass_len]='\0';
  printf("mass=%sGeV/c2\n",mass_str);
  if ( EOF==sscanf(mass_str, "%f", &M_Zp) ) return;

  TTree * tree = (TTree*) f->Get("T");
  //  RooRealVar mes("dimuon_mass","dimuon_mass",100,300) ;
  
  char factorysource[300];
  sprintf(factorysource,"RooBreitWigner::Shape(dimuon_mass[%f,%f],mean[%f,%f,%f],width[0.01,0,1])",M_Zp-10.,M_Zp+10.,M_Zp,M_Zp-10.,M_Zp+10.);
  RooWorkspace w("w",kTRUE);
  w.factory(factorysource);

  RooDataSet data("data","dataset with mass",tree,*w.var("dimuon_mass"));
 
  // --- Perform extended ML fit of composite PDF to toy data ---
  w.pdf("Shape")->fitTo(data);
  
  // --- Plot toy data and composite PDF overlaid ---
  RooPlot* mesframe = w.var("dimuon_mass")->frame();
  data.plotOn(mesframe);
  w.pdf("Shape")->plotOn(mesframe);
  mesframe->Draw();
  
  //  char result_filename[300];
  //  strncpy(result_filename,filename,strstr(filename,"_GenShape.root")-filename);
  //  strcpy(result_filename+strlen(filename)-14,".txt\0");
  //  ofstream RESULT(result_filename);
  //  RESULT<<M_Zp<<" "<<w.var("width")->getVal()<<endl;
  cout<<M_Zp<<"\t"<<w.var("width")->getVal()<<"\t"<<w.var("width")->getError()<<endl;
}
