#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>

void plotDphiForDifferentCut(char *infname)
{
   TFile *inf = new TFile(infname);
   TTree *t = (TTree*) inf->FindObjectAny("TrackletTree12");

   TH1F *hDeta = new TH1F*("hDeta","",0,1000);

}
