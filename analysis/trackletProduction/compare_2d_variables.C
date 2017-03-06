#include <TFile.h>
#include <TTree.h>
#include "TCut.h"
#include <TH2D.h>
#include <TCanvas.h>

#include <vector>
#include <string>

typedef struct var_t {
   std::string title;
   std::string variable[2];
   TCut selection;
   int nbins[2];
   float lrange[2];
   float hrange[2];
   std::string axis_label;
} var_t;

std::vector<var_t> options = {
   {"eta-phi", {"eta1", "phi1"}, "vz[1]>-99", {200, 200}, {-3.5, -3}, {3.5, 3}, ";#eta;#phi;"},
   {"vz-eta", {"vz[1]", "eta1"}, "vz[1]>-99", {200, 100}, {-3, -20}, {3, 20}, ";v_{z};#eta;"}
};

int TrackletType[3] = {12, 13, 23};

int compare_tracklet_variables(const char* input_file, const char* label, int opt) {
   TCut final_selection = options[opt].selection && "passHLT";
   final_selection = final_selection * "weight";

   for (int i=0; i<3; ++i) {
      TFile* finput = new TFile(input_file, "READ");
      TTree* tinput = (TTree*)finput->Get(Form("TrackletTree%i", TrackletType[i]));

      TH2D* hinput = new TH2D("hinput", options[opt].axis_label.c_str(), options[opt].nbins[0], options[opt].lrange[0], options[opt].hrange[0], options[opt].nbins[1], options[opt].lrange[1], options[opt].hrange[1]);
      tinput->Draw(Form("%s:%s>>hinput", options[opt].variable[0].c_str(), options[opt].variable[1].c_str()), final_selection, "goff");
      hinput->SetStats(0);

      TCanvas* c1 = new TCanvas(Form("c_%i", TrackletType[i]), "", 700, 600);
      hinput->Draw("colz");

      c1->SaveAs(Form("figs/tracklet-%s-%i-%s.png", options[opt].title.c_str(), TrackletType[i], label));
   }

   return 0;
}

int main(int argc, char* argv[]) {
   if (argc == 4)
      return compare_tracklet_variables(argv[1], argv[2], atoi(argv[3]));
   else
      return 1;
}
