#include <TFile.h>
#include <TTree.h>
#include "TCut.h"
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <vector>
#include <string>
#include <fstream>

typedef struct var_t {
   std::string title;
   std::string variable;
   TCut selection;
   int nbins;
   float lrange;
   float hrange;
   std::string label;
   int log_scale;
} var_t;

std::vector<var_t> options = {
   {"deta", "deta", "abs(deta)<0.1", 100, -0.1, 0.1, ";#Delta#eta;", 0},
   {"deta-wide", "deta", "abs(deta)<2", 100, -2, 2, ";#Delta#eta;", 1},
   {"dphi", "dphi", "abs(dphi)<0.1", 100, -0.1, 0.1, ";#Delta#phi;", 0},
   {"dphi-wide", "dphi", "abs(dphi)<2", 100, -2, 2, ";#Delta#phi;", 1},
   {"vz", "vz[1]", "vz[1]>-99", 100, -20, 20, ";v_{z};", 0},
   {"mult", "mult", "1", 100, 0, 200, ";Background-subtracted tracklets;", 0},
   {"ntracklets", "nTracklets", "1", 100, 0, 400, ";Raw tracklets;", 0},
   {"nhit1", "nhit1", "1", 100, 0, 400, ";First layer hits;", 0}
};

std::string data_legend;
std::string legends[2] = {"EPOS LHC ", "HIJING "};

int TrackletType[3] = {12, 13, 23};

int compare_variables(const char* data_file, std::string mc_list, int opt, int energy) {
   switch (energy) {
      case 5:
         data_legend = "Run 285090 Express";
         legends[0] += "5.02 TeV";
         break;
      case 8:
         data_legend = "Run 285832 Express";
         legends[0] += "8.16 TeV";
         legends[1] += "8.16 TeV";
         break;
      default:
         return 1;
   }

   std::vector<std::string> file_list;
   if (mc_list.substr(mc_list.find_last_of(".") + 1) == "root") {
      file_list.push_back(mc_list);
   } else {
      std::ifstream file_stream(mc_list);
      if (!file_stream) return 1;
      std::string line;
      while (std::getline(file_stream, line))
         file_list.push_back(line);
   }

   std::size_t nfiles = file_list.size();

   TCut final_selection = options[opt].selection && "passHLT";
   if (energy == 8)
      final_selection = final_selection && "nLumi>164";
   final_selection = final_selection * "weight";

   for (int i=0; i<3; ++i) {
      TFile* fmc[nfiles] = {0};
      TTree* tmc[nfiles] = {0};
      TH1D* hmc[nfiles] = {0};
      for (std::size_t j=0; j<nfiles; ++j) {
         fmc[j] = new TFile(file_list[j].c_str(), "READ");
         tmc[j] = (TTree*)fmc[j]->Get(Form("TrackletTree%i", TrackletType[i]));
         hmc[j] = new TH1D(Form("hmc_%zu", j), options[opt].label.c_str(), options[opt].nbins, options[opt].lrange, options[opt].hrange);
         tmc[j]->Draw(Form("%s>>hmc_%zu", options[opt].variable.c_str(), j), final_selection, "goff");

         hmc[j]->SetStats(0);
         hmc[j]->SetLineColor(4+j);
      }

      TFile* fdata = new TFile(data_file, "READ");
      TTree* tdata = (TTree*)fdata->Get(Form("TrackletTree%i", TrackletType[i]));
      TH1D* hdata = new TH1D("hdata", options[opt].label.c_str(), options[opt].nbins, options[opt].lrange, options[opt].hrange);
      tdata->Draw(Form("%s>>hdata", options[opt].variable.c_str()), final_selection, "goff");

      hdata->SetStats(0);
      hdata->SetMarkerStyle(21);

      TCanvas* c1 = new TCanvas(Form("c_%i", TrackletType[i]), "", 600, 600);
      if (options[opt].log_scale) c1->SetLogy();

      float max_y = hdata->GetMaximum();
      for (std::size_t j=0; j<nfiles; ++j)
         if (hmc[j]->GetMaximum() > max_y)
            max_y = hmc[j]->GetMaximum();

      float min_y = 0;
      if (options[opt].log_scale) {
         min_y = hdata->GetMinimum();
         for (std::size_t j=0; j<nfiles; ++j)
            if (hmc[j]->GetMinimum() < min_y)
               min_y = hmc[j]->GetMinimum();
         min_y /= 2;
      }

      hdata->SetAxisRange(min_y, max_y * 1.2 , "Y");

      hdata->DrawNormalized("p");
      for (std::size_t j=0; j<nfiles; ++j)
         hmc[j]->DrawNormalized("hist e same");
      hdata->DrawNormalized("p same");

      TLegend* l1 = new TLegend(0.54, 0.56, 0.9, 0.7);
      l1->SetBorderSize(0);
      for (std::size_t j=0; j<nfiles; ++j)
         l1->AddEntry(hmc[j], legends[j].c_str());
      l1->AddEntry(hdata, data_legend.c_str(), "p");
      l1->Draw();

      c1->SaveAs(Form("figs/%s-%i-%itev.png", options[opt].title.c_str(), TrackletType[i], energy));
   }

   return 0;
}

int main(int argc, char* argv[]) {
   if (argc == 5)
      return compare_variables(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
   else
      return 1;
}
