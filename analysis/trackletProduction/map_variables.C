#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"

#include <vector>
#include <string>

typedef struct twod_var_t {
   std::string title;
   std::string var_x[3];
   std::string var_y[3];
   int nbins_x;
   int nbins_y;
   float lrange_x;
   float lrange_y;
   float hrange_x;
   float hrange_y;
   std::string selection;
} twod_var_t;

std::vector<twod_var_t> options = {
   {"eta-phi", {"phi1", "phi2", "phi3"}, {"eta1", "eta2", "eta3"}, 200, 200, -3.5, -3, 3.5, 3, "(1)"}
};

std::string layer_cuts_5tev[3] = {
   "(1)", "(1)", "(1)"
};

std::string layer_cuts_8tev[3] = {
   "(fabs(phi1-1.78421)>0.00001 && fabs(eta1+1.16703)>0.00001 && fabs(phi1-1.43093)>0.00001 && fabs(eta1+0.66483)>0.00001 && fabs(phi1-1.42661)>0.00001 && fabs(eta1+0.66445)>0.00001 && fabs(phi1-1.40936)>0.00001 && fabs(eta1+0.66284)>0.00001 && fabs(phi1-1.40721)>0.00001 && fabs(eta1+0.66262)>0.00001 && fabs(phi1-2.62996)>0.00001 && fabs(eta1+2.10101)>0.00001 && fabs(phi1-1.35537)>0.00001 && fabs(eta1+2.48903)>0.00001 && fabs(phi1-1.35055)>0.00001 && fabs(eta1+2.48966)>0.00001 && fabs(phi1-1.40721)>0.00001 && fabs(eta1+0.66262)>0.00001 && fabs(phi1-1.40936)>0.00001 && fabs(eta1+0.66284)>0.00001 && fabs(phi1-1.41582)>0.00001 && fabs(eta1+0.66346)>0.00001 && fabs(phi1-1.42661)>0.00001 && fabs(eta1+0.66445)>0.00001 && fabs(phi1-1.43093)>0.00001 && fabs(eta1+0.66483)>0.00001 && fabs(phi1-3.09539)>0.00001 && fabs(eta1-2.08781)>0.00001 && fabs(phi1-1.15823)>0.00001 && fabs(eta1+2.09717)>0.00001 && fabs(phi1-1.16068)>0.00001 && fabs(eta1+2.09732)>0.00001)",
   "(fabs(phi2+2.02809)>0.00001 && fabs(eta2-0.34774)>0.00001 && fabs(phi2+2.03070)>0.00001 && fabs(eta2-0.34771)>0.00001 && fabs(phi2+2.29325)>0.00001 && fabs(eta2-0.70344)>0.00001 && fabs(phi2+2.29048)>0.00001 && fabs(eta2-0.70358)>0.00001 && fabs(phi2-1.50475)>0.00001 && fabs(eta2-1.41259)>0.00001 && fabs(phi2+1.91616)>0.00001 && fabs(eta2-1.12794)>0.00001 && fabs(phi2+0.96467)>0.00001 && fabs(eta2-0.65614)>0.00001 && fabs(phi2+0.96730)>0.00001 && fabs(eta2-0.65611)>0.00001 && fabs(phi2+0.34368)>0.00001 && fabs(eta2-1.84190)>0.00001 && fabs(phi2+2.50245)>0.00001 && fabs(eta2-0.28293)>0.00001 && fabs(phi2+2.50521)>0.00001 && fabs(eta2-0.28301)>0.00001 && fabs(phi2+0.88593)>0.00001 && fabs(eta2-0.68248)>0.00001 && fabs(phi2+0.96467)>0.00001 && fabs(eta2-0.65614)>0.00001 && fabs(phi2+0.96730)>0.00001 && fabs(eta2-0.65611)>0.00001 && fabs(phi2+2.12876)>0.00001 && fabs(eta2+0.34704)>0.00001 && fabs(phi2+2.13015)>0.00001 && fabs(eta2+0.34905)>0.00001 && fabs(phi2+2.13154)>0.00001 && fabs(eta2+0.34711)>0.00001 && fabs(phi2+0.96534)>0.00001 && fabs(eta2+0.24780)>0.00001)",
   "(fabs(phi3-0.01456)>0.00001 && fabs(eta3+0.06675)>0.00001 && fabs(phi3+0.72997)>0.00001 && fabs(eta3+0.62795)>0.00001 && fabs(phi3+1.97283)>0.00001 && fabs(eta3-0.73149)>0.00001)"
};

int map_variables(const char* fname, const char* label, int opt, int energy, int data = 0) {
   TFile* f = new TFile(fname, "r");
   TTree* t = (TTree*)f->Get("ana/PixelTree");

   std::string layer_cuts[3];
   if (data) {
      switch (energy) {
         case 5:
            std::copy(layer_cuts_5tev, layer_cuts_5tev + 3, layer_cuts);
            break;
         case 8:
            std::copy(layer_cuts_8tev, layer_cuts_8tev + 3, layer_cuts);
            break;
         default:
            printf("check sqrt(snn) energy");
            return 1;
      }
   }

   std::string layer_selections[3];
   for (std::size_t i=0; i<3; ++i) {
      if (data) {
         switch (energy) {
            case 5:
               layer_selections[i] = options[opt].selection + " && " + layer_cuts_5tev[i];
               break;
            case 8:
               layer_selections[i] = options[opt].selection + " && " + layer_cuts_8tev[i] + " && (nLumi>164)";
               break;
            default:
               printf("check sqrt(snn) energy!\n");
               return 1;
         }
      }
   }

   TH2D* h1 = new TH2D("h1", Form(";%s;%s", options[opt].var_x[0].c_str(), options[opt].var_y[0].c_str()), options[opt].nbins_x, options[opt].lrange_x, options[opt].hrange_x, options[opt].nbins_y, options[opt].lrange_y, options[opt].hrange_y);
   t->Draw(Form("%s:%s>>h1", options[opt].var_y[0].c_str(), options[opt].var_x[0].c_str()), layer_selections[0].c_str(), "colz goff");
   TH2D* h2 = new TH2D("h2", Form(";%s;%s", options[opt].var_x[1].c_str(), options[opt].var_y[1].c_str()), options[opt].nbins_x, options[opt].lrange_x, options[opt].hrange_x, options[opt].nbins_y, options[opt].lrange_y, options[opt].hrange_y);
   t->Draw(Form("%s:%s>>h2", options[opt].var_y[1].c_str(), options[opt].var_x[1].c_str()), layer_selections[1].c_str(), "colz goff");
   TH2D* h3 = new TH2D("h3", Form(";%s;%s", options[opt].var_x[2].c_str(), options[opt].var_y[2].c_str()), options[opt].nbins_x, options[opt].lrange_x, options[opt].hrange_x, options[opt].nbins_y, options[opt].lrange_y, options[opt].hrange_y);
   t->Draw(Form("%s:%s>>h3", options[opt].var_y[2].c_str(), options[opt].var_x[2].c_str()), layer_selections[2].c_str(), "colz goff");

   TCanvas* c1 = new TCanvas("c1", "", 700, 600);
   h1->SetStats(0);
   h1->Draw("colz");
   c1->SaveAs(Form("figs/%s-%s-l1.png", options[opt].title.c_str(), label));

   TCanvas* c2 = new TCanvas("c2", "", 700, 600);
   h2->SetStats(0);
   h2->Draw("colz");
   c2->SaveAs(Form("figs/%s-%s-l2.png", options[opt].title.c_str(), label));

   TCanvas* c3 = new TCanvas("c3", "", 700, 600);
   h3->SetStats(0);
   h3->Draw("colz");
   c3->SaveAs(Form("figs/%s-%s-l3.png", options[opt].title.c_str(), label));

   return 0;
}

int main(int argc, char* argv[]) {
   if (argc == 5)
      return map_variables(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
   else if (argc == 6)
      return map_variables(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
   else
      return 1;
}

