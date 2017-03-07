#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"

int draw_accep_ratio(int energy) {
    int TrackletType[3] = {12, 13, 23};

    for (int i=0; i<3; ++i) {
        TFile* finput = new TFile(Form("rootfiles/accep/%itev/100/acceptance-%i.root", energy, TrackletType[i]), "read");
        TH2F* haccep_data = (TH2F*)finput->Get("haccep_data");
        TH2F* haccep_mc = (TH2F*)finput->Get("haccep_mc");

        TH2F* hratio = (TH2F*)haccep_mc->Clone("hratio");
        hratio->Divide(haccep_data);

        TCanvas* c1 = new TCanvas(Form("c1_%i", TrackletType[i]), "", 600, 600);
        hratio->SetStats(0);
        hratio->SetAxisRange(0, 2, "Z");
        hratio->Draw("colz");
        c1->SaveAs(Form("figs/corrs/accep-ratio-%i-%itev.png", TrackletType[i], energy));

        finput->Close();
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 2)
        return draw_accep_ratio(atoi(argv[1]));
    else
        return 1;
}
