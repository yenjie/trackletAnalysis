#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"

void draw_legend(TH1F* h1, const char* label);

int draw_corr_plots(const char* label) {
    int TrackletType[3] = {12, 13, 23};

    for (int i=0; i<3; ++i) {
        TFile* finput = new TFile(Form("correction/correction-%i-%s.root", TrackletType[i], label), "read");

        TH1F* hxi = (TH1F*)((TH1F*)finput->Get("hTrigEff"))->Clone("hxi");
        TH1F* hsdfrac = (TH1F*)((TH1F*)finput->Get("hSDFrac"))->Clone("hsdfrac");
        TH1F* hempty = (TH1F*)((TH1F*)finput->Get("hEmptyEvtCorrection"))->Clone("hempty");

        TCanvas* c1 = new TCanvas(Form("c1_%i", TrackletType[i]), "", 600, 600);
        hxi->SetStats(0);
        hxi->SetAxisRange(0, 1.2, "Y");
        hxi->SetTitle("xi;multiplicity;#xi");
        hxi->Draw();
        draw_legend(hxi, label);
        c1->SaveAs(Form("figs/corrs/xi-%s-%i.png", label, TrackletType[i]));

        TCanvas* c2 = new TCanvas(Form("c2_%i", TrackletType[i]), "", 600, 600);
        hsdfrac->SetStats(0);
        hsdfrac->SetAxisRange(-0.01, 0.25, "Y");
        hsdfrac->SetTitle("SD fraction;multiplicity;SD_{frac}");
        hsdfrac->Draw();
        draw_legend(hsdfrac, label);
        c2->SaveAs(Form("figs/corrs/sdfrac-%s-%i.png", label, TrackletType[i]));

        TCanvas* c3 = new TCanvas(Form("c3_%i", TrackletType[i]), "", 600, 600);
        hempty->SetStats(0);
        hempty->SetAxisRange(0.8, 1.2, "Y");
        hempty->SetTitle(";#eta;Empty Event correction");
        hempty->Draw();
        draw_legend(hempty, label);
        c3->SaveAs(Form("figs/corrs/empty-event-%s-%i.png", label, TrackletType[i]));

        finput->Close();
    }

    return 0;
}

void draw_legend(TH1F* h1, const char* label) {
    TLegend* l1 = new TLegend(0.32, 0.32, 0.6, 0.4);
    l1->SetFillStyle(0);
    l1->SetBorderSize(0);
    l1->AddEntry(h1, label, "p");
    l1->Draw();
}

int main(int argc, char* argv[]) {
    if (argc == 2)
        return draw_corr_plots(argv[1]);
    else
        return 1;
}
