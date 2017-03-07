#include "TFile.h"
#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"

int get_scale_factor(TH2F* h1);

int draw_beta_alpha(const char* label) {
    int TrackletType[3] = {12, 13, 23};

    for (int i=0; i<3; ++i) {
        TFile* finput = new TFile(Form("correction/correction-%i-%s.root", TrackletType[i], label), "read");

        TH2F* haccep = (TH2F*)finput->Get("hAcceptance1D");

        TH3F* hbeta = (TH3F*)((TH3F*)finput->Get("hbeta"))->Clone();
        TH3F* halpha = (TH3F*)((TH3F*)finput->Get("halpha"))->Clone();

        TH2F* hbeta_xz = (TH2F*)hbeta->Project3D("zx");
        TH2F* halpha_xz = (TH2F*)halpha->Project3D("zx");

        hbeta_xz->Scale(1./20);
        halpha_xz->Scale(1./20);

        TH1F* hbeta_y = (TH1F*)hbeta->Project3D("y");
        TH1F* halpha_y = (TH1F*)halpha->Project3D("y");

        hbeta_y->Scale(1./get_scale_factor(haccep));
        halpha_y->Scale(1./get_scale_factor(haccep));

        TCanvas* c1 = new TCanvas(Form("c1_%i", TrackletType[i]), "", 600, 600);
        hbeta_xz->SetTitle("beta");
        hbeta_xz->SetStats(0);
        hbeta_xz->Draw("colz");
        c1->SaveAs(Form("figs/beta-alpha/beta-%i-%s.png", TrackletType[i], label));

        TCanvas* c2 = new TCanvas(Form("c2_%i", TrackletType[i]), "", 600, 600);
        halpha_xz->SetTitle("alpha");
        halpha_xz->SetStats(0);
        halpha_xz->Draw("colz");
        c2->SaveAs(Form("figs/beta-alpha/alpha-%i-%s.png", TrackletType[i], label));

        TCanvas* c3 = new TCanvas(Form("c3_%i", TrackletType[i]), "", 600, 600);
        hbeta_y->SetTitle("beta");
        hbeta_y->SetAxisRange(0, 0.4, "Y");
        hbeta_y->SetStats(0);
        hbeta_y->Draw();
        c3->SaveAs(Form("figs/beta-alpha/beta-mult-%i-%s.png", TrackletType[i], label));

        TCanvas* c4 = new TCanvas(Form("c4_%i", TrackletType[i]), "", 600, 600);
        halpha_y->SetTitle("alpha");
        halpha_y->SetAxisRange(0.5, 2, "Y");
        halpha_y->SetStats(0);
        halpha_y->Draw();
        c4->SaveAs(Form("figs/beta-alpha/alpha-mult-%i-%s.png", TrackletType[i], label));

        finput->Close();
    }

    return 0;
}

int get_scale_factor(TH2F* h1) {
    int nbins = 0;
    for (int i=0; i<h1->GetNcells(); ++i)
        if (h1->GetBinContent(i) != 0)
            ++nbins;

    return nbins;
}

int main(int argc, char* argv[]) {
    if (argc == 2)
        return draw_beta_alpha(argv[1]);
    else
        return 1;
}
