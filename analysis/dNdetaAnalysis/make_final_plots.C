#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"

#include <fstream>
#include <vector>
#include <string>

#include "alice.h"
#include "error_bands.h"

void set_style(TH1F* h1, int style, float size, int colour);
void draw_cms_prelim();

int make_final_plots(const char* list, const char* output_file) {
    TH1::AddDirectory(kFALSE);

    std::vector<std::string> file_list;
    std::ifstream file_stream(list);
    if (!file_stream) return 1;
    std::string line;
    while (std::getline(file_stream, line))
        file_list.push_back(line);

    std::size_t nfiles = file_list.size();
    if (nfiles != 4) {return 1;}

    TFile* f5tev = new TFile(file_list[0].c_str(), "read");
    TH1F* havg_5tev = (TH1F*)f5tev->Get("havg")->Clone("havg_5tev");
    f5tev->Close();

    TFile* f8tev = new TFile(file_list[1].c_str(), "read");
    TH1F* havg_8tev = (TH1F*)f8tev->Get("havg")->Clone("havg_8tev");
    f8tev->Close();

    TFile* f5tev_sys = new TFile(file_list[2].c_str(), "read");
    TH1F* hsys_5tev = (TH1F*)f5tev_sys->Get("havg_systematics")->Clone("hsys_5tev");
    f5tev_sys->Close();

    TFile* f8tev_sys = new TFile(file_list[3].c_str(), "read");
    TH1F* hsys_8tev = (TH1F*)f8tev_sys->Get("havg_systematics")->Clone("hsys_8tev");
    f8tev_sys->Close();

    set_style(havg_5tev, 21, 0.72, 38);
    set_style(havg_8tev, 21, 0.72, 42);

    TH1F* halice = get_alice_5tev();
    TH1F* halice_sys = get_alice_sys_5tev();

    set_style(halice, 20, 0.6, 30);

    TGraph* gr = new TGraph();
    gr->SetFillStyle(1001);

    TCanvas* c1 = new TCanvas("c1", "", 500, 500);
    gPad->SetTicky();
    havg_5tev->Draw();
    gr->SetFillColorAlpha(38, 0.7);
    draw_sys_unc(gr, havg_5tev, hsys_5tev);
    havg_5tev->Draw("p same");
    draw_cms_prelim();
    c1->SaveAs("5tev.png");

    TCanvas* c2 = new TCanvas("c2", "", 500, 500);
    gPad->SetTicky();
    havg_8tev->Draw();
    gr->SetFillColorAlpha(42, 0.7);
    draw_sys_unc(gr, havg_8tev, hsys_8tev);
    havg_8tev->Draw("p same");
    draw_cms_prelim();
    c2->SaveAs("8tev.png");

    TCanvas* c3 = new TCanvas("c3", "", 500, 500);
    gPad->SetTicky();
    havg_5tev->Draw();
    gr->SetFillColorAlpha(38, 0.7);
    draw_sys_unc(gr, havg_5tev, hsys_5tev);
    havg_5tev->Draw("p same");
    gr->SetFillColorAlpha(42, 0.7);
    draw_sys_unc(gr, havg_8tev, hsys_8tev);
    havg_8tev->Draw("p same");
    draw_cms_prelim();
    c3->SaveAs("58tev.png");

    TCanvas* c4 = new TCanvas("c4", "", 500, 500);
    gPad->SetTicky();
    havg_5tev->Draw();
    gr->SetFillColorAlpha(30, 0.7);
    draw_sys_unc(gr, halice, halice_sys);
    halice->Draw("p same");
    gr->SetFillColorAlpha(38, 0.7);
    draw_sys_unc(gr, havg_5tev, hsys_5tev);
    havg_5tev->Draw("p same");
    draw_cms_prelim();
    c4->SaveAs("5tevwalice.png");

    TCanvas* c5 = new TCanvas("c5", "", 500, 500);
    gPad->SetTicky();
    havg_5tev->Draw();
    gr->SetFillColorAlpha(30, 0.7);
    draw_sys_unc(gr, halice, halice_sys);
    halice->Draw("p same");
    gr->SetFillColorAlpha(38, 0.7);
    draw_sys_unc(gr, havg_5tev, hsys_5tev);
    havg_5tev->Draw("p same");
    gr->SetFillColorAlpha(42, 0.7);
    draw_sys_unc(gr, havg_8tev, hsys_8tev);
    havg_8tev->Draw("p same");
    draw_cms_prelim();
    c5->SaveAs("58tevwithalice.png");

    TFile* foutput = new TFile(output_file, "recreate");
    foutput->Write("", TObject::kOverwrite);

    return 0;
}

void set_style(TH1F* h1, int style, float size, int colour) {
    TAxis* x_axis = h1->GetXaxis();
    TAxis* y_axis = h1->GetYaxis();

    x_axis->SetLabelFont(43);
    x_axis->SetLabelSize(16);
    x_axis->SetLabelOffset(0.012);
    x_axis->SetTitleFont(43);
    x_axis->SetTitleSize(18);
    x_axis->SetTitleOffset(1.25);

    y_axis->SetLabelFont(43);
    y_axis->SetLabelSize(16);
    y_axis->SetLabelOffset(0.012);
    y_axis->SetTitleFont(43);
    y_axis->SetTitleSize(18);
    y_axis->SetTitleOffset(1.25);

    h1->SetLineColor(colour);
    h1->SetMarkerStyle(style);
    h1->SetMarkerSize(size);
    h1->SetMarkerColor(1);
}

void draw_cms_prelim() {
    TLatex* latexCMS = new TLatex();
    latexCMS->SetTextFont(63);
    latexCMS->SetTextSize(20);
    latexCMS->DrawLatexNDC(0.15, 0.84, "CMS");

    TLatex* latexPrelim = new TLatex();
    latexPrelim->SetTextFont(53);
    latexPrelim->SetTextSize(15);
    latexPrelim->DrawLatexNDC(0.15, 0.8, "Preliminary");
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return make_final_plots(argv[1], argv[2]);
    else
        return 1;
}
