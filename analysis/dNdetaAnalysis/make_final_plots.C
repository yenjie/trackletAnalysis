#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"

#include <stdarg.h>

#include <fstream>
#include <vector>
#include <string>

#include "alice.h"
#include "error_bands.h"

#include "predictions_hijing.h"
#include "predictions_cgc.h"

void set_style(TH1F* h1, int style, float size, int colour);
void set_mc_style(TH1F* h1, int colour);
void draw_legend(int nhists, ...);
void draw_predictions_legend(int ngraphs, ...);
void draw_cms_prelim();

int make_final_plots(const char* list, const char* output_file, const char* gen_list) {
    TH1::AddDirectory(kFALSE);

    std::vector<std::string> file_list;
    std::ifstream file_stream(list);
    if (!file_stream) return 1;
    std::string line;
    while (std::getline(file_stream, line))
        file_list.push_back(line);

    std::size_t nfiles = file_list.size();
    if (nfiles != 4) {return 1;}

    std::vector<std::string> gen_file_list;
    std::ifstream gen_file_stream(gen_list);
    if (!gen_file_stream) return 1;
    std::string gen_line;
    while (std::getline(gen_file_stream, gen_line))
        gen_file_list.push_back(gen_line);

    TFile* f5tev_gen = new TFile(gen_file_list[0].c_str(), "read");
    TH1F* hgen_5tev = (TH1F*)f5tev_gen->Get("hTruthWOSelection")->Clone("hgen_5tev");
    f5tev_gen->Close();

    TFile* f8tev_gen = new TFile(gen_file_list[1].c_str(), "read");
    TH1F* hgen_8tev = (TH1F*)f8tev_gen->Get("hTruthWOSelection")->Clone("hgen_8tev");
    f8tev_gen->Close();

    TFile* f8tev_hijing = new TFile(gen_file_list[2].c_str(), "read");
    TH1F* hhj_8tev = (TH1F*)f8tev_hijing->Get("hTruthWOSelection")->Clone("hhj_8tev");
    f8tev_hijing->Close();

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

    set_mc_style(hgen_5tev, 46);
    set_mc_style(hgen_8tev, 46);
    set_mc_style(hhj_8tev, 40);

    set_style(havg_5tev, 21, 0.7, 38);
    set_style(havg_8tev, 21, 0.7, 42);

    TH1F* halice = get_alice_5tev();
    TH1F* halice_sys = get_alice_sys_5tev();

    set_style(halice, 20, 0.6, 30);

    TGraph* gHIJING = get_predictions_hijing();
    gHIJING->SetLineColor(1);
    gHIJING->SetLineStyle(1);
    gHIJING->SetLineWidth(2);
    gHIJING->SetMarkerSize(0);

    TGraph* gHIJING_nosh = get_predictions_hijing_nosh();
    gHIJING_nosh->SetLineColor(1);
    gHIJING_nosh->SetLineStyle(2);
    gHIJING_nosh->SetLineWidth(2);
    gHIJING_nosh->SetMarkerSize(0);

    TGraph* gCGC = get_predictions_cgc();
    gCGC->SetLineColor(38);
    gCGC->SetLineStyle(1);
    gCGC->SetLineWidth(2);
    gCGC->SetMarkerSize(0);

    TGraph* gr = new TGraph();
    gr->SetFillStyle(1001);

    TCanvas* c1 = new TCanvas("c1", "", 500, 500);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    gPad->SetTicky();
    havg_5tev->Draw("p e x0");
    hgen_5tev->Draw("same hist c");
    gr->SetFillColorAlpha(38, 0.7);
    draw_sys_unc(gr, havg_5tev, hsys_5tev);
    havg_5tev->Draw("same p e x0");
    draw_legend(2, havg_5tev, hgen_5tev);
    draw_cms_prelim();
    c1->SaveAs("figs/final/CMS-5TeV.pdf");

    TCanvas* c2a = new TCanvas("c2a", "", 500, 500);
    c2a->SetLeftMargin(0.12);
    c2a->SetBottomMargin(0.12);
    gPad->SetTicky();
    havg_8tev->Draw("p e x0");
    hgen_8tev->Draw("same hist c");
    hhj_8tev->Draw("same hist c");
    gr->SetFillColorAlpha(42, 0.7);
    draw_sys_unc(gr, havg_8tev, hsys_8tev);
    havg_8tev->Draw("same p e x0");
    draw_legend(3, havg_8tev, hgen_8tev, hhj_8tev);
    draw_cms_prelim();
    c2a->SaveAs("figs/final/CMS-8TeV.pdf");

    havg_8tev->SetAxisRange(0, 35, "Y");

    TCanvas* c2 = new TCanvas("c2", "", 500, 500);
    c2->SetLeftMargin(0.12);
    c2->SetBottomMargin(0.12);
    gPad->SetTicky();
    havg_8tev->Draw("p e x0");
    hgen_8tev->Draw("same hist c");
    hhj_8tev->Draw("same hist c");
    gHIJING->Draw("same l z");
    gHIJING_nosh->Draw("same l z");
    gCGC->Draw("same l z");
    gr->SetFillColorAlpha(42, 0.7);
    draw_sys_unc(gr, havg_8tev, hsys_8tev);
    havg_8tev->Draw("same p e x0");
    draw_legend(3, havg_8tev, hgen_8tev, hhj_8tev);
    draw_predictions_legend(3, gHIJING, gHIJING_nosh, gCGC);
    draw_cms_prelim();
    c2->SaveAs("figs/final/CMS-8TeV-predictions.pdf");

    havg_8tev->SetAxisRange(0, 30, "Y");

    TCanvas* c3 = new TCanvas("c3", "", 500, 500);
    c3->SetLeftMargin(0.12);
    c3->SetBottomMargin(0.12);
    gPad->SetTicky();
    havg_5tev->Draw("p e x0");
    gr->SetFillColorAlpha(38, 0.7);
    draw_sys_unc(gr, havg_5tev, hsys_5tev);
    havg_5tev->Draw("same p e x0");
    gr->SetFillColorAlpha(42, 0.7);
    draw_sys_unc(gr, havg_8tev, hsys_8tev);
    havg_8tev->Draw("same p e x0");
    draw_legend(2, havg_8tev, havg_5tev);
    draw_cms_prelim();
    c3->SaveAs("figs/final/CMS-58TeV.pdf");

    TCanvas* c4 = new TCanvas("c4", "", 500, 500);
    c4->SetLeftMargin(0.12);
    c4->SetBottomMargin(0.12);
    gPad->SetTicky();
    havg_5tev->Draw("p e x0");
    hgen_5tev->Draw("same hist c");
    gr->SetFillColorAlpha(30, 0.7);
    draw_sys_unc(gr, halice, halice_sys);
    halice->Draw("same p e x0");
    gr->SetFillColorAlpha(38, 0.7);
    draw_sys_unc(gr, havg_5tev, hsys_5tev);
    havg_5tev->Draw("same p e x0");
    draw_legend(3, havg_5tev, halice, hgen_5tev);
    draw_cms_prelim();
    c4->SaveAs("figs/final/CMS-5TeV-ALICE.pdf");

    TCanvas* c5 = new TCanvas("c5", "", 500, 500);
    c5->SetLeftMargin(0.12);
    c5->SetBottomMargin(0.12);
    gPad->SetTicky();
    havg_5tev->Draw("p e x0");
    gr->SetFillColorAlpha(30, 0.7);
    draw_sys_unc(gr, halice, halice_sys);
    halice->Draw("same p e x0");
    gr->SetFillColorAlpha(38, 0.7);
    draw_sys_unc(gr, havg_5tev, hsys_5tev);
    havg_5tev->Draw("same p e x0");
    gr->SetFillColorAlpha(42, 0.7);
    draw_sys_unc(gr, havg_8tev, hsys_8tev);
    havg_8tev->Draw("same p e x0");
    draw_legend(3, havg_8tev, havg_5tev, halice);
    draw_cms_prelim();
    c5->SaveAs("figs/final/CMS-58TeV-ALICE.pdf");

    TFile* foutput = new TFile(output_file, "recreate");
    foutput->Write("", TObject::kOverwrite);

    return 0;
}

void set_style(TH1F* h1, int style, float size, int colour) {
    TAxis* x_axis = h1->GetXaxis();
    TAxis* y_axis = h1->GetYaxis();

    x_axis->SetLabelFont(43);
    x_axis->SetLabelSize(20);
    x_axis->SetLabelOffset(0.012);
    x_axis->SetTitleFont(43);
    x_axis->SetTitleSize(22);
    x_axis->SetTitleOffset(1.2);

    y_axis->SetLabelFont(43);
    y_axis->SetLabelSize(20);
    y_axis->SetLabelOffset(0.012);
    y_axis->SetTitleFont(43);
    y_axis->SetTitleSize(22);
    y_axis->SetTitleOffset(1.2);

    h1->SetLineColor(colour);
    h1->SetMarkerStyle(style);
    h1->SetMarkerSize(size);
    h1->SetMarkerColor(1);

    h1->SetFillStyle(1001);
    h1->SetFillColor(colour);
}

void set_mc_style(TH1F* h1, int colour) {
    h1->SetLineColor(colour);
    h1->SetLineWidth(2);
    h1->SetMarkerSize(0);
}

std::string get_label(std::string name) {
    if (name == "halice")
        return "ALICE 5.02 TeV pPb 2013";
    else if (name == "havg_5tev")
        return "CMS 5.02 TeV pPb";
    else if (name == "havg_8tev")
        return "CMS 8.16 TeV pPb";
    else if (name == "hgen_5tev")
        return "EPOS LHC 5.02 TeV";
    else if (name == "hgen_8tev")
        return "EPOS LHC 8.16 TeV";
    else if (name == "hhj_8tev")
        return "HIJING 1.3 8.16 TeV";
    else if (name == "gHIJING")
        return "HIJING 2.1 8.16 TeV";
    else if (name == "gHIJING_nosh")
        return "HIJING 2.1 (no shadowing)";
    else if (name == "gCGC")
        return "KLN model";
    else
        return {};
}

void draw_legend(int nhists, ...) {
    va_list hist_list;
    va_start(hist_list, nhists);

    TLegend* l1 = new TLegend(0.3, 0.45-nhists*0.05, 0.6, 0.45);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->SetTextFont(43);
    l1->SetTextSize(18);
    for (int i=0; i<nhists; ++i) {
        TH1F* h1 = va_arg(hist_list, TH1F*);
        if (strstr(h1->GetName(), "ha") != NULL)
            l1->AddEntry(h1, get_label(h1->GetName()).c_str(), "pf");
        else
            l1->AddEntry(h1, get_label(h1->GetName()).c_str(), "l");
    }
    l1->Draw();

    va_end(hist_list);
}

void draw_predictions_legend(int ngraphs, ...) {
    va_list graph_list;
    va_start(graph_list, ngraphs);

    TLegend* l1 = new TLegend(0.3, 0.3-ngraphs*0.05, 0.6, 0.3);
    l1->SetBorderSize(0);
    l1->SetTextFont(43);
    l1->SetTextSize(18);
    for (int i=0; i<ngraphs; ++i) {
        TGraph* g1 = va_arg(graph_list, TGraph*);
        l1->AddEntry(g1, get_label(g1->GetName()).c_str(), "l");
    }
    l1->Draw();

    va_end(graph_list);
}

void draw_cms_prelim() {
    TLatex* latexCMS = new TLatex();
    latexCMS->SetTextFont(63);
    latexCMS->SetTextSize(25);
    latexCMS->DrawLatexNDC(0.16, 0.84, "CMS");

    TLatex* latexPrelim = new TLatex();
    latexPrelim->SetTextFont(53);
    latexPrelim->SetTextSize(18);
    latexPrelim->DrawLatexNDC(0.16, 0.79, "Preliminary");
}

int main(int argc, char* argv[]) {
    if (argc == 4)
        return make_final_plots(argv[1], argv[2], argv[3]);
    else
        return 1;
}
