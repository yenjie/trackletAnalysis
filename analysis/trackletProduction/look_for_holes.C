#define _XLEDGE(x) (fabs(x_low_edge - x) < 0.01)
#define _YLEDGE(x) (fabs(y_low_edge - x) < 0.01)

#define _YLEDGE_RANGE(low, high) (y_low_edge > (low - 0.01) && y_low_edge < (high + 0.01))
#define _XLEDGE_RANGE(low, high) (x_low_edge > (low - 0.01) && x_low_edge < (high + 0.01))

int compass[4] = {1, -202, -1, 202};

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"

#include <vector>
#include <fstream>

int trace_hits(TH2D* hholes, int bin, int direction, std::vector<int>& trace, TH2D* htrace);
int get_n_neighbours(TH2D* hhist, int bin);

int look_for_holes(const char* mc_fname, const char* data_fname, int start, int end) {
    TFile* fmc = new TFile(mc_fname, "read");
    TTree* tmc = (TTree*)fmc->Get("ana/PixelTree");
    TH2D* hmc = new TH2D("hmc", "", 200, -3.5, 3.5, 200, -3, 3);
    printf("projecting mc hitmap\n");
    tmc->Draw("eta1:phi1>>hmc", "", "colz");

    TFile* fdata = new TFile(data_fname, "read");
    TTree* tdata = (TTree*)fdata->Get("ana/PixelTree");
    TH2D* hdata = new TH2D("hdata", "", 200, -3.5, 3.5, 200, -3, 3);
    printf("projecting data hitmap\n");
    tdata->Draw("-(eta1/abs(eta1)*r1/tan(2*atan(exp(-abs(eta1))))+0.551191)/abs((eta1/abs(eta1)*r1/tan(2*atan(exp(-abs(eta1))))+0.551191))*log(tan(atan(r1/abs(eta1/abs(eta1)*r1/tan(2*atan(exp(-abs(eta1))))+0.551191))/2)):phi1>>hdata", Form("Entry$>=%i && Entry$<%i", start, end), "colz");

    printf("calculating ratio\n");
    TH2D* hratio = (TH2D*)hdata->Clone("hratio");
    hratio->SetStats(0);
    hratio->Divide(hmc);

    TH2D* hholes = (TH2D*)hratio->Clone("hholes");
    hholes->SetStats(0);
    hholes->Scale(1/3.5);

    TH1D* hratio_dist = new TH1D("hratio_dist", "", 200, 0, 4);

    for (int i=0; i<hratio->GetNcells(); ++i) {
        if (hratio->GetBinContent(i) != 0)
            hratio_dist->Fill(hratio->GetBinContent(i));
        if (hratio->GetBinContent(i) > 4)
            hratio->SetBinContent(i, 4.01);
        if (hholes->GetBinContent(i) > 0.7)
            hholes->SetBinContent(i, 0);

        float y_low_edge = hholes->GetYaxis()->GetBinLowEdge(i / 202);
        if (hholes->GetBinContent(i) > 0.6 && _YLEDGE_RANGE(-2.16, -1.17))
            hholes->SetBinContent(i, 0);

        float x_low_edge = hholes->GetXaxis()->GetBinLowEdge(i % 202);
        if (_XLEDGE(1.400) || _XLEDGE(1.540) || _XLEDGE(1.575) || _XLEDGE(2.065) ||
            _XLEDGE(2.765) || _XLEDGE(2.415) || _XLEDGE(0.700) || _XLEDGE(0.000) ||
            _XLEDGE(-1.61) || _XLEDGE(-0.70) || _XLEDGE(-3.15) || _XLEDGE(-2.835) ||
            _XLEDGE(-1.75) || _XLEDGE(-2.45) || _XLEDGE(-1.40) || _XLEDGE(-1.575) ||
            _XLEDGE(-2.10) || _XLEDGE(1.715) || _XLEDGE(0.315) || _XLEDGE(-1.05) ||
            _XLEDGE(2.450)) {
            if (hholes->GetBinContent(i) > 0.45)
                hholes->SetBinContent(i, 0);
        }
    }

    for (int i=0; i<hholes->GetNcells(); ++i) {
        if (get_n_neighbours(hholes, i) == 0)
            hholes->SetBinContent(i, 0);
    }

    printf("tracing lines\n");
    TH2D* htrace = new TH2D("htrace", "", 200, -3.5, 3.5, 200, -3, 3);
    for (int i=0; i<hholes->GetNcells(); ++i) {
        if (hholes->GetBinContent(i) == 0 || htrace->GetBinContent(i) != 0)
            continue;

        float y_low_edge = hholes->GetYaxis()->GetBinLowEdge(i / 202);
        float x_low_edge = hholes->GetXaxis()->GetBinLowEdge(i % 202);
        if (_YLEDGE_RANGE(2.52, 2.55) || _YLEDGE(2.13) ||
            _YLEDGE_RANGE(2.25, 2.28) ||
            (_YLEDGE(2.430) && !_XLEDGE_RANGE(1.0, 1.6)) ||
            _YLEDGE_RANGE(1.74, 1.89) || _YLEDGE_RANGE(1.14, 1.26) ||
            _YLEDGE_RANGE(-0.06, 0.03) || _YLEDGE_RANGE(-2.58, -2.46) ||
            _YLEDGE_RANGE(-1.89, -1.77) || _YLEDGE_RANGE(-2.31, -2.16) ||
            _YLEDGE_RANGE(-1.29, -1.17)) {
            std::vector<int> trace;
            int trace_val = trace_hits(hholes, i, 0, trace, htrace);
            if (trace_val == get_n_neighbours(hholes, i)) {
                trace.push_back(i);
                for (std::size_t j=0; j<trace.size(); ++j)
                    hholes->SetBinContent(trace[j], 0);
            }
            trace.clear();
        }
    }

#define GET_MACRO(_1,_2,NAME,...) NAME
#define _ZERO(...) GET_MACRO(__VA_ARGS__, _ZERO2, _ZERO1)(__VA_ARGS__)
#define _ZERO1(x) hholes->SetBinContent(x, 0)
#define _ZERO2(low, high) {             \
    for (int i=low; i<=high; ++i)       \
        hholes->SetBinContent(i, 0);    \
}

    _ZERO(33000, 33005); _ZERO(33205, 33207);
    _ZERO(31077); _ZERO(31279); _ZERO(31481);
    _ZERO(3919, 3927); _ZERO(3717); _ZERO(3515); _ZERO(3310, 3313);
    _ZERO(5095, 5100); _ZERO(4890, 4893); _ZERO(4489); _ZERO(4287); _ZERO(4085);
    _ZERO(12101, 12106);
    _ZERO(11846, 11855);
    _ZERO(12443, 12451);
    _ZERO(12606, 12613); _ZERO(12815);
    _ZERO(11989, 11996); _ZERO(12191, 12192);
    _ZERO(12777, 12780); _ZERO(12577, 12579);
    _ZERO(12160, 12170); _ZERO(11963, 11966);
    _ZERO(12558, 12561); _ZERO(12758, 12766);
    _ZERO(12141, 12151);
    _ZERO(12537, 12542); _ZERO(12737, 12740); _ZERO(12742, 12746);
    _ZERO(36776, 36784); _ZERO(36986); _ZERO(36579, 36581);
    _ZERO(28775, 28783); _ZERO(28997, 29006);
    _ZERO(28401, 28409); _ZERO(28421, 28425);
    _ZERO(20766, 20773);
    _ZERO(7605, 7613);
    _ZERO(7620, 7625);

    for (int i=0; i<hholes->GetNcells(); ++i) {
        if (get_n_neighbours(hholes, i) == 0)
            hholes->SetBinContent(i, 0);
    }

    TFile* foutput = new TFile("holes.root", "recreate");
    hdata->Write("", TObject::kOverwrite);
    hmc->Write("", TObject::kOverwrite);
    hratio->Write("", TObject::kOverwrite);
    hholes->Write("", TObject::kOverwrite);
    hratio_dist->Write("", TObject::kOverwrite);

    std::ofstream ofs("hholes.h", std::ofstream::out);
    hholes->SavePrimitive(ofs);

    foutput->Close();

    return 0;
}

int trace_hits(TH2D* hholes, int bin, int direction, std::vector<int>& trace, TH2D* htrace) {
    bin += direction;
    if (hholes->GetBinContent(bin) == 0)
        return 0;
    if (htrace->GetBinContent(bin) != 0)
        return 2;

    htrace->SetBinContent(bin, 1);
    int trace_val = 0;
    for (int i=0; i<4; ++i) {
        if (direction != -compass[i])
            trace_val += trace_hits(hholes, bin, compass[i], trace, htrace);
    }
    if (trace_val < 2) {
        trace.push_back(bin);
        return 1;
    } else {
        return trace_val;
    }
}

int get_n_neighbours(TH2D* hhist, int bin) {
    int n = 0;

    int bin_x = bin % 202;
    int bin_y = bin / 202;
    for (int i=0; i<4; ++i)
        if (bin_x > 0 && bin_x < 201 && bin_y > 0 && bin_y < 201)
            if (hhist->GetBinContent(bin + compass[i]) != 0)
                ++n;

    return n;
}

int main(int argc, char* argv[]) {
    if (argc == 5)
        return look_for_holes(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
    else
        return 1;
}
