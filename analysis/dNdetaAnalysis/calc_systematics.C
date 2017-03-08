#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <fstream>
#include <vector>
#include <string>

#include "systematics.h"
#include "error_bands.h"

#define _NTYPES 4

std::string hist_labels[_NTYPES] = {
    "h12", "h13", "h23", "havg"
};

std::string sys_types[6] = {
    "split", "drop", "smear", "dphi", "mult1", "noise"
};

std::string fit_funcs[6] = {
    "pol4", "pol2", "pol4", "pol2", "pol2", "pol6"
};

int options[6] = {
    2, 2, 0, 2, 2, 2
};

int calc_systematics(const char* nominal_file, const char* list, const char* label) {
    TH1::AddDirectory(kFALSE);
    TH1::SetDefaultSumw2(kTRUE);

    TFile* fnominal = new TFile(nominal_file, "read");
    TH1F* hnominals[_NTYPES] = {0};
    for (int i=0; i<4; ++i)
        hnominals[i] = (TH1F*)fnominal->Get(hist_labels[i].c_str());

    std::vector<std::string> file_list;
    std::ifstream file_stream(list);
    if (!file_stream) return 1;
    std::string line;
    while (std::getline(file_stream, line))
        file_list.push_back(line);

    std::size_t nfiles = file_list.size();
    if (!nfiles) {printf("0 total files!\n"); return 1;}

    TFile* fsys[nfiles] = {0};
    for (std::size_t i=0; i<nfiles; ++i)
        fsys[i] = new TFile(file_list[i].c_str(), "read");

    TFile* fout = new TFile(Form("rootfiles/%s-systematics.root", label), "recreate");

    total_sys_var_t* total_sys_vars[_NTYPES] = {0};
    sys_var_t* sys_vars[_NTYPES][nfiles] = {0};
    for (int i=0; i<_NTYPES; ++i) {
        total_sys_vars[i] = new total_sys_var_t(hist_labels[i], hnominals[i]);

        for (std::size_t j=0; j<nfiles; ++j) {
            sys_vars[i][j] = new sys_var_t(hist_labels[i], sys_types[j], hnominals[i], (TH1F*)fsys[j]->Get(hist_labels[i].c_str()));
            sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), "pol2");
            sys_vars[i][j]->write();

            total_sys_vars[i]->add_sys_var(sys_vars[i][j], options[j]);
        }
        total_sys_vars[i]->write();
    }

    for (std::size_t i=0; i<nfiles; ++i)
        fsys[i]->Close();

    TCanvas* c1[_NTYPES] = {0};
    for (int i=0; i<_NTYPES; ++i) {
        c1[i] = new TCanvas(Form("sys_%s", hist_labels[i].c_str()), "", 900, 600);

        c1[i]->Divide(3, 2);
        for (std::size_t j=0; j<nfiles; ++j) {
            c1[i]->cd(j+1);
            sys_vars[i][j]->get_diff_abs()->Draw();
        }

        c1[i]->SaveAs(Form("figs/systematics/%s-%s.png", c1[i]->GetName(), label));
    }

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 4)
       return calc_systematics(argv[1], argv[2], argv[3]);
    else
        return 1;
}
