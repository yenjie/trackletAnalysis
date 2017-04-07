#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TCanvas.h"

int unboost(const char* input, const char* output) {
    TFile* finput = new TFile(input, "read");
    TTree* tinput = (TTree*)finput->Get("TrackletTree12");

    int npart = 0;
    tinput->SetBranchStatus("npart", 1);
    tinput->SetBranchAddress("npart", &npart);

#define _SET_BRANCH(tree, branch)               \
    tree->SetBranchStatus(#branch, 1);          \
    tree->SetBranchAddress(#branch, branch);

    float pt[720];
    float eta[720];
    float phi[720];
    _SET_BRANCH(tinput, pt);
    _SET_BRANCH(tinput, eta);
    _SET_BRANCH(tinput, phi);

    TFile* foutput = new TFile(output, "recreate");

    TH1D* hlab = new TH1D("hlab", "", 30, -3, 3);
    TH1D* hcm = new TH1D("hcm", "", 30, -3, 3);

    uint64_t nentries = tinput->GetEntries();
    for (uint64_t i=0; i<nentries; ++i) {
        tinput->GetEntry(i);

        for (int j=0; j<npart; ++j) {
            TLorentzVector part;
            part.SetPtEtaPhiM(pt[j], eta[j], phi[j], 0);
            part.Boost(0, 0, 0.434);

            hlab->Fill(eta[j]);
            hcm->Fill(part.Eta());
        }
    }

    hlab->SetMarkerSize(0);
    hcm->SetLineColor(2);
    hcm->SetMarkerSize(0);

    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    hlab->Draw("c");
    hcm->Draw("same c");
    c1->SaveAs("unboost.png");

    foutput->Write("", TObject::kOverwrite);
    foutput->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return unboost(argv[1], argv[2]);
    else
        return 1;
}
