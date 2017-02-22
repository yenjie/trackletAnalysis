#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include <fstream>

int calc_mult_weights(const char* fnMC, const char* fnData) {
    TFile* fMC = new TFile(fnMC, "read");
    TFile* fData = new TFile(fnData, "read");

    int TrackletType[3] = {12, 13, 23};
    for (int i=0; i<3; ++i) {
        TTree* tMC = (TTree*)fMC->Get(Form("TrackletTree%i", TrackletType[i]));
        TTree* tData = (TTree*)fData->Get(Form("TrackletTree%i", TrackletType[i]));

        TH1F* hMC = new TH1F("hMC", "", 61, -5, 300);
        TH1F* hData = new TH1F("hData", "", 61, -5, 300);

        tMC->Draw("mult>>hMC", "vz[1]>-99 && passHLT");
        tData->Draw("mult>>hData", "vz[1]>-99 && passHLT");

        hData->Divide(hMC);

        std::ofstream ofs(Form("mult_weights_%i.h", TrackletType[i]), std::ofstream::out);
        hData->SavePrimitive(ofs);

        hData->Delete();
        hMC->Delete();
    }

    return 0;
}
