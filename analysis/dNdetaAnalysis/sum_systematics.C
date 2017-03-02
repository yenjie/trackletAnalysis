#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include "error_bands.h"

int systematics(const char* fname) {
   return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 2)
       return systematics(argv[1]);
    else
        return 1;
}
