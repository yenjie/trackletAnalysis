#include <TCut.h>
#include <TString.h>
#include <iostream>

class selectionCut {
   public:
      selectionCut(bool isMC, int nLumiL = 69, int nLumiH = 145);
      ~selectionCut(){}

      TCut Cut;
      TString evtSelection;
      TCut CutWOVtxCut;
      TString vtxCut;
      TString myCut;
      int VzRangeL;
      int VzRangeH;
};

selectionCut::selectionCut(bool isMC, int nLumiL, int nLumiH) {
   VzRangeL = -20;
   VzRangeH = 20;
   vtxCut = Form("vz[1]<%d&&vz[1]>%d", VzRangeH, VzRangeL);

   evtSelection = "1";
   // evtSelection = "(nHFp>0 || nHFn>0)"; // HF OR
   // evtSelection = "(nHFp>0 && nHFn>0)"; // HF AND
   // evtSelection = "(nHFp>0 ^ nHFn>0)"; // HF XOR

   // evtSelection = "l1TBit[40]>0&&l1TBit[36]!=1&&l1TBit[37]!=1&&l1TBit[38]!=1&&l1TBit[39]!=1";
   // evtSelection = "(nHFp>=1&&nHFn>=1)&&l1TBit[34]!=0&&l1TBit[36]!=1&&l1TBit[37]!=1&&l1TBit[38]!=1&&l1TBit[39]!=1";
   // evtSelection = "(nHFp>=1&&nHFn>=1)&&l1ABit[124]!=0&&l1TBit[36]!=1&&l1TBit[37]!=1&&l1TBit[38]!=1&&l1TBit[39]!=1"; // old default
   // evtSelection = "(nHFp>=1&&nHFn>=1)&&l1TBit[40]!=0&&l1TBit[36]!=1&&l1TBit[37]!=1&&l1TBit[38]!=1&&l1TBit[39]!=1";
   // evtSelection = "(nHFp>=1&&nHFn==0)&&l1TBit[34]!=0&&l1TBit[36]!=1&&l1TBit[37]!=1&&l1TBit[38]!=1&&l1TBit[39]!=1";
   // evtSelection = "(nHFn>=1&&nHFp==0)&&l1TBit[34]!=0&&l1TBit[36]!=1&&l1TBit[37]!=1&&l1TBit[38]!=1&&l1TBit[39]!=1";

   // if (!isMC) evtSelection += Form("&&nLumi>=%d&&nLumi<=%d&&l1ABit[0]==1&&l1ABit[82]==1", nLumiL, nLumiH);
   // if (!isMC) {
   //    evtSelection += Form("&&mult2>=%d&&mult2<=%d&&l1ABit[0]==1&&l1ABit[82]==1", nLumiL, nLumiH);
   //    cout << "Selection on mult2 mode!" << endl;
   // }

   CutWOVtxCut = TCut(evtSelection);
   myCut = vtxCut + "&&" + evtSelection;
   cout << myCut << endl;
   Cut = TCut(myCut);
}
