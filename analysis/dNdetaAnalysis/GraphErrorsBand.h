#include "TError.h"
#include "TGraph.h"

TGraph* GetTGraphErrorBand(TH1F* hist, Double_t xoffset = 0) {
   // returns a band surrounding a TGraphErrors or TGraphAsymmErrors
   // xoffset: points at the edges can be moved such that band does
   //   not pass through the marker. It is assumed that points are
   //   in increasing order in x.
   TGraph tmp(hist);
   TGraph* ing = &tmp;
   if (!ing) return 0;
   Int_t n = ing->GetN();
   if (!n) return 0;
   TGraph* outg = new TGraph(2*n+1);
   Double_t* x = ing->GetX();
   Double_t* y = ing->GetY();
   Double_t* yerr = ing->GetEY();
   Double_t* yerrlow = ing->GetEYlow();
   // Double_t* yerrhi = ing->GetEYhigh();
   if (yerr && yerrlow) {
      Error("GetTGraphErrorBand", "both yerr and yerrlow are here, "
	    " GetTGraphErrorsBand is written with wrong assumptions");
      return 0;
   }
   for (Int_t i=0; i<n; i++) {
      Double_t shiftedx = x[i];
      if (i == 0) shiftedx -= xoffset;
      if (i == n-1) shiftedx += xoffset;
      outg->SetPoint(i, shiftedx, y[i]-hist->GetBinError(i+1));
      outg->SetPoint(2*n-1-i, shiftedx, y[i]+hist->GetBinError(i+1));
   }
   outg->SetPoint(2*n, x[0], y[0]-hist->GetBinError(1));
   outg->SetLineColor(ing->GetLineColor());
   outg->SetLineWidth(ing->GetLineWidth());
   outg->SetLineStyle(ing->GetLineStyle());
   outg->SetFillColor(ing->GetFillColor());
   outg->SetFillStyle(ing->GetFillStyle());
   return outg;
}

TGraph* GetErrorBand(TH1F* hist, Double_t ratio1, Double_t ratio2, Double_t xoffset) {
   // returns a band surrounding a TGraphErrors or TGraphAsymmErrors
   // xoffset: points at the edges can be moved such that band does
   //   not pass through the marker. It is assumed that points are
   //   in increasing order in x.
   TGraph tmp(hist);
   TGraph* ing = &tmp;
   if (!ing) return 0;
   Int_t n = ing->GetN();
   if (!n) return 0;
   TGraph* outg = new TGraph(2*n-11);
   Double_t* x = ing->GetX();
   Double_t* y = ing->GetY();
   Double_t* yerr = ing->GetEY();
   Double_t* yerrlow = ing->GetEYlow();
   // Double_t* yerrhi = ing->GetEYhigh();
   if (yerr && yerrlow) {
      Error("GetTGraphErrorBand", "both yerr and yerrlow are here, "
	    " GetTGraphErrorsBand is written with wrong assumptions");
      return 0;
   }
   for (Int_t i=3; i<n-3; i++) {
      Double_t shiftedx = x[i];
      if (i == 3) shiftedx -= xoffset;
      if (i == n-4) shiftedx += xoffset;
      outg->SetPoint(i-3, shiftedx, y[i]*(1-ratio2));
      outg->SetPoint(2*n-i-10, shiftedx, y[i]*(1+ratio1));
   }
   outg->SetPoint(2*n-12, x[3]-xoffset, y[3]*(1-ratio2));

   outg->SetLineColor(ing->GetLineColor());
   outg->SetLineWidth(ing->GetLineWidth());
   outg->SetLineStyle(ing->GetLineStyle());
   outg->SetFillColor(17);
   outg->SetFillStyle(1002);
   return outg;
}

TGraph* GetErrorBand(TH1F* hist, Double_t* ratio1, Double_t* ratio2, Double_t xoffset) {
   // returns a band surrounding a TGraphErrors or TGraphAsymmErrors
   // xoffset: points at the edges can be moved such that band does
   //   not pass through the marker. It is assumed that points are
   //   in increasing order in x.
   TGraph tmp(hist);
   TGraph* ing = &tmp;
   if (!ing)  return 0;
   Int_t n = ing->GetN();
   if (!n) return 0;
   TGraph* outg = new TGraph(2*n-13);
   Double_t* x = ing->GetX();
   Double_t* y = ing->GetY();
   Double_t* yerr = ing->GetEY();
   Double_t* yerrlow = ing->GetEYlow();
   // Double_t* yerrhi = ing->GetEYhigh();
   if (yerr && yerrlow) {
      Error("GetTGraphErrorBand", "both yerr and yerrlow are here, "
	    " GetTGraphErrorsBand is written with wrong assumptions");
      return 0;
   }
   for (Int_t i=4; i<n-4; i++) {
      Double_t shiftedx = x[i];
      if (i == 4) shiftedx -= xoffset;
      if (i == n-5) shiftedx += xoffset;
      outg->SetPoint(i-3, shiftedx, y[i]*(1-ratio2[i]));
      outg->SetPoint(2*n-i-12, shiftedx, y[i]*(1+ratio1[i]));
   }
   outg->SetPoint(2*n-15, x[4]-xoffset, y[4]*(1-ratio2[4]));

   outg->SetLineColor(ing->GetLineColor());
   outg->SetLineWidth(ing->GetLineWidth());
   outg->SetLineStyle(ing->GetLineStyle());
   outg->SetFillColor(17);
   outg->SetFillStyle(1002);
   outg->Print();
   return outg;
}

TGraph* GetErrorBand2(TH1F* hist, Double_t* ratio1, Double_t* ratio2, Double_t xoffset) {
   // returns a band surrounding a TGraphErrors or TGraphAsymmErrors
   // xoffset: points at the edges can be moved such that band does
   //   not pass through the marker. It is assumed that points are
   //   in increasing order in x.
   TGraph tmp(hist);
   TGraph* ing = &tmp;
   if (!ing)  return 0;
   Int_t n = ing->GetN();
   if (!n) return 0;
   TGraph* outg = new TGraph(2*n+1-8);
   Double_t* x = ing->GetX();
   Double_t* y = ing->GetY();
   Double_t* yerr = ing->GetEY();
   Double_t* yerrlow = ing->GetEYlow();
   // Double_t* yerrhi = ing->GetEYhigh();
   if (yerr && yerrlow) {
      Error("GetTGraphErrorBand", "both yerr and yerrlow are here, "
	    " GetTGraphErrorsBand is written with wrong assumptions");
      return 0;
   }
   for (Int_t i=2; i<n-2; i++) {
      Double_t shiftedx = x[i];
      cout << x[i] << endl;
      if (i == 2) shiftedx -= xoffset;
      if (i == n-3) shiftedx += xoffset;
      outg->SetPoint(i-1-1, shiftedx, y[i]*(1-ratio2[i]));
      outg->SetPoint(2*n-6-i-1, shiftedx, y[i]*(1+ratio1[i]));
      outg->SetPoint(2*n-6-1, x[1]-xoffset, y[i]*(1-ratio2[i]));
   }

   outg->SetLineColor(ing->GetLineColor());
   outg->SetLineWidth(ing->GetLineWidth());
   outg->SetLineStyle(ing->GetLineStyle());
   outg->SetFillColor(17);
   outg->SetFillStyle(1002);
   return outg;
}

TGraph* GetTGraphNoErrors(TGraph* ing) {
   // returns a TGraph with no errors
   if (!ing) return 0;
   Int_t n = ing->GetN();
   if (!n) return 0;
   TGraph* outg = new TGraph(n);
   Double_t* x = ing->GetX();
   Double_t* y = ing->GetY();

   for (Int_t i=0; i<n; i++)
      outg->SetPoint(i, x[i], y[i]);
   outg->SetLineColor(ing->GetLineColor());
   outg->SetLineWidth(ing->GetLineWidth());
   outg->SetLineStyle(ing->GetLineStyle());
   outg->SetFillColor(ing->GetFillColor());
   outg->SetFillStyle(ing->GetFillStyle());
   outg->SetMarkerColor(ing->GetMarkerColor());
   outg->SetMarkerSize(ing->GetMarkerSize());
   outg->SetMarkerStyle(ing->GetMarkerStyle());
   return outg;
}
