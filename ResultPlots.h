#include <TH1F.h>

TH1F* tracklet7TeV()
{
   TH1F *hAvg = new TH1F("hAvg","",12, -3, 3);
   hAvg->SetBinContent(2,6.30313);
   hAvg->SetBinContent(3,6.39303);
   hAvg->SetBinContent(4,6.22794);
   hAvg->SetBinContent(5,6.0133);
   hAvg->SetBinContent(6,5.75931);
   hAvg->SetBinContent(7,5.75931);
   hAvg->SetBinContent(8,6.0133);
   hAvg->SetBinContent(9,6.22794);
   hAvg->SetBinContent(10,6.39303);
   hAvg->SetBinContent(11,6.30313);
   hAvg->SetBinError(2,0.182791);
   hAvg->SetBinError(3,0.185398);
   hAvg->SetBinError(4,0.18061);
   hAvg->SetBinError(5,0.174386);
   hAvg->SetBinError(6,0.16702);
   hAvg->SetBinError(7,0.16702);
   hAvg->SetBinError(8,0.174386);
   hAvg->SetBinError(9,0.18061);
   hAvg->SetBinError(10,0.185398);
   hAvg->SetBinError(11,0.182791);
   hAvg->SetMinimum(0);
   hAvg->SetMaximum(7);
   hAvg->SetEntries(25);
   hAvg->SetFillColor(1);
   hAvg->SetFillStyle(0);
   hAvg->SetLineStyle(0);
   hAvg->SetMarkerStyle(4);
   hAvg->SetMarkerSize(1.7);
   hAvg->GetXaxis()->SetTitle("#eta");
   hAvg->GetXaxis()->CenterTitle(true);
   hAvg->GetXaxis()->SetLabelFont(42);
   hAvg->GetXaxis()->SetLabelOffset(0.01);
   hAvg->GetXaxis()->SetLabelSize(0.045);
   hAvg->GetXaxis()->SetTitleSize(0.055);
   hAvg->GetXaxis()->SetTitleFont(42);
   hAvg->GetYaxis()->SetTitle("dN/d#eta");
   hAvg->GetYaxis()->CenterTitle(true);
   hAvg->GetYaxis()->SetLabelFont(42);
   hAvg->GetYaxis()->SetLabelOffset(0.01);
   hAvg->GetYaxis()->SetLabelSize(0.045);
   hAvg->GetYaxis()->SetTitleSize(0.055);
   hAvg->GetYaxis()->SetTitleOffset(1.25);
   hAvg->GetYaxis()->SetTitleFont(42);
   hAvg->GetZaxis()->SetLabelFont(42);
   hAvg->GetZaxis()->SetLabelSize(0.045);
   hAvg->GetZaxis()->SetTitleFont(42);
   return hAvg;
}

TH1F* tracklet7TeVOld()
{
   TH1F *hAvg = new TH1F("hAvg","",12, -3, 3);
   hAvg->SetBinContent(2,6.302286);
   hAvg->SetBinContent(3,6.343485);
   hAvg->SetBinContent(4,6.233619);
   hAvg->SetBinContent(5,6.009152);
   hAvg->SetBinContent(6,5.785029);
   hAvg->SetBinContent(7,5.785029);
   hAvg->SetBinContent(8,6.009152);
   hAvg->SetBinContent(9,6.233619);
   hAvg->SetBinContent(10,6.343485);
   hAvg->SetBinContent(11,6.302286);
   hAvg->SetBinError(2,0.2394869);
   hAvg->SetBinError(3,0.2410524);
   hAvg->SetBinError(4,0.2368775);
   hAvg->SetBinError(5,0.2283478);
   hAvg->SetBinError(6,0.2198311);
   hAvg->SetBinError(7,0.2198311);
   hAvg->SetBinError(8,0.2283478);
   hAvg->SetBinError(9,0.2368775);
   hAvg->SetBinError(10,0.2410524);
   hAvg->SetBinError(11,0.2394869);
   hAvg->SetMinimum(0);
   hAvg->SetMaximum(7);
   hAvg->SetEntries(25);
   hAvg->SetFillColor(1);
   hAvg->SetFillStyle(0);
   hAvg->SetLineStyle(0);
   hAvg->SetMarkerStyle(4);
   hAvg->SetMarkerSize(1.7);
   hAvg->GetXaxis()->SetTitle("#eta");
   hAvg->GetXaxis()->CenterTitle(true);
   hAvg->GetXaxis()->SetLabelFont(42);
   hAvg->GetXaxis()->SetLabelOffset(0.01);
   hAvg->GetXaxis()->SetLabelSize(0.045);
   hAvg->GetXaxis()->SetTitleSize(0.055);
   hAvg->GetXaxis()->SetTitleFont(42);
   hAvg->GetYaxis()->SetTitle("dN/d#eta");
   hAvg->GetYaxis()->CenterTitle(true);
   hAvg->GetYaxis()->SetLabelFont(42);
   hAvg->GetYaxis()->SetLabelOffset(0.01);
   hAvg->GetYaxis()->SetLabelSize(0.045);
   hAvg->GetYaxis()->SetTitleSize(0.055);
   hAvg->GetYaxis()->SetTitleOffset(1.25);
   hAvg->GetYaxis()->SetTitleFont(42);
   hAvg->GetZaxis()->SetLabelFont(42);
   hAvg->GetZaxis()->SetLabelSize(0.045);
   hAvg->GetZaxis()->SetTitleFont(42);
   return hAvg;
}


TH1F* pixelCounting7TeV()
{
   TH1F *corr_result_1 = new TH1F("corr_result_1","corr_result_1",14,-3.5,3.5);
   corr_result_1->SetBinContent(4,6.53);
   corr_result_1->SetBinContent(5,6.25);
   corr_result_1->SetBinContent(6,6.09);
   corr_result_1->SetBinContent(7,5.73);
   corr_result_1->SetBinContent(8,5.73);
   corr_result_1->SetBinContent(9,6.11);
   corr_result_1->SetBinContent(10,6.29);
   corr_result_1->SetBinContent(11,6.55);
   corr_result_1->SetBinError(4,0.03429594);
   corr_result_1->SetBinError(5,0.03246259);
   corr_result_1->SetBinError(6,0.03176294);
   corr_result_1->SetBinError(7,0.03051352);
   corr_result_1->SetBinError(8,0.03040509);
   corr_result_1->SetBinError(9,0.03198735);
   corr_result_1->SetBinError(10,0.03306762);
   corr_result_1->SetBinError(11,0.03435959);
   corr_result_1->SetEntries(9.108065);
   corr_result_1->SetLineStyle(0);
   corr_result_1->SetMarkerStyle(4);
   corr_result_1->SetMarkerSize(1.8);
   corr_result_1->GetXaxis()->SetTitle("#eta");
   corr_result_1->GetXaxis()->SetLabelFont(42);
   corr_result_1->GetXaxis()->SetLabelOffset(0.007);
   corr_result_1->GetXaxis()->SetLabelSize(0.05);
   corr_result_1->GetXaxis()->SetTitleSize(0.06);
   corr_result_1->GetXaxis()->SetTitleOffset(0.9);
   corr_result_1->GetXaxis()->SetTitleFont(42);
   corr_result_1->GetYaxis()->SetTitle("dN/d#eta");
   corr_result_1->GetYaxis()->SetLabelFont(42);
   corr_result_1->GetYaxis()->SetLabelOffset(0.007);
   corr_result_1->GetYaxis()->SetLabelSize(0.05);
   corr_result_1->GetYaxis()->SetTitleSize(0.06);
   corr_result_1->GetYaxis()->SetTitleOffset(0.9);
   corr_result_1->GetYaxis()->SetTitleFont(42);
   corr_result_1->GetZaxis()->SetLabelFont(42);
   corr_result_1->GetZaxis()->SetLabelOffset(0.007);
   corr_result_1->GetZaxis()->SetLabelSize(0.05);
   corr_result_1->GetZaxis()->SetTitleSize(0.06);
   corr_result_1->GetZaxis()->SetTitleFont(42);
   return corr_result_1;
}
TH1F* tracklet900GeV()
{
   Double_t xAxis5[13] = {-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3};

   TH1F *hMeasuredFinal = new TH1F("hMeasuredFinal","",12, xAxis5);


   hMeasuredFinal->SetBinContent(7,(3.44304+3.40686)/2);
   hMeasuredFinal->SetBinContent(8,(3.61126+3.61123)/2);
   hMeasuredFinal->SetBinContent(9,(3.64045+3.77519)/2);
   hMeasuredFinal->SetBinContent(10,(3.69915+3.49038)/2);
   hMeasuredFinal->SetBinContent(6,(3.44304+3.40686)/2);
   hMeasuredFinal->SetBinContent(5,(3.61126+3.61123)/2);
   hMeasuredFinal->SetBinContent(4,(3.64045+3.77519)/2);
   hMeasuredFinal->SetBinContent(3,(3.69915+3.49038)/2);


   hMeasuredFinal->SetBinError(3,0.132634);
   hMeasuredFinal->SetBinError(4,0.143457);
   hMeasuredFinal->SetBinError(5,0.137227);
   hMeasuredFinal->SetBinError(6,0.129461);
   hMeasuredFinal->SetBinError(7,0.130835);
   hMeasuredFinal->SetBinError(8,0.137228);
   hMeasuredFinal->SetBinError(9,0.138337);
   hMeasuredFinal->SetBinError(10,0.140568);

   
   hMeasuredFinal->SetMarkerColor(kBlue);
   hMeasuredFinal->SetMarkerStyle(4);
   hMeasuredFinal->SetLineColor(4);
   hMeasuredFinal->SetMarkerSize(0.7); // use rootlogon size

   return hMeasuredFinal;
}


TH1F* tracklet2360GeVHF1()
{

   Double_t xAxis3[13] = {-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3}; 
   
   TH1F *hAvg = new TH1F("hAvg","",12, xAxis3);
   hAvg->SetBinContent(3,4.802082);
   hAvg->SetBinContent(4,4.658647);
   hAvg->SetBinContent(5,4.569177);
   hAvg->SetBinContent(6,4.405703);
   hAvg->SetBinContent(7,4.405703);
   hAvg->SetBinContent(8,4.569177);
   hAvg->SetBinContent(9,4.658647);
   hAvg->SetBinContent(10,4.802082);
   hAvg->SetBinError(3,0.1824791);
   hAvg->SetBinError(4,0.1770286);
   hAvg->SetBinError(5,0.1736287);
   hAvg->SetBinError(6,0.1674167);
   hAvg->SetBinError(7,0.1674167);
   hAvg->SetBinError(8,0.1736287);
   hAvg->SetBinError(9,0.1770286);
   hAvg->SetBinError(10,0.1824791);
   hAvg->SetMinimum(0);
   hAvg->SetMaximum(5.5);
   hAvg->SetEntries(194100.4);
   hAvg->SetFillColor(1);
   hAvg->SetFillStyle(0);
   hAvg->SetLineStyle(0);
   hAvg->SetMarkerStyle(20);
   hAvg->GetXaxis()->SetTitle("#eta");
   hAvg->GetXaxis()->CenterTitle(true);
   hAvg->GetXaxis()->SetLabelFont(42);
   hAvg->GetXaxis()->SetLabelOffset(0.01);
   hAvg->GetXaxis()->SetLabelSize(0.045);
   hAvg->GetXaxis()->SetTitleSize(0.055);
   hAvg->GetXaxis()->SetTitleFont(42);
   hAvg->GetYaxis()->SetTitle("dN/d#eta");
   hAvg->GetYaxis()->CenterTitle(true);
   hAvg->GetYaxis()->SetLabelFont(42);
   hAvg->GetYaxis()->SetLabelOffset(0.01);
   hAvg->GetYaxis()->SetLabelSize(0.045);
   hAvg->GetYaxis()->SetTitleSize(0.055);
   hAvg->GetYaxis()->SetTitleOffset(1.25);
   hAvg->GetYaxis()->SetTitleFont(42);
   hAvg->GetZaxis()->SetLabelFont(42);
   hAvg->GetZaxis()->SetLabelSize(0.045);
   hAvg->GetZaxis()->SetTitleFont(42);

   hAvg->SetMarkerColor(kBlue);
   hAvg->SetMarkerStyle(4);
   hAvg->SetLineColor(4);
   hAvg->SetMarkerSize(0.7); // use rootlogon size
   return hAvg;
}


TH1F* tracklet900GeVHF1()
{

   Double_t xAxis3[13] = {-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3}; 
   
   TH1F *hAvg = new TH1F("hAvg","",12, xAxis3);
   hAvg->SetBinContent(3,3.6899);
   hAvg->SetBinContent(4,3.617888);
   hAvg->SetBinContent(5,3.533942);
   hAvg->SetBinContent(6,3.441921);
   hAvg->SetBinContent(7,3.441921);
   hAvg->SetBinContent(8,3.533942);
   hAvg->SetBinContent(9,3.617888);
   hAvg->SetBinContent(10,3.6899);
   hAvg->SetBinError(3,0.1402162);
   hAvg->SetBinError(4,0.1374798);
   hAvg->SetBinError(5,0.1342898);
   hAvg->SetBinError(6,0.130793);
   hAvg->SetBinError(7,0.130793);
   hAvg->SetBinError(8,0.1342898);
   hAvg->SetBinError(9,0.1374798);
   hAvg->SetBinError(10,0.1402162);
   hAvg->SetMinimum(0);
   hAvg->SetMaximum(5.5);
   hAvg->SetEntries(194100.4);
   hAvg->SetFillColor(1);
   hAvg->SetFillStyle(0);
   hAvg->SetLineStyle(0);
   hAvg->SetMarkerStyle(20);
   hAvg->GetXaxis()->SetTitle("#eta");
   hAvg->GetXaxis()->CenterTitle(true);
   hAvg->GetXaxis()->SetLabelFont(42);
   hAvg->GetXaxis()->SetLabelOffset(0.01);
   hAvg->GetXaxis()->SetLabelSize(0.045);
   hAvg->GetXaxis()->SetTitleSize(0.055);
   hAvg->GetXaxis()->SetTitleFont(42);
   hAvg->GetYaxis()->SetTitle("dN/d#eta");
   hAvg->GetYaxis()->CenterTitle(true);
   hAvg->GetYaxis()->SetLabelFont(42);
   hAvg->GetYaxis()->SetLabelOffset(0.01);
   hAvg->GetYaxis()->SetLabelSize(0.045);
   hAvg->GetYaxis()->SetTitleSize(0.055);
   hAvg->GetYaxis()->SetTitleOffset(1.25);
   hAvg->GetYaxis()->SetTitleFont(42);
   hAvg->GetZaxis()->SetLabelFont(42);
   hAvg->GetZaxis()->SetLabelSize(0.045);
   hAvg->GetZaxis()->SetTitleFont(42);

   hAvg->SetMarkerColor(kBlue);
   hAvg->SetMarkerStyle(4);
   hAvg->SetLineColor(4);
   hAvg->SetMarkerSize(0.7); // use rootlogon size
   return hAvg;
}

TH1F* getPixelCounting2360GeV()
{
  TH1 *corr_result_1 = new TH1D("corr_result_1","corr_result_1",14,-3.5,3.5);
  corr_result_1->SetBinContent(4,5.1347);
  corr_result_1->SetBinContent(5,4.946645);
  corr_result_1->SetBinContent(6,4.756005);
  corr_result_1->SetBinContent(7,4.541892);
  corr_result_1->SetBinContent(8,4.505972);
  corr_result_1->SetBinContent(9,4.68336);
  corr_result_1->SetBinContent(10,4.872745);
  corr_result_1->SetBinContent(11,5.1347);
  corr_result_1->SetBinError(4,0.4279572);
  corr_result_1->SetBinError(5,0.4064926);
  corr_result_1->SetBinError(6,0.3911505);
  corr_result_1->SetBinError(7,0.3846572);
  corr_result_1->SetBinError(8,0.3827717);
  corr_result_1->SetBinError(9,0.3835248);
  corr_result_1->SetBinError(10,0.3989117);
  corr_result_1->SetBinError(11,0.4153993);
  corr_result_1->SetEntries(10);
  corr_result_1->SetLineStyle(0);
  corr_result_1->SetMarkerStyle(20);
  corr_result_1->SetMarkerSize(1.8);
  corr_result_1->GetXaxis()->SetTitle("#eta");
  corr_result_1->GetXaxis()->SetLabelFont(42);
  corr_result_1->GetXaxis()->SetLabelOffset(0.007);
  corr_result_1->GetXaxis()->SetLabelSize(0.05);
  corr_result_1->GetXaxis()->SetTitleSize(0.06);
  corr_result_1->GetXaxis()->SetTitleOffset(0.9);
  corr_result_1->GetXaxis()->SetTitleFont(42);
  corr_result_1->GetYaxis()->SetTitle("dN/d#eta");
  corr_result_1->GetYaxis()->SetLabelFont(42);
  corr_result_1->GetYaxis()->SetLabelOffset(0.007);
  corr_result_1->GetYaxis()->SetLabelSize(0.05);
  corr_result_1->GetYaxis()->SetTitleSize(0.06);
  corr_result_1->GetYaxis()->SetTitleOffset(0.9);
  corr_result_1->GetYaxis()->SetTitleFont(42);
  corr_result_1->GetZaxis()->SetLabelFont(42);
  corr_result_1->GetZaxis()->SetLabelOffset(0.007);
  corr_result_1->GetZaxis()->SetLabelSize(0.05);
  corr_result_1->GetZaxis()->SetTitleSize(0.06);
  corr_result_1->GetZaxis()->SetTitleFont(42);
//  corr_result_1->Draw("sameP");
  
  TH1 *corr_result_2 = new TH1D("corr_result_2","corr_result_2",14,-3.5,3.5);
  corr_result_2->SetBinContent(5,4.880269);
  corr_result_2->SetBinContent(6,4.706363);
  corr_result_2->SetBinContent(7,4.4455);
  corr_result_2->SetBinContent(8,4.446791);
  corr_result_2->SetBinContent(9,4.690773);
  corr_result_2->SetBinContent(10,4.84068);
  corr_result_2->SetBinError(5,0.401483);
  corr_result_2->SetBinError(6,0.3919279);
  corr_result_2->SetBinError(7,0.3663138);
  corr_result_2->SetBinError(8,0.3689177);
  corr_result_2->SetBinError(9,0.3864274);
  corr_result_2->SetBinError(10,0.4051236);
  corr_result_2->SetEntries(8);
  corr_result_2->SetLineColor(2);
  corr_result_2->SetLineStyle(0);
  corr_result_2->SetMarkerColor(2);
  corr_result_2->SetMarkerStyle(24);
  corr_result_2->SetMarkerSize(1.8);
  corr_result_2->GetXaxis()->SetTitle("#eta");
  corr_result_2->GetXaxis()->SetLabelFont(42);
  corr_result_2->GetXaxis()->SetLabelOffset(0.007);
  corr_result_2->GetXaxis()->SetLabelSize(0.05);
  corr_result_2->GetXaxis()->SetTitleSize(0.06);
  corr_result_2->GetXaxis()->SetTitleOffset(0.9);
  corr_result_2->GetXaxis()->SetTitleFont(42);
  corr_result_2->GetYaxis()->SetTitle("dN/d#eta");
  corr_result_2->GetYaxis()->SetLabelFont(42);
  corr_result_2->GetYaxis()->SetLabelOffset(0.007);
  corr_result_2->GetYaxis()->SetLabelSize(0.05);
  corr_result_2->GetYaxis()->SetTitleSize(0.06);
  corr_result_2->GetYaxis()->SetTitleOffset(0.9);
  corr_result_2->GetYaxis()->SetTitleFont(42);
  corr_result_2->GetZaxis()->SetLabelFont(42);
  corr_result_2->GetZaxis()->SetLabelOffset(0.007);
  corr_result_2->GetZaxis()->SetLabelSize(0.05);
  corr_result_2->GetZaxis()->SetTitleSize(0.06);
  corr_result_2->GetZaxis()->SetTitleFont(42);
//  corr_result_2->Draw("sameP");

  TH1 *corr_result_3 = new TH1D("corr_result_3","corr_result_3",14,-3.5,3.5);
  corr_result_3->SetBinContent(6,4.667683);
  corr_result_3->SetBinContent(7,4.347159);
  corr_result_3->SetBinContent(8,4.393262);
  corr_result_3->SetBinContent(9,4.604611);
  corr_result_3->SetBinError(6,0.38661);
  corr_result_3->SetBinError(7,0.3582893);
  corr_result_3->SetBinError(8,0.3647646);
  corr_result_3->SetBinError(9,0.3802615);
  corr_result_3->SetEntries(6);
  corr_result_3->SetLineColor(4);
  corr_result_3->SetLineStyle(0);
  corr_result_3->SetMarkerColor(4);
  corr_result_3->SetMarkerStyle(25);
  corr_result_3->SetMarkerSize(1.8);
  corr_result_3->GetXaxis()->SetTitle("#eta");
  corr_result_3->GetXaxis()->SetLabelFont(42);
  corr_result_3->GetXaxis()->SetLabelOffset(0.007);
  corr_result_3->GetXaxis()->SetLabelSize(0.05);
  corr_result_3->GetXaxis()->SetTitleSize(0.06);
  corr_result_3->GetXaxis()->SetTitleOffset(0.9);
  corr_result_3->GetXaxis()->SetTitleFont(42);
  corr_result_3->GetYaxis()->SetTitle("dN/d#eta");
  corr_result_3->GetYaxis()->SetLabelFont(42);
  corr_result_3->GetYaxis()->SetLabelOffset(0.007);
  corr_result_3->GetYaxis()->SetLabelSize(0.05);
  corr_result_3->GetYaxis()->SetTitleSize(0.06);
  corr_result_3->GetYaxis()->SetTitleOffset(0.9);
  corr_result_3->GetYaxis()->SetTitleFont(42);
  corr_result_3->GetZaxis()->SetLabelFont(42);
  corr_result_3->GetZaxis()->SetLabelOffset(0.007);
  corr_result_3->GetZaxis()->SetLabelSize(0.05);
  corr_result_3->GetZaxis()->SetTitleSize(0.06);
  corr_result_3->GetZaxis()->SetTitleFont(42);
//  corr_result_3->Draw("sameP");

  TH1F *avg = new TH1F("avg","avg",14,-3.5,3.5);

  for (int i=1;i<=14;i++)
  {
     double var1=corr_result_1->GetBinContent(i);
     double var2=corr_result_2->GetBinContent(i);
     double var3=corr_result_3->GetBinContent(i);
     
     int flag=0;
     
     if (var1!=0) flag++;
     if (var2!=0) flag++;
     if (var3!=0) flag++;
     
     double var=(var1+var2+var3);
     if (flag!=0) var/=flag;
     double varErr = var*0.04;
     avg->SetBinContent(i,var);
     avg->SetBinError(i,varErr);
  }
  avg->SetMarkerStyle(4);
  return avg;
}

TH1F *getBSC2360GeV()
{
   Double_t xAxis3[13] = {-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3}; 

   TH1F *hAvg = new TH1F("hBSC","",12, xAxis3);
   hAvg->SetBinContent(3,4.755457);
   hAvg->SetBinContent(4,4.634982);
   hAvg->SetBinContent(5,4.556972);
   hAvg->SetBinContent(6,4.410845);
   hAvg->SetBinContent(7,4.410845);
   hAvg->SetBinContent(8,4.556972);
   hAvg->SetBinContent(9,4.634982);
   hAvg->SetBinContent(10,4.755457);
   hAvg->SetBinError(3,0.2282619);
   hAvg->SetBinError(4,0.2224791);
   hAvg->SetBinError(5,0.2187346);
   hAvg->SetBinError(6,0.2117206);
   hAvg->SetBinError(7,0.2117206);
   hAvg->SetBinError(8,0.2187346);
   hAvg->SetBinError(9,0.2224791);
   hAvg->SetBinError(10,0.2282619);
   hAvg->SetMinimum(0);
   hAvg->SetMaximum(6);
   hAvg->SetEntries(25);
   hAvg->SetFillColor(1);
   hAvg->SetFillStyle(0);
   hAvg->SetLineStyle(0);
   hAvg->SetMarkerStyle(4);
   hAvg->SetMarkerSize(1.2);
   hAvg->SetLineColor(4);
   hAvg->SetMarkerColor(4);
   hAvg->GetXaxis()->SetTitle("#eta");
   hAvg->GetXaxis()->CenterTitle(true);
   hAvg->GetXaxis()->SetLabelFont(42);
   hAvg->GetXaxis()->SetLabelOffset(0.01);
   hAvg->GetXaxis()->SetLabelSize(0.045);
   hAvg->GetXaxis()->SetTitleSize(0.055);
   hAvg->GetXaxis()->SetTitleFont(42);
   hAvg->GetYaxis()->SetTitle("dN/d#eta");
   hAvg->GetYaxis()->CenterTitle(true);
   hAvg->GetYaxis()->SetLabelFont(42);
   hAvg->GetYaxis()->SetLabelOffset(0.01);
   hAvg->GetYaxis()->SetLabelSize(0.045);
   hAvg->GetYaxis()->SetTitleSize(0.055);
   hAvg->GetYaxis()->SetTitleOffset(1.25);
   hAvg->GetYaxis()->SetTitleFont(42);
   hAvg->GetZaxis()->SetLabelFont(42);
   hAvg->GetZaxis()->SetLabelSize(0.045);
   hAvg->GetZaxis()->SetTitleFont(42);
   
   return hAvg;
}


TH1F *getBSC900GeV()
{
   Double_t xAxis3[13] = {-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3}; 

   TH1F *hAvg = new TH1F("hBSC900","",12, xAxis3);
   hAvg->SetBinContent(3,3.801514);
   hAvg->SetBinContent(4,3.63721);
   hAvg->SetBinContent(5,3.52502);
   hAvg->SetBinContent(6,3.414644);
   hAvg->SetBinContent(7,3.414644);
   hAvg->SetBinContent(8,3.52502);
   hAvg->SetBinContent(9,3.63721);
   hAvg->SetBinContent(10,3.801514);
   hAvg->SetBinError(3,0.1824726);
   hAvg->SetBinError(4,0.1745861);
   hAvg->SetBinError(5,0.169201);
   hAvg->SetBinError(6,0.1639029);
   hAvg->SetBinError(7,0.1639029);
   hAvg->SetBinError(8,0.169201);
   hAvg->SetBinError(9,0.1745861);
   hAvg->SetBinError(10,0.1824726);
   hAvg->SetMinimum(0);
   hAvg->SetMaximum(6);
   hAvg->SetEntries(25);
   hAvg->SetFillColor(1);
   hAvg->SetFillStyle(0);
   hAvg->SetLineStyle(0);
   hAvg->SetLineColor(4);
   hAvg->SetMarkerStyle(4);
   hAvg->SetMarkerColor(4);
   hAvg->SetMarkerSize(1.2);
   hAvg->GetXaxis()->SetTitle("#eta");
   hAvg->GetXaxis()->CenterTitle(true);
   hAvg->GetXaxis()->SetLabelFont(42);
   hAvg->GetXaxis()->SetLabelOffset(0.01);
   hAvg->GetXaxis()->SetLabelSize(0.045);
   hAvg->GetXaxis()->SetTitleSize(0.055);
   hAvg->GetXaxis()->SetTitleFont(42);
   hAvg->GetYaxis()->SetTitle("dN/d#eta");
   hAvg->GetYaxis()->CenterTitle(true);
   hAvg->GetYaxis()->SetLabelFont(42);
   hAvg->GetYaxis()->SetLabelOffset(0.01);
   hAvg->GetYaxis()->SetLabelSize(0.045);
   hAvg->GetYaxis()->SetTitleSize(0.055);
   hAvg->GetYaxis()->SetTitleOffset(1.25);
   hAvg->GetYaxis()->SetTitleFont(42);
   hAvg->GetZaxis()->SetLabelFont(42);
   hAvg->GetZaxis()->SetLabelSize(0.045);
   hAvg->GetZaxis()->SetTitleFont(42);
   return hAvg;
}

TH1F *getHF900GeV()
{
   Double_t xAxis3[13] = {-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3}; 

   TH1F *hAvg = new TH1F("hBSC900","",12, xAxis3);
   hAvg->SetBinContent(3,3.701043);
   hAvg->SetBinContent(4,3.623568);
   hAvg->SetBinContent(5,3.538689);
   hAvg->SetBinContent(6,3.445783);
   hAvg->SetBinContent(7,3.445783);
   hAvg->SetBinContent(8,3.538689);
   hAvg->SetBinContent(9,3.623568);
   hAvg->SetBinContent(10,3.701043);
   hAvg->SetBinError(3,0.1776501);
   hAvg->SetBinError(4,0.1739312);
   hAvg->SetBinError(5,0.1698571);
   hAvg->SetBinError(6,0.1653976);
   hAvg->SetBinError(7,0.1653976);
   hAvg->SetBinError(8,0.1698571);
   hAvg->SetBinError(9,0.1739312);
   hAvg->SetBinError(10,0.1776501);
   hAvg->SetMinimum(0);
   hAvg->SetMaximum(6);
   hAvg->SetEntries(25);
   hAvg->SetFillColor(1);
   hAvg->SetFillStyle(0);
   hAvg->SetLineStyle(0);
   hAvg->SetLineColor(4);
   hAvg->SetMarkerStyle(4);
   hAvg->SetMarkerColor(4);
   hAvg->SetMarkerSize(1.2);
   hAvg->GetXaxis()->SetTitle("#eta");
   hAvg->GetXaxis()->CenterTitle(true);
   hAvg->GetXaxis()->SetLabelFont(42);
   hAvg->GetXaxis()->SetLabelOffset(0.01);
   hAvg->GetXaxis()->SetLabelSize(0.045);
   hAvg->GetXaxis()->SetTitleSize(0.055);
   hAvg->GetXaxis()->SetTitleFont(42);
   hAvg->GetYaxis()->SetTitle("dN/d#eta");
   hAvg->GetYaxis()->CenterTitle(true);
   hAvg->GetYaxis()->SetLabelFont(42);
   hAvg->GetYaxis()->SetLabelOffset(0.01);
   hAvg->GetYaxis()->SetLabelSize(0.045);
   hAvg->GetYaxis()->SetTitleSize(0.055);
   hAvg->GetYaxis()->SetTitleOffset(1.25);
   hAvg->GetYaxis()->SetTitleFont(42);
   hAvg->GetZaxis()->SetLabelFont(42);
   hAvg->GetZaxis()->SetLabelSize(0.045);
   hAvg->GetZaxis()->SetTitleFont(42);
   return hAvg;
}


TH1F* getPYTHIA900()
{
 
   TH1F *hTruthWOSelection = new TH1F("hTruthWOSelection","",12, -3,3);
   hTruthWOSelection->SetBinContent(1,3.209623);
   hTruthWOSelection->SetBinContent(2,3.278142);
   hTruthWOSelection->SetBinContent(3,3.242403);
   hTruthWOSelection->SetBinContent(4,3.095124);
   hTruthWOSelection->SetBinContent(5,2.949543);
   hTruthWOSelection->SetBinContent(6,2.804374);
   hTruthWOSelection->SetBinContent(7,2.810627);
   hTruthWOSelection->SetBinContent(8,2.944706);
   hTruthWOSelection->SetBinContent(9,3.118667);
   hTruthWOSelection->SetBinContent(10,3.227454);
   hTruthWOSelection->SetBinContent(11,3.274334);
   hTruthWOSelection->SetBinContent(12,3.201184);
   hTruthWOSelection->SetBinError(1,0.009087568);
   hTruthWOSelection->SetBinError(2,0.009184057);
   hTruthWOSelection->SetBinError(3,0.009133856);
   hTruthWOSelection->SetBinError(4,0.008924003);
   hTruthWOSelection->SetBinError(5,0.008711603);
   hTruthWOSelection->SetBinError(6,0.008494516);
   hTruthWOSelection->SetBinError(7,0.00850398);
   hTruthWOSelection->SetBinError(8,0.008704456);
   hTruthWOSelection->SetBinError(9,0.008957879);
   hTruthWOSelection->SetBinError(10,0.009112776);
   hTruthWOSelection->SetBinError(11,0.009178721);
   hTruthWOSelection->SetBinError(12,0.009075613);
   hTruthWOSelection->SetMinimum(0);
   hTruthWOSelection->SetMaximum(5.5);
   hTruthWOSelection->SetEntries(1444075);
   hTruthWOSelection->SetFillColor(1);
   hTruthWOSelection->SetFillStyle(0);
   hTruthWOSelection->SetLineColor(2);
   hTruthWOSelection->SetLineStyle(0);
   hTruthWOSelection->SetMarkerColor(2);
   hTruthWOSelection->SetMarkerStyle(20);
   hTruthWOSelection->SetMarkerSize(1.1);
   hTruthWOSelection->GetXaxis()->SetTitle("#eta");
   hTruthWOSelection->GetXaxis()->CenterTitle(true);
   hTruthWOSelection->GetXaxis()->SetLabelFont(42);
   hTruthWOSelection->GetXaxis()->SetLabelOffset(0.01);
   hTruthWOSelection->GetXaxis()->SetLabelSize(0.045);
   hTruthWOSelection->GetXaxis()->SetTitleSize(0.055);
   hTruthWOSelection->GetXaxis()->SetTitleFont(42);
   hTruthWOSelection->GetYaxis()->SetTitle("dN/d#eta");
   hTruthWOSelection->GetYaxis()->CenterTitle(true);
   hTruthWOSelection->GetYaxis()->SetLabelFont(42);
   hTruthWOSelection->GetYaxis()->SetLabelOffset(0.01);
   hTruthWOSelection->GetYaxis()->SetLabelSize(0.045);
   hTruthWOSelection->GetYaxis()->SetTitleSize(0.055);
   hTruthWOSelection->GetYaxis()->SetTitleOffset(1.25);
   hTruthWOSelection->GetYaxis()->SetTitleFont(42);
   hTruthWOSelection->GetZaxis()->SetLabelFont(42);
   hTruthWOSelection->GetZaxis()->SetLabelSize(0.045);
   hTruthWOSelection->GetZaxis()->SetTitleFont(42);
   return hTruthWOSelection;
}

TH1F* getPYTHIA2360()
{
   
   TH1F *hTruthWOSelection = new TH1F("hTruthWOSelection","",12, -3,3);
   hTruthWOSelection->SetBinContent(1,3.852361);
   hTruthWOSelection->SetBinContent(2,3.842731);
   hTruthWOSelection->SetBinContent(3,3.777578);
   hTruthWOSelection->SetBinContent(4,3.635311);
   hTruthWOSelection->SetBinContent(5,3.452649);
   hTruthWOSelection->SetBinContent(6,3.310205);
   hTruthWOSelection->SetBinContent(7,3.310737);
   hTruthWOSelection->SetBinContent(8,3.455538);
   hTruthWOSelection->SetBinContent(9,3.643952);
   hTruthWOSelection->SetBinContent(10,3.792327);
   hTruthWOSelection->SetBinContent(11,3.844986);
   hTruthWOSelection->SetBinContent(12,3.854641);
   hTruthWOSelection->SetBinError(1,0.009880511);
   hTruthWOSelection->SetBinError(2,0.009868154);
   hTruthWOSelection->SetBinError(3,0.00978414);
   hTruthWOSelection->SetBinError(4,0.009598132);
   hTruthWOSelection->SetBinError(5,0.009353889);
   hTruthWOSelection->SetBinError(6,0.009158902);
   hTruthWOSelection->SetBinError(7,0.009159638);
   hTruthWOSelection->SetBinError(8,0.009357801);
   hTruthWOSelection->SetBinError(9,0.009609533);
   hTruthWOSelection->SetBinError(10,0.009803222);
   hTruthWOSelection->SetBinError(11,0.00987105);
   hTruthWOSelection->SetBinError(12,0.009883436);
   hTruthWOSelection->SetMinimum(0);
   hTruthWOSelection->SetMaximum(5.5);
   hTruthWOSelection->SetEntries(1727327);
   hTruthWOSelection->SetFillColor(1);
   hTruthWOSelection->SetFillStyle(0);
   hTruthWOSelection->SetLineColor(2);
   hTruthWOSelection->SetLineStyle(0);
   hTruthWOSelection->SetMarkerColor(2);
   hTruthWOSelection->SetMarkerStyle(20);
   hTruthWOSelection->SetMarkerSize(1.1);
   hTruthWOSelection->GetXaxis()->SetTitle("#eta");
   hTruthWOSelection->GetXaxis()->CenterTitle(true);
   hTruthWOSelection->GetXaxis()->SetLabelFont(42);
   hTruthWOSelection->GetXaxis()->SetLabelOffset(0.01);
   hTruthWOSelection->GetXaxis()->SetLabelSize(0.045);
   hTruthWOSelection->GetXaxis()->SetTitleSize(0.055);
   hTruthWOSelection->GetXaxis()->SetTitleFont(42);
   hTruthWOSelection->GetYaxis()->SetTitle("dN/d#eta");
   hTruthWOSelection->GetYaxis()->CenterTitle(true);
   hTruthWOSelection->GetYaxis()->SetLabelFont(42);
   hTruthWOSelection->GetYaxis()->SetLabelOffset(0.01);
   hTruthWOSelection->GetYaxis()->SetLabelSize(0.045);
   hTruthWOSelection->GetYaxis()->SetTitleSize(0.055);
   hTruthWOSelection->GetYaxis()->SetTitleOffset(1.25);
   hTruthWOSelection->GetYaxis()->SetTitleFont(42);
   hTruthWOSelection->GetZaxis()->SetLabelFont(42);
   hTruthWOSelection->GetZaxis()->SetLabelSize(0.045);
   hTruthWOSelection->GetZaxis()->SetTitleFont(42);
   return hTruthWOSelection;
}

TH1F *getFinal900GeV()
{
   TH1F *hMeasuredFinal3 = new TH1F("hMeasuredFinal3","",12,-3,3);
   
   hMeasuredFinal3->GetXaxis()->SetRange(1,13);

   hMeasuredFinal3->SetMinimum(0);
   hMeasuredFinal3->SetMaximum(6.0);

   hMeasuredFinal3->SetBinContent(2,3.60278);  // -2.4 to -2.0    
   hMeasuredFinal3->SetBinContent(3,3.731399616); //-2 to -1.5 
   hMeasuredFinal3->SetBinContent(4,3.618227783);        // -1.5 to -1.0 
   hMeasuredFinal3->SetBinContent(5,3.538766473); // -1.0 to -0.5 
   hMeasuredFinal3->SetBinContent(6,3.477522368);        // -0.5 to 0

   hMeasuredFinal3->SetBinContent(7,hMeasuredFinal3->GetBinContent(6)); // 0 to 0.5 
   hMeasuredFinal3->SetBinContent(8,hMeasuredFinal3->GetBinContent(5)); // 0.5 to 1.0
   hMeasuredFinal3->SetBinContent(9,hMeasuredFinal3->GetBinContent(4)); // 1.0 to 1.5 
   hMeasuredFinal3->SetBinContent(10,hMeasuredFinal3->GetBinContent(3)); // 1.5 to 2.0 
   hMeasuredFinal3->SetBinContent(11,hMeasuredFinal3->GetBinContent(2)); // 2.0 to 2.4 


   // --- Stat error ---
   double statError= 0.01767767;
   

   for (int i=2; i<12; ++i) {
     hMeasuredFinal3->SetBinError(i,hMeasuredFinal3->GetBinContent(i)*statError);
     cout << "stat error (" << i << "): " << hMeasuredFinal3->GetBinError(i) << endl;
   }

   hMeasuredFinal3->GetYaxis()->SetTitle("dN_{ch}/d#eta");
   hMeasuredFinal3->GetXaxis()->SetTitle("#eta");
   hMeasuredFinal3->GetXaxis()->CenterTitle();
   hMeasuredFinal3->GetYaxis()->CenterTitle();


   hMeasuredFinal3->GetXaxis()->SetNdivisions(405);
   //hMeasuredFinal3->GetYaxis()->SetNdivisions(1005);
   hMeasuredFinal3->GetYaxis()->SetNdivisions(506);


   hMeasuredFinal3->SetMarkerColor(2);
   hMeasuredFinal3->SetMarkerStyle(20);
   //hMeasuredFinal3->SetMarkerSize(1.8); // use rootlogon size
   
   return hMeasuredFinal3;
}

TH1F* getFinal2360GeV()
{
   TH1F *hMeasuredFinal236 = new TH1F("hMeasuredFinal236","",12,-3,3);
   //hMeasuredFinal236->GetXaxis()->SetRange(3,10);
   hMeasuredFinal236->GetXaxis()->SetRange(1,13);
   hMeasuredFinal236->SetMarkerColor(kRed);
   hMeasuredFinal236->SetMarkerStyle(4);


   hMeasuredFinal236->SetBinContent(2,4.78296);  // -2.4 to -2.0  
   hMeasuredFinal236->SetBinContent(3,4.806189822); //-2 to -1.5
   hMeasuredFinal236->SetBinContent(4,4.662502625);        // -1.5 to -1.0 
   hMeasuredFinal236->SetBinContent(5,4.614742045); // -1.0 to -0.5 
   hMeasuredFinal236->SetBinContent(6,4.469256947);        // -0.5 to 0

   hMeasuredFinal236->SetBinContent(7,hMeasuredFinal236->GetBinContent(6)); // 0 to 0.5 
   hMeasuredFinal236->SetBinContent(8,hMeasuredFinal236->GetBinContent(5)); // 0.5 to 1.0
   hMeasuredFinal236->SetBinContent(9,hMeasuredFinal236->GetBinContent(4)); // 1.0 to 1.5 
   hMeasuredFinal236->SetBinContent(10,hMeasuredFinal236->GetBinContent(3)); // 1.5 to 2.0 
   hMeasuredFinal236->SetBinContent(11,hMeasuredFinal236->GetBinContent(2)); // 2.0 to 2.4 

   return hMeasuredFinal236;
}
