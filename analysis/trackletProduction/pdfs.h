TH1F* getMultRatio()
{
   TH1F*hData = new TH1F("hData","",61,-5,300);
   hData->SetBinContent(1,0.6450937);
   hData->SetBinContent(2,0.701944);
   hData->SetBinContent(3,0.8027075);
   hData->SetBinContent(4,0.8971947);
   hData->SetBinContent(5,1.126406);
   hData->SetBinContent(6,1.334744);
   hData->SetBinContent(7,1.471403);
   hData->SetBinContent(8,1.545206);
   hData->SetBinContent(9,1.578001);
   hData->SetBinContent(10,1.456554);
   hData->SetBinContent(11,1.269754);
   hData->SetBinContent(12,1.163332);
   hData->SetBinContent(13,1.003741);
   hData->SetBinContent(14,0.8560793);
   hData->SetBinContent(15,0.6800857);
   hData->SetBinContent(16,0.508435);
   hData->SetBinContent(17,0.4717374);
   hData->SetBinContent(18,0.3881982);
   hData->SetBinContent(19,0.3744298);
   hData->SetBinContent(20,0.3028831);
   hData->SetBinContent(21,0.3141549);
   hData->SetBinContent(22,0.2663495);
   hData->SetBinContent(23,0.2649249);
   hData->SetBinContent(24,0.251578);
   hData->SetBinContent(25,0.2151928);
   hData->SetBinContent(26,0.2914738);
   hData->SetBinContent(27,0.3616434);
   hData->SetBinError(1,0.203913);
   hData->SetBinError(2,0.0106297);
   hData->SetBinError(3,0.008664558);
   hData->SetBinError(4,0.01020254);
   hData->SetBinError(5,0.01490549);
   hData->SetBinError(6,0.02000016);
   hData->SetBinError(7,0.02454768);
   hData->SetBinError(8,0.02844274);
   hData->SetBinError(9,0.03213968);
   hData->SetBinError(10,0.0334598);
   hData->SetBinError(11,0.03322895);
   hData->SetBinError(12,0.03406088);
   hData->SetBinError(13,0.03373149);
   hData->SetBinError(14,0.03303213);
   hData->SetBinError(15,0.03167981);
   hData->SetBinError(16,0.02974478);
   hData->SetBinError(17,0.03160048);
   hData->SetBinError(18,0.03295892);
   hData->SetBinError(19,0.03688089);
   hData->SetBinError(20,0.03903735);
   hData->SetBinError(21,0.04905059);
   hData->SetBinError(22,0.05434879);
   hData->SetBinError(23,0.07223056);
   hData->SetBinError(24,0.09064875);
   hData->SetBinError(25,0.1093604);
   hData->SetBinError(26,0.1720088);
   hData->SetBinError(27,0.2627279);
   hData->SetEntries(55102);
   hData->SetFillColor(1);
   hData->SetFillStyle(0);
   hData->SetLineColor(2);
   hData->SetLineStyle(0);
   hData->SetMarkerColor(2);
   hData->SetMarkerStyle(20);
   hData->SetMarkerSize(1.25);
   hData->GetXaxis()->SetLabelFont(42);
   hData->GetXaxis()->SetLabelOffset(0.01);
   hData->GetXaxis()->SetLabelSize(0.045);
   hData->GetXaxis()->SetTitleSize(0.055);
   hData->GetXaxis()->SetTitleFont(42);
   hData->GetYaxis()->SetLabelFont(42);
   hData->GetYaxis()->SetLabelOffset(0.01);
   hData->GetYaxis()->SetLabelSize(0.045);
   hData->GetYaxis()->SetTitleSize(0.055);
   hData->GetYaxis()->SetTitleOffset(1.3);
   hData->GetYaxis()->SetTitleFont(42);
   hData->GetZaxis()->SetLabelFont(42);
   hData->GetZaxis()->SetLabelSize(0.045);
   hData->GetZaxis()->SetTitleFont(42);
   return hData;
}
