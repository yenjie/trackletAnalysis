void draw_sys_unc(TGraph* gr, TH1* h1, TH1* h1_sys) {
    for (int i=1; i<=h1->GetNbinsX(); ++i) {
        if (h1->GetBinError(i) == 0) continue;

        double x = h1->GetBinCenter(i);
        int sys_bin = h1_sys->FindBin(x);
        double bin_width = h1->GetBinLowEdge(i+1) - h1->GetBinLowEdge(i);

        double val = h1->GetBinContent(i);
        double error = TMath::Abs(h1_sys->GetBinContent(sys_bin));

        gr->SetPoint(0, x - (bin_width/2), std::max(val - error, h1->GetMinimum()));
        gr->SetPoint(1, x + (bin_width/2), std::max(val - error, h1->GetMinimum()));
        gr->SetPoint(2, x + (bin_width/2), std::min(val + error, h1->GetMaximum()));
        gr->SetPoint(3, x - (bin_width/2), std::min(val + error, h1->GetMaximum()));

        gr->DrawClone("f");
    }
}
