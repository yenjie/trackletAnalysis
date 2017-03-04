#ifndef _SYSTEMATICS_H
#define _SYSTEMATICS_H

#include "TH1.h"
#include "TF1.h"
#include "TMath.h"

#include <string>

void th1_abs(TH1F* h) {
    for (int i=1; i<=h->GetNbinsX(); ++i)
        h->SetBinContent(i, TMath::Abs(h->GetBinContent(i)));
}

void th1_ratio_abs(TH1F* h) {
    for (int i=1; i<=h->GetNbinsX(); ++i) {
        if (h->GetBinContent(i) != 0) {
            h->SetBinContent(i, TMath::Abs(h->GetBinContent(i) - 1));
            h->SetBinError(i, h->GetBinError(i));
        } else {
            h->SetBinContent(i, 0);
            h->SetBinError(i, 0);
        }
    }
}

void th1_sqrt_sum_squares(TH1F* h1, TH1F* h2) {
    for (int i=1; i<=h1->GetNbinsX(); ++i) {
        double s1 = h1->GetBinContent(i);
        double s2 = h2->GetBinContent(i);
        double s_total = TMath::Sqrt(s1 * s1 + s2 * s2);

        double e1 = h1->GetBinError(i);
        double e2 = h2->GetBinError(i);
        double e_total = TMath::Sqrt(e1 * e1 + e2 * e2);

        h1->SetBinContent(i, s_total);
        h1->SetBinError(i, e_total);
    }
}

void th1_from_tf1(TH1F* h, TF1* f) {
    for (int i=1; i<=h->GetNbinsX(); ++i)
        if (h->GetBinContent(i) != 0)
            h->SetBinContent(i, f->Eval(h->GetBinCenter(i)));
}

class sys_var_t {
friend class total_sys_var_t;

private:
    std::string label = "";
    std::string type = "";

    std::string hist_name = "";

    TH1F* hnominal = 0;
    TH1F* hvariation = 0;

    TH1F* hdiff = 0;
    TH1F* hdiff_abs = 0;
    TH1F* hratio = 0;
    TH1F* hratio_abs = 0;

    TF1* fdiff = 0;
    TF1* fratio = 0;
    TH1F* hdiff_abs_fit = 0;
    TH1F* hratio_abs_fit = 0;

    void calc_sys();

public:
    sys_var_t(const sys_var_t& sys_var);
    sys_var_t(std::string label, std::string type, TH1F* hnominal, TH1F* hvariation);
    ~sys_var_t();

    void fit_sys(std::string diff_fit_func, std::string ratio_fit_func);
    void write();

    TH1F* get_diff_abs() {return hdiff_abs;}
    TH1F* get_ratio_abs() {return hratio_abs;}
};

sys_var_t::sys_var_t(const sys_var_t& sys_var) {
    label = sys_var.label;
    type = sys_var.type;
}

sys_var_t::sys_var_t(std::string label, std::string type, TH1F* hnominal, TH1F* hvariation) {
    this->label = label;
    this->type = type;
    this->hist_name = label + "_" + type;
    this->hnominal = (TH1F*)hnominal->Clone(Form("%s_nominal", hist_name.c_str()));
    this->hvariation = (TH1F*)hvariation->Clone(Form("%s_variation", hist_name.c_str()));

    calc_sys();
}

sys_var_t::~sys_var_t() {};

void sys_var_t::calc_sys() {
    hdiff = (TH1F*)hvariation->Clone(Form("%s_diff", hist_name.c_str()));
    hdiff->Add(hnominal, -1);
    hdiff_abs = (TH1F*)hdiff->Clone(Form("%s_diff_abs", hist_name.c_str()));
    th1_abs(hdiff_abs);

    hratio = (TH1F*)hvariation->Clone(Form("%s_ratio", hist_name.c_str()));
    hratio->Divide(hvariation, hnominal);
    hratio_abs = (TH1F*)hratio->Clone(Form("%s_ratio_abs", hist_name.c_str()));
    th1_ratio_abs(hratio_abs);
}

void sys_var_t::fit_sys(std::string diff_fit_func, std::string ratio_fit_func) {
    double range_low = hnominal->GetBinLowEdge(hnominal->FindFirstBinAbove(0.1));
    double range_high = hnominal->GetBinLowEdge(hnominal->FindLastBinAbove(0.1) + 1);

    hdiff_abs_fit = (TH1F*)hdiff_abs->Clone(Form("%s_diff_abs_fit", hist_name.c_str()));
    TF1* diff_fit = new TF1(Form("%s_diff_fit_function", hist_name.c_str()), diff_fit_func.c_str());
    diff_fit->SetRange(range_low, range_high);

    hdiff_abs->Fit(Form("%s_diff_fit_function", hist_name.c_str()), "F Q", "", range_low, range_high);
    hdiff_abs->Fit(Form("%s_diff_fit_function", hist_name.c_str()), "F Q", "", range_low, range_high);
    hdiff_abs->Fit(Form("%s_diff_fit_function", hist_name.c_str()), "F M Q", "", range_low, range_high);
    fdiff = (TF1*)hdiff_abs->GetFunction(Form("%s_diff_fit_function", hist_name.c_str()))->Clone(Form("%s_diff_fit", hist_name.c_str()));
    th1_from_tf1(hdiff_abs_fit, fdiff);

    hratio_abs_fit = (TH1F*)hratio_abs->Clone(Form("%s_ratio_abs_fit", hist_name.c_str()));
    TF1* ratio_fit = new TF1(Form("%s_ratio_fit_function", hist_name.c_str()), ratio_fit_func.c_str());
    ratio_fit->SetRange(range_low, range_high);

    hratio_abs->Fit(Form("%s_ratio_fit_function", hist_name.c_str()), "F Q", "", range_low, range_high);
    hratio_abs->Fit(Form("%s_ratio_fit_function", hist_name.c_str()), "F Q", "", range_low, range_high);
    hratio_abs->Fit(Form("%s_ratio_fit_function", hist_name.c_str()), "F M Q", "", range_low, range_high);
    fratio = (TF1*)hratio_abs->GetFunction(Form("%s_ratio_fit_function", hist_name.c_str()))->Clone(Form("%s_ratio_fit", hist_name.c_str()));
    th1_from_tf1(hratio_abs_fit, fratio);
}

void sys_var_t::write() {
    hnominal->Write("", TObject::kOverwrite);
    hvariation->Write("", TObject::kOverwrite);

    hdiff->Write("", TObject::kOverwrite);
    hdiff_abs->Write("", TObject::kOverwrite);
    hratio->Write("", TObject::kOverwrite);
    hratio_abs->Write("", TObject::kOverwrite);

    fdiff->Write("", TObject::kOverwrite);
    fratio->Write("", TObject::kOverwrite);
    hdiff_abs_fit->Write("", TObject::kOverwrite);
    hratio_abs_fit->Write("", TObject::kOverwrite);
}

class total_sys_var_t {
private:
    std::string label = "";

    TH1F* hnominal = 0;
    TH1F* hsystematics = 0;

    void add_sqrt_sum_squares(TH1F* herr);

public:
    total_sys_var_t(const total_sys_var_t& total_sys_var);
    total_sys_var_t(std::string label, TH1F* hnominal);
    ~total_sys_var_t();

    void add_sys_var(sys_var_t* sys_var, int option);
    void write();

    TH1F* get_total() {return hsystematics;}
};

total_sys_var_t::total_sys_var_t(const total_sys_var_t& total_sys_var) {
    label = total_sys_var.label;
}

total_sys_var_t::total_sys_var_t(std::string label, TH1F* hnominal) {
    this->label = label;
    this->hnominal = (TH1F*)hnominal->Clone(Form("%s_nominal", label.c_str()));
    this->hsystematics = new TH1F(Form("%s_systematics", label.c_str()), "", hnominal->GetNbinsX(), hnominal->GetXaxis()->GetXbins()->GetArray());
}

total_sys_var_t::~total_sys_var_t() {};

void total_sys_var_t::add_sqrt_sum_squares(TH1F* herr) {
    th1_sqrt_sum_squares(hsystematics, herr);
}

void total_sys_var_t::add_sys_var(sys_var_t* sys_var, int option) {
    switch (option) {
        case 0:
            add_sqrt_sum_squares(sys_var->hdiff_abs);
            break;
        case 1: {
            TH1F* htmp = (TH1F*)sys_var->hratio_abs->Clone("htmp");
            htmp->Multiply(hnominal);
            add_sqrt_sum_squares(htmp);
            htmp->Delete();
            break;
        }
        case 2:
            if (!sys_var->hdiff_abs_fit) {printf("no fit found!\n"); return;}
            add_sqrt_sum_squares(sys_var->hdiff_abs_fit);
            break;
        case 3: {
            if (!sys_var->hratio_abs_fit) {printf("no fit found!\n"); return;}
            TH1F* htmp = (TH1F*)sys_var->hratio_abs_fit->Clone("htmp");
            htmp->Multiply(hnominal);
            add_sqrt_sum_squares(htmp);
            htmp->Delete();
            break;
        }
        case 4:
            break;
        default:
            return;
    }
}

void total_sys_var_t::write() {
    hnominal->Write("", TObject::kOverwrite);
    hsystematics->Write("", TObject::kOverwrite);
}

#endif
