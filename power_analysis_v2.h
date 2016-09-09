#include <math.h> 
#include <iostream>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom.h"
#include <vector>

class power_analysis
{
    public:
    TRandom* random1;
    TFile * my_TFile;

    TH1F * right_ratio;
    TH1F * left_ratio;
    TH1F * full_ratio;
    TH1F * poisson_right;
    TH1F * poisson_left;
    TH1F * poisson_full;
    TH1F * d_y;
    TH1F * w_jets;
    TH1F * single_top;
    TH1F * tt_semileptonic;
    TH1F * tt_dileptonic;
    
    TH1F * CLs_right_left_one;
    TH1F * CLs_right_left_two;
    TH1F * CLs_right_full_one;
    TH1F * CLs_right_full_two;
    TH1F * CLs_left_full_one;
    TH1F * CLs_left_full_two;

    TH1F * mc_right_output;
    TH1F * mc_left_output;
    TH1F * mc_full_output;

    TH1F * mc_right_test_statistic;
    TH1F * mc_left_test_statistic;
    TH1F * mc_full_test_statistic;


    TH2F * power_curve_left;
    TH2F * power_curve_right;
    TH2F * power_curve_left_fb;
    TH2F * power_curve_right_fb;

    TF1 * poisson_left_TF1;
    TF1 * poisson_right_TF1;
    TF1 * poisson_full_TF1;

    TH1F * right_mc_dist;
    TH1F * left_mc_dist;
    TH2F * TH2F_right_mc_dist;
    TH2F * TH2F_left_mc_dist;

    TLine * five_sigma_line;
    TCanvas * my_canvas;
    
    double count_signal_right, count_bkg_right;
    double count_signal_left, count_bkg_left;
    double count_signal_full, count_bkg_full;
    double alpha_left, alpha_right;
    double beta_left, beta_right;
    double alpha_temp, beta_temp;
    double cut_line; 
    double n_statistics, n_cross_section; 
    double fraction_event_poisson_left;    
    double fraction_event_poisson_right;    
    double fraction_event_poisson_full;    
    double five_sigma;
    double fb;
    double fb_temp;
    double n_left_tmp, n_right_tmp, n_full_tmp;
    double right_per_fb, left_per_fb, full_per_fb;
    double total_right_per_fb, total_left_per_fb, total_full_per_fb;
    double test_lumi;

    int cut_bin;
    int canvas_width, canvas_height;
    int bins;
    int max_sample;
    int n_sample;
    int right_entries, left_entries, full_entries;

    std::string temp;
    std::string file_name;
    std::string file_tag;
    std::string had_tag;
};

power_analysis::power_analysis(std::string file, int max_sample_l, std::string had, std::string bkg, double xs_scale, double test_lumi_local)
{
        had_tag = had;
        gROOT->ForceStyle();
        gStyle->SetPadTickY(1);
        gStyle->SetPadTickX(1);
        //gStyle->SetPadTickX(.03);
        //gErrorIgnoreLevel = kWarning;
      int n_loop = 6; 
    test_lumi = test_lumi_local;    
    double variable_bins[] = { 0,1,2,3,4,5,6, 7};
    //CLs_right_left_one = new TH1F("CLs_right_left_one", "CLs_right_left_one", n_loop,0,n_loop);
    //CLs_right_left_two = new TH1F("CLs_right_left_two", "CLs_right_left_two", n_loop,0,n_loop);
    //CLs_right_full_one = new TH1F("CLs_right_full_one", "CLs_right_full_one", n_loop,0,n_loop);
    //CLs_right_full_two = new TH1F("CLs_right_full_two", "CLs_right_full_two", n_loop,0,n_loop);
    //CLs_left_full_one = new TH1F("CLs_left_full_one", "CLs_left_full_one", n_loop,0,n_loop);
    //CLs_left_full_two = new TH1F("CLs_left_full_two", "CLs_left_full_two", n_loop,0,n_loop);
    CLs_right_left_one = new TH1F("CLs_right_left_one", "CLs_right_left_one", n_loop,variable_bins);
    CLs_right_left_two = new TH1F("CLs_right_left_two", "CLs_right_left_two", n_loop,variable_bins);
    CLs_right_full_one = new TH1F("CLs_right_full_one", "CLs_right_full_one", n_loop,variable_bins);
    CLs_right_full_two = new TH1F("CLs_right_full_two", "CLs_right_full_two", n_loop,variable_bins);
    CLs_left_full_one = new TH1F("CLs_left_full_one", "CLs_left_full_one", n_loop,variable_bins);
    CLs_left_full_two = new TH1F("CLs_left_full_two", "CLs_left_full_two", n_loop,variable_bins);

     random1 = new TRandom; 
    random1->SetSeed(500);
    bins = 100; 
    canvas_width = 900;
    canvas_height= 700;
    canvas_width = 630;
    canvas_height= 604;
    max_sample = max_sample_l;
    cut_line = .5;

    five_sigma = 0.0000005733;
    //five_sigma = .5;
    five_sigma_line = new TLine(five_sigma, -1, five_sigma, max_sample);
    five_sigma_line->SetLineColor(3);

    file_name = file + ".root";
    //TString exp;
    //exp.Form("%i", (int)xs_scale);
      char lum2[192];
      sprintf(lum2,"%.2d",124*xs_scale);
    file_tag = file + "_fb_" + lum2; //std::to_string((int)xs_scale);
    if(had == "lep")
    {
        file_tag = "leptonic_combined";
    }

    power_curve_left = new TH2F("power_curve_left", "alpha vs. count", bins, 0, 1, bins, 0, max_sample);
    power_curve_right = new TH2F("power_curve_right", "alpha vs. count", bins, 0, 1, bins, 0, max_sample);
    power_curve_left_fb = new TH2F("power_curve_left_fb", "alpha vs. count", bins, 0, 1, bins, 0, fb);
    power_curve_right_fb = new TH2F("power_curve_right_fb", "alpha vs. count", bins, 0, 1, bins, 0, fb);

    power_curve_left->SetLineColor(1);
    power_curve_left->SetMarkerColor(1);
    power_curve_left->SetMarkerStyle(2);
    power_curve_left->SetMarkerSize(3);
    power_curve_right->SetLineColor(2);
    power_curve_right->SetMarkerColor(2);
    power_curve_right->SetMarkerStyle(2);
    power_curve_right->SetMarkerSize(3);

    my_canvas = new TCanvas("name", "title", canvas_width, canvas_height);
    my_canvas->Range(-2.25,-0.009716237,10.25,0.05435566);
    my_canvas->SetFillColor(0);
    my_canvas->SetBorderMode(0);
    my_canvas->SetBorderSize(2);
    my_canvas->SetTickx(1);
    my_canvas->SetTicky(1);
    my_canvas->SetLeftMargin(0.18);
    my_canvas->SetBottomMargin(0.15);
    my_canvas->SetFrameFillStyle(0);
    my_canvas->SetFrameBorderMode(0);
    my_canvas->SetFrameFillStyle(0);
    my_canvas->SetFrameBorderMode(0); 

    my_TFile = new TFile(file_name.c_str());
    
    TDirectory * right_TDirectory = (TDirectory*)my_TFile->Get("monotop_right");
    TDirectory * left_TDirectory = (TDirectory*)my_TFile->Get("monotop_left");
    TDirectory * full_TDirectory = (TDirectory*)my_TFile->Get("monotop_full");
    
    if(had == "had")
    {
        right_ratio = (TH1F*)right_TDirectory->Get("Hbjet_ratio_pT_lumi");
        left_ratio = (TH1F*)left_TDirectory->Get("Hbjet_ratio_pT_lumi");
        full_ratio = (TH1F*)full_TDirectory->Get("Hbjet_ratio_pT_lumi");
        right_ratio->Scale(xs_scale);
        left_ratio->Scale(xs_scale);
        full_ratio->Scale(xs_scale);

        if( bkg == "bkg")
        {
            TDirectory * d_y_TDirectory = (TDirectory*)my_TFile->Get("Drell-yan");
            d_y = (TH1F*)d_y_TDirectory->Get("Hbjet_ratio_pT_lumi");
            TDirectory * single_top_TDirectory = (TDirectory*)my_TFile->Get("single-top");
            single_top = (TH1F*)single_top_TDirectory->Get("Hbjet_ratio_pT_lumi");
            //double error; 
            //std::cout << single_top->IntegralAndError(0,52,error , "") << " " << error << std::endl;
            right_ratio->Add(d_y);
            right_ratio->Add(single_top);
            left_ratio->Add(d_y);
            left_ratio->Add(single_top);
            full_ratio->Add(d_y);
            full_ratio->Add(single_top);

            
        }

        
    } else {

        right_ratio = (TH1F*)right_TDirectory->Get("bjet_ratio_pT_lumi");
        left_ratio = (TH1F*)left_TDirectory->Get("bjet_ratio_pT_lumi");
        full_ratio = (TH1F*)full_TDirectory->Get("bjet_ratio_pT_lumi");
        if (had == "lep")
        {
            
            my_TFile2 = new TFile("Leptonic_histograms_el.root");
            TDirectory * right_TDirectory2 = (TDirectory*)my_TFile2->Get("monotop_right");
            right_ratio->Add((TH1F*)right_TDirectory2->Get("bjet_ratio_pT_lumi"));

            TDirectory * left_TDirectory2 = (TDirectory*)my_TFile2->Get("monotop_left");
            left_ratio->Add((TH1F*)left_TDirectory2->Get("bjet_ratio_pT_lumi"));

            TDirectory * full_TDirectory2 = (TDirectory*)my_TFile2->Get("monotop_full");
            full_ratio->Add((TH1F*)full_TDirectory2->Get("bjet_ratio_pT_lumi"));

        }
        right_ratio->Scale(xs_scale);
        left_ratio->Scale(xs_scale);
        full_ratio->Scale(xs_scale);

        if( bkg == "bkg")
        {
            TDirectory * d_y_TDirectory = (TDirectory*)my_TFile->Get("Drell-yan");
            d_y = (TH1F*)d_y_TDirectory->Get("bjet_ratio_pT_lumi");

            TDirectory * single_top_TDirectory = (TDirectory*)my_TFile->Get("single-top");
            single_top = (TH1F*)single_top_TDirectory->Get("bjet_ratio_pT_lumi");

            TDirectory * w_jets_TDirectory = (TDirectory*)my_TFile->Get("W+jets");
            w_jets = (TH1F*)w_jets_TDirectory->Get("bjet_ratio_pT_lumi");

            TDirectory * tt_semileptonic_TDirectory = (TDirectory*)my_TFile->Get("ttbar-semileptonic");
            tt_semileptonic = (TH1F*)tt_semileptonic_TDirectory->Get("bjet_ratio_pT_lumi");

            TDirectory * tt_dileptonic_TDirectory = (TDirectory*)my_TFile->Get("ttbar-dileptonic");
            tt_dileptonic = (TH1F*)tt_dileptonic_TDirectory->Get("bjet_ratio_pT_lumi");

            right_ratio->Add(d_y);
            right_ratio->Add(single_top);
            right_ratio->Add(w_jets);
            right_ratio->Add(tt_semileptonic);
            right_ratio->Add(tt_dileptonic);

            left_ratio->Add(d_y);
            left_ratio->Add(single_top);
            left_ratio->Add(w_jets);
            left_ratio->Add(tt_semileptonic);
            left_ratio->Add(tt_dileptonic);

            full_ratio->Add(d_y);
            full_ratio->Add(single_top);
            full_ratio->Add(w_jets);
            full_ratio->Add(tt_semileptonic);
            full_ratio->Add(tt_dileptonic);
        }

    }
        //right_ratio->Rebin(5);
        //left_ratio->Rebin(5);
        //full_ratio->Rebin(5);
        
        
        

    mc_right_output = (TH1F*)right_ratio->Clone("mc_right_output");
    mc_left_output = (TH1F*)left_ratio->Clone("mc_left_output");
    mc_full_output = (TH1F*)full_ratio->Clone("mc_full_output");
    double start_mc = -.5;
    double end_mc = .5;
    int bins_mc = 1000;
    mc_right_test_statistic = new TH1F("mc_right_test_statistic", "mc_right_test_statistic", bins_mc, start_mc,end_mc);
    mc_right_test_statistic->SetLineColor(9);
    mc_left_test_statistic = new TH1F("mc_left_test_statistic", "mc_left_test_statistic", bins_mc, start_mc,end_mc);
    mc_left_test_statistic->SetLineColor(6);
    mc_full_test_statistic = new TH1F("mc_full_test_statistic", "mc_full_test_statistic", bins_mc, start_mc,end_mc);
    mc_full_test_statistic->SetLineColor(2);


    left_entries =  left_ratio->Integral(0,left_ratio->GetSize());//left_ratio->GetEntries();
    right_entries =  right_ratio->Integral(0,right_ratio->GetSize());//right_ratio->GetEntries();
    full_entries =  full_ratio->Integral(0,full_ratio->GetSize());//full_ratio->GetEntries();
    fb = 100;

    cut_bin =  cut_line*right_ratio->GetSize(); 

    right_ratio->SetLineColor(9);
    left_ratio->SetLineColor(6);
    full_ratio->SetLineColor(2);

    temp = "left_right_distributions_" + file_tag + ".png";
    right_ratio->Draw();
    left_ratio->Draw("same");
    full_ratio->Draw("same");
    my_canvas->SaveAs(temp.c_str());

   if(had =="had")
   {
    fraction_event_poisson_right = right_ratio->Integral(cut_bin, right_ratio->GetSize())/right_ratio->Integral(0,right_ratio->GetSize());
    fraction_event_poisson_left = left_ratio->Integral(cut_bin, left_ratio->GetSize())/left_ratio->Integral(0,left_ratio->GetSize());
    fraction_event_poisson_full = full_ratio->Integral(cut_bin, full_ratio->GetSize())/full_ratio->Integral(0,full_ratio->GetSize());
   } else
   {
    fraction_event_poisson_right = right_ratio->Integral(0,cut_bin)/right_ratio->Integral(0,right_ratio->GetSize());
    fraction_event_poisson_left = left_ratio->Integral(0,cut_bin)/left_ratio->Integral(0,left_ratio->GetSize());
    fraction_event_poisson_full = full_ratio->Integral(0,cut_bin)/full_ratio->Integral(0,full_ratio->GetSize());

   } 

    right_per_fb = fraction_event_poisson_right*(double)right_entries/fb;
    full_per_fb = fraction_event_poisson_full*(double)full_entries/fb;
    left_per_fb = fraction_event_poisson_left*(double)left_entries/fb;
    total_right_per_fb = (double)right_entries/fb; 
    total_full_per_fb = (double)full_entries/fb; 
    total_left_per_fb = (double)left_entries/fb; 
    std::cout << fraction_event_poisson_right << " " << right_entries <<std::endl;
    std::cout << fraction_event_poisson_left << " " << left_entries <<std::endl;
    std::cout << fraction_event_poisson_full << " " << full_entries <<std::endl;
    //std::cout << right_per_fb << " " << fraction_event_poisson_right << std::endl;
    //std::cout << left_per_fb << " " << fraction_event_poisson_left << std::endl;

   //poisson_left_TF1 = new TF1("poisson",  "Poisson Distribution 100 fb^-1; Probability per Count;Count in Sample", 0, max_sample);

   //poisson_left_TF1 = new TF1("poisson_left_pink",  "TMath::PoissonI(x*[1],[0])/[2]", 0, max_sample);
   //poisson_right_TF1 = new TF1("poisson_right_blue",  "TMath::PoissonI(x*[1],[0])/[2]", 0, max_sample);
   //poisson_full_TF1 = new TF1("poisson_full_red",  "TMath::PoissonI(x*[1],[0])/[2]", 0, max_sample);
   poisson_left_TF1 = new TF1("poisson_left_pink",  "TMath::Poisson(x*[1],[0])/[2]", 0, max_sample);
   poisson_right_TF1 = new TF1("poisson_right_blue",  "TMath::Poisson(x*[1],[0])/[2]", 0, max_sample);
   poisson_full_TF1 = new TF1("poisson_full_red",  "TMath::Poisson(x*[1],[0])/[2]", 0, max_sample);

    poisson_left_TF1->SetLineColor(6);   
    poisson_right_TF1->SetLineColor(9);   
    poisson_full_TF1->SetLineColor(2);   
}



void power_analysis::make_poisson(double femto_barn, std::string femto)
{
    fb_temp = femto_barn; 
    n_left_tmp = left_per_fb*fb_temp;
    n_right_tmp = right_per_fb*fb_temp;
    n_full_tmp = full_per_fb*fb_temp;

    poisson_left_TF1->SetParameters(n_left_tmp,1,1);
    poisson_right_TF1->SetParameters(n_right_tmp,1,1);
    poisson_full_TF1->SetParameters(n_full_tmp,1,1);
    //std::cout << n_left_tmp << " "  << n_right_tmp << std::endl;
    //poisson_left_TF1->Draw();
    //poisson_right_TF1->Draw("same");
    //my_canvas->SaveAs("temp.png");
}


void power_analysis::draw_poisson(double femto_barn)
{
     gROOT->ProcessLine(".x tdrstyle.cc");
     gStyle->SetPadLeftMargin(.16);
     double ptLow= 0.1, ptHigh = 0.35;
      ptLow= 0.17;
      ptHigh = 0.45;
         //TPaveText *pt = new TPaveText(ptLow,0.94,ptHigh,0.98,"NDC");
         TPaveText *pt = new TPaveText(ptLow,0.88,ptHigh,0.96,"brNDC");
    
          pt->SetFillColor(0);
          pt->SetTextFont(42);
          char lum2[192];
          sprintf(lum2,"DELPHES simulation                               #sqrt{s} = 13 TeV, L = %.2i fb^{-1}",100);
          pt->AddText(lum2);
            pt->SetBorderSize(0);
            pt->SetFillStyle(0);
            pt->SetTextAlign(12);
            pt->SetTextFont(42);
            pt->SetTextSize(0.03);
    
        double offset = .47;
          TPaveText *pt3 = new TPaveText(ptLow+ offset,0.88,ptHigh+ offset,0.96,"NDC");
          pt3->SetFillColor(0);
          pt3->SetTextFont(42);
            pt3->SetBorderSize(0);
            pt3->SetFillStyle(0);
            pt3->SetTextAlign(12);
            pt3->SetTextFont(42);
          sprintf(lum2," #sqrt{s} = 13 TeV, L = %.2i fb^{-1}",100);
          //pt3->AddText(lum2);  




        make_poisson(femto_barn, "fb");
        temp = "poisson_distribution_" + file_tag + ".pdf";
        //temp = "poisson_distribution_" + file_tag + "_fb_" + std::to_string(femto_barn)  + ".png";
        double max_y = .03;
        if(poisson_left_TF1->GetParameter(0) > 500)
        {
            TH1F *temphisto = new TH1F("temphisto", "Probability vs. Count (100 fb^{-1}); Count; Probability per Count", 100, 150, 500); //
        } else 
        {
            double param_one = std::min(poisson_left_TF1->GetParameter(0),poisson_right_TF1->GetParameter(0))*.6;
            double param_two = std::max(poisson_left_TF1->GetParameter(0),poisson_right_TF1->GetParameter(0))*1.5;

            TH1F *temphisto = new TH1F("temphisto", "Probability vs. Count (100 fb^{-1}); Count; Probability per Count", 100, param_one, param_two); //
        }
        temphisto->SetTitle("");
        
        max_y = std::max(poisson_right_TF1->GetMaximum()*1.2,poisson_left_TF1->GetMaximum()*1.2);
        //`std::cout << poisson_left_TF1->GetMaximum(); 
        temphisto->GetYaxis()->SetRangeUser(0,max_y);
        temphisto->GetYaxis()->SetTitleSize(0.05);
        temphisto->GetXaxis()->SetTitleSize(0.05);
        temphisto->GetXaxis()->SetLabelFont(42);
        temphisto->GetXaxis()->SetLabelOffset(0.007);
        temphisto->GetXaxis()->SetTitleSize(0.05);
        temphisto->GetXaxis()->SetTitleOffset(1.5);
        temphisto->GetXaxis()->SetTitleFont(42);
        temphisto->GetXaxis()->SetLabelSize(0.045);
        temphisto->GetYaxis()->SetTitleOffset(1.8);
        temphisto->GetYaxis()->SetLabelOffset(0.007);
        temphisto->GetYaxis()->SetLabelFont(42);
        temphisto->GetYaxis()->SetLabelOffset(0.007);
        temphisto->GetYaxis()->SetLabelSize(0.045);
        temphisto->GetYaxis()->SetTitleSize(0.05);
        temphisto->GetYaxis()->SetTitleOffset(1.8);
        temphisto->GetYaxis()->SetTitleFont(42);

        poisson_left_TF1->SetLineColor(6);
        poisson_left_TF1->SetLineStyle(2);
        poisson_left_TF1->SetLineWidth(3);
        poisson_left_TF1->SetFillColor(6);
        poisson_left_TF1->SetFillStyle(3004);

        poisson_right_TF1->SetLineColor(9);
        poisson_right_TF1->SetLineStyle(2);
        poisson_right_TF1->SetLineWidth(3);
        poisson_right_TF1->SetFillColor(9);
        poisson_right_TF1->SetFillStyle(3005);

        poisson_full_TF1->SetLineColor(2);
        poisson_full_TF1->SetLineStyle(2);
        poisson_full_TF1->SetLineWidth(3);
        poisson_full_TF1->SetFillColor(2);
        poisson_full_TF1->SetFillStyle(3006);
        


        temphisto->SetStats(0);
        temphisto->Draw();
        pt->Draw();
        pt3->Draw();

        poisson_left_TF1->Draw("same");
        poisson_right_TF1->Draw("same");
        //poisson_full_TF1->Draw("same");

        TLine * right_center =  new TLine(poisson_right_TF1->GetParameter(0), 0, poisson_right_TF1->GetParameter(0), max_y);
        right_center->SetLineColor(9);
        right_center->SetLineStyle(2);
        right_center->SetLineWidth(3);

        TLine * full_center =  new TLine(poisson_full_TF1->GetParameter(0), 0, poisson_full_TF1->GetParameter(0), max_y);
        full_center->SetLineColor(2);
        full_center->SetLineStyle(2);
        full_center->SetLineWidth(3);

        TLine * left_center =  new TLine(poisson_left_TF1->GetParameter(0), 0, poisson_left_TF1->GetParameter(0), max_y);
        left_center->SetLineColor(6);
        left_center->SetLineStyle(2);
        left_center->SetLineWidth(3);

        right_center->Draw("same");
        //full_center->Draw("same");
        left_center->Draw("same");

         TLegend * leg = new TLegend(0.65,0.7,0.89,0.85);
        leg->SetLineColorAlpha(1,0);
        leg->SetFillColorAlpha(0,1);
        leg->AddEntry(poisson_left_TF1, "LH Model", "LF");
        leg->AddEntry(poisson_right_TF1, "RH Model", "LF");
        //leg->AddEntry(poisson_full_TF1, "Signal-full 1-jet", "LF");
        //leg->Draw("same");
        leg->Draw();
       // my_canvas->SetXTitle("Count"); 
       //
       //
        my_canvas->SaveAs(temp.c_str());
        full_center->Draw("same");
        poisson_full_TF1->Draw("same");
        leg->AddEntry(poisson_full_TF1, "M3", "LF");
        leg->Draw();

        temp = "poisson_distribution_full_" + file_tag + ".pdf";
        my_canvas->SaveAs(temp.c_str());


}

void power_analysis::return_power(std::string femto, TF1 * dist_one, TF1 * dist_two)
{
    //double count_one = (int)dist_one->GetParameter(0); 
    //double count_two = (int)dist_two->GetParameter(0); 
    double count_one = dist_one->GetParameter(0); 
    double count_two = dist_two->GetParameter(0); 

    if((count_one > count_two))
    {
        beta_left = ROOT::Math::poisson_cdf(count_two,count_two);
        alpha_right = ROOT::Math::poisson_cdf(count_two,count_one);
        beta_right = 1-  ROOT::Math::poisson_cdf(count_one-1,count_one);
        alpha_left =  1 - ROOT::Math::poisson_cdf(count_one-1,count_two);
        //std::cout << dist_two->GetName() << ":\t"<<  count_two << "\t" << dist_one->GetName() << ":\t" << count_one << std::endl;
    } else if ((count_two > count_one))
    {
        beta_left = ROOT::Math::poisson_cdf(count_one,count_one);
        alpha_right = ROOT::Math::poisson_cdf(count_one,count_two);
        beta_right = 1-  ROOT::Math::poisson_cdf(count_two-1,count_two);
        alpha_left =  1 - ROOT::Math::poisson_cdf(count_two-1,count_one);
        //std::cout << dist_one->GetName() << "\t"<<  count_one << "\t" << dist_two->GetName() << ":\t" << count_two << std::endl;
    } else 
    {
        std::cout << "0" <<  std::endl;
    }

    //beta_left = ROOT::Math::poisson_cdf((int)n_left_tmp,(int)n_left_tmp);
    //alpha_right = ROOT::Math::poisson_cdf((int)n_left_tmp,(int)n_right_tmp);
    //beta_right = 1-  ROOT::Math::poisson_cdf((int)n_right_tmp-1,(int)n_right_tmp);
    //alpha_left =  1 - ROOT::Math::poisson_cdf((int)n_right_tmp-1,(int)n_left_tmp);
    //if( alpha_right > beta_left)
    //{
    //    beta_left = 1- ROOT::Math::poisson_cdf((int)n_left_tmp-1,(int)n_left_tmp);
    //    alpha_right = 1- ROOT::Math::poisson_cdf((int)n_left_tmp-1 ,(int)n_right_tmp);

    //    beta_right = ROOT::Math::poisson_cdf((int)n_right_tmp,(int)n_right_tmp);
    //    alpha_left = ROOT::Math::poisson_cdf((int)n_right_tmp,(int)n_left_tmp);
    //} 
    //std::cout <<fb_temp << "\t" << count_one << "\t" << count_two << "\t" << alpha_right << "\t" << beta_left << "\t" << alpha_left << "\t" << beta_right <<std::endl;
    std::cout <<fb_temp << " " << count_one <<" " << count_two << "\t" << dist_one->GetName()<< "\t" << dist_two->GetName() << "\t" << alpha_right << "\t" << beta_left << "\t" << alpha_left << "\t" << beta_right <<std::endl;

}


void power_analysis::loop(int n_steps, std::string femto)
{
    //double step_size = 10;
    double step_size = pow(2.71828,log(fb)/n_steps);
    step_size = 20;
    //std::cout << step_size << std::endl;
    for(int step = 0; step < n_steps; step++)
    {
        //make_poisson(pow(step_size,step+1), "fb");
        make_poisson(40+step_size*step, "fb");
        return_power("fb", poisson_left_TF1, poisson_right_TF1);
        CLs_right_left_one->Fill(step, (1-alpha_right/beta_left)*100);
        CLs_right_left_two->Fill(step, (1-alpha_left/beta_right)*100);
    } 

        make_poisson(124, "fb");
        return_power("fb", poisson_left_TF1, poisson_right_TF1);
        CLs_right_left_one->Fill(step, (1-alpha_right/beta_left)*100);
        CLs_right_left_two->Fill(step, (1-alpha_left/beta_right)*100);
        step = step +1;
        make_poisson(300, "fb");
        return_power("fb", poisson_left_TF1, poisson_right_TF1);
        CLs_right_left_one->Fill(step, (1-alpha_right/beta_left)*100);
        CLs_right_left_two->Fill(step, (1-alpha_left/beta_right)*100);

    for(int step = 0; step < n_steps; step++)
    {
        //make_poisson(pow(step_size,step+1), "fb");
        make_poisson(40+step_size*step, "fb");
        return_power("fb", poisson_right_TF1, poisson_full_TF1);
        CLs_right_full_one->Fill(step, (1-alpha_right/beta_left)*100);
        CLs_right_full_two->Fill(step, (1-alpha_left/beta_right)*100);
    } 
        make_poisson(124, "fb");
        return_power("fb", poisson_right_TF1, poisson_full_TF1);
        CLs_right_full_one->Fill(step, (1-alpha_right/beta_left)*100);
        CLs_right_full_two->Fill(step, (1-alpha_left/beta_right)*100);
        step = step +1;
        make_poisson(300, "fb");
        return_power("fb", poisson_right_TF1, poisson_full_TF1);
        CLs_right_full_one->Fill(step, (1-alpha_right/beta_left)*100);
        CLs_right_full_two->Fill(step, (1-alpha_left/beta_right)*100);

    for(int step = 0; step < n_steps; step++)
    {
        //make_poisson(pow(step_size,step+1), "fb");
        make_poisson(40+step_size*step, "fb");
        return_power("fb", poisson_left_TF1, poisson_full_TF1);
        CLs_left_full_one->Fill(step, (1-alpha_right/beta_left)*100);
        CLs_left_full_two->Fill(step, (1-alpha_left/beta_right)*100);
    } 
        make_poisson(124, "fb");
        return_power("fb", poisson_left_TF1, poisson_full_TF1);
        CLs_left_full_one->Fill(step, (1-alpha_right/beta_left)*100);
        CLs_left_full_two->Fill(step, (1-alpha_left/beta_right)*100);

        step = step +1;
        make_poisson(300, "fb");
        return_power("fb", poisson_left_TF1, poisson_full_TF1);
        CLs_left_full_one->Fill(step, (1-alpha_right/beta_left)*100);
        CLs_left_full_two->Fill(step, (1-alpha_left/beta_right)*100);
        std::cout << " _____________" << std::endl;
        double power = 0;
        double power1 = 0;
        double power2 = 0;
        test_lumi = 1;
        while (power < 95){
            make_poisson(test_lumi, "fb");
            return_power("fb", poisson_left_TF1, poisson_right_TF1);
            power1 = (1-alpha_right/beta_left)*100;
            power2 = (1-alpha_left/beta_right)*100;
            power = min(power1,power2);
            std::cout << power << " " <<test_lumi <<std::endl;
            test_lumi = test_lumi + 1;
        }



    temp = "power_curve_" + file_tag + ".png";    
    power_curve_right->Draw();
    power_curve_left->Draw("same");
    five_sigma_line->Draw("same");
    //my_canvas->SetLogy();
    //my_canvas->SetLogx();
    my_canvas->SaveAs(temp.c_str());


    temp = "poisson_" + file_tag + ".png";    
    poisson_left_TF1->Draw();
    poisson_right_TF1->Draw("same");
    my_canvas->SaveAs(temp.c_str()); 

    draw_CLs(CLs_right_left_one, CLs_right_left_two, "right_left_", "LH Model" , "RH Model");
    draw_CLs(CLs_right_full_one, CLs_right_full_two, "right_full_", "M3" , "LH Model");
    draw_CLs(CLs_left_full_one, CLs_left_full_two, "left_full_", "RH Model" , "M3");

}

void power_analysis::draw_CLs(TH1F * plot_one, TH1F * plot_two, std::string save_name, std::string plot_string_one, std::string plot_string_two)
{
    //my_canvas->SetLogy();
     gROOT->ProcessLine(".x tdrstyle.cc");
    gStyle->SetOptTitle(0);//
     gStyle->SetPadLeftMargin(0.16);
    //gStyle->SetErrorX(1)
     double ptLow= 0.1, ptHigh = 0.35;
     ptLow= 0.17;
     ptHigh = 0.45;
      TPaveText *pt = new TPaveText(ptLow,0.88,ptHigh,0.96,"brNDC");
              pt->SetFillColor(0);
          pt->SetTextFont(42);
          pt->AddText("DELPHES simulation                                                    #sqrt{s} = 13 TeV");
            pt->SetBorderSize(0);
            pt->SetFillStyle(0);
            pt->SetTextAlign(12);
            pt->SetTextFont(42);
            pt->SetTextSize(0.03);

        double offset = .5;
          TPaveText *pt3 = new TPaveText(ptLow+ offset,0.88,ptHigh+ offset,0.96,"NDC");
          pt3->SetFillColor(0);
          pt3->SetTextFont(42);
            pt3->SetBorderSize(0);
            pt3->SetFillStyle(0);
            pt3->SetTextAlign(12);
            pt3->SetTextFont(42);
          char lum2[192];
          //sprintf(lum2," #sqrt{s} = 13 TeV, L = %.2i fb^{-1}",100);
         temp = " #sqrt{s} = 13 TeV"; // , " + plot_string_one + " vs. " + plot_string_two; 
          sprintf(lum2,temp.c_str(),100);
    //      pt3->AddText(lum2);
    plot_two->GetXaxis()->SetTitle("Luminosity [ fb^{-1} ]");
    plot_two->GetYaxis()->SetTitle("CL_{S}");
    plot_two->GetYaxis()->SetTitleSize(0.05);
    plot_two->GetXaxis()->SetTitleSize(0.05);
    plot_two->GetXaxis()->SetTitleOffset(1.);
    plot_two->GetYaxis()->SetTitleOffset(1.3);
    plot_two->GetXaxis()->SetTitleFont(42);
    plot_two->GetYaxis()->SetTitleFont(42);
    plot_two->GetXaxis()->SetLabelSize(0.045);
    plot_two->GetYaxis()->SetLabelSize(0.045);

    plot_one->GetXaxis()->SetTitle("Luminosity [ fb^{-1} ]");
    plot_one->GetYaxis()->SetTitle("CL_{S}");
    plot_one->GetYaxis()->SetTitleSize(0.05);
    plot_one->GetXaxis()->SetTitleSize(0.05);
    plot_one->GetXaxis()->SetTitleOffset(1.);
    plot_one->GetYaxis()->SetTitleOffset(1.3);
    plot_one->GetXaxis()->SetTitleFont(42);
    plot_one->GetYaxis()->SetTitleFont(42);
    plot_one->GetXaxis()->SetLabelSize(0.045);
    plot_one->GetYaxis()->SetLabelSize(0.045);
    //std::cout <<  plot_one->GetBinContent(1) << " jsjhjsj " << plot_two->GetBinContent(2) << std::endl;
    if ( plot_one->GetBinContent(1) > plot_two->GetBinContent(2))
    {
        plot_one->SetMinimum(plot_two->GetBinContent(2)-2);
    }
    if(plot_one->GetBinContent(1) > 95)
    {

    //    plot_one->SetMinimum(93);
    }    
    TLine * CLS_95 =  new TLine(0, 95, 6, 95);
    CLS_95->SetLineColor(1);
    CLS_95->SetLineWidth(3);
    CLS_95->SetLineStyle(2);

    TLine * CLS_68 =  new TLine(0, 68, 6, 68);
    CLS_68->SetLineColor(1);
    CLS_68->SetLineWidth(3);
    CLS_68->SetLineStyle(3);

    TLine * CLS_99_7 =  new TLine(0, 99.7, 6, 99.7);
    CLS_99_7->SetLineColor(1);
    CLS_99_7->SetLineWidth(3);
    CLS_99_7->SetLineStyle(2);


    plot_one->SetStats(0);
    plot_two->SetMarkerColor(2);
    plot_two->SetMarkerColor(2);
    plot_two->SetMarkerColor(2);
    plot_one->SetLineColor(4);
    plot_one->SetMarkerColor(4);
    plot_one->SetMarkerStyle(23);
    plot_one->SetMarkerSize(1.5);
    plot_one->SetBinError(1,0);
    plot_one->SetBinError(2,0);
    plot_one->SetBinError(3,0);
    plot_one->SetBinError(4,0);
    plot_one->SetBinError(5,0);
    plot_one->SetBinError(6,.00000000000001);

    plot_two->SetStats(0);
    plot_two->SetLineColor(2);
    plot_two->SetMarkerColor(2);
    plot_two->SetMarkerStyle(22);
    plot_two->SetMarkerSize(1.5);
    plot_two->SetBinError(1,0);
    plot_two->SetBinError(2,0);
    plot_two->SetBinError(3,0);
    plot_two->SetBinError(4,0);
    plot_two->SetBinError(5,0);
    plot_two->SetBinError(6,.00000000000000001);

    plot_two->GetXaxis()->SetBinLabel(1,"40");
    plot_two->GetXaxis()->SetBinLabel(2,"60");
    plot_two->GetXaxis()->SetBinLabel(3,"80");
    plot_two->GetXaxis()->SetBinLabel(4,"100");
    plot_two->GetXaxis()->SetBinLabel(5,"124");
    plot_two->GetXaxis()->SetBinLabel(6,"300");
    
    plot_one->GetXaxis()->SetBinLabel(1,"40");
    plot_one->GetXaxis()->SetBinLabel(2,"60");
    plot_one->GetXaxis()->SetBinLabel(3,"80");
    plot_one->GetXaxis()->SetBinLabel(4,"100");
    plot_one->GetXaxis()->SetBinLabel(5,"124");
    plot_one->GetXaxis()->SetBinLabel(6,"300");
    //plot_one->GetXaxis()->SetBinLabel(6,"test");

    if(plot_one->GetBinContent(5) < 60){

         TLegend * leg = new TLegend(0.23,0.65,0.53,0.85);
    } else
    {

         TLegend * leg = new TLegend(0.55,0.3,0.85,0.5);
    }
        leg->SetLineColorAlpha(1,0);
        leg->SetFillColorAlpha(0,1);
        temp = "CL_{S} " + plot_string_one + " from " + plot_string_two;
        //leg->AddEntry(plot_one, temp.c_str(), "P");
        temp = "CL_{S} " + plot_string_two + " from " + plot_string_one;
        leg->AddEntry(plot_two, temp.c_str(), "P");
        if(plot_one->GetBinContent(1) < 95){

        leg->AddEntry(CLS_95, "95\% CL_{S}", "L");
        } else {

        leg->AddEntry(CLS_99_7, "99.7\% CL_{S}", "L");
        }
    
    //    leg->AddEntry(CLS_68, "68\% CL_{S}", "L");
    

    temp =  save_name + file_tag + ".pdf";    
    //plot_one->Draw("X1");
    plot_two->Draw("X1");
        if(plot_one->GetBinContent(1) < 95){

        CLS_95->Draw("");
        } else {

        CLS_99_7->Draw("");
        }


    if(plot_one->GetBinContent(1)< 70)
    {

    //CLS_68->Draw("same");
    }

    leg->Draw();
    pt->Draw();
   pt3->Draw("same"); 
    //plot_two->Draw("same");
    my_canvas->SaveAs(temp.c_str()); 

}

double power_analysis::integral_test_one(TH1F * sample_dis)
{
    TH1F * temp_right = (TH1F*)sample_dis->Clone(); 
    TH1F * temp_left = (TH1F*)sample_dis->Clone(); 
    TH1F * temp_right = (TH1F*)sample_dis->Clone(); 
    TH1F * temp_left = (TH1F*)sample_dis->Clone(); 

    temp_right->Multiply(right_ratio);
    temp_left->Multiply(left_ratio);
    
//    std::cout << temp_left->Integral()/left_ratio->Integral() << " " << temp_right->Integral()/right_ratio->Integral() << std::endl;
    return (temp_left->Integral()/left_ratio->Integral());
}

double power_analysis::integral_test_two(TH1F * sample_dis)
{
    TH1F * temp_right = (TH1F*)sample_dis->Clone(); 
    TH1F * temp_left = (TH1F*)sample_dis->Clone(); 
    TH1F * temp_right = (TH1F*)sample_dis->Clone(); 
    TH1F * temp_left = (TH1F*)sample_dis->Clone(); 

    temp_right->Multiply(right_ratio);
    temp_left->Multiply(left_ratio);
    
    return (temp_right->Integral()/right_ratio->Integral());
}

double power_analysis::integral_test(TH1F * sample_dis)
{
    TH1F * temp_right = (TH1F*)sample_dis->Clone(); 
    TH1F * temp_left = (TH1F*)sample_dis->Clone(); 
    TH1F * temp_right = (TH1F*)sample_dis->Clone(); 
    TH1F * temp_left = (TH1F*)sample_dis->Clone(); 

    temp_right->Multiply(right_ratio);
    temp_left->Multiply(left_ratio);
    double  left, right;
    left = temp_left->Integral()/left_ratio->Integral();
    right = temp_right->Integral()/right_ratio->Integral();

    std::cout << left - right << std::endl; 
    return left - right;
}
void power_analysis::shoot_MC(TH1F * template_dist, TH1F * output_dist, double femto_barn)
{

    double n_expectation_per_fb = template_dist->Integral()/100;
    double lambda = femto_barn*n_expectation_per_fb; 
    double count = random1->Poisson(lambda);

    //std::cout << n_expectation_per_fb << " " << lambda << " " << count << std::endl;


    output_dist->Reset();

    for(int i = 0; i < count ; i++)
    {
        output_dist->Fill(template_dist->GetRandom());
    }

    //std::cout << output_dist->Integral() <<std::endl;
    

}

double power_analysis::test_statistic(TH1F * output_dist)
{

    double statistic; 

    //TH1F * temp_template_one = template_dist_one->Clone("temp_template_one");
    //TH1F * temp_template_two = template_dist_two->Clone("temp_template_two");
    TH1F * temp_template_one = right_ratio->Clone("temp_template_one");
    TH1F * temp_template_one_2 = right_ratio->Clone("temp_template_one_2");
    TH1F * temp_template_two = left_ratio->Clone("temp_template_two");
    TH1F * temp_template_two_2 = left_ratio->Clone("temp_template_two_2");

    TH1F * temp_output_one = output_dist->Clone("temp_output_one");
    TH1F * temp_output_one_2 = output_dist->Clone("temp_output_one_2");
    TH1F * temp_output_two = output_dist->Clone("temp_output_two");
    TH1F * temp_output_two_2 = output_dist->Clone("temp_output_two_2");

//    TH1F * temp_output_one = right_ratio->Clone("temp_output_one");
//    TH1F * temp_output_one_2 = right_ratio->Clone("temp_output_one_2");
//    TH1F * temp_output_two = right_ratio->Clone("temp_output_two");
//    TH1F * temp_output_two_2 = right_ratio->Clone("temp_output_two_2");

    temp_template_one_2->Multiply(temp_template_one_2);
    temp_template_two_2->Multiply(temp_template_two_2);
    temp_output_one_2->Multiply(temp_output_one_2);
    temp_output_two_2->Multiply(temp_output_two_2);

    temp_template_one->Scale(1.0/sqrt(temp_template_one_2->Integral()));
    temp_template_two->Scale(1.0/sqrt(temp_template_two_2->Integral()));
    temp_output_one->Scale(1.0/sqrt(temp_output_one_2->Integral()));
    temp_output_two->Scale(1.0/sqrt(temp_output_two_2->Integral()));

    //temp_template_one->Add(temp_template_one->Integral(),-1);
    //temp_template_two->Scale(1.0/temp_template_two->Integral());
    //temp_output_one->Scale(1.0/temp_output_one->Integral());
    //temp_output_two->Scale(1.0/temp_output_two->Integral());
    
  //  temp_template_one->Multiply(temp_template_one);  
   // temp_template_two->Multiply(temp_template_two);  
//    std::cout << "norm: " << temp_template_one->Integral() << " " << temp_template_two->Integral() << std::endl;
    

    double one, two;
    temp_output_one->Multiply(temp_template_one); 
    temp_output_two->Multiply(temp_template_two); 


    one = temp_output_one->Integral();
    two = temp_output_two->Integral();


    //temp_template_one->Multiply(temp_template_one);    
    //temp_template_two->Multiply(temp_template_two);    
    //one = one/temp_template_one_2->Integral();
    //two = two/temp_template_two_2->Integral();

    std::cout << "output: " << one << " " << two <<std::endl;

    statistic = one - two;
    std::cout << "output: " << statistic << std::endl;
    return statistic;
}

void power_analysis::fill_mc_test_statistic(TH1F * template_dist_one, TH1F * output_dist, TH1F * test_statistic, double femto_barn, int n_runs)
{
    for(int i = 0; i < n_runs; i++)
    {
        shoot_MC(template_dist_one,output_dist,femto_barn);
        test_statistic->Fill(test_statistic(output_dist));

    }
}





void power_analysis::shoot_MC(double femto_barn, int n_population)
{
    //TRandom* random1 = new TRandom; 
    //random1->SetSeed(500);
   
    double right_lambda = femto_barn*total_right_per_fb; 
    double left_lambda = femto_barn*total_left_per_fb; 
    double right_n, left_n;
    double right_statistic, left_statistic;

    TH1F * left_toy_mc = (TH1F*)left_ratio->Clone("left_toy_mc");  
    TH1F * right_toy_mc = (TH1F*)right_ratio->Clone("right_toy_mc");  

    //right_mc_dist = new TH1F("right_mc_dist", "right_mc_dist", 100, -2,2);
    right_mc_dist = new TH1F("right_mc_dist", "right_mc_dist", 1000, -.6,.6);
    right_mc_dist->SetLineColor(2);
    //left_mc_dist = new TH1F("left_mc_dist", "left_mc_dist", 100, -2,2);
    left_mc_dist = new TH1F("left_mc_dist", "left_mc_dist", 1000, -.6,.6);

    TH2F_right_mc_dist = new TH2F("TH2F_right_mc_dist", "TH2F_right_mc_dist", 100, 8, 15, 100, 8, 15);
    TH2F_right_mc_dist->SetMarkerColor(2);
    TH2F_left_mc_dist = new TH2F("TH2F_left_mc_dist", "TH2F_left_mc_dist", 100, 8, 15, 100, 8, 15);

    for(int i = 0; i< n_population; i++)
    {
        right_n = random1->Poisson(right_lambda);
        left_n = random1->Poisson(left_lambda);
        std::cout << right_n << " " << left_n << std::endl;
        std::cout << right_lambda << " " << left_lambda << std::endl;

        right_toy_mc->Reset();
        for(int j_right = 0; j_right < right_n; j_right++)
        {
            right_toy_mc->Fill(right_ratio->GetRandom());   
        }

        left_toy_mc->Reset();
        for(int j_left = 0; j_left < left_n; j_left++)
        {
            left_toy_mc->Fill(left_ratio->GetRandom());   
        }
        //std::cout << right_toy_mc->Chi2Test(right_ratio,"UW") << " " << right_toy_mc->Chi2Test(left_ratio,"UW") << std::endl; 
        //std::cout << left_toy_mc->Chi2Test(left_ratio,"UW") << " " << left_toy_mc->Chi2Test(right_ratio,"UW") << std::endl; 
//        right_statistic =  atan(log(right_toy_mc->Chi2Test(right_ratio,"UW CHI2/NDF")/right_toy_mc->Chi2Test(left_ratio,"UW CHI2/NDF")));
//        left_statistic=  atan(log(left_toy_mc->Chi2Test(right_ratio,"UW CHI2/NDF")/left_toy_mc->Chi2Test(left_ratio,"UW CHI2/NDF")));
        //right_statistic =  atan((right_toy_mc->Chi2Test(right_ratio,"UW")/right_toy_mc->Chi2Test(left_ratio,"UW")));
        //left_statistic=  atan((left_toy_mc->Chi2Test(right_ratio,"UW")/left_toy_mc->Chi2Test(left_ratio,"UW")));

    left_statistic = integral_test(left_toy_mc);
    right_statistic = integral_test(right_toy_mc); 
        right_mc_dist->Fill(right_statistic);
        left_mc_dist->Fill(left_statistic);
        //TH2F_right_mc_dist->Fill(right_toy_mc->Chi2Test(right_ratio,"UW CHI2/NDF"),right_toy_mc->Chi2Test(left_ratio,"UW CHI2/NDF"));
        TH2F_left_mc_dist->Fill(integral_test_one(left_toy_mc),integral_test_two(left_toy_mc));
        TH2F_right_mc_dist->Fill(integral_test_one(right_toy_mc),integral_test_two(right_toy_mc));
        //TH2F_left_mc_dist->Fill(left_toy_mc->Chi2Test(right_ratio,"UW CHI2/NDF"), left_toy_mc->Chi2Test(left_ratio,"UW CHI2/NDF"));


    }

    //right_toy_mc->Draw();
    //left_toy_mc->Draw("same");
    //my_canvas->SaveAs("toy_mc.png");

    //left_toy_mc->Draw();
    //my_canvas->SaveAs("test_clone.png");

}



double power_analysis::Median(const TH1F * h1) { 

   int n = h1->GetXaxis()->GetNbins();  
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray(); 
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]); 
}

void power_analysis::mc_power()
{
    right_mc_dist->Draw();
    left_mc_dist->Draw("same");
    my_canvas->SaveAs("toy_mc_dist.png"); 
    
    int bin = right_mc_dist->GetSize();
    double mean_one = 0;
    double mean_two = 0;
    int bin_one;
    int bin_two;

    for (int i = 0; i < bin; i++)
    {
        if(mean_one < .5)
        {
            mean_one = right_mc_dist->Integral(0,i)/right_mc_dist->Integral();
            bin_one = i;
        }
        if(mean_two < .5)
        {
            mean_two = left_mc_dist->Integral(0,i)/left_mc_dist->Integral();
            bin_two = i;
        }
    } 

    double mean_right = right_mc_dist->GetMean();
    double mean_left = left_mc_dist->GetMean();
    std::cout << mean_right << " " <<mean_left <<std::endl;
    int bin_right = (int)(mean_right*(double)bin/2/.6) + (int)(bin/2);
    int bin_left = (int)(mean_left*(double)bin/2/.6) + (int)(bin/2);
    bin_right = bin_one;
    bin_left = bin_two;
    std::cout << bin_right << " " << bin_left <<std::endl;

    std::cout << right_mc_dist->Integral(0,bin_right)/right_mc_dist->Integral() << "\t" << left_mc_dist->Integral(0,bin_right)/left_mc_dist->Integral() << std::endl;
    std::cout << left_mc_dist->Integral(bin_left,bin)/left_mc_dist->Integral() << "\t" << right_mc_dist->Integral(bin_left,bin)/right_mc_dist->Integral() << std::endl;
    TH2F_right_mc_dist->Draw();
    TH2F_left_mc_dist->Draw("same");
    my_canvas->SaveAs("toy_mc_dist_th2f.png"); 


}

void power_analysis::mc_power(TH1F * histo_one, TH1F * histo_two)
{
    histo_one->Draw();
    histo_two->Draw("same");
    my_canvas->SaveAs("toy_mc_dist_method2.png"); 
    int bin = histo_one->GetSize();
    double mean_one = 1;
    double mean_two = 1;
    int bin_one;
    int bin_two;

    //std::cout << "should be about 1000: "  << histo_one->Integral() << " " <<  histo_two->Integral() << std::endl;
    for (int i = bin; i > 0; i--)
    {
        if(mean_one > .5)
        {
            mean_one = histo_one->Integral(0,i)/histo_one->Integral();
            bin_one = i;
        //std::cout << "mean one: " << mean_one <<std::endl; 
        }
        if(mean_two > .5)
        {
            mean_two = histo_two->Integral(0,i)/histo_two->Integral();
            bin_two = i;
        std::cout << "mean two: " << mean_two << " " << bin_two << std::endl; 
        }
    } 
    //std::cout << "should be about .5" << histo_one->Integral(0,bin_one)<< std::endl;

    double mean_right = histo_one->GetMean();
    double mean_left = histo_two->GetMean();
    

    std::cout << mean_right << " " <<mean_left <<std::endl;
    int bin_right = (int)(mean_right*(double)bin/2/.6) + (int)(bin/2);
    int bin_left = (int)(mean_left*(double)bin/2/.6) + (int)(bin/2);
    bin_right = bin_one;
    bin_left = bin_two;
   std::cout << bin_right << " " << bin_left <<std::endl;

    if(mean_right > mean_left)
    {
    std::cout << histo_one->Integral(0,bin_right)/histo_one->Integral() << "\t" << histo_two->Integral(0,bin_right)/histo_two->Integral() << std::endl;
    std::cout << histo_two->Integral(bin_left,bin)/histo_two->Integral() << "\t" << histo_one->Integral(bin_left,bin)/histo_one->Integral() << std::endl;
    }
    if(mean_left > mean_right)
    {
    std::cout << histo_one->Integral(0,bin_right)/histo_one->Integral() << "\t" << histo_two->Integral(0,bin_right)/histo_two->Integral() << std::endl;
    std::cout << histo_two->Integral(bin_left,bin)/histo_two->Integral() << "\t" << histo_one->Integral(bin_left,bin)/histo_one->Integral() << std::endl;
        
    }

}

void power_analysis::chi2_test( double femto_barn, int n_population)
{


   
    double right_lambda = femto_barn*total_right_per_fb; 
    double left_lambda = femto_barn*total_left_per_fb; 
    double right_n, left_n;
    double right_statistic, left_statistic;

    TH1F * left_toy_mc = (TH1F*)left_ratio->Clone("left_toy_mc");  
    TH1F * right_toy_mc = (TH1F*)right_ratio->Clone("right_toy_mc");  

    //right_mc_dist = new TH1F("right_mc_dist", "right_mc_dist", 100, -2,2);
    right_mc_dist = new TH1F("right_mc_dist", "right_mc_dist", 100, -1,1);
    right_mc_dist->SetLineColor(2);
    //left_mc_dist = new TH1F("left_mc_dist", "left_mc_dist", 100, -2,2);
    left_mc_dist = new TH1F("left_mc_dist", "left_mc_dist", 100, -1,1);

    TH2F_right_mc_dist = new TH2F("TH2F_right_mc_dist", "TH2F_right_mc_dist", 100, 8, 15, 100, 8, 15);
    TH2F_right_mc_dist->SetMarkerColor(2);
    TH2F_left_mc_dist = new TH2F("TH2F_left_mc_dist", "TH2F_left_mc_dist", 100, 8, 15, 100, 8, 15);

    for(int i = 0; i< n_population; i++)
    {
        right_n = random1->Poisson(right_lambda);
        left_n = random1->Poisson(left_lambda);
        //std::cout << "n for MC: " << right_n << " " << left_n << std::endl;
        //std::cout << "lambda for MC: " << right_lambda << " " << left_lambda << std::endl;

        right_toy_mc->Reset();
        for(int j_right = 0; j_right < right_n; j_right++)
        {
            right_toy_mc->Fill(right_ratio->GetRandom());   
        }

        left_toy_mc->Reset();
        for(int j_left = 0; j_left < left_n; j_left++)
        {
            left_toy_mc->Fill(left_ratio->GetRandom());   
        }

        //std::cout << "right chi2 test: " << right_toy_mc->Chi2Test(right_ratio,"UW") << " " << right_toy_mc->Chi2Test(left_ratio,"UW") << std::endl; 
        //std::cout << "left chi2 test: " << left_toy_mc->Chi2Test(right_ratio,"UW") << " " << left_toy_mc->Chi2Test(left_ratio,"UW") << std::endl; 
        //
        //std::cout << "right chi2 test: " << right_toy_mc->Chi2Test(right_ratio,"UW")-right_toy_mc->Chi2Test(left_ratio,"UW") << std::endl; 
        //std::cout << "left chi2 test: " << left_toy_mc->Chi2Test(right_ratio,"UW")-left_toy_mc->Chi2Test(left_ratio,"UW") << std::endl; 
//        std::cout << "right chi2 test: " << (atan(right_toy_mc->Chi2Test(right_ratio,"UW")/right_toy_mc->Chi2Test(left_ratio,"UW"))-atan(3.1415/4)) << std::endl; 
//        std::cout << "left chi2 test: " << (atan(left_toy_mc->Chi2Test(right_ratio,"UW")/left_toy_mc->Chi2Test(left_ratio,"UW"))-atan(3.1415/4)) << std::endl; 

right_toy_mc->Scale(1.0/right_toy_mc->Integral());
right_ratio->Scale(1.0/right_ratio->Integral());
left_toy_mc->Scale(1.0/left_toy_mc->Integral());
left_ratio->Scale(1.0/left_ratio->Integral());


        //right_statistic = (atan(right_toy_mc->Chi2Test(right_ratio,"WW")/right_toy_mc->Chi2Test(left_ratio,"WW"))-atan(3.1415/4));
        //left_statistic = (atan(left_toy_mc->Chi2Test(right_ratio,"WW")/left_toy_mc->Chi2Test(left_ratio,"WW"))-atan(3.1415/4));
        right_statistic = (right_toy_mc->Chi2Test(right_ratio,"WW")-right_toy_mc->Chi2Test(left_ratio,"WW"));
        left_statistic = (left_toy_mc->Chi2Test(right_ratio,"WW")-left_toy_mc->Chi2Test(left_ratio,"WW"));
        //right_statistic = right_toy_mc->Chi2Test(right_ratio,"UW")-right_toy_mc->Chi2Test(left_ratio,"UW");
        //left_statistic = left_toy_mc->Chi2Test(right_ratio,"UW")-left_toy_mc->Chi2Test(left_ratio,"UW");
    //    right_statistic = right_toy_mc->AndersonDarlingTest(right_ratio)-right_toy_mc->AndersonDarlingTest(left_ratio);
     //   left_statistic = left_toy_mc->AndersonDarlingTest(right_ratio)-left_toy_mc->AndersonDarlingTest(left_ratio);
//        right_statistic = right_toy_mc->KolmogorovTest(right_ratio)-right_toy_mc->KolmogorovTest(left_ratio);
//        left_statistic = left_toy_mc->KolmogorovTest(right_ratio)-left_toy_mc->KolmogorovTest(left_ratio);
//        right_statistic = right_toy_mc->GetMean();
//        left_statistic = left_toy_mc->GetMean();

//        right_statistic =  atan(log(right_toy_mc->Chi2Test(right_ratio,"UW CHI2/NDF")/right_toy_mc->Chi2Test(left_ratio,"UW CHI2/NDF")));
//        left_statistic=  atan(log(left_toy_mc->Chi2Test(right_ratio,"UW CHI2/NDF")/left_toy_mc->Chi2Test(left_ratio,"UW CHI2/NDF")));
        //right_statistic =  atan((right_toy_mc->Chi2Test(right_ratio,"UW")/right_toy_mc->Chi2Test(left_ratio,"UW")));
        //left_statistic=  atan((left_toy_mc->Chi2Test(right_ratio,"UW")/left_toy_mc->Chi2Test(left_ratio,"UW")));

  //  left_statistic = integral_test(left_toy_mc);
  //  right_statistic = integral_test(right_toy_mc); 
        right_mc_dist->Fill(right_statistic);
        left_mc_dist->Fill(left_statistic);
  //      //TH2F_right_mc_dist->Fill(right_toy_mc->Chi2Test(right_ratio,"UW CHI2/NDF"),right_toy_mc->Chi2Test(left_ratio,"UW CHI2/NDF"));
  //      TH2F_left_mc_dist->Fill(integral_test_one(left_toy_mc),integral_test_two(left_toy_mc));
  //      TH2F_right_mc_dist->Fill(integral_test_one(right_toy_mc),integral_test_two(right_toy_mc));
  //      //TH2F_left_mc_dist->Fill(left_toy_mc->Chi2Test(right_ratio,"UW CHI2/NDF"), left_toy_mc->Chi2Test(left_ratio,"UW CHI2/NDF"));


    }

    std::cout << "middle split: " << right_mc_dist->Integral(87, 100) << " " << left_mc_dist->Integral(87, 100) <<std::endl;

//    right_toy_mc->Draw();
//    my_canvas->SaveAs("test.png");
}
