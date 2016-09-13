{
//=========Macro generated from canvas: name/title
//=========  (Tue Sep 13 23:02:33 2016) by ROOT version5.34/36
   TCanvas *name = new TCanvas("name", "title",0,0,630,604);
   gStyle->SetOptFit(1);
   gStyle->SetOptTitle(0);
   name->Range(-250,-0.008995262,1138.889,0.05097315);
   name->SetFillColor(0);
   name->SetBorderMode(0);
   name->SetBorderSize(2);
   name->SetTickx(1);
   name->SetTicky(1);
   name->SetLeftMargin(0.18);
   name->SetBottomMargin(0.15);
   name->SetFrameFillStyle(0);
   name->SetFrameBorderMode(0);
   Double_t xAxis1[7] = {0, 1, 2, 3, 4, 5, 6}; 
   
   TH1F *CLs_right_left_two = new TH1F("CLs_right_left_two","CLs_right_left_two",6, xAxis1);
   CLs_right_left_two->SetBinContent(1,77.56848);
   CLs_right_left_two->SetBinContent(2,87.56239);
   CLs_right_left_two->SetBinContent(3,92.90997);
   CLs_right_left_two->SetBinContent(4,95.07765);
   CLs_right_left_two->SetBinContent(5,97.57918);
   CLs_right_left_two->SetBinContent(6,99.9559);
    
   CLs_right_left_two->SetBinError(1,1e-17);
   CLs_right_left_two->SetBinError(2,1e-17);
   CLs_right_left_two->SetBinError(3,1e-17);
   CLs_right_left_two->SetBinError(4,1e-17);
   CLs_right_left_two->SetBinError(5,1e-17);
   CLs_right_left_two->SetBinError(6,1e-17);
   CLs_right_left_two->SetEntries(6);
   CLs_right_left_two->SetStats(0);
   CLs_right_left_two->SetLineColor(2);
   CLs_right_left_two->SetMarkerColor(2);
   CLs_right_left_two->SetMarkerStyle(22);
   CLs_right_left_two->SetMarkerSize(1.5);
   CLs_right_left_two->GetXaxis()->SetTitle("Luminosity [ fb^{-1} ]");
   CLs_right_left_two->GetXaxis()->SetBinLabel(1,"40");
   CLs_right_left_two->GetXaxis()->SetBinLabel(2,"60");
   CLs_right_left_two->GetXaxis()->SetBinLabel(3,"80");
   CLs_right_left_two->GetXaxis()->SetBinLabel(4,"100");
   CLs_right_left_two->GetXaxis()->SetBinLabel(5,"124");
   CLs_right_left_two->GetXaxis()->SetBinLabel(6,"300");
   CLs_right_left_two->GetXaxis()->SetLabelFont(42);
   CLs_right_left_two->GetXaxis()->SetLabelSize(0.045);
   CLs_right_left_two->GetXaxis()->SetTitleSize(0.05);
   CLs_right_left_two->GetXaxis()->SetTitleFont(42);
   CLs_right_left_two->GetYaxis()->SetTitle("CL_{S}");
   CLs_right_left_two->GetYaxis()->SetLabelFont(42);
   CLs_right_left_two->GetYaxis()->SetLabelSize(0.045);
   CLs_right_left_two->GetYaxis()->SetTitleSize(0.05);
   CLs_right_left_two->GetYaxis()->SetTitleOffset(1.3);
   CLs_right_left_two->GetYaxis()->SetTitleFont(42);
   CLs_right_left_two->GetZaxis()->SetLabelFont(42);
   CLs_right_left_two->GetZaxis()->SetLabelSize(0.035);
   CLs_right_left_two->GetZaxis()->SetTitleSize(0.035);
   CLs_right_left_two->GetZaxis()->SetTitleFont(42);
    gStyle->SetErrorX(0);
   CLs_right_left_two->Draw();
   TLine *line = new TLine(0,95,6,95);
   line->SetLineStyle(2);
   line->SetLineWidth(3);
   line->Draw();
   
   TLegend * leg = new TLegend(0.45,0.3,0.85,0.5);
   //TLegend * leg = new TLegend(0.45,0.2,0.85,0.5);
  // TLegend *leg = new TLegend(-2.353437e-185,-2.353437e-185,-2.353437e-185,-2.353437e-185,NULL,"brNDC");
   leg->SetBorderSize(1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = 924;
   color = new TColor(ci, 0, 0, 0, " ", 0);
   leg->SetLineColor(ci);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);

   ci = TColor::GetColor("#ffffff");
   leg->SetFillColor(ci);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("CLs_right_left_two","CL_{S} RH Model from LH Model","P");
   //TLegendEntry *entry=leg->AddEntry("CLs_right_left_two","CL_{S} RH from LH","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("TLine","95% CL_{S}","L");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   
   TPaveText *pt = new TPaveText(0.17, 0.88, 0.45, 0.96,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.03);
   TText *text = pt->AddText("DELPHES simulation                                                    #sqrt{s} = 13 TeV");
   pt->Draw();
   
   pt = new TPaveText(-2.353437e-185,-2.353437e-185,-2.353437e-185,-2.353437e-185,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->Draw();
   name->Modified();
   name->cd();
   name->SetSelected(name);

    name->SaveAs("CLs_had_20fb.pdf");
}
