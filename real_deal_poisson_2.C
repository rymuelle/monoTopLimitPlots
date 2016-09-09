#include <math.h> 
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom.h"
#include "power_analysis_v2.h"

void real_deal_poisson_2()
{
    power_analysis pow_analysis("Leptonic_histograms_mu",500, "", "bkg",  105.0/124.0, 100); 
//    power_analysis pow_analysis("Leptonic_histograms_mu",500, "lep", "bkg",  105.0/124.0, 100); 
    pow_analysis.loop(4, "fb");
    pow_analysis.draw_poisson(100);
//
//    power_analysis el("Leptonic_histograms_el",1000, "", "bkg", 105/124.0, 100); 
//    el.loop(4, "fb");
//    el.draw_poisson(100);

//    power_analysis had("Hadronic_histograms",1000, "had", "bkg", 105.0/124.0, 100);     
//    had.loop(4, "fb");
//    had.draw_poisson(100);
}

