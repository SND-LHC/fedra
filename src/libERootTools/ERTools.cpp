#include <vector>
#include "TSpectrum.h"
#include "ERTools.h"
#include "EdbLog.h"
NamespaceImp(ERTools);

using namespace ERTools;

std::vector<PeakInfo> ERTools::FindPeaksWithIntegral(TH1* hist, int window_size, double threshold = 0.1, 
					    int from_bin=0, int to_bin=0)
{
    std::vector<PeakInfo> peaks;
    
    if (!hist || window_size < 1) return peaks;
    
    const int nBins = hist->GetNbinsX();
    window_size = std::min(window_size, nBins); // Ensure window isn't larger than histogram

    threshold *= hist->GetBinContent(hist->GetMaximumBin());
    //printf("threshold = %f ,   maximum %f     maxbin %f\n", 	   threshold, hist->GetMaximum(), hist->GetBinContent(hist->GetMaximumBin()) );
    
    if(from_bin==to_bin&&from_bin==0) {from_bin=1; to_bin=nBins;}
    
    for (int i = from_bin; i <= to_bin; ++i) {
        double centerContent = hist->GetBinContent(i);
        if (centerContent <= threshold) continue;
        
        bool isPeak = true;
        int left = std::max(1, i - window_size/2);
        int right = std::min(nBins, i + window_size/2);
        
        // Check if current bin is maximum in the window
        for (int j = left; j <= right; ++j) {
            if (j == i) continue;
            if (hist->GetBinContent(j) >= centerContent) {
                isPeak = false;
                break;
            }
        }
        
        if (isPeak) {
            // Calculate window statistics
            double sum = 0.0;
            double sum2 = 0.0;
	    double xsum = 0.0;
            int count = 0;
            
            for (int j = left; j <= right; ++j) {
                double val = hist->GetBinContent(j);
                sum += val;
                sum2 += val * val;
		double x = hist->GetBinCenter(j);
		xsum += x*val;
                count++;
            }
            
            double mean = sum / count;
            double rms = TMath::Sqrt(sum2/count - mean*mean);

	    double xmean = xsum/sum;
            
            peaks.push_back({
                i,                          // bin_id
		xmean,			    // peak_position
                centerContent,              // peak_height
                sum,                        // window_integral
                mean,                       // window_mean
                rms                         // window_rms
            });
        }
    }
    return peaks;
}

//----------------------------------------------------------------------------------------
void ERTools::DiffProfile2D(const TProfile2D* prof1, const TProfile2D* prof2, TH2D *hDiff)
{
  for (int binx = 1; binx <= prof1->GetNbinsX(); ++binx) {
    for (int biny = 1; biny <= prof1->GetNbinsY(); ++biny) {
      int bin = prof1->GetBin(binx, biny);
      // Skip empty bins in either histogram
      if (prof1->GetBinEntries(bin) == 0 || prof2->GetBinEntries(bin) == 0) {
	hDiff->SetBinContent(binx, biny, 0);  // Mark as 0 (or kNaN)
	continue;
      }
      double diff = prof1->GetBinContent(bin) - prof2->GetBinContent(bin);
      hDiff->SetBinContent(binx, biny, diff);
    }
  }
}

//----------------------------------------------------------------------------------------
TH1D* ERTools::get_h_var( TTree *tree, const char *var, const char *hname, double bin, const char *cut  )
{
  Log(3,"ERTools::get_h_var","%s %s %f %s", var,hname,bin,cut);
  if(!tree) return 0;
  
  // First draw with cut to get selected entries
  Long64_t nentries = tree->Draw(var, cut, "goff");
  
  if (nentries == 0) {
    Log(3,"ERTools::get_h_var","No entries pass the cut");
    return 0;
  }
  
  // Get the drawn values array
  Double_t *v = tree->GetV1();
  
  // Calculate min/max from the selected entries
  double min_ = 1e30;
  double max_ = -1e30;
  for (Long64_t i = 0; i < nentries; i++) {
    if (v[i] < min_) min_ = v[i];
    if (v[i] > max_) max_ = v[i];
  }
  
  int n = (max_-min_)/bin;
  Log(3,"ERTools::get_h_var","%d %f %f", n,min_,max_);
  
  TH1D *h = new TH1D(hname,var,n,min_,max_);
  tree->Draw(Form("%s>>%s",var,hname), cut, "goff");
  return h;
}

//----------------------------------------------------------------------------------------
TH2D* ERTools::get_h2_var( TTree *tree, const char *var1, const char *var2, const char *hname, double bin1, double bin2 )
{
  if(!tree) return 0;
  tree->Draw(var1, "", "goff");
  double  min1_ = tree->GetMinimum(var1);
  double  max1_ = tree->GetMaximum(var1);
  tree->Draw(var2, "", "goff");
  double  min2_ = tree->GetMinimum(var2);
  double  max2_ = tree->GetMaximum(var2);
  int n1 = (max1_-min1_)/bin1;
  int n2 = (max2_-min2_)/bin2;
  TH2D *h = new TH2D(hname,Form("%s vs %s",var1,var2),n1,min1_,max1_, n2,min2_,max2_);
  tree->Draw(Form("%s:%s>>%s",var1,var2,hname), "", "goff");
  return h;
}

double ERTools::GetMaxBinHeight(const TH1* hist) {
  if (!Internal::ValidateHistogram(hist)) return 0;
  return hist->GetBinContent(hist->GetMaximumBin());
}

TH1D* ERTools::RebinHistogram(const TH1* hist, int group) {
  if (!Internal::ValidateHistogram(hist)) return nullptr;
  TH1D* hnew = dynamic_cast<TH1D*>(hist->Clone());
  hnew->Rebin(group);
  return hnew;
}

std::vector<double> ERTools::FindPeaksTSpectrum(TH1* hist, int npeaks) {
  std::vector<double> peaks;
  TSpectrum spec(npeaks);
  spec.Search(hist, 2, "goff", 0.1);
    
  for (int i=0; i<spec.GetNPeaks(); i++) {
    peaks.push_back(spec.GetPositionX()[i]);
  }
  return peaks;
}

bool ERTools::Internal::ValidateHistogram(const TH1* hist) {
  return hist && hist->GetEntries() > 0;
}
