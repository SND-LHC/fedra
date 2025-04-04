#ifndef EDB_ROOT_TOOLS_H
#define EDB_ROOT_TOOLS_H

//#include "TH1.h"
#include "TMath.h"
#include "TProfile2D.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

namespace ERTools {
  // peaks search
  struct PeakInfo {
    int bin_id;
    double peak_position;
    double peak_height;
    double window_integral;
    double window_mean;
    double window_rms;
  };
  std::vector<PeakInfo> FindPeaksWithIntegral(TH1* hist, int window_size, double threshold, int from_bin, int to_bin);
  
  TH1D* get_h_var( TTree *tree, const char *var, const char *hname, double bin, const char *cut="" );
  TH2D* get_h2_var( TTree *tree, const char *var1, const char *var2, const char *hname, double bin1, double bin2 );

  // Histogram utilities
  double GetMaxBinHeight(const TH1* hist);
  TH1D*  RebinHistogram(const TH1* hist, int group);
  void   DiffProfile2D(const TProfile2D* prof1, const TProfile2D* prof2, TH2D *hDiff);
    
  // Statistical functions
  double ComputeAsymmetry(TH1* hist1, TH1* hist2);
    
  // Peak finding
  std::vector<double> FindPeaksTSpectrum(TH1* hist, int npeaks=5);
    
  // Helper functions (not exposed to interpreter)
  namespace Internal {
    bool ValidateHistogram(const TH1* hist);
  }
}

#endif
#include "TH2D.h"
