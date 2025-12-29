// emtang.cpp
// Read tracks from tracks tree and estimate angle resolution

#include <iostream>
#include <cstdio>
#include <cstdint>
#include <stdexcept>
#include <cmath>
#include "TTree.h"
#include "TCut.h"
#include "TNtuple.h"
#include "EdbLog.h"
#include "EdbScanProc.h"
#include "EdbTrackFitter.h"

EdbPVRec rec;
EdbScanProc sproc;
TNtuple ntuple("ntuple", "angle resolution ntuple", "dtx:dty:dtheta");

void procTrack(EdbTrackP &t)
{
  EdbTrackFitter fitter;
  EdbTrackP t2;
  if(!fitter.SplitTrack(t,t2, 5))
  {
    std::cerr << "Track splitting failed!" << std::endl;
    return;
  }
  fitter.Fit3Pos(t);
  fitter.Fit3Pos(t2);
  float dtx = t.TX() - t2.TX();
  float dty = t.TY() - t2.TY();
  float dtheta = sqrt(dtx*dtx + dty*dty);

  ntuple.Fill(dtx, dty, dtheta);
}

void ProcTracks(EdbPVRec &rec)
{
  TFile file("emtang.root", "RECREATE");
  for(int i=0; i<rec.Ntracks(); i++)
  {
    EdbTrackP *track = rec.GetTrack(i);
    procTrack(*track);
  }
  ntuple.Write();
  file.Close();
}

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 1;
  }

  try
  {
    TCut cut="nseg==10"; // No cut, read all tracks
    int ntr = sproc.ReadTracksTree(argv[1], rec, cut);
    if(ntr <= 0)
    {
      std::cerr << "No tracks read from file " << argv[1] << std::endl;
      return 1;
    }
    std::cout << "Successfully read " << ntr << " tracks" << std::endl;
    ProcTracks(rec);
  }
  catch (const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
