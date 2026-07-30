//-- Author :  Valeri Tioukov   11/06/2008

#include <string.h>
#include <iostream>
#include <TRint.h>
#include <TEnv.h>
#include "EdbLog.h"
#include "EdbScanProc.h"

using namespace std;

void print_help_message()
{
  cout << "\nUsage: \n\t  emalign -A=idA  -B=idB[-p=NPRE -f=NFULL -o=DATA_DIRECTORY -v=DEBUG] \n";
  cout << "\t  emalign  -set=ID [-p=NPRE -f=NFULL -o=DATA_DIRECTORY -v=DEBUG -m -env=PARFILE] \n\n";
  cout << "\t\t  idA      - id of the first piece formed as BRICK.PLATE.MAJOR.MINOR \n";
  cout << "\t\t  idB      - id of the second piece formed as BRICK.PLATE.MAJOR.MINOR \n";
  cout << "\t\t  ID       - id of the dataset formed as BRICK.PLATE.MAJOR.MINOR \n";
  cout << "\t\t  NPRE     - number of the prealignments (default is 0)\n";
  cout << "\t\t  NFULL    - number of the fullalignments (default is 0)\n";
  cout << "\t\t  DEBUG    - verbosity level: 0-print nothing, 1-errors only, 2-normal, 3-print all messages\n";
  cout << "\t\t  PARFILE  - for the new alignment: take parameters from here (default: align.rootrc)\n";
  cout << "\t\t  -m       - make the affine files starting from EdbScanSet\n";
  cout << "\t\t  -new     - use the new alignment\n";
  cout << "\t\t  -check   - check the alignment\n";
  cout << "\t\t  -readaff - read par files and update set file\n";
  cout << "\nExample: \n";
  cout << "\t  emalign -set=4554.0.1.1 -o=/scratch/BRICKS -new -v2\n";
  cout << "\n If the data location directory if not explicitly defined\n";
  cout << " the current directory will be assumed to be the brick directory \n";
  cout << "\n If the parameters file (align.rootrc) is not presented - the default \n";
  cout << " parameters will be used. After the execution them are saved into align.save.rootrc file\n";
  cout << endl;
}

void set_default(TEnv &cenv)
{
  // default parameters for the new alignment
  cenv.SetValue("fedra.align.OffsetMax", 1000.);
  cenv.SetValue("fedra.align.DZ", 250.);
  cenv.SetValue("fedra.align.DPHI", 0.02);
  cenv.SetValue("fedra.align.SigmaR", 25.);
  cenv.SetValue("fedra.align.SigmaT", 0.012);
  cenv.SetValue("fedra.align.DoFine", 1);
  cenv.SetValue("fedra.readCPcut", "eCHI2P<2.0&&s.eW>10&&eN1==1&&eN2==1&&s.Theta()>0.05&&s.Theta()<0.99");
  cenv.SetValue("fedra.readCPcutA", "1");
  cenv.SetValue("fedra.readCPcutB", "1");
  cenv.SetValue("fedra.align.SaveCouples", 1);

  cenv.SetValue("emalign.outdir", "..");
  cenv.SetValue("emalign.reportfileformat", "pdf");
  cenv.SetValue("emalign.env", "align.rootrc");
  cenv.SetValue("emalign.EdbDebugLevel", 1);
}

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    print_help_message();
    return 0;
  }

  TEnv cenv("alignenv");
  set_default(cenv);
  gEDBDEBUGLEVEL = cenv.GetValue("emalign.EdbDebugLevel", 1);
  const char *env = cenv.GetValue("emalign.env", "align.rootrc");
  const char *outdir = cenv.GetValue("emalign.outdir", "..");
  const char *reportfileformat = cenv.GetValue("emalign.reportfileformat", "pdf");

  bool do_ida = false;
  bool do_idb = false;
  bool do_new = false; // apply new alignment algorithm
  bool do_set = false;
  bool do_check = false;
  bool do_makeAff = false;
  bool do_readAff = false;
  Int_t brick = 0, plate = 0, major = 0, minor = 0;
  Int_t npre = 0, nfull = 0;
  EdbID idA, idB;

  for (int i = 1; i < argc; i++)
  {
    char *key = argv[i];

    if (!strncmp(key, "-A=", 3))
    {
      if (strlen(key) > 3)
        sscanf(key + 3, "%d.%d.%d.%d", &brick, &plate, &major, &minor);
      idA.Set(brick, plate, major, minor);
      do_ida = true;
    }
    else if (!strncmp(key, "-B=", 3))
    {
      if (strlen(key) > 3)
        sscanf(key + 3, "%d.%d.%d.%d", &brick, &plate, &major, &minor);
      idB.Set(brick, plate, major, minor);
      do_idb = true;
    }
    else if (!strncmp(key, "-set=", 5))
    {
      if (strlen(key) > 5)
        sscanf(key + 5, "%d.%d.%d.%d", &brick, &plate, &major, &minor);
      do_set = true;
    }
    else if (!strncmp(key, "-new", 4))
    {
      do_new = true;
    }
    else if (!strncmp(key, "-check", 6))
    {
      do_check = true;
    }
    else if (!strncmp(key, "-env=", 5))
    {
      if (strlen(key) > 5)
        env = key + 5;
    }
    else if (!strncmp(key, "-o=", 3))
    {
      if (strlen(key) > 3)
        outdir = key + 3;
    }
    else if (!strncmp(key, "-p=", 3))
    {
      if (strlen(key) > 3)
        npre = atoi(key + 3);
    }
    else if (!strncmp(key, "-f=", 3))
    {
      if (strlen(key) > 3)
        nfull = atoi(key + 3);
    }
    else if (!strncmp(key, "-v=", 3))
    {
      if (strlen(key) > 3)
        gEDBDEBUGLEVEL = atoi(key + 3);
    }
    else if (!strncmp(key, "-m", 2))
    {
      do_makeAff = true;
    }
    else if (!strncmp(key, "-readaff", 8))
    {
      do_readAff = true;
    }
  }

  if (!((do_ida && do_idb) || do_set))
  {
    print_help_message();
    return 0;
  }

  cenv.SetValue("emalign.env", env);
  cenv.ReadFile(cenv.GetValue("emalign.env", "align.rootrc"), kEnvLocal);
  cenv.SetValue("emalign.outdir", outdir);

  EdbScanProc sproc;
  sproc.eProcDirClient = cenv.GetValue("emalign.outdir", "..");

  if (do_ida && do_idb)
  {
    printf("\n----------------------------------------------------------------------------\n");
    printf("align  %d.%d.%d.%d and  %d.%d.%d.%d\n", idA.eBrick, idA.ePlate, idA.eMajor, idA.eMinor, idB.eBrick, idB.ePlate, idB.eMajor, idB.eMinor);
    printf("----------------------------------------------------------------------------\n\n");

    if (do_new)
    {
      EdbID id0 = idA;
      id0.ePlate = 0;
      EdbScanSet *ss = sproc.ReadScanSet(id0);
      if (ss)
      {
        EdbAffine2D aff;
        float dz = -1300;
        if (ss->GetAffP2P(idA.ePlate, idB.ePlate, aff))
          dz = ss->GetDZP2P(idA.ePlate, idB.ePlate);
        sproc.AlignNewNopar(idA, idB, cenv, &aff, dz);
      }
    }
    else
      sproc.AlignAll(idA, idB, npre, nfull);
  }
  else if (do_set)
  {
    printf("\n----------------------------------------------------------------------------\n");
    printf("align set %d.%d.%d.%d\n", brick, plate, major, minor);
    printf("----------------------------------------------------------------------------\n\n");
    cenv.WriteFile("align.save.rootrc");

    EdbID id(brick, plate, major, minor);
    EdbScanSet *ss = sproc.ReadScanSet(id);
    if (!ss)
      return 0;
    ss->Brick().SetID(brick);
    ss->MakePIDList();

    if (do_makeAff)
      sproc.MakeAFFSet(*ss);
    else if (do_readAff)
    {
      sproc.AssembleScanSet(*ss);
      sproc.WriteScanSet(id, *ss);
    }
    else if (do_new)
    {
      sproc.AlignSetNewNopar(*ss, cenv);
      sproc.MakeAlignSetSummary(id, reportfileformat);
    }
    else if (do_check)
    {
      sproc.MakeAlignSetSummary(id, reportfileformat);
    }
    else
      sproc.AlignSet(*ss, npre, nfull);
  }

  cenv.WriteFile("align.save.rootrc");
  return 1;
}
