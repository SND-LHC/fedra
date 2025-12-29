// em_hts.cpp
// Read HTS-like ASCII segment file (one segment per line) and write a .raw.root

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdint>
#include "EdbRun.h"
#include "EdbView.h"
#include "EdbSegment.h"

using namespace std;

struct DataPoint
{
  int pos;
  int zone;
  long long shotid;
  int isg;
  int ph;
  float tx;
  float ty;
  float x;
  float y;
  float z;
  float z1;
  float z2;
  float px;
  float py;
  int id0;
  int id1;
  int id2;
};

struct HTSView
{
  int view_id = 0;
  int nsegments = 0;
  int plate = -1;
  int side = -1; // 2: top, 1: bottom
  int col = -1;
  int row = -1;
  float x0, y0, z0, z1, z2;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  std::vector<DataPoint> segments;
};

void decodeShotID(long long shotid, uint32_t &ViewID, uint32_t &ImagerID, int16_t &col, int16_t &row)
{
  // Convert to uint32_t (assuming shotid stores the uint32_t value)
  uint32_t shotid_uint32 = static_cast<uint32_t>(shotid);

  // Extract col and row
  col = static_cast<int16_t>(shotid_uint32 & 0x0000FFFF);
  row = static_cast<int16_t>((shotid_uint32 & 0xFFFF0000) >> 16);

  // Reconstruct ShotID from col and row
  uint32_t reconstructed_shotid = ((uint32_t)(uint16_t)row << 16) | ((uint32_t)(uint16_t)col);

  // Decode ViewID and ImagerID
  const uint32_t NumberOfImager = 72;
  ViewID = reconstructed_shotid / NumberOfImager;
  ImagerID = reconstructed_shotid % NumberOfImager;
}

void calculate_HTSView_Parameters(HTSView &view)
{
  if (view.segments.empty())
    return;

  float xsum = 0.0f, ysum = 0.0f, zsum = 0.0f, z1sum = 0.0f, z2sum = 0.0f;
  float xmin = view.segments[0].x, xmax = view.segments[0].x;
  float ymin = view.segments[0].y, ymax = view.segments[0].y;
  view.plate = view.segments[0].pos / 10;
  view.side = view.segments[0].pos % 10;
  view.col = view.segments[0].id0;
  view.row = view.segments[0].id1;

  for (const auto &dp : view.segments)
  {
    xsum += dp.x;
    ysum += dp.y;
    zsum += dp.z;
    z1sum += dp.z1;
    z2sum += dp.z2;

    if (dp.x < xmin)
      xmin = dp.x;
    if (dp.x > xmax)
      xmax = dp.x;
    if (dp.y < ymin)
      ymin = dp.y;
    if (dp.y > ymax)
      ymax = dp.y;
  }

  int n = view.segments.size();
  view.nsegments = n;
  view.x0 = xsum / n;
  view.y0 = ysum / n;
  view.z0 = zsum / n;
  view.z1 = z1sum / n;
  view.z2 = z2sum / n;

  view.xmin = xmin;
  view.xmax = xmax;
  view.ymin = ymin;
  view.ymax = ymax;
  view.zmin = min(view.z1, view.z2);
  view.zmax = max(view.z1, view.z2);
}

void ConvertViewToEdbView(HTSView &view, EdbRun &run)
{
  // Here you would convert the HTSView to an EdbView and add it to the EdbRun
  // This is a placeholder function
  std::cout << "Converting View " << view.view_id << " with " << view.nsegments << " segments to EdbView." << std::endl;

  calculate_HTSView_Parameters(view);

  float zbase = (view.side == 1) ? view.zmax : view.zmin;
  EdbView *edbView = run.GetView();
  edbView->GetHeader()->SetViewID(view.view_id);
  edbView->GetHeader()->SetCoordXY(view.x0, view.y0);
  EdbAffine2D aa(*(edbView->GetHeader()->GetAffine()));
  edbView->GetHeader()->SetAffine(aa.A11(), aa.A12(), aa.A21(), aa.A22(), view.x0, view.y0);
  if (view.side == 2)
  { // top
    edbView->GetHeader()->SetCoordZ(view.zmax, view.zmin, 0, 0);
    edbView->GetHeader()->SetNframes(16, 0); // assuming 16 frames for top
  }
  else
  { // bottom
    edbView->GetHeader()->SetCoordZ(0, 0, view.zmax, view.zmin);
    edbView->GetHeader()->SetNframes(0, 16); // assuming 16 frames for bottom
  }
  edbView->GetHeader()->SetColRow(view.col, view.row);

  int cnt = 0;
  EdbSegment seg;
  for (const auto &dp : view.segments)
  {
    float x = dp.x - view.x0 + (zbase - dp.z) * dp.tx; // propagate to zbase
    float y = dp.y - view.y0 + (zbase - dp.z) * dp.ty; // propagate to zbase
    seg.Set(x, y, zbase, dp.tx, dp.ty, dp.z2 - dp.z1, view.side, dp.ph, dp.isg);
    edbView->AddSegment(&seg);
    cnt++;
  }
  run.AddView(edbView);
}

int readDataPoint(const std::string &line, DataPoint &point)
{
  std::istringstream iss(line);
  iss >> point.pos >> point.zone >> point.shotid >> point.isg >> point.ph >>
      point.tx >> point.ty >> point.x >> point.y >> point.z >> point.z1 >> point.z2 >>
      point.px >> point.py >> point.id0 >> point.id1 >> point.id2;
  if (!iss)
  {
    return -1; // Error in parsing
  }
  point.z = (point.z1 + point.z2) / 2; // mt x,y corresponds to 8-th layer z
  return 0;                            // Success
}

void CheckNewView(const DataPoint &point, const DataPoint &last_point, HTSView &current_view, EdbRun &run, int line_number)
{
  if (last_point.isg == -1)
    return; // First point, nothing to compare

  if (point.id0 != last_point.id0 || point.id1 != last_point.id1)
  {
    // New view detected
    std::cout << "line=" << line_number
              << " last_isg/isg =" << last_point.isg << " / " << point.isg
              << " last_id0/id0 = " << last_point.id0 << " / " << point.id0
              << " last_id1/id1 = " << last_point.id1 << " / " << point.id1 << std::endl;
    if (current_view.nsegments > 0)
    {
      std::cout << "View " << current_view.view_id << " has " << current_view.nsegments << " segments." << std::endl;
      ConvertViewToEdbView(current_view, run);
      current_view.view_id++;
      current_view.nsegments = 0;
      current_view.segments.clear();
    }
  }
}

int readHTSFileTGZ(const std::string &filename, EdbRun &run)
{
  FILE *input_stream = nullptr;
  bool is_pipe = false;

  // Determine input source
  if (filename.find(".tgz") != std::string::npos ||
      filename.find(".tar.gz") != std::string::npos)
  {
    // For .tgz files: use pipe
    std::string command = "tar -xzOf '" + filename + "' 2>/dev/null";
    input_stream = popen(command.c_str(), "r");
    is_pipe = true;
    if (!input_stream)
    {
      throw std::runtime_error("Failed to open pipe for: " + filename);
    }
  }
  else if (filename.find(".gz") != std::string::npos)
  {
    // For .gz files: use pipe
    std::string command = "gunzip -c '" + filename + "' 2>/dev/null";
    input_stream = popen(command.c_str(), "r");
    is_pipe = true;
    if (!input_stream)
    {
      throw std::runtime_error("Failed to open pipe for: " + filename);
    }
  }
  else
  {
    // For uncompressed files: use regular FILE*
    input_stream = fopen(filename.c_str(), "r");
    if (!input_stream)
    {
      throw std::runtime_error("Cannot open file: " + filename);
    }
  }

  // Now read from the stream (pipe or regular file)
  char *line_buffer = nullptr;
  size_t buffer_size = 0;
  ssize_t bytes_read;
  int total_lines = 0;

  HTSView current_view;
  DataPoint last_point;
  last_point.isg = -1;
  last_point.id0 = 12345678;
  last_point.id1 = 12345678;

  try
  {
    while ((bytes_read = getline(&line_buffer, &buffer_size, input_stream)) != -1)
    {
      std::string line(line_buffer);

      // Remove trailing newline if present
      if (!line.empty() && line.back() == '\n')
      {
        line.pop_back();
      }

      // Remove trailing carriage return if present
      if (!line.empty() && line.back() == '\r')
      {
        line.pop_back();
      }

      total_lines++;
      if (line.empty())
        continue;

      DataPoint point;
      if (readDataPoint(line, point) != 0)
      {
        std::cerr << "Warning: could not parse line " << total_lines
                  << " in file " << filename
                  << ": " << line << std::endl;
        continue;
      }

      CheckNewView(point, last_point, current_view, run, total_lines);

      last_point = point;
      current_view.nsegments++;
      current_view.segments.push_back(point);
    }

    // Check for read errors
    if (ferror(input_stream))
    {
      throw std::runtime_error("Error reading from " + filename);
    }
  }
  catch (...)
  {
    free(line_buffer);
    if (is_pipe)
    {
      pclose(input_stream);
    }
    else
    {
      fclose(input_stream);
    }
    throw;
  }

  free(line_buffer);

  // Finalize last view
  if (current_view.nsegments > 0)
  {
    std::cout << "View " << current_view.view_id << " has "
              << current_view.nsegments << " segments." << std::endl;
    ConvertViewToEdbView(current_view, run);
  }

  // Close the input stream properly
  int status = 0;
  if (is_pipe)
  {
    status = pclose(input_stream);
    if (status != 0)
    {
      std::cerr << "Warning: decompression command exited with status "
                << status << " for file: " << filename << std::endl;
    }
  }
  else
  {
    fclose(input_stream);
  }

  return total_lines;
}

int main(int argc, char *argv[])
{
  std::vector<std::string> input_files;
  std::string output_file;

  // Parse command line arguments
  for (int i = 1; i < argc; i++)
  {
    std::string arg = argv[i];

    if (arg == "-o" || arg == "--output")
    {
      if (i + 1 < argc)
      {
        output_file = argv[++i];
      }
      else
      {
        std::cerr << "Error: -o option requires an output filename" << std::endl;
        return 1;
      }
    }
    else if (arg == "-h" || arg == "--help")
    {
      std::cerr << "Usage: " << argv[0] << " [options] <input_files...>" << std::endl;
      std::cerr << "Options:" << std::endl;
      std::cerr << "  -o, --output <file>   Specify output ROOT file (default: auto-generated)" << std::endl;
      std::cerr << "  -h, --help            Show this help message" << std::endl;
      std::cerr << std::endl;
      std::cerr << "Examples:" << std::endl;
      std::cerr << "  " << argv[0] << " side1.tgz side2.tgz -o output.root" << std::endl;
      std::cerr << "  " << argv[0] << " data1.gz data2.gz data3.gz" << std::endl;
      std::cerr << "  " << argv[0] << " *.tgz -o combined.root" << std::endl;
      return 0;
    }
    else if (arg[0] == '-')
    {
      std::cerr << "Error: Unknown option " << arg << std::endl;
      return 1;
    }
    else
    {
      // It's an input file
      input_files.push_back(arg);
    }
  }

  // Check if we have at least one input file
  if (input_files.empty())
  {
    std::cerr << "Error: No input files specified" << std::endl;
    std::cerr << "Usage: " << argv[0] << " [options] <input_files...>" << std::endl;
    std::cerr << "Use " << argv[0] << " -h for more help" << std::endl;
    return 1;
  }

  // Generate output filename if not specified
  if (output_file.empty())
  {
    if (input_files.size() == 1)
    {
      // Single input file: use its base name
      std::string infile = input_files[0];
      auto pos = infile.find_last_of("/");
      std::string base = (pos == std::string::npos) ? infile : infile.substr(pos + 1);

      // Remove extension
      auto dot = base.find_last_of('.');
      while (dot != std::string::npos &&
             (base.substr(dot) == ".tgz" ||
              base.substr(dot) == ".gz" ||
              base.substr(dot) == ".tar.gz" ||
              base.substr(dot) == ".txt"))
      {
        base = base.substr(0, dot);
        dot = base.find_last_of('.');
      }

      output_file = base + ".raw.root";
    }
    else
    {
      // Multiple input files: use generic name
      output_file = "combined.raw.root";
    }
  }

  // Optional: Check if output file has .root extension
  if (output_file.size() < 5 || output_file.substr(output_file.size() - 5) != ".root")
  {
    std::cout << "Note: Output filename doesn't end with .root, adding it" << std::endl;
    output_file += ".root";
  }

  std::cout << "Processing " << input_files.size() << " input file(s):" << std::endl;
  for (size_t i = 0; i < input_files.size(); i++)
  {
    std::cout << "  [" << i + 1 << "] " << input_files[i] << std::endl;
  }
  std::cout << "Output file: " << output_file << std::endl;

  // Create ROOT run object
  EdbRun run(output_file.c_str(), "RECREATE");

  try
  {
    int total_records = 0;

    // Process all input files
    for (size_t i = 0; i < input_files.size(); i++)
    {
      std::cout << "\n=== Processing file " << i + 1 << "/" << input_files.size()
                << ": " << input_files[i] << " ===" << std::endl;

      int numRecords = readHTSFileTGZ(input_files[i], run);
      std::cout << "Successfully read " << numRecords << " records from "
                << input_files[i] << std::endl;

      total_records += numRecords;
    }

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Total records processed: " << total_records << std::endl;
    std::cout << "Output saved to: " << output_file << std::endl;
  }
  catch (const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    run.Close();
    return 1;
  }

  run.Close();

  return 0;
}
