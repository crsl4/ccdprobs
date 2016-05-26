#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>

#include "parameter.h"

using namespace std;

void usage(ostream& f)
{
  f << "Usage: bl [options]" << endl;
  f << "Options:" << endl;
  f << "    -f sequence-file-name          |  the FASTA format file name with the sequence data" << endl;
  f << "    -t topology-string             |  the quoted string with the tree topology" << endl;
  f << "    -p stationary-distribution     |  four relative probabilities for A,C,G,T, comma-separated, no spaces" << endl;
  f << "    -q symmetric-q-parameters      |  six relative values for AC,AG,AT,CG,CT,GT, comma-separated, no spaces" << endl;
  f << "    -s seed                        |  positive integer random seed (machine chosen if not provided)" << endl;
  f << "    -n sample-size                 |  positive integer sample size (1000 if not provided)" << endl;
  f << "    --verbose                      |  print comments" << endl;
  f << "    --mvnormal                     |  use multivariate normal instead of gamma/beta" << endl;
  f << "    -h || --help                   |  print this help message" << endl;
  exit(1);
}

vector<double> getDoublesFromCommaList (string x)
{
  stringstream s(x);
  vector<double> v;
  double d;
  s >> d;
  v.push_back(d);
  do
  {
    char c;
    s >> c;
    if ( c==EOF )
      break;
    if ( c != ',' )
    {
      cerr << "Error: expected `,' when reading in list of numbers, read `" << c << "' instead." << endl;
      exit(1);
    }
    s >> d;
    v.push_back(d);
  } while (s.good());
  return v;
}

void Parameter::processCommandLine(int argc,char* argv[])
{
  if ( argc==1 )
    usage(cerr);

  int k=0;
  verbose = false; //false by default
  mvnormal = false;
  while ( ++k < argc )
  {
    if ( strcmp(argv[k],"-h") == 0 || strcmp(argv[k],"--help") == 0 )
      usage(cerr);

    if ( strcmp(argv[k],"-f") == 0 )
    {
      if ( ++k < argc )
      {
        stringstream s;
        s << argv[k];
        s >> sequenceFileName;
      }
      else {
        cerr << "Error: flag `-f' not followed by a filename" << endl;
        usage(cerr);
      }
    }
    else if ( strcmp(argv[k],"-t") == 0 )
    {
      if ( ++k < argc )
      {
        stringstream s;
        s << argv[k];
        s >> topology;
      }
      else
      {
        cerr << "Error: flag `-t' not followed by a topology string" << endl;
        usage(cerr);
      }
    }
    else if ( strcmp(argv[k],"-s") == 0 )
    {
      if ( ++k < argc )
      {
        stringstream s;
        s << argv[k];
        s >> seed;
      }
      else
      {
        cerr << "Error: flag `-s' not followed by a int seed" << endl;
        usage(cerr);
      }
    }
    else if ( strcmp(argv[k],"-n") == 0 )
    {
      if ( ++k < argc )
      {
        stringstream s;
        s << argv[k];
        s >> sampleSize;
      }
      else
      {
        cerr << "Error: flag `-n' not followed by a int sample size" << endl;
        usage(cerr);
      }
    }
    else if ( strcmp(argv[k],"--verbose") == 0 )
    {
      verbose = true;
    }
    else if ( strcmp(argv[k],"--mvnormal") == 0 )
    {
      mvnormal = true;
    }
    else if ( strcmp(argv[k],"-p") == 0 )
    {
      if ( ++k < argc )
      {
        vector<double> v = getDoublesFromCommaList(argv[k]);
        if ( v.size() != 4 )
        {
          cerr << "Error: expected 4 values in comma-separated list, no spaces, read " << v.size() << " numbers instead." << endl;
          exit(1);
        }
        double sum=0;
        for ( vector<double>::iterator p=v.begin(); p != v.end(); ++p )
        {
          if ( *p < 0 )
          {
            cerr << "Error: values after -p must be positive. Read " << *p << "." << endl;
            exit(1);
          }
          sum += *p;
        }
        for ( int i=0; i<4; ++i )
          stationaryP[i] = v[i]/sum;
      }
      else
      {
        cerr << "Error: flag `-p' not followed by comma-separated no-space list of four (relative) probabilities" << endl;
        usage(cerr);
      }
    }
    else if ( strcmp(argv[k],"-q") == 0 )
    {
      if ( ++k < argc )
      {
        vector<double> v = getDoublesFromCommaList(argv[k]);
        if ( v.size() != 6 )
        {
          cerr << "Error: expected 6 values in comma-separated list, no spaces, read " << v.size() << " numbers instead." << endl;
          exit(1);
        }
        double sum=0;
        for ( vector<double>::iterator p=v.begin(); p != v.end(); ++p )
        {
          if ( *p < 0 )
          {
            cerr << "Error: values after -q must be positive. Read " << *p << "." << endl;
            exit(1);
          }
          sum += *p;
        }
        for ( int i=0; i<6; ++i )
          symmetricQP[i] = v[i]/sum;
      }
      else
      {
        cerr << "Error: flag `-q' not followed by comma-separated no-space list of six (unweighted) values" << endl;
        usage(cerr);
      }
    }
    else // does not match a flag
    {
      cerr << "Option `" << argv[k] << "' does not match a valid program flag." << endl;
      usage(cerr);
    }
  }
}
