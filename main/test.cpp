#include "options.h"
#include <string>
#include <complex.h>

#include "meshreader.h"

using namespace MainNamespace;

int main(int argc, char **argv)
{
  Options options(argc,argv);
  // MeshReader reader;
  auto tmp = parseGmsh(meshFile());
  
}
