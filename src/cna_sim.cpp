/*
*   Copyright (C) 2020 Rishvanth Prabakar
*
*   Authors: Rish Prabakar
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <unistd.h>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

static string
print_usage(const string &name) {
  std::ostringstream oss;
  oss << name
      << " -i [in_file.txt]"
      << " -n [number of reads]"
      << " -o [out_file.txt]"
      << " -v verbose" << endl;
  return oss.str();
}

int
main (int argc, char *argv[]) {
  try {
    
    string in_file;
    string out_file;
    size_t n_reads{0};
    bool VERBOSE{false};    

    int opt;
    while ((opt = getopt(argc, argv, "i:n:o:v")) != -1) {
      if (opt == 'i')
        in_file = optarg;
      else if (opt == 'n')
        n_reads = atoi(optarg);
      else if (opt == 'o')
        out_file = optarg;
      else if (opt == 'v')
        VERBOSE = true;
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (in_file.empty() || out_file.empty() || !n_reads) 
      throw std::runtime_error(print_usage(argv[0]));

  }
  catch (std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  } 
  return EXIT_SUCCESS;
}
