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

#include<algorithm>
#include<iostream>
#include<fstream>

#include "Genome.hpp"

using std::cout;
using std::cerr;
using std::endl;

Genome::Genome (const string &in_file, const bool V) {
  VERBOSE = V;
  read_genome(in_file);
}

void
Genome::read_genome (const string &in_file) {
  std::ifstream in{in_file};
  if (!in)
    throw std::runtime_error("Cannot open " + in_file);

  string line;
  while (getline(in, line)){
    if (line[0] == '>') {
      chr_name.push_back(line.substr(1));
      chr_seq.push_back(string{});
      ++n_chr;
    }
    else
      chr_seq.back() += line;
  }

  chr_abs_pos.push_back(0);
  for (size_t i = 0; i < n_chr; ++i) {
    std::transform(chr_seq[i].begin(), chr_seq[i].end(),
                    chr_seq[i].begin(), toupper);
    genome_size += chr_seq[i].length();
    chr_abs_pos.push_back(genome_size);
  }

  if (VERBOSE) {
    cerr << "[GENOME STATS]" << endl;
    cerr << "\tGenome size: " << genome_size << endl;
    cerr << "\tchr\tsize\tabs_pos" << endl;
    for (size_t i = 0; i < n_chr; ++i) {
      cerr << "\t" << chr_name[i]
           << "\t" << chr_seq[i].length()
           << "\t" << chr_abs_pos[i] << endl;
    }
  }

  in.close();
}


string
Genome::chr_substr (const size_t chr, const size_t pos,
                     const size_t len) const {
  return chr_seq[chr].substr(pos, len);
}
