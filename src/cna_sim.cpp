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
#include <fstream>
#include <unistd.h>

#include "Genome.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

void
split_string (const string &in, vector<string> &tokens,
              const char delim = '\t'){

  tokens.clear();
  size_t pos = 0;
  size_t next_pos = in.find(delim);
  while (next_pos != string::npos) {
    tokens.push_back(in.substr(pos, next_pos - pos));
    pos = ++next_pos;
    next_pos = in.find(delim, pos);
  }
  tokens.push_back(in.substr(pos));
}

class CnaSim {
public:
  CnaSim (const Genome &genome) : CnaSim(genome, VERBOSE) {}
  CnaSim (const Genome &genome, const bool V) {
    VERBOSE = V;
    set_default_cna(genome);
  }
  CnaSim (const Genome &genome, const string& cna_regions) :
            CnaSim(genome, cna_regions, VERBOSE) {}
  CnaSim (const Genome &genome, const string& cna_regions, const bool V) {
    VERBOSE = V;
    set_cna(genome, cna_regions);
  }

private:
  vector<size_t> chr;
  vector<size_t> start;
  vector<size_t> end; // half-open intervals
  vector<size_t> cn;
  vector<size_t> abs_pos;
  bool VERBOSE = false;

  void set_default_cna (const Genome &genome);
  void set_cna (const Genome &genome, const string &cna_regions);
  string print_cna_list (const Genome &genome);
};

void
CnaSim::set_cna (const Genome &genome, const string &cna_regions) {
  std::ifstream in{cna_regions};
  if (!in)
    throw std::runtime_error("cannot open: " + cna_regions);

  const size_t normal_cn = 2;

  string line;
  vector<string> tokens;
  string cna_chr;
  size_t cna_start, cna_end, cna_cn;
  size_t curr_chr = 0;
  while (getline(in, line)) {
    split_string(line, tokens);
    cna_chr = tokens[0];
    cna_start = atoi(tokens[1].c_str());
    cna_end = atoi(tokens[2].c_str());
    cna_cn = atoi(tokens[3].c_str());

    if (chr.size() != 0 && genome.chr_tag(curr_chr) != cna_chr){
      if (end.back() != genome.chr_len(curr_chr)) {
        chr.push_back(curr_chr);
        start.push_back(end.back());
        end.push_back(genome.chr_len(curr_chr));
        cn.push_back(normal_cn);
      }
      ++curr_chr;
    }

    while (genome.chr_tag(curr_chr) != cna_chr) {
      chr.push_back(curr_chr);
      start.push_back(0);
      end.push_back(genome.chr_len(curr_chr));
      cn.push_back(normal_cn);
      ++curr_chr;
    }

    if (cna_start != 0) {
      if (genome.chr_tag(chr.back()) != cna_chr)
        start.push_back(0);
      else
        start.push_back(end.back());
      chr.push_back(curr_chr);
      end.push_back(cna_start);
      cn.push_back(normal_cn);
    }

    chr.push_back(curr_chr);
    start.push_back(cna_start);
    end.push_back(cna_end);
    cn.push_back(cna_cn);
  }

  if (end.back() != genome.chr_len(curr_chr)) {
    chr.push_back(curr_chr);
    start.push_back(end.back());
    end.push_back(genome.chr_len(curr_chr));
    cn.push_back(normal_cn);
  }

  while (++curr_chr < genome.chr_count()) {
    chr.push_back(curr_chr);
    start.push_back(0);
    end.push_back(genome.chr_len(curr_chr));
    cn.push_back(normal_cn);
  }

  abs_pos.push_back(0);
  for (size_t i = 1; i < chr.size(); ++i)
    abs_pos.push_back(abs_pos[i-1] + ((end[i-1] - start[i-1])*cn[i-1]));

  if (VERBOSE)
    cerr << print_cna_list(genome) << endl;

  in.close();
}

void
CnaSim::set_default_cna (const Genome &genome) {

  const size_t normal_cn = 2;
  const size_t n_chr = genome.chr_count();
  for (size_t i = 0; i < n_chr; ++i) {
    chr.push_back(i);
    start.push_back(0);
    end.push_back(genome.chr_len(i));
    cn.push_back(normal_cn);
  }

  abs_pos.push_back(0);
  for (size_t i = 1; i < chr.size(); ++i)
    abs_pos.push_back(abs_pos[i-1] + ((end[i-1] - start[i-1])*cn[i-1]));

  if (VERBOSE)
    cerr << print_cna_list(genome) << endl;
}

string
CnaSim::print_cna_list (const Genome &genome) {
  std::ostringstream oss;
  oss << "[CNA REGIONS]" << endl;
  oss << "\tchr\tstart\tend\tcn\tabs_pos" << endl;
  for (size_t i = 0; i < chr.size(); ++i) {
    oss << "\t" << genome.chr_tag(chr[i])
        << "\t" << start[i]
        << "\t" << end[i]
        << "\t" << cn[i]
        << "\t" << abs_pos[i] << endl;
  }
  return oss.str();
}

static string
print_usage(const string &name) {
  std::ostringstream oss;
  oss << name
      << " -i [in_file.txt]"
      << " -c [cna_regions.bed]"
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
    string cna_regions;
    size_t n_reads{0};
    bool VERBOSE{false};

    int opt;
    while ((opt = getopt(argc, argv, "i:c:n:o:v")) != -1) {
      if (opt == 'i')
        in_file = optarg;
      else if (opt == 'c')
        cna_regions = optarg;
      else if (opt == 'n')
        n_reads = atoi(optarg);
      else if (opt == 'o')
        out_file = optarg;
      else if (opt == 'v')
        VERBOSE = true;
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (in_file.empty() || cna_regions.empty() || out_file.empty()
        || !n_reads)
      throw std::runtime_error(print_usage(argv[0]));

    if (VERBOSE)
      cerr << "[READING GENOME]" << endl;
    Genome genome(in_file, VERBOSE);

    if (VERBOSE)
      cerr << "[SETTING TUMOR CNA REGIONS]" << endl;
    CnaSim tumor_cna(genome, cna_regions, VERBOSE);

    if (VERBOSE)
      cerr << "[SETTING NORMAL CNA REGIONS]" << endl;
    CnaSim normal_cna(genome, VERBOSE);

  }
  catch (std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
