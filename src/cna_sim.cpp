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
#include <random>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <algorithm>

#include "Genome.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

struct ReadInfo {
  string read;
  string chr_name;
  size_t chr_pos;
};

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
    abs_pos_dist = std::uniform_int_distribution<size_t>(0, cn_genome_size);
  }

  CnaSim (const Genome &genome, const string& cna_regions) :
            CnaSim(genome, cna_regions, VERBOSE) {}
  CnaSim (const Genome &genome, const string& cna_regions, const bool V) {
    VERBOSE = V;
    set_cna(genome, cna_regions);
    abs_pos_dist = std::uniform_int_distribution<size_t>(0, cn_genome_size);
  }

  void sample_genome(const Genome &genome, const size_t readlen,
                      ReadInfo &info);

private:
  vector<size_t> chr;
  vector<size_t> start;
  vector<size_t> end; // half-open intervals
  vector<size_t> cn;
  vector<size_t> abs_pos;
  size_t cn_genome_size = 0;

  std::mt19937 rng{std::random_device()()};
  std::uniform_int_distribution<size_t> abs_pos_dist;
  bool VERBOSE = false;

  void set_default_cna (const Genome &genome);
  void set_cna (const Genome &genome, const string &cna_regions);
  void add_cna_segment (const size_t cna_chr, const size_t cna_start,
                        const size_t cna_end, const size_t cna_cn);
  string print_cna_list (const Genome &genome) const;
};

void
CnaSim::sample_genome (const Genome &genome, const size_t readlen,
                        ReadInfo &info) {

  size_t cn_abs_pos, index, normal_pos;
  string read;
  bool done = false;
  while (!done) {
    cn_abs_pos = abs_pos_dist(rng);
    index = std::upper_bound(abs_pos.begin(), abs_pos.end(), cn_abs_pos)
              - abs_pos.begin() - 1;
    normal_pos = ((cn_abs_pos - abs_pos[index])
                  % (end[index] - start[index])) + start[index];
    read = genome.chr_substr(chr[index], normal_pos, readlen);

    if (read.find('N') == string::npos && read.length() == readlen) {
      info.read = read;
      info.chr_name = genome.chr_tag(chr[index]);
      info.chr_pos = normal_pos;
      done = true;
    }
  }
}


void
CnaSim::add_cna_segment (const size_t cna_chr, const size_t cna_start,
                         const size_t cna_end, const size_t cna_cn) {

  chr.push_back(cna_chr);
  start.push_back(cna_start);
  end.push_back(cna_end);
  cn.push_back(cna_cn);
  abs_pos.push_back(cn_genome_size);

  cn_genome_size += (end.back() - start.back())*cn.back();
}


void
CnaSim::set_cna (const Genome &genome, const string &cna_regions) {
  std::ifstream in{cna_regions};
  if (!in)
    throw std::runtime_error("cannot open: " + cna_regions);

  const size_t normal_cn = 2;

  string line;
  vector<string> tokens;
  size_t curr_chr = 0;
  while (getline(in, line)) {
    split_string(line, tokens);
    const string cna_chr = tokens[0];
    const size_t cna_start = atoi(tokens[1].c_str());
    const size_t cna_end = atoi(tokens[2].c_str());
    const size_t cna_cn = atoi(tokens[3].c_str());

    if (chr.size() != 0 && genome.chr_tag(curr_chr) != cna_chr){
      if (end.back() != genome.chr_len(curr_chr))
        add_cna_segment(curr_chr, end.back(), genome.chr_len(curr_chr),
                        normal_cn);
      ++curr_chr;
    }

    while (genome.chr_tag(curr_chr) != cna_chr) {
      add_cna_segment(curr_chr, 0, genome.chr_len(curr_chr), normal_cn);
      ++curr_chr;
    }

    if (cna_start != 0) {
      size_t seg_start = 0;
      if (chr.size() != 0 && genome.chr_tag(chr.back()) == cna_chr)
        seg_start = end.back();
      if (cna_start - seg_start > 0)
        add_cna_segment(curr_chr, seg_start, cna_start, normal_cn);
    }

    add_cna_segment(curr_chr, cna_start, cna_end, cna_cn);
  }

  if (end.back() != genome.chr_len(curr_chr))
    add_cna_segment(curr_chr, end.back(), genome.chr_len(curr_chr), normal_cn);

  while (++curr_chr < genome.chr_count())
    add_cna_segment(curr_chr, 0, genome.chr_len(curr_chr), normal_cn);

  if (VERBOSE)
    cerr << print_cna_list(genome);

  in.close();
}

void
CnaSim::set_default_cna (const Genome &genome) {

  const size_t normal_cn = 2;
  const size_t n_chr = genome.chr_count();
  for (size_t i = 0; i < n_chr; ++i)
    add_cna_segment(i, 0, genome.chr_len(i), normal_cn);

  if (VERBOSE)
    cerr << print_cna_list(genome);
}

string
CnaSim::print_cna_list (const Genome &genome) const {
  std::ostringstream oss;
  oss << "[CNA REGIONS]" << endl;
  oss << "\tAltered genome size: " << cn_genome_size << endl;
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
format_fasta (const string &name, const ReadInfo &info) {
  std::ostringstream oss;
  oss << ">" << name
      << " " << info.chr_name
      << ":" << info.chr_pos << endl
      << info.read << endl;
  return oss.str();
}

static string
print_usage(const string &name) {
  std::ostringstream oss;
  oss << name << " [options]" << endl
      << "\t-i in_file.txt" << endl
      << "\t-c cna_regions.bed (required if TFx > 0.0)" << endl
      << "\t-t tumor fraction (default: 100)" << endl
      << "\t-n number of reads (default: 1M)" << endl
      << "\t-l read length (default: 100)" << endl
      << "\t-o out_file.txt" << endl
      << "\t-v verbose (default: false)" << endl;
  return oss.str();
}

int
main (int argc, char *argv[]) {
  try {

    string in_file;
    string out_file;
    string cna_regions;
    size_t n_reads{1000000};
    size_t read_len{100};
    float tfx{100};
    bool VERBOSE{false};

    int opt;
    while ((opt = getopt(argc, argv, "i:c:t:n:l:o:v")) != -1) {
      if (opt == 'i')
        in_file = optarg;
      else if (opt == 'c')
        cna_regions = optarg;
      else if (opt == 't')
        tfx = atof(optarg);
      else if (opt == 'n')
        n_reads = atoi(optarg);
      else if (opt == 'l')
        read_len = atoi(optarg);
      else if (opt == 'o')
        out_file = optarg;
      else if (opt == 'v')
        VERBOSE = true;
      else
        throw std::runtime_error(print_usage(argv[0]));
    }

    if (in_file.empty() || (tfx > 0 && cna_regions.empty())
        || tfx > 100 || out_file.empty())
      throw std::runtime_error(print_usage(argv[0]));

    if (VERBOSE)
      cerr << "[READING GENOME]" << endl;
    Genome genome(in_file, VERBOSE);

    std::ofstream out{out_file};
    if (!out)
      throw std::runtime_error("cannot open " + out_file);


    ReadInfo info;
    size_t read_count = 0;
    const size_t n_tumor_reads = (tfx/100) * n_reads;
    if (n_tumor_reads > 0) {
      if (VERBOSE)
        cerr << "[SETTING TUMOR CNA REGIONS]" << endl;
      CnaSim tumor_cna(genome, cna_regions, VERBOSE);

      if (VERBOSE) {
        cerr << "[GENERATING TUMOR READS]" << endl;
        cerr << "\tTumor reads: " << n_tumor_reads << endl;
      }
      while (++read_count <= n_tumor_reads) {
        tumor_cna.sample_genome(genome, read_len, info);
        out << format_fasta(string("tumor_" + std::to_string(read_count)),
                  info);
      }
    }

    const size_t n_normal_reads = n_reads - n_tumor_reads;
    if (n_normal_reads > 0) {
      if (VERBOSE)
        cerr << "[SETTING NORMAL CNA REGIONS]" << endl;
      CnaSim normal_cna(genome, VERBOSE);

      read_count = 0;
      if (VERBOSE){
        cerr << "[GENERATING NORMAL READS]" << endl;
        cerr << "\tNormal reads: " << n_normal_reads << endl;
      }
      while (++read_count <= n_normal_reads) {
        normal_cna.sample_genome(genome, read_len, info);
        out << format_fasta(string("normal_" + std::to_string(read_count)),
                  info);
      }
    }

    out.close();

  }
  catch (std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
