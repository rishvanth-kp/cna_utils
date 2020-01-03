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

#ifndef GENOME_HPP
#define GENOME_HPP

#include <vector>
#include <string>

using std::vector;
using std::string;

class Genome {
public:
  Genome (const string &in_file);
  Genome (const string &in_file, const bool V);

private:
  vector<string> chr_name;
  vector<string> chr_seq;
  vector<size_t> chr_abs_pos;
  size_t n_chr;
  size_t genome_size;
  bool VERBOSE = false;

  void read_genome (const string &in_file);
};


#endif
