#pragma once

#include "common.h"

#include <string>
#include <iostream>
#include <fstream>
#include <array>

void print(const DArray2 &data, int i_start, int i_end, int j_start, int j_end, std::ofstream &out);
void output(const std::string &header, const DArray2 &data, int i_start, int i_end, int j_start, int j_end, std::ofstream &out);
void output(const std::string &header, const DArray2 &data, const std::array<int, 4> &range, std::ofstream &out);
void output(const std::string &header, const DArray2 &data, std::ofstream &out);
