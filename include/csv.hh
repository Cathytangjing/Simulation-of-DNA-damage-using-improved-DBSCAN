//
// Created by larry on 2020/4/6.
//

#ifndef CLUSTERING_CHEMICAL_INCLUDE_CSV_H_
#define CLUSTERING_CHEMICAL_INCLUDE_CSV_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

class CSVWriter {
 public:
  CSVWriter(std::string file_name, std::vector<std::string> header) : col_index_(0) {
    s_.open(file_name, std::ios::out | std::ios::trunc);

    for (int i = 0; i < header.size(); i++) {
      if (i != 0) s_ << ", ";
      s_ << header[i];
    }
    s_ << std::endl;

    s_.close();

    file_name_ = file_name;
    header_ = header;
  }

  template <typename T>
  void WriteRow(std::vector<T> content) {
    if (col_index_ != 0) {
      std::cerr << "Unfinished row!" << std::endl;
      exit(233);
    }

    s_.open(file_name_, std::ios::app);
    for (int i = 0; i < content.size(); i++) {
      if (i != 0) s_ << ", ";
      s_ << content[i];
    }
    s_ << std::endl;

    s_.close();
  }

  template <typename T>
  void WriteSingle(T value) {
    s_.open(file_name_, std::ios::app);

    if (col_index_ != 0) s_ << ", ";
    s_ << value;

    col_index_++;
    if (col_index_ == header_.size()) {
      s_ << std::endl;
      col_index_ = 0;
    }

    s_.close();
  }

 private:
  std::string file_name_;
  std::vector<std::string> header_;
  std::ofstream s_;

  size_t col_index_;
};

#endif //CLUSTERING_CHEMICAL_INCLUDE_CSV_H_
