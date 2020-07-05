//
// Created by nuccathy on 18-12-21.
//

#ifndef CLUSTERING_CHEMICAL_SIZEDISTRIBUTION_H
#define CLUSTERING_CHEMICAL_SIZEDISTRIBUTION_H

#include <vector>
#include <iostream>

class SizeDistribution {
 public:
  SizeDistribution(std::vector<int> vec) {
    this->distribution_ = vec;
  }

  SizeDistribution(int n) { this->distribution_.push_back(n); }

  SizeDistribution() = default;

  inline std::vector<int>::iterator Begin() {
    return this->distribution_.begin();
  }

  inline std::vector<int>::iterator End() {
    return this->distribution_.end();
  }

  inline std::vector<int>::const_iterator ConstBegin() const {
    return this->distribution_.cbegin();
  }

  inline std::vector<int>::const_iterator ConstEnd() const {
    return this->distribution_.cend();
  }

  inline int Get(G4int key) {
    return (this->distribution_)[key];
  }

  inline void Set(G4int key, G4int value) {
    this->distribution_[key] = value;
  }

  inline std::vector<int> &GetContent() {
    return this->distribution_;
  };

  inline void Push(int n) {
    this->distribution_.push_back(n);
  }

  SizeDistribution &operator+=(const SizeDistribution &other) {
    for (auto it = other.ConstBegin(); it != other.ConstEnd(); it++) {
      this->distribution_.push_back(*it);
    }

    return *this;
  }

//  SizeDistribution& operator+= ( SizeDistribution &other) {
//    for (auto it = other. Begin(); it != other. End(); it++) {
//      this->distribution_.push_back(*it);
//    }
//
//    return *this;
//  }

  SizeDistribution &operator+(const SizeDistribution &a) {
    SizeDistribution temp;

    for (auto it = a.ConstBegin(); it != a.ConstEnd(); it++) {
      temp.Push(*it);
    }

    for (auto it = this->ConstBegin(); it != this->ConstEnd(); it++) {
      temp.Push(*it);
    }

    return temp;
  }

  SizeDistribution &operator+=(int n) {
    this->distribution_.push_back(n);

    return *this;
  }

//  SizeDistribution& operator+= (const G4int &value) {
//    if (this->key_set_) {
//      this->distribution_[this->key_]=this->distribution_[this->key_]+value;
//      this->key_set_=false;
//    } else {
//      std::cerr << "Key not set in size distribution!" << std::endl;
//      exit(1);
//    }
//
//    return *this;
//  }

  SizeDistribution &operator/(double denominator) {
    for (auto entry : this->distribution_) {
      entry /= denominator;
    }

    return *this;
  }

 private:
//  G4int key_;
//  G4bool key_set_;

  std::vector<int> distribution_;
};

#endif //CLUSTERING_CHEMICAL_SIZEDISTRIBUTION_H
