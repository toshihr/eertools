#ifndef EERTOOLS_H
#define EERTOOLS_H

#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <string>
#include <fstream>
#include <map>
#include <random>
#include <memory>
#include <functional>

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T> >
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };
}

namespace eertools
{
typedef std::vector<double> vec_double;
typedef std::vector< vec_double > vmat_double;
typedef std::vector<std::string> vec_string;
typedef std::vector< vec_string > vmat_string;
typedef std::vector<int> vec_int;
typedef std::vector< vec_int > vmat_int;
typedef std::list<double> list_double;
typedef std::vector<long double> vec_ldouble;
typedef std::shared_ptr<std::string> pString;
typedef std::vector<pString> vec_pString;

const int PROB_TABLE_SIZE = 0x7f + 1;
const double ln2 = log(2);

bool calcComplexProb(vmat_double& complexProb, const std::string& s1, const std::string& s2);
bool calcProb(vec_double& prob, const std::string s);
double getEntropy(const vec_double& prob);
double getMutualEntropy(const vec_double& prob1, const vec_double& prob2, const vmat_double& complexProb);
double getEER(const vec_double& prob1, const vec_double& prob2, const vmat_double& complexProb);
double getEER2(const vec_double& prob1, const vec_double& prob2, const vmat_double& complexProb);
double getECD(const vec_double& prob1, const vec_double& prob2, const vmat_double& complexProb);
double getECD(const vec_double& nodes, const int& nStart, const int& L, const int& M);
double getECD(const std::string& seq);
double getLYAP(const vec_double& nodes, const int& nStart, const int& L, const std::function<double (double t)> & df);
}

#endif
