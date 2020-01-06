//
// Created by mario on 31/07/19.
//

#ifndef CPPMODEL_STIRLING_NUMBERS_HPP
#define CPPMODEL_STIRLING_NUMBERS_HPP

#include "options.hpp"
#include <map>

typedef std::tuple<uint, uint> keyType;

class memoizer {
 private:
   std::map<keyType, unsigned long int> memo;

 public:
  template <class F>
  const unsigned long int& operator()(F f, int n, int m) {
    keyType key = std::make_tuple<uint, uint>(n, m);
    std::map<keyType, unsigned long int>::const_iterator it = memo.find(std::tie(n, m));
    if (it == memo.end()) {
      it = memo.insert(std::make_pair(key, f(n, m))).first;
    }
    return it->second;
  }
};

unsigned long int stirling(int n, int m);

namespace {
  unsigned long int stirling_(int n, int m) {
    if (((n == 0) & (m == 0)) ||( (n == 1) & (m == 1)))
      return 1;
    else if ((n > 0) & (m == 0))
      return 0;
    else if (m > n)
      return 0;
    else
      return stirling(n-1, m-1) + (n-1) * stirling(n-1, m);
  }
}


unsigned long int stirling(int n, int m) {
  static memoizer memo;
  return memo(stirling_, n, m);
}


#endif  // CPPMODEL_STIRLING_NUMBERS_HPP
