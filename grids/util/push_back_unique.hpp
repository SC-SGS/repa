#pragma once

namespace repa {
namespace util {

template <typename T>
void push_back_unique(std::vector<T>& v, T val)
{
  if (std::find(std::begin(v), std::end(v), val) == std::end(v))
    v.push_back(val);
}

}
}

