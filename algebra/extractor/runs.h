// MIT License
//
// Copyright (c) 2020 Jonas Ellert, jonas.ellert@tu-dortmund.de
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <vector>

namespace linear_time_runs {

template<typename index_type>
struct run_type {
  index_type start;
  index_type period;
  index_type length;
};

namespace internal {

template<bool increasing, typename value_type, typename index_type>
void unidirectional_runs(value_type const *const text, size_t const n,
                         std::vector <run_type<index_type>> &result) {

  auto compare = [&text](index_type const &i, index_type const &j,
                         index_type const &lce) {
      if constexpr (increasing) {
        return j > 0 && text[i + lce] > text[j + lce];
      } else {
        return text[i + lce] < text[j + lce];
      }
  };

  auto get_rlce = [&](index_type const &i, index_type const &j,
                      index_type lce = 0) {
      while (text[i + lce] == text[j + lce])
        ++lce;
      return lce;
  };

  auto get_llce = [&](index_type const &i, index_type const &j,
                      index_type llce = 0) {
      while (text[i - llce] == text[j - llce])
        ++llce;
      return llce;
  };

  struct distance_edge_type {
    index_type distance, lce;

    distance_edge_type(index_type const d, index_type const l)
        : distance(d), lce(l) {}
  };

  struct target_edge_type {
    index_type target, lce;

    target_edge_type(index_type const t, index_type const l)
        : target(t), lce(l) {}
  };

  std::deque <target_edge_type> q;
  auto top_j = [&q]() { return q.back().target; };
  auto top_lce = [&q]() { return q.back().lce; };
  auto pop = [&q]() { q.pop_back(); };
  auto push = [&q](index_type const j, index_type const lce) {
      q.emplace_back(j, lce);
  };

  push(0, 0);
  push(1, 0);

  std::vector <index_type> first_edge_of_node(n + 1);
  std::vector <distance_edge_type> edges;
  edges.reserve(2 * n);

  first_edge_of_node[1] = 0;
  edges.emplace_back(1, 0);

  index_type distance = 1;
  index_type rhs = 1;
  for (index_type i = 2; i < n - 1; ++i) {
    first_edge_of_node[i] = edges.size();

    index_type const copy_from = i - distance;
    index_type const stop_edge = first_edge_of_node[copy_from + 1];
    index_type e = first_edge_of_node[copy_from];

    for (; e < stop_edge; ++e) {
      if (i + edges[e].lce < rhs) {
        edges.emplace_back(edges[e]);
      } else
        break;
    }

    if (e == stop_edge) {
      index_type const target = i - edges.back().distance;

      while (top_j() > target)
        pop();
      push(i, edges.back().lce);
      continue;
    }

    index_type j = i - edges[e].distance;
    index_type lce = get_rlce(i, j, (rhs > i) ? (rhs - i) : (index_type) 0);
    rhs = i + lce;
    distance = i - j;
    edges.emplace_back(distance, lce);

    while (top_j() > j)
      pop();

    while (compare(i, j, lce)) {
      if (top_lce() < lce) {
        lce = top_lce();
        pop();
        edges.emplace_back(i - top_j(), lce);
        break;
      }

      pop();
      j = top_j();
      lce = get_rlce(i, j, lce);
      rhs = i + lce;
      distance = i - j;
      edges.emplace_back(distance, lce);
    }

    push(i, lce);
  }

  first_edge_of_node[n - 1] = edges.size();
  while (top_j() > 0) {
    edges.emplace_back(n - 1 - top_j(), 0);
    pop();
  }
  first_edge_of_node[n] = edges.size() + 1;

  struct nss_edge_type {
    index_type distance, llce, rlce;
  };

  std::vector <nss_edge_type> nss_edges(n);
  for (index_type i = 1; i < n; ++i) {
    index_type e = first_edge_of_node[i];
    index_type end_edge = first_edge_of_node[i + 1] - 1; // last edge = pss
    for (; e < end_edge; ++e) {
      nss_edges[i - edges[e].distance] =
          nss_edge_type{edges[e].distance, 0, edges[e].lce};
    }
  }

  nss_edges[n - 2].llce = 0;
  index_type lhs = n - 2;
  distance = 1;

  for (index_type i = n - 3; i > 0; --i) {
    if (i > lhs + nss_edges[i + distance].llce) {
      nss_edges[i].llce = nss_edges[i + distance].llce;
      continue;
    }

    nss_edges[i].llce = get_llce(i, i + nss_edges[i].distance,
                                 (lhs < i) ? (i - lhs) : (index_type) 0);
    lhs = i - nss_edges[i].llce;
    distance = nss_edges[i].distance;
  }

  std::vector <index_type> count_runs_at_idx(n);
  for (index_type i = 1; i < n - 1; ++i) {
    if (nss_edges[i].distance < nss_edges[i].llce + nss_edges[i].rlce) {
      ++count_runs_at_idx[i - nss_edges[i].llce + 1];
    }
  }

  index_type runs_with_duplicates = 0;
  for (index_type i = 1; i < n - 1; ++i) {
    index_type const gsize = count_runs_at_idx[i];
    count_runs_at_idx[i] = runs_with_duplicates;
    runs_with_duplicates += gsize;
  }

  if (runs_with_duplicates > 0) {
    auto prev_res_size = result.size();
    result.resize(prev_res_size + runs_with_duplicates);
    auto runs = &result[prev_res_size];
    for (index_type i = 1; i < n - 1; ++i) {
      if (nss_edges[i].distance < nss_edges[i].llce + nss_edges[i].rlce) {
        auto &border = count_runs_at_idx[i - nss_edges[i].llce + 1];
        runs[border++] = run_type<index_type>{
            i - nss_edges[i].llce + 1, nss_edges[i].distance,
            nss_edges[i].distance + nss_edges[i].llce + nss_edges[i].rlce - 1};
      }
    }

    index_type runs_without_duplicates = 1;
    for (index_type i = 1; i < runs_with_duplicates; ++i) {
      if ((runs[i].length != runs[i - 1].length) ||
          (runs[i].period != runs[i - 1].period) ||
          (runs[i].start != runs[i - 1].start)) {
        runs[runs_without_duplicates++] = runs[i];
      }
    }
    result.resize(prev_res_size + runs_without_duplicates);
  }
}

} // namespace internal

template<typename value_type, typename index_type = uint32_t>
auto compute_all_runs(value_type const *const text, size_t const n) {
  using run_t = run_type<index_type>;
  std::vector <run_t> result;
  internal::unidirectional_runs<true>(text, n, result);
  auto increasing = result.size();
  internal::unidirectional_runs<false>(text, n, result);

  std::inplace_merge(result.begin(), result.begin() + increasing, result.end(),
                     [](run_t const &a, run_t const &b) {
                         return a.start < b.start ||
                                (a.start == b.start && a.period < b.period);
                     });

  result.shrink_to_fit();
  return result;
}

} // namespace linear_time_runs
