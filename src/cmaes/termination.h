#pragma once

#include <fmt/core.h>

#include <functional>
#include <initializer_list>
#include <range/v3/view/all.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/view/map.hpp>
#include <string>
#include <unordered_map>

#include "parameter.h"
#include "solution.h"

namespace ew_cmaes {

using termination_name = std::string;

using termination_callback = std::function<bool(const parameters&, const solutions&)>;

struct termination_status {
  bool terminate_{false};
  std::string msg_;
};

struct termination_criteria {
  auto check(const auto& params, const auto& sols) -> termination_status {
    for (auto&& [name, terminate_cb] : criteria_) {
      if (terminate_cb(params, sols)) {
        return termination_status{.terminate_ = true,
                                  .msg_ = fmt::format("Solver is terminated with the following "
                                                      "termination criteria: {}",
                                                      name)};
      }
    }
    return {};
  }

  auto add_critiera(const termination_name& name, const termination_callback& cb) -> void { criteria_.try_emplace(name, cb); }

  auto remove_critiera(const termination_name& name) {
    if (criteria_.contains(name)) {
      criteria_.erase(name);
    }
  }

  std::unordered_map<termination_name, termination_callback> criteria_;
};

namespace predefined {
const auto predefined_termination_criteria = termination_criteria{

    .criteria_ = {{"max_fevals", [](const auto& p, const auto& s) { return s.fevals_ > p.max_fevals_; }},
                  {"max_iter", [](const auto& p, const auto& s) { return s.fevals_ > p.max_fevals_; }},
                  {"stop_fitness", [](const auto& p, const auto& s) { return s.fevals_ > p.max_fevals_; }}}};
}

}  // namespace ew_cmaes
