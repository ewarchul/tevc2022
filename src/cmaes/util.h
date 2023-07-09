#pragma once

namespace ew_cmaes {

constexpr auto to_underlying(auto input)
    -> std::underlying_type_t<decltype(input)> {
  return static_cast<std::underlying_type_t<decltype(input)>>(input);
}

}  // namespace ew_cmaes
