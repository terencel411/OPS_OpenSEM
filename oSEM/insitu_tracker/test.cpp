#include "3D_slice_tracker.hh"
#include <array>
#include <cstddef>
#include <vector>

int main() {
  std::vector<double> data;

  for (int i = 0; i < 27; i++)
    data.push_back(i);

  std::array<std::size_t, 3> dims = {3, 3, 3};

  auto tracker = Insitu::SliceTracker3D<double, 0>(dims, 0, "foo", {0, 0}, 1);

  tracker.generate_graph(data.data());
}
