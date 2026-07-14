#include "3D_slice_tracker.hh"
#include <array>
#include <cstddef>
#include <vector>

int main() {
  std::vector<double> data;

  std::array<std::size_t, 3> dims = {3, 4, 5};

  for (int i = 0; i < dims[0] * dims[1] * dims[2]; i++)
    data.push_back(i);

  auto tracker = Insitu::SliceTracker3D<double, 2>(dims, 0, "foo", {0, 0}, 1);

  tracker.generate_graph(data.data());
}
