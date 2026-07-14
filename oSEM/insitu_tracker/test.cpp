#include "3D_slice_tracker.hh"
#include <array>
#include <cstddef>
#include <vector>

int main() {
  std::vector<double> data;

  std::array<std::size_t, 3> dims = {100, 400, 700};

  for (int i = 0; i < dims[0] * dims[1] * dims[2]; i++)
    data.push_back(i);

  auto tracker_x =
      Insitu::SliceTracker3D<double, 0>(dims, 0, "x.bmp", {0, 0}, 1);
  auto tracker_y =
      Insitu::SliceTracker3D<double, 1>(dims, 0, "y.bmp", {0, 0}, 1);
  auto tracker_z =
      Insitu::SliceTracker3D<double, 2>(dims, 0, "z.bmp", {0, 0}, 1);

  tracker_x.generate_graph(data.data());
  tracker_y.generate_graph(data.data());
  tracker_z.generate_graph(data.data());
}
