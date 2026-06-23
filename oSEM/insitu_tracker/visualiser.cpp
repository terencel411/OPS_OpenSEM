#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

void println(std::string string) { std::cout << string << std::endl; }

void load_new_image(int num) {
  auto file = std::ifstream("test_output" + std::to_string(num));

  if (!file.is_open()) {
    return;
  }

  system("clear");

  auto curr_line = std::string{};
  while (std::getline(file, curr_line)) {
    std::cout << curr_line << std::endl;
  }

  std::cout << "File Name: " << "test_output" + std::to_string(num)
            << std::endl;

  file.close();
}

int main() {
  bool moving_forward = true;

  println("|===| Visualiser |===|");
  println("- To navigate through slides, use enter with no input");
  println("'f', Enter -> Enter moves Forward (default)");
  println("'b', Enter -> Enter moves Backward (default)");
  println("'start', Enter -> Show first slide");
  println("exit, Enter -> finish program");

  std::string input;

  int current_image = 0;
  int min = 0;
  int max = 1999;

  while (true) {
    std::getline(std::cin, input);

    if (input == "exit")
      return EXIT_SUCCESS;

    if (input == "f") {
      moving_forward = true;
    } else if (input == "b") {
      moving_forward = false;
    } else if (input == "start") {
      current_image = 0;
      load_new_image(current_image);
    } else if (input == "") {
      current_image = moving_forward ? current_image + 1 : current_image - 1;

      if (current_image < min)
        current_image = max;

      if (current_image > max)
        current_image = min;

      load_new_image(current_image);
    }

    input = "";
  }
}
