#include <SDL2/SDL.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

void println(std::string string) { std::cout << string << std::endl; }

auto new_image_surface(SDL_Renderer *renderer, int image_num) -> SDL_Texture * {
  auto file_name = "bmp_output" + std::to_string(image_num) + ".bmp";
  SDL_Surface *surface = SDL_LoadBMP(file_name.c_str());
  assert(surface);

  SDL_Texture *texture = SDL_CreateTextureFromSurface(renderer, surface);

  SDL_FreeSurface(surface);

  return texture;
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

  if (SDL_Init(SDL_INIT_VIDEO) != 0) {
    std::cerr << SDL_GetError() << "\n";
    return EXIT_FAILURE;
  }

  SDL_Window *window =
      SDL_CreateWindow("BMP Viewer", SDL_WINDOWPOS_CENTERED,
                       SDL_WINDOWPOS_CENTERED, 1000, 1000, SDL_WINDOW_SHOWN);
  SDL_Renderer *renderer =
      SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

  if (!renderer)
    std::cerr << SDL_GetError() << std::endl;

  SDL_Event event;

  while (true) {
    std::getline(std::cin, input);

    SDL_Texture *texture;
    if (input == "exit")
      return EXIT_SUCCESS;

    if (input == "f") {
      moving_forward = true;
    } else if (input == "b") {
      moving_forward = false;
    } else if (input == "start") {
      current_image = 0;
      texture = new_image_surface(renderer, current_image);
    } else if (input == "") {
      current_image = moving_forward ? current_image + 1 : current_image - 1;

      if (current_image < min)
        current_image = max;

      if (current_image > max)
        current_image = min;

      texture = new_image_surface(renderer, current_image);
    }

    SDL_RenderClear(renderer);
    SDL_RenderCopy(renderer, texture, nullptr, nullptr);
    SDL_RenderPresent(renderer);

    SDL_DestroyTexture(texture);

    while (SDL_PollEvent(&event)) {
      if (event.type == SDL_QUIT)
        break;
    }

    input = "";
  }

  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();
}
