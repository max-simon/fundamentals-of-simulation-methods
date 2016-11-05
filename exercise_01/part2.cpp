#include <iostream>

int main() {

  ///////////////////
  // Part 2
  //////////////////
  std::cout << "--- Part 2 --- \n" << std::endl;

  double a = 1.0e17;
  double b = -1.0e17;
  double c = 1.0;
  double x = (a + b) + c;
  double y = a + (b + c);

  std::cout << "x = " << x << std::endl;
  std::cout << "y = " << y << std::endl;

  return 0;
}
