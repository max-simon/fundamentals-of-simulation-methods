#include <iostream>
#include <limits>
#include <iomanip>
#include <string>

using namespace std;

template <typename T>
int print_epsilon(string name) {
  // generate epsilon of Type T
  T epsilon = 1;
  int precession = 0;
  // while there is a difference
  while(1 + epsilon != 1) {
    // divide epsilon per 2
    epsilon /= 2;
    precession++;
  }

  T one = 1;

  // Output
  cout << "Machine-Epsilon of " << name << ": " << 2*epsilon << endl;
  // 2 times epsilon, because per definition the last, where there is a difference
  cout << "Precession (bits of mantissa): " << precession - 1 << endl;
  // - 1, because per definiton the last, where there is a difference

  cout << "1 + epsilon = " << setprecision(30) << one + (2*epsilon) << endl;

  return epsilon;

}

int main() {

  ///////////////////
  // Part 1
  //////////////////
  cout << "--- Part 1 --- \n" << endl;

  // machine precession of float
  print_epsilon<float>("float");

  // machine precession of double
  print_epsilon<double>("double");

  // machine precession of double
  print_epsilon<long double>("long double");



  return 0;

}
