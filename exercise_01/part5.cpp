#import <iostream>
#import <cstdio>
#import <fstream>
#import <algorithm>
#import <cmath>
#import <vector>
#import <string>

// For ordering the numbers
bool comparism(double a, double b) {
  return std::log10(std::abs(b)) < std::log10(std::abs(a));
}

// template type refering to variable type for summation
template <typename T>
void read_and_sum(std::string path, bool backwards, bool sorting) {

  // open file
  std::filebuf fb;
  fb.open(path, std::ios_base::in | std::ios_base::binary);

  // open file and read number of numbers
  int size = 0;
  fb.sgetn((char*)&size, sizeof(size));
  std::cout << "File " << path << " successfully read! There are " << size << " numbers." << std::endl;
  // read numbers into vector
  std::vector<double> numbers(size);
  fb.sgetn((char*)&numbers[0], numbers.size() * sizeof(numbers[0]));



  if(sorting) {
     // Sort with comparism-Function
     std::sort(numbers.begin(), numbers.end(), comparism);
     std::cout << "\tArray sorted." << std::endl;
  }

  T result = 0;
  // summation
  if(backwards) {
    std::cout << "\tSummation backwards." << std::endl;
    for(int i = size - 1; i >= 0; i--) {
      result += numbers[i];
    }
  }
  else {
    std::cout << "\tSummation forward." << std::endl;
    for(int i = 0; i < size; i++) {
      result += numbers[i];
    }
  }

  std::cout << "\tResult: " << std::scientific << result << std::endl;

}


int main() {

  read_and_sum<double>("numbers.dat", false, false); // forward without sorting

  read_and_sum<double>("numbers.dat", true, false); // backward without sorting

  read_and_sum<double>("numbers.dat", false, true); // forward with sorting

  read_and_sum<double>("numbers.dat", true, true); // backward with sorting

  read_and_sum<long double>("numbers.dat", false, true); // forward with sorting

  read_and_sum<long double>("numbers.dat", true, true); // backward with sorting

  return 0;

}
