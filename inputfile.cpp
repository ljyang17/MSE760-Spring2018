#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <math.h>

int main(int argc, char const *argv[]) {
  const std::string demo_dir = "/Users/lijingyang/Documents/MATLAB/chrono-calc";
  std::ifstream input;
    input.open((demo_dir + "/analytic-pflow").c_str());
    char str[100];
    double analytical[100];
    for (int i = 0; i < 21; i++) {
      input.getline(str,100);
      analytical[i] = atof(str);
      printf("%f\n", analytical[i]);
    }
  return 0;
}
