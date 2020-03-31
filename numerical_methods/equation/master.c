#include <stdio.h>
#include "tinyexpr.h"
int main(int argc, char const *argv[]) {
  char operator;
  printf("######  ####  #    #   ##   ##### #  ####  #    #     ####   ####  #      #    # ###### #####  \n");
  printf("#      #    # #    #  #  #    #   # #    # ##   #    #      #    # #      #    # #      #    # \n");
  printf("#####  #    # #    # #    #   #   # #    # # #  #     ####  #    # #      #    # #####  #    # \n");
  printf("#      #  # # #    # ######   #   # #    # #  # #         # #    # #      #    # #      #####  \n");
  printf("#      #   #  #    # #    #   #   # #    # #   ##    #    # #    # #       #  #  #      #   #  \n");
  printf("######  ### #  ####  #    #   #   #  ####  #    #     ####   ####  ######   ##   ###### #    # \n");
  printf("This program finds roots using diffrenet numerical methods:\n"
          "1- Graphical Method             - https://wikipedia.org/wiki/List_of_graphical_methods\n"
          "2- Bisection method             - https://wikipedia.org/wiki/Bisection_method\n"
          "3- Regula-falsi method          - https://wikipedia.org/wiki/Regula_falsi\n"
          "4- Fixed-point iteration method - https://wikipedia.org/wiki/Fixed-point_iteration\n"
          "5- Newton-raphson method        - https://wikipedia.org/wiki/Newton's_method\n"
          "6- Secant method                - https://wikipedia.org/wiki/Secant_method\n"
          "Choose one of the methods above (1-6):");scanf("%c", &operator);
  switch (operator) {
    case '1':
    printf("case 1\n");
    break;
  }
  return 0;
}
