#include <stdio.h>
#include <math.h>
#define METHODS 4
void inverse(double a[][100],int size);
void gauss(double a[][100],int size);
void diagonal(double a[][100],int size);
void jacobi(double a[][100],int size,double x0[][100],double epsilon);
void seidel(double a[][100],int size,double x0[][100],double epsilon);
int main(int argc, char const *argv[]) {
  char operator;
  int size;
  static double a[100][100],x0[2][100],epsilon;
  printf("This program includes several equation solving methods involving matrices.\n"
          "1- The inverse of a matrix             - https://en.wikipedia.org/wiki/Gaussian_elimination#Finding_the_inverse_of_a_matrix\n"
          "2- Gauss Elimination                   - https://en.wikipedia.org/wiki/Gaussian_elimination\n"
          "3- Jacobi method                       - https://en.wikipedia.org/wiki/Jacobi_method\n"
          "4- Gauss Seidel                        - https://en.wikipedia.org/wiki/Gauss–Seidel_method\n"
          "Choose one of the methods above (1-%d): ",METHODS);scanf("%c", &operator);
  switch (operator) {
    case '1':
    printf("Enter the size of the matrix nXn (n): ");scanf("%d",&size);
    for (size_t i = 0; i < size; i++) {
      a[i][size+i]=1;
    }
    for (size_t i = 0; i < size; i++) {
      for (size_t j = 0; j < size; j++) {
        printf("a[%ld][%ld]: ",i+1,j+1);scanf("%lf",&a[i][j]);
      }
    }
    printf("\n");
    inverse(a,size);
    for (size_t i = 0; i < size; i++) {
      for (size_t j = size; j < 2*size; j++) {
        printf("%4.4lf\t",a[i][j]);
      }
      printf("\n");
    }
    break;
    case '2':
    printf("Given the equations are represented as a[1]*X1+a[2]*X2+...+a[n]*Xn=C\nEnter the number of equations / variables: ");scanf("%d",&size);
    for (size_t i = 0; i < size; i++) {
      for (size_t j = 0; j < size; j++) {
        printf("In the %ld.row: Enter the coefficient of x%ld: ",i+1,j+1);scanf("%lf",&a[i][j]);
      }
    printf("Enter the %ld. constant (c%ld): ",i+1,i+1);scanf("%lf",&a[i][size]);
    }
    printf("\n");
    gauss(a,size);
    for (size_t i = 0; i < size; i++) {
      printf("x%ld=%lf\n",i+1,a[i][size]);
    }
    break;
    case '3':
    printf("Before starting, enter an epsilon value (i.e. max error accepted): ");scanf("%lf",&epsilon);
    printf("Enter the number of equations / variables: ");scanf("%d",&size);
    for (size_t i = 0; i < size; i++) {
      for (size_t j = 0; j < size; j++) {
        printf("In the %ld.row: Enter the coefficient of x%ld: ",i+1,j+1);scanf("%lf",&a[i][j]);
      }
    printf("Enter the %ld. constant (c%ld): ",i+1,i+1);scanf("%lf",&a[i][size]);
    }
    printf("\n");
    printf("Enter the starting value for:\n");
    for (size_t i = 0; i < size; i++) {
      printf("x%ld: ",i+1);scanf("%lf",&x0[0][i]);
    }
    diagonal(a,size);
    jacobi(a,size,x0,epsilon);
    for (size_t i = 0; i < size; i++) {
      printf("x%ld=%lf\n",i+1,x0[1][i]);
    }
    break;
    case '4':
    printf("Before starting, enter an epsilon value (i.e. max error accepted): ");scanf("%lf",&epsilon);
    printf("Enter the number of equations / variables: ");scanf("%d",&size);
    for (size_t i = 0; i < size; i++) {
      for (size_t j = 0; j < size; j++) {
        printf("In the %ld.row: Enter the coefficient of x%ld: ",i+1,j+1);scanf("%lf",&a[i][j]);
      }
    printf("Enter the %ld. constant (c%ld): ",i+1,i+1);scanf("%lf",&a[i][size]);
    }
    printf("\n");
    printf("Enter the starting value for:\n");
    for (size_t i = 0; i < size; i++) {
      printf("x%ld: ",i+1);scanf("%lf",&x0[0][i]);
    }
    diagonal(a,size);
    seidel(a,size,x0,epsilon);
    for (size_t i = 0; i < size; i++) {
      printf("x%ld=%lf\n",i+1,x0[1][i]);
    }
    break;
  }
  return 0;
}
void inverse(double a[][100],int size) {
  double div;
  for (size_t i = 0; i < size; i++) {
    div=a[i][i];
    for (size_t j = i; j < size*2; j++) {
      a[i][j]/=div;
    }
    for (size_t j = i+1; j < size; j++) {
      div=a[j][i];
      for (size_t k = 0; k < size*2; k++) {
        a[j][k]-=div*a[i][k];
      }
    }
  }
  for (size_t i = 1; i < size; i++) {
    for (size_t j = i; j < size; j++) {
      div=a[size-j-1][size-i];
      a[size-j-1][size-i]-=div*a[size-i][size-i];
      for (size_t k = 0; k < size; k++) {
        a[size-j-1][size+k]-=div*a[size-i][size+k];
      }
    }
  }
}
void gauss(double a[][100],int size) {
  double div;
  for (size_t i = 0; i < size; i++) {
    div=a[i][i];
    for (size_t j = i; j < size+1; j++) {
      a[i][j]/=div;
    }
    for (size_t j = i+1; j < size; j++) {
      div=a[j][i];
      for (size_t k = 0; k < size+1; k++) {
        a[j][k]-=div*a[i][k];
      }
    }
  }
  for (size_t i = 1; i < size; i++) {
    for (size_t j = i; j < size; j++) {
      div=a[size-j-1][size-i];
      a[size-j-1][size-i]-=div*a[size-i][size-i];
      for (size_t k = 0; k < size; k++) {
        a[size-j-1][size+k]-=div*a[size-i][size+k];
      }
    }
  }
}
void diagonal(double a[][100],int size) {
  int t,k,counter,elements=0;
  static double max[3][1000],rows[100],copy[100][100];
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      counter=0;
      while (fabs(a[i][j])>max[0][elements-counter] && counter<=elements) {
        max[0][elements-counter+1]=max[0][elements-counter];
        max[1][elements-counter+1]=max[1][elements-counter];
        max[2][elements-counter+1]=max[2][elements-counter];
        counter++;
      }
        max[0][elements-counter+1]=fabs(a[i][j]);
        max[1][elements-counter+1]=i;
        max[2][elements-counter+1]=j;
      elements++;
    }
  }
  for (size_t i = 0; i < size; i++) {
    rows[i]=-1;
  }
  counter=0;k=0;
  while (counter<size) {
    if (rows[(int)max[1][k]]==-1) {
      t=0;
      while (max[2][k]!=rows[t] && t<size) {
        t++;
      }
      if (t==size) {
        rows[(int)max[1][k]]=max[2][k];
        counter++;
      }
    }
    k++;
  }
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size+1; j++) {
      copy[i][j]=a[i][j];
    }
  }
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size+1; j++) {
      a[(int)rows[i]][j]=copy[i][j];
    }
  }
}
void jacobi(double a[][100],int size,double x0[][100],double epsilon) {
  int positive=0;
  while (positive<size) {
    for (size_t i = 0; i < size; i++) {
      x0[1][i]=a[i][size]+x0[0][i]*a[i][i];
      for (size_t j = 0; j < size; j++) {
        x0[1][i]-=a[i][j]*x0[0][j];
      }
      x0[1][i]/=a[i][i];
    }
    positive=0;
    for (size_t i = 0; i < size; i++) {
      if (fabs(x0[1][i]-x0[0][i])<epsilon) {
        positive++;
      }
    }
    for (size_t i = 0; i < size; i++) {
      x0[0][i]=x0[1][i];
    }
  }
}
void seidel(double a[][100],int size,double x0[][100],double epsilon) {
  int positive=0;
  static double swap;
  while (positive<size) {
    for (size_t i = 0; i < size; i++) {
      x0[1][i]=a[i][size]+x0[0][i]*a[i][i];
      for (size_t j = 0; j < size; j++) {
        x0[1][i]-=a[i][j]*x0[0][j];
      }
      x0[1][i]/=a[i][i];
      swap=x0[0][i];
      x0[0][i]=x0[1][i];
      x0[1][i]=swap;
    }
    positive=0;
    for (size_t i = 0; i < size; i++) {
      if (fabs(x0[1][i]-x0[0][i])<epsilon) {
        positive++;
      }
    }
    for (size_t i = 0; i < size; i++) {
      x0[1][i]=x0[0][i];
    }
  }
}
