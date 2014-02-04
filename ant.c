/* Siddharth Vishwakarma 
   10MA20041
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

float** allocate(int x,int y)                                   
{ int i;
  float **c=(float **)malloc(x*sizeof(float*));  
  for( i=0;i<x;i++)
    { c[i] = (float *)malloc(y*sizeof(float));
    }
  return c;
}
void thomas(float *coeff,int k)
{
  float A,B,C,D,E,step,b1,b2;
  A=coeff[0];  B=coeff[1];  C=coeff[2];  D=coeff[3];  E=coeff[4];  b1=coeff[5];  b2=coeff[6];
  step =  1/((float)k);
  int n = k-1,i;
  float **mat = allocate(n,5);

  // Creating the matrice for the equation

  // for the first and last row
  mat[0][0]=0;
  mat[0][1]=-1*(4*A-(C*2*step*step));
  mat[0][2]=(2*A+B*step);
  mat[0][3]=2*step*step*(D*step+E*sin(step))-(2*A+B*step)*b1;
  
  mat[n-1][0]=2*A-B*step;
  mat[n-1][1]=-1*(4*A-C*2*step*step);
  mat[n-1][2]=0;
  mat[n-1][3]=2*step*step*(D*n*step+E*sin(n*step))-(2*A+B*step)*b2;
  // other rows
  for( i=1;i<n-1;i++)
    {
      mat[i][0]=2*A-B*step;
      mat[i][1]=-1*(4*A-C*2*step*step);
      mat[i][2]=(2*A+B*step);
      mat[i][3]=2*step*step*(D*(i+1)*step+E*sin((i+1)*step));
    }
  //Thomas algorithm
  for ( i = 1; i < n; i++)
    {
      float  m = mat[i][0]/mat[i-1][1];
      mat[i][1] = mat[i][1] - m * mat[i - 1][2];
      mat[i][3] = mat[i][3] - m*mat[i-1][3];
    }
  //Back substitution
  mat[n-1][4] = mat[n-1][3]/mat[n-1][1];
  for ( i = n - 2; i >= 0; i--)
    mat[i][4] = (mat[i][3] - mat[i][2] * mat[i+1][4]) / mat[i][1];

  printf("For step size h = %.4f \n",step);
  for(i=0;i<n;i++)
    {
      printf("y(%.2f) = %f \n",(i+1)*step,mat[i][4]);
    }
  printf("\n");
}
int main()
{
  float coeff[7];
  int i,n,p,q,r;

  printf("Enter the coefficients of y\", y',y,x,sin(x) \n");
  scanf("%f %f %f %f %f",&coeff[0],&coeff[1],&coeff[2],&coeff[3],&coeff[4]);
  printf("Enter  y(0),y(1) and three choices of number of steps \n");
  scanf("%f %f %d %d %d",&coeff[5],&coeff[6],&p,&q,&r);
  printf(" The equation is :- \n");
  printf("%.2fy\" + %.2fy'+ %.2fy = %.2fx + %.2fsin(x)\n",coeff[0],coeff[1],coeff[2],coeff[3],coeff[4]);
  printf(" y(0) = %.2f  y(1) = %.2f \n",coeff[5],coeff[6]);
  
  thomas(coeff,p);
  thomas(coeff,q);
  thomas(coeff,r);
    
  return 0;
}

