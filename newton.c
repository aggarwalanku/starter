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

float* thomas(float **mat,int n)
{
  float *X=(float *)malloc(n*sizeof(float));
  int i;
  //Thomas algorithm
  for ( i = 1; i < n; i++)
    {
      float  m = mat[i][0]/mat[i-1][1];
      mat[i][1] = mat[i][1] - m * mat[i - 1][2];
      mat[i][3] = mat[i][3] - m*mat[i-1][3];
    }
  //Back substitution
  X[n-1] = mat[n-1][3]/mat[n-1][1];
  for ( i = n - 2; i >= 0; i--)
   X[i] = (mat[i][3] - mat[i][2] * X[i+1]) / mat[i][1];

  /*  printf("For step size h = %.4f \n",step);
  for(i=0;i<n;i++)
    {
      printf("y(%.2f) = %f \n",(i+1)*step,mat[i][4]);
    }
  printf("\n");*/
  return X;
}

void newton(float *coeff,int k)
{
  float A,B,step,b1,b2;
  A=coeff[0];B=coeff[1];b1=coeff[2];b2=coeff[3];
  step =  1/((float)k);
  int n = k-1,i;
  
  // initial solution
  float *y = (float *)malloc((k+1)*sizeof(float));
  float **mat = allocate(n,4);
  float *dy;

  y[0]=b1;y[n+1]=b2;
  for(i=1;i<n+1;i++)
    {
      y[i]=b1+i*(b2-b1)*step;
    }
  int j;
  for( j=0;j<5;j++)
    {
      //Creating the matrice for the equation
      
      for( i=1;i<n+1;i++)
	{
	  mat[i-1][0]=(A*y[i]+0.5*B*y[i-1]-0.5*B*y[i+1])/step*step;
	  mat[i-1][1]=A*(y[i+1]-4*y[i]+y[i-1])/step*step;
	  mat[i-1][2]=(A*y[i]+0.5*B*y[i+1]-0.5*B*y[i-1])/step*step;
	  mat[i-1][3]=-1*(A*y[i]*(y[i+1]-2*y[i]+y[i-1])+B*0.25*(y[i+1]-y[i-1])*(y[i+1]-y[i-1]))/step*step;
	}
      mat[0][0]=0;
      mat[n-1][2]=0;
      dy=thomas(mat,n);
      for(i=1;i<n+1;i++)
	{
	  y[i]=y[i]+dy[i-1];
	  printf("%.2f %f \n",(i)*step,y[i]);
	}
      printf("\n");
    }
}
      int main()
{
  float coeff[4];
  int i,n,p,q,r;

  printf("Enter the coefficients of y*y\"and (y')^2 \n");
  scanf("%f %f ",&coeff[0],&coeff[1]);
  printf("Enter  y(0),y(1) and three choices of number of steps \n");
  scanf("%f %f %d %d %d",&coeff[2],&coeff[3],&p,&q,&r);
  printf(" The equation is :- \n");
  printf("%.2fy*y\" + %.2f(y')^2 = 0 \n",coeff[0],coeff[1]);
  printf(" y(0) = %.2f  y(1) = %.2f \n",coeff[2],coeff[3]);
  newton(coeff,p);
  newton(coeff,q);
  newton(coeff,r);
  return 0;
}

