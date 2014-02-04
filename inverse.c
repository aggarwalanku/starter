
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
float f(float *mat,int x,int y,int n)
{
  return mat[x*n+y];
}
*/

	
void printmat(float **mat,int m,int n)
{
  int i,j;
  for(i=0;i<m;i++)
    {
      for(j=0;j<n;j++)
	{
	  printf("%.2f ",*(*(mat+i)+j));
	}
      printf("\n");
    }
}
float** allot2d(int m,int n)
     {
  float **mat=(float **)malloc(m*sizeof(float *)); 
  int i;
  for(i=0;i<m;i++)
    {
      *(mat+i)=(float *)malloc(n*sizeof(float));
    }
  return mat;
}

float **copy(float **a,int n)
{
  int i,j;
  float **ret = allot2d(n,n);
  for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
	{
	  ret[i][j]= a[i][j];
	}
    }
  return ret;
}

/*
float* allotcate(int m,int n)
{
  float *mat=(float *)malloc(m*n*sizeof(float)); 
  return mat;
}*/

void identity(float **mat,int n)
{
  int i,j;
  for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
	{
	  if(i == j)
	    mat[i][j]=1;
	  else
	    mat[i][j]=0;
	}
    }
}

void rowoperation(float *dest,float*a,float *b,float coeff,int n)
{
  int i;
  for(i=0;i<n;i++)
    {
      dest[i]=a[i]+( b[i]*coeff);
    }
}
void inverse(float** mat,int n)
{
  float **inv=allot2d(n,n);
  identity(inv,n);
  int i,j;
  for(j=0;j<n;j++)
    {
      for(i=j+1;i<n;i++)
	{
	  float k = mat[i][j]/mat[j][j]*-1;
	  rowoperation(mat[i],mat[i],mat[j],k,n);
	  rowoperation(inv[i],inv[i],inv[j],k,n);
	}
    }
  for(j=n-1;j>=0;j--)
    {
      for(i=j-1;i>=0;i--)
	{
	  float k = mat[i][j]/mat[j][j]*-1;
	  rowoperation(mat[i],mat[i],mat[j],k,n);
	  rowoperation(inv[i],inv[i],inv[j],k,n);
	}
    }
  for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
	inv[i][j]= inv[i][j]/mat[i][i];
      mat[i][i]=1;
    }
  printf("\n");
  printmat(inv,n,n);
  printf("\n");
	
}
float ** create(float**mat,int k,int n)
{
  float** ret =allot2d(n,n);
  int i;
  identity(ret,n) ;
  for(i=0;i<n;i++)
    {
      ret[k][i]=mat[k][i];
    }
  return ret;
}
void mqq(float **mat,int n)
{
  int i=0;
  for(i=0;i<n;i++)
    {
      float ** co=copy(mat,n);
      float ** temp = create(co,i,n);
      printf("E%i is \n",i);
      inverse(temp,n);
    }
}
  
int main()
{
  int m,n,i,j;
  printf("Enter n ");
  scanf("%d",&n);
  float **mat = allot2d(n,n);
  printf("Enter the matrix\n");
  for(i=0;i<n;i++)
    {
      printf("Enter Row number %d\n",i+1);
      for(j=0;j<n;j++)
	{
	  scanf("%f",(j+*(mat+i)));
	}
    }
  printf("The matrix you have entered is \n");
  printmat(mat,n,n);
  printf("\n");
  // inverse(mat,n);
  // printmat(mat,m,n);
  mqq(mat,n);
  printf("Final M is \n");
  inverse(mat,n);
  return 0;
}



