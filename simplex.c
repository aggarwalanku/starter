#include <stdio.h>
#include <stdlib.h>
#include <math.h>
	
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
void rowoperation(float *dest,float*a,float *b,float coeff,int n)
{
  int i;
  for(i=0;i<n;i++)
    {
      dest[i]=a[i]+( b[i]*coeff);
    }
}
int optimality(float *a,int n)
{
  int i=0;
  for(i=0;i<n;i++)
    {
      if(fabs(a[i]) > 0.001 && a[i] < 0)
	return 0;
    }
  return 1;
}

void simplex(float ** table,int m,int n)
{
  int ent =1,lev =1,i,j,index;
  float min=0,minr=1000000;
  int *basic=(int *)malloc((m+1)*sizeof(int));
  basic[0]=0;
  for(i=0;i<m;i++)
    {
      basic[i+1]=n+1+i;
    }
  while(!optimality(table[0],m+n+1))
    {
      min =0;
      minr=1000000;
      for(i=0;i<m+n;i++)
	{
	  if(table[0][i+1] < min)
	    {
	      min = table[0][i+1];
	      ent = i+1;
	    }
	}
      printf("\n %d \n",ent);
      for(i=0;i<m;i++)
	{
	  if(fabs(table[i+1][ent]) > 0.001 && table[i+1][ent] > 0 )
	    {
	      float ratio = table[i+1][m+n+1]/table[i+1][ent];
	      if(ratio < minr)
		{
		  if( ratio > 0 )
		    {
		      minr = ratio;
		      lev = basic[i+1];
		      index=i+1;
		    }
		}
	    }
	}	  
      basic[index]=ent;
      float q = table[index][ent];
      for(i=0;i<m+n+2;i++)
	{
	  table[index][i]=table[index][i]/q;
	}
      for(i=0;i<m+1;i++)
	{
	  if(i != index)
	    {
	      rowoperation(table[i],table[i],table[index],-1*table[i][ent],m+n+2);
	    }
	}
      printf("\n");
    }
  printf("The optimal value of Z is %.3f at \n",table[0][m+n+1]);
  for(i=0;i<m;i++)
    {
      printf("x[%d] = %.3f\n",basic[i+1],table[i+1][m+n+1]);
    }
}	      	  	      
   	      
int main()
{
  int m,n,i,j;
  printf("Enter the number of equations and variables ");
  scanf("%d %d",&m,&n);
  float *z=(float *)malloc(n*sizeof(float));
  float **mat = allot2d(m,n+1);
  printf("Input the equation of Z\n");
  for(i=0;i<n;i++)
    {
      scanf("%f",&z[i]);
    }
  printf("Enter the original equations\n");
  for(i=0;i<m;i++)
    {
      printf("Enter the %d equation\n",i+1);
      for(j=0;j<n+1;j++)
	{
	  scanf("%f",(j+*(mat+i)));
	}
    }
  float **table=allot2d(m+1,m+n+2);
  for(i=0;i<n;i++)
    table[0][i+1]=-1*z[i];
  table[0][0]=1;
  //generating table
  for(i=1;i<=m;i++)
    {
      for(j=0;j<n;j++)
	table[i][j+1]=mat[i-1][j];
      table[i][n+i]=1;
      table[i][m+n+1]=mat[i-1][n];
    }
  printf("\n");
  simplex(table,m,n);
  return 0;
}

