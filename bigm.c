#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M 100000
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
int optimality(float *a,int n,float c)
{
  int i=0;
  for(i=0;i<n;i++)
    {
      if(fabs(a[i]) > 0.001 && c*a[i] < 0)
	return 0;
    }
  return 1;
}

void simplex(float ** table,int m,int n,int art,int slack,int * type)
{
  int ent =1,lev =1,i,j,index,c1=0,c2=0;
  float min=0,minr=1000000;
  int *basic=(int *)malloc((m+1)*sizeof(int));
  int vars = n+art+slack;
  basic[0]=0;
  float c=1;
  if(*type == 1)
    c=-1;

  for(i=0;i<m;i++)
    {
      if(type[i+1] != 2)
	{
	  basic[i+1] = n+c1+1;
	  c1++;
	  if(type[i+1]==1)
	    c2++;
	}
      else
	{
	  basic[i+1]=n+art+c2+1;
	  c2++;
	}	  
    }
  while(!optimality(table[0],vars+1,c))
    {
      min =0;
      minr=1000000;
      for(i=0;i<vars;i++)
	{
	  if(c*table[0][i+1] < min)
	    {
	      min = c*table[0][i+1];
	      ent = i+1;
	    }
	}
      for(i=0;i<m;i++)
	{
	  if(fabs(table[i+1][ent]) > 0.001 )
	    {
	      float ratio = table[i+1][vars+1]/table[i+1][ent];
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
      for(i=0;i<vars+2;i++)
	{
	  table[index][i]=table[index][i]/q;
	}
      for(i=0;i<m+1;i++)
	{
	  if(i != index)
	    {
	      rowoperation(table[i],table[i],table[index],-1*table[i][ent],vars+2);
	    }
	}
    }

  printf("The optimal value of Z is %.3f at \n",-1*table[0][vars+1]);
  for(i=0;i<m;i++)
    {
      printf("x[%d] = %.3f\n",basic[i+1],table[i+1][vars+1]);
    }
}	      	  	      
   	      
int main()
{
  int m,n,i,j,slack,art,vars,c1,c2;
  slack=0;art=0;c1=0,c2=0;
  printf("Enter the number of equations and variables ");
  scanf("%d %d",&m,&n);
  float *z=(float *)malloc(n*sizeof(float));
  float **mat = allot2d(m,n+1);
  int *type=(int *)malloc((m+1)*sizeof(int));
  printf("Input the equation of Z\n");
  for(i=0;i<n;i++)
    {
      scanf("%f",&z[i]);
    }
  printf("Enter 1 to Maximize 2 to minimize\n");
  scanf("%d",type);
  printf("Enter the original equations\n");
  for(i=0;i<m;i++)
    {
      printf("Enter the %d equation\n",i+1);
      for(j=0;j<n+1;j++)
	{
	  scanf("%f",(j+*(mat+i)));
	}      
      printf("Enter 1 for >= 2 for <= 3 for = \n");
      scanf("%d",type+i+1);
      if(type[i+1] != 2)
	art++;
      if(type[i+1] != 3)
	slack++;
    }
  vars=n+slack+art+1;
  float **table=allot2d(m+1,vars+1);

  //generating table
  float c=1;
  
  if(*type == 1)
    c=-1;
  for(i=0;i<n;i++)
    table[0][i+1]=z[i];

  for(i=0;i<art;i++)
    table[0][i+n+1]=c*M;

  for(i=1;i<=m;i++)
    {
      for(j=0;j<n;j++)
	table[i][j+1]=mat[i-1][j];
      if(type[i]!=2)
	{
	  table[i][n+c1+1]=1;
	  c1++;
	}
      if(type[i]!=3)
	{
	  if(type[i]== 1)
	    table[i][n+art+c2+1]=-1;
	  else
	    table[i][n+art+c2+1]=1;
	  
	  c2++;
	}  
      table[i][vars]=mat[i-1][n];
    }
  printf("\n");
  printmat(table,m+1,vars+1);
  printf("\n");

  for(i=0;i<m;i++)
    {
      if(type[i+1]!=2)
	{
	  rowoperation(table[0],table[0],table[i+1],-c*M,vars+1);
	}
    }
  printmat(table,m+1,vars+1);
  simplex(table,m,n,art,slack,type);
  return 0;
}

