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

float ** identity(int n)
{
  float **mat =allot2d(n,n);
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
float ** inverse(float** cmat,int n)
{
  float**mat=copy(cmat,n);
  float **inv= identity(n);
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
  return inv;
	
}/*
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
*/
void subs(float **B,int col_num,float * col,int n)
{
  int i;
  for(i=0;i<n;i++)
    {
      B[i][col_num]=col[i];
    }
}
float * col_Pj(float **table,int col_num,int n)
{
  float * result = (float *)malloc(n*sizeof(float));
  int i;
  for(i=0;i<n;i++)
    {
      result[i]=table[i][col_num-1];
    }
  return result;
}

float * colmat(float * col,float ** mat,int n,int flag)
{
  float * result = (float *)malloc(n*sizeof(float));
  int i,j,k;
  for(i=0;i<n;i++)
    {
      float sum=0;
      for(j=0;j<n;j++)
	{
	  if(flag)
	    sum += col[j]*mat[j][i];
	  else
	    sum += col[j]*mat[i][j];
	}
      result[i]=sum;
    }
  return result;
}
float colcol(float *a,float *b, int n)
{
  int i;float sum =0;
  for(i=0;i<n;i++)
    sum+=a[i]*b[i];
  return sum;
}
void revised_simplex(float ** table,int m,int n,float *z)
{
  int *nbasic = (int *)malloc(n*sizeof(int));
  int *basic = (int *)malloc(m*sizeof(int));
  float *cnbasic = (float *)malloc(n*sizeof(float));
  float *cbasic = (float *)malloc(m*sizeof(float));
  float **B=identity(m);
  int i,j,temp;
  float* xB;
  float *b = (float *)malloc(m*sizeof(float));
  for(i=0;i<m;i++)
    b[i]=table[i][m+n];
  //Setting initial basic and non basic variables
  for(i=0;i<n;i++)
    nbasic[i]=i+1;
  for(i=0;i<m;i++)
    basic[i]=m+i+1;
  for(i=0;i<n;i++)
    cnbasic[i]=z[i];
  for(i=0;i<m;i++)
    cbasic[i]=0;  
  j=0;
  while(1)
    {
      printf("\nITERATION NUMBER %d\n",j+1);
      printf("\n The basic variables are \n");
      for(i=0;i<m;i++)
	{
	  printf("%d ",basic[i]);
	}
      printf("\n The nonbasic variables are \n");
        for(i=0;i<n;i++)
	{
	  printf("%d ",nbasic[i]);
	}
      
      float **iB=inverse(B,m);

      // Selecting entering variable
      float min=0;int ent =1;
      float *Y=colmat(cbasic,iB,m,1);
      for(i=0;i<n;i++)
	{
	  float *Pi=col_Pj(table,nbasic[i],m);
	  float YPi=colcol(Y,Pi,m);
	  float s= YPi-cnbasic[i];
	  if( s < min)
	    {
	      min = s;
	      ent =i;
	    }	  
	}
      xB= colmat(b,iB,m,0);
      if (min == 0) break;
      printf("\n The entering variable is x%d ",nbasic[ent]);
      
      //Selecting Leaving variable
      float *Pj=col_Pj(table,nbasic[ent],m);
      float * Aj=colmat(Pj,iB,m,0);
      min =100000;
      int leave =1;
      for(i=0;i<m;i++)
	{
	  float s= xB[i]/Aj[i];
	  if(s < min && s > 0)
	    {
	      min = s;
	      leave =i;
	    }
	}
      printf("\n The leaving  variable is x%d ",basic[leave]);

      //Computing new B
      temp=nbasic[ent];
      nbasic[ent]=basic[leave];
      basic[leave]=temp;
      Pj = col_Pj(table,temp,m);
      temp=cbasic[leave];
      cbasic[leave]=cnbasic[ent];
      cnbasic[ent]=temp;
      subs(B,leave,Pj,m);
      printf("\n New B is\n");
      printmat(B,m,m);
      printf("\n");
      j++;
    }
  printf("\nWe stop here.\n");
  printf("\n The optimal value of Z is %f\nat \n",colcol(xB,cbasic,m));
  for(i=0;i<m;i++)
    {
      printf("x%d = %f\n",basic[i],xB[i]);
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
  float **table=allot2d(m,m+n+1);
  //generating table
  for(i=0;i<m;i++)
    {
      for(j=0;j<n;j++)
	table[i][j]=mat[i][j];
      table[i][n+i]=1;
      table[i][m+n]=mat[i][n];
    }
  revised_simplex(table,m,n,z);
}

