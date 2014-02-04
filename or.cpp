#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

int count =0;

void printmatrix(float** p,int x,int y)            
{ int i,j; 
  for( i=0;i<x;i++)
   { for(j=0;j<y;j++)
       { printf(" %.2f ",p[i][j]);
       }
     printf("\n");
   }
}
float** allocate(int x,int y)                                   
{ int i;
  float **c=(float **)malloc(x*sizeof(float*));  
  for( i=0;i<x;i++)
    { c[i] = (float *)malloc(y*sizeof(float));
    }
  return c;
}

float arraymult(float *a,float *b,int n)
{ 
  float sum=0;
  int i;
  for(i=0;i<n;i++)
    { sum += a[i]*b[i];
    }
 return sum;
}

int factorial(int x)
{
  if(x==1||x==0)
    return(1);
  else
    return(x*factorial(x-1));
}
int ncr(int n,int r)
{ 
  return factorial(n)/((factorial(r))*(factorial(n-r)));    
}
void fill(float *a,int index,int num,int n,float ** store)
{
  float *copy = (float *)malloc(n*sizeof(float));
  memcpy(copy,a,n*sizeof(float));

  copy[index]=1;
  if(num == 1)
    {
      int j;
      for(j=0;j<n;j++)
	{
	  //printf(" %.2f ",copy[j]);
	  store[count][j]=copy[j];
       }
      count++;
      //printf("\n");
    }
  else
    {
      fill(copy,index+1,num-1,n,store);
    }
  
  if(n-index-1 >= num)
    {
      fill(a,index+1,num,n,store);
    }

}

int gauss(float **a,int n,float *x)
{
   int i,j,k,max;
   double tmp;

   for (i=0;i<n;i++) {
      max = i;
      for (j=i+1;j<n;j++) {
         if (fabs(a[j][i]) > fabs(a[max][i]))
            max = j;
      }
      
      for (k=i;k<n+1;k++) {
         tmp = a[i][k];
         a[i][k] = a[max][k];
         a[max][k] = tmp;
      }
       if (fabs(a[i][i]) < 0.001)
         return 0;
      
      for (j=i+1;j<n;j++) {
         for (k=n;k>=i;k--) {
            a[j][k] -= a[i][k] * a[j][i] / a[i][i];
         }
      }
   }
   //printmatrix(a,n,n);
   for (j=n-1;j>=0;j--) {
      tmp = 0;
      for (k=j+1;k<n;k++)
         tmp += a[j][k] * x[k];
      x[j] = (a[j][n] - tmp) / a[j][j];
   }

   return 1;
}

void printsolution(float **sols,int m,int n,int *type,float *optim,int max)
{
  int i =0,j,maxc=0,num=0;
  float zz=0;
  int comb = ncr(n,n-m);
  printf("\n");
  for(;i<n;i++)
    printf("x%d\t",i+1);
  printf(" Z      Type \n");
  for(i=0;i<comb;i++)
    {
      int flag=0;
      if(type[i] == 1)
	{
	  num++;
	  for(j=0;j<n;j++)
	    {
	      if(sols[i][j] < 0)
		flag =1;
	      printf("%.2f\t",sols[i][j]);
	    }
	  float z = arraymult(sols[i],optim,n);
	  if(z > zz)
	    {
	      zz= z;
	      maxc =i;
	    }
	  printf("%.2f\t",z);
	  if(flag == 1)
	    printf("Not feasible\n");
	  else
	    printf("Feasible\n");
	}
    }
  if(num != 0)
    {
      printf("The optimal solution is \n");
      for(i=0;i<n;i++)
	printf("x%d=%f\n",i+1,sols[maxc][i]);
      printf("Z=%.2f\n",zz);
    }
  else
    {
      printf("Given system is either inconsistent or has infinite solutions\n");
    }
}
int main()
{
  int  m,n,p,i,j,comb,max;
  printf("Enter number of equations and variables\n");
  scanf("%d %d",&m,&n);
  p=n+1;
  float *optim=(float *)malloc(n*sizeof(float));
  float *blank=(float *)malloc(n*sizeof(float));
  float **mat = allocate(m,p);
  printf("Enter the coefficients of the equations\n");
  scanf("%*c");
  for(i=0;i<m;i++)
    {
      printf("Enter equation %d \n",i+1);
      for(j=0;j<p;j++)
	{
	  scanf("%f",(j+*(mat+i)));
	}
    }
  printf("The coefficient matrix is \n");
  printmatrix(mat,m,p);
  printf("\n");
  printf("Enter the function to be optimized \n");
  for(j=0;j<n;j++)
    {
      scanf("%f",optim+j);
    }
  printf("Enter 1 to maximize, 0 to minimize \n");
  scanf("%d",&max);
  
  comb =ncr(n,n-m);
  float **sols = allocate(comb,n);
  float **store = allocate(comb,n);
  float **temp = allocate(m,m+1);
  int *type=(int *)malloc(comb*sizeof(int));
  fill(blank,0,m,n,store);
  for( i=0;i<comb;i++)
    {
      int k=0,l;
      for(j=0;j<n;j++)
	{
	  if(fabs(store[i][j] - 1) < 0.001)
	    {
	     
	      for(l=0;l<m;l++)
		{
		  temp[l][k]=mat[l][j];
		}
	      k++;
	    }
	}
      for(l=0;l<m;l++)
	{
	  temp[l][k]=mat[l][n];
	}
      //printmatrix(temp,m,m+1);
      //  printf("\n");
      float * fs = (float *)malloc(m*sizeof(float));
     type[i]= gauss(temp,m,fs);
      k =0;
      for(j=0;j<n;j++)
	{
	  if(store[i][j] < 0.001)
	    {
	      sols[i][j]=0;
	    }
	  else
	    {
	      sols[i][j]=fs[k];
	      k++;
	    }
	}
    }
  //  printmatrix(sols,comb,n);
  printsolution(sols,m,n,type,optim,max);
  return 0;
}

	      
	      
	  
	  
      
