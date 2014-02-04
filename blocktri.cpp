#include<iostream>
#include<vector>
using namespace std;
typedef std::vector<float> row;
typedef std::vector<row > dat;

class matrix {
public:
  dat data;
  int rows,col;
  friend ostream& operator<<(ostream& out,matrix &A);
  friend istream& operator>>(istream& in,matrix &A);
  matrix operator+(const matrix &a);
  matrix operator*(const matrix &a);
  matrix(int a,int b) : data(a, row(b)), rows(a), col(b){}
  row& operator[](int x) { return data[x];}
  const row& operator[](int x) const { return data[x];}
  matrix operator*(float k);

};
matrix matrix::operator+(const matrix &a)
{
  matrix c=a;
  for(int i=0;i<c.rows;i++)
    {
      for(int j=0;j<c.col;j++)
	{
	  c[i][j]+=(*this)[i][j];
	}
    }
  return c;
}
matrix matrix::operator*(const matrix &a)
{
  matrix c(this->rows,a.col);
  for(int i=0;i<this->rows;i++)
    for(int j=0;j<a.col;j++)
      for(int k=0;k<a.rows;k++)
	c[i][j]+=(*this)[i][k] * (a[k][j]);
  return c;
} 
matrix matrix::operator*(float k)
{
  matrix c=*this;
  for(int i=0;i<c.rows;i++)
    for(int j=0;j<c.col;j++)
      c[i][j]=c[i][j]*k;
  return c;
}
ostream& operator<<(ostream& out ,matrix &a)
{
  int m,n;
  m=a.rows;n=a.col;
  for(int i=0;i<m;i++)
    {
      for(int j=0;j<n;j++)
	{
	  out<<a[i][j]<<" ";
	}
      out<<"\n";
    }
  return out;
}
istream& operator>>(istream& in ,matrix &a)
{
  for(int i=0;i<a.rows;i++)
    for(int j=0;j<a.col;j++)
      in>>a[i][j];
  return in;
}

void rowoperation(row &dest,row a,row b,float coeff,int n)
{
  int i;
  for(i=0;i<n;i++)
    {
      dest[i]=a[i]+( b[i]*coeff);
    }
}

matrix identity(int n)
{
  matrix mat(n,n);
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

matrix inverse(matrix cmat)
{
  matrix mat =cmat;
  int n = mat.rows;
  matrix inv = identity(n);
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
}
	
void block_thomas(float * coeff,int n)
{
  float h = 1/(float)n;
  // generating matrices from the equation coefficients
  std::vector<matrix> A(n-2,matrix(2,2));
  std::vector<matrix> B(n-1,matrix(2,2));
  std::vector<matrix> C(n-2,matrix(2,2));
  
  matrix d(2,1);
  d[0][0]=coeff[4];d[1][0]=0;
  std::vector<matrix> D(n-1,d);
  std::vector<matrix> X(n-1,matrix(2,1));

  for(int i=0;i<n-2;i++)
    {
      A[i][0][0]= (coeff[0]/(h*h))-(coeff[1]/(2*h));
      A[i][0][1]=0;
      A[i][1][0]=0;
      A[i][1][1]=-0.5/h;
    }
  for(int i=0;i<n-1;i++)
    {
      B[i][0][0]= (-2*coeff[0]/(h*h))+ coeff[2];
      B[i][0][1]=coeff[3];
      B[i][1][0]= -1;
      B[i][1][1]= 0;
    }
  for(int i=0;i<n-2;i++)
    {
      C[i][0][0]= (coeff[0]/(h*h))+(coeff[1]/(2*h));
      C[i][0][1]=0;
      C[i][1][0]=0;
      C[i][1][1]=0.5/h;
    }
  D[0][0][0]=coeff[4]-((coeff[0]/(h*h))-(coeff[1]/(2*h)))*coeff[6];
  D[0][1][0]=0.5/h*coeff[5];

  D[n-2][0][0]=coeff[4]-((coeff[0]/(h*h))+(coeff[1]/(2*h)))*coeff[8];
  D[n-2][1][0]=-0.5/h*coeff[7];
  //thomas algo for block tri diagonal
  matrix temp = inverse(B[0]);
  C[0]=temp*C[0];
  D[0]=temp*D[0];
  for(int i=1;i<n-2;i++)
    {
      B[i]=B[i]+((A[i-1]*C[i-1])*(-1));
      temp=inverse(B[i]);
      C[i]=temp*C[i];
      D[i]=temp*(D[i]+A[i-1]*D[i-1]*-1);
    }
  B[n-2]=B[n-2]+((A[n-3]*C[n-3])*(-1));
  temp=inverse(B[n-2]);
  D[n-2]=temp*(D[n-2]+A[n-3]*D[n-3]*-1);
  X[n-2]=D[n-2];
  for(int i=n-3;i>=0;i--)
    {
      X[i]= D[i]+C[i]*X[i+1]*-1;
    }

  for(int i=0;i<n-1;i++)
    {
      cout<<(i+1)*h<<" "<<X[i][1][0]<<"\n";
    }

} 
      
      
  
int main()
{
  /* int i,j,k;
  cin>>i>>j>>k;
  matrix a(i,j);
  //matrix b(j,k);
  cin>>a;
    matrix c = a*-2;
  cout<<c;
  return 0;
  */
  float coeff[9];
  int i,n,p,q,r;

  cout<<"Enter the coefficients of y\"', y\",y',y,x^2 \n";
  cin>>coeff[0]>>coeff[1]>>coeff[2]>>coeff[3]>>coeff[4];
  cout<<"Enter  y(0),y'(0),y(1),y'(1) and three choices of number of steps \n";
  cin>>coeff[5]>>coeff[6]>>coeff[7]>>coeff[8]>>p>>q>>r;
  cout<<" The equation is :- \n";
  cout<<coeff[0]<<"y\"' + "<<coeff[1]<<"y\"+ "<<coeff[2]<<"y'+"<<coeff[3]<<"y = "<<coeff[4]<<"x^2 \n";
  cout<<" y(0) = "<<coeff[5]<<"  y'(0) = "<<coeff[6]<<" y(1) = "<<coeff[7]<<" y'(1) = "<<coeff[8]<<" \n";
  block_thomas(coeff,p);
  cout<<"\n";
  block_thomas(coeff,q);
  cout<<"\n";
  block_thomas(coeff,r);

  return 0;
}
      
