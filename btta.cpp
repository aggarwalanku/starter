#include <iostream>
#include <vector>
#include <cstdlib>
#include "matrix.h"

using namespace std;


int main()
{
	matrix V, W, X, Y, Z;
	int i, j, k, row, col, n;//n -> dimension of the block tridiagonal matrix
	float h;//h -> step size
	cout << "enter the step size: ";
	cin >> h;
	n = (1/h)-1;
	row = 2; col = 2;
	
	matrix_3d A, B, C, D, B1, C1, D1, X1;
	A.input(n, row, col, 1, h);
	B.input(n, row, col, 2, h);
	C.input(n, row, col, 3, h);
	D.input(n, row, 1, 4, h);
	
	W = B.transform(0);
	W = ~W;
	X = C.transform(0);
	X = W * X;
	C1.mat.push_back(X.mat);C1.m = 1; C1.n = row; C1.k = col;
	X = D.transform(0);
	X = W * X;
	D1.mat.push_back(X.mat);D1.m = 1; D1.n = row; D1.k = 1;
	for(i=1; i<n; i++)
	{
		//1st equation
		W = A.transform(i);
		X = C1.transform(i-1);
		Y = W * X;
		Z = B.transform(i);
		Z = Z - Y;
		B1.mat.push_back(Z.mat);B1.m = i; B1.n = row; B1.k = col;
		//2nd equation
		W = B1.transform(i-1);//'i-1' as it is getting pushed for the first time
		W = ~W;
		X = C.transform(i);
		Y = W * X;
		C1.mat.push_back(Y.mat);C1.m = i+1; C1.n = row; C1.k = col;
		//3rd equation
		X = D1.transform(i-1);
		Y = A.transform(i);
		Z = Y * X;
		V = D.transform(i);
		V = V - Z;
		W = W * V;
		D1.mat.push_back(W.mat);D1.m = i+1; D1.n = row; D1.k = 1;
	}
	
	cout << "result\n";
	//displaying the result
	X = D1.transform(n-1);
	X1.mat.push_back(X.mat);X1.m = 1; X1.n = row; X1.k = 1;
	for(i=n-2; i>=0; i--)
	{
		Y = C1.transform(i);
		Z = Y * X;
		W = D1.transform(i);
		W = W - Z;
		X1.mat.push_back(W.mat);X1.m = n-i; X1.n = row; X1.k = 1;
		X = X1.transform(n-i-1);
	}

	for(i=n-1; i>=0; i--)
	{
		X = X1.transform(i);
		cout << h*(n-i) << " -> " << X.mat[0][0] << endl;
	}
	return 0;
}
