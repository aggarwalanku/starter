#include<stdio.h>
#include<iostream>
long int binary(int n) {
	long int rem,i=1,bin=0;
	do {
		rem=(n%2);
		bin+=(i*rem);
		n=n/2;
		i=i*10;
	}while(n>0);
	std::cout<<"BIN = "<<bin<<"\n";
	return bin;
}
int main()
{ 
 binary(1);
	return 0;
}
