#include<iostream>

#include"Array.h"
#include"ArrayCross9.h"
#include"ArrayStaggered.h"

using namespace std;

int main()
{
	Array phi(2, 3, 2, 1.0);

	ArrayStaggered u(64, 32, 1.0);

	ArrayCross9 phin(102, 3, 2, 1.0);
}