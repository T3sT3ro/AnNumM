#include <stdio.h>

int main(){
	double x = 1.0;
	double exact = 1.0;

	while(x < 2.0) {
	volatile double div = 1.0/x;
	volatile double ret = x*div;
		if (*(long long*)(&ret) != *(long long*)(&exact))
			printf("%.64lf\n", x);
		(*(long long*)(&x))++;
	}
	return 0;
}
