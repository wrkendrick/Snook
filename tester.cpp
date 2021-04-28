#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <numeric>
#include <vector>

std::vector<int> vec_1;

int main() {
	int a[] = {1,2,3};
	int b[] = {4,5,6};
	vec_1.assign(a,a+3);
	
	std::cout << "HERE\n";
	std::vector<int>* vecPtr = &vec_1;
	for (int j=0; j<vecPtr->size(); j++) {
		std::cout << vecPtr->at(j) << "\n";
	}
	std::cout << "RUN FINISHED" << "\n";
}