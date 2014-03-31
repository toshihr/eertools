#include <iostream>
#include <string>
#include <CppUTest/CommandLineTestRunner.h>
#include "eertools.h"

using namespace std;
using namespace eertools;

TEST_GROUP(EntropyTestGroup) {
	TEST_SETUP() {
	}
 
	TEST_TEARDOWN() {
	}
};

TEST(EntropyTestGroup, TestProbabilityDistribution) {
	const double eps = __FLT_EPSILON__;
	const string s = "abab";
	vec_double p1;

	calcProb(p1, s);
	const double e = getEntropy(p1);

	DOUBLES_EQUAL(e, 1.0, eps);  
} 

int main(int argc, char** argv) {
  return RUN_ALL_TESTS(argc, argv);
}