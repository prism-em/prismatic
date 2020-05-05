#include "potentialTests.h"
#include <boost/test/unit_test.hpp>

void test_case1(){
    BOOST_TEST(1==0);
};

void test_case2(){
    int a = 1;
    int b = 2;
    int c = a + b;
};
