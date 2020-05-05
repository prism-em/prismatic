#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

namespace bt = boost::unit_test;

BOOST_AUTO_TEST_CASE( filler )
{
    BOOST_TEST( 1== 1);
}