#include <boost/test/included/unit_test.hpp>
#include "potentialTests.h"

namespace bt = boost::unit_test;

bt::test_suite* init_unit_test_suite( int /*argc*/, char* /*argv*/[] )
{
    bt::framework::master_test_suite().p_name.value = "Prismatic Test Suite";

    bt::test_suite* ts1 = BOOST_TEST_SUITE( "potentialTests" );
    ts1->add( BOOST_TEST_CASE( &test_case1 ) );
    ts1->add( BOOST_TEST_CASE( &test_case2 ) );

    // bt::test_suite* ts2 = BOOST_TEST_SUITE( "test_suite2" );
    // ts2->add( BOOST_TEST_CASE( &test_case3 ) );
    // ts2->add( BOOST_TEST_CASE( &test_case4 ) );

    bt::framework::master_test_suite().add( ts1 );
    // bt::framework::master_test_suite().add( ts2 );

    return 0;
}