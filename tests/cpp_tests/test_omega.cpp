#include <catch.hpp>

// Include headers for the functionality you want to test
// #include "your_header.h"

TEST_CASE("omega basic functionality", "[omega]") {
    // Add your test cases here
    SECTION("Test section 1") {
        // Test code here
        REQUIRE(true);
    }

    SECTION("Test section 2") {
        // More test code here
        REQUIRE(1 == 1);
    }
}
