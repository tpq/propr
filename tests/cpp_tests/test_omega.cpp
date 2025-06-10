#include "catch2/catch.hpp"

TEST_CASE("omega basic functionality", "[omega]") {
    SECTION("Test section 1") {
        REQUIRE(true);
    }

    SECTION("Test section 2") {
        REQUIRE(1 == 1);
    }
}
