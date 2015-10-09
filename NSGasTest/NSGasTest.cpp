#include <gtest/gtest.h>

using namespace std;
using namespace ::testing;


TEST(nsgas, main_test)
{
}

int main(int ac, char* av [])
{
	testing::InitGoogleTest(&ac, av);
	return RUN_ALL_TESTS();
}
