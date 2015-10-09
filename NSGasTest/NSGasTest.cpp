#include <gtest/gtest.h>
#include <iostream>

using namespace std;
using namespace ::testing;

int factorian(int n)
{
	if (n == 0 || n == 1)
		return 1;
	else
		return n * factorian(n - 1);
}

TEST(factorianTest, simpleTest)
{
	ASSERT_EQ(1, factorian(0));
	ASSERT_EQ(1, factorian(1));
	ASSERT_EQ(2, factorian(2));
}

int main(int ac, char* av [])
{
	testing::InitGoogleTest(&ac, av);
	return RUN_ALL_TESTS();
}
