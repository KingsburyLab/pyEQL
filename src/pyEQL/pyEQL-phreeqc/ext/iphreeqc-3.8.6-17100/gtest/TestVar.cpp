#include <gtest/gtest.h>

#include "../src/Var.h"

TEST(TestVar, VarInit)
{
	VAR v;
	::VarInit(&v);
	ASSERT_EQ(TT_EMPTY, v.type);
}
