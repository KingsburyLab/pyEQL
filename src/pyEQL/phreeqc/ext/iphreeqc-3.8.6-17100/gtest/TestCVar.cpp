#include <gtest/gtest.h>

#include "../src/CVar.hxx"

TEST(TestCVar, CVarCtor)
{
	CVar v;
	ASSERT_EQ(TT_EMPTY, v.type);
}
