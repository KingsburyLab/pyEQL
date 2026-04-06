#include <gtest/gtest.h>

#include "../src/CSelectedOutput.hxx"
#include "IPhreeqc.hpp"
#include "Phreeqc.h"

#if defined(_WIN32)
#define strdup _strdup
#endif

TEST(TestSelectedOutput, TestEmpty)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());
	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());
}

TEST(TestSelectedOutput, TestSinglePushBack)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());
	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	CVar v(7.0);
	ASSERT_EQ(0, co.PushBack("pH", v));

	ASSERT_EQ((size_t)1, co.GetColCount());
	// row count doesn't change until EndRow is called
	ASSERT_EQ((size_t)1, co.GetRowCount());
	ASSERT_EQ(0, co.EndRow());

	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());
#if defined(_DEBUG)
	co.Dump("TestSinglePushBack");
#endif
}

TEST(TestSelectedOutput, TestMultiplePushBack)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());
	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	CVar v1(7.0);
	ASSERT_EQ(0, co.PushBack("pH", v1));

	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)1, co.GetRowCount());

	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	CVar v2(8.0);
	ASSERT_EQ(0, co.PushBack("pH", v2));

	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)3, co.GetRowCount());
#if defined(_DEBUG)
	co.Dump("TestMultiplePushBack");
#endif
}

TEST(TestSelectedOutput, TestNewHeadingsPushBack)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	CVar v1(7.0);
	ASSERT_EQ(0, co.PushBack("pH", v1));

	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)1, co.GetRowCount());

	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	CVar v2(8.0);
	ASSERT_EQ(0, co.PushBack("pH", v2));

	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	CVar v3(9.0);
	ASSERT_EQ(0, co.PushBack("user_pH", v3));

	ASSERT_EQ((size_t)2, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());


	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)2, co.GetColCount());
	ASSERT_EQ((size_t)3, co.GetRowCount());
#if defined(_DEBUG)
	co.Dump("TestNewHeadingsPushBack");
#endif
}

TEST(TestSelectedOutput, TestPushBackDouble)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());


	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	ASSERT_EQ(0, co.PushBackDouble("pH", 7.0));

	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)1, co.GetRowCount()); // heading

	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	CVar v;
	ASSERT_EQ(TT_EMPTY, v.type);
	ASSERT_EQ(VR_OK, co.Get(0, 0, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("pH"), std::string(v.sVal));

	CVar vval;
	ASSERT_EQ(TT_EMPTY, vval.type);
	ASSERT_EQ(VR_OK, co.Get(1, 0, &vval));
	ASSERT_EQ(TT_DOUBLE, vval.type);
	ASSERT_EQ(7.0, vval.dVal);
}

TEST(TestSelectedOutput, TestPushBackLong)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	ASSERT_EQ(0, co.PushBackLong("Sim", 2));

	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)1, co.GetRowCount()); // heading plus first row

	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	CVar v;
	ASSERT_EQ(TT_EMPTY, v.type);
	ASSERT_EQ(VR_OK, co.Get(0, 0, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Sim"), std::string(v.sVal));

	CVar vval;
	ASSERT_EQ(TT_EMPTY, vval.type);
	ASSERT_EQ(VR_OK, co.Get(1, 0, &vval));
	ASSERT_EQ(TT_LONG, vval.type);
	ASSERT_EQ(2l, vval.lVal);
}

TEST(TestSelectedOutput, TestPushBackString)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	ASSERT_EQ(0, co.PushBackString("state", "i_soln"));

	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)1, co.GetRowCount()); // heading

	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	CVar v;
	ASSERT_EQ(TT_EMPTY, v.type);
	ASSERT_EQ(VR_OK, co.Get(0, 0, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("state"), std::string(v.sVal));

	CVar vval;
	ASSERT_EQ(TT_EMPTY, vval.type);
	ASSERT_EQ(VR_OK, co.Get(1, 0, &vval));
	ASSERT_EQ(TT_STRING, vval.type);
	ASSERT_EQ(std::string("i_soln"), std::string(vval.sVal));
}

TEST(TestSelectedOutput, TestPushBackEmpty)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	ASSERT_EQ(0, co.PushBackEmpty("Empty"));

	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)1, co.GetRowCount()); // heading

	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	CVar v;
	ASSERT_EQ(TT_EMPTY, v.type);
	ASSERT_EQ(VR_OK, co.Get(0, 0, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("Empty"), std::string(v.sVal));

	CVar vval;
	ASSERT_EQ(TT_EMPTY, vval.type);
	ASSERT_EQ(VR_OK, co.Get(1, 0, &vval));
	ASSERT_EQ(TT_EMPTY, vval.type);
}

TEST(TestSelectedOutput, TestDuplicateHeadings)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	ASSERT_EQ(0, co.PushBackDouble("pH", 7.0));

	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)1, co.GetRowCount()); // heading

	// overwrite pH with 8.0
	//
	ASSERT_EQ(0, co.PushBackDouble("pH", 8.0));

	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)1, co.GetRowCount()); // heading

	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	CVar v;
	ASSERT_EQ(TT_EMPTY, v.type);
	ASSERT_EQ(VR_OK, co.Get(0, 0, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("pH"), std::string(v.sVal));


	CVar vval;
	ASSERT_EQ(TT_EMPTY, vval.type);
	ASSERT_EQ(VR_OK, co.Get(1, 0, &vval));
	ASSERT_EQ(TT_DOUBLE, vval.type);
	ASSERT_EQ(8.0, vval.dVal);
}

TEST(TestSelectedOutput, TestEndRow)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	ASSERT_EQ(0, co.PushBackDouble("pH", 7.0));

	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)1, co.GetRowCount()); // heading

	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	CVar v;
	ASSERT_EQ(TT_EMPTY, v.type);
	ASSERT_EQ(VR_OK, co.Get(0, 0, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("pH"), std::string(v.sVal));

	CVar vval;
	ASSERT_EQ(TT_EMPTY, vval.type);
	ASSERT_EQ(VR_OK, co.Get(1, 0, &vval));
	ASSERT_EQ(TT_DOUBLE, vval.type);
	ASSERT_EQ(7.0, vval.dVal);

	ASSERT_EQ(0, co.PushBackDouble("pH", 8.0));
	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)3, co.GetRowCount());

	CVar vval3;
	ASSERT_EQ(TT_EMPTY, vval3.type);
	ASSERT_EQ(VR_OK, co.Get(1, 0, &vval3));
	ASSERT_EQ(TT_DOUBLE, vval3.type);
	ASSERT_EQ(7.0, vval3.dVal);

	CVar vval2;
	ASSERT_EQ(TT_EMPTY, vval2.type);
	ASSERT_EQ(VR_OK, co.Get(2, 0, &vval2));
	ASSERT_EQ(TT_DOUBLE, vval2.type);
	ASSERT_EQ(8.0, vval2.dVal);
}

TEST(TestSelectedOutput, TestEndRow2)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	ASSERT_EQ(0, co.PushBackDouble("pH", 6.0));
	ASSERT_EQ(0, co.PushBackDouble("pH", 7.0));
	ASSERT_EQ(0, co.PushBackDouble("pH", 8.0));
	ASSERT_EQ(0, co.PushBackDouble("pH", 9.0));

	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	CVar v;
	ASSERT_EQ(TT_EMPTY, v.type);
	ASSERT_EQ(VR_OK, co.Get(0, 0, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("pH"), std::string(v.sVal));

	CVar vval;
	ASSERT_EQ(TT_EMPTY, vval.type);
	ASSERT_EQ(VR_OK, co.Get(1, 0, &vval));
	ASSERT_EQ(TT_DOUBLE, vval.type);
	ASSERT_EQ(9.0, vval.dVal);   // dups get overwritten

	ASSERT_EQ(0, co.PushBackDouble("pH", 8.0));
	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)3, co.GetRowCount());

	CVar vval3;
	ASSERT_EQ(TT_EMPTY, vval3.type);
	ASSERT_EQ(VR_OK, co.Get(1, 0, &vval3));
	ASSERT_EQ(TT_DOUBLE, vval3.type);
	ASSERT_EQ(9.0, vval3.dVal);

	CVar vval2;
	ASSERT_EQ(TT_EMPTY, vval2.type);
	ASSERT_EQ(VR_OK, co.Get(2, 0, &vval2));
	ASSERT_EQ(TT_DOUBLE, vval2.type);
	ASSERT_EQ(8.0, vval2.dVal);
}


TEST(TestSelectedOutput, TestTooManyHeadings)
{
	// COMMENT: {8/26/2013 4:12:03 PM}	IPhreeqc p;
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ((size_t)0, p.PtrSelectedOutput->GetColCount());
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ((size_t)0, p.PtrSelectedOutput->GetRowCount());
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	p.PtrSelectedOutput->Clear();
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ((size_t)0, p.PtrSelectedOutput->GetColCount());
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ((size_t)0, p.PtrSelectedOutput->GetRowCount());
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	// USER_PUNCH
	// COMMENT: {8/26/2013 4:12:03 PM}	// -headings 1.name 1.type 1.moles
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	p.PhreeqcPtr->n_user_punch_index        = 0;
	// COMMENT: {8/26/2013 4:12:03 PM}	p.PhreeqcPtr->UserPunch_map[1]          = UserPunch();
	// COMMENT: {8/26/2013 4:12:03 PM}	p.PhreeqcPtr->current_user_punch        = &(p.PhreeqcPtr->UserPunch_map[1]);
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	std::vector< std::string > headings;
	// COMMENT: {8/26/2013 4:12:03 PM}	headings.push_back("1.name");
	// COMMENT: {8/26/2013 4:12:03 PM}	headings.push_back("1.type");
	// COMMENT: {8/26/2013 4:12:03 PM}	headings.push_back("1.moles");
	// COMMENT: {8/26/2013 4:12:03 PM}	p.PhreeqcPtr->UserPunch_map[1].Set_headings(headings);
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(0,   p.EndRow());
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ((size_t)3, p.PtrSelectedOutput->GetColCount());
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ((size_t)2, p.PtrSelectedOutput->GetRowCount());
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}#if defined(_DEBUG)
	// COMMENT: {8/26/2013 4:12:03 PM}	p.PtrSelectedOutput->Dump("TestTooManyHeadings");
	// COMMENT: {8/26/2013 4:12:03 PM}#endif
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	// clean up headings
	// COMMENT: {8/26/2013 4:12:03 PM}	p.PhreeqcPtr->UserPunch_map[1].Get_headings().empty();
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	CVar head0, head1, head2;
	// COMMENT: {8/26/2013 4:12:03 PM}	CVar val0, val1, val2;
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(VR_OK, p.PtrSelectedOutput->Get(0, 0, &head0));
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(VR_OK, p.PtrSelectedOutput->Get(0, 1, &head1));
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(VR_OK, p.PtrSelectedOutput->Get(0, 2, &head2));
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(VR_OK, p.PtrSelectedOutput->Get(1, 0, &val0));
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(VR_OK, p.PtrSelectedOutput->Get(1, 1, &val1));
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(VR_OK, p.PtrSelectedOutput->Get(1, 2, &val2));
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(TT_STRING, head0.type);
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(TT_STRING, head1.type);
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(TT_STRING, head2.type);
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(TT_EMPTY, val0.type);
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(TT_EMPTY, val1.type);
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(TT_EMPTY, val2.type);
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(std::string("1.name"), std::string(head0.sVal));
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(std::string("1.type"), std::string(head1.sVal));
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(std::string("1.moles"), std::string(head2.sVal));
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(0, p.PtrSelectedOutput->PushBackLong("sim", 1));
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(0, p.PtrSelectedOutput->PushBackString("state", "i_soln"));
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(0, p.PtrSelectedOutput->PushBackLong("soln", 22));
	// COMMENT: {8/26/2013 4:12:03 PM}
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ(0,  p.PtrSelectedOutput->EndRow());
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ((size_t)6, p.PtrSelectedOutput->GetColCount());
	// COMMENT: {8/26/2013 4:12:03 PM}	ASSERT_EQ((size_t)3, p.PtrSelectedOutput->GetRowCount());
	// COMMENT: {8/26/2013 4:12:03 PM}#if defined(_DEBUG)
	// COMMENT: {8/26/2013 4:12:03 PM}	p.PtrSelectedOutput->Dump("TestTooManyHeadings");
	// COMMENT: {8/26/2013 4:12:03 PM}#endif
}

TEST(TestSelectedOutput, TestNotEnoughHeadings)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	// USER_PUNCH
	// -headings 1.name 1.type 1.moles

	ASSERT_EQ(0, co.PushBackLong("sim", 1));
	ASSERT_EQ(0, co.PushBackString("state", "i_soln"));
	ASSERT_EQ(0, co.PushBackLong("soln", 22));

	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)3, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());
#if defined(_DEBUG)
	co.Dump("TestNotEnoughHeadings");
#endif

	ASSERT_EQ(0, co.PushBackLong("sim", 2));
	ASSERT_EQ(0, co.PushBackString("state", "react"));
	ASSERT_EQ(0, co.PushBackLong("soln", 23));

	ASSERT_EQ(0, co.PushBackEmpty("no_heading_1"));
	ASSERT_EQ(0, co.PushBackEmpty("no_heading_2"));
	ASSERT_EQ(0, co.PushBackEmpty("no_heading_3"));

#if defined(_DEBUG)
	co.Dump("TestNotEnoughHeadings");
#endif

	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)6, co.GetColCount());
	ASSERT_EQ((size_t)3, co.GetRowCount());

	CVar head0, head1, head2, head3, head4, head5;
	CVar val0, val1, val2, val3, val4, val5;

	ASSERT_EQ(VR_OK, co.Get(0, 0, &head0));
	ASSERT_EQ(VR_OK, co.Get(0, 1, &head1));
	ASSERT_EQ(VR_OK, co.Get(0, 2, &head2));
	ASSERT_EQ(VR_OK, co.Get(0, 3, &head3));
	ASSERT_EQ(VR_OK, co.Get(0, 4, &head4));
	ASSERT_EQ(VR_OK, co.Get(0, 5, &head5));

	ASSERT_EQ(VR_OK, co.Get(1, 0, &val0));
	ASSERT_EQ(VR_OK, co.Get(1, 1, &val1));
	ASSERT_EQ(VR_OK, co.Get(1, 2, &val2));
	ASSERT_EQ(VR_OK, co.Get(1, 3, &val3));
	ASSERT_EQ(VR_OK, co.Get(1, 4, &val4));
	ASSERT_EQ(VR_OK, co.Get(1, 5, &val5));

	ASSERT_EQ(TT_STRING, head0.type);
	ASSERT_EQ(TT_STRING, head1.type);
	ASSERT_EQ(TT_STRING, head2.type);
	ASSERT_EQ(TT_STRING, head3.type);
	ASSERT_EQ(TT_STRING, head4.type);
	ASSERT_EQ(TT_STRING, head5.type);

	ASSERT_EQ(TT_LONG, val0.type);
	ASSERT_EQ(TT_STRING, val1.type);
	ASSERT_EQ(TT_LONG, val2.type);
	ASSERT_EQ(TT_EMPTY, val3.type);
	ASSERT_EQ(TT_EMPTY, val4.type);
	ASSERT_EQ(TT_EMPTY, val5.type);

	ASSERT_EQ(std::string("sim"), std::string(head0.sVal));
	ASSERT_EQ(std::string("state"), std::string(head1.sVal));
	ASSERT_EQ(std::string("soln"), std::string(head2.sVal));
	ASSERT_EQ(std::string("no_heading_1"), std::string(head3.sVal));
	ASSERT_EQ(std::string("no_heading_2"), std::string(head4.sVal));
	ASSERT_EQ(std::string("no_heading_3"), std::string(head5.sVal));

	ASSERT_EQ(1l, val0.lVal);
	ASSERT_EQ(std::string("i_soln"), std::string(val1.sVal));
	ASSERT_EQ(22l, val2.lVal);
}

TEST(TestSelectedOutput, TestInvalidRow)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	CVar v;
	ASSERT_EQ(VR_INVALIDROW, co.Get(0, 0, &v));
	ASSERT_EQ(TT_ERROR, v.type);
	ASSERT_EQ(VR_INVALIDROW, v.vresult);


	ASSERT_EQ(VR_INVALIDROW, co.Get(-1, -1, &v));
	ASSERT_EQ(TT_ERROR, v.type);
	ASSERT_EQ(VR_INVALIDROW, v.vresult);

	ASSERT_EQ(0, co.PushBackEmpty("heading"));
	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	ASSERT_EQ(VR_OK, co.Get(0, 0, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("heading"), std::string(v.sVal));

	ASSERT_EQ(VR_OK, co.Get(1, 0, &v));
	ASSERT_EQ(TT_EMPTY, v.type);


	ASSERT_EQ(VR_INVALIDROW, co.Get(2, 0, &v));
	ASSERT_EQ(TT_ERROR, v.type);
	ASSERT_EQ(VR_INVALIDROW, v.vresult);
}

TEST(TestSelectedOutput, TestInvalidCol)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	CVar v;
	ASSERT_EQ(VR_INVALIDROW, co.Get(0, 0, &v));
	ASSERT_EQ(TT_ERROR, v.type);
	ASSERT_EQ(VR_INVALIDROW, v.vresult);


	ASSERT_EQ(VR_INVALIDROW, co.Get(-1, -1, &v));
	ASSERT_EQ(TT_ERROR, v.type);
	ASSERT_EQ(VR_INVALIDROW, v.vresult);

	ASSERT_EQ(0, co.PushBackEmpty("heading"));
	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	ASSERT_EQ(VR_OK, co.Get(0, 0, &v));
	ASSERT_EQ(TT_STRING, v.type);
	ASSERT_EQ(std::string("heading"), std::string(v.sVal));

	ASSERT_EQ(VR_OK, co.Get(1, 0, &v));
	ASSERT_EQ(TT_EMPTY, v.type);


	ASSERT_EQ(VR_INVALIDCOL, co.Get(0, 1, &v));
	ASSERT_EQ(TT_ERROR, v.type);
	ASSERT_EQ(VR_INVALIDCOL, v.vresult);

	ASSERT_EQ(VR_INVALIDCOL, co.Get(0, -1, &v));
	ASSERT_EQ(TT_ERROR, v.type);
	ASSERT_EQ(VR_INVALIDCOL, v.vresult);
}

TEST(TestSelectedOutput, TestGet)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	ASSERT_EQ(0, co.PushBackEmpty("heading"));
	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());


	CVar v0 = co.Get(0, 0);
	ASSERT_EQ(TT_STRING, v0.type);
	ASSERT_EQ(std::string("heading"), std::string(v0.sVal));

	CVar v1 = co.Get(1, 0);
	ASSERT_EQ(TT_EMPTY, v1.type);
}

TEST(TestSelectedOutput, TestLongHeadings)
{
	CSelectedOutput co;
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	co.Clear();
	ASSERT_EQ((size_t)0, co.GetColCount());
	ASSERT_EQ((size_t)0, co.GetRowCount());

	ASSERT_EQ(0, co.PushBackEmpty("heading890123456789012345678901234567890123456789"));
	ASSERT_EQ(0, co.EndRow());
	ASSERT_EQ((size_t)1, co.GetColCount());
	ASSERT_EQ((size_t)2, co.GetRowCount());

	CVar v0 = co.Get(0, 0);
	ASSERT_EQ(TT_STRING, v0.type);
	ASSERT_EQ(std::string("heading890123456789012345678901234567890123456789"), std::string(v0.sVal));

	CVar v1 = co.Get(1, 0);
	ASSERT_EQ(TT_EMPTY, v1.type);
}
