#pragma once

#ifdef DISABLE_ASSERT_TDS_CONDTROL
#define TDS_CONTROL_PRECONDITION(cond, message)
#else
#include <cassert>
#define TDS_CONTROL_PRECONDITION(cond, message) assert(cond)
#endif
