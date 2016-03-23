#include "Log.h"

const int Log::LEVEL_DEBUG = 0;
const int Log::LEVEL_NORMAL = 1;

Log::Log() 
{
	enabled = false;
	out = &std::cerr;
}

Log& Log::operator<< (std::ostream& (*pf)(std::ostream&))
{
	if (enabled) { *this->out << pf;  }
	return *this;
}

Log& Log::operator<< (std::ios& (*pf)(std::ios&))
{
	if (enabled) { *this->out << pf;  }
	return *this;
}

Log& Log::operator<< (std::ios_base& (*pf)(std::ios_base&))
{
	if (enabled) { *this->out << pf;  }
	return *this;
}

