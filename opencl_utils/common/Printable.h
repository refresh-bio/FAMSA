#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <ostream>

/// <summary>
/// Base class for all printable objects.
/// </summary>
class Printable
{
public:
	/// <summary>
	/// Gets string representation of an object.
	/// </summary>
	/// <return>String representation.</return>
	virtual std::string toString()
	{ 
		return "";
	}
};

/// <summary>
/// Directs string representation of object to an output stream.
/// </summary>
/// <param name ="s">Output stream.</param>
/// <param name ="printable">Printable object.</param>
/// <return>Reference to an output stream.</return>
inline std::ostream &operator<<(std::ostream &s, Printable &printable)
{
	return s << printable.toString();
}