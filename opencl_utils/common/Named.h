#pragma once
#include <string>

/// <summary>
/// Base class for all object having names.
/// </summary>
class Named
{
public:
	/// <summary>
	/// Gets object name.
	/// </summary>
	virtual std::string getName()
	{
		return "";
	}
};
