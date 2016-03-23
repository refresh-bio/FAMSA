#pragma once
#include <CL/cl.h>
#include <string>
#include <iostream>

namespace  clex
{
	class EventHandler
	{
	public:

		EventHandler(cl::Event& event, std::string message) : event(event), message(message) 
		{
			int code = event.setCallback(CL_COMPLETE, handle, this);
		}

		static void CL_CALLBACK handle(cl_event id, cl_int status, void * data) 
		{
			auto ths = static_cast<EventHandler*>(data);
			std::cout << ths->message << std::endl;
		}

	protected:
		cl::Event& event;
		const std::string message;
	};
}