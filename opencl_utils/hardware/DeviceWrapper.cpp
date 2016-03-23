#include "DeviceWrapper.h"


clex::DeviceWrapper::DeviceWrapper(const cl::Device& ref, const cl::Context& context,  bool useKernelProfiling)
{
	this->device = ref;
	this->info = std::make_shared<clex::DeviceInfo>(ref);

	int props = 0;
	if (useKernelProfiling) {
		props |= CL_QUEUE_PROFILING_ENABLE;
	}

	// create five command queues
	int code;
	for (int i = 0; i < 10; ++i) {
		queues.push_back(std::make_shared<cl::CommandQueue>(context, device, props, &code));
	//	clCall(code);
	}

	mainQueue = queues[queues.size() - 1];
}
