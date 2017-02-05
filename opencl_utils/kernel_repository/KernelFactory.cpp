#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "KernelFactory.h"
#include "KernelRepository.h"

#include "../opencl_utils/common/dbgnew.h"
#include "../opencl_utils/common/Log.h"

using namespace std;

/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<clex::Kernel> KernelFactory::create(
	std::vector<std::string> filenames,
	std::string kernelName,
	std::vector<std::string> defines,
	int maxRegisters)
{
	return nullptr;
}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<clex::Kernel> KernelFactory::create(
	std::vector<int> identifiers,
	std::string kernelName,
	std::vector<std::string> defines,
	int maxRegisters)
{
	// generate a key
	
	std::string key = openCl->mainDevice->info->deviceName + kernelName + std::to_string(this->fastMath);
	key = accumulate( defines.begin(), defines.end(), key );
	
	// if prototype with given key does not exists...
//	if (prototypes.find(key) == prototypes.end()) {

		auto numericKey =  std::hash<std::string>()(key);
		std::string binaryFilename = std::to_string(numericKey) + ".bin";
		std::unique_ptr<cl::Program> program;
		bool binaryExists = false;// boost::filesystem::exists(binaryFilename);
		// fixme:
#ifdef _DEBUG
		binaryExists = false;
#endif
		//binaryExists = false;

		// check if there is a file with compiled binaries
		if (openCl->mainDevice->info->vendor == clex::AMD && binaryExists) {
			LOG_DEBUG << "Loading kernel: " << kernelName << "...";
			program = loadProgram(binaryFilename, defines, maxRegisters);
			LOG_DEBUG << "OK";
		} else {
			LOG_DEBUG << "Compiling kernel: " << kernelName << "...";

			program = loadProgram(identifiers, defines, maxRegisters);
			if (openCl->mainDevice->info->vendor == clex::AMD) {
				saveProgram(*program, binaryFilename);
			}
#ifdef DEBUG
			ofstream reportFile("kernel-report-" + kernelName + ".txt");
			reportFile << buildLog;
#endif
			LOG_DEBUG << "OK";
		}

		auto kernel = std::unique_ptr<clex::Kernel>(new clex::Kernel(openCl, *program, kernelName.c_str()));
	//	prototypes[key] = std::unique_ptr<clex::Kernel>(new clex::Kernel(openCl, *program, kernelName.c_str()));
		LOG_DEBUG << " (fastMath = " << fastMath << ", " << kernel->privateMemSize / 4 << " registers)" << endl;
//	}

	// clone a prototype and return it
	//auto copy = std::unique_ptr<clex::Kernel>(new clex::Kernel(*prototypes[key]));
	return kernel;
}

/// <summary>
/// See declaration for all the details.
/// </summary>
void KernelFactory::saveProgram(const cl::Program& program, std::string binaryFilename)
{
	cl_int code;
	const std::vector< ::size_t> sizes = program.getInfo<CL_PROGRAM_BINARY_SIZES>(&code);
	clCall(code);
	const auto binaries = program.getInfo<CL_PROGRAM_BINARIES>(&code);
	clCall(code);
	ofstream file(binaryFilename, ios::binary);

	for (int i = 0; i < binaries.size(); i++) {
		::size_t binarySize = sizes[i];
		file.write((const char*)&binarySize, sizeof(::size_t));
		file.write(binaries[i], binarySize);
		delete binaries[i];
	}

}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<cl::Program> KernelFactory::loadProgram(
	std::vector<std::string> sourceFilenames,
	std::vector<std::string> defines,
	int maxRegisters)
{
	char *buffer = new char[MAX_SOURCE_SIZE];
	char *position = buffer;

	// join all files together
	for (int i = 0; i < sourceFilenames.size(); i++) {
		ifstream file;
		file.open(sourceFilenames[i].c_str());

		if (file.is_open()) {
			file.read(position, MAX_SOURCE_SIZE);
			position += file.gcount();
			*position = '\n';
			position++;
			file.close();
		}
		else
		{
			throw "Unable to open file!";
		}
	}

	*position = 0;

	int code;
	cl::Program::Sources sources;
	sources.push_back(std::make_pair(buffer, strlen(buffer) + 1));
	auto program = std::unique_ptr<cl::Program>(new cl::Program(*openCl->context, sources, &code));
	clCall(code);

	compileProgram(*program, defines, maxRegisters);

	delete [] buffer;
	return program;
}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<cl::Program> KernelFactory::loadProgram(
	std::vector<int> resourceIds,
	std::vector<std::string> defines,
	int maxRegisters)
{
	cl::Program::Sources sources;

	for (int i = 0; i < resourceIds.size(); i++)
	{
		const char* data = KernelRepository::kernels[resourceIds[i]];
		::size_t size = strlen(data);
		sources.push_back(std::make_pair(data, size));
	}

	int code;
	auto program = std::unique_ptr<cl::Program>(new cl::Program(*openCl->context, sources, &code));
	clCall(code);

	compileProgram(*program, defines, maxRegisters);

	return program;
}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<cl::Program> KernelFactory::loadProgram(
	std::string binaryFilename,
	std::vector<std::string> defines,
	int maxRegisters)
{
	char* buffer = new char[MAX_SOURCE_SIZE];
	char* position = buffer;
	ifstream file(binaryFilename, ios::binary);
	cl::Program::Binaries binaries;

	while (!file.eof())
	{
		::size_t binarySize;
		file.read((char*)&binarySize, sizeof(::size_t));
		file.read(position, binarySize);
		binaries.push_back(std::make_pair(position, binarySize));
		position += binarySize;
		file.peek();
	}

	int code;
	std::vector<cl_int> stats(binaries.size());
	std::vector<cl::Device> devices(1, openCl->mainDevice->device);

	auto program = std::unique_ptr<cl::Program>(new cl::Program(*openCl->context, devices, binaries, &stats, &code));
	clCall(code);

	compileProgram(*program, defines, maxRegisters);

	delete [] buffer;
	return program;
}

/// <summary>
/// See declaration for all the details.
/// </summary>
void KernelFactory::compileProgram(
	cl::Program &program,
	std::vector<std::string> defines,
	int maxRegisters)
{
	std::vector<std::string> options;

	if (this->fastMath) {
		options.push_back("-cl-mad-enable");
		options.push_back("-cl-fast-relaxed-math");
	//	options.push_back("-cl-unsafe-math-optimizations");
	//	options.push_back("-cl-denorms-are-zero");
	//	options.push_back("-cl-finite-math-only");
	}

	//options.push_back("-save-temps=d:/cl/");
	//options.push_back("-cl-opt-disable");
	//options.push_back("-cl-std=CL1.1");
	//options.push_back("-cl-strict-aliasing");
	//options.push_back("-cl-no-signed-zeros");
	//options.push_back("-cl-single-precision-constant");

	if (openCl->mainDevice->info->vendor == clex::NVidia) {
		options.push_back("-cl-nv-verbose");

		if (maxRegisters > 0) {
			std::string option = "-cl-nv-maxrregcount=" + std::to_string(maxRegisters);
			options.push_back(option);
		}
	}

	std::ostringstream oss;
    copy(options.begin(), options.end(), ostream_iterator<std::string>(oss, " "));
	oss << " -D ";
	std::copy(defines.begin(), defines.end(), ostream_iterator<std::string>(oss, " -D "));
	oss << "KERNEL";

	std::string commandLine = oss.str();
	std::vector<cl::Device> devices(1, openCl->mainDevice->device);

	int code = program.build(devices, commandLine.c_str());
	int status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(openCl->mainDevice->device);
	this->buildLog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(openCl->mainDevice->device);

	if ((openCl->mainDevice->info->vendor == clex::NVidia && buildLog.length() > 0) || code != CL_SUCCESS) {
		cout << endl << buildLog << endl;
	}

	clCall(code);
}
