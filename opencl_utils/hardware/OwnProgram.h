#pragma once
#include <OpenCl.h>

class OwnProgram : public cl::Program
{
public:
	/*static std::unique_ptr<cl::Program> load(std::shared_ptr<OpenCL> openCl, std::vector<std::string> sourceFiles);
	static std::unique_ptr<cl::Program> load(std::shared_ptr<OpenCL> openCl, std::vector<int> resourceIds);
	static std::unique_ptr<cl::Program> load(std::shared_ptr<OpenCL> openCl, std::string binaryFile);*/
	
	void compile(std::string kernelName, std::vector<std::string> defines);
	void save(std::string binaryFilename);

protected:
	std::shared_ptr<OpenCL> openCl;
};