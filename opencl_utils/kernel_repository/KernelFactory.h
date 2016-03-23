#pragma once
#include <memory>
#include <vector>
#include "../opencl_utils/hardware/Kernel.h"
#include "../opencl_utils/hardware/OpenCl.h"
#include "KernelRepository.h"

class KernelFactory
{
public:
	std::string buildLog;
	
	bool getFastMath() const { return this->fastMath; } 
	void setFastMath(bool fastMath) { this->fastMath = fastMath; } 

	/// Singleton access method.
	static KernelFactory& instance(std::shared_ptr<clex::OpenCL> cl) 
	{
		static std::map<std::shared_ptr<clex::OpenCL>, std::shared_ptr<KernelFactory>> devices2factories;
		
		if (devices2factories.find(cl) == devices2factories.end())
			devices2factories[cl] = std::shared_ptr<KernelFactory>(new KernelFactory(cl));

		return *devices2factories[cl];
	}
	
	/// <summary>
	/// </summary>
	/// <param name=openCl></param>
	/// <param name=filenames></param>
	/// <param name=kernelName></param>
	/// <param name=defines></param>
	/// <returns></returns>
	std::unique_ptr<clex::Kernel> create(
		std::vector<std::string> filenames, 
		std::string kernelName,
		std::vector<std::string> defines = std::vector<std::string>(0),
		int maxRegisters = 0);

	/// <summary>
	/// </summary>
	/// <param name=identifiers></param>
	/// <param name=kernelName></param>
	/// <param name=defines></param>
	/// <returns></returns>
	std::unique_ptr<clex::Kernel> create(
		std::vector<int> resourceIds, 
		std::string kernelName,
		std::vector<std::string> defines = std::vector<std::string>(0),
		int maxRegisters = 0);

protected:
	/// Maximum size of OpenCL source in bytes.
	static const int MAX_SOURCE_SIZE = 10000000;
	
	std::shared_ptr<clex::OpenCL> openCl;

	bool fastMath;

	/// Map of prototypes.
	std::map<std::string, std::unique_ptr<clex::Kernel>> prototypes;

	/// Protected constructor needed by singleton.
	KernelFactory(std::shared_ptr<clex::OpenCL> openCl) : openCl(openCl), fastMath(false) {}

	/// <summary>
	/// </summary>
	/// <param name=program></param>
	/// <param name=kernelName></param>
	/// <param name=defines></param>
	/// <returns></returns>
	void compileProgram(
		cl::Program &program,  
		std::vector<std::string> defines,
		int maxRegisters);

	/// <summary>
	/// </summary>
	/// <param name=binaryFilename></param>
	/// <returns></returns>
	void saveProgram(const cl::Program& program, std::string binaryFilename);

	/// <summary>
	/// </summary> 
	/// <param name=sourceFiles></param>
	/// <returns></returns>
	std::unique_ptr<cl::Program> loadProgram(
		std::vector<std::string> sourceFiles, 
		std::vector<std::string> defines,
		int maxRegisters);

	/// <summary>
	/// </summary>
	/// <param name=identifiers></param>
	/// <returns></returns>
	std::unique_ptr<cl::Program> loadProgram(
		std::vector<int> resourceIds, 
		std::vector<std::string> defines,
		int maxRegisters);

	/// <summary>
	/// </summary>
	/// <param name=binaryFile></param>
	/// <returns></returns>
	std::unique_ptr<cl::Program> loadProgram(
		std::string binaryFile,
		std::vector<std::string> defines,
		int maxRegisters);
};