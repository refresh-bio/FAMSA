//#include "OwnProgram.h"

/*
OwnProgram::OwnProgram(std::shared_ptr<OpenCL> openCl, std::vector<std::string> sourceFilenames)
	: openCl(openCl)
{
	char *buffer = new char[MAX_SOURCE_SIZE];
	char *position = buffer;

	// join all files together
	for (int i = 0; i < sourceFilenames.size(); i++)
	{
		ifstream file;
		file.open(sourceFilenames[i].c_str());

		if (file.is_open())
		{
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
	Program::Sources sources;
	sources.push_back(std::make_pair(buffer, strlen(buffer) + 1));
	auto program = std::unique_ptr<cl::Program>(new cl::Program(openCl->context, sources, &code));
	clCall(code);

	delete [] buffer;
	return program;
}

OwnProgram::OwnProgram(std::shared_ptr<OpenCL> openCl, std::vector<int> resourceIds)
	: openCl(openCl)
{
	Program::Sources sources;

	for (int i = 0; i < resourceIds.size(); i++)
	{
		::size_t size;
		const char* data = ResourceManager::load(resourceIds[i], "KERNEL_FILE", size);
		sources.push_back(std::make_pair(data, size));
	}

	int code;
	auto program = std::unique_ptr<cl::Program>(new cl::Program(openCl->context, sources, &code));
	clCall(code);

	return program;
}

OwnProgram::OwnProgram(std::shared_ptr<OpenCL> openCl, std::string binaryFilename)
	: openCl(openCl)
{
	char* buffer = new char[MAX_SOURCE_SIZE];
	char* position = buffer;
	ifstream file(binaryFilename, ios::binary);
	Program::Binaries binaries;

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
	auto program = std::unique_ptr<cl::Program>(
		new cl::Program(openCl->context, std::vector<cl::Device>(1, openCl->device), binaries, nullptr, &code));
	clCall(code);

	delete [] buffer;
}


void OwnProgram::compile(std::string kernelName, std::vector<std::string> defines)
{
	cerr << "Compiling kernel: " << kernelName << "...";

	std::ostringstream oss;
	std::copy(defines.begin(), defines.end(), ostream_iterator<std::string>(oss, " -D "));
	std::string options = " -D " + oss.str() + "KERNEL";

	int code = this->build(std::vector<cl::Device>(1, openCl->device), options.c_str());
	int status = this->getBuildInfo<CL_PROGRAM_BUILD_STATUS>(openCl->device);
	std::string buildLog = this->getBuildInfo<CL_PROGRAM_BUILD_LOG>(openCl->device);

#ifdef DEBUG
	ofstream reportFile("kernel-report-" + kernelName + ".txt");
	reportFile << buildLog;
#endif
	clCall(code);
	cerr << "OK!" << endl;
}

void OwnProgram::save(std::string binaryFilename)
{
	cl_int code;
	const std::vector<::size_t> sizes = this->getInfo<CL_PROGRAM_BINARY_SIZES>(&code);
	clCall(code);
	const auto binaries = this->getInfo<CL_PROGRAM_BINARIES>(&code);
	clCall(code);
	ofstream file(binaryFilename, ios::binary);

	for (int i = 0; i < binaries.size(); i++)
	{
		::size_t binarySize = sizes[i];
		file.write((const char*)&binarySize, sizeof(::size_t));
		file.write(binaries[i], binarySize);
	}
}
*/