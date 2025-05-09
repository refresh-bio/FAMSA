#ifndef _FILE_WRAPPER_H
#define _FILE_WRAPPER_H

// *** History of updates
// *** v. 1.0.1 (2024-03-11) - bug fix (wrong zlib initialization)
// *** v. 1.0.2 (2024-05-01) - bug fix (wrong reading from file)
// ***

#include <cstdint>
#include <cstring>
#include <cstdio>
#include <string>
#include <utility>
#include <memory>
#include <tuple>
#include <vector>
#include <set>
#include <array>
#include <map>
#include <cassert>
#include <algorithm>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#if defined(ARCH_X64)
#if defined(REFRESH_USE_IGZIP)
#define REFRESH_STREAM_DECOMPRESSION_ENABLE_IGZIP
#else
#define REFRESH_STREAM_DECOMPRESSION_ENABLE_ZLIB
#endif
#else
#define REFRESH_STREAM_DECOMPRESSION_ENABLE_ZLIB
#endif

#ifdef REFRESH_STREAM_DECOMPRESSION_ENABLE_IGZIP
#include <igzip_lib.h>
#endif

#ifdef REFRESH_STREAM_DECOMPRESSION_ENABLE_ZLIB
#include <zlib.h>
#endif

#ifdef REFRESH_STREAM_DECOMPRESSION_ENABLE_ZSTD
#include <zstd.h>
#endif


namespace refresh
{
	// **********************************************************************************
	// Base class for basic I/O operations, like reading from file, stdin, ...
	// **********************************************************************************
	class stream_in_base
	{
	public:
		stream_in_base() = default;
		virtual ~stream_in_base() = default;

		virtual std::pair<char*, size_t> read() = 0;
		virtual void release(char*) = 0;
		virtual std::string get_file_name() const = 0;
	};

	// **********************************************************************************
	// Base class for buffered reading
	// **********************************************************************************
	class stream_in_buffered : public stream_in_base
	{
	protected:
		size_t buffer_size;
		char* buffer;
		size_t buffer_filled;
		bool buffer_released;

	public:
		stream_in_buffered(size_t buffer_size) : 
			stream_in_base(),
			buffer_size(buffer_size)
		{
			buffer = new char[buffer_size];
			buffer_released = true;
			buffer_filled = 0;
		}

		virtual ~stream_in_buffered()
		{
			delete[] buffer;
		}
	};

	// **********************************************************************************
	// Low-level reading from stdin
	// **********************************************************************************
	class stream_in_stdin : public stream_in_buffered
	{
	public:
		stream_in_stdin(size_t buffer_size = 8 << 20) :
			stream_in_buffered(buffer_size)
		{
#ifdef _WIN32
			_setmode(_fileno(stdin), _O_BINARY);
#endif
		}

		virtual ~stream_in_stdin()
		{}

		virtual std::pair<char*, size_t> read()
		{
			if (!buffer_released)
				assert(0);

			if (!feof(stdin))
				buffer_filled = fread(buffer, 1, buffer_size, stdin);
			else
				buffer_filled = 0;

			buffer_released = false;

			return std::make_pair(buffer, buffer_filled);
		}

		virtual void release(char *ptr)
		{
			buffer_released = true;
		}

		virtual std::string get_file_name()  const
		{
			return "";
		}
	};

	// **********************************************************************************
	// Low-level reading from file
	// **********************************************************************************
	class stream_in_file : public stream_in_buffered
	{
		std::string file_name;
		size_t io_buffer_size;
		FILE* file = nullptr;
		bool test_extension = true;

	public:
		stream_in_file(const std::string& file_name, size_t io_buffer_size = 16 << 20, size_t buffer_size = 8 << 20, bool test_extension = true) :
			stream_in_buffered(buffer_size),
			file_name(file_name),
			io_buffer_size(io_buffer_size),
			test_extension(test_extension)
		{
			open(file_name);
		}

		virtual ~stream_in_file()
		{
			close();
		}

		virtual bool open(const std::string& _file_name)
		{
			file_name = _file_name;

			if (file)
				close();

			file = fopen(file_name.c_str(), "rb");

			if (!file)
				return false;

			setvbuf(file, nullptr, _IOFBF, io_buffer_size);

			buffer_released = true;

			return true;
		}

		virtual bool close()
		{
			if (file)
			{
				fclose(file);
				file = nullptr;

				return true;
			}

			return false;
		}

		virtual bool restart()
		{
			if (!file)
				return false;

#ifdef _WIN32
			_fseeki64(file, 0, SEEK_SET);
#else
			fseek(file, 0, SEEK_SET);
#endif

			buffer_filled = 0;
			buffer_released = true;

			return true;
		}

		virtual std::pair<char*, size_t> read()
		{
			if (!buffer_released)
				assert(0);

			buffer_filled = fread(buffer, 1, buffer_size, file);
			buffer_released = false;

			return std::make_pair(buffer, buffer_filled);
		}

		virtual void release(char *ptr)
		{
			buffer_released = true;
		}

		virtual std::string get_file_name()  const
		{
			return test_extension ? file_name : "";
		}

		bool is_open() const
		{
			return file != nullptr;
		}
	};

	// **********************************************************************************
	// Base class of decompression engines, like zlib, igzip, plain-text, zstd, ...
	// **********************************************************************************
	class stream_decompression_engine
	{
	protected:
		stream_in_base* stream_in;
		size_t out_buffer_size;
		char* in_buffer_data = nullptr;
		size_t in_buffer_filled = 0;
		size_t in_buffer_pos = 0;

	public:
		stream_decompression_engine(stream_in_base* stream_in, size_t out_buffer_size, char *in_buffer_data, size_t in_buffer_filled) :
			stream_in(stream_in),
			out_buffer_size(out_buffer_size),
			in_buffer_data(in_buffer_data),
			in_buffer_filled(in_buffer_filled)
		{}

		virtual ~stream_decompression_engine()
		{}

		static std::string ext(const std::string& fn, const size_t len)
		{
			if (fn.length() <= len)
				return "";

			return fn.substr(fn.length() - len, std::string::npos);
		}

		virtual int read(char* ptr, size_t &readed) = 0;
	};

// **********************************************************************************
// Decompression engine passing text files without any decompression
// **********************************************************************************
	class stream_decompression_engine_text : public stream_decompression_engine
	{

	public:
		stream_decompression_engine_text(stream_in_base* stream_in, size_t out_buffer_size, char* in_buffer_data, size_t in_buffer_filled)
			: stream_decompression_engine(stream_in, out_buffer_size, in_buffer_data, in_buffer_filled)
		{}

		virtual ~stream_decompression_engine_text()
		{
		}

		static bool knows_it(const std::string& file_name, char* data, size_t size)
		{
			return true;
		}

		virtual int read(char* ptr, size_t& readed)
		{
			readed = 0;

			while (true)
			{
				size_t to_copy = std::min<size_t>(in_buffer_filled - in_buffer_pos, out_buffer_size - readed);
				memcpy(ptr + readed, in_buffer_data + in_buffer_pos, to_copy);
				in_buffer_pos += to_copy;
				readed += to_copy;

				if (in_buffer_pos == in_buffer_filled)
				{
					stream_in->release(in_buffer_data);
					in_buffer_pos = 0;
					std::tie(in_buffer_data, in_buffer_filled) = stream_in->read();

					if (in_buffer_filled == 0)
						return readed ? 0 : -1;
				}

				if (readed == out_buffer_size)
					return 0;
			}
		}
	};

#ifdef REFRESH_STREAM_DECOMPRESSION_ENABLE_IGZIP
	// **********************************************************************************
	// Decompression engine using igzip library (for .gz files)
	// **********************************************************************************
	class stream_decompression_engine_igzip : public stream_decompression_engine
	{
		constexpr static std::array<uint8_t, 2> magic_numbers = {0x1f, 0x8b};

		struct isal_gzip_header gz_hdr;
		struct inflate_state state;

		enum class internal_state_t {none, file_start, before_header, inside};
		internal_state_t internal_state = internal_state_t::none;

		void initialize()
		{
			isal_gzip_header_init(&gz_hdr);
			isal_inflate_init(&state);
//			state.crc_flag = ISAL_GZIP_NO_HDR_VER;
			state.crc_flag = ISAL_GZIP;

			state.next_in = (uint8_t*) in_buffer_data;
			state.avail_in = (uint32_t) in_buffer_filled;

			internal_state = internal_state_t::file_start;
		}

		bool check_header()
		{
			if (state.avail_in >= magic_numbers.size())
				return std::equal(magic_numbers.begin(), magic_numbers.end(), state.next_in);

			uint32_t in_state_len = state.avail_in;

			if (!std::equal(magic_numbers.begin(), magic_numbers.begin() + in_state_len, state.next_in))
				return false;

			if (isal_inflate(&state) != ISAL_DECOMP_OK)
				return false;

			if (!load_new_part())
				return false;

			return std::equal(magic_numbers.begin() + in_state_len, magic_numbers.end(), state.next_in);
		}

		bool load_new_part() 
		{
			if(in_buffer_data)
				stream_in->release(in_buffer_data);
			std::tie(in_buffer_data, in_buffer_filled) = stream_in->read();

			state.next_in = (uint8_t*)in_buffer_data;
			state.avail_in = (uint32_t)in_buffer_filled;

			return in_buffer_filled != 0;
		}

	public:
		stream_decompression_engine_igzip(stream_in_base *stream_in, size_t out_buffer_size, char *in_buffer_data, size_t in_buffer_filled)
			: stream_decompression_engine(stream_in, out_buffer_size, in_buffer_data, in_buffer_filled)
		{
			isal_inflate_init(&state);
			isal_gzip_header_init(&gz_hdr);
		}

		virtual ~stream_decompression_engine_igzip()
		{
			if (in_buffer_data)
				stream_in->release(in_buffer_data);
		}

		static bool knows_it(const std::string& file_name, const char* data, const size_t size)
		{
			if (file_name.size() > 3 && ext(file_name, 3) != std::string(".gz"))
				return false;

			if (size < magic_numbers.size())
				return false;

			return std::equal(magic_numbers.begin(), magic_numbers.end(), (uint8_t*) data);
		}

		virtual int read(char* ptr, size_t& readed)
		{
			if (internal_state == internal_state_t::none)
				initialize();

			readed = 0;

			state.next_out = (uint8_t*) ptr;
			state.avail_out = (uint32_t) out_buffer_size;

			if (internal_state == internal_state_t::before_header || internal_state == internal_state_t::file_start)
			{
				isal_inflate_reset(&state);

				if (!check_header())
				{
					if (internal_state == internal_state_t::file_start)
						return -2;													// Error, no data in file
					else
						return -1;													// Just end-of-file, maybe some junk data are still present, but do not care about them
				}

				internal_state = internal_state_t::inside;
			}

			do
			{
				if (state.avail_in == 0)
					if (!load_new_part())
						return -1;

				int ret = isal_inflate(&state);

				if (ret != ISAL_DECOMP_OK)
					return -3;														// Error, broken gzip file

				readed = state.next_out - (uint8_t*)ptr;
				state.avail_in = (uint32_t)(in_buffer_filled - (state.next_in - (uint8_t*)in_buffer_data));
				state.avail_out = (uint32_t)(out_buffer_size - readed);

				if (state.block_state == ISAL_BLOCK_FINISH)
				{
					internal_state = internal_state_t::before_header;
					break;
				}

				if (readed == out_buffer_size)
					return 0;

			} while (true);			

			return 0;
		}
	};
#endif
	 
#ifdef REFRESH_STREAM_DECOMPRESSION_ENABLE_ZLIB
	// **********************************************************************************
	// Decompression engine using zlib library (for .gz files)
	// **********************************************************************************
	class stream_decompression_engine_zlib : public stream_decompression_engine
	{
		constexpr static std::array<uint8_t, 2> magic_numbers = { 0x1f, 0x8b };

		z_stream state;

		enum class internal_state_t { none, file_start, before_header, inside };
		internal_state_t internal_state = internal_state_t::none;

		void initialize()
		{
			state.next_in = (uint8_t*)in_buffer_data;
			state.avail_in = (uint32_t)in_buffer_filled;

			state.zalloc = NULL;
			state.zfree = NULL;
			state.opaque = NULL;

			inflateInit2(&state, 15 + 32);

			internal_state = internal_state_t::file_start;
		}

		bool check_header()
		{
			if (state.avail_in >= magic_numbers.size())
				return std::equal(magic_numbers.begin(), magic_numbers.end(), state.next_in);

			uint32_t in_state_len = state.avail_in;

			if (!std::equal(magic_numbers.begin(), magic_numbers.begin() + in_state_len, state.next_in))
				return false;

			if (inflate(&state, Z_NO_FLUSH) != Z_OK)
				return false;

			if (!load_new_part())
				return false;

			return std::equal(magic_numbers.begin() + in_state_len, magic_numbers.end(), state.next_in);
		}

		bool load_new_part()
		{
			if (in_buffer_data)
				stream_in->release(in_buffer_data);
			std::tie(in_buffer_data, in_buffer_filled) = stream_in->read();

			state.next_in = (uint8_t*)in_buffer_data;
			state.avail_in = (uint32_t)in_buffer_filled;

			return in_buffer_filled != 0;
		}

	public:
		stream_decompression_engine_zlib(stream_in_base *stream_in, size_t buffer_size, char* in_buffer_data, size_t in_buffer_filled)
			: stream_decompression_engine(stream_in, buffer_size, in_buffer_data, in_buffer_filled)
		{
		}

		virtual ~stream_decompression_engine_zlib()
		{
			if (in_buffer_data)
				stream_in->release(in_buffer_data);

			if(internal_state != internal_state_t::none)
				inflateEnd(&state);
		}

		static bool knows_it(const std::string& file_name, const char* data, const size_t size)
		{
			if (file_name.size() > 3 && ext(file_name, 3) != std::string(".gz"))
				return false;

			if (size < magic_numbers.size())
				return false;

			return std::equal(magic_numbers.begin(), magic_numbers.end(), (uint8_t*) data);
		}

		virtual int read(char* ptr, size_t& readed)
		{
			if (internal_state == internal_state_t::none)
				initialize();

			readed = 0;

			state.next_out = (uint8_t*)ptr;
			state.avail_out = out_buffer_size;

			if (internal_state == internal_state_t::before_header || internal_state == internal_state_t::file_start)
			{
				if(internal_state == internal_state_t::before_header)
					inflateReset(&state);

				if (!check_header())
				{
					if (internal_state == internal_state_t::file_start)
						return -2;													// Error, no data in file
					else
						return -1;													// Just end-of-file, maybe some junk data are still present, but do not care about them
				}

				internal_state = internal_state_t::inside;
			}

			do
			{
				if (state.avail_in == 0)
					if (!load_new_part())
						return -1;

				int ret = inflate(&state, Z_NO_FLUSH);

				if (ret == Z_DATA_ERROR)
					return -3;														// Error, broken gzip file

				readed = state.next_out - (uint8_t*)ptr;
				state.avail_in = (uint32_t)in_buffer_filled - (state.next_in - (uint8_t*)in_buffer_data);
				state.avail_out = (uint32_t)(out_buffer_size - readed);

				if (ret == Z_STREAM_END)
				{
					internal_state = internal_state_t::before_header;
					break;
				}

				if (readed == out_buffer_size)
					return 0;

			} while (true);

			return 0;
		}
	};
#endif
	 
#ifdef REFRESH_STREAM_DECOMPRESSION_ENABLE_ZSTD
	// **********************************************************************************
	// Decompression engine using zstd library (for .zst files)
	// **********************************************************************************
	class stream_decompression_engine_zst : public stream_decompression_engine
	{
		constexpr static std::array<uint8_t, 4> magic_numbers = { 0x28, 0xb5, 0x2f, 0xfd };

		ZSTD_DStream *state;

		ZSTD_inBuffer in_buffer;
		ZSTD_outBuffer out_buffer;

		enum class internal_state_t { none, file_start, before_header, inside };
		internal_state_t internal_state = internal_state_t::none;

		void initialize()
		{
			ZSTD_initDStream(state);

			in_buffer.src = in_buffer_data;
			in_buffer.size = in_buffer_filled;
			in_buffer.pos = 0;

			internal_state = internal_state_t::file_start;
		}

		bool check_header()
		{
			if (in_buffer.size - in_buffer.pos >= magic_numbers.size())
				return std::equal(magic_numbers.begin(), magic_numbers.end(), (uint8_t*) in_buffer.src + in_buffer.pos);

			uint32_t in_state_len = in_buffer.size - in_buffer.pos;

			if (!std::equal(magic_numbers.begin(), magic_numbers.begin() + in_state_len, (uint8_t*)in_buffer.src + in_buffer.pos))
				return false;

			if(ZSTD_decompressStream(state, &out_buffer, &in_buffer) < 0)
				return false;

			if (!load_new_part())
				return false;

			return std::equal(magic_numbers.begin() + in_state_len, magic_numbers.end(), (uint8_t*) in_buffer.src);
		}

		bool load_new_part()
		{
			if (in_buffer_data)
				stream_in->release(in_buffer_data);
			std::tie(in_buffer_data, in_buffer_filled) = stream_in->read();

			in_buffer.src = in_buffer_data;
			in_buffer.size = in_buffer_filled;
			in_buffer.pos = 0;

			return in_buffer_filled != 0;
		}

	public:
		stream_decompression_engine_zst(stream_in_base *stream_in, size_t buffer_size, char* in_buffer_data, size_t in_buffer_filled)
			: stream_decompression_engine(stream_in, buffer_size, in_buffer_data, in_buffer_filled)
		{
			state = ZSTD_createDStream();
		}

		virtual ~stream_decompression_engine_zst()
		{
			ZSTD_freeDStream(state);
		}

		static bool knows_it(const std::string& file_name, const char* data, const size_t size)
		{
			if (file_name.size() > 4 && ext(file_name, 4) != std::string(".zst"))
				return false;

			if (size < magic_numbers.size())
				return false;

			return std::equal(magic_numbers.begin(), magic_numbers.end(), (uint8_t*) data);
		}

		virtual int read(char* ptr, size_t& readed)
		{
			if (internal_state == internal_state_t::none)
				initialize();

			readed = 0;

			out_buffer.dst = ptr;
			out_buffer.size = out_buffer_size;
			out_buffer.pos = 0;

			if (internal_state == internal_state_t::before_header || internal_state == internal_state_t::file_start)
			{
				if (internal_state == internal_state_t::before_header)
					ZSTD_initDStream(state);

				if (!check_header())
				{
					if (internal_state == internal_state_t::file_start)
						return -2;													// Error, no data in file
					else
						return -1;													// Just end-of-file, maybe some junk data are still present, but do not care about them
				}

				internal_state = internal_state_t::inside;
			}

			do
			{
				if (in_buffer.pos == in_buffer.size)
					if (!load_new_part())
						return -1;

				int ret = ZSTD_decompressStream(state, &out_buffer, &in_buffer);

				if (ret < 0)
					return -3;														// Error, broken gzip file

				readed = out_buffer.pos;

				if (ret == 0)
				{
					internal_state = internal_state_t::before_header;
					break;
				}

				if (readed == out_buffer_size)
					return 0;

			} while (true);

			return 0;
		}
	};
#endif
	 
#if defined(REFRESH_STREAM_DECOMPRESSION_ENABLE_IGZIP) && defined (REFRESH_STREAM_DECOMPRESSION_ENABLE_ZLIB)
	using stream_decompression_engine_gz = stream_decompression_engine_igzip;
#else
#if defined(REFRESH_STREAM_DECOMPRESSION_ENABLE_IGZIP)
	using stream_decompression_engine_gz = stream_decompression_engine_igzip;
#else
#if defined(REFRESH_STREAM_DECOMPRESSION_ENABLE_ZLIB)
	using stream_decompression_engine_gz = stream_decompression_engine_zlib;
#endif
#endif
#endif

	// **********************************************************************************
	// Main class for decompression of stream data (file, stdin, ...) in some compressed format (.gz, .zstd, ...)
	// **********************************************************************************
	class stream_decompression
	{
	public:
		enum class format_t { unknown, text, gzip, zstd };

	private:
		stream_decompression_engine* engine = nullptr;
		format_t format = format_t::unknown;
		size_t engine_part_size;

		char* buffer;
		size_t size;
		size_t filled;
		size_t pos;
		bool eof_marker = true;

		bool determine_format(stream_in_base* stream_in)
		{
			format = format_t::unknown;
			if (engine)
			{
				delete engine;
				engine = nullptr;
			}

			if (!stream_in)
				return false;

			char* ptr;
			std::tie(ptr, filled) = stream_in->read();

#if defined(REFRESH_STREAM_DECOMPRESSION_ENABLE_IGZIP) || defined(REFRESH_STREAM_DECOMPRESSION_ENABLE_ZLIB)
			if (stream_decompression_engine_gz::knows_it(stream_in->get_file_name(), ptr, filled))
			{
				format = format_t::gzip;
				engine = new stream_decompression_engine_gz(stream_in, engine_part_size, ptr, filled);
			}
			else 
#endif
#ifdef REFRESH_STREAM_DECOMPRESSION_ENABLE_ZSTD
			if (stream_decompression_engine_zst::knows_it(stream_in->get_file_name(), ptr, filled))
			{
				format = format_t::zstd;
				engine = new stream_decompression_engine_zst(stream_in, engine_part_size, ptr, filled);
			}
			else 
#endif
			if (stream_decompression_engine_text::knows_it(stream_in->get_file_name(), ptr, filled))
			{
				format = format_t::text;
				engine = new stream_decompression_engine_text(stream_in, engine_part_size, ptr, filled);
			}

			return format != format_t::unknown;
		}

		int fill_buffer()
		{
			int ret = engine->read(buffer, filled);
			pos = 0;

			return ret;
		}

	public:
		stream_decompression(stream_in_base *stream_in, size_t engine_part_size = 16 << 20) :
			engine_part_size(engine_part_size)
		{
			determine_format(stream_in);
			eof_marker = false;
			filled = 0;

			buffer = new char[engine_part_size];
			pos = 0;
		}

		~stream_decompression()
		{
			if (engine)
				delete engine;

			delete[] buffer;
		}

		bool restart(stream_in_base* stream_in)
		{
			auto ret = determine_format(stream_in);

			filled = 0;
			pos = 0;

			eof_marker = false;

			return ret;
		}

		bool release()
		{
			if (!engine)
				return false;

			delete engine;
			engine = nullptr;

			return true;
		}

		int getc()
		{
			if (pos < filled)
				return buffer[pos++];

			auto ret = fill_buffer();

			if (ret < 0)
			{
				eof_marker = true;
				return ret;
			}

			if (filled == 0)
				return -1;

			return buffer[pos++];
		}

		int read(char* dest, size_t requested, size_t &readed)
		{
			char* p = dest;
			readed = 0;
			int ret = 0;

			while (readed < requested)
			{
				size_t to_copy = std::min<size_t>(requested - readed, filled - pos);
				memcpy(p, buffer + pos, to_copy);
				readed += to_copy;
				pos += to_copy;
				p += to_copy;

				if (pos == filled)
				{
					ret = fill_buffer();

					if(ret < 0)
					{
						eof_marker = true;
						return ret;
					}
				}
			}

			return 0;
		}

		int getline(std::string& str)
		{
			str.clear();

			int ret = 0;

			while (true)
			{
				auto q = std::find(buffer + pos, buffer + filled, (char) 0x0a);

				str.append(buffer + pos, q);

				pos = q - buffer;

				if (pos != filled)
				{
					++pos;
					break;
				}

				ret = fill_buffer();

				if (ret < 0)
				{
					eof_marker = true;
					break;
				}
			}

			if (!str.empty() && str.back() == 0x0d)
				str.pop_back();

			return ret;
		}

		bool eof() const
		{
			return eof_marker;
		}

		format_t get_format() const
		{
			return format;
		}
	};
}

#endif
