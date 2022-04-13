#pragma once
#include <string>
#include <map>
#include <sstream>
#include <cstdint>
#include <memory>

#define GET_STATS

#ifdef GET_STATS
	#define STATS_WRITE(key,value) writeStats(key, value)
	#define STATS_ADD(key,value) addStats(key, value)
#else
	#define STATS_WRITE(key,value) 
	#define STATS_ADD(key,value) 
#endif


class IStat {
public:
	virtual ~IStat() {}
	virtual std::string toString() const = 0;
	virtual void add(const IStat& other) = 0;

	virtual std::shared_ptr<IStat> clone() = 0;
};

template <class T>
class Stat : public IStat {
public:
	typedef T value_type;
	
	Stat(T value) : value(value) {}
	
	std::shared_ptr<IStat> clone() {
		auto copy = std::make_shared<Stat<T>>(this->value);
		return copy;
	}
	
	virtual void add(T other) {
		value += other;
	}

	virtual void add(const IStat& other) {
		auto casted = dynamic_cast<const Stat<T>&>(other);
		value += casted.value;
	}

	virtual std::string toString() const { return std::to_string(value); }
protected:
	T value;
}; 

inline std::ostream& operator<<(std::ostream& os, IStat& stats) {
	os << stats.toString();
	return os;
}

class Statistics
{
public:
	
	virtual ~Statistics() {}

	template<class T>
	void put(const std::string& key, const T& value) {
		statistics[key] = std::make_shared<Stat<T>>(value);
	}

	template <class T>
	void add(const std::string& key, const T& other) {
		auto casted = std::dynamic_pointer_cast<Stat<T>>(statistics[key]);
		casted->add(other);
	}

	void clear() { statistics.clear(); }

	std::string toString() const
	{
		std::ostringstream oss;

		for (auto s = statistics.begin(); s != statistics.end(); s++) {
			oss << s->first << "=" << s->second->toString() << std::endl;
		}

		return oss.str();
	}


protected:
	std::map<std::string, std::shared_ptr<IStat>> statistics;
};