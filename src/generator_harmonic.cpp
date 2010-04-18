#include <stdexcept>
#include <sstream>
#include <string>
#include "generator.h"

using namespace std;

struct HarmonicGenerator: public Generator
{
	HarmonicGenerator(): Generator()
	{
	}

	~HarmonicGenerator()
	{
	}

	void set_surface(const string &p)
	{
		if (p != "sphere")
		{
			std::stringstream str;
			str << "unknown surface type ";
			str << p;
			str << " for harmonic method!";
			throw logic_error(str.str());
		}
	}

	void set_approximation(const string &p)
	{
		if (p != "time")
		{
			std::stringstream str;
			str << "unknown approximation type ";
			str << p;
			str << " for harmonic method!";
			throw logic_error(str.str());
		}
	}
};

Generator * harmonic_generator_creator()
{
	return new HarmonicGenerator();
}
