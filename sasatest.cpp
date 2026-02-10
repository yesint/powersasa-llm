#include <cmath>
#include <fstream>
#include <numeric>
#include <print>
#include <sstream>
#include <string>
#include <vector>

#include <power_sasa.h>
#include "vec3.h"

int main()
{
	using Scalar = float;
	using Coord = Vec3<Scalar>;
	std::println("Precision: float");

	std::vector<Coord> coords; // Vector of coordinates
	std::vector<Scalar> weights; // Vector of radii including probe

	// Read coordinates from protein_coords.txt
	const std::string coordFile = "protein_coords.txt";
	std::ifstream infile(coordFile);
	
	if (!infile.is_open()) {
		std::cerr << "Error: Could not open " << coordFile << std::endl;
		return 1;
	}

	std::string line;
	int atomCount = 0;
	
	while (std::getline(infile, line)) {
		// Skip comments and empty lines
		if (line.empty() || line[0] == '#') {
			continue;
		}

		std::istringstream iss(line);
		Scalar x, y, z, radius;
		std::string element, atomName;

		if (iss >> x >> y >> z >> radius >> element >> atomName) {
			coords.push_back(Coord(x, y, z));
			weights.push_back(radius + 0.14f);
			atomCount++;
		}
	}

	std::println("Loaded {} atoms from {}", atomCount, coordFile);

	if (coords.empty()) {
		std::cerr << "Error: No atoms loaded from " << coordFile << std::endl;
		return 1;
	}

	std::println("Computing SASA...");
	POWERSASA::PowerSasa<Scalar,Coord> ps = POWERSASA::PowerSasa<Scalar,Coord>(coords, weights, 1, 0, 1, 0);
	ps.calc_sasa_all();

	// Compute total SASA and volume
	const std::vector<Scalar>& sasa = ps.getSasa();
	const std::vector<Scalar>& vol = ps.getVol();

	const Scalar totalSasa = std::accumulate(sasa.begin(), sasa.end(), Scalar(0));
	const Scalar totalVol = std::accumulate(vol.begin(), vol.end(), Scalar(0));

	std::println("Total SASA: {} nm^2", totalSasa);
	std::println("Total Volume: {} nm^3", totalVol);

	// Regression test: Check against Golden Values
	constexpr Scalar goldenSasa = 144.812f;
	constexpr Scalar goldenVol = 56.747f;
	constexpr Scalar epsilon = 1e-3f;

	if (std::abs(totalSasa - goldenSasa) > epsilon) {
		std::cerr << "REGRESSION FAILURE: SASA value changed! Expected: " << goldenSasa << ", Got: " << totalSasa << std::endl;
		return 1;
	}

	if (std::abs(totalVol - goldenVol) > epsilon) {
		std::cerr << "REGRESSION FAILURE: Volume value changed! Expected: " << goldenVol << ", Got: " << totalVol << std::endl;
		return 1;
	}

	std::println("REGRESSION PASS: SASA and Volume match golden values.");

	return 0;
}
