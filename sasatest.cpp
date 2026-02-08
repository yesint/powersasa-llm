#include <stdio.h>
#include <fstream>
#include <sstream>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

#include <power_sasa.h>
#include <testincludes.h>

// g++ -o sasatest -I. sasatest.cpp

int main(int argc, char* argv[])
{
	typedef float Scalar;
	typedef Eigen::Vector3f Coord;
	std::cout << "Precision: float" << std::endl;

	std::vector<Coord> coords; //Vector of coordinates
	std::vector<Scalar> weights; //Vector of weights (VdW radii)

	// Read coordinates from protein_coords.txt
	std::string coordFile = "protein_coords.txt";
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
			weights.push_back(radius+0.14);
			atomCount++;
		}
	}

	infile.close();

	std::cout << "Loaded " << atomCount << " atoms from " << coordFile << std::endl;

	if (coords.empty()) {
		std::cerr << "Error: No atoms loaded from " << coordFile << std::endl;
		return 1;
	}

	std::cout << "Computing SASA..." << std::endl;
	POWERSASA::PowerSasa<Scalar,Coord> ps = POWERSASA::PowerSasa<Scalar,Coord>(coords, weights, 1, 0, 1, 0);
	ps.calc_sasa_all();

	// Compute total SASA and volume
	const std::vector<Scalar>& sasa = ps.getSasa();
	const std::vector<Scalar>& vol = ps.getVol();

	Scalar totalSasa = 0.0;
	Scalar totalVol = 0.0;

	for (const auto& s : sasa) {
		totalSasa += s;
	}

	for (const auto& v : vol) {
		totalVol += v;
	}

	std::cout << "Total SASA: " << totalSasa << " nm^2" << std::endl;
	std::cout << "Total Volume: " << totalVol << " nm^3" << std::endl;

	// Regression test: Check against Golden Values
	const Scalar goldenSasa = 144.812f;
	const Scalar goldenVol = 56.747f;
	const Scalar epsilon = 1e-3f; // Tolerance

	if (std::abs(totalSasa - goldenSasa) > epsilon) {
		std::cerr << "REGRESSION FAILURE: SASA value changed! Expected: " << goldenSasa << ", Got: " << totalSasa << std::endl;
		return 1;
	}

	if (std::abs(totalVol - goldenVol) > epsilon) {
		std::cerr << "REGRESSION FAILURE: Volume value changed! Expected: " << goldenVol << ", Got: " << totalVol << std::endl;
		return 1;
	}

	std::cout << "REGRESSION PASS: SASA and Volume match golden values." << std::endl;

	return 0;
}
