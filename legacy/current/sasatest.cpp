#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <print>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <power_sasa.h>
#include "vec3.h"

namespace {
namespace fs = std::filesystem;

template <class Scalar, class Coord>
bool load_case_file(const std::string& path, std::vector<Coord>& coords, std::vector<Scalar>& weights)
{
	coords.clear();
	weights.clear();

	std::ifstream infile(path);
	if (!infile.is_open()) {
		std::cerr << "Error: Could not open " << path << '\n';
		return false;
	}

	std::string line;
	while (std::getline(infile, line)) {
		if (line.empty() || line[0] == '#') {
			continue;
		}
		std::istringstream iss(line);
		Scalar x, y, z, radius;
		if (!(iss >> x >> y >> z >> radius)) {
			continue;
		}
		coords.push_back(Coord(x, y, z));
		weights.push_back(radius + static_cast<Scalar>(0.14));
	}

	return !coords.empty();
}

struct GoldenCase {
	std::string filePath;
	float expectedSasa = 0.0f;
	float expectedVol = 0.0f;
};

bool load_golden_manifest(const fs::path& manifestPath, std::vector<GoldenCase>& cases)
{
	cases.clear();
	std::ifstream infile(manifestPath);
	if (!infile.is_open()) {
		std::cerr << "Error: Could not open manifest " << manifestPath << '\n';
		return false;
	}

	const fs::path manifestDir = manifestPath.parent_path();
	const fs::path repoRoot = manifestDir.parent_path().parent_path();
	std::string line;
	while (std::getline(infile, line)) {
		if (line.empty() || line[0] == '#') {
			continue;
		}
		std::istringstream iss(line);
		GoldenCase c;
		if (!(iss >> c.filePath >> c.expectedSasa >> c.expectedVol)) {
			std::cerr << "Error: Malformed manifest entry: " << line << '\n';
			return false;
		}
		const fs::path listed = c.filePath;
		if (!listed.is_absolute()) {
			if (fs::exists(listed)) {
				c.filePath = listed.lexically_normal().string();
			} else if (fs::exists(repoRoot / listed)) {
				c.filePath = (repoRoot / listed).lexically_normal().string();
			} else if (fs::exists(manifestDir / listed)) {
				c.filePath = (manifestDir / listed).lexically_normal().string();
			} else {
				c.filePath = (manifestDir / listed).lexically_normal().string();
			}
		}
		cases.push_back(c);
	}

	return !cases.empty();
}

fs::path find_manifest()
{
	constexpr std::string_view kRel = "tests/testdata/sasa_cases/golden_912340c.txt";
	const fs::path candidates[] = {
	    fs::path(kRel),
	    fs::path("..") / kRel,
	    fs::path("../..") / kRel,
	    fs::path("../../..") / kRel,
	};
	for (const fs::path& candidate : candidates) {
		if (fs::exists(candidate)) {
			return candidate;
		}
	}
	return fs::path(kRel);
}

}  // namespace

int main()
{
	using Scalar = float;
	using Coord = Vec3<Scalar>;
	constexpr Scalar epsilon = 1e-3f;
	const fs::path goldenManifest = find_manifest();
	std::println("Precision: {}", std::string_view{"float"});
	std::println("Manifest: {}", goldenManifest.string());

	std::vector<GoldenCase> cases;
	if (!load_golden_manifest(goldenManifest, cases)) {
		return 1;
	}

	bool allPassed = true;
	std::vector<Coord> coords;
	std::vector<Scalar> weights;

	for (const GoldenCase& c : cases) {
		if (!load_case_file<Scalar, Coord>(c.filePath, coords, weights)) {
			std::cerr << "Error: Failed to load atoms from " << c.filePath << '\n';
			return 1;
		}

		POWERSASA::PowerSasa<Scalar, Coord> ps(coords, weights, 1, 0, 1, 0);
		ps.calc_sasa_all();

		const std::vector<Scalar>& sasa = ps.getSasa();
		const std::vector<Scalar>& vol = ps.getVol();

		const Scalar totalSasa = std::accumulate(sasa.begin(), sasa.end(), Scalar(0));
		const Scalar totalVol = std::accumulate(vol.begin(), vol.end(), Scalar(0));

		const bool sasaOk = std::abs(totalSasa - c.expectedSasa) <= epsilon;
		const bool volOk = std::abs(totalVol - c.expectedVol) <= epsilon;
		allPassed = allPassed && sasaOk && volOk;

		std::println("Case: {}", c.filePath);
		std::println("  Atoms: {}", coords.size());
		std::println("  SASA:   {:>12.6f}  expected {:>12.6f}  [{}]",
		             totalSasa, c.expectedSasa, sasaOk ? "OK" : "FAIL");
		std::println("  Volume: {:>12.6f}  expected {:>12.6f}  [{}]",
		             totalVol, c.expectedVol, volOk ? "OK" : "FAIL");
	}

	if (!allPassed) {
		std::cerr << "REGRESSION FAILURE: one or more cases differ from 912340c golden values." << '\n';
		return 1;
	}

	std::println("REGRESSION PASS: all cases match 912340c golden values.");
	return 0;
}
