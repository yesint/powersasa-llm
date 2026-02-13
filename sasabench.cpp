#include <chrono>
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

using Clock = std::chrono::steady_clock;

struct Frame {
  std::vector<Vec3<float>> coords;
  std::vector<float> weights;
};

bool load_protein(std::vector<Vec3<float>>& coords, std::vector<float>& weights) {
  const fs::path candidates[] = {
      fs::path("protein_coords.txt"),
      fs::path("../protein_coords.txt"),
      fs::path("../../protein_coords.txt"),
  };

  fs::path path;
  for (const auto& c : candidates) {
    if (fs::exists(c)) {
      path = c;
      break;
    }
  }

  if (path.empty()) {
    std::cerr << "Error: Could not find protein_coords.txt\n";
    return false;
  }

  std::ifstream infile(path);
  if (!infile.is_open()) {
    std::cerr << "Error: Could not open " << path << '\n';
    return false;
  }

  coords.clear();
  weights.clear();
  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream iss(line);
    float x = 0.0f, y = 0.0f, z = 0.0f, radius = 0.0f;
    if (!(iss >> x >> y >> z >> radius)) {
      continue;
    }
    coords.emplace_back(x, y, z);
    weights.push_back(radius + 0.14f);
  }

  if (coords.empty()) {
    std::cerr << "Error: No atoms loaded from " << path << '\n';
    return false;
  }

  std::println("Loaded {} atoms from {}", coords.size(), path.string());
  return true;
}

std::vector<Frame> make_frames(const std::vector<Vec3<float>>& baseCoords,
                               const std::vector<float>& baseWeights,
                               int iterations,
                               float maxCoordVariationFraction) {
  std::vector<Frame> frames;
  frames.reserve(static_cast<std::size_t>(iterations));

  for (int it = 0; it < iterations; ++it) {
    Frame frame;
    frame.coords = baseCoords;
    frame.weights = baseWeights;

    const float phase = 0.07f * static_cast<float>(it + 1);
    for (std::size_t i = 0; i < frame.coords.size(); ++i) {
      const float fi = static_cast<float>(i + 1);
      const float dx =
          maxCoordVariationFraction * baseCoords[i].x() * std::sin(0.013f * fi + phase);
      const float dy =
          maxCoordVariationFraction * baseCoords[i].y() * std::cos(0.017f * fi + 1.3f * phase);
      const float dz =
          maxCoordVariationFraction * baseCoords[i].z() * std::sin(0.011f * fi - 0.7f * phase);
      frame.coords[i].x() += dx;
      frame.coords[i].y() += dy;
      frame.coords[i].z() += dz;
    }

    frames.push_back(std::move(frame));
  }

  return frames;
}

struct RunResult {
  double milliseconds = 0.0;
  double checksum = 0.0;
};

RunResult run_reconstruct_each_iter(const std::vector<Frame>& frames) {
  using Scalar = float;
  using Coord = Vec3<Scalar>;

  const auto t0 = Clock::now();
  double checksum = 0.0;
  for (const Frame& frame : frames) {
    POWERSASA::PowerSasa<Scalar, Coord> ps(frame.coords, frame.weights, 1, 0, 1, 0);
    ps.calc_sasa_all();

    const auto& sasa = ps.getSasa();
    const auto& vol = ps.getVol();
    const Scalar totalSasa = std::accumulate(sasa.begin(), sasa.end(), Scalar(0));
    const Scalar totalVol = std::accumulate(vol.begin(), vol.end(), Scalar(0));
    checksum += static_cast<double>(totalSasa + totalVol);
  }
  const auto t1 = Clock::now();

  const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  return RunResult{ms, checksum};
}

RunResult run_update_coords_each_iter(const std::vector<Frame>& frames) {
  using Scalar = float;
  using Coord = Vec3<Scalar>;

  POWERSASA::PowerSasa<Scalar, Coord> ps(frames.front().coords, frames.front().weights, 1, 0, 1, 0);

  const auto t0 = Clock::now();
  double checksum = 0.0;
  for (const Frame& frame : frames) {
    ps.update_coords(frame.coords, frame.weights);
    ps.calc_sasa_all();

    const auto& sasa = ps.getSasa();
    const auto& vol = ps.getVol();
    const Scalar totalSasa = std::accumulate(sasa.begin(), sasa.end(), Scalar(0));
    const Scalar totalVol = std::accumulate(vol.begin(), vol.end(), Scalar(0));
    checksum += static_cast<double>(totalSasa + totalVol);
  }
  const auto t1 = Clock::now();

  const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  return RunResult{ms, checksum};
}

void report(const std::string_view label,
            const RunResult& reconstruct,
            const RunResult& update,
            int iterations) {
  const double perIterReconstruct = reconstruct.milliseconds / static_cast<double>(iterations);
  const double perIterUpdate = update.milliseconds / static_cast<double>(iterations);
  const double speedup = reconstruct.milliseconds / update.milliseconds;
  const double pctFaster = (1.0 - (update.milliseconds / reconstruct.milliseconds)) * 100.0;

  std::println("{}:", label);
  std::println("Reconstruct per iteration:");
  std::println("  total:    {:>10.3f} ms", reconstruct.milliseconds);
  std::println("  per-iter: {:>10.3f} ms", perIterReconstruct);
  std::println("  checksum: {:>14.6f}", reconstruct.checksum);

  std::println("update_coords per iteration:");
  std::println("  total:    {:>10.3f} ms", update.milliseconds);
  std::println("  per-iter: {:>10.3f} ms", perIterUpdate);
  std::println("  checksum: {:>14.6f}", update.checksum);

  if (speedup >= 1.0) {
    std::println("Result: update_coords is faster by {:.2f}% ({:.2f}x).", pctFaster, speedup);
  } else {
    std::println("Result: reconstruct is faster by {:.2f}% ({:.2f}x).",
                 (1.0 - speedup) * 100.0, 1.0 / speedup);
  }
}

}  // namespace

int main(int argc, char** argv) {
  int iterations = 20;
  if (argc >= 2) {
    iterations = std::max(1, std::atoi(argv[1]));
  }

  std::vector<Vec3<float>> baseCoords;
  std::vector<float> baseWeights;
  if (!load_protein(baseCoords, baseWeights)) {
    return 1;
  }

  std::println("Benchmark iterations: {}", iterations);
  std::println("Note: random frame generation is precomputed and excluded from timed sections.");
  const std::vector<Frame> smallFrames = make_frames(baseCoords, baseWeights, iterations, 0.0003f);
  const std::vector<Frame> largeFrames = make_frames(baseCoords, baseWeights, iterations, 0.30f);

  const RunResult reconstructSmall = run_reconstruct_each_iter(smallFrames);
  const RunResult updateSmall = run_update_coords_each_iter(smallFrames);
  report("Scenario A (0.03% coordinate variation; update_coords included)", reconstructSmall,
         updateSmall, iterations);

  const RunResult reconstructLarge = run_reconstruct_each_iter(largeFrames);
  const RunResult updateLarge = run_update_coords_each_iter(largeFrames);
  report("Scenario B (up to 30% coordinate variation; update_coords included)",
         reconstructLarge, updateLarge, iterations);

  return 0;
}
