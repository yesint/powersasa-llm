/*
You are only allow to use this software if you have signed the following license.

ACADEMIC/NON-PROFIT
SASA SOFTWARE LICENSE AGREEMENT
IMPORTANT: This SASA license Agreement is a legal agreement between you, the end user (either an individual or an entity), and the Karlsruhe Institute of Technology.
SASA Software License
GRANT OF LICENSE. Karlsruhe Institute of Technology grants, and you hereby accept, a non-exclusive license to use the SASA software product of the version specified above ("Software") to the extent of its rights and in accordance with the terms of this Agreement. This licensed copy of the Software may only be used on computers at your site by you and members of your organization at your site who have read and agreed to this license. You may install the Software on computers at your site for your own use or use by members of your organization at your site. You may not distribute copies of the Software to others outside of your site. You may make only those copies of the Software which are necessary to install and use it as permitted by this Agreement, or are for purposes of backup and archival records.
OWNERSHIP. This ownership is protected by the copyright laws of the Federal Republic of Germany and by international treaty provisions. Upon expiration or termination of this Agreement, you shall promptly return all copies of the Software and accompanying written materials to the Karlsruhe Institute of Technology.
MODIFICATIONS AND DERIVATIVE WORKS. You may modify the software, and use it to create derivative works, for your internal use at the site covered by this license. You may not distribute such modified or derivative software to others outside of your site without written permission. You may distribute the modifications themselves (e.g. as "patches") under terms of your choice. We encourage users to contribute modifications back into the Software, but you are under no obligation to do so.
REPORTS OF PUBLICATIONS. You agree to acknowledge use of the Software in any reports or publications of results obtained with the Software and cite the publications listed on the download page where you obtained the software in any report or publication in which the Software was used. If you fail to properly acknowledge the use of the software you agree to pay the industrial license fee.
ASSIGNMENT RESTRICTIONS. You shall not use the Software (or any part thereof) in connection with the provision of consultancy, modeling or other services, whether for value or otherwise, on behalf of any third party who does not hold a current valid SASA  Software License Agreement. You shall not use the Software to write other software that duplicates the functionality of the Software. You shall not rent, lease, or otherwise sublet the Software or any part thereof. You may transfer on a permanent basis the rights granted under this license provided you transfer this Agreement and all copies of the Software, including prior versions, and all accompanying written materials. The recipient must agree to the terms of this Agreement in full and register this transfer with the Karlsruhe Institute of Technology.
LIMITED WARRANTY. LICENSEE acknowledges that LICENSORS make no warranty, expressed or implied, that the program will function without error, or in any particular hardware environment, or so as to generate any particular function or result, and further excluding any other warranty, as to the condition of the program, its merchantability, or its fitness for a particular purpose. LICENSORS shall not be liable for any direct, consequential, or other damages suffered by the LICENSEE or any others as a result of their use of the program, whether or not the same could have been foreseen by LICENSORS prior to granting this License. In no event shall LICENSORS liability for any breach of this agreement exceed the fee paid for the license.
KARLSRUHE INSTITUTE OF TECHNOLOGY'S LIABILITY. In no event shall the Karlsruhe Institute of Technology be liable for any indirect, special, or consequential damages, such as, but not limited to, loss of anticipated profits or other economic loss in connection with or arising out of the use of the software by you or the services provided for in this Agreement, even if the Karlsruhe Institute of Technology has been advised of the possibility of such damages. The Karlsruhe Institute of Technology's entire liability and your exclusive remedy shall be, at the Karlsruhe Institute of Technology's discretion, to return the Software and proof of purchase to the Karlsruhe Institute of Technology for either (a) return of any license fee, or (b) correction or replacement of Software that does not meet the terms of this limited warranty.
NO OTHER WARRANTIES. The Karlsruhe Institute of Technology disclaims other implied warranties, including, but not limited to, implied warranties of merchantability or fitness for any purpose, and implied warranties arising by usage of trade, course of dealing, or course of performance. Some states do not allow the limitation of the duration or liability of implied warranties, so the above restrictions might not apply to you.
LICENSE FEE. The software is free for non-profit, government and academic organizations. For-profit and commercial organizations wishing to license SASA shall contact:
SASA-support@kit.edu
and request a quote for a commercial license.
Supercomputer centers can license SASA under the same conditions and make it available to their users from non-profit organizations as executable code. However, for-profit organizations who want to use the program at supercomputer centers must have signed a separate license agreement.


If you have no license please contact SASA-support@kit.edu
*/
#ifndef POWER_DIAGRAM_H_
#define POWER_DIAGRAM_H_

#define __power_diagram_internal_timing__ 0
#ifndef PD_ENABLE_TOPOLOGY_ASSERTS
#define PD_ENABLE_TOPOLOGY_ASSERTS 0
#endif
#include <array>
#include <iostream>
#include <print>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <deque>
#include <algorithm>
#include <cstddef>
#include <limits>
#include <ctime>
// #include "basic_vector_calc.h"

namespace POWER_DIAGRAM
{

	class MyException : public std::exception
	{
	public:
		MyException() {}
	};
	class IdenticalPointException : public std::exception {};
	class VerticesFullException : public std::exception {};

	// Return the nth natural number of a sequence where one index ("without") is skipped.
	// Used to compare generator triples while omitting one generator position.
	inline int nth(const int n, const int without)
	{
		return n + (without <= n);
	}

	template <typename T>
	inline std::string pd_to_string(const T &value)
	{
		std::ostringstream oss;
		oss << value;
		return oss.str();
	}

	template <class PDCoord, class PDFloat, class Pos_iterator, class Strength_iterator, const int dimension>
	void getBoundingBox(PDCoord &lowestCorner, PDCoord &highestCorner, const unsigned int &size, const Pos_iterator pos_begin, const Strength_iterator strength_begin, const PDFloat additionalCubeSize = pow(2.0, 1.0 / dimension) - 1)
	{
		// Build a radius-aware axis-aligned box covering all weighted points.
		// The optional inflation factor keeps the initial clipping cube safely outside the data.
		if (size > 0)
		{
			lowestCorner = *pos_begin;
			highestCorner = *pos_begin;

			Pos_iterator pos_end = pos_begin;
			Strength_iterator strength_end = strength_begin;

			for (unsigned int i = 0; i < size; i++)
			{
				for (int g = dimension - 1; g >= 0; g--)
				{
					if ((*pos_end)[g] - (*strength_end) < lowestCorner[g])
						lowestCorner[g] = (*pos_end)[g] - (*strength_end);
					if ((*pos_end)[g] + (*strength_end) > highestCorner[g])
						highestCorner[g] = (*pos_end)[g] + (*strength_end);
				}
				++pos_end;
				++strength_end;
			}

			const PDCoord center = 0.5 * (lowestCorner + highestCorner);
			lowestCorner += /*PowerDiagram<PDFloat,PDCoord,dimension>::error((highestCorner+lowestCorner).norm())**/ (lowestCorner - center) * additionalCubeSize;
			highestCorner += /*PowerDiagram<PDFloat,PDCoord,dimension>::error((highestCorner+lowestCorner).norm())**/ (highestCorner - center) * additionalCubeSize;
		}
	}

	struct PowerDiagramRuntimeParams
	{
		bool radiiGiven;	  // otherwise the second input array is interpreted as power
		bool fill_myVertices; // otherwise no cell information is constructed, but only a diagram of vertices
		bool fill_neighbours; // otherwise no neighbourhood-information is extracted from myVertices
		bool fill_zeroPoints; // otherwise no information where power diagram is zero is created
		bool with_warnings;	  // otherwise powerdiagram is in silent mode and will not tell when reducing powers
		bool without_check;
		PowerDiagramRuntimeParams(const bool _radiiGiven, const bool _fill_myVertices, const bool _fill_neighbours, const bool _fill_zeroPoints, const bool _with_warnings, bool _without_check)
		{
			radiiGiven = _radiiGiven;
			fill_myVertices = _fill_myVertices;
			fill_neighbours = _fill_neighbours;
			fill_zeroPoints = _fill_zeroPoints;
			with_warnings = _with_warnings;
			without_check = _without_check;
		}
	};

	template <class PDFloat, class PDCoord, class Pos_iterator, class Strength_iterator, class BondTo_iterator>
	struct PowerDiagramParams
	{
		bool create_vertices; // otherwise the object is only constructed, but no calculation is done
		PowerDiagramRuntimeParams runpar;

		const PDCoord lowestCorner;
		const PDCoord highestCorner;

		const unsigned int size;
		Pos_iterator pos_begin;
		Strength_iterator strength_begin;
		BondTo_iterator bondTo_begin;

		PowerDiagramParams with_radiiGiven(const bool yon)
		{
			runpar.radiiGiven = yon;
			return *this;
		};
		PowerDiagramParams with_calculate(const bool yon)
		{
			create_vertices = yon;
			return *this;
		};
		PowerDiagramParams with_myVertices(const bool yon)
		{
			runpar.fill_myVertices = yon;
			return *this;
		};
		PowerDiagramParams with_cells(const bool yon)
		{
			runpar.fill_neighbours = yon;
			runpar.fill_myVertices = yon;
			return *this;
		};
		PowerDiagramParams with_zeroPoints(const bool yon)
		{
			runpar.fill_zeroPoints = yon;
			runpar.fill_myVertices = yon;
			return *this;
		};
		PowerDiagramParams with_Warnings(const bool yon)
		{
			runpar.with_warnings = yon;
			return *this;
		};
		PowerDiagramParams without_Check(const bool yon)
		{
			runpar.without_check = yon;
			return *this;
		};
		PowerDiagramParams(const unsigned int size_, Pos_iterator &pos_begin_, Strength_iterator &strength_begin_, BondTo_iterator bondTo_begin_, PDCoord &lc, PDCoord &hc,
						   bool radiiGiven_ = 1, bool create_vertices_ = 1, bool fill_cellVertices = 1, bool fill_cellNeighbours = 1, bool fill_zeroPoints = 1, bool withWarnings = 1, bool withoutChecks = 1) 
						   : create_vertices(create_vertices_), runpar(PowerDiagramRuntimeParams(radiiGiven_, fill_cellVertices, fill_cellNeighbours, fill_zeroPoints, withWarnings, withoutChecks)), lowestCorner(lc), highestCorner(hc),
							 size(size_), pos_begin(pos_begin_), strength_begin(strength_begin_), bondTo_begin(bondTo_begin_)
		{
		}
	};

	template <class PDFloat, class PDCoord, const int dimension>
	class PowerDiagram
	{
	public:
		// Forward declarations for cell and vertex
		struct zeroPoint;
		struct vertex;
		struct cell;
		struct EdgeEnds;

		using CellId = std::size_t;
		using VertexId = std::size_t;
		using ZeroId = std::size_t;
		inline static constexpr std::size_t kInvalidId = std::numeric_limits<std::size_t>::max();
		enum class GeneratorKind
		{
			point,
			side
		};

		struct GeneratorRef
		{
			GeneratorKind kind;
			std::size_t index;
			constexpr GeneratorRef() : kind(GeneratorKind::point), index(kInvalidId) {}
			constexpr GeneratorRef(const GeneratorKind kind_, const std::size_t index_) : kind(kind_), index(index_) {}
			constexpr bool is_valid() const { return index != kInvalidId; }
		};

		PDCoord center;
		PDFloat maxr2;

	private:
		PowerDiagramRuntimeParams params;
		unsigned int _nVertices;
		unsigned int _nUnused;

		std::vector<VertexId> unused;
		std::vector<cell> points;
		std::vector<vertex> vertices; // we use the plane space of a vector because of speed. never push_back or it will cause a realloc !!!
		std::vector<zeroPoint> zeros; // this is a vector of all zero points, but maybe some are inactive (inactives do not appear in myZeros!)
		std::vector<cell> sideGenerators;

		unsigned int nRevertVertices;
		unsigned int nRevertZeros;
		unsigned int nRevertPoints;
		std::array<CellId, 1 << dimension> cornerOwners{};

		PDFloat powerErr;
		PDFloat insertionErrorScale;
		std::vector<VertexId> ReplacedIds;		// old vertices that are removed (by ID)
		std::vector<VertexId> Invalids;			// old vertices that are removed reloaded
		std::vector<GeneratorRef> InvolvedRefs; // all cells which are involved (point or side refs)
		std::vector<EdgeEnds> planes;			// used for connecting new vertices, stores "open ends"

		enum class ReplaceState
		{
			persisting,
			replaced,
			ambiguous
		};

	public:
		struct cell
		{
			int visitedAs;
			PDCoord position;
			PDFloat r;
			PDFloat r2;
			CellId bondToId;
			std::vector<CellId> neighboursIds;
			std::vector<VertexId> myVerticesIds;
			std::vector<int> myZeroPoints;

			inline cell(PDCoord const &pos, PDFloat const &root, PDFloat const &power) : visitedAs(0), position(pos), r(root), r2(power), bondToId(kInvalidId)
			{
				myVerticesIds.reserve(12);
			}
			inline cell(PDCoord const &pos, const PDFloat &str) : visitedAs(0), position(pos), r(str), r2(str * str), bondToId(kInvalidId)
			{
				myVerticesIds.reserve(12);
			}
			inline PDFloat power(const PDCoord &coord) const
			{
				// Power distance to this weighted site.
				return (position - coord).squaredNorm() - r2;
			}

		private:
		};
		inline void AddToInvolved(const GeneratorRef &ref)
		{
			// Geometry:
			//   Tag one generator cell as part of the active local insertion neighborhood.
			// Code structure:
			//   Converts ref->cell, stores traversal mark in visitedAs, and appends ref to InvolvedRefs.
			cell &thit = cell_from_ref(ref);
			thit.visitedAs = involved_size();
			push_involved(ref);
		}
		clock_t t1, t2, t3, t4, t5, t6;
		cell const &get_point(const int n) const { return points[n]; }
		std::vector<cell> const &get_points() const { return points; }
		inline cell &cell_at(const CellId id)
		{
			return points[id];
		}
		inline const cell &cell_at(const CellId id) const
		{
			return points[id];
		}
		inline vertex &vertex_at(const VertexId id)
		{
			return vertices[id];
		}
		inline const vertex &vertex_at(const VertexId id) const
		{
			return vertices[id];
		}
		inline bool ref_is_real_point(const GeneratorRef &ref) const
		{
			// Geometry:
			//   Distinguish true input generators from synthetic side generators.
			// Code structure:
			//   Checks GeneratorRef validity, kind, and bounds against points[].
			return ref.is_valid() && ref.kind == GeneratorKind::point && ref.index < points.size();
		}
		inline void sync_cell_link_mirrors(cell &a_cell)
		{
			// Geometry:
			//   No geometric change; clamp stale adjacency/ownership ids after edits.
			// Code structure:
			//   Bounds-checks bondToId, neighboursIds, and myVerticesIds against current vector sizes.
			if (a_cell.bondToId >= points.size())
				a_cell.bondToId = kInvalidId;
			for (std::size_t i = 0; i < a_cell.neighboursIds.size(); ++i)
				if (a_cell.neighboursIds[i] >= points.size())
					a_cell.neighboursIds[i] = kInvalidId;
			for (std::size_t i = 0; i < a_cell.myVerticesIds.size(); ++i)
				if (a_cell.myVerticesIds[i] >= vertices.size())
					a_cell.myVerticesIds[i] = kInvalidId;
		}
		inline void sync_vertex_link_mirrors(vertex &a_vertex) const
		{
			// Geometry:
			//   No geometric change; sanitize vertex generator/endpoint references.
			// Code structure:
			//   Invalidates out-of-range GeneratorRefs and endpoint ids.
			for (int g = 0; g <= dimension; ++g)
			{
				const GeneratorRef &ref = a_vertex.generatorRefs[g];
				if (ref.kind == GeneratorKind::point && ref.index >= points.size())
					a_vertex.generatorRefs[g] = GeneratorRef();
				if (ref.kind == GeneratorKind::side && ref.index >= sideGenerators.size())
					a_vertex.generatorRefs[g] = GeneratorRef();
				if (a_vertex.endPointIds[g] >= vertices.size())
					a_vertex.endPointIds[g] = kInvalidId;
			}
		}
		inline void sync_zero_link_mirrors(zeroPoint &a_zero) const
		{
			// Geometry:
			//   No geometric change; sanitize zero-point origin and generator refs.
			// Code structure:
			//   Invalidates fromId/ref entries that no longer point into owner arrays.
			if (a_zero.fromId >= vertices.size())
				a_zero.fromId = kInvalidId;
			for (int g = 0; g < dimension; ++g)
			{
				const GeneratorRef &ref = a_zero.generatorRefs[g];
				if (ref.kind == GeneratorKind::point && ref.index >= points.size())
					a_zero.generatorRefs[g] = GeneratorRef();
				if (ref.kind == GeneratorKind::side && ref.index >= sideGenerators.size())
					a_zero.generatorRefs[g] = GeneratorRef();
			}
		}
		inline void push_zero_from_edge(const VertexId source_id, const int branch, const PDFloat sol)
		{
			// Geometry:
			//   Record one edge-parameterized zero crossing event.
			// Code structure:
			//   Copies the 3 defining generator refs from source vertex excluding branch slot.
			const vertex &source_vertex = vertex_at(source_id);
			zeroPoint zp(
				sol,
				branch,
				source_vertex.generatorRefs[nth(0, branch)],
				source_vertex.generatorRefs[nth(1, branch)],
				source_vertex.generatorRefs[nth(2, branch)],
				source_id);
			zeros.push_back(zp);
		}
		inline void sync_all_link_mirrors()
		{
			// Geometry:
			//   No geometric change; global ID/ref sanitation pass.
			// Code structure:
			//   Applies sync_* helpers to all cells, active vertices, and zero points.
			for (cell &a_cell : points)
				sync_cell_link_mirrors(a_cell);
			for (unsigned int vi = 0; vi < _nVertices; ++vi)
				sync_vertex_link_mirrors(vertices[vi]);
			for (zeroPoint &a_zero : zeros)
				sync_zero_link_mirrors(a_zero);
		}
		inline void validate_transient_mirror_invariants() const
		{
			// Geometry:
			//   Debug-only integrity assertion for transient insertion containers.
			// Code structure:
			//   Ensures involved refs and replaced vertex ids are valid when asserts are enabled.
#if PD_ENABLE_TOPOLOGY_ASSERTS
			for (std::size_t i = 0; i < involved_size(); ++i)
			{
				if (!valid_generator_ref(InvolvedRefs[i]))
					throw MyException();
			}
			for (std::size_t i = 0; i < ReplacedIds.size(); ++i)
			{
				if (ReplacedIds[i] == kInvalidId)
					throw MyException();
				if (!valid_vertex_id(ReplacedIds[i]))
					throw MyException();
			}
#endif
		}
		inline void validate_cell_mirror_invariants(const cell &a_cell) const
		{
			// Geometry:
			//   Debug-only integrity assertion for one cell link mirror.
			// Code structure:
			//   Verifies bounds of myVerticesIds/neighboursIds/bondToId under assert guard.
#if PD_ENABLE_TOPOLOGY_ASSERTS
			for (std::size_t i = 0; i < a_cell.myVerticesIds.size(); ++i)
				if (a_cell.myVerticesIds[i] != kInvalidId && a_cell.myVerticesIds[i] >= vertices.size())
					throw MyException();
			for (std::size_t i = 0; i < a_cell.neighboursIds.size(); ++i)
				if (a_cell.neighboursIds[i] != kInvalidId && a_cell.neighboursIds[i] >= points.size())
					throw MyException();
			if (a_cell.bondToId != kInvalidId && a_cell.bondToId >= points.size())
				throw MyException();
#endif
		}
		inline void validate_all_cell_mirror_invariants() const
		{
#if PD_ENABLE_TOPOLOGY_ASSERTS
			for (const cell &a_cell : points)
				validate_cell_mirror_invariants(a_cell);
#endif
		}
		inline void validate_phase_mirror_invariants() const
		{
#if PD_ENABLE_TOPOLOGY_ASSERTS
			validate_transient_mirror_invariants();
			validate_all_cell_mirror_invariants();
#endif
		}
		inline bool valid_cell_id(const CellId id) const { return id != kInvalidId && id < points.size(); }
		inline bool valid_vertex_id(const VertexId id) const { return id != kInvalidId && id < vertices.size(); }
		inline bool valid_generator_ref(const GeneratorRef &ref) const
		{
			// Geometry:
			//   Determine whether generator ref maps to a currently existing real/side cell.
			// Code structure:
			//   Checks kind-specific bounds in points[] or sideGenerators[].
			if (!ref.is_valid())
				return false;
			if (ref.kind == GeneratorKind::point)
				return ref.index < points.size();
			if (ref.kind == GeneratorKind::side)
				return ref.index < sideGenerators.size();
			return false;
		}
		inline const cell &cell_from_ref(const GeneratorRef &ref) const
		{
			if (!valid_generator_ref(ref))
				throw MyException();
			return (ref.kind == GeneratorKind::point) ? points[ref.index] : sideGenerators[ref.index];
		}
		inline cell &cell_from_ref(const GeneratorRef &ref)
		{
			if (!valid_generator_ref(ref))
				throw MyException();
			return (ref.kind == GeneratorKind::point) ? points[ref.index] : sideGenerators[ref.index];
		}
		inline void set_vertex_generator(vertex &a_vertex, const int slot, const GeneratorRef &ref)
		{
			// Geometry:
			//   Assign one generator in vertex-defining tuple.
			// Code structure:
			//   Direct write to generatorRefs[slot].
			a_vertex.generatorRefs[slot] = ref;
		}
		inline void set_vertex_endpoint(vertex &a_vertex, const int slot, const VertexId id)
		{
			// Geometry:
			//   Assign one endpoint link of a vertex edge branch.
			// Code structure:
			//   Direct write to endPointIds[slot].
			a_vertex.endPointIds[slot] = id;
		}
		inline void set_vertex_endpoint_deferred(vertex &a_vertex, const int slot, const VertexId id)
		{
			a_vertex.endPointIds[slot] = id;
		}
		inline void swap_vertex_link_slots(vertex &a_vertex, const int a, const int b)
		{
			// Geometry:
			//   Permute local generator/endpoint slot ordering without changing represented vertex.
			// Code structure:
			//   Swaps generatorRefs and endPointIds in lockstep.
			std::swap(a_vertex.generatorRefs[a], a_vertex.generatorRefs[b]);
			std::swap(a_vertex.endPointIds[a], a_vertex.endPointIds[b]);
		}
		inline void clear_bond_to(cell &a_cell)
		{
			a_cell.bondToId = kInvalidId;
			validate_cell_mirror_invariants(a_cell);
		}
		inline void set_bond_to_id(cell &a_cell, const CellId id)
		{
			a_cell.bondToId = (id < points.size()) ? id : kInvalidId;
			validate_cell_mirror_invariants(a_cell);
		}
		inline void clear_replaced()
		{
			// Geometry:
			//   Reset replaced-vertex frontier container.
			// Code structure:
			//   Clears ReplacedIds and validates transient invariants.
			ReplacedIds.clear();
			validate_transient_mirror_invariants();
		}
		inline void push_replaced_id(const VertexId id)
		{
			ReplacedIds.push_back(id);
			validate_transient_mirror_invariants();
		}
		inline void set_replaced_id(const std::size_t index, const VertexId id)
		{
			ReplacedIds[index] = id;
			validate_transient_mirror_invariants();
		}
		inline void pop_replaced()
		{
			ReplacedIds.pop_back();
			validate_transient_mirror_invariants();
		}
		inline VertexId replaced_id_at(const std::size_t index) const
		{
			if (index >= ReplacedIds.size())
				return kInvalidId;
			return ReplacedIds[index];
		}
		inline std::size_t replaced_size() const
		{
			return ReplacedIds.size();
		}
		inline void clear_involved()
		{
			InvolvedRefs.clear();
			validate_transient_mirror_invariants();
		}
		inline void sort_involved_by_ref()
		{
			// Geometry:
			//   Canonicalize involved-generator ordering for deterministic local processing.
			// Code structure:
			//   Stable rank by (kind, index) through index permutation then container rewrite.
			std::vector<std::size_t> order(InvolvedRefs.size());
			for (std::size_t i = 0; i < order.size(); ++i)
				order[i] = i;
			auto ref_rank = [](const GeneratorRef &ref) -> int
			{
				if (!ref.is_valid())
					return 2;
				return (ref.kind == GeneratorKind::point) ? 0 : 1;
			};
			std::sort(order.begin(), order.end(),
					  [&](const std::size_t a, const std::size_t b)
					  {
						  const GeneratorRef &ra = InvolvedRefs[a];
						  const GeneratorRef &rb = InvolvedRefs[b];
						  const int ka = ref_rank(ra);
						  const int kb = ref_rank(rb);
						  if (ka != kb)
							  return ka < kb;
						  return ra.index < rb.index;
					  });
			std::vector<GeneratorRef> sorted_refs;
			sorted_refs.reserve(InvolvedRefs.size());
			for (const std::size_t idx : order)
			{
				sorted_refs.push_back(InvolvedRefs[idx]);
			}
			InvolvedRefs.swap(sorted_refs);
			validate_transient_mirror_invariants();
		}
		inline void push_involved(const GeneratorRef &ref)
		{
			// Geometry:
			//   Add one generator to active local neighborhood.
			// Code structure:
			//   Pushes into InvolvedRefs and revalidates transient invariants.
			InvolvedRefs.push_back(ref);
			validate_transient_mirror_invariants();
		}
		inline CellId involved_id_at(const std::size_t index) const
		{
			if (index >= InvolvedRefs.size())
				return kInvalidId;
			const GeneratorRef &ref = InvolvedRefs[index];
			if (ref.kind != GeneratorKind::point)
				return kInvalidId;
			if (ref.index >= points.size())
				return kInvalidId;
			return static_cast<CellId>(ref.index);
		}
		inline std::size_t involved_size() const
		{
			return InvolvedRefs.size();
		}
		inline void clear_cell_my_vertices(cell &a_cell)
		{
			// Geometry:
			//   Clear inverse cell->vertex incidence for one cell.
			// Code structure:
			//   Empties myVerticesIds and validates cell invariants.
			a_cell.myVerticesIds.clear();
			validate_cell_mirror_invariants(a_cell);
		}
		inline void push_cell_my_vertex(cell &a_cell, const VertexId id)
		{
			a_cell.myVerticesIds.push_back(id);
			validate_cell_mirror_invariants(a_cell);
		}
		inline void pop_cell_my_vertex(cell &a_cell)
		{
			a_cell.myVerticesIds.pop_back();
			validate_cell_mirror_invariants(a_cell);
		}
		inline void set_cell_my_vertex(cell &a_cell, const std::size_t index, const VertexId id)
		{
			a_cell.myVerticesIds[index] = id;
			validate_cell_mirror_invariants(a_cell);
		}
		inline void erase_cell_my_vertex(cell &a_cell, const std::size_t index)
		{
			a_cell.myVerticesIds.erase(a_cell.myVerticesIds.begin() + index);
			validate_cell_mirror_invariants(a_cell);
		}
		inline void erase_cell_my_vertex_by_id(cell &a_cell, const VertexId id)
		{
			for (std::size_t i = 0; i < a_cell.myVerticesIds.size(); ++i)
			{
				if (a_cell.myVerticesIds[i] != id)
					continue;
				erase_cell_my_vertex(a_cell, i);
				return;
			}
		}
		inline void clear_cell_neighbours(cell &a_cell)
		{
			// Geometry:
			//   Clear cell-cell adjacency list.
			// Code structure:
			//   Empties neighboursIds and validates cell invariants.
			a_cell.neighboursIds.clear();
			validate_cell_mirror_invariants(a_cell);
		}
		inline void push_cell_neighbour(cell &a_cell, const CellId id)
		{
			a_cell.neighboursIds.push_back(id);
			validate_cell_mirror_invariants(a_cell);
		}
		inline void set_cell_neighbour(cell &a_cell, const std::size_t index, const CellId id)
		{
			a_cell.neighboursIds[index] = id;
			validate_cell_mirror_invariants(a_cell);
		}
		inline void assign_cell_neighbours_from_ids(cell &a_cell, const std::vector<CellId> &neighbour_ids)
		{
			a_cell.neighboursIds = neighbour_ids;
			validate_cell_mirror_invariants(a_cell);
		}
		inline void erase_cell_neighbour(cell &a_cell, const std::size_t index)
		{
			a_cell.neighboursIds.erase(a_cell.neighboursIds.begin() + index);
			validate_cell_mirror_invariants(a_cell);
		}
		inline zeroPoint const &get_zero(const ZeroId n) const { return zeros[n]; }
		inline zeroPoint &get_zero(const ZeroId n) { return zeros[n]; }
		inline vertex const &get_vertex(const VertexId n) const { return vertex_at(n); }
		inline vertex &get_vertex(const VertexId n) { return vertex_at(n); }
		inline cell const &get_cell(const CellId n) const { return cell_at(n); }
		inline cell &get_cell(const CellId n) { return cell_at(n); }
		inline cell const &get_generator(const GeneratorRef &ref) const
		{
			return (ref.kind == GeneratorKind::point) ? points[ref.index] : sideGenerators[ref.index];
		}
		inline cell &get_generator(const GeneratorRef &ref)
		{
			return (ref.kind == GeneratorKind::point) ? points[ref.index] : sideGenerators[ref.index];
		}
		inline GeneratorRef get_point_ref(const CellId index) const { return GeneratorRef(GeneratorKind::point, index); }
		inline GeneratorRef get_side_ref(const std::size_t index) const { return GeneratorRef(GeneratorKind::side, index); }
		unsigned int get_point_num(const CellId my_cell_id)
		{
			return static_cast<unsigned int>(my_cell_id);
		}
		std::vector<vertex> const &get_vertices() const { return vertices; }

		std::vector<zeroPoint> const &get_zeroPoints() const { return zeros; }
		inline VertexId zeroPointFromId(const zeroPoint &zp) const
		{
			return (zp.fromId != kInvalidId && zp.fromId < _nVertices) ? zp.fromId : kInvalidId;
		}
		inline bool zeroPointValid(const zeroPoint &zp) const
		{
			// Geometry:
			//   Test whether a stored zero-crossing still references two valid active edge endpoints.
			// Code structure:
			//   Validates fromId, branch endpoint id, and endpoint connectivity flags.
			const VertexId from_id = zeroPointFromId(zp);
			if (from_id == kInvalidId)
				return false;
			const vertex &from = vertex_at(from_id);
			VertexId to_id = from.resolved_endpoint_id(*this, zp.branch);
			if (to_id == kInvalidId || to_id >= _nVertices)
				return false;
			const vertex &to = vertex_at(to_id);
			return ((!from.invalid) && (!to.invalid));
		}
		inline PDCoord zeroPointPos(const zeroPoint &zp) const
		{
			// Geometry:
			//   Evaluate 3D position of zero crossing by affine interpolation along its source edge.
			// Code structure:
			//   Resolves source/target endpoints and computes linear blend via zp.pos.
			const VertexId from_id = zeroPointFromId(zp);
			if (from_id == kInvalidId)
				throw MyException();
			const vertex &from = vertex_at(from_id);
			VertexId to_id = from.resolved_endpoint_id(*this, zp.branch);
			if (to_id == kInvalidId || to_id >= _nVertices)
				throw MyException();
			const vertex &to = vertex_at(to_id);
			return to.position * zp.pos - from.position * (zp.pos - static_cast<PDFloat>(1.0));
		}

		inline static PDFloat error(const PDFloat &f)
		{
			// Scale machine epsilon to the magnitude of f with a non-zero floor.
			if (f > std::numeric_limits<PDFloat>::min() / std::numeric_limits<PDFloat>::epsilon())
				return f * (std::numeric_limits<PDFloat>::epsilon());
			else if (f < -std::numeric_limits<PDFloat>::min() / std::numeric_limits<PDFloat>::epsilon())
				return -f * (std::numeric_limits<PDFloat>::epsilon());
			else
				return std::numeric_limits<PDFloat>::min() / std::numeric_limits<PDFloat>::epsilon();
		}
		inline bool above_power_err(const PDFloat value) const { return value > powerErr; }
		inline bool below_neg_power_err(const PDFloat value) const { return value < -powerErr; }
		inline bool within_power_err(const PDFloat value) const { return std::abs(value) <= powerErr; }
		inline bool below_power_err(const PDFloat value) const { return value < powerErr; }
		inline PDFloat power_err_scaled_epsilon() const { return powerErr * std::numeric_limits<PDFloat>::epsilon(); }
		template <class Pos_iterator, class Strength_iterator, class BondTo_iterator>
		static PowerDiagramParams<PDFloat, PDCoord, Pos_iterator, Strength_iterator, BondTo_iterator> create(unsigned int size, Pos_iterator pos_begin, Strength_iterator strength_begin, BondTo_iterator bondTo_begin)
		{
			// Geometry:
			//   Build initial global clipping bounds from weighted input sites.
			// Code structure:
			//   Computes bounding box and returns deferred parameter object for fluent construction flags.
			PDCoord highestCorner;
			PDCoord lowestCorner;
			{
				if (size >= 1)
					getBoundingBox<PDCoord, PDFloat, Pos_iterator, Strength_iterator, 3>(lowestCorner, highestCorner, size, pos_begin, strength_begin);
				else
				{
					std::println("create empty PD(not implemented, yet)");
					throw MyException();
				}
				return PowerDiagramParams<PDFloat, PDCoord, Pos_iterator, Strength_iterator, BondTo_iterator>(size, pos_begin, strength_begin, bondTo_begin, lowestCorner, highestCorner);
			}
		}
		template <typename Pos_iterator, typename Strength_iterator, typename BondTo_iterator>
		PowerDiagram(PowerDiagramParams<PDFloat, PDCoord, Pos_iterator, Strength_iterator, BondTo_iterator> _params) : center(0.5 * (_params.highestCorner + _params.lowestCorner)), params(_params.runpar)
		{
			// Build initial clipping cube, insert all points, then derive optional adjacency/zero-point caches.
			_nUnused = 0;
			insertionErrorScale = 0;
			cornerOwners.fill(kInvalidId);
			nRevertPoints = 0;
			nRevertZeros = 0;
			nRevertVertices = (1 << dimension);
			planes.resize(64 * 64);
			points.reserve(_params.size);
			vertices.reserve(_params.size * 32 + (1 << dimension));
			vertices.resize(vertices.capacity());

			PDCoord lowest = _params.lowestCorner;
			PDCoord highest = _params.highestCorner;
			buildCube(PDCoord(lowest - center), PDCoord(highest - center));
			// build and connect corners

			Pos_iterator pos_it = _params.pos_begin;
			Strength_iterator strength_it = _params.strength_begin;
			BondTo_iterator bondTo_it = _params.bondTo_begin;

			if (_params.size > 0)
			{
				if (params.radiiGiven)
				{
					points.push_back(cell(*pos_it - center, *strength_it));
					clear_bond_to(points.back());
					for (unsigned int i = 1; i < _params.size; i++)
					{
						++pos_it;
						++strength_it;
						++bondTo_it;
						const unsigned int bond_to = static_cast<unsigned int>(*bondTo_it);
						points.push_back(cell(*pos_it - center, *strength_it));
						set_bond_to_id(points.back(), bond_to);
					}
				}
				else
				{
					points.push_back(cell(*pos_it - center, sqrt(*strength_it), *strength_it));
					clear_bond_to(points.back());
					for (unsigned int i = 1; i < _params.size; i++)
					{
						++pos_it;
						++strength_it;
						++bondTo_it;
						const unsigned int bond_to = static_cast<unsigned int>(*bondTo_it);
						points.push_back(cell(*pos_it - center, sqrt(*strength_it), *strength_it));
						set_bond_to_id(points.back(), bond_to);
					}
				}
			}

			if (_params.create_vertices)
				buildVertices(points.size());

			if (params.fill_myVertices)
				FillAllMyVertices();
			if (params.fill_neighbours)
				FillAllNeighbours();
			if (params.fill_zeroPoints)
				FillAllZeroPoints();
			sync_all_link_mirrors();
		}
		const std::vector<cell> &getPoints() const { return points; }

		void revert()
		{
			// Geometry:
			//   Undo incremental insertions and restore the exact previously valid power-diagram state
			//   (corner ownership, finite-vertex connectivity, and derived adjacency/zero sets).
			// Code structure:
			//   Uses snapshot fields captured in addMore() (nRevert* and cornerOwners), rebuilds mirror IDs,
			//   restores corner generators, reconnects preserved invalidated vertices, and trims inserted points/zeros.
			for (unsigned int vi = nRevertVertices; vi < _nVertices; ++vi)
				for (int gi = 1; gi <= dimension; ++gi)
					if (vertices[vi].generatorRefs[gi].kind == GeneratorKind::point && vertices[vi].generatorRefs[gi].index < points.size())
						pop_cell_my_vertex(points[vertices[vi].generatorRefs[gi].index]);
			_nVertices = nRevertVertices;
			clear_involved();
			if (points.size() > nRevertPoints)
				points.erase(points.begin() + nRevertPoints, points.end());
			for (int c = 0; c < (1 << dimension); c++)
			{
				const CellId owner_id = cornerOwners[c];
				if (owner_id == kInvalidId || owner_id >= points.size())
					throw MyException();
				cell &owner = cell_at(owner_id);
				set_vertex_generator(vertices[c], 0, GeneratorRef(GeneratorKind::point, owner_id));
				vertices[c].powerValue = cell_from_ref(vertices[c].generatorRefs[0]).power(vertices[c].position);
			}
			for (const VertexId invalid_id : Invalids)
			{
				if (invalid_id == kInvalidId || invalid_id >= vertices.size())
					continue;
				vertex &it = vertex_at(invalid_id);
				it.invalid = 0;
				it.rrv = 0;
				for (int endpoint_idx = it.isCorner(); endpoint_idx <= dimension; ++endpoint_idx)
				{
					VertexId endpoint_id = it.resolved_endpoint_id(*this, endpoint_idx);
					if (endpoint_id == kInvalidId || endpoint_id >= vertices.size())
						continue;
					vertex &endpoint = vertex_at(endpoint_id);
					for (int g1 = it.isCorner(); g1 <= dimension; g1++)
						for (int g2 = endpoint.isCorner(); g2 <= dimension; g2++)
						{
							if (it.generatorRefs[nth(0, g1)].kind == endpoint.generatorRefs[nth(0, g2)].kind && it.generatorRefs[nth(0, g1)].index == endpoint.generatorRefs[nth(0, g2)].index && it.generatorRefs[nth(1, g1)].kind == endpoint.generatorRefs[nth(1, g2)].kind && it.generatorRefs[nth(1, g1)].index == endpoint.generatorRefs[nth(1, g2)].index && it.generatorRefs[nth(2, g1)].kind == endpoint.generatorRefs[nth(2, g2)].kind && it.generatorRefs[nth(2, g1)].index == endpoint.generatorRefs[nth(2, g2)].index)
							{
								//						connection is already set
								set_vertex_endpoint(endpoint, g2, invalid_id);
							}
						}
				}
				for (int g = 0; g <= dimension; ++g)
					if (it.generatorRefs[g].kind == GeneratorKind::point && it.generatorRefs[g].index < points.size())
					{
						cell &generator = points[it.generatorRefs[g].index];
						push_cell_my_vertex(generator, invalid_id);
						if (generator.visitedAs == 0)
						{
							generator.visitedAs = -1;
							push_involved(GeneratorRef(GeneratorKind::point, it.generatorRefs[g].index));
						}
					}
			}
			Invalids.clear();
			if (params.fill_neighbours)
				FillAllNeighboursOfInvolved();

			if (params.fill_zeroPoints)
			{
				for (std::size_t involved_idx = 0; involved_idx < involved_size(); ++involved_idx)
				{
					const CellId involved_id = involved_id_at(involved_idx);
					if (involved_id == kInvalidId)
						continue;
					cell &involved_cell = cell_at(involved_id);
					while (!involved_cell.myZeroPoints.empty() && involved_cell.myZeroPoints.back() > static_cast<int>(nRevertZeros))
						involved_cell.myZeroPoints.pop_back();
				}
				zeros.erase(zeros.begin() + nRevertZeros, zeros.end());
			}
			nRevertPoints = 0;
			nRevertVertices = 0;
			for (int c = 0; c < (1 << dimension); c++)
				cornerOwners[c] = kInvalidId;
			sync_all_link_mirrors();
			validate_phase_mirror_invariants();
		}
		template <class Pos_iterator, class Strength_iterator>
		inline void recalculate(const Pos_iterator pos_it, const Strength_iterator strength_it, const unsigned int size)
		// does standard deletion and calculation of Vertices, neighbour information, ...
		{
			// Geometry:
			//   Recompute the complete clipped power diagram for a new full point set.
			// Code structure:
			//   Resets transient and per-cell caches, rebuilds the cube, updates all sites, then recomputes
			//   vertices / myVertices / neighbours / zero points in that strict dependency order.
			clearAllmyVertices();
			clear_interna();
			if (size > points.size())
			{
				std::vector<CellId> old_bond_ids(points.size(), kInvalidId);
				for (std::size_t i = 1; i < points.size(); ++i)
				{
					old_bond_ids[i] = points[i].bondToId;
				}
				const std::size_t old_capacity = points.capacity();
				points.reserve(size);
				if (points.capacity() != old_capacity)
					for (std::size_t i = 1; i < points.size(); ++i)
					{
						if (old_bond_ids[i] != kInvalidId)
							set_bond_to_id(points[i], old_bond_ids[i]);
					}
			}
			zeros.clear();
			PDCoord lowest, highest;
			getBoundingBox<PDCoord, PDFloat, PDCoord const *, PDFloat const *, dimension>(lowest, highest, size, &(*pos_it), &(*strength_it));
			buildCube(lowest - center, highest - center);

			if (params.radiiGiven)
			{
				for (std::size_t i = 0; i < points.size(); i++)
				{
					points[i].position = (*(pos_it + i)) - center;
					points[i].r = (*(strength_it + i));
					points[i].r2 = points[i].r * points[i].r;
					points[i].visitedAs = 0;
				}

				for (std::size_t i = points.size(); i < static_cast<std::size_t>(size); i++)
				{
					points.push_back(cell((*(pos_it + i)) - center, (*(strength_it + i))));
					if (points.size() > 1)
						set_bond_to_id(points.back(), points.size() - 2);
					else
						clear_bond_to(points.back());
				}
			}
			else
			{
				for (std::size_t i = 0; i < points.size(); i++)
				{
					points[i].position = (*(pos_it + i)) - center;
					points[i].r2 = (*(strength_it + i));
					points[i].r = sqrt(points[i].r2);
					points[i].visitedAs = 0;
				}
				for (std::size_t i = points.size(); i < static_cast<std::size_t>(size); i++)
				{
					points.push_back(cell((*(pos_it + i)) - center, sqrt(*(strength_it + i)), *(strength_it + i)));
					if (points.size() > 1)
						set_bond_to_id(points.back(), points.size() - 2);
					else
						clear_bond_to(points.back());
				}
			}

			buildVertices(size);
			if (params.fill_myVertices || params.fill_neighbours)
				FillAllMyVertices();
			if (params.fill_neighbours)
				FillAllNeighbours();
			if (params.fill_zeroPoints)
				FillAllZeroPoints();
			sync_all_link_mirrors();
		}

		template <class Pos_iterator, class Strength_iterator>
		inline void addMore(const Pos_iterator pos_it, const Strength_iterator strength_it, const int _newSize)
		// does standard deletion and calculation of Vertices, neighbour information, ...
		{
			// Geometry:
			//   Insert additional weighted sites incrementally into the existing tessellation.
			// Code structure:
			//   Captures a full revert snapshot, appends points, optionally rebuilds the outer cube if needed,
			//   updates generator refs after reallocations, then processes only new insertions and touched caches.
			const unsigned int gap = _newSize < points.size() ? points.size() - _newSize : 1;
			const unsigned int newSize = _newSize < points.size() ? points.size() + 1 : _newSize;
			nRevertVertices = _nVertices;
			nRevertZeros = zeros.size();
			nRevertPoints = points.size();
			std::vector<CellId> old_bond_ids(nRevertPoints, kInvalidId);
			std::vector<std::vector<CellId>> old_neighbour_ids(nRevertPoints);
			std::vector<std::array<GeneratorRef, dimension + 1>> old_vertex_generator_refs(_nVertices);
			for (std::size_t i = 0; i < nRevertPoints; ++i)
			{
				old_bond_ids[i] = points[i].bondToId;
				old_neighbour_ids[i].resize(points[i].neighboursIds.size(), kInvalidId);
				for (std::size_t j = 0; j < points[i].neighboursIds.size(); ++j)
				{
					old_neighbour_ids[i][j] = points[i].neighboursIds[j];
				}
			}
			for (unsigned int vi = 0; vi < _nVertices; ++vi)
				for (int g = 0; g <= dimension; ++g)
					old_vertex_generator_refs[vi][g] = vertices[vi].generatorRefs[g];
			for (int c = 0; c < (1 << dimension); c++)
			{
				const GeneratorRef &ref = vertices[c].generatorRefs[0];
				cornerOwners[c] = (ref.kind == GeneratorKind::point && ref.index < points.size()) ? static_cast<CellId>(ref.index) : kInvalidId;
			}
			const std::size_t old_points_capacity = points.capacity();
			points.reserve(newSize);
			const bool points_reallocated = (points.capacity() != old_points_capacity);
			if (params.radiiGiven)
			{
				const std::size_t add_count = newSize - nRevertPoints;
				for (std::size_t i = 0; i < add_count; i++)
				{
					points.push_back(cell((*(pos_it + i)) - center, (*(strength_it + i))));
					const std::size_t new_idx = points.size() - 1;
					if (new_idx >= gap)
						set_bond_to_id(points.back(), new_idx - gap);
					else
						clear_bond_to(points.back());
				}
			}
			else
			{
				const std::size_t add_count = newSize - nRevertPoints;
				for (std::size_t i = 0; i < add_count; i++)
				{
					points.push_back(cell((*(pos_it + i)) - center, sqrt(*(strength_it + i)), *(strength_it + i)));
					const std::size_t new_idx = points.size() - 1;
					if (new_idx >= gap)
						set_bond_to_id(points.back(), new_idx - gap);
					else
						clear_bond_to(points.back());
				}
			}
			{
				PDCoord lowest;
				PDCoord highest;
				PDCoord rebuild(0, 0, 0);
				getBoundingBox<PDCoord, PDFloat, PDCoord const *, PDFloat const *, dimension>(lowest, highest, newSize - nRevertPoints, &(*pos_it), &(*strength_it), 0.0);

				for (int d = 0; d < dimension; d++)
				{
					if (vertices.begin()->position[d] + center[d] - lowest[d] > rebuild[d])
						rebuild[d] = vertices.begin()->position[d] + center[d] - lowest[d];
					if (vertices[(1 << dimension) - 1].position[d] + center[d] - highest[d] < rebuild[d])
						rebuild[d] = -(vertices[(1 << dimension) - 1].position[d] + center[d] - highest[d]);
				}

				if (rebuild.squaredNorm() > 0)
				{
					clearAllmyVertices();
					for (unsigned int vi = 0; vi < _nVertices; ++vi)
					{
						vertices[vi].invalid = 0;
						vertices[vi].rrv = 0;
					}

					{ // createlike
						for (unsigned int i = 0; i < nRevertPoints; ++i)
							if (old_bond_ids[i] != kInvalidId)
								set_bond_to_id(points[i], old_bond_ids[i]);
						buildCube(vertices.begin()->position - 2 * rebuild, vertices[(1 << dimension) - 1].position + 2 * rebuild);

						buildVertices(nRevertPoints);

						if (params.fill_myVertices || params.fill_neighbours)
							FillAllMyVertices(0, (1 << dimension));
						if (params.fill_neighbours)
							FillAllNeighbours();
						if (params.fill_zeroPoints)
							FillAllZeroPoints(nRevertZeros);
					}
					nRevertVertices = _nVertices;
				}
				else
				{
					if (points_reallocated)
					{
						for (unsigned int vi = 0; vi < _nVertices; ++vi)
							for (int g = 0; g <= dimension; ++g)
							{
								const GeneratorRef &ref = old_vertex_generator_refs[vi][g];
								if (!ref.is_valid())
									continue;
								if (ref.kind == GeneratorKind::point)
								{
									if (ref.index >= points.size())
										continue;
									set_vertex_generator(vertices[vi], g, ref);
								}
								else
								{
									if (ref.index >= sideGenerators.size())
										continue;
									set_vertex_generator(vertices[vi], g, ref);
								}
							}

						for (std::size_t point_idx = 0; point_idx < nRevertPoints; ++point_idx)
						{
							cell &point = points[point_idx];
							if (point_idx < old_bond_ids.size() && old_bond_ids[point_idx] != kInvalidId)
								set_bond_to_id(point, old_bond_ids[point_idx]);
							assign_cell_neighbours_from_ids(point, old_neighbour_ids[point_idx]);
						}
					}
				}
			}

			buildVertices(newSize, nRevertPoints);

			if (params.fill_myVertices || params.fill_neighbours)
				FillAllMyVertices(nRevertPoints, nRevertVertices);
			if (params.fill_neighbours)
				FillAllNeighboursOfInvolved();
			if (params.fill_zeroPoints)
				FillAllZeroPoints(nRevertVertices, nRevertZeros);
			sync_all_link_mirrors();
		}

		void addMore(const PDCoord &pos, const PDFloat &radius, const int near)
		{
			// Single-point incremental insertion convenience wrapper.
			addMore(&pos, &radius, near);
		}

		static void make_inputfile(std::vector<PDCoord> const &position, std::vector<PDFloat> const &power)
		{
			for (unsigned int i = 0; i < position.size(); ++i)
			{
				std::println("{} {} {} {}", position[i].x(), position[i].y(), position[i].z(), power[i]);
			}
		}

		void clearAllmyVertices()
		{
			// Drop per-cell cached vertex/zero-point ownership lists.
			for (cell &point : points)
			{
				clear_cell_my_vertices(point);
				point.myZeroPoints.clear();
			}
		}
		void buildCube(const PDCoord &lowest, const PDCoord &highest)
		{
			// Geometry:
			//   Construct the artificial bounding polytope (axis-aligned cube) that clips the unbounded diagram.
			// Code structure:
			//   Creates 2^d corner vertices, assigns side generators, wires corner-corner edges by dimension bits,
			//   and normalizes local generator slot ordering for deterministic downstream matching.
			_nVertices = 1 << dimension;
			sideGenerators.clear();
			for (int i = 0; i < 2 * dimension; i++)
				sideGenerators.push_back(cell(PDCoord(0, 0, 0), 0));
			PDCoord lhc = lowest;
			vertices[0].setTo(lowest);
			for (int j = dimension - 1; j >= 0; j--)
				set_vertex_generator(vertices[0], j + 1, GeneratorRef(GeneratorKind::side, static_cast<std::size_t>(j)));
			for (int i = 0; i < (1 << dimension); i++)
				vertices[i].rrv = 0;
			for (int i = 0; i < (1 << dimension); i++)
				vertices[i].invalid = 0;
			for (int i = 1; i < (1 << dimension); i++)
			{
				int j = 0;
				while (lhc[j] == highest[j])
				{
					lhc[j] = lowest[j];
					j++;
				}
				lhc[j] = highest[j];
				vertices[i].setTo(lhc);
				for (j = dimension - 1; j >= 0; j--)
				{
					const std::size_t side_idx = static_cast<std::size_t>((lhc[j] == lowest[j]) ? j : (j + dimension));
					set_vertex_generator(vertices[i], j + 1, GeneratorRef(GeneratorKind::side, side_idx));
				}
			}
			for (int i = 0; i < (1 << dimension); i++)
				for (int d = 0; d < dimension; d++)
				{
					const int ii = i;
					const int j = (ii >> d) % 2 ? ii - (1 << d) : ii + (1 << d);
					set_vertex_endpoint(vertices[i], d + 1, static_cast<VertexId>(j));
				}
			for (int i = 0; i < (1 << dimension); i++)
			{ //: TODO: direct sort or faster?
				for (int g = dimension - 1; g > 0; g--)
					for (int j = g; j > 0; j--)
					{
						// Adjacent elements in the same array always have increasing addresses.
						swap_vertex_link_slots(vertices[i], j, j + 1);
					}
			}
			for (int d = 1; d <= dimension; d++) // n
			{
				const GeneratorRef &ref0 = vertices[0].generatorRefs[d];
				if (ref0.kind == GeneratorKind::point && ref0.index < points.size())
					push_cell_my_vertex(points[ref0.index], static_cast<VertexId>(0));
				const GeneratorRef &ref1 = vertices[(1 << dimension) - 1].generatorRefs[d];
				if (ref1.kind == GeneratorKind::point && ref1.index < points.size())
					push_cell_my_vertex(points[ref1.index], static_cast<VertexId>(0));
			}
		}
		void buildVertices(const unsigned int &nPoints, const int from = 0)
		{
			// Geometry:
			//   Perform incremental regular triangulation/power-cell update by inserting points one by one.
			// Code structure:
			//   For each insertion: prepare replacement frontier, attempt local rebuild, fallback on numerical
			//   dampening when ambiguous, then run consistency checks and optional cache recomputation.
			//	try
			{

				if (points.size() > 0)
				{
					maxr2 = points[0].r2;
					for (int i = from; i < static_cast<int>(nPoints); i++)
					{
						if (maxr2 < points[i].r2)
							maxr2 = points[i].r2;
					}
					powerErr = 1000 * error(maxr2);
					if (from == 0)
					{
						insertFirst();
					}

					for (unsigned int i = (from == 0) ? 1 : from; i < nPoints; i++)
					{
						unsigned int done = 1;
						while (1)
						{
							if (doInsertion(prepareInsertion(i)))
								break;
							const PDFloat errorScale = insertionErrorScale;
							if (__power_diagram_internal_timing__)
							{
								t3 += clock();
								t4 += clock();
							}

							CellId identicalPointId = kInvalidId;
							done++;
							if (done > 100)
								throw MyException();
							cell &insertion_cell = points[i];
							// is this an Identical point problem?
							{
								CellId closest_id = kInvalidId;
								PDFloat mindist = 0;
								for (std::size_t involved_idx = 1; involved_idx < involved_size(); ++involved_idx)
								{
									const CellId candidate_id = involved_id_at(involved_idx);
									if (candidate_id == kInvalidId)
										continue;
									const cell &candidate = cell_at(candidate_id);
									const PDFloat dist = (candidate.position - insertion_cell.position).squaredNorm();
									if (closest_id == kInvalidId || dist < mindist)
									{
										mindist = dist;
										closest_id = candidate_id;
									}
								}
								if (closest_id != kInvalidId && error(insertion_cell.r) > sqrt(mindist))
								{
									identicalPointId = closest_id;
									if (params.with_warnings)
									{
										std::println("numerical similar point to {} found. {} is ignored", closest_id + 1, i + 1);
									}
								}
							}
							// delete new vertices built directly into vertices (when unused has been empty)
							if (_nUnused == 0)
							{
								// here comes the deletion ;)
								_nVertices -= insertion_cell.myVerticesIds.size() - unused.size() + _nUnused;
								for (const VertexId involved_vid : insertion_cell.myVerticesIds)
								{
									if (involved_vid == kInvalidId || involved_vid >= vertices.size())
										continue;
									vertex &involved_vertex = vertex_at(involved_vid);
									if (involved_vid < (1u << dimension))
									{
										_nVertices++;
									}
									else
										involved_vertex.disconnect();
								}
							}

							// delete new vertices built on unused (not needed because endpoints not constructed,yet)
							for (std::size_t unused_idx = _nUnused; unused_idx < unused.size(); ++unused_idx)
							{
								const VertexId unused_id = unused[unused_idx];
								if (unused_id == kInvalidId || unused_id >= vertices.size())
									continue;
								vertex_at(unused_id).disconnect();
							}
							_nUnused = unused.size();

							// reconnect replaced with persisting
							for (std::size_t replaced_idx = 0; replaced_idx < replaced_size(); ++replaced_idx)
							{
								const VertexId replaced_id = replaced_id_at(replaced_idx);
								if (replaced_id == kInvalidId || replaced_id >= vertices.size())
									continue;
								vertex &replaced_vertex = vertex_at(replaced_id);
								for (int endpoint_idx = replaced_vertex.isCorner(); endpoint_idx <= dimension; ++endpoint_idx)
								{
									VertexId endpoint_id = replaced_vertex.resolved_endpoint_id(*this, endpoint_idx);
									if (endpoint_id == kInvalidId || endpoint_id >= vertices.size())
										continue;
									if (vertex_at(endpoint_id).rrv <= 0)
									{
										vertex &endpoint = vertex_at(endpoint_id);
										endpoint.rrv = 0;
										for (int g1 = replaced_vertex.isCorner(); g1 <= dimension; g1++)
											for (int g2 = endpoint.isCorner(); g2 <= dimension; g2++)
											{
												if (replaced_vertex.generatorRefs[nth(0, g1)].kind == endpoint.generatorRefs[nth(0, g2)].kind && replaced_vertex.generatorRefs[nth(0, g1)].index == endpoint.generatorRefs[nth(0, g2)].index && replaced_vertex.generatorRefs[nth(1, g1)].kind == endpoint.generatorRefs[nth(1, g2)].kind && replaced_vertex.generatorRefs[nth(1, g1)].index == endpoint.generatorRefs[nth(1, g2)].index && replaced_vertex.generatorRefs[nth(2, g1)].kind == endpoint.generatorRefs[nth(2, g2)].kind && replaced_vertex.generatorRefs[nth(2, g1)].index == endpoint.generatorRefs[nth(2, g2)].index)
												{
													set_vertex_endpoint_deferred(replaced_vertex, g1, endpoint_id);
													set_vertex_endpoint_deferred(endpoint, g2, replaced_id);
												}
											}
									}
								}
							}

							const VertexId fallback_replaced_id = replaced_id_at(0);
							for (std::size_t involved_idx = 1; involved_idx < involved_size(); ++involved_idx)
							{
								const CellId involved_id = involved_id_at(involved_idx);
								if (involved_id == kInvalidId)
									continue;
								cell &involved_cell = cell_at(involved_id);
								if (involved_cell.myVerticesIds.empty())
									continue;
								const VertexId representative_id = involved_cell.myVerticesIds[0];
								if (representative_id == kInvalidId || representative_id >= vertices.size())
									continue;
								if (!vertex_at(representative_id).isConnected() && fallback_replaced_id != kInvalidId && fallback_replaced_id < vertices.size())
									set_cell_my_vertex(involved_cell, 0, fallback_replaced_id);
							}

							// set everything zero again
							SetInvolvedPersistingVisitedToZero();
							clear_cell_my_vertices(insertion_cell);
							for (std::size_t replaced_idx = 0; replaced_idx < replaced_size(); ++replaced_idx)
							{
								const VertexId replaced_id = replaced_id_at(replaced_idx);
								if (replaced_id == kInvalidId || replaced_id >= vertices.size())
									continue;
								vertex &replaced_vertex = vertex_at(replaced_id);
								replaced_vertex.rrv = 0;
								replaced_vertex.invalid = 0;
							}

							clear_replaced();
							if (identicalPointId != kInvalidId)
							{
								cell &identicalPoint = cell_at(identicalPointId);
								if (!identicalPoint.myVerticesIds.empty())
								{
									const VertexId identical_vid = identicalPoint.myVerticesIds.front();
									if (identical_vid != kInvalidId && identical_vid < vertices.size())
										push_cell_my_vertex(insertion_cell, identical_vid);
								}
								break;
							}
							else
							{
								const PDFloat oldr2 = insertion_cell.r2;
								if (params.with_warnings)
									std::print(" Numerical Zero Warning: Power of {} is reduced from {}", i + 1, insertion_cell.r2);
								insertion_cell.r2 -= pow(2.0, done) * (errorScale);
								if (insertion_cell.r2 >= 0)
									insertion_cell.r = sqrt(insertion_cell.r2);
								else
									insertion_cell.r = -sqrt(-insertion_cell.r2);
								if (params.with_warnings)
									std::println(" to {} ( Change was {} )", insertion_cell.r2, insertion_cell.r2 - oldr2);
							}

							if (done > 100)
							{
								std::println("Exception : cannot get stable results");
								throw MyException();
							}
						}
					}

					if (!params.without_check)
					{
						PDFloat checkconst = 0;
						//				std::println("checking diagram");
						for (unsigned int vi = 0; vi < _nVertices; ++vi)
							if (vertices[vi].isConnected())
								for (std::size_t point_id = 0; point_id < points.size(); ++point_id)
								{
									cell &point = points[point_id];
									if (point.power(vertices[vi].position) - vertices[vi].powerValue < -checkconst)
									{
										const auto &v = vertices[vi];
										const auto is_point_ref = [point_id](const GeneratorRef &ref) -> bool
										{
											return ref.kind == GeneratorKind::point && ref.index == point_id;
										};
										if (!is_point_ref(v.generatorRefs[0]) && !is_point_ref(v.generatorRefs[1]) && !is_point_ref(v.generatorRefs[2]) && !is_point_ref(v.generatorRefs[3]))
										{
											checkconst = point.power(vertices[vi].position) - vertices[vi].powerValue;

											cell &g0 = cell_from_ref(vertices[vi].generatorRefs[0]);
											cell &g1 = cell_from_ref(vertices[vi].generatorRefs[1]);
											VertexId e0_id = vertices[vi].resolved_endpoint_id(*this, 0);
											VertexId e1_id = vertices[vi].resolved_endpoint_id(*this, 1);
											if (e0_id == kInvalidId || e0_id >= vertices.size() || e1_id == kInvalidId || e1_id >= vertices.size())
												throw MyException();
											const vertex &e0 = vertex_at(e0_id);
											const vertex &e1 = vertex_at(e1_id);
											std::println("totaly wrong are {} {} {} {}          {} {} {}", point_id, point.power(vertices[vi].position), g0.power(vertices[vi].position), g1.power(vertices[vi].position), g0.position[0], g0.position[1], g0.position[2]);
											std::println("{} {} {} {}",
														 ((vertices[vi].generatorRefs[0].kind == GeneratorKind::point) ? static_cast<long long>(vertices[vi].generatorRefs[0].index) : -1),
														 ((vertices[vi].generatorRefs[1].kind == GeneratorKind::point) ? static_cast<long long>(vertices[vi].generatorRefs[1].index) : -1),
														 ((vertices[vi].generatorRefs[2].kind == GeneratorKind::point) ? static_cast<long long>(vertices[vi].generatorRefs[2].index) : -1),
														 ((vertices[vi].generatorRefs[3].kind == GeneratorKind::point) ? static_cast<long long>(vertices[vi].generatorRefs[3].index) : -1));
											std::println("{} {}", point.power(e0.position), point.power(e1.position));
											std::println("{} {}", g0.power(e0.position), g0.power(e1.position));
											std::println("{} {}", g1.power(e0.position), g1.power(e1.position));
											std::println("{}\n", pd_to_string(vertices[vi].position));
											std::println("{}\n", pd_to_string(e0.position));
											std::println("{}", pd_to_string(e1.position));
											throw MyException();
										}
									}
								}

						if (std::abs(checkconst) > 0.001)
							std::println("the error of the worst vertex is around {}", checkconst);
						// dump_vertices();
					}
				}
				else
				{ /*no vertices*/
				}
			}

			// delete waste that was produced, but during cleanup, never move waste (unvalidating "unused" pointer)
			std::sort(unused.begin(), unused.end());
			for (int i = static_cast<int>(unused.size()) - 1; i >= 0; --i)
			{
				const VertexId unused_id = unused[i];
				if (unused_id != (_nVertices - 1))
				{
					if (unused_id == kInvalidId || unused_id >= vertices.size())
						throw MyException();
					const VertexId source_id = _nVertices - 1;
					--_nVertices;
					vertices[_nVertices].moveAddressNetworkUpdateOnly(*this, vertex_at(unused_id), source_id, unused_id);
				}
				else
				{
					--_nVertices;
				}
			}

			unused.clear();
			_nUnused = 0;
		}
		inline void findReplacedVertex(VertexId &this_id, PDFloat &value, const cell &insertionPoint)
		{
			// Geometry:
			//   Gradient-like walk on power values to find a seed vertex likely inside the replaced region.
			// Code structure:
			//   Uses 1-hop/2-hop/3-hop neighborhood scans with tie-handling; falls back to full scan only under
			//   near-degenerate numeric conditions.
			if (value < 0)
				return;
			vertex &This = vertex_at(this_id);
			PDFloat newValue;
			PDFloat smallVal = std::numeric_limits<PDFloat>::max(); // something larger than value
			// look at the neighbours
			for (int idx = This.isCorner(); idx <= dimension; ++idx)
			{
				// each powerdiff value defines a plane between insertionPoint and current cell (generator0) approach the direction perpendicular to that plane in direction of insertion point !
				VertexId endpoint_id = This.resolved_endpoint_id(*this, idx);
				if (endpoint_id == kInvalidId || endpoint_id >= vertices.size())
					continue;
				vertex &endpoint = vertex_at(endpoint_id);
				newValue = endpoint.powerdiff3D(cell_from_ref(endpoint.generatorRefs[0]), insertionPoint);
				if (newValue < value)
				{
					value = newValue;
					this_id = endpoint_id;
					if (value < 0)
						return;
					idx = This.isCorner() - 1;
				}
				else if (newValue == value)
					smallVal = newValue;
			}
			if ((smallVal != value))
				return;
			smallVal = std::numeric_limits<PDFloat>::max();
			// found value is not definitely the best one
			// try hard to be sure not beeing in a local minimum (second neighbour)
			for (int g = This.isCorner(); g <= dimension; ++g)
			{
				VertexId ep1_id = This.resolved_endpoint_id(*this, g);
				if (ep1_id == kInvalidId || ep1_id >= vertices.size())
					continue;
				vertex &ep1 = vertex_at(ep1_id);
				for (int g2 = ep1.isCorner(); g2 <= dimension; ++g2)
				{
					VertexId candidate_id = ep1.resolved_endpoint_id(*this, g2);
					if (candidate_id == kInvalidId || candidate_id >= vertices.size() || candidate_id == this_id)
						continue;
					vertex &candidate = vertex_at(candidate_id);
					newValue = candidate.powerdiff3D(cell_from_ref(candidate.generatorRefs[0]), insertionPoint);
					if (newValue < value)
					{
						value = newValue;
						this_id = candidate_id;
						return findReplacedVertex(this_id, value, insertionPoint);
					}
					else if (newValue == value)
						smallVal = newValue;
				}
			}
			if ((smallVal != value))
				return;
			smallVal = std::numeric_limits<PDFloat>::max();
			// second was also close... third neighbour...
			for (int g = This.isCorner(); g <= dimension; ++g)
			{
				VertexId ep1_id = This.resolved_endpoint_id(*this, g);
				if (ep1_id == kInvalidId || ep1_id >= vertices.size())
					continue;
				vertex &ep1 = vertex_at(ep1_id);
				for (int g2 = ep1.isCorner(); g2 <= dimension; ++g2)
				{
					VertexId ep2_id = ep1.resolved_endpoint_id(*this, g2);
					if (ep2_id == kInvalidId || ep2_id >= vertices.size() || ep2_id == this_id)
						continue;
					vertex &ep2 = vertex_at(ep2_id);
					for (int g3 = ep2.isCorner(); g3 <= dimension; ++g3)
					{
						VertexId candidate_id = ep2.resolved_endpoint_id(*this, g3);
						if (candidate_id == kInvalidId || candidate_id >= vertices.size() || candidate_id == ep1_id || candidate_id == this_id)
							continue;
						vertex &candidate = vertex_at(candidate_id);
						newValue = candidate.powerdiff3D(cell_from_ref(candidate.generatorRefs[0]), insertionPoint);
						if (newValue < value)
						{
							value = candidate.powerdiff3D(cell_from_ref(candidate.generatorRefs[0]), insertionPoint);
							this_id = candidate_id;
							return findReplacedVertex(this_id, value, insertionPoint);
						}
						else if (newValue == value)
							smallVal = newValue;
					}
				}
			}
			if (smallVal != value)
				return;
			if (params.with_warnings)
				std::println("warning : program slowed down because of too small accuracy");
			//...so the numerical problem wants to be tough? A fat lot we care!
			for (unsigned int vi = 0; vi < nVertices(); ++vi)
				if (vertices[vi].isConnected())
				{
					if (vertices[vi].powerdiff3D(cell_from_ref(vertices[vi].generatorRefs[0]), insertionPoint) < value)
					{
						value = vertices[vi].powerdiff3D(cell_from_ref(vertices[vi].generatorRefs[0]), insertionPoint);
						this_id = vi;
					}
				}
		}
		inline ReplaceState finiteReplaced(vertex &This, const CellId cell_id)
		{
			if (cell_id == kInvalidId)
				return ReplaceState::ambiguous;
			// Geometry:
			//   Compare power of candidate insertion site against current owner at this vertex.
			// Code structure:
			//   Writes signed replacement margin into rrv and classifies using robust tolerance thresholds.
			This.rrv = This.powerdiff3D(cell_at(cell_id), cell_from_ref(This.generatorRefs[0]));
			if (above_power_err(This.rrv))
				return ReplaceState::replaced;
			if (below_neg_power_err(This.rrv))
				return ReplaceState::persisting;
			This.rrv = 0;
			return ReplaceState::ambiguous;
		}
		void dump_vertices(std::ostream &out = std::cout)
		{
			// Geometry:
			//   Diagnostic dump of current local combinatorics and vertex power consistency.
			// Code structure:
			//   Emits one line per vertex plus neighbor endpoint lines using GeneratorRef-derived IDs.
			out << "vertices, generators and neighbours " << std::endl;
			const auto generator_idx = [this](const GeneratorRef &ref) -> long long
			{
				if (ref.kind != GeneratorKind::point || ref.index >= points.size())
					return -1;
				const CellId id = static_cast<CellId>(ref.index);
				return (id == kInvalidId) ? -1 : static_cast<long long>(id);
			};

			for (unsigned int vi = 0; vi < _nVertices; ++vi)
			{
				vertex &v = vertices[vi];
				if (v.isConnected() && (!v.isCorner()))
				{
					out << generator_idx(v.generatorRefs[0]) << " " << std::flush;
					out << generator_idx(v.generatorRefs[1]) << " " << std::flush;
					out << generator_idx(v.generatorRefs[2]) << " " << std::flush;
					out << generator_idx(v.generatorRefs[3]) << "    " << std::flush;
					out << v.position[0] << " " << std::flush;
					out << v.position[1] << " " << std::flush;
					out << v.position[2] << " " << std::flush;
					out << v.powerValue << "   " << std::flush;
					out << cell_from_ref(v.generatorRefs[0]).power(v.position) << std::endl;

					for (int eg = 0; eg <= dimension; ++eg)
					{
						VertexId endpoint_id = v.resolved_endpoint_id(*this, eg);
						if (endpoint_id == kInvalidId || endpoint_id >= vertices.size())
							continue;
						const vertex &endpoint = vertex_at(endpoint_id);
						out << " " << generator_idx(endpoint.generatorRefs[0]) << " " << generator_idx(endpoint.generatorRefs[1]) << " " << generator_idx(endpoint.generatorRefs[2]) << " " << generator_idx(endpoint.generatorRefs[3]) << std::endl;
					}
				}
				else if (v.isCorner())
				{
					out << generator_idx(v.generatorRefs[0]) << " " << generator_idx(v.generatorRefs[1]) << " " << generator_idx(v.generatorRefs[2]) << " " << generator_idx(v.generatorRefs[3]) << "    " << v.position[0] << " " << v.position[1] << " " << v.position[2] << " " << v.powerValue << "   " << cell_from_ref(v.generatorRefs[0]).power(v.position) << std::endl;
				}
				else
				{
					out << "outtake" << std::endl;
				}
				out << std::endl;
			}
		}

		//	inline const int dimension()const{return dimension;}
		inline int nPoints() const { return points.size(); }
		inline unsigned int const &nVertices() const { return _nVertices; }

		bool hasVirtualGenerators(const vertex &that) const
		{
			GeneratorRef ref = that.resolved_generator_ref(*this, dimension);
			return !ref_is_real_point(ref);
		}
		int nVirtualGenerators(const vertex &that) const
		{
			GeneratorRef ref = that.resolved_generator_ref(*this, dimension);
			const bool dim_is_real = ref_is_real_point(ref);
			if (dim_is_real)
				return 0;
			return that.isCorner() ? 3 : 2;
		}
		CellId findCellInsideCube(const PDCoord &pos, CellId hint_id = kInvalidId)
		{
			// Greedy neighbor walk to find the cell with minimal power at pos inside the current cube.
			if (points.empty())
				return kInvalidId;
			if (hint_id == kInvalidId)
				hint_id = points.size() / 2;
			const cell &hint_cell = points[hint_id];
			for (const CellId neighbour_id : hint_cell.neighboursIds)
			{
				if (neighbour_id == kInvalidId || neighbour_id >= points.size())
					continue;
				const cell &neighbour = points[neighbour_id];
				if (neighbour.power(pos) < hint_cell.power(pos))
					return findCellInsideCube(pos, neighbour_id);
			}
			return hint_id;
		}

	private:
		VertexId getRepresentative(const CellId start_id)
		{
			// Retrieve a connected representative vertex following bondToId chain by IDs.
			CellId current_id = start_id;
			while (current_id != kInvalidId && current_id < points.size())
			{
				cell &current = cell_at(current_id);
				for (const VertexId vid : current.myVerticesIds)
				{
					if (vid == kInvalidId || vid >= vertices.size())
						continue;
					vertex &candidate = vertex_at(vid);
					if (!candidate.isConnected())
						continue;
					for (int g = 0; g <= dimension; ++g)
					{
						const GeneratorRef &ref = candidate.generatorRefs[g];
						if (ref.kind == GeneratorKind::point && ref.index == current_id)
							return vid;
					}
				}
				if (current.bondToId == current_id)
					break;
				current_id = current.bondToId;
			}
			return 0;
		}
		VertexId prepareInsertion(const CellId this_id, VertexId hint_id = kInvalidId)
		{
			// Geometry:
			//   Compute the local conflict region for a new site and stabilize it if near degeneracy appears.
			// Code structure:
			//   Picks a representative hint, runs replacement flood-fill, and if ambiguous repeatedly perturbs
			//   radius/power until a stable replaced/persisting partition is obtained or throws.
			if (this_id == kInvalidId || this_id >= points.size())
				throw MyException();
			cell &This = cell_at(this_id);
			if (__power_diagram_internal_timing__)
				t2 -= clock();
			// there is a power of new cell that is so low, that only one vertex would be replaced. *hint will be the one
			hint_id = getRepresentative(This.bondToId);
			if (hint_id == kInvalidId || hint_id >= vertices.size())
				hint_id = 0;
			PDFloat value = vertex_at(hint_id).powerdiff3D(cell_from_ref(vertex_at(hint_id).generatorRefs[0]), This);
			findReplacedVertex(hint_id, value, This);
			if (__power_diagram_internal_timing__)
				t2 += clock();

			unsigned int done = 1;
			while (1)
			{
				if (done != 1)
				{
					value = vertex_at(hint_id).powerdiff3D(cell_from_ref(vertex_at(hint_id).generatorRefs[0]), This);
					findReplacedVertex(hint_id, value, This);
				}

				if (FillReplacedPersistingAndInvolved(this_id, hint_id))
					break;

				const PDFloat oldr2 = This.r2;
				if (params.with_warnings)
					std::print("Numerical Warning: Power of {} is reduced from {}", this_id + 1, This.r2);
				SetInvolvedPersistingVisitedToZero();
				clear_cell_my_vertices(This);
				for (std::size_t replaced_idx = 0; replaced_idx < replaced_size(); ++replaced_idx)
				{
					const VertexId replaced_id = replaced_id_at(replaced_idx);
					if (replaced_id == kInvalidId || replaced_id >= vertices.size())
						continue;
					vertex &replaced_vertex = vertex_at(replaced_id);
					replaced_vertex.rrv = 0;
					for (int g = replaced_vertex.isCorner(); g <= dimension; g++)
					{
						VertexId endpoint_id = replaced_vertex.resolved_endpoint_id(*this, g);
						if (endpoint_id != kInvalidId && endpoint_id < vertices.size())
							vertex_at(endpoint_id).rrv = 0;
					}
				}
				clear_replaced();
				This.r2 -= pow(2.0, done) * (PowerDiagram<PDFloat, PDCoord, dimension>::powerErr);
				if (This.r2 > 0)
					This.r = sqrt(This.r2);
				else
					This.r = -sqrt(-This.r2);
				if (params.with_warnings)
					std::println(" to {} ( Change was {} )", This.r2, This.r2 - oldr2);
				done++;
				if (done > 100)
				{
					std::println("exception : cannot get stable results with atom {} {}", this_id, pd_to_string(This.position + center));
					throw MyException();
				}
			}
			return hint_id;
		}

		bool doInsertion(const VertexId hint_id)
		{
			// Geometry:
			//   Materialize the topological transition after the conflict region is known.
			// Code structure:
			//   Creates all required new finite vertices on cut edges, connects them, compacts invalid vertices,
			//   then refreshes representative ownership and visited flags.
			//			if(hint exists)
			{
				if (__power_diagram_internal_timing__)
				{
					const unsigned int zeit = clock();
					t3 -= zeit;
					t4 -= zeit;
				}
				if (!CreateFiniteVerticesFromReplaced())
					return false;
				if (__power_diagram_internal_timing__)
				{
					const unsigned int zeit = clock();
					t4 += zeit;
					t5 -= zeit;
				}
				ConnectNewFinitesAmongThemselves3D();
				if (__power_diagram_internal_timing__)
					t5 += clock();
				UpdateUnused();
				AssignRepresentativeVerticesToCells(hint_id);
				SetInvolvedPersistingVisitedToZero();
				validate_phase_mirror_invariants();
				if (__power_diagram_internal_timing__)
					t3 += clock();
			}
			return true;
		}

		void clear_interna()
		{
			// Reset per-insertion transient containers.
			clear_replaced();
			clear_involved();
			//		planes.clear();//should always be clean
		}

		inline void insertFirst()
		{
			// Geometry:
			//   Initialize first finite cell against the clipping cube.
			// Code structure:
			//   Sets corner power values from the first site and marks all corners as owned by that site.
			clear_interna();
			for (int i = 0; i < (1 << dimension); i++)
				vertices[i].setPowerData(points[0]);
			for (int i = 0; i < (1 << dimension); i++)
				set_vertex_generator(vertices[i], 0, GeneratorRef(GeneratorKind::point, 0));
			for (int i = 0; i < (1 << dimension); i++)
				push_cell_my_vertex(points[0], static_cast<VertexId>(i));
		}
		void FillAllMyVertices(const int fromPoint = 0, const int fromVertex = 1 << dimension)
		{
			// Geometry:
			//   Build inverse incidence: for each real cell, list all finite boundary vertices.
			// Code structure:
			//   Scans vertices, filters corners/virtual generators based on mode, and populates myVerticesIds;
			//   optionally marks involved cells for incremental neighbour refresh.
			{
				if (fromPoint > 0)
					clear_involved();
				for (std::size_t point_idx = static_cast<std::size_t>(fromPoint); point_idx < points.size(); ++point_idx)
					clear_cell_my_vertices(points[point_idx]);

				for (unsigned int vi = fromVertex; vi < _nVertices; ++vi)
					if (!(vertices[vi].invalid))
					{
						if (!(hasVirtualGenerators(vertices[vi])))
							for (int g = 0; g <= dimension; ++g)
							{
								GeneratorRef ref = vertices[vi].resolved_generator_ref(*this, g);
								if (ref.kind != GeneratorKind::point || ref.index >= points.size())
									continue;
								cell &generator = points[ref.index];
								push_cell_my_vertex(generator, static_cast<VertexId>(vi));
								if (fromPoint > 0)
									if (generator.visitedAs == 0)
									{
										generator.visitedAs = -1;
										push_involved(GeneratorRef(GeneratorKind::point, ref.index));
									}
							}
						else if (!vertices[vi].isCorner())
							for (int g = 0; g <= dimension; ++g)
							{
								GeneratorRef ref = vertices[vi].resolved_generator_ref(*this, g);
								if (ref.kind == GeneratorKind::side)
									break;
								if (ref.kind != GeneratorKind::point || ref.index >= points.size())
									continue;
								cell &generator = points[ref.index];
								push_cell_my_vertex(generator, static_cast<VertexId>(vi));
								if (fromPoint > 0)
									if (generator.visitedAs == 0)
									{
										generator.visitedAs = -1;
										push_involved(GeneratorRef(GeneratorKind::point, ref.index));
									}
							}
						else
						{ /*dont give corners to sasa code, it doesnt check for it... so we define myVertices as not holding corners!*/
							std::println("wrong internal order, SASA stopped");
							throw MyException();
						}
					}
			}
		}

		void FillAllNeighbours()
		{
			// Geometry:
			//   Derive cell-cell adjacency graph from shared finite vertices.
			// Code structure:
			//   Clears neighbour lists, then for each cell traverses incident vertices and collects distinct
			//   point generators using visitedAs as a temporary duplicate filter.
			for (cell &point : points)
			{
				clear_cell_neighbours(point);
				point.visitedAs = -1;
			}
			for (std::size_t point_idx = 0; point_idx < points.size(); ++point_idx)
			{
				cell &point = points[point_idx];
				const int current_cell_order = static_cast<int>(point_idx);
				for (const VertexId vid : point.myVerticesIds)
				{
					if (vid == kInvalidId || vid >= vertices.size())
						continue;
					vertex &vtx = vertex_at(vid);
					if (vtx.isCorner())
						continue;
					for (int g = dimension; g >= 0; g--)
					{
						GeneratorRef ref = vtx.resolved_generator_ref(*this, g);
						if (ref.kind == GeneratorKind::point && ref.index < points.size())
						{
							cell &neighbour = points[ref.index];
							if (neighbour.visitedAs < current_cell_order && ref.index != point_idx)
							{
								push_cell_neighbour(point, static_cast<CellId>(ref.index));
								neighbour.visitedAs = current_cell_order;
							}
						}
					}
				}
			}

			for (cell &point : points)
				point.visitedAs = 0;
		}

		void FillAllNeighboursOfInvolved()
		{
			// Geometry:
			//   Incremental adjacency update restricted to recently touched cells.
			// Code structure:
			//   Removes stale neighbours, adds new ones from local vertex incidence, and finally clears
			//   temporary visited markers on the involved subset.
			sort_involved_by_ref();
			for (std::size_t involved_idx = 0; involved_idx < involved_size(); ++involved_idx)
			{
				const CellId involved_id = involved_id_at(involved_idx);
				if (involved_id == kInvalidId)
					continue;
				cell &involved_cell = cell_at(involved_id);
				std::size_t neighbour_idx = 0;
				while (neighbour_idx < involved_cell.neighboursIds.size())
				{
					const CellId neighbour_id = involved_cell.neighboursIds[neighbour_idx];
					if (neighbour_id == kInvalidId || neighbour_id >= points.size() || points[neighbour_id].visitedAs == -1)
						erase_cell_neighbour(involved_cell, neighbour_idx);
					else
						++neighbour_idx;
				}
			}
			for (std::size_t involved_idx = 0; involved_idx < involved_size(); ++involved_idx)
			{
				const CellId current_cell_id = involved_id_at(involved_idx);
				if (current_cell_id == kInvalidId)
					continue;
				cell &involved_cell = cell_at(current_cell_id);
				const int current_cell_order = static_cast<int>(current_cell_id);
				for (const VertexId vid : involved_cell.myVerticesIds)
				{
					if (vid == kInvalidId || vid >= vertices.size())
						continue;
					vertex &vtx = vertex_at(vid);
					for (int g = dimension; g >= 0; g--)
					{
						GeneratorRef ref = vtx.resolved_generator_ref(*this, g);
						if (ref.kind == GeneratorKind::point && ref.index < points.size())
						{
							cell &neighbour = points[ref.index];
							if (neighbour.visitedAs != 0 && neighbour.visitedAs <= current_cell_order && ref.index != current_cell_id)
							{
								push_cell_neighbour(involved_cell, static_cast<CellId>(ref.index));
								neighbour.visitedAs = current_cell_order + 1;
							}
						}
					}
				}
			}

			for (std::size_t involved_idx = 0; involved_idx < involved_size(); ++involved_idx)
			{
				const CellId involved_id = involved_id_at(involved_idx);
				if (involved_id == kInvalidId)
					continue;
				cell &involved_cell = cell_at(involved_id);
				involved_cell.visitedAs = 0;
			}
		}
		void FillAllZeroPoints(unsigned int fromVertex = (1 << dimension), const unsigned int fromZero = 0)
		{
			// Geometry:
			//   Locate all edge parameters where interpolated power crosses zero (surface intersection seeds).
			// Code structure:
			//   Iterates eligible edges, solves quadratic/degenerate crossing equations, stores zeroPoint records,
			//   then back-links them to real cells via myZeroPoints.
			zeros.erase(zeros.begin() + fromZero, zeros.end());
			for (unsigned int vertex_index = fromVertex; vertex_index < this->_nVertices; ++vertex_index)
			{
				const vertex &current = vertices[vertex_index];
				if (!(current.invalid))
				{
					GeneratorRef boundary_ref = current.resolved_generator_ref(*this, dimension - 1);
					if (boundary_ref.kind == GeneratorKind::point && boundary_ref.index < points.size())
						for (int endpoint_idx = (hasVirtualGenerators(current)) * 3; endpoint_idx <= dimension; ++endpoint_idx)
						{
							VertexId endpoint_id = current.resolved_endpoint_id(*this, endpoint_idx);
							if (endpoint_id != kInvalidId && endpoint_id > vertex_index)
							{
								const vertex &endpoint = vertex_at(endpoint_id);
								if (current.powerValue > 0)
								{
									const int branch = endpoint_idx;
									const PDFloat &v3 = endpoint.powerValue;
									const PDFloat &v2 = current.powerValue;
									GeneratorRef ref = current.resolved_generator_ref(*this, branch == 0);
									const PDFloat &v1 = cell_from_ref(ref).power(2 * current.position - endpoint.position);
									const PDFloat quot = 2 * (v1 + v3 - 2 * v2);
									// const PDFloat rootsq=sqr(v1-v3)-4*quot*v2;
									const PDFloat rootsq = (v1 - v3) * (v1 - v3) - 4 * quot * v2;
									if (rootsq <= 0)
										continue;
									if (below_power_err(quot))
										if (v1 >= 0 && v2 >= 0 && v3 >= 0)
											continue;
									const PDFloat rootquot = sqrt(rootsq) / quot;
									const PDFloat min = (v1 - v3) / quot;
									const PDFloat sol1 = min + rootquot;
									const PDFloat sol2 = min - rootquot;
									if (sol1 > 0 && sol1 < 1)
										if (endpoint.powerValue > 0)
										{
											push_zero_from_edge(vertex_index, branch, sol1);
											push_zero_from_edge(vertex_index, branch, sol2);
										}
										else
											push_zero_from_edge(vertex_index, branch, sol1);
									else if (sol2 > 0 && sol2 < 1)
										push_zero_from_edge(vertex_index, branch, sol2);
									else
									{ // the covered zeros
										push_zero_from_edge(vertex_index, branch, sol1);
										push_zero_from_edge(vertex_index, branch, sol2);
									}
								}
								else if (endpoint.powerValue > 0)
								{
									const int branch = endpoint_idx;
									const PDFloat &v3 = endpoint.powerValue;
									const PDFloat &v2 = current.powerValue;
									GeneratorRef ref = current.resolved_generator_ref(*this, branch == 0);
									const PDFloat &v1 = cell_from_ref(ref).power(2 * current.position - endpoint.position);
									const PDFloat quot = 2 * (v1 + v3 - 2 * v2);
									// const PDFloat rootsq=sqr(v1-v3)-4*quot*v2;
									const PDFloat rootsq = (v1 - v3) * (v1 - v3) - 4 * quot * v2;
									if (rootsq <= 0)
										continue;
									if (below_power_err(quot))
										if (v1 >= 0 && v2 >= 0 && v3 >= 0)
											continue;
									const PDFloat rootquot = sqrt(rootsq) / quot;
									const PDFloat min = (v1 - v3) / quot;
									const PDFloat sol1 = min + rootquot;
									const PDFloat sol2 = min - rootquot;
									if (sol1 > 0 && sol1 < 1)
										push_zero_from_edge(vertex_index, branch, sol1);
									else
										push_zero_from_edge(vertex_index, branch, sol2);
								}
							}
						}
				}
				/*
				else //vertex on edge of cube dont have connections through spheres
				{

				}
				*/
			}
			for (unsigned int i = fromZero; i < zeros.size(); i++)
			{
				for (int g = 0; g < dimension; ++g)
				{
					const GeneratorRef &ref = zeros[i].generatorRefs[g];
					if (ref.is_valid() && ref.kind == GeneratorKind::point && ref.index < points.size())
					{
						points[ref.index].myZeroPoints.push_back(i);
					}
				}
			}
		}

		inline bool tryToBuildVertexOnEdge(const vertex &This, const int &here, const VertexId this_id)
		{
			// Geometry:
			//   Spawn one new finite vertex at the cut point on a replaced->persisting edge.
			// Code structure:
			//   Takes either free slot from unused or appends new vertex, initializes endpoint and generator
			//   tuple, then validates finite power consistency.
			// edge between This (replaced and finite) and that defined by generators s1,s2,s3 will get a vertex (of newest,s1,s2,s3)
			{
				VertexId persisting_id = This.resolved_endpoint_id(*this, here);
				if (persisting_id == kInvalidId || persisting_id >= vertices.size())
					throw MyException();
				vertex &persisting = vertex_at(persisting_id);
				VertexId builtVertexId = kInvalidId;
				if (_nUnused == 0)
				{
					if (nVertices() == vertices.capacity())
						throw MyException();
					vertices[_nVertices].endPointsAndPositionOverwrite(persisting_id, This.getPowerPointOnLine2(persisting));
					builtVertexId = ++_nVertices - 1;
				}
				else
				{
					const VertexId unused_id = unused[_nUnused - 1];
					if (unused_id == kInvalidId || unused_id >= vertices.size())
						throw MyException();
					vertex &unused_vertex = vertex_at(unused_id);
					unused_vertex.endPointsAndPositionOverwrite(persisting_id, This.getPowerPointOnLine2(persisting));
					builtVertexId = unused_id;
					--_nUnused;
				}
				if (!vertex_at(builtVertexId).Init(This, here, *this, builtVertexId, this_id))
				{
					insertionErrorScale = powerErr;
					return false;
				}
			}

			return true;
		}

		//  void replace_a_vertex(vertex& old_vertex,const Cell& newGenerator);
		//  void checkvertex(vertexIter& myvertex,Cell& newGenerator,std::vector<vertexIter>& replaced,std::vector <vertexIter>&surroundings,const vertexIter former);

		bool FillReplacedPersistingAndInvolved(const CellId this_id, const VertexId start_id)
		{
			// Geometry:
			//   Classify local subgraph into replaced/persisting regions around insertion site.
			// Code structure:
			//   Seeds with start vertex, pushes insertion site into involved refs, then recursively propagates
			//   replacement state using vertex-level checks.
			if (this_id == kInvalidId || this_id >= points.size())
				return false;
			cell &This = cell_at(this_id);
			clear_interna();
			if (start_id == kInvalidId || start_id >= vertices.size())
				return false;
			vertex &start = vertex_at(start_id);
			push_involved(GeneratorRef(GeneratorKind::point, this_id));
			const ReplaceState startState = finiteReplaced(start, this_id);
			if (startState == ReplaceState::ambiguous)
				return false;
			if (startState == ReplaceState::replaced)
			{
				if (start.isCorner())
				{
					return start.cornerToReplacedAndGo(*this, start_id);
				}
				else
					return start.finiteToReplacedAndGo(*this, start_id);
			}
			return true;
		}
		bool CreateFiniteVerticesFromReplaced()
		{
			// Geometry:
			//   Generate the complete set of new finite vertices induced by replaced/persisting interface edges.
			// Code structure:
			//   Counts required edge cuts per replaced vertex, reserves memory if needed, then dispatches to
			//   corner/non-corner build paths.
			// best procedure for new vertices : knowledge : each new (finite) vertex MUST lie on
			// exactly one old EXISTING edge which is NOT disappearing totally
			// all possible edges are the ones coming out our "replaced" vertices
			// so we only try to create if an endPoint of a replaced vertex is not replaced (visitedAs ==-1)
			for (std::size_t replaced_idx = 0; replaced_idx < replaced_size(); ++replaced_idx)
			{
				const VertexId replaced_id = replaced_id_at(replaced_idx);
				if (replaced_id == kInvalidId || replaced_id >= vertices.size())
					continue;
				vertex &replaced_vertex = vertex_at(replaced_id);
				int needed = 0;
				for (int g = dimension; g >= replaced_vertex.isCorner(); g--)
				{
					VertexId endpoint_id = replaced_vertex.resolved_endpoint_id(*this, g);
					if (endpoint_id == kInvalidId || endpoint_id >= vertices.size())
						continue;
					if (vertex_at(endpoint_id).rrv <= 0)
						needed++;
				}
				const int additionalNeeded = needed - static_cast<int>(_nUnused);
				if (additionalNeeded > 0 && nVertices() + additionalNeeded > vertices.capacity())
				{
					this->ReserveNewVertices();
				}
				bool ok = true;
				if (!replaced_vertex.isCorner())
					ok = replaced_vertex.template buildIn<0>(*this, replaced_id);
				else
					ok = replaced_vertex.template buildIn<1>(*this, replaced_id);
				if (!ok)
				{
					return false;
				}
			}
			return true;
		}

		inline void ConnectNewFinitesAmongThemselves3D()
		{
			// Geometry:
			//   Connect newly created vertices into valid polyhedral topology.
			// Code structure:
			//   Uses pair-of-old-generator keys into planes[] as a sparse matching table and stitches both edge
			//   endpoints when second endpoint for a key arrives.
			// we use an InvolvedSize*InvolvedSize-matrix (sparse) and fill in all new edges
			// these are generated by the new Cell and two older Cells
			// we identify vertices to be connected over an new edge by only the two old cells! (new one is everywhere)
			// the earlier dumps itself into the array at spot [a,b] (with a<b) the later simply picks it up!
			if (involved_size() * involved_size() > planes.size())
				planes.resize(involved_size() * involved_size());

			const CellId involved_front_id = involved_id_at(0);
			if (involved_front_id == kInvalidId)
				return;
			cell &involved_front = cell_at(involved_front_id);
			for (unsigned int vertex_idx = 0; vertex_idx < (1u << dimension); ++vertex_idx)
				if (vertices[vertex_idx].rrv > 0)
				{
					vertices[vertex_idx].setPowerData(involved_front);
					set_vertex_generator(vertices[vertex_idx], 0, GeneratorRef(GeneratorKind::point, involved_front_id));
				}

			for (const VertexId vid : involved_front.myVerticesIds)
			// if(!(*it)->isCorner())
			{
				if (vid == kInvalidId || vid >= vertices.size())
					continue;
				vertex_at(vid).registerForConnection3D(*this, vid);
			}
		}

		// const int GoAlongFace(const vertexIter& former, const vertexIter& current,const vertexIter& finish,const CellPtr& Generator1,const Ptr& Generator2,const int function(const int&))const;
		void ReserveNewVertices()
		{
			// Geometry:
			//   Capacity management only; no geometric state change.
			// Code structure:
			//   Reserves larger vertex storage, refreshes index mirrors after reallocation, and restores all
			//   transient containers (replaced/involved/front representatives).
			std::vector<VertexId> _replaced;
			std::vector<VertexId> _currentmyVertices;
			std::vector<VertexId> _first;
			_replaced.reserve(replaced_size());
			_first.reserve(vertices.capacity());

			for (std::size_t replaced_idx = 0; replaced_idx < replaced_size(); ++replaced_idx)
			{
				const VertexId replaced_id = replaced_id_at(replaced_idx);
				if (replaced_id == kInvalidId)
					throw MyException();
				_replaced.push_back(replaced_id);
			}
			const CellId involved_front_id = involved_id_at(0);
			if (involved_front_id == kInvalidId)
				throw MyException();
			cell &involved_front = cell_at(involved_front_id);
			_currentmyVertices.reserve(involved_front.myVerticesIds.size());
			for (const VertexId vid : involved_front.myVerticesIds)
				_currentmyVertices.push_back(vid);
			const std::size_t involved_prefix = involved_front_id;
			for (std::size_t point_idx = 0; point_idx < involved_prefix; ++point_idx)
			{
				const cell &point = points[point_idx];
				_first.push_back(point.myVerticesIds.empty() ? kInvalidId : point.myVerticesIds.front());
			}

			vertices.reserve(2 * vertices.capacity() + 1);

			for (unsigned int i = 0; i < _nVertices; ++i)
				vertices[i].refreshAfterRealloc(*this, i);
			vertices.resize(vertices.capacity());

			for (std::size_t i = 0; i < replaced_size(); ++i)
			{
				set_replaced_id(i, _replaced[i]);
			}
			for (std::size_t i = 0; i < unused.size(); ++i)
				if (!valid_vertex_id(unused[i]))
					throw MyException();
			for (std::size_t i = 0; i < involved_front.myVerticesIds.size(); ++i)
			{
				const VertexId restored_id = _currentmyVertices[i];
				if (restored_id == kInvalidId)
					throw MyException();
				if (restored_id >= vertices.size())
					throw MyException();
				set_cell_my_vertex(involved_front, i, restored_id);
			}
			for (std::size_t i = 0; i < involved_prefix; i++)
			{
				const VertexId restored_id = _first[i];
				if (restored_id == kInvalidId)
					continue;
				if (restored_id >= vertices.size())
					throw MyException();
				set_cell_my_vertex(points[i], 0, restored_id);
			}
		}
		void UpdateUnused()
		{
			// Geometry:
			//   Finalize removed-region cleanup after an insertion.
			// Code structure:
			//   Keeps corners alive, marks finite replaced vertices disconnected, collects reusable slots,
			//   and erases stale myVertices references from non-new cells.
			unused.resize(_nUnused);
			// mark replaced as unused
			for (std::size_t replaced_idx = 0; replaced_idx < replaced_size();)
			{
				const VertexId replaced_id = replaced_id_at(replaced_idx);
				if (replaced_id == kInvalidId || replaced_id >= vertices.size())
				{
					const VertexId last_id = replaced_id_at(replaced_size() - 1);
					set_replaced_id(replaced_idx, last_id);
					pop_replaced();
					continue;
				}
				vertex &replaced_vertex = vertex_at(replaced_id);
				if (replaced_vertex.isCorner())
				{
					replaced_vertex.rrv = 0;
					// Replaced.erase(it);--it;//the corners are always part of diagram
					const VertexId last_id = replaced_id_at(replaced_size() - 1);
					set_replaced_id(replaced_idx, last_id);
					pop_replaced();
				}
				else
				{
					if (replaced_id < nRevertVertices)
						Invalids.push_back(replaced_id);
					replaced_vertex.disconnect(); // vertex has no connection any more.  we delete it later
					++replaced_idx;
				}
			}
			if (nRevertVertices == 0)
				for (std::size_t replaced_idx = 0; replaced_idx < replaced_size(); ++replaced_idx)
				{
					const VertexId replaced_id = replaced_id_at(replaced_idx);
					if (replaced_id == kInvalidId || replaced_id >= vertices.size())
						continue;
					unused.push_back(replaced_id);
				}

			else
				for (std::size_t replaced_idx = 0; replaced_idx < replaced_size(); ++replaced_idx)
				{
					const VertexId replaced_id = replaced_id_at(replaced_idx);
					if (replaced_id == kInvalidId || replaced_id >= vertices.size())
						continue;
					vertex &replaced_vertex = vertex_at(replaced_id);
					if (replaced_id >= nRevertVertices)
						unused.push_back(replaced_id);
					else
						for (const GeneratorRef &ref : replaced_vertex.generatorRefs)
							if (valid_generator_ref(ref))
								erase_cell_my_vertex_by_id(cell_from_ref(ref), replaced_id);
				}
			_nUnused = unused.size();
		}
		void AssignRepresentativeVerticesToCells(const VertexId default_id)
		{
			// Geometry:
			//   Maintain one stable connected representative vertex per involved real cell.
			// Code structure:
			//   Uses new cell representative as fallback and patches cells whose previous representative became invalid.
			const CellId involved_front_id = involved_id_at(0);
			if (involved_front_id == kInvalidId)
				return;
			cell &involved_front = cell_at(involved_front_id);
			// if there are no new vertices, the new cell is covered
			if (involved_front.myVerticesIds.empty())
			{
				if (default_id != kInvalidId && default_id < vertices.size())
					push_cell_my_vertex(involved_front, default_id);
			}
			else // we need one existing vertex close to each cell => we give every cell without representativ a new vertex
			{
				const VertexId new_representative_id = involved_front.myVerticesIds.front();
				if (new_representative_id == kInvalidId || new_representative_id >= vertices.size())
					return;
				for (std::size_t involved_idx = 1; involved_idx < involved_size(); ++involved_idx)
				{
					const CellId involved_id = involved_id_at(involved_idx);
					if (involved_id == kInvalidId)
						continue;
					cell &involved_cell = cell_at(involved_id);
					if (involved_cell.myVerticesIds.empty())
						continue;
					const VertexId representative_id = involved_cell.myVerticesIds.front();
					if (representative_id == kInvalidId || representative_id >= vertices.size())
						continue;
					if (!vertex_at(representative_id).isConnected())				 // if representing vertex has been erased
						set_cell_my_vertex(involved_cell, 0, new_representative_id); // we assign representative of new also to this one
				}
			}
		}

		void SetInvolvedPersistingVisitedToZero()
		{
			// Geometry:
			//   Reset temporary traversal state after insertion/revert local operations.
			// Code structure:
			//   Clears visitedAs on involved refs and resets rrv on involved-front vertices and their incident endpoints.
			for (std::size_t involved_idx = 1; involved_idx < involved_size(); ++involved_idx)
			{
				if (involved_idx >= InvolvedRefs.size())
					continue;
				const GeneratorRef &ref = InvolvedRefs[involved_idx];
				if (!valid_generator_ref(ref))
					continue;
				cell &involved_cell = cell_from_ref(ref);
				involved_cell.visitedAs = 0;
			}
			const CellId involved_front_id = involved_id_at(0);
			if (involved_front_id == kInvalidId)
				return;
			cell &involved_front = cell_at(involved_front_id);
			for (const VertexId vid : involved_front.myVerticesIds)
			{
				if (vid == kInvalidId || vid >= vertices.size())
					continue;
				vertex &vtx = vertex_at(vid);
				vtx.rrv = 0;
				if (!vtx.isCorner())
				{
					VertexId endpoint_id = vtx.resolved_endpoint_id(*this, 0);
					if (endpoint_id != kInvalidId && endpoint_id < vertices.size())
						vertex_at(endpoint_id).rrv = 0;
				}
				else
				{
					for (int g = 1; g <= 3; ++g)
					{
						VertexId endpoint_id = vtx.resolved_endpoint_id(*this, g);
						if (endpoint_id != kInvalidId && endpoint_id < vertices.size())
							vertex_at(endpoint_id).rrv = 0;
					}
				}
			}
		}

	public:
		template <class VectorSubtraction>
		inline /*static*/ PDCoord intersectionOfLineAndPlane3D(const VectorSubtraction &direction, const VectorSubtraction &supportVec, const VectorSubtraction &normal, const PDFloat &planeVal)
		{
			// Compute intersection of a line and power-bisector plane in shifted coordinates.
			const PDFloat tmp = (normal.dot(direction));
			//		const PDFloat sqr=normal.squaredNorm();

			if (tmp > power_err_scaled_epsilon())
				return direction * ((0.5 * (planeVal + normal.squaredNorm()) - normal.dot(supportVec)) / tmp);
			throw MyException();
		}

		struct zeroPoint
		{
			PDFloat pos;
			VertexId fromId;
			int branch;
			std::array<GeneratorRef, dimension> generatorRefs;

			zeroPoint(
				const PDFloat &position,
				const int &way,
				const GeneratorRef &aref = GeneratorRef(),
				const GeneratorRef &bref = GeneratorRef(),
				const GeneratorRef &cref = GeneratorRef(),
				const VertexId origin_id = kInvalidId) : pos(position), fromId(origin_id), branch(way)
			{
				generatorRefs[0] = aref;
				generatorRefs[1] = bref;
				generatorRefs[2] = cref;
			}
		};
		struct vertex
		{
			PDFloat rrv; // relative replace value (power difference)
			bool invalid;
			std::array<GeneratorRef, dimension + 1> generatorRefs;
			PDCoord position;
			PDFloat powerValue;
			std::array<VertexId, dimension + 1> endPointIds;

			friend class PowerDiagram<PDFloat, PDCoord, dimension>;
			inline bool isCorner() const { return endPointIds[0] == kInvalidId; }
			inline bool isOnEdge(const PowerDiagram<PDFloat, PDCoord, dimension> &This)
			{
				GeneratorRef ref = this->resolved_generator_ref(This, dimension - 2);
				return !This.ref_is_real_point(ref);
			}
			inline bool isOnSurface(const PowerDiagram<PDFloat, PDCoord, dimension> &This)
			{
				GeneratorRef ref = this->resolved_generator_ref(This, dimension - 1);
				return !This.ref_is_real_point(ref);
			}
			//	inline int hasVirtualGenerators()const {return (generatorRefs[dimension].is_valid());}
			//	inline int isFinite() const { return generatorRefs[0].is_valid(); }
			inline void disconnect() { invalid = 1; }
			inline int isConnected() const { return !invalid; }

			//  vertex(const vertex& copy);
			inline vertex() : invalid(1)
			{
				for (int g = 0; g <= dimension; ++g)
				{
					endPointIds[g] = kInvalidId;
					generatorRefs[g] = GeneratorRef();
				}
			}

			inline bool Init(
				const vertex &This,
				const int &keep,
				PowerDiagram<PDFloat, PDCoord, dimension> &owner,
				const VertexId self_id,
				const VertexId this_id)
			{
				// Geometry:
				//   Define one new finite vertex from parent replaced vertex by substituting one generator.
				// Code structure:
				//   Copies generator tuple with one dropped slot, assigns insertion generator in slot 0,
				//   rewires reciprocal endpoint in neighbor, and rejects near-zero power unstable cases.
				const CellId involved_front_id = owner.involved_id_at(0);
				if (involved_front_id == kInvalidId)
					return false;
				cell &involved_front = owner.cell_at(involved_front_id);
				this->setPowerData(involved_front);

				for (int g = dimension; g > 0; g--)
					owner.set_vertex_generator(*this, g, This.generatorRefs[g - (g <= keep)]);
				owner.set_vertex_generator(*this, 0, GeneratorRef(GeneratorKind::point, involved_front_id));
				owner.push_cell_my_vertex(involved_front, self_id);
				{
					VertexId endpoint_id = resolved_endpoint_id(owner, 0);
					if (endpoint_id == kInvalidId || endpoint_id >= owner.vertices.size())
						return false;
					vertex &endpoint = owner.vertex_at(endpoint_id);
					const int slot = endpoint.endpoint_slot_to(this_id);
					owner.set_vertex_endpoint_deferred(endpoint, slot, self_id);
				}

				if (owner.within_power_err(powerValue))
				{
					return false;
				}
				return true;
			}

			inline PDCoord getPowerPointOnLine2(const vertex &persist) const
			{
				// Geometry:
				//   Parametric edge interpolation at rrv=0 crossing between replaced and persisting endpoint.
				// Code structure:
				//   Uses already computed endpoint rrv values to avoid recomputing power differences.
				//	const PDCoord PlaneNormal=(b->position-a->position)/*/(a->position-b->position).norm()*/;
				//	const PDFloat PlaneValue=0.5*(PlaneNormal.squaredNorm()+(a->r2-b->r2)/*(a->position-b->position).norm()*/);
				// PlaneNormal and PlaneValue are a factor of (a->position-b->position).norm() too big but they cancel each other out
				//	const PDFloat lowPower=newOne->power(replaced->position);
				//	const PDFloat persistPower=newOne->power(persist->position);
				return ((rrv) / ((rrv) - (persist.rrv))) * (persist.position - position) + position;
			}

			void operator=(const vertex &that)
			{
				generatorRefs = that.generatorRefs;
				position = that.position;
				powerValue = that.powerValue;
				endPointIds = that.endPointIds;
				invalid = that.invalid;
				rrv = that.rrv;
			}

		private:
			inline void setPowerData(const cell &aCell)
			{
				powerValue = (aCell.position - position).squaredNorm() - aCell.r2;
			}
			inline void setTo(const PDCoord pos)
			{
				position = pos;
				for (int g = dimension; g >= 0; g--)
				{
					endPointIds[g] = kInvalidId;
					generatorRefs[g] = GeneratorRef();
				}
			}

			inline PDFloat powerdiff3D(const cell &aCell, const cell &bCell) const
			{
				// Geometry:
				//   Signed power(bCell)-power(aCell) variant evaluated at this vertex.
				// Code structure:
				//   Uses numerically stable algebraic rearrangement to reduce cancellation when cells are close.
				// this has best accuracy when vertex is far away and atoms are close. something similar for a far atom is :
				// -bCall->r2+aCell->r2-(closeCell->position-position).squaredNorm()+2.0*((closeCell->position-position).dot(position-farCell->position)
				return -bCell.r2 + aCell.r2 - (aCell.position - bCell.position).squaredNorm() + 2.0 * ((aCell.position - bCell.position).dot(position - bCell.position));
			}

			template <class PDCalc>
			inline void endPointsAndPositionOverwrite(const VertexId endPointId, const PDCalc &pos)
			{
				// Geometry:
				//   Reinitialize a vertex slot as a provisional cut-vertex with one known endpoint.
				// Code structure:
				//   Writes minimal fields required before full Init(): endpoint 0, position, rrv/invalid reset.
				endPointIds[0] = endPointId;
				rrv = 0;
				invalid = 0;
				position = pos;
			}

			void refreshAfterRealloc(PowerDiagram<PDFloat, PDCoord, dimension> &owner, const VertexId copy_id)
			{
				// Geometry:
				//   No geometric transformation; only preserves incidence consistency after storage move.
				// Code structure:
				//   Revalidates point-generator ownership front entries that track this vertex ID.
				for (int g = dimension; g >= 0; g--)
				{
					GeneratorRef ref = resolved_generator_ref(owner, g);
					if (ref.kind != GeneratorKind::point || ref.index >= owner.points.size())
						continue;
					cell &generator = owner.points[ref.index];
					if ((!generator.myVerticesIds.empty()) && generator.myVerticesIds.front() == copy_id)
						generator.myVerticesIds.front() = copy_id;
				}
			}

			inline void moveAddressNetworkUpdateOnly(
				PowerDiagram<PDFloat, PDCoord, dimension> &owner,
				vertex &whereTo,
				const VertexId this_id,
				const VertexId where_to_id)
			{
				// Geometry:
				//   Topology-preserving move of a vertex record into a recycled slot.
				// Code structure:
				//   Rewrites reciprocal endpoint links in all incident neighbors from this_id to where_to_id,
				//   then copies full payload.
				{
					const VertexId endpoint_id = resolved_endpoint_id(owner, dimension);
					if (endpoint_id != kInvalidId && endpoint_id < owner.vertices.size())
					{
						vertex &endpoint = owner.vertex_at(endpoint_id);
						const int slot = endpoint.endpoint_slot_to(this_id);
						owner.set_vertex_endpoint_deferred(endpoint, slot, where_to_id);
					}
				}
				for (int g = 0; g < dimension; ++g)
				{
					const VertexId endpoint_id = resolved_endpoint_id(owner, g);
					if (endpoint_id == kInvalidId || endpoint_id >= owner.vertices.size())
						continue;
					vertex &endpoint = owner.vertex_at(endpoint_id);
					const int slot = endpoint.endpoint_slot_to(this_id);
					owner.set_vertex_endpoint_deferred(endpoint, slot, where_to_id);
				}

				whereTo = *this;
			}
			inline int endpoint_slot_to(const VertexId comp_id) const
			{
				// Geometry:
				//   Find local edge branch index that points to a given adjacent vertex.
				// Code structure:
				//   Linear scan over endpoint slots; slot 0 is corner sentinel and scanned last by caller contracts.
				for (int g = dimension; g > 0; --g)
					if (endPointIds[g] == comp_id)
						return g;
				return 0;
			}
			inline VertexId resolved_endpoint_id(const PowerDiagram<PDFloat, PDCoord, dimension> &owner, const int g) const
			{
				(void)owner;
				return endPointIds[g];
			}
			inline GeneratorRef resolved_generator_ref(const PowerDiagram<PDFloat, PDCoord, dimension> &owner, const int g) const
			{
				(void)owner;
				return generatorRefs[g];
			}
			inline vertex &resolved_endpoint(PowerDiagram<PDFloat, PDCoord, dimension> &owner, const int g) const
			{
				// Geometry:
				//   Access adjacent vertex at branch g.
				// Code structure:
				//   Validates id range before dereference and throws on broken topology.
				const VertexId id = resolved_endpoint_id(owner, g);
				if (id == kInvalidId || id >= owner.vertices.size())
					throw MyException();
				return owner.vertex_at(id);
			}
			inline cell &resolved_generator(PowerDiagram<PDFloat, PDCoord, dimension> &owner, const int g) const
			{
				// Geometry:
				//   Access generator cell defining this vertex at slot g.
				// Code structure:
				//   Validates GeneratorRef against owner arrays and throws on stale refs.
				GeneratorRef ref = resolved_generator_ref(owner, g);
				if (!owner.valid_generator_ref(ref))
					throw MyException();
				return owner.cell_from_ref(ref);
			}

			bool cornerToReplacedAndGo(PowerDiagram<PDFloat, PDCoord, dimension> &owner, const VertexId self_id)
			{
				// Geometry:
				//   Replacement flood-fill branch for corner vertex.
				// Code structure:
				//   Marks self replaced, accumulates involved generators, assigns corner to insertion front,
				//   and recursively checks incident endpoints.
				owner.push_replaced_id(self_id);
				for (int g = 0; g <= dimension; ++g)
				{
					GeneratorRef ref = resolved_generator_ref(owner, g);
					if (!owner.valid_generator_ref(ref))
						throw MyException();
					cell &generator = owner.cell_from_ref(ref);
					if (generator.visitedAs == 0)
						owner.AddToInvolved(ref);
				}
				const CellId involved_front_id = owner.involved_id_at(0);
				if (involved_front_id != kInvalidId)
				{
					cell &involved_front = owner.cell_at(involved_front_id);
					owner.push_cell_my_vertex(involved_front, self_id); // although replaced it will be part of the new cell!its a corner!
				}

				for (int g = dimension; g > 0; --g)
				{
					vertex &endpoint = resolved_endpoint(owner, g);
					if (endpoint.rrv == 0)
					{
						const VertexId endpoint_id = resolved_endpoint_id(owner, g);
						if (!endpoint.replaceCheck(owner, endpoint_id))
							return false;
					}
				}

				return true;
			}
			bool finiteToReplacedAndGo(PowerDiagram<PDFloat, PDCoord, dimension> &owner, const VertexId self_id)
			{
				// Geometry:
				//   Replacement flood-fill branch for finite (non-corner) vertex.
				// Code structure:
				//   Marks self, accumulates involved generators, and recursively propagates replace checks.
				owner.push_replaced_id(self_id);
				for (int g = 0; g <= dimension; ++g)
				{
					GeneratorRef ref = resolved_generator_ref(owner, g);
					if (!owner.valid_generator_ref(ref))
						throw MyException();
					cell &generator = owner.cell_from_ref(ref);
					if (generator.visitedAs == 0)
						owner.AddToInvolved(ref);
				}

				for (int g = 0; g <= dimension; ++g)
				{
					vertex &endpoint = resolved_endpoint(owner, g);
					if (endpoint.rrv == 0)
					{
						const VertexId endpoint_id = resolved_endpoint_id(owner, g);
						if (!endpoint.replaceCheck(owner, endpoint_id))
							return false;
					}
				}

				return true;
			}
			inline bool replaceCheck(PowerDiagram<PDFloat, PDCoord, dimension> &owner, const VertexId self_id)
			{
				// Geometry:
				//   Decide if this vertex is replaced by insertion site.
				// Code structure:
				//   Dispatches to corner/non-corner logic because endpoint traversal domains differ.
				if (this->isCorner())
					return this->cornerReplaceCheck(owner, self_id);
				return this->finiteReplaceCheck(owner, self_id);
			}

			bool finiteReplaceCheck(PowerDiagram<PDFloat, PDCoord, dimension> &owner, const VertexId self_id)
			{
				// Geometry:
				//   Replacement predicate for finite vertex.
				// Code structure:
				//   Calls finiteReplaced() and executes replacement expansion only on positive classification.
				const CellId involved_front_id = owner.involved_id_at(0);
				if (involved_front_id == kInvalidId)
					return false;
				const typename PowerDiagram<PDFloat, PDCoord, dimension>::ReplaceState state =
					owner.finiteReplaced(*this, involved_front_id);
				if (state == PowerDiagram<PDFloat, PDCoord, dimension>::ReplaceState::ambiguous)
					return false;
				if (state == PowerDiagram<PDFloat, PDCoord, dimension>::ReplaceState::replaced)
					return this->finiteToReplacedAndGo(owner, self_id);
				return true;
			}
			bool cornerReplaceCheck(PowerDiagram<PDFloat, PDCoord, dimension> &owner, const VertexId self_id)
			{
				// Geometry:
				//   Replacement predicate for corner vertex.
				// Code structure:
				//   Same classifier as finite case, but expansion path preserves corner-specific endpoint handling.
				const CellId involved_front_id = owner.involved_id_at(0);
				if (involved_front_id == kInvalidId)
					return false;
				const typename PowerDiagram<PDFloat, PDCoord, dimension>::ReplaceState state =
					owner.finiteReplaced(*this, involved_front_id);
				if (state == PowerDiagram<PDFloat, PDCoord, dimension>::ReplaceState::ambiguous)
					return false;
				if (state == PowerDiagram<PDFloat, PDCoord, dimension>::ReplaceState::replaced)
					return this->cornerToReplacedAndGo(owner, self_id);
				return true;
			}
			template <const int cornerInfo>
			bool buildIn(PowerDiagram<PDFloat, PDCoord, dimension> &pd, const VertexId self_id) const
			{
				// Geometry:
				//   For each surviving outgoing edge, create one intersection vertex with insertion boundary.
				// Code structure:
				//   Scans valid endpoint branches and delegates actual creation to tryToBuildVertexOnEdge().
				for (int g = dimension; g >= cornerInfo; g--)
				{
					VertexId endpoint_id = this->resolved_endpoint_id(pd, g);
					if (endpoint_id == kInvalidId || endpoint_id >= pd.vertices.size())
						continue;
					if (pd.vertex_at(endpoint_id).rrv <= 0)
					{
						if (!pd.tryToBuildVertexOnEdge(*this, g, self_id))
							return false;
					}
				}
				return true;
			}

			inline void registerForConnection3D(PowerDiagram<PDFloat, PDCoord, dimension> &owner, const VertexId self_id)
			{
				// Geometry:
				//   Register a newly built vertex into the edge-matching map for 3D cell-facet stitching.
				// Code structure:
				//   Uses visitedAs indices of old generators to compute deterministic key pairs in planes[].
				const int g1 = resolved_generator(owner, 1).visitedAs;
				const int g2 = resolved_generator(owner, 2).visitedAs;
				const int g3 = resolved_generator(owner, 3).visitedAs;
				owner.planes[g2 * owner.involved_size() + g1].storeOrConnect(owner, self_id, 3);
				owner.planes[g3 * owner.involved_size() + g1].storeOrConnect(owner, self_id, 2);
				owner.planes[g3 * owner.involved_size() + g2].storeOrConnect(owner, self_id, 1);
			}
		};
		struct EdgeEnds
		{
			VertexId aId;
			int aSlot;
			inline EdgeEnds() : aId(kInvalidId), aSlot(-1) {}
			inline void storeOrConnect(PowerDiagram<PDFloat, PDCoord, dimension> &owner, const VertexId pvertex_id, const int slot)
			{
				// First endpoint stores itself; second endpoint closes the pair and connects both vertices.
				if (this->aId == kInvalidId)
				{
					this->aId = pvertex_id;
					this->aSlot = slot;
				}
				else
				{
					vertex &pvertex = owner.vertex_at(pvertex_id);
					vertex &other = owner.vertex_at(this->aId);
					owner.set_vertex_endpoint_deferred(pvertex, slot, this->aId);
					owner.set_vertex_endpoint_deferred(other, this->aSlot, pvertex_id);
					this->aId = kInvalidId;
					this->aSlot = -1;
				}
			}
			inline void connect(PowerDiagram<PDFloat, PDCoord, dimension> &owner, const VertexId pvertex_id, const int slot)
			{
				// Connect against a previously stored endpoint if one exists for this key.
				if (this->aId != kInvalidId)
				{
					vertex &pvertex = owner.vertex_at(pvertex_id);
					vertex &other = owner.vertex_at(this->aId);
					owner.set_vertex_endpoint_deferred(pvertex, slot, this->aId);
					owner.set_vertex_endpoint_deferred(other, this->aSlot, pvertex_id);
					this->aId = kInvalidId;
					this->aSlot = -1;
				}
			}
		};
	};
	// static initializations
	// template <class PDFloat, class PDCoord>
	// std::vector<cell<PDFloat,PDCoord> > PowerDiagram<PDFloat,PDCoord>::sideGenerators; //outside

} // namespace POWER_DIAGRAM
#endif /* POWER_DIAGRAM_H_ */
