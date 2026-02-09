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
#include <array>
#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <algorithm>
#include <cstddef>
#include <limits>
#include <ctime>
//#include "basic_vector_calc.h"


namespace POWER_DIAGRAM
{

class MyException : public std::exception
{
public:
	MyException() {	}
};
class IdenticalPointException : public std::exception {};
class VerticesFullException : public std::exception {};




// Return the nth natural number of a sequence where one index ("without") is skipped.
// Used to compare generator triples while omitting one generator position.
inline int nth(const int n,const int without)
{
	return n+(without<=n);
}
template <class PDCoord, class PDFloat, class Pos_iterator, class Strength_iterator, const int dimension>
void getBoundingBox(PDCoord& lowestCorner,PDCoord& highestCorner,const unsigned int& size,const Pos_iterator pos_begin,const Strength_iterator  strength_begin,const PDFloat additionalCubeSize=pow(2.0,1.0/dimension)-1)
{
	// Build a radius-aware axis-aligned box covering all weighted points.
	// The optional inflation factor keeps the initial clipping cube safely outside the data.
	if(size>0)
	{
		lowestCorner=*pos_begin;
		highestCorner=*pos_begin;

		Pos_iterator pos_end=pos_begin;
		Strength_iterator strength_end=strength_begin;

		for(unsigned int i=0;i<size;i++)
		{
			for(int g=dimension-1;g>=0;g--)
			{
				if((*pos_end)[g]-(*strength_end)<lowestCorner[g])lowestCorner[g]=(*pos_end)[g]-(*strength_end);
				if((*pos_end)[g]+(*strength_end)>highestCorner[g])highestCorner[g]=(*pos_end)[g]+(*strength_end);
			}
			++pos_end;
			++strength_end;
		}


		const PDCoord center=0.5*(lowestCorner+highestCorner);
		lowestCorner+=/*PowerDiagram<PDFloat,PDCoord,dimension>::error((highestCorner+lowestCorner).norm())**/(lowestCorner-center)*additionalCubeSize;
		highestCorner+=/*PowerDiagram<PDFloat,PDCoord,dimension>::error((highestCorner+lowestCorner).norm())**/(highestCorner-center)*additionalCubeSize;
	}
}




struct PowerDiagramRuntimeParams
{
	bool radiiGiven;//otherwise the second input array is interpreted as power
	bool fill_myVertices;//otherwise no cell information is constructed, but only a diagram of vertices
	bool fill_neighbours;//otherwise no neighbourhood-information is extracted from myVertices
	bool fill_zeroPoints;//otherwise no information where power diagram is zero is created
	bool with_warnings;//otherwise powerdiagram is in silent mode and will not tell when reducing powers
	bool without_check;
	PowerDiagramRuntimeParams(const bool _radiiGiven,const bool _fill_myVertices,const bool _fill_neighbours,const bool _fill_zeroPoints,const bool _with_warnings,bool _without_check)
	{
		radiiGiven=_radiiGiven;
		fill_myVertices=_fill_myVertices;
		fill_neighbours=_fill_neighbours;
		fill_zeroPoints=_fill_zeroPoints;
		with_warnings=_with_warnings;
		without_check=_without_check;
	}
};
template <class PDFloat, class PDCoord,class Pos_iterator, class Strength_iterator, class BondTo_iterator>
struct PowerDiagramParams
{
	bool create_vertices;//otherwise the object is only constructed, but no calculation is done
	PowerDiagramRuntimeParams runpar;

	const PDCoord lowestCorner;
	const PDCoord highestCorner;

	const unsigned int size;
	Pos_iterator pos_begin;
	Strength_iterator strength_begin;
	BondTo_iterator bondTo_begin;
	
	PowerDiagramParams with_radiiGiven(const bool yon){runpar.radiiGiven=yon; return *this;};
	PowerDiagramParams with_calculate(const bool yon){create_vertices=yon; return *this;};
	PowerDiagramParams with_myVertices(const bool yon){runpar.fill_myVertices=yon; return *this;};
	PowerDiagramParams with_cells(const bool yon){runpar.fill_neighbours=yon;runpar.fill_myVertices=yon; return *this;};
	PowerDiagramParams with_zeroPoints(const bool yon){runpar.fill_zeroPoints=yon;runpar.fill_myVertices=yon; return *this;};
	PowerDiagramParams with_Warnings(const bool yon){runpar.with_warnings=yon;return *this;};
	PowerDiagramParams without_Check(const bool yon){runpar.without_check=yon;return *this;};
	PowerDiagramParams(const unsigned int size_, Pos_iterator& pos_begin_, Strength_iterator& strength_begin_, BondTo_iterator bondTo_begin_, PDCoord& lc, PDCoord& hc,
				 bool radiiGiven_=1, bool create_vertices_=1, bool fill_cellVertices=1,bool fill_cellNeighbours=1,bool fill_zeroPoints=1,bool withWarnings=1,bool withoutChecks=1):
		create_vertices(create_vertices_),runpar(PowerDiagramRuntimeParams(radiiGiven_,fill_cellVertices,fill_cellNeighbours,fill_zeroPoints,withWarnings,withoutChecks)),lowestCorner(lc),highestCorner(hc),
		size(size_),pos_begin(pos_begin_),strength_begin(strength_begin_),bondTo_begin(bondTo_begin_)
	{
	}
};



template <class PDFloat, class PDCoord,const int dimension>
class PowerDiagram
{
public:
	//Forward declarations for cell and vertex
	struct zeroPoint;
	struct vertex;
	struct cell;
	struct EdgeEnds;
	typedef cell* cellPtr;
	typedef cell const* const_cellPtr;
	typedef vertex* vertexPtr;
	typedef vertex const* const_vertexPtr;
	using CellId = std::size_t;
	using VertexId = std::size_t;
	using ZeroId = std::size_t;
	inline static constexpr std::size_t kInvalidId = std::numeric_limits<std::size_t>::max();
	enum class GeneratorKind { point, side };
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

	std::vector<vertexPtr> unused;
	std::vector< cell > points;
	std::vector< vertex > vertices;//we use the plane space of a vector because of speed. never push_back or it will cause a realloc !!!
	std::vector< zeroPoint> zeros;//this is a vector of all zero points, but maybe some are inactive (inactives do not appear in myZeros!)
	std::vector< cell > sideGenerators;


	unsigned int nRevertVertices;
	unsigned int nRevertZeros;
	unsigned int nRevertPoints;
	int cornerOwners[1<<dimension];

	PDFloat powerErr;
	PDFloat insertionErrorScale;
	std::vector<vertexPtr> Replaced;//old vertices that are removed
	std::vector<vertexPtr> Invalids;//old vertices that are removed reloaded
	std::vector<cellPtr> Involved; // all cells which are involved
	std::vector< EdgeEnds> planes;//used for connecting new vertices, stores "open ends"
	enum class ReplaceState
	{
		persisting,
		replaced,
		ambiguous
	};
public :
		struct cell
		{
			int visitedAs;
			PDCoord position;
			PDFloat r;
			PDFloat r2;
			typedef zeroPoint* zeroPointPtr;
			cellPtr bondTo;
			std::vector<cellPtr> neighbours;
			std::vector<vertexPtr> myVertices;
			CellId bondToId;
			std::vector<CellId> neighboursIds;
			std::vector<VertexId> myVerticesIds;
			std::vector<int> myZeroPoints;

			inline cell(PDCoord const& pos, PDFloat const& root,PDFloat const& power,const cellPtr bond):visitedAs(0),position(pos),r(root),r2(power),bondTo(bond),bondToId(kInvalidId)
			{
				myVertices.reserve(12);
				myVerticesIds.reserve(12);
			}
			inline cell(PDCoord const& pos,const PDFloat& str,const cellPtr bond):visitedAs(0),position(pos),r(str),r2(str*str),bondTo(bond),bondToId(kInvalidId)
			{
				myVertices.reserve(12);
				myVerticesIds.reserve(12);
			}
		inline PDFloat power(const PDCoord& coord) const
		{
			// Power distance to this weighted site.
			return (position-coord).squaredNorm()-r2;
		}
		inline bool isReal(const PowerDiagram<PDFloat,PDCoord,dimension>& This)
		{
			// Distinguish true input sites from temporary cube side generators.
			return (this>=&This.getPoints()[0]&&this<&This.getPoints()[This.getPoints().size()]);
		}
	private:

	};
	inline void AddToInvolved(cell& thit)
	{
		// Mark a cell as part of the local insertion neighborhood.
		thit.visitedAs=Involved.size();
		Involved.push_back(&thit);
	}
	clock_t t1,t2,t3,t4,t5,t6;
	cell const & get_point(const int n) const { return points[n]; }
	std::vector< cell > const & get_points() const { return points; }
	inline CellId get_cell_id(const cell& my_cell) const
	{
		return static_cast<CellId>(&my_cell - points.data());
	}
	inline VertexId get_vertex_id(const vertex& my_vertex) const
	{
		return static_cast<VertexId>(&my_vertex - vertices.data());
	}
	inline CellId cell_id_or_invalid(const cellPtr ptr) const
	{
		if (ptr == NULL || points.empty()) return kInvalidId;
		const cell* const first = points.data();
		const cell* const last = first + points.size();
		if (ptr < first || ptr >= last) return kInvalidId;
		return static_cast<CellId>(ptr - first);
	}
	inline CellId cell_id_or_invalid(const const_cellPtr ptr) const
	{
		if (ptr == NULL || points.empty()) return kInvalidId;
		const cell* const first = points.data();
		const cell* const last = first + points.size();
		if (ptr < first || ptr >= last) return kInvalidId;
		return static_cast<CellId>(ptr - first);
	}
	inline VertexId vertex_id_or_invalid(const vertexPtr ptr) const
	{
		if (ptr == NULL || vertices.empty()) return kInvalidId;
		const vertex* const first = vertices.data();
		const vertex* const last = first + vertices.size();
		if (ptr < first || ptr >= last) return kInvalidId;
		return static_cast<VertexId>(ptr - first);
	}
	inline VertexId vertex_id_or_invalid(const const_vertexPtr ptr) const
	{
		if (ptr == NULL || vertices.empty()) return kInvalidId;
		const vertex* const first = vertices.data();
		const vertex* const last = first + vertices.size();
		if (ptr < first || ptr >= last) return kInvalidId;
		return static_cast<VertexId>(ptr - first);
	}
	inline GeneratorRef generator_ref_or_invalid(const cellPtr ptr) const
	{
		if (ptr == NULL) return GeneratorRef();
		if (!points.empty())
		{
			const cell* const pfirst = points.data();
			const cell* const plast = pfirst + points.size();
			if (ptr >= pfirst && ptr < plast)
			{
				return GeneratorRef(GeneratorKind::point, static_cast<std::size_t>(ptr - pfirst));
			}
		}
		if (!sideGenerators.empty())
		{
			const cell* const sfirst = sideGenerators.data();
			const cell* const slast = sfirst + sideGenerators.size();
			if (ptr >= sfirst && ptr < slast)
			{
				return GeneratorRef(GeneratorKind::side, static_cast<std::size_t>(ptr - sfirst));
			}
		}
		return GeneratorRef();
	}
	inline void sync_cell_link_mirrors(cell& a_cell) const
	{
		a_cell.bondToId = cell_id_or_invalid(a_cell.bondTo);
		a_cell.neighboursIds.resize(a_cell.neighbours.size(), kInvalidId);
		for (std::size_t i = 0; i < a_cell.neighbours.size(); ++i)
		{
			a_cell.neighboursIds[i] = cell_id_or_invalid(a_cell.neighbours[i]);
		}
		a_cell.myVerticesIds.resize(a_cell.myVertices.size(), kInvalidId);
		for (std::size_t i = 0; i < a_cell.myVertices.size(); ++i)
		{
			a_cell.myVerticesIds[i] = vertex_id_or_invalid(a_cell.myVertices[i]);
		}
	}
	inline void sync_vertex_link_mirrors(vertex& a_vertex) const
	{
		for (int g = 0; g <= dimension; ++g)
		{
			a_vertex.generatorRefs[g] = generator_ref_or_invalid(a_vertex.generators[g]);
			a_vertex.endPointIds[g] = vertex_id_or_invalid(a_vertex.endPoints[g]);
		}
	}
	inline void sync_zero_link_mirrors(zeroPoint& a_zero) const
	{
		a_zero.fromId = vertex_id_or_invalid(a_zero.from);
		for (int g = 0; g < dimension; ++g)
		{
			a_zero.generatorRefs[g] = generator_ref_or_invalid(a_zero.generators[g]);
		}
	}
	inline void sync_all_link_mirrors()
	{
		for (cell& a_cell : points) sync_cell_link_mirrors(a_cell);
		for (unsigned int vi = 0; vi < _nVertices; ++vi) sync_vertex_link_mirrors(vertices[vi]);
		for (zeroPoint& a_zero : zeros) sync_zero_link_mirrors(a_zero);
	}
	inline cellPtr cell_ptr_from_id(const CellId id)
	{
		return (id == kInvalidId || id >= points.size()) ? NULL : &points[id];
	}
	inline const_cellPtr cell_ptr_from_id_const(const CellId id) const
	{
		return (id == kInvalidId || id >= points.size()) ? NULL : &points[id];
	}
	inline const_cellPtr cell_ptr_from_ref_const(const GeneratorRef& ref) const
	{
		if (!ref.is_valid()) return NULL;
		if (ref.kind == GeneratorKind::point && ref.index < points.size()) return &points[ref.index];
		if (ref.kind == GeneratorKind::side && ref.index < sideGenerators.size()) return &sideGenerators[ref.index];
		return NULL;
	}
	inline vertexPtr vertex_ptr_from_id(const VertexId id)
	{
		return (id == kInvalidId || id >= vertices.size()) ? NULL : &vertices[id];
	}
	inline void set_bond_to(cell& a_cell, const cellPtr ptr)
	{
		a_cell.bondTo = ptr;
		a_cell.bondToId = cell_id_or_invalid(ptr);
	}
	inline void set_bond_to_id(cell& a_cell, const CellId id)
	{
		a_cell.bondToId = id;
		a_cell.bondTo = cell_ptr_from_id(id);
	}
	inline void clear_cell_my_vertices(cell& a_cell)
	{
		a_cell.myVertices.clear();
		a_cell.myVerticesIds.clear();
	}
	inline void push_cell_my_vertex(cell& a_cell, const vertexPtr ptr)
	{
		a_cell.myVertices.push_back(ptr);
		a_cell.myVerticesIds.push_back(vertex_id_or_invalid(ptr));
	}
	inline void pop_cell_my_vertex(cell& a_cell)
	{
		a_cell.myVertices.pop_back();
		a_cell.myVerticesIds.pop_back();
	}
	inline void set_cell_my_vertex(cell& a_cell, const std::size_t index, const vertexPtr ptr)
	{
		a_cell.myVertices[index] = ptr;
		a_cell.myVerticesIds[index] = vertex_id_or_invalid(ptr);
	}
	inline void erase_cell_my_vertex(cell& a_cell, const std::size_t index)
	{
		a_cell.myVertices.erase(a_cell.myVertices.begin() + index);
		a_cell.myVerticesIds.erase(a_cell.myVerticesIds.begin() + index);
	}
	inline void clear_cell_neighbours(cell& a_cell)
	{
		a_cell.neighbours.clear();
		a_cell.neighboursIds.clear();
	}
	inline void push_cell_neighbour(cell& a_cell, const cellPtr ptr)
	{
		a_cell.neighbours.push_back(ptr);
		a_cell.neighboursIds.push_back(cell_id_or_invalid(ptr));
	}
	inline void erase_cell_neighbour(cell& a_cell, const std::size_t index)
	{
		a_cell.neighbours.erase(a_cell.neighbours.begin() + index);
		a_cell.neighboursIds.erase(a_cell.neighboursIds.begin() + index);
	}
	inline zeroPoint const& get_zero(const ZeroId n) const { return zeros[n]; }
	inline zeroPoint& get_zero(const ZeroId n) { return zeros[n]; }
	inline vertex const& get_vertex(const VertexId n) const { return vertices[n]; }
	inline vertex& get_vertex(const VertexId n) { return vertices[n]; }
	inline cell const& get_cell(const CellId n) const { return points[n]; }
	inline cell& get_cell(const CellId n) { return points[n]; }
	inline cell const& get_generator(const GeneratorRef& ref) const
	{
		return (ref.kind == GeneratorKind::point) ? points[ref.index] : sideGenerators[ref.index];
	}
	inline cell& get_generator(const GeneratorRef& ref)
	{
		return (ref.kind == GeneratorKind::point) ? points[ref.index] : sideGenerators[ref.index];
	}
	inline GeneratorRef get_point_ref(const CellId index) const { return GeneratorRef(GeneratorKind::point, index); }
	inline GeneratorRef get_side_ref(const std::size_t index) const { return GeneratorRef(GeneratorKind::side, index); }
		        unsigned int get_point_num ( cell const & my_cell )
		        {
		            return static_cast<unsigned int>(get_cell_id(my_cell));
		        }
	std::vector< vertex> const & get_vertices() const { return vertices; }

	std::vector< zeroPoint> const & get_zeroPoints() const { return zeros; }
	inline const_vertexPtr zeroPointFrom(const zeroPoint& zp) const
	{
		if (zp.fromId != kInvalidId && zp.fromId < _nVertices) return &vertices[zp.fromId];
		return zp.from;
	}
	inline const_cellPtr zeroPointGenerator(const zeroPoint& zp, const int g) const
	{
		if (g < 0 || g >= dimension) return NULL;
		const GeneratorRef& ref = zp.generatorRefs[g];
		if (ref.is_valid())
		{
			if (ref.kind == GeneratorKind::point && ref.index < points.size()) return &points[ref.index];
			if (ref.kind == GeneratorKind::side && ref.index < sideGenerators.size()) return &sideGenerators[ref.index];
		}
		return zp.generators[g];
	}
	inline bool zeroPointValid(const zeroPoint& zp) const
	{
		const const_vertexPtr from_ptr = zeroPointFrom(zp);
		if (from_ptr == NULL) return false;
		const vertexPtr to_ptr = from_ptr->endPoints[zp.branch];
		if (to_ptr == NULL) return false;
		return ((!from_ptr->invalid) && (!to_ptr->invalid));
	}
	inline PDCoord zeroPointPos(const zeroPoint& zp) const
	{
		const const_vertexPtr from_ptr = zeroPointFrom(zp);
		if (from_ptr == NULL) throw MyException();
		const vertexPtr to_ptr = from_ptr->endPoints[zp.branch];
		if (to_ptr == NULL) throw MyException();
		return (to_ptr->position) * zp.pos - from_ptr->position * (zp.pos - static_cast<PDFloat>(1.0));
	}


	inline static PDFloat error(const PDFloat &f) 
	{
		// Scale machine epsilon to the magnitude of f with a non-zero floor.
		if(f>std::numeric_limits<PDFloat>::min()/std::numeric_limits<PDFloat>::epsilon())
			return f*(std::numeric_limits<PDFloat>::epsilon());
		else if(f<-std::numeric_limits<PDFloat>::min()/std::numeric_limits<PDFloat>::epsilon())
			return -f*(std::numeric_limits<PDFloat>::epsilon());
		else return std::numeric_limits<PDFloat>::min()/std::numeric_limits<PDFloat>::epsilon();
	}
	inline bool above_power_err(const PDFloat value) const { return value > powerErr; }
	inline bool below_neg_power_err(const PDFloat value) const { return value < -powerErr; }
	inline bool within_power_err(const PDFloat value) const { return std::abs(value) <= powerErr; }
	inline bool below_power_err(const PDFloat value) const { return value < powerErr; }
	inline PDFloat power_err_scaled_epsilon() const { return powerErr * std::numeric_limits<PDFloat>::epsilon(); }
template <class Pos_iterator, class Strength_iterator, class BondTo_iterator>
static PowerDiagramParams<PDFloat,PDCoord,Pos_iterator,Strength_iterator,BondTo_iterator> create(unsigned int size, Pos_iterator pos_begin, Strength_iterator strength_begin, BondTo_iterator bondTo_begin)
{
	// Factory helper: derive bounding cube and return a fully parameterized construction object.
	PDCoord highestCorner;
	PDCoord lowestCorner;
	{
		if(size>=1)
			getBoundingBox<PDCoord,PDFloat,Pos_iterator,Strength_iterator,3>(lowestCorner,highestCorner,size,pos_begin,strength_begin);
		else{std::cout<<"create empty PD(not implemented, yet)"<<std::endl;throw MyException();}
		return PowerDiagramParams<PDFloat,PDCoord,Pos_iterator,Strength_iterator,BondTo_iterator>	(size,pos_begin,strength_begin,bondTo_begin,lowestCorner,highestCorner);
	}
}
template <typename Pos_iterator, typename Strength_iterator, typename BondTo_iterator>
		PowerDiagram(PowerDiagramParams<PDFloat,PDCoord,Pos_iterator,Strength_iterator,BondTo_iterator> _params):center(0.5*(_params.highestCorner+_params.lowestCorner)),params(_params.runpar)
		{
			// Build initial clipping cube, insert all points, then derive optional adjacency/zero-point caches.
			_nUnused=0;
			insertionErrorScale=0;
			nRevertPoints=0;
		nRevertZeros=0;
		nRevertVertices=(1<<dimension);
		Replaced.reserve(64);
		planes.resize(64*64);
		points.reserve(_params.size);
		vertices.reserve(_params.size*32+(1<<dimension));
		vertices.resize(vertices.capacity());

		PDCoord lowest=_params.lowestCorner;
		PDCoord highest=_params.highestCorner;
		buildCube(PDCoord(lowest-center),PDCoord(highest-center));
		//build and connect corners
		
		Pos_iterator pos_it=_params.pos_begin;
		Strength_iterator strength_it=_params.strength_begin;
		BondTo_iterator bondTo_it=_params.bondTo_begin;

			if(_params.size>0)
			{
				if(params.radiiGiven)
				{
					points.push_back(cell(*pos_it-center,*strength_it,NULL));
					set_bond_to(points.back(), NULL);
					for(unsigned int i=1;i<_params.size;i++)
					{
						points.push_back(cell(*(++pos_it)-center,*(++strength_it),&points[*(++bondTo_it)]));
						set_bond_to(points.back(), &points[*(bondTo_it)]);
					}
				}
				else
				{
					points.push_back(cell(*pos_it-center,sqrt(*strength_it),*strength_it,NULL));
					set_bond_to(points.back(), NULL);
					for(unsigned int i=1;i<_params.size;i++)
					{
						points.push_back(cell(*(++pos_it)-center,sqrt(*(++strength_it)),*(strength_it),&points[*(++bondTo_it)]));
						set_bond_to(points.back(), &points[*(bondTo_it)]);
					}
				}
			}

		if(_params.create_vertices)
			buildVertices(points.size());


		if(params.fill_myVertices)
			FillAllMyVertices();
		if(params.fill_neighbours)
			FillAllNeighbours();
			if(params.fill_zeroPoints)
				FillAllZeroPoints();
			sync_all_link_mirrors();


		}
	const std::vector<cell >& getPoints()const {return points;}

		void revert()
		{
			//if addmore was used this function can revert to the diagram without the added atoms
			// Roll back topology and cached adjacency/zero-point state to the snapshot taken before addMore().
				for(unsigned int vi=nRevertVertices;vi<_nVertices;++vi)
					for(int gi=1;gi<=dimension;++gi)
						if(vertices[vi].generators[gi]->isReal(*this))
							pop_cell_my_vertex(*vertices[vi].generators[gi]);
		_nVertices=nRevertVertices;
		Involved.clear();
		if(points.size()>nRevertPoints)
			points.erase(points.begin()+nRevertPoints,points.end());
			for(int c=0;c<(1<<dimension);c++)
			{
				const CellId owner_id = (cornerOwners[c] < 0) ? kInvalidId : static_cast<CellId>(cornerOwners[c]);
				cellPtr owner_ptr = cell_ptr_from_id(owner_id);
				if(owner_ptr == NULL) owner_ptr = points.data();
				vertices[c].generators[0]=owner_ptr;
				vertices[c].powerValue=vertices[c].generators[0]->power(vertices[c].position);
			}
			for(vertexPtr it : Invalids)
			{
				it->invalid=0;
				it->rrv=0;
				for(int endpoint_idx=it->isCorner();endpoint_idx<=dimension;++endpoint_idx)
				{
					vertexPtr endpoint = it->endPoints[endpoint_idx];
					for(int g1=it->isCorner();g1<=dimension;g1++)
					for(int g2=endpoint->isCorner();g2<=dimension;g2++)
					{
						if(it->generators[nth(0,g1)]==endpoint->generators[nth(0,g2)]&&it->generators[nth(1,g1)]==endpoint->generators[nth(1,g2)]&&it->generators[nth(2,g1)]==endpoint->generators[nth(2,g2)])
						{
	//						(*it)->endPoints[g1]=*it2;		//is already set 
							endpoint->endPoints[g2]=&(*it);
						}
					}
				}
					for(int g=0;g<=dimension;++g)
						if(it->generators[g]->isReal(*this))
						{
							push_cell_my_vertex(*it->generators[g], &(*it));
							if(it->generators[g]->visitedAs==0)
							{
							it->generators[g]->visitedAs=-1;
							Involved.push_back(it->generators[g]);
						}
					}
			}
		Invalids.clear();
		if(params.fill_neighbours)
			FillAllNeighboursOfInvolved();

			if(params.fill_zeroPoints)
			{
				for(const cellPtr involved_cell : Involved)
					while(!involved_cell->myZeroPoints.empty()&&involved_cell->myZeroPoints.back()>static_cast<int>(nRevertZeros))
						involved_cell->myZeroPoints.pop_back();
				zeros.erase(zeros.begin()+nRevertZeros,zeros.end());
			}
		nRevertPoints=0;
		nRevertVertices=0;
			for(int c=0;c<(1<<dimension);c++)
				cornerOwners[c]=0;
			sync_all_link_mirrors();
		}
	template <class Pos_iterator, class Strength_iterator>
	inline void recalculate(const Pos_iterator pos_it,const Strength_iterator strength_it,const unsigned int size)
	//does standard deletion and calculation of Vertices, neighbour information, ...
	{
		// Rebuild the full diagram from scratch for updated coordinates/radii.
		clearAllmyVertices();
		clear_interna();
			if(size>points.size())
			{
				std::vector<CellId> old_bond_ids(points.size(), kInvalidId);
				for(std::size_t i=1;i<points.size();++i)
				{
					old_bond_ids[i] = (points[i].bondToId != kInvalidId)
						? points[i].bondToId
						: cell_id_or_invalid(points[i].bondTo);
				}
				const cellPtr old = points.empty() ? NULL : points.data();
				points.reserve(size);
					if(old != NULL && old != points.data())
						for(std::size_t i=1;i<points.size();++i)
						{
							if(old_bond_ids[i] != kInvalidId) set_bond_to_id(points[i], old_bond_ids[i]);
						}
				}
		zeros.clear();
		PDCoord lowest,highest;
		getBoundingBox<PDCoord,PDFloat,PDCoord const*,PDFloat const*,dimension>(lowest,highest,size,&(*pos_it),&(*strength_it));
		buildCube(lowest-center,highest-center);

		if(params.radiiGiven)
			{
				for(unsigned int i=0;i<points.size();i++)
				{
					points[i].position=(*(pos_it+i))-center;
					points[i].r=(*(strength_it+i));
					points[i].r2=points[i].r*points[i].r;
					points[i].visitedAs=0;
				}

						for(unsigned int i=points.size();i<size;i++)
						{
							points.push_back(cell((*(pos_it+i))-center,(*(strength_it+i)),&points.back()));
							set_bond_to(points.back(), points.size() > 1 ? &points[points.size()-2] : NULL);
						}

			}
			else
			{
				for(unsigned int i=0;i<points.size();i++)
				{
					points[i].position=(*(pos_it+i))-center;
					points[i].r2=(*(strength_it+i));
					points[i].r=sqrt(points[i].r2);
					points[i].visitedAs=0;
				}
						for(unsigned int i=points.size();i<size;i++)
						{
							points.push_back(cell((*(pos_it+i))-center,sqrt(*(strength_it+i)),*(strength_it+i),&points.back()));
							set_bond_to(points.back(), points.size() > 1 ? &points[points.size()-2] : NULL);
						}
				}

		buildVertices(size);
		if(params.fill_myVertices||params.fill_neighbours)
			FillAllMyVertices();
		if(params.fill_neighbours)
			FillAllNeighbours();
			if(params.fill_zeroPoints)
				FillAllZeroPoints();
			sync_all_link_mirrors();
		}


	template <class Pos_iterator, class Strength_iterator>
	inline void addMore(const Pos_iterator pos_it,const Strength_iterator strength_it,const int _newSize)
	//does standard deletion and calculation of Vertices, neighbour information, ...
	{
		// Incrementally insert additional points while preserving a revert snapshot.
		const unsigned int gap=_newSize<points.size()?points.size()-_newSize:1;
		const unsigned int newSize=_newSize<points.size()?points.size()+1:_newSize;
			nRevertVertices=_nVertices;
			nRevertZeros=zeros.size();
			nRevertPoints=points.size();
			std::vector<CellId> old_bond_ids(nRevertPoints, kInvalidId);
			std::vector<std::vector<CellId> > old_neighbour_ids(nRevertPoints);
			std::vector<std::array<GeneratorRef,dimension+1> > old_vertex_generator_refs(_nVertices);
			for(std::size_t i=0;i<nRevertPoints;++i)
			{
				old_bond_ids[i] = (points[i].bondToId != kInvalidId) ? points[i].bondToId : cell_id_or_invalid(points[i].bondTo);
				old_neighbour_ids[i].resize(points[i].neighbours.size(), kInvalidId);
				for(std::size_t j=0;j<points[i].neighbours.size();++j)
				{
					old_neighbour_ids[i][j] = (!points[i].neighboursIds.empty() && j < points[i].neighboursIds.size())
						? points[i].neighboursIds[j]
						: cell_id_or_invalid(points[i].neighbours[j]);
				}
			}
			for(unsigned int vi=0;vi<_nVertices;++vi)
				for(int g=0;g<=dimension;++g)
					old_vertex_generator_refs[vi][g] = vertices[vi].generatorRefs[g].is_valid()
						? vertices[vi].generatorRefs[g]
						: generator_ref_or_invalid(vertices[vi].generators[g]);
			for(int c=0;c<(1<<dimension);c++)
			{
				const CellId owner_id = cell_id_or_invalid(vertices[c].generators[0]);
				cornerOwners[c]=(owner_id == kInvalidId) ? 0 : static_cast<int>(owner_id);
			}
			const const_cellPtr oldstorage=points.data();
		points.reserve(newSize);
		if(params.radiiGiven)
			{
					for(unsigned int i=0;i<newSize-nRevertPoints;i++)
					{
						points.push_back(cell((*(pos_it+i))-center,(*(strength_it+i)),NULL));
						set_bond_to(points.back(), (&points.back())-gap);
					}
				}
				else
				{
					for(unsigned int i=0;i<newSize-nRevertPoints;i++)
					{
						points.push_back(cell((*(pos_it+i))-center,sqrt(*(strength_it+i)),*(strength_it+i),&points.back()));
						set_bond_to(points.back(), (&points.back())-gap);
					}
				}
		{
			PDCoord lowest;
			PDCoord highest;
			PDCoord rebuild(0,0,0);
			getBoundingBox<PDCoord,PDFloat,PDCoord const*,PDFloat const*,dimension>(lowest,highest,newSize-nRevertPoints,&(*pos_it),&(*strength_it),0.0);


			for(int d=0;d<dimension;d++)
			{
				if(vertices.begin()->position[d]+center[d]-lowest[d]>rebuild[d])rebuild[d]=vertices.begin()->position[d]+center[d]-lowest[d];
				if(vertices[(1<<dimension)-1].position[d]+center[d]-highest[d]<rebuild[d])rebuild[d]=-(vertices[(1<<dimension)-1].position[d]+center[d]-highest[d]);
			}

			if(rebuild.squaredNorm()>0)
			{
				clearAllmyVertices();
					for(unsigned int vi=0;vi<_nVertices;++vi)
					{
						vertices[vi].invalid=0;
						vertices[vi].rrv=0;
					}

				{//createlike
						for(unsigned int i=0;i<nRevertPoints;++i)
							if(old_bond_ids[i] != kInvalidId) set_bond_to_id(points[i], old_bond_ids[i]);
						buildCube(vertices.begin()->position-2*rebuild,vertices[(1<<dimension)-1].position+2*rebuild);



				buildVertices(nRevertPoints);

				if(params.fill_myVertices||params.fill_neighbours)
					FillAllMyVertices(0,(1<<dimension));
				if(params.fill_neighbours)
					FillAllNeighbours();
				if(params.fill_zeroPoints)
					FillAllZeroPoints(nRevertZeros);
			}
				nRevertVertices=_nVertices;

			}
			else
			{
				if(points.data()!= oldstorage)
					{
							for(unsigned int vi=0;vi<_nVertices;++vi)
								for(int g=0;g<=dimension;++g)
								{
									const GeneratorRef& ref = old_vertex_generator_refs[vi][g];
									if(ref.is_valid())
									{
										if(ref.kind == GeneratorKind::point && ref.index < points.size())
										{
											vertices[vi].generators[g] = &points[ref.index];
											continue;
										}
										if(ref.kind == GeneratorKind::side && ref.index < sideGenerators.size())
										{
											vertices[vi].generators[g] = &sideGenerators[ref.index];
											continue;
										}
									}
								}

							for(std::size_t point_idx=0;point_idx<nRevertPoints;++point_idx)
							{
									cell& point = points[point_idx];
									if(point_idx < old_bond_ids.size() && old_bond_ids[point_idx] != kInvalidId) set_bond_to_id(point, old_bond_ids[point_idx]);
									point.neighbours.assign(old_neighbour_ids[point_idx].size(), NULL);
									point.neighboursIds = old_neighbour_ids[point_idx];
									for(std::size_t neighbour_idx=0;neighbour_idx<point.neighboursIds.size();++neighbour_idx)
										point.neighbours[neighbour_idx] = cell_ptr_from_id(point.neighboursIds[neighbour_idx]);
							}
					}
				}
		}


		buildVertices(newSize,nRevertPoints);

		if(params.fill_myVertices||params.fill_neighbours)
			FillAllMyVertices(nRevertPoints,nRevertVertices);
		if(params.fill_neighbours)
			FillAllNeighboursOfInvolved();
			if(params.fill_zeroPoints)
				FillAllZeroPoints(nRevertVertices,nRevertZeros);
			sync_all_link_mirrors();
		}

	void addMore(const PDCoord& pos,const PDFloat& radius,const int near)
	{
		// Single-point incremental insertion convenience wrapper.
		addMore(&pos,&radius,near);
	}

	static void make_inputfile(std::vector<PDCoord> const& position,std::vector<PDFloat> const& power)
	{
		for (unsigned int i = 0; i < position.size(); ++i)
		{
			std::cout << position[i].x() << " " <<  position[i].y() << " " << position[i].z() << " " <<power[i] << std::endl;
		}
	}

		void clearAllmyVertices()
		{
			// Drop per-cell cached vertex/zero-point ownership lists.
			for(cell& point : points)
			{
				clear_cell_my_vertices(point);
				point.myZeroPoints.clear();
			}
		}
	void buildCube(const PDCoord& lowest,const PDCoord& highest)
	{
		// Initialize the outer clipping cube and its connectivity as the starting polytope.
		_nVertices=1<<dimension;
		sideGenerators.clear();
		for(int i=0;i<2*dimension;i++)
			sideGenerators.push_back(cell(PDCoord(0,0,0),0,NULL));
		PDCoord lhc=lowest;
		vertices[0].setTo(lowest);
		for(int j=dimension-1;j>=0;j--)
			vertices[0].generators[j+1]=&sideGenerators[j];
		for(int i=0;i<(1<<dimension);i++)
			vertices[i].rrv=0;
		for(int i=0;i<(1<<dimension);i++)
			vertices[i].invalid=0;
		for(int i=1;i<(1<<dimension);i++)
		{
			int j=0;
			while(lhc[j]==highest[j])
			{
				lhc[j]=lowest[j];
				j++;
			}
			lhc[j]=highest[j];
			vertices[i].setTo(lhc);
			for(j=dimension-1;j>=0;j--)
				vertices[i].generators[j+1]=((lhc[j]==lowest[j])?&sideGenerators[j]:&sideGenerators[j+dimension]);
		}
		for(int i=0;i<(1<<dimension);i++)
			for(int d=0;d<dimension;d++)
			{
				const int ii=i;
				const int j=(ii>>d)%2?ii-(1<<d):ii+(1<<d);
				vertices[i].endPoints[d+1]=&vertices[j];
			}
		for(int i=0;i<(1<<dimension);i++)
		{//:TODO: direct sort or faster?
			for(int g=dimension-1;g>0;g--)
				for(int j=g;j>0;j--)
					if(&vertices[i].generators[j+1]-&vertices[i].generators[j]>0)
					{
						std::swap(vertices[i].generators[j],vertices[i].generators[j+1]);
						std::swap(vertices[i].endPoints[j],vertices[i].endPoints[j+1]);
					}
		}
		for(int d=1;d<=dimension;d++)//n
		{
			push_cell_my_vertex(*vertices[0].generators[d], vertices.data());
			push_cell_my_vertex(*vertices[(1<<dimension)-1].generators[d], vertices.data());
		}
	}
	void buildVertices(const unsigned int& nPoints,const int from=0)
	{
		// Insert points one-by-one and maintain a consistent power-diagram vertex network.
		//	try
		{

			if(points.size()>0)
			{
				maxr2=points[0].r2;
				for(int i=from;i< static_cast<int>(nPoints) ;i++)
				{
					if(maxr2<points[i].r2)
						maxr2=points[i].r2;
				}
				powerErr=1000*error(maxr2);
				if(from==0)
				{
					insertFirst();
				}


					for(unsigned int i=(from==0)?1:from;i<nPoints;i++)
					{
						unsigned int done=1;
						while(1)
						{
							if(doInsertion(prepareInsertion(points[i])))
								break;
							const PDFloat errorScale=insertionErrorScale;
							if(__power_diagram_internal_timing__){t3+=clock();t4+=clock();}

							cellPtr identicalPoint=NULL;
							done++;
							if(done>100)
								throw MyException();
							//is this an Identical point problem?
								{
									PDFloat mindist=(Involved[1]->position-Involved.front()->position).squaredNorm();
									cellPtr closest=Involved[1];
									for(std::size_t involved_idx=2;involved_idx<Involved.size();++involved_idx)
										if((Involved[involved_idx]->position-Involved.front()->position).squaredNorm()<mindist)
										{
											mindist=(Involved[involved_idx]->position-Involved.front()->position).squaredNorm();
											closest=Involved[involved_idx];
										}
										if(error(Involved.front()->r)>sqrt(mindist))
										{
										identicalPoint=closest;
										if(params.with_warnings)
										{
											std::cout<<"numerical similar point to "<<get_cell_id(*closest)+1<<" found. ";
											std::cout<<get_cell_id(*Involved.front())+1<<" is ignored"<<std::endl;
										}
									}
							}
							//delete new vertices built directly into vertices (when unused has been empty)
								if(_nUnused==0)
								{
									//here comes the deletion ;)
									_nVertices-=Involved.front()->myVertices.size()-unused.size()+_nUnused;
									for(vertexPtr involved_vertex : Involved.front()->myVertices)
										if(get_vertex_id(*involved_vertex)<(1<<dimension))
										{
											_nVertices++;
										}else involved_vertex->disconnect();
								}

								//delete new vertices built on unused (not needed because endpoints not constructed,yet)
								for(std::size_t unused_idx=_nUnused;unused_idx<unused.size();++unused_idx)
									unused[unused_idx]->disconnect();
								_nUnused=unused.size();	

							//reconnect replaced with persisting
								for(const vertexPtr replaced_vertex : Replaced)
								for(int endpoint_idx=replaced_vertex->isCorner();endpoint_idx<=dimension;++endpoint_idx)
				 				if(replaced_vertex->endPoints[endpoint_idx]->rrv<=0)
								{
									vertexPtr endpoint = replaced_vertex->endPoints[endpoint_idx];
									endpoint->rrv=0;
									for(int g1=replaced_vertex->isCorner();g1<=dimension;g1++)
									for(int g2=endpoint->isCorner();g2<=dimension;g2++)
									{
										if(replaced_vertex->generators[nth(0,g1)]==endpoint->generators[nth(0,g2)]&&replaced_vertex->generators[nth(1,g1)]==endpoint->generators[nth(1,g2)]&&replaced_vertex->generators[nth(2,g1)]==endpoint->generators[nth(2,g2)])	
										{
											replaced_vertex->endPoints[g1]=endpoint;
											endpoint->endPoints[g2]=replaced_vertex;
										}
									}
								}


								for(std::size_t involved_idx=1;involved_idx<Involved.size();++involved_idx)
									if(!Involved[involved_idx]->myVertices[0]->isConnected())
										Involved[involved_idx]->myVertices[0]=Replaced[0];


								//set everything zero again
								SetInvolvedPersistingVisitedToZero();
								clear_cell_my_vertices(*Involved.front());
							 	for(vertexPtr replaced_vertex : Replaced)
								{
				  					replaced_vertex->rrv=0;
									replaced_vertex->invalid=0;
								}

							Replaced.clear();
							if(identicalPoint!=NULL)
							{	
								push_cell_my_vertex(*Involved.front(), identicalPoint->myVertices.front());
								break;
							}
								else
								{
									const PDFloat oldr2=Involved.front()->r2;
									if(params.with_warnings)
										std::cout<<" Numerical Zero Warning: Power of "<<get_cell_id(*Involved.front())+1<<" is reduced from "<<Involved.front()->r2;
									Involved.front()->r2-=pow(2.0,done)*(errorScale);
								if(Involved.front()->r2>=0)	Involved.front()->r= sqrt(Involved.front()->r2);
								else 		Involved.front()->r=-sqrt(-Involved.front()->r2);
								if(params.with_warnings)
									std::cout<<" to "<<Involved.front()->r2<<" ( Change was "<<Involved.front()->r2-oldr2<<" )"<<std::endl;
							}




							if(done>100){std::cout<<"Exception : cannot get stable results"<<std::endl; throw MyException();}
						}
					}









if(!params.without_check){
					PDFloat checkconst=0;
//				std::cout<<"checking diagram"<<std::endl;
					for(unsigned int vi=0;vi<_nVertices;++vi)
						if(vertices[vi].isConnected())
							for(cell& point : points)
								if(point.power(vertices[vi].position)-vertices[vi].powerValue<-checkconst)
			if((&point!=vertices[vi].generators[0])&&(&point!=vertices[vi].generators[1])&&(&point!=vertices[vi].generators[2])&&(&point!=vertices[vi].generators[3]))
					{checkconst=point.power(vertices[vi].position)-vertices[vi].powerValue;
		
						std::cout<<"totaly wrong are "<<get_cell_id(point)<<" "<<point.power(vertices[vi].position)<<" "<<vertices[vi].generators[0]->power(vertices[vi].position)<<" "<<vertices[vi].generators[1]->power(vertices[vi].position)<<"          "<<vertices[vi].generators[0]->position[0]<<" "<<vertices[vi].generators[0]->position[1]<<" "<<vertices[vi].generators[0]->position[2]<<std::endl;
std::cout<<cell_id_or_invalid(vertices[vi].generators[0])<<" "<<cell_id_or_invalid(vertices[vi].generators[1])<<" "<<cell_id_or_invalid(vertices[vi].generators[2])<<" "<<cell_id_or_invalid(vertices[vi].generators[3])<<std::endl;
std::cout<<point.power(vertices[vi].endPoints[0]->position)<<" "<<point.power(vertices[vi].endPoints[1]->position)<<std::endl;
std::cout<<vertices[vi].generators[0]->power(vertices[vi].endPoints[0]->position)<<" "<<vertices[vi].generators[0]->power(vertices[vi].endPoints[1]->position)<<std::endl;
std::cout<<vertices[vi].generators[1]->power(vertices[vi].endPoints[0]->position)<<" "<<vertices[vi].generators[1]->power(vertices[vi].endPoints[1]->position)<<std::endl;
std::cout<<vertices[vi].position<<std::endl<<std::endl;
std::cout<<vertices[vi].endPoints[0]->position<<std::endl<<std::endl;
std::cout<<vertices[vi].endPoints[1]->position<<std::endl;
	throw MyException();

	}

if(std::abs(checkconst)>0.001)
	std::cout<<"the error of the worst vertex is around "<<checkconst<<std::endl;
//dump_vertices();
}










			}
			else{/*no vertices*/}
		}

			//delete waste that was produced, but during cleanup, never move waste (unvalidating "unused" pointer)
			std::sort(unused.begin(),unused.end());
			for(int i=static_cast<int>(unused.size())-1;i>=0;--i)
			{
				if(unused[i]!=&vertices[_nVertices-1])
					vertices[--_nVertices].moveAddressNetworkUpdateOnly(unused[i]);
				else {--_nVertices;}
			}

		unused.clear();
		_nUnused=0;

	}
	inline void findReplacedVertex(vertexPtr& This, PDFloat& value,const cell& insertionPoint)
	{
		// Starting from a hint vertex, descend the local graph to a vertex most clearly dominated by insertionPoint.
		if(value<0)
			return;
		PDFloat newValue;
		PDFloat smallVal=std::numeric_limits<PDFloat>::max();//something larger than value
		//look at the neighbours
			for(int idx=This->isCorner();idx<=dimension;++idx)
			{
				//each powerdiff value defines a plane between insertionPoint and current cell (generator0) approach the direction perpendicular to that plane in direction of insertion point !
				vertexPtr endpoint=This->endPoints[idx];
				newValue=endpoint->powerdiff3D(endpoint->generators[0],&insertionPoint);
			 	if(newValue<value)
				{
					value=newValue;
					This=endpoint;
					if(value<0)return;
					idx=This->isCorner()-1;
				}
				else if(newValue==value)
					smallVal=newValue;
			}
		if((smallVal!=value))
			return;
		smallVal=std::numeric_limits<PDFloat>::max();
		//found value is not definitely the best one
		//try hard to be sure not beeing in a local minimum (second neighbour)
			for(int g=This->isCorner();g<=dimension;++g)
				for(int g2=This->endPoints[g]->isCorner();g2<=dimension;++g2)
					if(This->endPoints[g]->endPoints[g2]!=This)
					{
						vertexPtr candidate=This->endPoints[g]->endPoints[g2];
						newValue=candidate->powerdiff3D(candidate->generators[0],&insertionPoint);
						if(newValue<value)
						{
							value=newValue;
							This=candidate;
							return findReplacedVertex(This,value,insertionPoint);
						}
						else if(newValue==value)
						smallVal=newValue;
				}
		if((smallVal!=value))
			return;
		smallVal=std::numeric_limits<PDFloat>::max();
		//second was also close... third neighbour...
			for(int g=This->isCorner();g<=dimension;++g)
				for(int g2=This->endPoints[g]->isCorner();g2<=dimension;++g2)
					if(This->endPoints[g]->endPoints[g2]!=This)
						for(int g3=This->endPoints[g]->endPoints[g2]->isCorner();g3<=dimension;++g3)
							if(This->endPoints[g]->endPoints[g2]->endPoints[g3]!=This->endPoints[g]&&This->endPoints[g]->endPoints[g2]->endPoints[g3]!=This)
							{
								vertexPtr candidate=This->endPoints[g]->endPoints[g2]->endPoints[g3];
								newValue=candidate->powerdiff3D(candidate->generators[0],&insertionPoint);
								if(newValue<value)
								{
									value=candidate->powerdiff3D(candidate->generators[0],&insertionPoint);
									This=candidate;
									return findReplacedVertex(This,value,insertionPoint);
								}else if(newValue==value)
									smallVal=newValue;
						}
		if(smallVal!=value)
			return;
		if(params.with_warnings)
			std::cout<<"warning : program slowed down because of too small accuracy"<<std::endl;
		//...so the numerical problem wants to be tough? A fat lot we care!
				for(unsigned int vi=0;vi<nVertices();++vi)
				if(vertices[vi].isConnected())
				{
					if(vertices[vi].powerdiff3D(vertices[vi].generators[0],&insertionPoint)<value)
				{
					value=vertices[vi].powerdiff3D(vertices[vi].generators[0],&insertionPoint);
						This=&vertices[vi];
					}
				}

		}
	inline ReplaceState finiteReplaced(vertex& This,const_cellPtr const& aCell)
	{
		// Classify whether vertex This is replaced by aCell, persists, or is numerically ambiguous.
		This.rrv=This.powerdiff3D(aCell,This.generators[0]);
			if(above_power_err(This.rrv)) return ReplaceState::replaced;
			if(below_neg_power_err(This.rrv)) return ReplaceState::persisting;
			This.rrv=0;
			return ReplaceState::ambiguous;
		}
	void dump_vertices(std::ostream& out=std::cout)
	{
		std::cout<<"vertices, generators and neighbours "<<std::endl;
		const auto generator_idx = [this](const_cellPtr generator) -> long long
		{
			const CellId id = cell_id_or_invalid(generator);
			return (id == kInvalidId) ? -1 : static_cast<long long>(id);
		};

		for(unsigned int vi=0;vi<_nVertices;++vi)
		{
			vertex& v = vertices[vi];
			if(v.isConnected() && (!v.isCorner()))
			{
				out<<generator_idx(v.generators[0])<<" "<<std::flush;
				out<<generator_idx(v.generators[1])<<" "<<std::flush;
				out<<generator_idx(v.generators[2])<<" "<<std::flush;
				out<<generator_idx(v.generators[3])<<"    "<<std::flush;
				out<<v.position[0]<<" "<<std::flush;
				out<<v.position[1]<<" "<<std::flush;
				out<<v.position[2]<<" "<<std::flush;
				out<<v.powerValue<<"   "<<std::flush;
				out<<v.generators[0]->power(v.position)<<std::endl;

				out<<" "<<generator_idx(v.endPoints[0]->generators[0])<<" "<<generator_idx(v.endPoints[0]->generators[1])<<" "<<generator_idx(v.endPoints[0]->generators[2])<<" "<<generator_idx(v.endPoints[0]->generators[3])<<std::endl;
				out<<" "<<generator_idx(v.endPoints[1]->generators[0])<<" "<<generator_idx(v.endPoints[1]->generators[1])<<" "<<generator_idx(v.endPoints[1]->generators[2])<<" "<<generator_idx(v.endPoints[1]->generators[3])<<std::endl;
				out<<" "<<generator_idx(v.endPoints[2]->generators[0])<<" "<<generator_idx(v.endPoints[2]->generators[1])<<" "<<generator_idx(v.endPoints[2]->generators[2])<<" "<<generator_idx(v.endPoints[2]->generators[3])<<std::endl;
				out<<" "<<generator_idx(v.endPoints[3]->generators[0])<<" "<<generator_idx(v.endPoints[3]->generators[1])<<" "<<generator_idx(v.endPoints[3]->generators[2])<<" "<<generator_idx(v.endPoints[3]->generators[3])<<std::endl;
			}
			else if(v.isCorner())
			{
				out<<generator_idx(v.generators[0])<<" "<<generator_idx(v.generators[1])<<" "<<generator_idx(v.generators[2])<<" "<<generator_idx(v.generators[3])<<"    "<<v.position[0]<<" "<<v.position[1]<<" "<<v.position[2]<<" "<<v.powerValue<<"   "<<v.generators[0]->power(v.position)<<std::endl;
			}
			else
			{
				std::cout<<"outtake"<<std::endl;
			}
			std::cout<<std::endl;
		}
	}

//	inline const int dimension()const{return dimension;}
	inline int nPoints()const{return points.size();}
	inline unsigned int const& nVertices() const {return _nVertices;}

		bool hasVirtualGenerators(const const_vertexPtr& that)const
		{
			const GeneratorRef& ref = that->generatorRefs[dimension];
			const const_cellPtr gptr = that->generators[dimension];
			if (ref.is_valid())
		{
			const const_cellPtr ref_ptr = cell_ptr_from_ref_const(ref);
			if (ref_ptr == gptr)
			{
					return !(ref.kind == GeneratorKind::point && ref.index < points.size());
				}
			}
			return (cell_id_or_invalid(gptr) == kInvalidId);
		}
		int nVirtualGenerators(const const_vertexPtr& that)const
		{
			const GeneratorRef& ref = that->generatorRefs[dimension];
			const const_cellPtr gptr = that->generators[dimension];
			bool dim_is_real = (cell_id_or_invalid(gptr) != kInvalidId);
			if (ref.is_valid())
			{
				const const_cellPtr ref_ptr = cell_ptr_from_ref_const(ref);
				if (ref_ptr == gptr)
			{
				dim_is_real = (ref.kind == GeneratorKind::point && ref.index < points.size());
			}
		}
		if(!dim_is_real)
			if(!dim_is_real)
				if(that->isCorner())
					return 3;
				else return 2;
			else	return 1;
		else	return 0;
	}
	cell const *findCellInsideCube(const PDCoord& pos,cell const * hint=NULL)
	{
		// Greedy neighbor walk to find the cell with minimal power at pos inside the current cube.
		if(hint==NULL)hint =&points[points.size()/2];
			if(!hint->neighboursIds.empty())
			{
				for(const CellId neighbour_id : hint->neighboursIds)
				{
					const_cellPtr neighbour = cell_ptr_from_id_const(neighbour_id);
					if(neighbour != NULL && neighbour->power(pos)<hint->power(pos))return findCellInsideCube(pos,neighbour);
				}
			}
			else
			{
				for(const cellPtr neighbour : hint->neighbours)
				{
					if(neighbour->power(pos)<hint->power(pos))return findCellInsideCube(pos,neighbour);
				}
			}
			return hint;
		}
private:
	vertexPtr getRepresentative(const const_cellPtr This)
	{
		// Retrieve a connected representative vertex for This, falling back along bondTo if needed.
		if(This == NULL)
		{
			return vertices.data();
		}
		const std::size_t n_vertices = This->myVerticesIds.empty() ? This->myVertices.size() : This->myVerticesIds.size();
		for(std::size_t i = 0; i < n_vertices; ++i)
		{
			vertexPtr candidate = NULL;
			if(!This->myVerticesIds.empty())
			{
				candidate = vertex_ptr_from_id(This->myVerticesIds[i]);
			}
			if(candidate == NULL && i < This->myVertices.size())
			{
				candidate = This->myVertices[i];
			}
			if(candidate != NULL && candidate->isConnected() && candidate->hasGenerator(This))
				{
					return candidate;
				}
		}
		if(This->bondToId != kInvalidId)
		{
			cellPtr const bond = cell_ptr_from_id(This->bondToId);
			if(bond != NULL && bond != This) return getRepresentative(bond);
		}
		const CellId this_id = cell_id_or_invalid(This);
		if(this_id != 0 && This->bondTo != NULL && This->bondTo != This) return getRepresentative(This->bondTo);
		else
		{
			return vertices.data();
		}
	}
	vertexPtr prepareInsertion(cell & This,vertexPtr hint=NULL)
	{
		// Build replaced/persisting/involved sets for inserting This, including numerical fallback reductions.
		if(__power_diagram_internal_timing__)t2-=clock();
			//there is a power of new cell that is so low, that only one vertex would be replaced. *hint will be the one
			hint=getRepresentative((This.bondToId != kInvalidId) ? cell_ptr_from_id(This.bondToId) : This.bondTo);
			PDFloat value=hint->powerdiff3D(hint->generators[0],&This);
			findReplacedVertex(hint,value,This);
			if(__power_diagram_internal_timing__)t2+=clock();

			unsigned int done=1;
			while(1)
			{
				if(done!=1)
				{
					value=hint->powerdiff3D(hint->generators[0],&This);
					findReplacedVertex(hint,value,This);
				}

				if(FillReplacedPersistingAndInvolved(This,hint))
					break;

				const PDFloat oldr2=Involved.front()->r2;
					if(params.with_warnings)
						std::cout<<"Numerical Warning: Power of "<<get_cell_id(*Involved.front())+1<<" is reduced from "<<Involved.front()->r2;
					SetInvolvedPersistingVisitedToZero();
					clear_cell_my_vertices(This);
						for(const vertexPtr replaced_vertex : Replaced)
						{
						replaced_vertex->rrv=0;
						for(int g=replaced_vertex->isCorner();g<=dimension;g++)
							replaced_vertex->endPoints[g]->rrv=0;
					}
				Replaced.clear();
				This.r2-=pow(2.0,done)*(PowerDiagram<PDFloat,PDCoord,dimension>::powerErr);
				if(This.r2>0)	This.r= sqrt(This.r2);
				else 		This.r=-sqrt(-This.r2);
				if(params.with_warnings)
					std::cout<<" to "<<Involved.front()->r2<<" ( Change was "<<Involved.front()->r2-oldr2<<" )"<<std::endl;
				done++;
				if(done>100){std::cout<<"exception : cannot get stable results with atom "<<get_cell_id(*Involved.front())<<" "<<Involved.front()->position+center<<std::endl;
					throw MyException();}
			}
			return hint;
		}

	bool doInsertion(const vertexPtr& hint)
	{
		// Materialize insertion after prepareInsertion(): create new finite vertices, connect, and update caches.
//			if(hint!=NULL)
				{
						if(__power_diagram_internal_timing__){const unsigned int zeit=clock();t3-=zeit;t4-=zeit;}
					if(!CreateFiniteVerticesFromReplaced())
						return false;
						if(__power_diagram_internal_timing__){const unsigned int zeit=clock();t4+=zeit;t5-=zeit;}
					ConnectNewFinitesAmongThemselves3D();
						if(__power_diagram_internal_timing__)t5+=clock();
					UpdateUnused();
					AssignRepresentativeVerticesToCells(hint);
					SetInvolvedPersistingVisitedToZero();
						if(__power_diagram_internal_timing__)t3+=clock();
				}
			return true;
		}

	void clear_interna()
	{
		// Reset per-insertion transient containers.
		Replaced.clear();
		Involved.clear();
//		planes.clear();//should always be clean
	}

		inline void insertFirst()
		{
		// Seed the diagram: assign the first real generator to all cube corners.
		clear_interna();
		for(int i=0;i<(1<<dimension);i++)
			vertices[i].setPowerData(points.data());
		for(int i=0;i<(1<<dimension);i++)
			vertices[i].generators[0]=points.data();
			for(int i=0;i<(1<<dimension);i++)
				push_cell_my_vertex(points[0], &vertices[i]);
		}
	void FillAllMyVertices(const int fromPoint=0,const int fromVertex=1<<dimension)
	{
		// Recompute per-cell vertex ownership lists from the current global vertex array.
		{
				if(fromPoint>0)
					Involved.clear();
				for(std::size_t point_idx=static_cast<std::size_t>(fromPoint);point_idx<points.size();++point_idx)
					clear_cell_my_vertices(points[point_idx]);


			for(unsigned int vi=fromVertex;vi<_nVertices;++vi)
			if(!(vertices[vi].invalid))
			{
					if(!(hasVirtualGenerators(&vertices[vi])))
						for(int g=0;g<=dimension;++g)
						{
							cellPtr generator = vertices[vi].generators[g];
							push_cell_my_vertex(*generator, &vertices[vi]);
							if(fromPoint>0)
								if(generator->visitedAs==0)
							{
								generator->visitedAs=-1;
								Involved.push_back(generator);
							}
						}
						else if(!vertices[vi].isCorner())
							for(int g=0;g<=dimension;++g)
							{
									cellPtr generator = vertices[vi].generators[g];
									if(generator_ref_or_invalid(generator).kind == GeneratorKind::side) break;
									push_cell_my_vertex(*generator, &vertices[vi]);
									if(fromPoint>0)
										if(generator->visitedAs==0)
								{
									generator->visitedAs=-1;
									Involved.push_back(generator);
								}
						}
					else
					{/*dont give corners to sasa code, it doesnt check for it... so we define myVertices as not holding corners!*/
						std::cout<<"wrong internal order, SASA stopped"<<std::endl;
						throw MyException();
					}

		}
		}
	}

		void FillAllNeighbours()
		{
			// Recompute full cell adjacency from shared finite vertices.
			for(cell& point : points)
			{
				clear_cell_neighbours(point);
				point.visitedAs=-1;
			}
			for(std::size_t point_idx=0;point_idx<points.size();++point_idx)
			{
				cell& point = points[point_idx];
				int current_cell_order = static_cast<int>(point_idx);
				const CellId current_cell_id = cell_id_or_invalid(&point);
				if(current_cell_id != kInvalidId) current_cell_order = static_cast<int>(current_cell_id);
				if(!point.myVerticesIds.empty())
				{
					for(const VertexId vid : point.myVerticesIds)
					{
						vertexPtr vtx = vertex_ptr_from_id(vid);
						if(vtx == NULL || vtx->isCorner()) continue;
							for(int g=dimension;g>=0;g--)
								if(vtx->generators[g]->isReal(*this))
								{
									if(vtx->generators[g]->visitedAs<current_cell_order&&vtx->generators[g]!=&point)
									{
										push_cell_neighbour(point, vtx->generators[g]);
										vtx->generators[g]->visitedAs=current_cell_order;
									}
								}
						}
				}
				else
				{
					for(typename std::vector<vertexPtr>::const_iterator it2=point.myVertices.begin();it2!=point.myVertices.end();++it2)
						if(!(*it2)->isCorner())
								for(int g=dimension;g>=0;g--)
									if((*it2)->generators[g]->isReal(*this))
									{
										if((*it2)->generators[g]->visitedAs<current_cell_order&&(*it2)->generators[g]!=&point)
											{
												push_cell_neighbour(point, (*it2)->generators[g]);
												(*it2)->generators[g]->visitedAs=current_cell_order;
											}
									}
					}
			}

		for(cell& point : points)
			point.visitedAs=0;
	}

			void FillAllNeighboursOfInvolved()
			{
				// Incrementally refresh adjacency only for cells touched by recent insertion/revert operations.
				sort(Involved.begin(),Involved.end(),
					[this](const cellPtr a, const cellPtr b)
					{
						const CellId a_id = cell_id_or_invalid(a);
						const CellId b_id = cell_id_or_invalid(b);
						if(a_id != kInvalidId && b_id != kInvalidId) return a_id < b_id;
						return a < b;
					});
				for(cellPtr involved_cell : Involved)
				{
				std::size_t neighbour_idx=0;
					while(neighbour_idx<involved_cell->neighbours.size())
					{
						cellPtr neighbour = involved_cell->neighbours[neighbour_idx];
						if(!involved_cell->neighboursIds.empty() && neighbour_idx < involved_cell->neighboursIds.size())
						{
							cellPtr id_neighbour = cell_ptr_from_id(involved_cell->neighboursIds[neighbour_idx]);
							if(id_neighbour != NULL) neighbour = id_neighbour;
						}
						if(neighbour == NULL || neighbour->visitedAs==-1)
							erase_cell_neighbour(*involved_cell, neighbour_idx);
						else
							++neighbour_idx;
					}
			}
			for(const cellPtr involved_cell : Involved)
			{
				int current_cell_order = 0;
				const CellId current_cell_id = cell_id_or_invalid(involved_cell);
				if(current_cell_id != kInvalidId)
					current_cell_order = static_cast<int>(current_cell_id);
				else
				{
					for(std::size_t point_idx=0;point_idx<points.size();++point_idx)
						if(&points[point_idx]==involved_cell)
						{
							current_cell_order = static_cast<int>(point_idx);
							break;
						}
				}
				if(!involved_cell->myVerticesIds.empty())
				{
					for(const VertexId vid : involved_cell->myVerticesIds)
					{
						vertexPtr vtx = vertex_ptr_from_id(vid);
						if(vtx == NULL) continue;
						for(int g=dimension;g>=0;g--)
							if(vtx->generators[g]->isReal(*this))
							{
								if(vtx->generators[g]->visitedAs!=0&&vtx->generators[g]->visitedAs<=current_cell_order&&vtx->generators[g]!=involved_cell)
								{
									push_cell_neighbour(*involved_cell, vtx->generators[g]);
									vtx->generators[g]->visitedAs=current_cell_order+1;
								}
							}
					}
				}
				else
				{
					for(typename std::vector<vertexPtr>::const_iterator it2=involved_cell->myVertices.begin();it2!=involved_cell->myVertices.end();++it2)
					{
						for(int g=dimension;g>=0;g--)
							if((*it2)->generators[g]->isReal(*this))
							{
									if((*it2)->generators[g]->visitedAs!=0&&(*it2)->generators[g]->visitedAs<=current_cell_order&&(*it2)->generators[g]!=involved_cell)
									{
										push_cell_neighbour(*involved_cell, (*it2)->generators[g]);
										(*it2)->generators[g]->visitedAs=current_cell_order+1;
									}
							}
					}
				}
			}

		for(cellPtr involved_cell : Involved)
			involved_cell->visitedAs=0;
	}
	void FillAllZeroPoints(	unsigned int fromVertex=(1<<dimension),const unsigned int fromZero=0)
	{
		// Recompute intersections where power value crosses zero along diagram edges.
		zeros.erase(zeros.begin()+fromZero,zeros.end());
		for(unsigned int vertex_index=fromVertex;vertex_index<this->_nVertices;++vertex_index)
		{
			const_vertexPtr it=&vertices[vertex_index];
			if(!(it->invalid))
				if(it->generators[dimension-1]->isReal(*this))
					for(typename std::array <vertexPtr,dimension+1>::const_iterator it2=it->endPoints.begin()+(hasVirtualGenerators(it))*3;it2!=it->endPoints.end();++it2)
					{
						if((*it2)-&(*it)>0)
						{
							if(it->powerValue>0)
							{
								const int branch=it2-it->endPoints.begin();
								const PDFloat& v3=(*it2)->powerValue;
								const PDFloat& v2=it->powerValue;
								const PDFloat& v1=it->generators[branch==0]->power(2*it->position-(*it2)->position);
								const PDFloat quot=2*(v1+v3-2*v2);
								//const PDFloat rootsq=sqr(v1-v3)-4*quot*v2;
								const PDFloat rootsq=(v1-v3)*(v1-v3)-4*quot*v2;
								if(rootsq<=0)
									continue;
								if(below_power_err(quot))
									if(v1>=0&&v2>=0&&v3>=0)continue;
								const PDFloat rootquot=sqrt(rootsq)/quot;
								const PDFloat min=(v1-v3)/quot;
								const PDFloat sol1=min+rootquot;
								const PDFloat sol2=min-rootquot;
								if(sol1>0&&sol1<1)
									if((*it2)->powerValue>0)
									{
											zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol1,&vertices[vertex_index],branch));
											zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol2,&vertices[vertex_index],branch));
										}
										else
											zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol1,&vertices[vertex_index],branch));
									else if(sol2>0&&sol2<1)
										zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol2,&vertices[vertex_index],branch));
									else
									{//the covered zeros
										zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol1,&vertices[vertex_index],branch));
										zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol2,&vertices[vertex_index],branch));
									}
								}else if((*it2)->powerValue>0)
							{
								const int branch=it2-it->endPoints.begin();
								const PDFloat& v3=(*it2)->powerValue;
								const PDFloat& v2=it->powerValue;
								const PDFloat& v1=it->generators[branch==0]->power(2*it->position-(*it2)->position);
								const PDFloat quot=2*(v1+v3-2*v2);
								//const PDFloat rootsq=sqr(v1-v3)-4*quot*v2;
								const PDFloat rootsq=(v1-v3)*(v1-v3)-4*quot*v2;
								if(rootsq<=0)
									continue;
								if(below_power_err(quot))
									if(v1>=0&&v2>=0&&v3>=0)continue;
								const PDFloat rootquot=sqrt(rootsq)/quot;
								const PDFloat min=(v1-v3)/quot;
								const PDFloat sol1=min+rootquot;
								const PDFloat sol2=min-rootquot;
								if(sol1>0&&sol1<1)
										zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol1,&vertices[vertex_index],branch));
									else
										zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol2,&vertices[vertex_index],branch));
								}
							}
						}
                    /*
					else //vertex on edge of cube dont have connections through spheres
					{

					}
                    */

		}
			for(unsigned int i=fromZero;i<zeros.size();i++)
			{
				sync_zero_link_mirrors(zeros[i]);
				for(int g = 0; g < dimension; ++g)
				{
					const GeneratorRef& ref = zeros[i].generatorRefs[g];
					if(ref.is_valid() && ref.kind == GeneratorKind::point && ref.index < points.size())
					{
						points[ref.index].myZeroPoints.push_back(i);
					}
					else if(zeros[i].generators[g] != NULL)
					{
						zeros[i].generators[g]->myZeroPoints.push_back(i);
					}
				}
			}
		}

		inline bool tryToBuildVertexOnEdge(const const_vertexPtr& This,const int& here)//,const cellPtr s1, const cellPtr s2,const cellPtr s3,const PDCoord& direction);
		{
			// Create one new finite vertex on a surviving edge between replaced and persisting regions.
			//edge between This (replaced and finite) and that defined by generators s1,s2,s3 will get a vertex (of newest,s1,s2,s3)
			{
				vertexPtr builtVertex = NULL;
				if(_nUnused==0)
				{
					if(nVertices()==vertices.capacity())
						throw MyException();
					vertices[_nVertices].endPointsAndPositionOverwrite(This->endPoints[here],This->getPowerPointOnLine2(This->endPoints[here]));
					builtVertex = &vertices[++_nVertices-1];
				}
				else
				{
					unused[_nUnused-1]->endPointsAndPositionOverwrite(This->endPoints[here],This->getPowerPointOnLine2(This->endPoints[here]));
					builtVertex = unused[--_nUnused];
				}
				if(!builtVertex->Init(This,here,*this))
				{
					insertionErrorScale=powerErr;
					return false;
				}
			}

			return true;
		}


	//  void replace_a_vertex(vertex& old_vertex,const Cell& newGenerator);
	//  void checkvertex(vertexIter& myvertex,Cell& newGenerator,std::vector<vertexIter>& replaced,std::vector <vertexIter>&surroundings,const vertexIter former);

	bool FillReplacedPersistingAndInvolved(cell& This,vertexPtr start)
	{
		// Flood-fill from start to identify replaced vertices and all cells involved by insertion of This.
		clear_interna();
			Involved.push_back(&This);
			const ReplaceState startState = finiteReplaced(*start,Involved.front());
			if(startState==ReplaceState::ambiguous) return false;
			if(startState==ReplaceState::replaced)
			{
				if(start->isCorner())
				{
					return start->cornerToReplacedAndGo(*this);
				}
				else
					return start->finiteToReplacedAndGo(*this);
			}
			return true;
		}
	bool CreateFiniteVerticesFromReplaced()
	{
		// For every replaced vertex, build the new finite vertices that lie on cut edges.
		//best procedure for new vertices : knowledge : each new (finite) vertex MUST lie on
		//exactly one old EXISTING edge which is NOT disappearing totally
		//all possible edges are the ones coming out our "replaced" vertices
		//so we only try to create if an endPoint of a replaced vertex is not replaced (visitedAs ==-1)
				for(const vertexPtr replaced_vertex : Replaced)
				{
					int needed=0;
					for(int g=dimension;g>=replaced_vertex->isCorner();g--)
						if(replaced_vertex->endPoints[g]->rrv<=0)
							needed++;
					const int additionalNeeded=needed-static_cast<int>(_nUnused);
					if(additionalNeeded>0&&nVertices()+additionalNeeded>vertices.capacity())
				{
					this->ReserveNewVertices();
				}
				bool ok=true;
				if(!replaced_vertex->isCorner())
					ok=replaced_vertex->template buildIn<0>(this);
				else
					ok=replaced_vertex->template buildIn<1>(this);
				if(!ok)
				{
					return false;
				}
			}
			return true;
		}

	inline void ConnectNewFinitesAmongThemselves3D()
	{
		// Wire newly created finite vertices using shared old-generator pairs as edge keys.
		//we use an InvolvedSize*InvolvedSize-matrix (sparse) and fill in all new edges
		//these are generated by the new Cell and two older Cells
		//we identify vertices to be connected over an new edge by only the two old cells! (new one is everywhere)
		//the earlier dumps itself into the array at spot [a,b] (with a<b) the later simply picks it up!
		if(Involved.size()*Involved.size()>planes.size())
			planes.resize(Involved.size()*Involved.size());

			for(unsigned int vertex_idx=0;vertex_idx<(1u<<dimension);++vertex_idx)
				if(vertices[vertex_idx].rrv>0)
				{
					vertices[vertex_idx].setPowerData(Involved.front());
					vertices[vertex_idx].generators[0]=(Involved.front());
				}


			for(const vertexPtr involved_vertex : Involved.front()->myVertices)
			//if(!(*it)->isCorner())
			{	
				involved_vertex->registerForConnection3D(this);
			}

	}

	


	// const int GoAlongFace(const vertexIter& former, const vertexIter& current,const vertexIter& finish,const CellPtr& Generator1,const Ptr& Generator2,const int function(const int&))const;
	void ReserveNewVertices()
	{
		// Grow vertex storage while preserving all pointer-based topology links.
		std::vector<int> _replaced;
		std::vector<int> _unused;
		std::vector<int> _currentmyVertices;
		std::vector<int> _first;
		const_vertexPtr const oldMemoryPointer=vertices.data();
			_replaced.reserve(Replaced.size());
			_unused.reserve(unused.size());
			_currentmyVertices.reserve(Involved.front()->myVertices.size());
			_first.reserve(vertices.capacity());

			for(const vertexPtr replaced_vertex : Replaced)
				_replaced.push_back(static_cast<int>(get_vertex_id(*replaced_vertex)));
			for(const vertexPtr unused_vertex : unused)
				_unused.push_back(static_cast<int>(get_vertex_id(*unused_vertex)));
				if(!Involved.front()->myVerticesIds.empty())
				{
					for(const VertexId vid : Involved.front()->myVerticesIds)
						_currentmyVertices.push_back(static_cast<int>(vid));
				}
				else
				{
					for(const vertexPtr current_vertex : Involved.front()->myVertices)
						_currentmyVertices.push_back(static_cast<int>(get_vertex_id(*current_vertex)));
				}
				std::size_t involved_prefix = 0;
				const CellId involved_front_id = cell_id_or_invalid(Involved.front());
				if(involved_front_id != kInvalidId) involved_prefix = involved_front_id;
				else
					for(std::size_t point_idx=0;point_idx<points.size();++point_idx)
						if(&points[point_idx]==Involved.front())
						{
							involved_prefix = point_idx;
							break;
						}
				for(std::size_t point_idx=0;point_idx<involved_prefix;++point_idx)
				{
					const cell& point = points[point_idx];
					if(!point.myVerticesIds.empty()) _first.push_back(static_cast<int>(point.myVerticesIds.front()));
					else _first.push_back(static_cast<int>(get_vertex_id(*point.myVertices.front())));
				}

		vertices.reserve(2*vertices.capacity()+1);

		for(unsigned int i=0;i<_nVertices;++i)
			vertices[i].refreshAfterRealloc(oldMemoryPointer+i);
			vertices.resize(vertices.capacity());

			for(unsigned int i=0;i<Replaced.size();i++)
			{
				vertexPtr restored = NULL;
				if(_replaced[i] >= 0) restored = vertex_ptr_from_id(static_cast<VertexId>(_replaced[i]));
				if(restored == NULL) restored = vertices.data()+_replaced[i];
				Replaced[i]=restored;
			}
			for(unsigned int i=0;i<unused.size();i++)
			{
				vertexPtr restored = NULL;
				if(_unused[i] >= 0) restored = vertex_ptr_from_id(static_cast<VertexId>(_unused[i]));
				if(restored == NULL) restored = vertices.data()+_unused[i];
				unused[i]=restored;
			}
				for(unsigned int i=0;i<Involved.front()->myVertices.size();i++)
				{
					vertexPtr restored = NULL;
					if(_currentmyVertices[i] >= 0) restored = vertex_ptr_from_id(static_cast<VertexId>(_currentmyVertices[i]));
					if(restored == NULL) restored = vertices.data()+_currentmyVertices[i];
					set_cell_my_vertex(*Involved.front(), i, restored);
				}
				for(std::size_t i=0;i<involved_prefix;i++)
				{
					vertexPtr restored = NULL;
					if(_first[i] >= 0) restored = vertex_ptr_from_id(static_cast<VertexId>(_first[i]));
					if(restored == NULL) restored = vertices.data()+_first[i];
					set_cell_my_vertex(points[i], 0, restored);
				}

	}
	void UpdateUnused()
	{
		// Mark replaced vertices as invalid/unused and update per-cell myVertices ownership lists.
		unused.resize(_nUnused);
		//mark replaced as unused
		for(std::size_t replaced_idx=0;replaced_idx<Replaced.size();)
		{
			vertexPtr replaced_vertex=Replaced[replaced_idx];
			if(replaced_vertex->isCorner())
			{
				replaced_vertex->rrv=0;
				//Replaced.erase(it);--it;//the corners are always part of diagram	
				Replaced[replaced_idx]=Replaced.back();
				Replaced.pop_back();
			}
			else
			{
				if(get_vertex_id(*replaced_vertex)<nRevertVertices)
					Invalids.push_back(replaced_vertex);
				replaced_vertex->disconnect();//vertex has no connection any more.  we delete it later
				++replaced_idx;
			}
		}
		if(nRevertVertices==0)
			for(vertexPtr replaced_vertex : Replaced)
				unused.push_back(replaced_vertex);

		else
			for(vertexPtr replaced_vertex : Replaced)
				if(get_vertex_id(*replaced_vertex)>=nRevertVertices)
					unused.push_back(replaced_vertex);
					else
						for(const cellPtr generator : replaced_vertex->generators)
							for(unsigned int i=0;i<generator->myVertices.size();i++)
								if(generator->myVertices[i]==replaced_vertex)
								{
									erase_cell_my_vertex(*generator, i);
									break;
								}
				_nUnused=unused.size();
		}
			void AssignRepresentativeVerticesToCells(const vertexPtr& aDefault)
		{
			// Ensure each involved cell keeps at least one representative connected vertex after insertion.
				//if there are no new vertices, the new cell is covered
				if(Involved.front()->myVertices.empty())
					push_cell_my_vertex(*Involved.front(), aDefault);
				else// we need one existing vertex close to each cell => we give every cell without representativ a new vertex
				{
					vertexPtr new_representative = NULL;
					if(!Involved.front()->myVerticesIds.empty())
					{
						new_representative = vertex_ptr_from_id(Involved.front()->myVerticesIds.front());
					}
					if(new_representative == NULL)
					{
						new_representative = Involved.front()->myVertices.front();
					}
					for(std::size_t involved_idx=1;involved_idx<Involved.size();++involved_idx)
					{
						cellPtr involved_cell = Involved[involved_idx];
						if(cell_id_or_invalid(involved_cell) == kInvalidId) continue;
						vertexPtr representative = NULL;
						if(!involved_cell->myVerticesIds.empty())
						{
							representative = vertex_ptr_from_id(involved_cell->myVerticesIds.front());
						}
						if(representative == NULL)
						{
							if(involved_cell->myVertices.empty()) continue;
							representative = involved_cell->myVertices.front();
						}
						if(!representative->isConnected())//if representing vertex has been erased
							set_cell_my_vertex(*involved_cell, 0, new_representative);//we assign representative of new also to this one
					}
				}

		}

		void SetInvolvedPersistingVisitedToZero() 
		{
			// Clear temporary rrv/visited marks used during local insertion traversal.
			for(std::size_t involved_idx=1;involved_idx<Involved.size();++involved_idx)
				Involved[involved_idx]->visitedAs=0;
			if(!Involved.front()->myVerticesIds.empty())
			{
				for(const VertexId vid : Involved.front()->myVerticesIds)
				{
					vertexPtr vtx = vertex_ptr_from_id(vid);
					if(vtx == NULL) continue;
					vtx->rrv=0;
					if(!vtx->isCorner())
						vtx->endPoints[0]->rrv=0;
					else
					{
						vtx->endPoints[1]->rrv=0;
						vtx->endPoints[2]->rrv=0;
						vtx->endPoints[3]->rrv=0;
					}
				}
			}
			else
			{
				for(const vertexPtr vtx : Involved.front()->myVertices)
				{
					vtx->rrv=0;
					if(!vtx->isCorner())
						vtx->endPoints[0]->rrv=0;
					else
					{
						vtx->endPoints[1]->rrv=0;
						vtx->endPoints[2]->rrv=0;
						vtx->endPoints[3]->rrv=0;
					}
				}
			}
		}

public:
template <class VectorSubtraction>
    inline /*static*/ PDCoord intersectionOfLineAndPlane3D(const VectorSubtraction& direction,const VectorSubtraction& supportVec,const VectorSubtraction& normal,const PDFloat& planeVal)
	{
		// Compute intersection of a line and power-bisector plane in shifted coordinates.
		const PDFloat tmp=(normal.dot(direction));
//		const PDFloat sqr=normal.squaredNorm();

        if(tmp>power_err_scaled_epsilon())
            return direction*((0.5*(planeVal+normal.squaredNorm())-normal.dot(supportVec))/tmp);
        throw MyException();
		}


inline const PDCoord getPowerCenterOf2(const cell *const g0,const cell *const g1)
{
	// Weighted midpoint on the bisector of two generators in power metric.
	//0.5*(1+(Ra²-Rb²)/dist2)*(a_pos-b_pos)+a_pos
	return (0.5*(1.+(g0->r2-g1->r2)/(g1->position-g0->position).squaredNorm()))*(g1->position-g0->position)+g0->position;
}

inline PDCoord getPowerPointOnLine(const PDCoord& direction,const PDCoord& supportVector,cell const* const& a, cell const* const& b)
{
	// Intersect a line with the power bisector plane of generators a and b.
//	const PDCoord PlaneNormal=(b->position-a->position)/*/(a->position-b->position).norm()*/;
//	const PDFloat PlaneValue=0.5*(PlaneNormal.squaredNorm()+(a->r2-b->r2)/*(a->position-b->position).norm()*/);
	//PlaneNormal and PlaneValue are a factor of (a->position-b->position).norm() too big but they cancel each other out
	return (PowerDiagram<PDFloat,PDCoord,3>::intersectionOfLineAndPlane3D(direction-supportVector,supportVector-a->position,(b->position-a->position),a->r2-b->r2)+supportVector);
}



	struct zeroPoint
	{
		PDFloat pos;
		vertexPtr from;
		VertexId fromId;
		int branch;
		std::array<cellPtr,dimension> generators;
		std::array<GeneratorRef,dimension> generatorRefs;

		zeroPoint(const cellPtr& a,const cellPtr& b,const cellPtr& c,const PDFloat& position,const vertexPtr& origin,const int& way):
			pos(position),from(origin),fromId(kInvalidId),branch(way)
		{
			generators[0]=a;
//		generators[0]->myZeroPoints.push_back(this);
		generators[1]=b;
//		generators[1]->myZeroPoints.push_back(this);
		generators[2]=c;
//		generators[2]->myZeroPoints.push_back(this);

	}
	PDCoord getPos()const
	{
		// Evaluate the 3D position of this zero crossing along its source edge.
		return (from->endPoints[branch]->position)*pos-from->position*(pos-1);
	}
	bool isValid()const{return ((!from->invalid)&&(!from->endPoints[branch]->invalid));}
};
	struct vertex
	{
		PDFloat rrv;//relative replace value (power difference)
		bool invalid;
		std::array<cellPtr,dimension+1> generators;
		std::array<GeneratorRef,dimension+1> generatorRefs;
		PDCoord position;
		PDFloat powerValue;
		std::array <vertexPtr,dimension+1> endPoints;
		std::array<VertexId,dimension+1> endPointIds;


	friend class PowerDiagram <PDFloat, PDCoord,dimension>;
	inline bool isCorner() const {return endPoints[0]==NULL; }
	inline bool isOnEdge(const PowerDiagram<PDFloat,PDCoord,dimension>& This)
	{
		return (!this->generators[dimension-2]->isReal(This));
	}
	inline bool isOnSurface(const PowerDiagram<PDFloat,PDCoord,dimension>& This)
	{
		return (!this->generators[dimension-1]->isReal(This));
	}
//	inline int hasVirtualGenerators()const {return (generators[dimension]->id<0);}
	inline int hasGenerator(const const_cellPtr& that)const {return (generators[0]==that||generators[1]==that||generators[2]==that||generators[3]==that);}
//	inline int isFinite() const { return generators[0]!=NULL; }
	inline void disconnect(){invalid=1;}
	inline int isConnected()const{return !invalid;}

	//  vertex(const vertex& copy);
		inline vertex():invalid(1)/*,generators(dimension+1,NULL),endPoints(dimension+1,NULL)*/
		{
			for (int g = 0; g <= dimension; ++g) endPointIds[g] = kInvalidId;
		}

			inline bool Init(const const_vertexPtr& This,const int& keep,PowerDiagram<PDFloat,PDCoord,dimension>& owner)
		{
			// Initialize a newly created finite vertex generated by cutting one edge.
				this->setPowerData(owner.Involved.front());

			for(int g=dimension;g>0;g--)
				generators[g]=This->generators[g-(g<=keep)];
			generators[0]=owner.Involved.front();
			owner.push_cell_my_vertex(*owner.Involved.front(), this);
			endPoints[0]->fastWhichis(This)=this;

			if(owner.within_power_err(powerValue))
			{
				return false;
			}
				return true;
		}

inline PDCoord getPowerPointOnLine2(vertex const* const& persist)const
{
	// Linear interpolation point where the replacement score crosses zero on an edge.
//	const PDCoord PlaneNormal=(b->position-a->position)/*/(a->position-b->position).norm()*/;
//	const PDFloat PlaneValue=0.5*(PlaneNormal.squaredNorm()+(a->r2-b->r2)/*(a->position-b->position).norm()*/);
	//PlaneNormal and PlaneValue are a factor of (a->position-b->position).norm() too big but they cancel each other out
//	const PDFloat lowPower=newOne->power(replaced->position);
//	const PDFloat persistPower=newOne->power(persist->position);
	return ((rrv)/((rrv)-(persist->rrv)))*(persist->position-position)+position;
}

		void operator=(const vertex& that)
		{
			generators=that.generators;
			generatorRefs=that.generatorRefs;
			position=that.position;
			powerValue=that.powerValue;
			endPoints=that.endPoints;
			endPointIds=that.endPointIds;
			invalid=that.invalid;
			rrv=that.rrv;
		}
private :
	inline void setPowerData(const const_cellPtr& aCell)
	{
		powerValue=(aCell->position-position).squaredNorm()-aCell->r2;
	}
	inline void setTo(const PDCoord pos)
	{
		position=pos;
		for(int g=dimension;g>=0;g--)
		{
			endPoints[g]=NULL;
			generators[g]=NULL;
		}
	}

	inline	PDFloat powerdiff3D(const_cellPtr const& aCell,const_cellPtr const& bCell)const
	{
		//this has best accuracy when vertex is far away and atoms are close. something similar for a far atom is :
		// -bCall->r2+aCell->r2-(closeCell->position-position).squaredNorm()+2.0*((closeCell->position-position).dot(position-farCell->position)
		return -bCell->r2+aCell->r2-(aCell->position-bCell->position).squaredNorm()+2.0*((aCell->position-bCell->position).dot(position-bCell->position));
	}

	template<class PDCalc>
	inline void endPointsAndPositionOverwrite(const vertexPtr& endPoint,const PDCalc& pos)
	{
		endPoints[0]=endPoint;
		rrv=0;
		invalid=0;
		position=pos;
	}

	void refreshAfterRealloc(const vertex*const& copy)
	{
		// Rebase endpoint pointers after vertex vector reallocation.
		for(int g=dimension;g>=0;g--)
			if(this->generators[g]!=NULL&&(!this->generators[g]->myVertices.empty())&&this->generators[g]->myVertices.front()==copy)
				this->generators[g]->myVertices.front()=this;
		for(typename std::array<vertexPtr,dimension+1>::iterator it=endPoints.begin();it!=endPoints.end();++it)
		{
				if((*it)!=NULL)
					*it=this+(*it-copy);
		}
	}

	inline void moveAddressNetworkUpdateOnly(const vertexPtr& whereTo)
	{
		// Move one vertex object to a new address and patch neighboring endpoint links.
		if(this->endPoints[dimension]!=NULL)
			(this->endPoints[dimension]->fastWhichis(this))=whereTo;
		for(typename std::array<vertexPtr,dimension+1>::iterator it=endPoints.begin();it!=endPoints.begin()+dimension;++it)
			((*it)->fastWhichis(this))=whereTo;

		*whereTo=*this;
	}
	inline vertexPtr& fastWhichis (const const_vertexPtr& comp)
	{
		for(typename std::array<vertexPtr,dimension+1>::iterator it=endPoints.begin()+dimension;it!=endPoints.begin();--it)
			if(*it==comp)return *it;
		return endPoints[0];
	}
	inline vertexPtr& persistingWhichis3D (const const_vertexPtr& newOne)
	{
		if(generators[2]==newOne->generators[2])
			if(generators[1]==newOne->generators[1])
				return endPoints[0];
			else
				return endPoints[1];
		else
			if(generators[2]!=newOne->generators[3])
				return endPoints[2];
			else
				return endPoints[3];
	}




		bool cornerToReplacedAndGo(PowerDiagram<PDFloat,PDCoord,dimension>& owner)
		{
			// Mark a corner vertex as replaced and continue replacement flood-fill through neighbors.
			owner.Replaced.push_back(this);
				for(typename std::array<cellPtr,dimension+1>::const_iterator it=this->generators.begin();it!=this->generators.end();++it)
					if((*it)->visitedAs==0)
						owner.AddToInvolved(*(*it));
				owner.push_cell_my_vertex(*owner.Involved.front(), this);//although replaced it will be part of the new cell!its a corner!


			for(typename std::array<vertexPtr,dimension+1>::const_iterator it=this->endPoints.begin()+dimension;it!=this->endPoints.begin();--it)
				if((*it)->rrv==0)
				{
					if(!(*it)->replaceCheck(owner))
						return false;
				}else{}

			return true;
		}
		bool finiteToReplacedAndGo(PowerDiagram<PDFloat,PDCoord,dimension>& owner)
		{
			// Mark a finite vertex as replaced and continue replacement flood-fill through neighbors.
			owner.Replaced.push_back(this);
			for(typename std::array<cellPtr,dimension+1>::const_iterator it=this->generators.begin();it!=this->generators.end();++it)
				if((*it)->visitedAs==0)
					owner.AddToInvolved(*(*it));

			for(typename std::array<vertexPtr,dimension+1>::const_iterator it=this->endPoints.begin();it!=this->endPoints.end();++it)
				if((*it)->rrv==0)
					if(!(*it)->replaceCheck(owner))
						return false;

			return true;
		}
		inline bool replaceCheck( PowerDiagram<PDFloat,PDCoord,dimension>& owner)
		{
			// Dispatch replacement test according to corner/non-corner vertex type.
			if(this->isCorner())
				return this->cornerReplaceCheck(owner);
			return this->finiteReplaceCheck(owner);
		}

		bool finiteReplaceCheck(PowerDiagram<PDFloat,PDCoord,dimension>& owner)
		{
			const typename PowerDiagram<PDFloat,PDCoord,dimension>::ReplaceState state=
				owner.finiteReplaced(*this,owner.Involved.front());
			if(state==PowerDiagram<PDFloat,PDCoord,dimension>::ReplaceState::ambiguous)
				return false;
			if(state==PowerDiagram<PDFloat,PDCoord,dimension>::ReplaceState::replaced)
				return this->finiteToReplacedAndGo(owner);
			return true;
		}
		bool cornerReplaceCheck(PowerDiagram<PDFloat,PDCoord,dimension>& owner)
		{
			const typename PowerDiagram<PDFloat,PDCoord,dimension>::ReplaceState state=
				owner.finiteReplaced(*this,owner.Involved.front());
			if(state==PowerDiagram<PDFloat,PDCoord,dimension>::ReplaceState::ambiguous)
				return false;
			if(state==PowerDiagram<PDFloat,PDCoord,dimension>::ReplaceState::replaced)
				return this->cornerToReplacedAndGo(owner);
			return true;
		}
		template <const int cornerInfo>
		bool buildIn(PowerDiagram<PDFloat,PDCoord,dimension>*const& pd)const
		{
			// Build all required new vertices on outgoing edges whose far endpoint persists.
			for(int g=dimension;g>=cornerInfo;g--)
			{
				if(this->endPoints[g]->rrv<=0)
				{
					if(!pd->tryToBuildVertexOnEdge(this,g))
						return false;
				}
			}
			return true;
		}


	inline void registerForConnection3D(PowerDiagram<PDFloat,PDCoord,dimension>*const& owner)
	{
		// Register candidate edge endpoints so matching generator pairs can be connected.
		owner->planes[generators[2]->visitedAs*owner->Involved.size()+generators[1]->visitedAs].storeOrConnect(this,endPoints[3]);
		owner->planes[generators[3]->visitedAs*owner->Involved.size()+generators[1]->visitedAs].storeOrConnect(this,endPoints[2]);
		owner->planes[generators[3]->visitedAs*owner->Involved.size()+generators[2]->visitedAs].storeOrConnect(this,endPoints[1]);
	}


};
struct EdgeEnds
{
	vertexPtr a;
	vertexPtr* b;
	inline void storeOrConnect(const vertexPtr& pvertex, vertexPtr& itsEndPointStorage)
	{
		// First endpoint stores itself; second endpoint closes the pair and connects both vertices.
		if(this->a==NULL)
		{
			this->a=pvertex;								//we store ourself
			this->b=&itsEndPointStorage; //and where the other should write itself into
		}
		else
		{
			itsEndPointStorage=this->a;//we connect ourself to the other
			*(this->b)=pvertex;								//and the other to us
			this->a=NULL;
		}
	}
	inline void connect(const vertexPtr pvertex, vertexPtr& itsEndPointStorage)
	{
		// Connect against a previously stored endpoint if one exists for this key.
		if(this->a!=NULL)
		{
			itsEndPointStorage=this->a;//we connect ourself to the other
			*(this->b)=pvertex;								//and the other to us
			this->a=NULL;
		}
	}
};

};
//static initializations
//template <class PDFloat, class PDCoord>
//std::vector<cell<PDFloat,PDCoord> > PowerDiagram<PDFloat,PDCoord>::sideGenerators; //outside 



} //namespace POWER_DIAGRAM
#endif /* POWER_DIAGRAM_H_ */
