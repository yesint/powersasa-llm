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


#ifndef POWERSASA_H_
#define POWERSASA_H_

#include <cstdlib>
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include "power_diagram.h"
#include <memory>
#include <vector>

//=============================================================

namespace POWERSASA
{


class PowerSasaException : public std::exception {};

/*!
 * \class PowerSasa
 * \brief Sasa calculation class
 *
 * This class calculates the Sasa and stores it by using a power diagram.
 *
 * \endcode
 */

template <class PDFloat, class PDCoord>      // floating type and vector type
class PowerSasa 
{
public:
	using PowerDiagramT = POWER_DIAGRAM::PowerDiagram<PDFloat, PDCoord, 3>;
	static constexpr unsigned int kMaxNb = 20;
	static constexpr unsigned int kMaxVx = 12;
	static constexpr unsigned int kMaxPnt = 4;
	static constexpr int kMaxCount = 100;

	// Relative tolerance used when deciding whether two radii/powers should be treated as equal.
	inline static PDFloat DRAD2() { return static_cast<PDFloat>(1000.0) * std::numeric_limits<PDFloat>::epsilon() ; }
	// Angular ordering tolerance used when ranking contour vertices on a circle.
	inline static PDFloat DANG() { return static_cast<PDFloat>(1000.0) * std::numeric_limits<PDFloat>::epsilon() ; }
	// Cosine threshold around |1| where cross-product based angle reconstruction is more stable.
	inline static PDFloat near_one_cosine_threshold() { return static_cast<PDFloat>(0.999); }
	// Minimum axis component considered numerically reliable in signed-angle orientation tests.
	inline static PDFloat axis_component_threshold() { return static_cast<PDFloat>(0.001); }
	// Clamp cosine-like values before acos() to avoid NaNs from tiny floating-point overshoots.
	inline static PDFloat clamp_unit_interval(const PDFloat value)
	{
		if (value < static_cast<PDFloat>(-1.0)) return static_cast<PDFloat>(-1.0);
		if (value > static_cast<PDFloat>(1.0)) return static_cast<PDFloat>(1.0);
		return value;
	}

	const std::vector<PDFloat>& getSasa() const {return Sasa;}
	const std::vector<PDFloat>& getVol() const {return Vol;}
	const std::vector<PDCoord>& getDVol() const {return DVol;}
	const std::vector<PDCoord>& getDSasa() const {return DSasa;}

	// Compute all geometric contributions for one atom: SASA, derivatives, volume and volume derivatives.
	void calc_sasa_single(const unsigned int iatom);
	// Compute contributions for every atom using the current power diagram topology.
	void calc_sasa_all();

	template<class Coordcontainer, class Floatcontainer>
	// Rebuild the underlying power diagram after coordinate/radius updates and resize output buffers if atom count grows.
	void update_coords(Coordcontainer const& coords, Floatcontainer const& radii)
	{
		const std::size_t n_old = power_diagram->get_points().size();
		power_diagram->recalculate(coords.begin(),radii.begin(),coords.size());
		if (n_old < power_diagram->get_points().size()) Resize_NA();
	}

	template <class Pos_iterator, class Strength_iterator>
	// Incrementally add one or more points to the power diagram and refresh atom-sized output storage.
	void add_more(const Pos_iterator pos_it, const Strength_iterator strength_it, const unsigned int newSize) {
		power_diagram->addMore(pos_it, strength_it, newSize);
		Resize_NA();
	}

	// Convenience single-point incremental insertion.
	void add_more(const PDCoord& position, const PDFloat& radius,const int near) {
		power_diagram->addMore(position,radius,near);
		Resize_NA();
	}
	// Undo the latest incremental insertion sequence performed through add_more().
	void revert() {power_diagram->revert();}
	// Access to the underlying power-diagram data structure (cells, vertices, zero points).
	POWER_DIAGRAM::PowerDiagram<PDFloat,PDCoord,3 > & get_power_diagram()  { return *power_diagram; }

	// Number of neighboring cells for atom iatom in the current power diagram.
	unsigned int NumOfNeighbours(unsigned int iatom) const
	{
		return static_cast<unsigned int>(power_diagram->get_points()[iatom].neighboursIds.size());
	}
	// Convert local neighbour index to global atom index.
	unsigned int AtomNo(unsigned int i_atom, unsigned int i_neighbour) const
	{
		const auto& points = power_diagram->get_points();
		const auto& atom = points[i_atom];
		return static_cast<unsigned int>(atom.neighboursIds[i_neighbour]);
	}

	template<class Coordcontainer, class Floatcontainer, class Intcontainer>
	// Construct from coordinates, radii and explicit bond-to hints used by incremental insertion heuristics.
	PowerSasa(Coordcontainer const& coords, Floatcontainer const& radii, Intcontainer const& bond_to,
		const bool with_Sasa, const bool with_dSasa, const bool with_Vol, const bool with_dVol) :
		power_diagram(std::make_unique<PowerDiagramT>(
			PowerDiagramT::create(coords.size(), coords.begin(), radii.begin(), bond_to.begin())
			.with_radiiGiven(1).with_calculate(1).with_cells(1).with_zeroPoints(1).with_Warnings(0).withoutCheck(1))),
		withSasa(with_Sasa), withDSasa(with_dSasa), withVol(with_Vol), withDVol(with_dVol)
	{
		Init();
	}
	template<class Coordcontainer, class Floatcontainer>
	// Construct from coordinates/radii only; bond-to hints are generated as a simple chain.
	PowerSasa(Coordcontainer const& coords, Floatcontainer const& radii,
				const bool with_Sasa=1, const bool with_dSasa=0, const bool with_Vol=0, const bool with_dVol=0) :
			power_diagram(nullptr), withSasa(with_Sasa), withDSasa(with_dSasa), withVol(with_Vol), withDVol(with_dVol)
	{
		std::vector<int> bond_to;
		bond_to.reserve(coords.size());
		bond_to.push_back(0);
		for (unsigned int i = 1; i < coords.size(); ++i)
		{
			bond_to.push_back(i-1);
		}
		power_diagram = std::make_unique<PowerDiagramT>(
			PowerDiagramT::create(coords.size(), coords.begin(), radii.begin(), bond_to.begin())
			.with_radiiGiven(1).with_calculate(1).with_cells(1).with_zeroPoints(1).with_Warnings(0).without_Check(1));
		Init();
	}

	virtual ~PowerSasa() = default;
private:
	// Allocate all per-atom/per-neighbour scratch buffers to default capacities.
	inline void Init()
	{
		Resize_NA();
		Resize_NB(kMaxNb);
		Resize_VX(kMaxVx);
		Resize_PNT(kMaxPnt);
	}
	// Ensure neighbour-indexed scratch arrays are large enough for current atom degree.
	inline void Resize_NB(unsigned int nnb)
	{
		np.resize(nnb);
		nt.resize(nnb);
		e.resize(nnb);
		sintheta.resize(nnb);
		costheta.resize(nnb);
		nb_RAD2.resize(nnb);
		nb_dist.resize(nnb);
		volnb.resize(nnb);
		knot.resize(nnb);
		fknot.resize(nnb);
		std::size_t npnt = kMaxPnt;
		if (!p.empty()) npnt = p[0].size();
		const std::size_t nnb_old = p.size();
		next.resize(nnb);
		p.resize(nnb);
		ang.resize(nnb);
		for (std::size_t i = nnb_old; i < nnb; ++i)
		{
			next[i].resize(npnt);
			p[i].resize(npnt);
			ang[i].resize(npnt);
		}
		if (withDSasa)
		{
			for (std::size_t i = 0; i < DSasa_parts.size(); ++i) DSasa_parts[i].resize(nnb);
		}
	}
	// Ensure global surface-vertex storage arrays are large enough.
        inline void Resize_VX(unsigned int nvx)
	{
		off.resize(nvx);
		vx.resize(nvx);
		br_c.resize(nvx);
		br_p.resize(nvx);
		for (unsigned int i = 0; i < nvx; ++i)
		{
			br_c[i].resize(2);
			br_p[i].resize(2);
		}
	}
	// Ensure per-neighbour contour-point arrays are large enough.
	inline void Resize_PNT(unsigned int npnt)
	{
		for (std::size_t i = 0; i < p.size(); ++i)
		{
			next[i].resize(npnt);
			p[i].resize(npnt);
			ang[i].resize(npnt);
		}
	        rang.resize(npnt);
		pos.resize(npnt);
	}
	// Resize per-atom outputs and refresh global tolerance derived from current maximum radius^2.
	inline void Resize_NA()
	{
		const std::size_t n = power_diagram->get_points().size();
		if (withSasa) Sasa.resize(n,0);
		if(withDSasa)
		{
			std::size_t nnb = kMaxNb;
			if (!DSasa_parts.empty() && !DSasa_parts.front().empty()) nnb = DSasa_parts.front().size();
			const std::size_t n_old = DSasa.size();
			DSasa.resize(n);
			DSasa_parts.resize(n);
			for (std::size_t i = n_old; i < n; ++i) DSasa_parts[i].resize(nnb);
			
		}
		if (withVol) Vol.resize(n,0);
		if (withDVol) DVol.resize(n,PDCoord(0,0,0));

		PDFloat maxr2 = 0.0;
		for (std::size_t i = 0; i < n; ++i)
		{
			if (power_diagram->get_points()[i].r2 > maxr2)
			  maxr2 = power_diagram->get_points()[i].r2;
		}
		tol_pow = maxr2 * DRAD2();
	}

	// Signed twist angle from a to b around axis c in [0, 2*pi).
 	inline PDFloat Ang_About(PDCoord const& a, PDCoord const& b, PDCoord const& c);
	// Compute per-vertex polar angles on a neighbour intersection circle.
	inline void Get_Ang(const int & np, const std::vector<int> & p, const PDCoord & e,
			const PDFloat & sintheta, const PDFloat & costheta, std::vector<PDFloat> & ang);
	// Build cyclic next-vertex links by ordering circle angles with degeneracy handling.
        inline void Get_Next(int n, std::vector<PDFloat> & ang, std::vector<int> & next,
			const std::vector<int> & p, const PDCoord & e);
	std::unique_ptr<PowerDiagramT> power_diagram;
	
	const bool withSasa;
	const bool withDSasa;
	const bool withVol;
	const bool withDVol;
	std::vector<PDFloat>                  Sasa;
	std::vector< std::vector<PDCoord> >   DSasa_parts;
	std::vector<PDCoord>		      DSasa;
	std::vector<PDFloat>                  Vol;
	std::vector<PDCoord>                  DVol;
	
	PDFloat                               tol_pow;

	// ----- for calc_sasa_single --------------------

	std::vector<int> np;            // number of points (registered vertices) of i-th atom
	std::vector<int> nt;            // counts situations that exclude "single cirle"
	std::vector<PDCoord> e;         // direction to neighbour
        std::vector<PDFloat> sintheta;
        std::vector<PDFloat> costheta;
	std::vector<PDFloat> nb_RAD2;
	std::vector<PDFloat> nb_dist;

	std::vector<int> off;           // off[i] != 0 if vertix i is already taken into account
	std::vector<PDCoord> vx;        // all surface vertices
	std::vector< std::vector<int> > br_c;       // bridge between circles (intersections with neighbors)
	std::vector< std::vector<int> > br_p;       // bridge between point numbers
  
	std::vector< std::vector<PDFloat> > ang; // angles
	std::vector< std::vector<int> >next;     // next[][i] is the number of angle that follows the i-th one
	std::vector< std::vector<int> > p;       // number of surface vertices

	std::vector<PDFloat> volnb;
	std::vector<PDCoord> knot;
	std::vector<bool>    fknot;
        //----- for Get_Next ----------------------------

	std::vector<int> rang;
	std::vector<int> pos;
	// special thing
public:
	template <const int n>
	// Evaluate SASA/volume response of each atom to n probe-radius perturbations via temporary add/revert.
	void calc_sasa_all(const PDFloat steps[n],PDFloat resultS[][n+1],PDFloat resultV[][n+1])
	{
		for (unsigned int i = 0; i < power_diagram->getPoints().size(); ++i)
		{
			calc_sasa_single(i);
			resultS[i][0]=Sasa[i];
			resultV[i][0]=Vol[i];
		}
		const int& size=power_diagram->getPoints().size();
		for(int a=0;a<size;a++)
		{
			for(int s=0;s<n;s++)
			{
				add_more(power_diagram->getPoints()[a].position+power_diagram->center,power_diagram->getPoints()[a].r+steps[s],a);
				calc_sasa_single(size);
				resultS[a][s]=Sasa[size];
				resultV[a][s]=Vol[size];
				revert();
			}
		}
	};
};

//Implementation:

//-----------------------------------------------------------------------------
// Retuns twist angle (in interval [0, 2*pi()) ) between Vectors a and b about c.
// a and b should be of unit length and perpendicular to c.

template <class PDFloat, class PDCoord> PDFloat PowerSasa<PDFloat, PDCoord>::
Ang_About(PDCoord const& a, PDCoord const& b, PDCoord const& c)
{
	PDFloat ang, co, vp;
	PDCoord  v;

	co = a.dot(b);
	if (co <= -near_one_cosine_threshold())
	{
		v = a.cross(b);
		ang = std::numbers::pi_v<PDFloat> - std::asin( v.norm() );
	}
	else if (co >= near_one_cosine_threshold())
	{
		v = a.cross(b);
		ang = std::asin( v.norm() );
	}
	else ang = std::acos(co);

	if (std::fabs(c[0]) > axis_component_threshold())
	{
		vp = a[1]*b[2] - a[2]*b[1];
		if ( (vp < 0.0) != (c[0] < 0.0) ) ang = - ang;
	}
	else if (std::fabs(c[1]) > axis_component_threshold())
	{
		vp = a[2]*b[0] - a[0]*b[2];
		if ( (vp < 0.0) != (c[1] < 0.0) ) ang = - ang;
	}
	else if (std::fabs(c[2]) > axis_component_threshold())
	{
		vp = a[0]*b[1] - a[1]*b[0];
		if ( (vp < 0.0) != (c[2] < 0.0) ) ang = - ang;
	}
	else
	{
		std::cerr << "PowerSasa: Axis too short" << std::endl;
		throw PowerSasaException();
	}
	if (ang < 0) ang += static_cast<PDFloat>(2.0) * std::numbers::pi_v<PDFloat>;
	return ang;
}
//-----------------------------------------------------------------------------
//Finds phi-angles of s-vertices

template <class PDFloat, class PDCoord> void PowerSasa<PDFloat, PDCoord>::
Get_Ang(const int & np,                              // total number of s-vertices for given circle
        const std::vector<int> & p,                  // s-vertice number for given circle
	const PDCoord & e,                           // direction to neighbour
        const PDFloat & sintheta, const PDFloat & costheta,  // theta angle
	std::vector<PDFloat> & ang)                          // output
	
{
//	if (np == 0) return;//is tested outside
	static PDCoord pu0, pu;
	ang[0] = 0.0;
	pu0 = (vx[p[0]] - costheta * e) / sintheta;
	for (int j = 1; j < np; j++)
	{
		pu = (vx[p[j]] - costheta * e) / sintheta;
		ang[j] =  Ang_About (pu0, pu, e);
	}
}


//-----------------------------------------------------------------------------
// given array ang of length n founds number next[i] of ang element that is next in magnitude to element i
// can be optimized (???)

template <class PDFloat, class PDCoord> void PowerSasa<PDFloat, PDCoord>::
Get_Next(int n, std::vector<PDFloat> & ang, std::vector<int> & next,
	 const std::vector<int> & p, const PDCoord & e)

{
	int j, k, m;
	//static int rang[MAX_PNT], pos[MAX_PNT];
	if (n == 0) return;

	for (j = 1; j < n; ++j)
	{
		if (ang[j] <= static_cast<PDFloat>(2.0) * std::numbers::pi_v<PDFloat> - DANG()) continue;
		PDFloat dp = vx[p[j]].cross(vx[p[0]]).dot(e);
		if (dp < 0.0) ang[j] = 0.0;
		else if (dp == 0.0)
		{
			std::cout << "PowerSasa: Precision insufficient to resolve angles" << std::endl;
			throw PowerSasaException();
		}
	}

	rang[0] = 0;
	rang[1] = 1;

	for (j = 2; j < n; j++)
	{
		m = j;
		for (k = 1; k < j; k++)
		{
			if (ang[k] > ang[j] + DANG())
			{
				(rang[k])++;
				m--;
			}
			else if (ang[k] > ang[j] - DANG())
			{

			        PDFloat dp = vx[p[j]].cross(vx[p[k]]).dot(e);
				if (dp > 0.0)
				{
					(rang[k])++;
					m--;
			        }
			        else if (dp == 0.0)
				{
					std::cout << "PowerSasa: Precision insufficient to resolve angles" << std::endl;
					throw PowerSasaException();
				}
			}
		}
		rang[j] = m;
	}

	for (j = 0; j < n; j++) pos[rang[j]] = j;

	for (j = 0; j < n; j++)
	{
		if (rang[j] == n-1) next[j] = pos[0];
		else next[j] = pos[ rang[j]+1 ];
	}
}
//===================================================================================

template <class PDFloat, class PDCoord> void PowerSasa<PDFloat, PDCoord>::
calc_sasa_single(const unsigned int iatom)
{

// ----- some initialisations  --------------------

	int i, j, kn;
	bool ok;
	using PowerDiagram3D = typename POWER_DIAGRAM::PowerDiagram<PDFloat, PDCoord, 3>;
	using CellId = typename PowerDiagram3D::CellId;
	using VertexId = typename PowerDiagram3D::VertexId;
	using GeneratorKind = typename PowerDiagram3D::GeneratorKind;

	std::vector<typename PowerDiagram3D::cell> const& atoms = power_diagram->get_points();
	const typename PowerDiagram3D::cell& atom = atoms[iatom];
	const int nnb = static_cast<int>(atom.neighboursIds.size());               // number of neighbors
	std::vector<CellId> neighbour_ids(static_cast<std::size_t>(nnb), PowerDiagram3D::kInvalidId);
	for (int nb_idx = 0; nb_idx < nnb; ++nb_idx)
	{
		const std::size_t idx = static_cast<std::size_t>(nb_idx);
		if (idx < atom.neighboursIds.size()) neighbour_ids[idx] = atom.neighboursIds[idx];
		if (neighbour_ids[idx] == PowerDiagram3D::kInvalidId)
		{
			std::cerr << "PowerSasa: Invalid neighbour link" << std::endl;
			throw PowerSasaException();
		}
	}
	auto neighbour = [&](const int idx) -> const typename PowerDiagram3D::cell&
	{
		return power_diagram->get_cell(neighbour_ids[static_cast<std::size_t>(idx)]);
	};
	auto neighbour_mut = [&](const int idx) -> typename PowerDiagram3D::cell&
	{
		return power_diagram->get_cell(neighbour_ids[static_cast<std::size_t>(idx)]);
	};
	PDCoord const &pos = atom.position;       // my coordinates

	if (static_cast<std::size_t>(nnb) >= np.size())
	{
		Resize_NB(nnb);
		//std::cerr << "Power Sasa: number of neighbors exeeds MAX_NB=" << MAX_NB << std::endl;
		//throw PowerSasaException();
	}

	if (withSasa) Sasa[iatom] = 0.0;
	if (withDSasa){DSasa[iatom].setZero(); for (i = 0; i < nnb; ++i) DSasa_parts[iatom][i].setZero();}
	if (withVol) Vol[iatom] = 0.0;
	if (withDVol) DVol[iatom].setZero();
	const bool do_sasa = withSasa || withVol;

	const PDFloat RAD  = atom.r;
	const PDFloat RAD2 = atom.r2;

	if (nnb == 0)
	{
		const auto& first_vertex = power_diagram->get_vertices()[0];
		bool is_owner = false;
		const auto& gref = first_vertex.generatorRefs[0];
		if (gref.is_valid() && gref.kind == GeneratorKind::point && gref.index == iatom)
		{
			is_owner = true;
		}
		if (is_owner)
		{
			if (withSasa)Sasa[iatom] = 4*3.1415926535897932384626433832795*RAD2;
			if(withVol)Vol[iatom] = 4*3.1415926535897932384626433832795*0.33333333333333333333333333333333*RAD*RAD2;
		}
		return;
	}

// check that the Power values are not close to zero

	bool covered = true;
	ok = true;
	
	const auto& pd_vertices = power_diagram->get_vertices();
	const std::size_t my_vertex_count = atom.myVerticesIds.size();
	for (std::size_t n = 0; n < my_vertex_count; ++n)
	{
		const std::size_t vid = atom.myVerticesIds[n];
		if (vid >= pd_vertices.size()) continue;
		const auto& atom_vertex = pd_vertices[vid];
		if (std::fabs(atom_vertex.powerValue) < tol_pow)
		{
			std::cerr << "PowerSasa: Vertex power value is too close to zero" << std::endl;
			throw PowerSasaException();
		}
		if (atom_vertex.powerValue > 0.0) covered = false;
	}

	if (covered && !withVol) return;

//----- get angles and other properties of neighbors -----

	int n_apart = 0;                // number of non-contributing neighbours
	PDFloat dist, nb_RAD;           // distance to neighbor, RADius of neighbor
	PDCoord  rel_pos;		// vector to neighbor
	for (i = 0; i < nnb; ++i)       // over neighbors
	{	  
		if (withVol)
		{
			volnb[i] = 0.0;
			fknot[i] = 0;
		}
		neighbour_mut(i).visitedAs = i;
		rel_pos = neighbour(i).position - pos;       // vector to neighbor
		dist = rel_pos.norm();                  // distance to neighbor
		if (dist == 0)
		{
			std::cerr << "PowerSasa: Invalid distance to neighbour" << std::endl;
			throw PowerSasaException();
		}
		costheta[i] = 1.0;                      // will be overwritten
		nb_RAD = neighbour(i).r;                     // neighbor RADius
		nb_RAD2[i] = neighbour(i).r2;                // neighbor RADius^2
		np[i] = 0;				// initialize number of points
		nt[i] = 0;                              // ?????

		if (dist <= nb_RAD - RAD) // totally covered
		{
			return;          
		}

		if (dist >= RAD + nb_RAD || dist <= RAD - nb_RAD)
		{
			++n_apart;
			np[i] = -1;
			continue;         // neighbour i does not contribute
		}

		costheta[i] = 0.5 * ( dist + (RAD2 - nb_RAD2[i]) / dist ) / RAD;

		if (costheta[i] <= -1.0)  // totally covered (impossible but still...)
		{
			return;          
		}
		if (costheta[i] >= 1.0)   // impossible but still...
		{
			++n_apart;
			np[i] = -1;
			continue;         // neighbour i does not contribute
		}

		sintheta[i] = std::sqrt(1.0 - costheta[i] * costheta[i]);
		e[i] = rel_pos / dist;
		nb_dist[i] = dist;
	}

	if (n_apart == nnb) // no contributing neighbours
	{
		if (withSasa) Sasa[iatom] = static_cast<PDFloat>(4.0) * std::numbers::pi_v<PDFloat> * RAD2;
		if (withVol)  Vol[iatom] = (static_cast<PDFloat>(4.0) / static_cast<PDFloat>(3.0)) * std::numbers::pi_v<PDFloat> * RAD * RAD2;
		return; 
	}
	

//------ register surface vertices ------------------------

        int partner[2], ptn, ptn0, ptn1;
	PDCoord zp_pos;

	int nvx = 0;
	for (unsigned int k = 0; k < atom.myZeroPoints.size(); ++k)
	{
		const std::size_t zp_id = static_cast<std::size_t>(atom.myZeroPoints[k]);
		if (zp_id >= power_diagram->get_zeroPoints().size()) continue;
		const auto& zp = power_diagram->get_zeroPoints()[zp_id];
		if(power_diagram->zeroPointValid(zp))
	{
		zp_pos = power_diagram->zeroPointPos(zp);
		ptn = 0;
		for (int kg = 0; kg < 3; ++kg)
		{
			const auto& zp_generator_ref = zp.generatorRefs[kg];
			if (!zp_generator_ref.is_valid()) continue;
			if (zp_generator_ref.kind != GeneratorKind::point || zp_generator_ref.index != iatom)
			{
				partner[ptn] = power_diagram->get_generator(zp_generator_ref).visitedAs;
				++ptn;
			}
		}
		ptn0 = partner[0];
		ptn1 = partner[1];

		if (zp.pos < 0.0 || 1.0 < zp.pos)
		{
			++nt[ptn0];
			++nt[ptn1];
			continue;
		}

		if ((np[ptn0] < 0) || (np[ptn1] < 0))
		{
			std::cerr << "PowerSasa: Invalid surface vertex" << std::endl;
			throw PowerSasaException();
		}
		if ((np[ptn0] >= static_cast<int>(rang.size())) || (np[ptn1] >= static_cast<int>(rang.size())))
		{
			Resize_PNT( (np[ptn0] > np[ptn1]) ? np[ptn0]+1 : np[ptn1]+1 );
			//std::cerr << "PowerSasa: Number of surface vertices for single neighbor exeeds MAX_PNT=" << MAX_PNT << std::endl;
			//throw PowerSasaException();
		}
		if (nvx >= static_cast<int>(vx.size()))
		{
		        Resize_VX(nvx+1);
			//std::cerr << "PowerSasa: number of surface vertices exeeds MAX_VX=" << MAX_VX << std::endl;
			//throw PowerSasaException();
		}

		vx[nvx] = (zp_pos - pos) / RAD;               // faster
		//vx[nvx] = (zp_pos - pos).normalized();      // more accurate?
		p[ptn0][np[ptn0]] = p[ptn1][np[ptn1]] = nvx;  // s-vertex number

		br_c[nvx][0] = ptn0;
		br_c[nvx][1] = ptn1;
		br_p[nvx][0] = np[ptn0];
		br_p[nvx][1] = np[ptn1];
		++nvx;
		++np[ptn0];
		++np[ptn1];

			const VertexId node1_id = zp.fromId;
			if (node1_id == PowerDiagram3D::kInvalidId || node1_id >= pd_vertices.size())
			{
				std::cerr << "PowerSasa: Invalid zeroPoint source vertex" << std::endl;
				throw PowerSasaException();
			}
			const auto& node1 = pd_vertices[node1_id];
			VertexId node2_id = node1.endPointIds[zp.branch];
			if (node2_id == PowerDiagram3D::kInvalidId)
			{
				const auto* endpoint_ptr = node1.endPoints[zp.branch];
				if (endpoint_ptr != nullptr) node2_id = power_diagram->get_vertex_id(*endpoint_ptr);
			}
			if (node2_id == PowerDiagram3D::kInvalidId || node2_id >= pd_vertices.size())
			{
				std::cerr << "PowerSasa: Invalid zeroPoint endpoint" << std::endl;
				throw PowerSasaException();
			}
			const auto& node2 = pd_vertices[node2_id];
			if (node1.powerValue < 0.0 && node2.powerValue > 0.0)
			{
			if (fknot[ptn0] == 0)
			{
				fknot[ptn0] = 1;
				knot[ptn0] = node1.position;
			}
			else volnb[ptn0] +=
			  std::abs((node1.position - knot[ptn0]).cross(zp_pos - knot[ptn0]).dot(e[ptn0]));
			if (fknot[ptn1] == 0)
			{
				fknot[ptn1] = 1;
				knot[ptn1] = node1.position;
			}
			else volnb[ptn1] +=
			  std::abs((node1.position - knot[ptn1]).cross(zp_pos - knot[ptn1]).dot(e[ptn1]));
		}
		else if (node1.powerValue > 0.0 && node2.powerValue < 0.0)
		{
			if (fknot[ptn0] == 0)
			{
				fknot[ptn0] = 1;
				knot[ptn0] = node2.position;
			}
			else volnb[ptn0] +=
			  std::abs((node2.position - knot[ptn0]).cross(zp_pos - knot[ptn0]).dot(e[ptn0]));
			if (fknot[ptn1] == 0)
			{
				fknot[ptn1] = 1;
				knot[ptn1] = node2.position;
			}
			else volnb[ptn1] +=
			  std::abs((node2.position - knot[ptn1]).cross(zp_pos - knot[ptn1]).dot(e[ptn1]));
		}
		else if (node1.powerValue > 0.0 && node2.powerValue > 0.0)
		{
			PDFloat dpos = node1.powerValue*(1.0 - zp.pos) /
			  (node2.powerValue*zp.pos + node1.powerValue*(1.0 - zp.pos)) - zp.pos;
			if (fknot[ptn0] == 0)
			{
				fknot[ptn0] = 1;
				knot[ptn0] = zp_pos;
			}
			else volnb[ptn0] += 0.5 *
			  std::abs(dpos*(zp_pos - knot[ptn0]).cross(node2.position - node1.position).dot(e[ptn0]));
			if (fknot[ptn1] == 0)
			{
				fknot[ptn1] = 1;
				knot[ptn1] = zp_pos;
			}
			else volnb[ptn1] += 0.5 *
			  std::abs(dpos*(zp_pos - knot[ptn1]).cross(node2.position - node1.position).dot(e[ptn1]));
		}
		else
		{
			std::cerr << "PowerSasa: Impossible zeroPoint" << std::endl;
			throw PowerSasaException();
		}
	}
	}

	const std::size_t node_count = atom.myVerticesIds.size();
	const auto generators_match = [](const typename PowerDiagram3D::GeneratorRef& a, const typename PowerDiagram3D::GeneratorRef& b) -> bool
	{
		return a.is_valid() && b.is_valid() && a.kind == b.kind && a.index == b.index;
	};

	for (j = 0; j < static_cast<int>(node_count); ++j)
	{
		const VertexId node1_id = atom.myVerticesIds[j];
		if (node1_id == PowerDiagram3D::kInvalidId || node1_id >= pd_vertices.size()) continue;
		const auto& node1 = pd_vertices[node1_id];
		for (kn = 0; kn < 4; ++kn)
		{
					VertexId node2_id = node1.endPointIds[kn];
					if (node2_id == PowerDiagram3D::kInvalidId)
					{
						const auto* node2_ptr = node1.endPoints[kn];
						if (node2_ptr == nullptr) continue;
						node2_id = power_diagram->get_vertex_id(*node2_ptr);
					}
					if (node2_id == PowerDiagram3D::kInvalidId || node1_id == PowerDiagram3D::kInvalidId) continue;
					if (node2_id >= pd_vertices.size()) continue;
					const auto& node2 = pd_vertices[node2_id];
					if (node2_id > node1_id) continue;
					if (node1.powerValue > 0.0 || node2.powerValue > 0.0) continue;
				bool node2_contains_atom = false;
				for (int kg = 0; kg < 4; ++kg)
				{
					const auto& g2ref = node2.generatorRefs[kg];
					if (g2ref.is_valid())
					{
						if (g2ref.kind == GeneratorKind::point && g2ref.index == iatom)
						{
							node2_contains_atom = true;
							break;
						}
					}
				}
				if (!node2_contains_atom) continue;
				ptn = 0;
				for (int kg = 0; kg < 4; ++kg)
				{
					const auto& g1ref = node1.generatorRefs[kg];
					if (!g1ref.is_valid()) continue;
					if (g1ref.kind == GeneratorKind::point && g1ref.index == iatom) continue;
					bool shared_generator = false;
					for (int kh = 0; kh < 4; ++kh)
					{
						const auto& g2ref = node2.generatorRefs[kh];
						if (generators_match(g1ref, g2ref))
						{
							shared_generator = true;
							break;
						}
					}
					if (shared_generator)
					{
						partner[ptn] = power_diagram->get_generator(g1ref).visitedAs;
						++ptn;
					}
				}
			ptn0 = partner[0];
			ptn1 = partner[1];
			++nt[ptn0];
			++nt[ptn1];
			if (fknot[ptn0] == 0)
			{
				fknot[ptn0] = 1;
				knot[ptn0] = node1.position;
			}
			else volnb[ptn0] +=
			  std::abs((node1.position - knot[ptn0]).cross(node2.position - knot[ptn0]).dot(e[ptn0]));
			if (fknot[ptn1] == 0)
			{
				fknot[ptn1] = 1;
				knot[ptn1] = node1.position;
			}
			else volnb[ptn1] +=
			  std::abs((node1.position - knot[ptn1]).cross(node2.position - knot[ptn1]).dot(e[ptn1]));
		}
	}

// ----- get angles of s-vertices within each circle ---------

	for (i = 0; i < nnb; i++)
	{
		if (np[i] <= 0) continue;
		if (np[i] % 2 != 0)
		{
			std::cerr << "PowerSasa: odd number of crossing between circles of " << iatom
				<< " and " << AtomNo(iatom, static_cast<unsigned int>(i)) << std::endl;
			throw PowerSasaException();
		}
		Get_Ang(np[i], p[i], e[i], sintheta[i], costheta[i], ang[i]);
		Get_Next(np[i], ang[i], next[i], p[i], e[i]);
	}

// ---------- get sasa from contours ---------------------------

        int ic1, ic2, ic_0, ic_1, ip2, ip_next, ivx, count;
	PDFloat phi, co, dirdet, ds1, ds2, scone;
        int p_ini_idx, pt_idx, pt0_idx;
	PDCoord vv;
	PDFloat vol2 = 0.0, vol3 = 0.0;

	PDFloat sasa_ia = 0.0;
	for (int iv = 0; iv < nvx; iv++) off[iv] = 0;

	for (int iv = 0; iv < nvx; iv++)
	{
		if (off[iv]) continue;
		p_ini_idx = iv;
		ic_0 = br_c[iv][0];
		ic_1 = br_c[iv][1];

		// determine which directoin
		dirdet = (e[ic_1].cross(e[ic_0])).dot(vx[p_ini_idx]);
		if (dirdet == 0.0)
		{
		    std::cerr << "PowerSasa: dirdet == 0.0" << std::endl;
		    throw PowerSasaException();
		}
		//assert(dirdet != 0.0);
		if ( dirdet > 0.0 )
		{
			ic1 = ic_0;
//			ip1 = br_p[iv][0];
			ic2 = ic_1;
			ip2 = br_p[iv][1];
		}
		else
		{
			ic1 = ic_1;
//			ip1 = br_p[iv][1];
			ic2 = ic_0;
			ip2 = br_p[iv][0];
		}

		pt_idx = p_ini_idx;
		if (do_sasa) sasa_ia += static_cast<PDFloat>(2.0) * std::numbers::pi_v<PDFloat>;
		count = 0;

		do                                 // main loop
		{
			++count;
			if (count > kMaxCount)
			{
				std::cerr << "PowerSasa: Wrong contour" << std::endl;
				throw PowerSasaException();
			}
			//assert(count <= MAX_PNT);   // else wrong contour
				ip_next = next[ic2][ip2];
				phi = ang[ic2][ip_next] - ang[ic2][ip2];
				if (phi < 0.0) phi += static_cast<PDFloat>(2.0) * std::numbers::pi_v<PDFloat>;
				
				if (do_sasa)
				{
					co = (e[ic1].dot(e[ic2]) - costheta[ic1]*costheta[ic2])
						/ (sintheta[ic1] * sintheta[ic2]);

					co = clamp_unit_interval(co);
					sasa_ia += phi * costheta[ic2] - acos(co);
				}
			
			off[ p[ic2][ip2] ] = 1;
			ic1 = ic2;

			ivx = p[ic1][ip_next];
			pt0_idx = pt_idx;
			pt_idx = ivx;
			
			if (br_c[ivx][0] == ic1)
			{
				ic2 = br_c[ivx][1];
				ip2 = br_p[ivx][1];
			}
			else
			{
				ic2 = br_c[ivx][0];
				ip2 = br_p[ivx][0];
			}

				if (withDSasa) {
					ds1 = 0.5 * RAD * phi * 
						( 1.0 + (nb_RAD2[ic1] - RAD2)/(nb_dist[ic1] * nb_dist[ic1]) );
					ds2 = RAD2 / nb_dist[ic1];
					DSasa_parts[iatom][ic1]+= ds1 * e[ic1] - ds2 * (vx[pt_idx] - vx[pt0_idx]).cross(e[ic1]);
				}

				if (withVol || withDVol) {
					vv = (pos + RAD * vx[pt0_idx]).cross(pos + RAD * vx[pt_idx]);
					scone = sintheta[ic1]*sintheta[ic1]*(phi - sin(phi));
				}
				
				if (withVol)
				{
					vol2 += costheta[ic1]*scone;
				if (fknot[ic1] == 0)
				{
				  fknot[ic1] = 1;
				  knot[ic1] = pos+RAD*vx[pt0_idx];
				}
				else volnb[ic1] +=
				  std::abs((pos+RAD*vx[pt0_idx] - knot[ic1]).cross(pos+RAD*vx[pt_idx] - knot[ic1]).dot(e[ic1]));
			}
			
				if (withDVol)
				{
					DVol[iatom]-= 0.5 * (vv + (RAD2*scone) * e[ic1]);
				}


		} while (pt_idx != p_ini_idx);
		
		if (do_sasa && sasa_ia > static_cast<PDFloat>(4.0) * std::numbers::pi_v<PDFloat>) {
			sasa_ia -= static_cast<PDFloat>(4.0) * std::numbers::pi_v<PDFloat>;
		}
	}

// ---------- sasa from single circles ----------------------------

	PDFloat pw_i, pw_j;
	PDCoord  cc;

	for (i = 0; i < nnb; ++i)
	{
		if (np[i] != 0 || nt[i] != 0) continue;

			ok = true;
		cc = pos + (RAD*costheta[i])*e[i];
		pw_i = -sintheta[i]*sintheta[i]*RAD2;

		for (j = 0; j < nnb; ++j)
		{
			if (j == i) continue;
				pw_j = (neighbour(j).position - cc).squaredNorm() - nb_RAD2[j];
			        if (pw_j <= pw_i)
				{
					ok = false;
					break;
				}
		}
		if (ok)
		{
			if (do_sasa)
			{  
				sasa_ia += static_cast<PDFloat>(2.0) * std::numbers::pi_v<PDFloat> * (static_cast<PDFloat>(1.0) + costheta[i]);
				if (sasa_ia > static_cast<PDFloat>(4.0) * std::numbers::pi_v<PDFloat>) {
					sasa_ia -= static_cast<PDFloat>(4.0) * std::numbers::pi_v<PDFloat>;
				}
			}
				if (withDSasa)
				{
					DSasa_parts[iatom][i] = (RAD * std::numbers::pi_v<PDFloat> *
						(static_cast<PDFloat>(1.0) + (nb_RAD2[i] - RAD2) / (nb_dist[i] * nb_dist[i]))) * e[i];
				}
				if (withVol || withDVol) {
					scone = sintheta[i] * sintheta[i] * static_cast<PDFloat>(2.0) * std::numbers::pi_v<PDFloat>;
				}
				if (withVol)   vol2 += costheta[i] * scone;
				if (withDVol)  DVol[iatom] = DVol[iatom] - (0.5*RAD2*scone) * e[i];
			}

		}

	if (withSasa) Sasa[iatom] = RAD2 * sasa_ia;
	if (withDSasa)
	{
		for (i = 0; i < nnb; ++i)
		{
			DSasa[iatom] -= DSasa_parts[iatom][i];
		}
	}

	if (withVol)
	{
		for (i = 0; i < nnb; ++i)
		{
			if (fknot[i]  != 0) vol3 += RAD * volnb[i] * costheta[i];
		}
	}
	if (withVol)   Vol[iatom] = RAD*RAD2*sasa_ia/3.0 + RAD*RAD2*vol2/6.0 + vol3/6.0;

	for (i = 0; i < nnb; ++i)       // over neighbors, set zero again (only necessary if an expansion of power diagram is planned)
	{
		neighbour_mut(i).visitedAs = 0;
	}
	return;
}

//=============================================================


template <class PDFloat, class PDCoord> void PowerSasa<PDFloat, PDCoord>::
calc_sasa_all()
{
	for (std::size_t i = 0; i < power_diagram->get_points().size(); ++i)
	{
		calc_sasa_single(static_cast<unsigned int>(i));
	}
}


}
#endif /* SASA_H_ */
