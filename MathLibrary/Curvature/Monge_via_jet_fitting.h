#pragma once

#include <iterator>
#include <math.h>
#include "SurfaceMesh/Vector.h"
#include "Macros.h"

#include "Eigen/LU"
#include "Eigen/Geometry"
#include "Eigen/Eigenvalues"

typedef double											FT;
typedef Vector<FT,3>									Point_3;
typedef Point_3											Vector3;
typedef std::vector<Point_3>::iterator					InputIterator;
typedef Eigen::Transform<FT, 3, Eigen::AffineCompact>	Aff_transformation;
typedef Eigen::VectorXd									LAVector;
typedef Eigen::MatrixXd									LAMatrix;

////////////////////// CLASS Monge_via_jet_fitting ////////////////////////
class Monge_via_jet_fitting {
public:
	//////////////////////begin nested CLASS Monge_form ///////////////////
	class Monge_form {
	protected:
		//point on the fitted surface where diff quantities are computed
		Vector3 m_origin_pt;
		//the monge trihedron (d1,d2,n) is orthonormal direct
		Vector3 m_d1;   //maximal ppal dir
		Vector3 m_d2;   //minimal ppal dir
		Vector3 m_n;    //normal direction
		//coeff = (k1, k2, //ppal curv
		//         b0, b1, b2, b3, //third order
		//         c0, c1, c2, c3, c4) //fourth order
		//     if (degree==1) no coeff needed
		std::vector<FT> m_coefficients;

	public:
		//constructor
		Monge_form() {
			m_origin_pt  = Point_3(0.,0.,0.); 
			m_d1 = Vector3(0.,0.,0.);
			m_d2 = Vector3(0.,0.,0.);
			m_n = Vector3(0.,0.,0.);
			m_coefficients = std::vector<FT>();
		}
		~Monge_form() {}
		//access
		const Vector3 origin() const { return m_origin_pt; }
		Vector3& origin() { return m_origin_pt; }
		Vector3 maximal_principal_direction() const { return m_d1; }
		Vector3& maximal_principal_direction() { return m_d1; }
		const Vector3 minimal_principal_direction() const { return m_d2; }
		Vector3& minimal_principal_direction() { return m_d2; }
		const Vector3 Normaldirection() const { return m_n; }
		Vector3& Normaldirection() { return m_n; }
		const std::vector<FT> coefficients() const { return m_coefficients; }
		std::vector<FT>& coefficients() { return m_coefficients; }

		const FT principal_curvatures(size_t i) const {
			return coefficients()[i];}

		const FT third_order_coefficients(size_t i) const {
			return coefficients()[i+2]; }
		const FT fourth_order_coefficients(size_t i) const {
			return coefficients()[i+6]; }

		//if d>=2, number of coeffs = (d+1)(d+2)/2 -4. 
		//we remove cst, linear and the xy coeff which vanish
		void set_up(std::size_t degree);
		//switch min-max ppal curv/dir wrt a given normal orientation.
		// if given_normal.monge_normal < 0 then change the orientation
		// if z=g(x,y) in the basis (d1,d2,n) then in the basis (d2,d1,-n)
		// z=h(x,y)=-g(y,x)
		void comply_wrt_given_normal(const Vector3 given_normal);
		void dump_verbose(std::ostream& out_stream) const;
		void dump_4ogl(std::ostream& out_stream, const FT scale);
	};
	//////////////////////end nested CLASS Monge_form /////////////////////

	//continue main class Monge_via_jet_fitting ////////
public:
	Monge_via_jet_fitting(); 
	Monge_form operator()(InputIterator begin, InputIterator end,  size_t d, size_t dprime);
	const FT inverse_condition_number() const {return inverse_condition_nb;}
	const std::pair<FT, Vector3> pca_basis(size_t i) const {
		return m_pca_basis[i];}

protected:
	int deg;
	int deg_monge;
	int nb_d_jet_coeff;
	int nb_input_pts;
	FT preconditionning;
	FT inverse_condition_nb;

	std::vector< std::pair<FT, Vector3> > m_pca_basis;

	//translate_p0 changes the origin of the world to p0 the first point 
	//  of the input data points
	//change_world2fitting (coord of a vector in world) = coord of this 
	//  vector in fitting. The matrix transform has as lines the coord of
	//  the basis vectors of fitting in the world coord. 
	//idem for change_fitting2monge
	Aff_transformation translate_p0, change_world2fitting, change_fitting2monge;

	//eigen val and vector stored in m_pca_basis
	// change_world2fitting is computed 
	void compute_PCA(InputIterator begin, InputIterator end); 

	//Coordinates of input points are computed in the fitting basis with 
	//  p0 as origin.
	//Pre-conditioning is computed, M and Z are filled
	void fill_matrix(InputIterator begin, InputIterator end, std::size_t d, LAMatrix& M, LAVector& Z);
	//A is computed, solving MA=Z in the ls sense, the solution A is stored in Z
	//Pre-conditioning is needed
	void solve_linear_system(LAMatrix &M, LAVector& Z);

	//Classical differential geometric calculus
	//change_fitting2monge is computed
	//if deg_monge =1 only 1st order info
	//if deg_monge >= 2 2nd order info are computed
	void compute_Monge_basis(const FT* A, Monge_form& monge_form);

	//if deg_monge >=3 then 3rd (and 4th) order info are computed
	void compute_Monge_coefficients(FT* A, std::size_t dprime, 
		Monge_form& monge_form);

	//for a trihedron (v1,v2,v3) switches v1 to -v1 if det(v1,v2,v3) < 0
	void switch_to_direct_orientation(Vector3& v1, const Vector3& v2,
		const Vector3& v3);

	//Convert Vector3 to Eigen::Vector3d back and forth
	Vector3 E2V_converter(Eigen::Vector3d vec);
	Eigen::Vector3d V2E_converter(Vector3 vec);

	friend std::ostream& operator<< (std::ostream& out_stream, const Monge_via_jet_fitting::Monge_form& monge)
	{
		monge.dump_verbose(out_stream);
		return out_stream;
	}

};


