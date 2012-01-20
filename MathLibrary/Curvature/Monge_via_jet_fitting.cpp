#include "Monge_via_jet_fitting.h"

using namespace std;

void Monge_via_jet_fitting::Monge_form::set_up(std::size_t degree) {
	if ( degree >= 2 ) std::fill_n(std::back_inserter(m_coefficients),
		(degree+1)*(degree+2)/2-4, 0.);
}

void Monge_via_jet_fitting::Monge_form::comply_wrt_given_normal(const Vector3 given_normal)
{
	if ( dot(given_normal,this->Normaldirection()) < 0 )
	{
		Normaldirection() = -Normaldirection();
		std::swap(maximal_principal_direction(), minimal_principal_direction());
		if ( coefficients().size() >= 2) 
			std::swap(coefficients()[0],coefficients()[1]);
		if ( coefficients().size() >= 6) {
			std::swap(coefficients()[2],coefficients()[5]);
			std::swap(coefficients()[3],coefficients()[4]);}
		if ( coefficients().size() >= 11) {
			std::swap(coefficients()[6],coefficients()[10]);
			std::swap(coefficients()[7],coefficients()[9]);}
		 std::vector<FT>::iterator itb = coefficients().begin(),
			ite = coefficients().end();
		for (;itb!=ite;itb++) { *itb = -(*itb); }
	}
}


void Monge_via_jet_fitting::Monge_form::dump_verbose(std::ostream& out_stream) const
{
	out_stream << "origin : " << origin() << std::endl
		<< "n : " << Normaldirection() << std::endl;
	if ( coefficients().size() >= 2) 
		out_stream << "d1 : " << maximal_principal_direction() << std::endl 
		<< "d2 : " << minimal_principal_direction() << std::endl
		<< "k1 : " << coefficients()[0] << std::endl 
		<< "k2 : " << coefficients()[1] << std::endl;	      
	if ( coefficients().size() >= 6) 
		out_stream << "b0 : " << coefficients()[2] << std::endl 
		<< "b1 : " << coefficients()[3] << std::endl
		<< "b2 : " << coefficients()[4] << std::endl
		<< "b3 : " << coefficients()[5] << std::endl;
	if ( coefficients().size() >= 11) 
		out_stream << "c0 : " << coefficients()[6] << std::endl 
		<< "c1 : " << coefficients()[7] << std::endl
		<< "c2 : " << coefficients()[8] << std::endl
		<< "c3 : " << coefficients()[9] << std::endl 
		<< "c4 : " << coefficients()[10] << std::endl
		<< std::endl; 
}


void Monge_via_jet_fitting::Monge_form::dump_4ogl(std::ostream& out_stream, const FT scale)
{
	out_stream << origin()  << " "
		<< maximal_principal_direction() * scale << " "
		<< minimal_principal_direction() * scale << " "
		<< coefficients()[0] << " "
		<< coefficients()[1] << " "
		<< std::endl;
}
//////////////////////////////////////////////////////////////
// Implementation main Monge_via_jet_fiting

Vector3 Monge_via_jet_fitting::E2V_converter(Eigen::Vector3d vec)
{
	return Vector3(vec[0], vec[1], vec[2]);
}

Eigen::Vector3d  Monge_via_jet_fitting::V2E_converter(Vector3 vec)
{
	return Eigen::Vector3d(vec[0], vec[1], vec[2]);
}

Monge_via_jet_fitting::	Monge_via_jet_fitting()
{
	m_pca_basis = std::vector< std::pair<FT, Vector3> >(3);
} 

Monge_via_jet_fitting::Monge_form	Monge_via_jet_fitting::operator()(InputIterator begin, InputIterator end, 
	size_t d, size_t dprime)
{
	this->deg = static_cast<int>(d);
	this->deg_monge = static_cast<int>(dprime);
	this->nb_d_jet_coeff = static_cast<int>((d+1)*(d+2)/2);
	this->nb_input_pts = static_cast<int>(end - begin);

	//Initialize
	Monge_form monge_form;
	monge_form.set_up(dprime);
	//for the system MA=Z
	LAMatrix M(nb_input_pts, nb_d_jet_coeff);
	LAVector Z = LAVector::Zero(nb_input_pts);

	compute_PCA(begin, end);
	fill_matrix(begin, end, d, M, Z);//with precond
	solve_linear_system(M, Z);  //correct with precond
	compute_Monge_basis(Z.data(), monge_form);
	if ( dprime >= 3) compute_Monge_coefficients(Z.data(), dprime, monge_form);

	//monge_form.dump_verbose(std::cout);
	//std::cout << endl;

	return monge_form;
}


void Monge_via_jet_fitting::compute_PCA(InputIterator begin, InputIterator end)
{
	int n = this->nb_input_pts;
	FT x, y, z,
		sumX = 0., sumY = 0., sumZ = 0.,
		sumX2 = 0., sumY2 = 0., sumZ2 = 0.,
		sumXY = 0., sumXZ = 0., sumYZ = 0., 
		xx, yy, zz, xy, xz, yz;

	for (; begin != end; begin++)
	{
		Point_3 lp = *begin;
		x = lp.x();
		y = lp.y();
		z = lp.z();   
		sumX += x / n;
		sumY += y / n;
		sumZ += z / n;
		sumX2 += x * x / n;
		sumY2 += y * y / n;
		sumZ2 += z * z / n;
		sumXY += x * y / n;
		sumXZ += x * z / n;
		sumYZ += y * z / n;
	}
	xx = sumX2 - sumX * sumX;
	yy = sumY2 - sumY * sumY;
	zz = sumZ2 - sumZ * sumZ;
	xy = sumXY - sumX * sumY;
	xz = sumXZ - sumX * sumZ;
	yz = sumYZ - sumY * sumZ;

	// assemble covariance matrix as a
	// semi-definite matrix. 
	// Matrix numbering:
	// 0
	// 1 2
	// 3 4 5
	Eigen::Matrix3d covariance;
	covariance << xx, xy, xz,
		xy, yy, yz,
		xz, yz, zz;

	// solve for eigenvalues and eigenvectors.
	// eigen values are sorted in descending order, 
	// eigen vectors are sorted in accordance.
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(covariance);

	Eigen::Vector3d e_values = Eigen::Reverse<Eigen::Vector3d, Eigen::Vertical>(es.eigenvalues());
	Eigen::Matrix3d e_vectors = Eigen::Reverse<Eigen::Matrix3d, Eigen::Horizontal>(es.eigenvectors());

	//store in m_pca_basis
	for (int i=0; i<3; i++)
	{
		m_pca_basis[i].first =  e_values[i];
		m_pca_basis[i].second = E2V_converter(e_vectors.col(i));	
	}

	switch_to_direct_orientation(m_pca_basis[0].second,	m_pca_basis[1].second, m_pca_basis[2].second);

	//Store the change of basis W->F
	Eigen::Matrix3d tM;
	tM << m_pca_basis[0].second[0], m_pca_basis[0].second[1], m_pca_basis[0].second[2], 
		m_pca_basis[1].second[0], m_pca_basis[1].second[1], m_pca_basis[1].second[2],
		m_pca_basis[2].second[0], m_pca_basis[2].second[1], m_pca_basis[2].second[2];

	Aff_transformation change_basis (tM.transpose());

	this->change_world2fitting = change_basis; 
}


void Monge_via_jet_fitting::fill_matrix(InputIterator begin, InputIterator end,
	std::size_t d, LAMatrix &M, LAVector& Z)
{
	//origin of fitting coord system = first input data point
	//transform coordinates of sample points with a
	//translation ($-p$) and multiplication by $ P_{W\rightarrow F}$.
	Point_3 origin = *begin;

	this->translate_p0 = Aff_transformation(Eigen::Translation3d(V2E_converter(-origin)));

	// Rotate then translate
	Aff_transformation transf_points = this->change_world2fitting.inverse() * translate_p0;

	//compute and store transformed points
	std::vector<Point_3> pts_in_fitting_basis;
	for(;begin != end; begin++)
	{
		Vector3 p_orig = *begin;
		Vector3 p = E2V_converter(transf_points * V2E_converter(*begin));

		pts_in_fitting_basis.push_back(p);

		/*cout << "From: " << p_orig[0] << " " << p_orig[1] << " " << p_orig[2] << endl;
		cout << "To: " << p[0] << " " << p[1] << " " << p[2] << endl;
		cout << endl;*/
	}

	//Compute preconditionning
	FT precond = 0.;
	 std::vector<Point_3>::iterator itb = pts_in_fitting_basis.begin(),
		ite = pts_in_fitting_basis.end();

	for(;itb != ite; itb++){ 
		precond += abs(itb->x()) + abs(itb->y());
	}

	precond /= 2*this->nb_input_pts;
	this->preconditionning = precond;

	//fill matrices M and Z
	itb = pts_in_fitting_basis.begin();
	int line_count = 0;
	FT x, y;
	for(;itb != ite; itb++)
	{
		x = itb->x();
		y = itb->y();

		Z(line_count) = itb->z();

		for (std::size_t k=0; k <= d; k++) {
			for (std::size_t i=0; i<=k; i++) {

				FT g = std::pow(x,static_cast<double>(k-i)) * std::pow(y,static_cast<double>(i))
					/( fact(i) * fact(k-i) *std::pow(this->preconditionning,static_cast<double>(k)));

				size_t c = k*(k+1)/2+i;

				M(line_count, c) = g;
			}
		}

		line_count++;
	}
}

void Monge_via_jet_fitting:: solve_linear_system(LAMatrix &M, LAVector& Z)
{
	Eigen::VectorXd A = M.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Z);

	for (int k=0; k <= this->deg; k++) 
	{
		for (int i=0; i<=k; i++)
		{
			int id = k*(k+1)/2+i;

			FT z_entry = A(id) / std::pow(this->preconditionning,k);
			Z(id) = z_entry;
		}
	}

	Eigen::JacobiSVD<Eigen::MatrixXd::PlainObject>::SingularValuesType sing_vals;
	sing_vals = M.jacobiSvd().singularValues();

	this->inverse_condition_nb = sing_vals(sing_vals.size()-1) / sing_vals(0);
}

void Monge_via_jet_fitting::compute_Monge_basis(const FT* A, Monge_form& monge_form)
{
	// only 1st order info.
	if ( this->deg_monge == 1 ) {  
		/*Eigen::Vector3d orig_monge(0., 0., A[0]);
		Eigen::Vector3d  normal(-A[1], -A[2], 1.);
		FT norm2 = normal.dot(normal);
		normal = normal / sqrt(norm2);
		monge_form.origin() = E2V_converter((translate_p0.inverse() * change_world2fitting.inverse()) * (orig_monge)) ;
		monge_form.Normaldirection() = E2V_converter(this->change_world2fitting.inverse() * (normal));*/
	}
	// else (deg_monge >= 2) then 2nd order info are computed
	else 
	{
		//bi-index to uni-index conversion : A(i,j)=A[(i+j)(i+j+1)/2+j]
		Eigen::Vector3d orig_monge(0., 0., A[0]);

		//normal = Xu cross product Xv
		Eigen::Vector3d Xu(1, 0, A[1]), Xv(0, 1, A[2]);
		Eigen::Vector3d normal(-A[1], -A[2], 1);

		FT norm2 = normal.dot(normal);
		normal = normal / sqrt(norm2);

		//Surface in fitting_basis : X(u,v)=(u,v,J_A(u,v))
		//in the basis Xu=(1,0,A[1]), Xv=(0,1,A[2]), Weingarten=-I^{-1}II
		//first fond form I=(e,f,f,g)
		//                 =(Xu.Xu, Xu.Xv, Xu.Xv, Xv.Xv)
		//second fond form II=(l,m,m,n)/norm2^(1/2)
		//                   =(n.Xuu, n.Xuv, n.Xuv, n.Xvv)
		//ppal curv are the opposite of the eigenvalues of Weingarten or the
		//  eigenvalues of weingarten = -Weingarten = I^{-1}II

		FT e = 1+A[1]*A[1], f = A[1]*A[2], g = 1+A[2]*A[2], l = A[3], m = A[4], n = A[5];
		Eigen::Matrix2d  weingarten;
		weingarten << (g*l-f*m)/ (sqrt(norm2)*norm2), (g*m-f*n)/ (sqrt(norm2)*norm2),
					  (e*m-f*l)/ (sqrt(norm2)*norm2), (e*n-f*m)/ (sqrt(norm2)*norm2);
		// Y, Z are normalized GramSchmidt of Xu, Xv
		// Xu->Y=Xu/||Xu||;
		// Xv->Z=Xv-(Xu.Xv)Xu/||Xu||^2;
		// Z-> Z/||Z||
		Eigen::Vector3d Y, Z;
		FT normXu = Xu.norm();
		Y = Xu / normXu;
		FT XudotXv = Xu.dot(Xv);
		Z = Xv - XudotXv * Xu / (normXu*normXu);
		FT normZ = Z.norm();
		Z = Z / normZ;
		Eigen::Matrix2d change_XuXv2YZ;
		change_XuXv2YZ(0,0) = 1 / normXu;
		change_XuXv2YZ(0,1) = -XudotXv / (normXu * normXu * normZ);
		change_XuXv2YZ(1,0) = 0;
		change_XuXv2YZ(1,1) = 1 / normZ;

		Eigen::Matrix2d inv = change_XuXv2YZ.inverse();

		//in the new orthonormal basis (Y,Z) of the tangent plane :
		weingarten = inv * weingarten * change_XuXv2YZ;


		Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(weingarten);

		Eigen::Vector2d e_values = Eigen::Reverse<Eigen::Vector2d, Eigen::Vertical>(es.eigenvalues());
		Eigen::Matrix2d e_vectors = Eigen::Reverse<Eigen::Matrix2d, Eigen::Horizontal>(es.eigenvectors());

		const FT * eval = e_values.data();
		const FT * evec = e_vectors.data();		//eval in decreasing order


		Vector3 d_max = E2V_converter(evec[0]*Y + evec[1]*Z);
		Vector3 d_min = E2V_converter(evec[2]*Y + evec[3]*Z);

		switch_to_direct_orientation(d_max, d_min, E2V_converter(normal));
		
		Eigen::Matrix3d tM;
		tM <<	d_max[0], d_max[1], d_max[2], 
				d_min[0], d_min[1], d_min[2],
				normal[0], normal[1], normal[2];

		Aff_transformation change_basis (tM.transpose());
		
		this->change_fitting2monge = change_basis;

		//store the monge basis origin and vectors with their world coord
		//store ppal curv

		monge_form.origin() = E2V_converter((translate_p0.inverse() * change_world2fitting) * orig_monge);

		monge_form.maximal_principal_direction() = E2V_converter(change_world2fitting * V2E_converter(d_max));
		monge_form.minimal_principal_direction() = E2V_converter(change_world2fitting * V2E_converter(d_min));
		
		monge_form.Normaldirection() = E2V_converter(change_world2fitting * normal);
		
		monge_form.coefficients()[0] = eval[0];
		monge_form.coefficients()[1] = eval[1];
	}
	//end else
}


void Monge_via_jet_fitting::compute_Monge_coefficients(FT* A, std::size_t dprime, Monge_form& monge_form)
{
	// One has the equation w=J_A(u,v) of the fitted surface S in the fitting_basis
	// Substituing (u,v,w)=change_fitting2monge^{-1}(x,y,z)
	// One has the equation f(x,y,z)=0 on this surface S in the monge basis
	// The monge form of the surface at the origin is the bivariate fct
	//   g(x,y) s.t. f(x,y,g(x,y))=0
	// voir les calculs Maple dans monge.mws
	// Notations are f123= d^3f/dxdydz
	//              g(x,y)=sum (gij x^i y^j/ i!j!) with 
	//              g00=g10=g01=g11=0, g20=kmax, g02=kmin
	//
	// g(x,y)= 1/2*(k1x^2 +k2y^2) 
	//       +1/6*(b0x^3 +3b1x^2y +3b2xy^2 +b3y^3) 
	//       +1/24*(c0x^4 +4c1x^3y +6c2x^2y^2 +4c3xy^3 +c4y^4)
	//       +...
	// p stores change_fitting2monge^{-1}=change_fitting2monge^{T}
	FT p[3][3];

	Aff_transformation t = change_fitting2monge.inverse();

	p[0][0] = t(0,0);
	p[1][0] = t(0,1);
	p[2][0] = t(0,2);
	p[0][1] = t(1,0);
	p[1][1] = t(1,1);
	p[2][1] = t(1,2);
	p[0][2] = t(2,0);
	p[1][2] = t(2,1);
	p[2][2] = t(2,2);

	// formula are designed for w=sum( Aij ui vj), but we have J_A = sum( Aij/i!j! ui vj)
	for (int k=0; k <= this->deg; k++) 
		for (int i=0; i<=k; i++)
			A[k*(k+1)/2+i] /= fact(k-i)*fact(i);//this is A(k-i;i)

	/*   //debug */
	/*   std::cout << "coeff of A" << std::endl */
	/* 	    << A[0] << " "<< A[1] << " "<< A[2] << std::endl */
	/* 	    << A[3] << " "<< A[4] << " "<< A[5] << std::endl */
	/* 	    << A[6] << " "<< A[7] << " "<< A[8] << " "<< A[9]<< std::endl */
	/* 	    << A[10] << " "<< A[11] << " "<< A[12] << " "<< A[13]<< " " << A[14] << std::endl; */



	//     note f1 = f2 = f12 = 0 
	//     FT f1 = A[1] * p[0][0] + A[2] * p[1][0] - p[2][0];
	//     FT f2 = A[2] * p[1][1] + A[1] * p[0][1] - p[2][1];
	//     FT f12 = 
	//     2 * A[3] * p[0][0] * p[0][1]
	//     + 2 * A[5] * p[1][0] * p[1][1]
	//     + A[4] * p[0][1] * p[1][0] 
	//     + A[4] * p[0][0] * p[1][1];
	//         -f11 / f3 = kmax
	//         -f22 / f3 = kmin 

	FT f3 = A[1] * p[0][2] + A[2] * p[1][2] - p[2][2];
	FT f11 =
		2 * A[4] * p[0][0] * p[1][0]
	+ 2 * A[5] * p[1][0] * p[1][0]
	+ 2 * A[3] * p[0][0] * p[0][0];
	FT f13 =
		A[4] * p[0][0] * p[1][2]
	+ A[4] * p[0][2] * p[1][0]
	+ 2 * A[5] * p[1][0] * p[1][2]
	+ 2 * A[3] * p[0][0] * p[0][2];
	FT f22 =
		2 * A[4] * p[0][1] * p[1][1]
	+ 2 * A[5] * p[1][1] * p[1][1]
	+ 2 * A[3] * p[0][1] * p[0][1];
	FT f23 =
		A[4] * p[0][1] * p[1][2]
	+ 2 * A[5] * p[1][1] * p[1][2]
	+ A[4] * p[0][2] * p[1][1]
	+ 2 * A[3] * p[0][1] * p[0][2];
	FT f33 =
		2 * A[5] * p[1][2] * p[1][2]
	+ 2 * A[3] * p[0][2] * p[0][2]
	+ 2 * A[4] * p[0][2] * p[1][2];
	FT f111 =
		6 * A[8] * p[0][0] * p[1][0] * p[1][0]
	+ 6 * A[7] * p[0][0] * p[0][0] * p[1][0]
	+ 6 * A[6] * p[0][0] * p[0][0] * p[0][0]
	+ 6 * A[9] * p[1][0] * p[1][0] * p[1][0];
	FT f222 =
		6 * A[7] * p[0][1] * p[0][1] * p[1][1]
	+ 6 * A[8] * p[0][1] * p[1][1] * p[1][1]
	+ 6 * A[9] * p[1][1] * p[1][1] * p[1][1]
	+ 6 * A[6] * p[0][1] * p[0][1] * p[0][1];
	FT f112 =
		2 * A[7] * p[0][0] * p[0][0] * p[1][1]
	+ 6 * A[6] * p[0][0] * p[0][0] * p[0][1]
	+ 2 * A[8] * p[0][1] * p[1][0] * p[1][0]
	+ 4 * A[8] * p[0][0] * p[1][0] * p[1][1]
	+ 6 * A[9] * p[1][0] * p[1][0] * p[1][1]
	+ 4 * A[7] * p[0][0] * p[0][1] * p[1][0];
	FT f122 =
		4 * A[8] * p[0][1] * p[1][0] * p[1][1]
	+ 2 * A[8] * p[0][0] * p[1][1] * p[1][1]
	+ 6 * A[6] * p[0][0] * p[0][1] * p[0][1]
	+ 2 * A[7] * p[0][1] * p[0][1] * p[1][0]
	+ 4 * A[7] * p[0][0] * p[0][1] * p[1][1]
	+ 6 * A[9] * p[1][0] * p[1][1] * p[1][1];
	FT f113 = 
		6*A[6]*p[0][0]*p[0][0]*p[0][2]
	+6*A[9]*p[1][0]*p[1][0]*p[1][2]
	+2*A[7]*p[0][0]*p[0][0]*p[1][2]
	+2*A[8]*p[0][2]*p[1][0]*p[1][0]
	+4*A[7]*p[0][0]*p[0][2]*p[1][0]
	+4*A[8]*p[0][0]*p[1][0]*p[1][2];
	FT f223 = 
		2*A[8]*p[0][2]*p[1][1]*p[1][1]
	+6*A[6]*p[0][1]*p[0][1]*p[0][2]
	+6*A[9]*p[1][1]*p[1][1]*p[1][2]
	+2*A[7]*p[0][1]*p[0][1]*p[1][2]
	+4*A[7]*p[0][1]*p[0][2]*p[1][1]
	+4*A[8]*p[0][1]*p[1][1]*p[1][2];
	FT f123 = 
		2*A[8]*p[0][2]*p[1][0]*p[1][1]
	+2*A[7]*p[0][0]*p[0][1]*p[1][2]
	+2*A[7]*p[0][0]*p[0][2]*p[1][1]
	+6*A[9]*p[1][0]*p[1][1]*p[1][2]
	+2*A[7]*p[0][1]*p[0][2]*p[1][0]
	+6*A[6]*p[0][0]*p[0][1]*p[0][2]
	+2*A[8]*p[0][0]*p[1][1]*p[1][2]
	+2*A[8]*p[0][1]*p[1][0]*p[1][2];

	FT b0 = 1/(f3*f3)*(-f111*f3+3*f13*f11);
	FT b1 = 1/(f3*f3)*(-f112*f3+f23*f11);
	FT b2 = 1/(f3*f3)*(-f122*f3+f13*f22);
	FT b3 = -1/(f3*f3)*(f222*f3-3*f23*f22);

	monge_form.coefficients()[2] = b0;
	monge_form.coefficients()[3] = b1;
	monge_form.coefficients()[4] = b2;
	monge_form.coefficients()[5] = b3;

	if ( dprime == 4 )
	{
		FT f1111 = 
			24*A[13]*p[0][0]*p[1][0]*p[1][0]*p[1][0]
		+24*A[12]*p[0][0]*p[0][0]*p[1][0]*p[1][0]
		+24*A[11]*p[0][0]*p[0][0]*p[0][0]*p[1][0]
		+24*A[14]*p[1][0]*p[1][0]*p[1][0]*p[1][0]
		+24*A[10]*p[0][0]*p[0][0]*p[0][0]*p[0][0];
		FT f1112 = 
			6*A[13]*p[0][1]*p[1][0]*p[1][0]*p[1][0]
		+18*A[13]*p[0][0]*p[1][0]*p[1][0]*p[1][1]
		+24*A[10]*p[0][0]*p[0][0]*p[0][0]*p[0][1]
		+12*A[12]*p[0][0]*p[0][1]*p[1][0]*p[1][0]
		+18*A[11]*p[0][0]*p[0][0]*p[0][1]*p[1][0]
		+24*A[14]*p[1][0]*p[1][0]*p[1][0]*p[1][1]
		+6*A[11]*p[0][0]*p[0][0]*p[0][0]*p[1][1]
		+12*A[12]*p[0][0]*p[0][0]*p[1][0]*p[1][1];
		FT f1122 = 
			12*A[11]*p[0][0]*p[0][0]*p[0][1]*p[1][1]
		+12*A[13]*p[0][0]*p[1][0]*p[1][1]*p[1][1]
		+12*A[13]*p[0][1]*p[1][0]*p[1][0]*p[1][1]
		+16*A[12]*p[0][0]*p[0][1]*p[1][0]*p[1][1]
		+12*A[11]*p[0][0]*p[0][1]*p[0][1]*p[1][0]
		+24*A[10]*p[0][0]*p[0][0]*p[0][1]*p[0][1]
		+4*A[12]*p[0][1]*p[0][1]*p[1][0]*p[1][0]
		+4*A[12]*p[0][0]*p[0][0]*p[1][1]*p[1][1]
		+24*A[14]*p[1][0]*p[1][0]*p[1][1]*p[1][1];
		FT f1222 = 
			6*A[13]*p[0][0]*p[1][1]*p[1][1]*p[1][1]
		+24*A[10]*p[0][0]*p[0][1]*p[0][1]*p[0][1]
		+24*A[14]*p[1][0]*p[1][1]*p[1][1]*p[1][1]
		+6*A[11]*p[0][1]*p[0][1]*p[0][1]*p[1][0]
		+18*A[11]*p[0][0]*p[0][1]*p[0][1]*p[1][1]
		+12*A[12]*p[0][0]*p[0][1]*p[1][1]*p[1][1]
		+12*A[12]*p[0][1]*p[0][1]*p[1][0]*p[1][1]
		+18*A[13]*p[0][1]*p[1][0]*p[1][1]*p[1][1];
		FT f2222 =
			24*A[13]*p[0][1]*p[1][1]*p[1][1]*p[1][1]
		+24*A[11]*p[0][1]*p[0][1]*p[0][1]*p[1][1]
		+24*A[12]*p[0][1]*p[0][1]*p[1][1]*p[1][1]
		+24*A[10]*p[0][1]*p[0][1]*p[0][1]*p[0][1]
		+24*A[14]*p[1][1]*p[1][1]*p[1][1]*p[1][1];

		FT c0 =	-1/(f3*f3*f3)*(f1111*(f3*f3)-4*f13*f3*f111+12*f13*f13*f11-6*f113*f3*f11+3*f33*f11*f11);
		FT c1 =	1/(f3*f3*f3)*(f23*f3*f111+3*f3*f123*f11+3*f13*f3*f112-f1112*(f3*f3)-6*f13*f23*f11); 
		FT c2 =	1/(f3*f3*f3)*(-f33*f22*f11+f113*f3*f22+2*f13*f3*f122-2*f13*f13*f22+f223*f3*f11+2*f23*f3*f112-2*f23*f23*f11-f1122*(f3*f3)); 
		FT c3 =	1/(f3*f3*f3)*(-f1222*(f3*f3)-6*f13*f23*f22+3*f123*f3*f22+f13*f3*f222+3*f23*f3*f122); 
		FT c4 =	-1/(f3*f3*f3)*(f2222*(f3*f3)+3*f33*f22*f22-6*f223*f3*f22-4*f23*f3*f222+12*f23*f23*f22) ; 

		monge_form.coefficients()[6] = c0;
		monge_form.coefficients()[7] = c1;
		monge_form.coefficients()[8] = c2;
		monge_form.coefficients()[9] = c3;
		monge_form.coefficients()[10] = c4;
	}
}

void Monge_via_jet_fitting::switch_to_direct_orientation(Vector3& v1, const Vector3& v2,	const Vector3& v3) 
{
	Eigen::Matrix3d M;

	for (int i=0; i<3; i++) M(i,0) = v1[i];
	for (int i=0; i<3; i++) M(i,1) = v2[i];
	for (int i=0; i<3; i++) M(i,2) = v3[i];

	//if(dot(cross(v1, v2), v3))
	//	v1 = -v1;

	if (M.determinant() < 0) 
		v1 = -v1;
}
