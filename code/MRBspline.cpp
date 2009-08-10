//#include <stdio.h>
//#include <stdlib.h>
#include "stdafx.h" 
#include <iostream>
//#include "globals.h"
#include <cmath>
#include "MRBspline.h"
#include "globals.h"
//#include "TIP.h"
//#include "TIPDoc.h"
//#include "TIPView.h"
//#include "MRBspline.h"

////////////////////////////////////////////////////////////////
/// IMPORTANT!!!!!
/// This "using" directive should "follow" all the "include"'s!!!
using namespace std; // Resolves any conflicts among 'include' file definitions
//////////////////////////////////////////////////////////////

//#define MAX_PTS	500
//#define MAX_CNT_PTS	500
#define MAX_PTS	1000
#define MAX_CNT_PTS	1000

//typedef double		Vec[MAX_CNT_PTS][3]; /* level n, 2^n + 3 pnts possible */
//typedef double		Mat[MAX_CNT_PTS][MAX_CNT_PTS];

#define dist3(x, y, z) sqrt(((double)x)*((double)x) + ((double)y)*((double)y) + ((double)z)*((double)z))
inline double dist2(double x1, double y1, double x2, double y2)
{
	return ( sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) ) );
}

void vector::print()
{
	for (int i = 0; i < N; i++) {
		TRACE("[%d] = %.2f  ", i, p[i]);
	}
	TRACE("\n");
}

void vector::copy(vector& b)
{
	delete[] p;
	N = b.N;
	p = new double[N];
	for (int i = 0; i < N; i++) {
		p[i] = b[i];
	}
}

double operator*(const vector& a, const vector& b)
{
	if (a.N != b.N) {
		TRACE("VECTOR INNER PRODUCT: a.N != b.N\n");
		exit(1);
	}
	double c = 0.0;
	for (int i = 0; i < a.N; i++) {
		c += a[i] * b[i];
	}
	return c;
}

void VecVecAdd( const vector& a, const vector& b, vector& c )
{
	c.init(a.N);

	for (int i = 0; i < a.N; i++) {
		c[i] = a[i] + b[i];
	}
}

matrix matrix::copy()
{
	matrix b(Nr, Nc);
	for (int i = 0; i < Nr; i++) {
		for (int j = 0; j < Nc; j++) {
			b[i][j] = p[i][j];
		}
	}

    return b;
}

void matrix::copy(matrix& b)
{
	init(b.Nr, b.Nc);
	for (int i = 0; i < Nr; i++) 
		for (int j = 0; j < Nc; j++) 
			p[i][j] = b.p[i][j];

	//p[0][0] = 0;
	/*
	for (int i = 0; i < Nr; i++) {
		for (int j = 0; j < Nc; j++) {
			//p[i][j] = b[i][j];
			//get(i, j) = b[i][j];
			//p[i][j] = b.get(i, j);
			//p[i][j] = 1;
		}
	}
	*/
}

vector operator*( matrix const& a, vector const& b )
{
    vector c(a.Nr);

	if (a.Nc != b.N) {
		TRACE("matrix.Nc = %d, vector.N = %d\n", a.Nc, b.N);
		TRACE("Matrix-Vector Mult Error!\n");
		//cerr << "Matrix-Vector Mult Error!\n";
		exit(1);
	}
	
	for (int i = 0; i < a.Nr; i++) {
		c[i] = 0.0;
		for (int j = 0; j < a.Nc; j++) {
			c[i] += a.get(i,j) * b.get(j);
			///////////////////////////////////////////////////////////
			//c[i] += a[i][j] * b[j]; // this doesn't work, why?
		}
		//TRACE("c[%d] = %.3f\n", i, c[i]);
	}
	
	//c.print();

    return c;
}

void MatVecMult( matrix const& a, vector const& b, vector& c)
{
	c.init(a.Nr);

	if (a.Nc != b.N) {
		TRACE("matrix.Nc = %d, vector.N = %d\n", a.Nc, b.N);
		TRACE("Matrix-Vector Mult Error!\n");
		//cerr << "Matrix-Vector Mult Error!\n";
		exit(1);
	}
	
	for (int i = 0; i < a.Nr; i++) {
		c[i] = 0.0;
		for (int j = 0; j < a.Nc; j++) {
			c[i] += a.get(i,j) * b.get(j);
			///////////////////////////////////////////////////////////
			//c[i] += a[i][j] * b[j]; // this doesn't work, why?
		}
		//TRACE("c[%d] = %.3f\n", i, c[i]);
	}
	
	//c.print();
    //return c;
	
}

void MatMatMult( matrix const& a, matrix const& b, matrix& c)
{
	c.init(a.Nr, b.Nc);

	if (a.Nc != b.Nr) {
		TRACE("a.Nc = %d, b.Nr = %d\n", a.Nc, b.Nr);
		TRACE("Matrix-Matrix Mult Error!\n");
		//cerr << "Matrix-Vector Mult Error!\n";
		exit(1);
	}
	
	for (int i = 0; i < a.Nr; i++) {
		for (int j = 0; j < b.Nc; j++) {
			c[i][j] = 0.0;
			for (int k = 0; k < a.Nc; k++) {
				c[i][j] += a.get(i,k) * b.get(k,j);
				///////////////////////////////////////////////////////////
				//c[i] += a[i][j] * b[j]; // this doesn't work, why?
			}
		}
		//TRACE("c[%d] = %.3f\n", i, c[i]);
	}
	
	//c.print();
    //return c;
	
}

void MatMatAdd( matrix const& a, matrix const& b, matrix& c)
{
	c.init(a.Nr, b.Nc);

	if (a.Nc != b.Nr) {
		TRACE("a.Nr = %d, b.Nr = %d\n", a.Nr, b.Nr);
		TRACE("a.Nc = %d, b.Nc = %d\n", a.Nc, b.Nc);
		TRACE("Matrix-Matrix Add Error!\n");
		exit(1);
	}
	
	for (int i = 0; i < a.Nr; i++) {
		for (int j = 0; j < a.Nc; j++) {
			c[i][j] = a.get(i,j) + b.get(i,j);
		}
		//TRACE("c[%d] = %.3f\n", i, c[i]);
	}
	
	//c.print();
    //return c;
	
}


void MatVecMultOverwrite( matrix const& a, vector& b)
// Overwrite the result in vector b
{
	vector c(a.Nr);

	if (a.Nc != b.N) {
		TRACE("matrix.Nc = %d, vector.N = %d\n", a.Nc, b.N);
		TRACE("Matrix-Vector Mult Overwrite Error!\n");
		//cerr << "Matrix-Vector Mult Error!\n";
		exit(1);
	}
	
	for (int i = 0; i < a.Nr; i++) {
		c[i] = 0.0;
		for (int j = 0; j < a.Nc; j++) {
			c[i] += a.get(i,j) * b.get(j);
			///////////////////////////////////////////////////////////
			//c[i] += a[i][j] * b[j]; // this doesn't work, why?
		}
		//TRACE("c[%d] = %.3f\n", i, c[i]);
	}
	
	b.copy(c);
    //return c;
}

matrix operator*( matrix const& a, matrix const& b)
/////////////////////////////////////////////////////
// IMPORTANT!!!
// The reason why A = B * C didn't work in my case was because
// my matrix included "dynamically allocated storage!!!".
// But, in Jehee's matrix definition, there was no pointers or dynamically allocated storage
// in the matrix definition. So the assignment statement like "A = B * C" works without any problem,
// because the compiler automatically provides copy constructor and assignment operator for 
// static-member-only classes!!! 
// So to solve this problem, YOU NEED TO DEFINE YOUR OWN VERSION OF ASSIGNMENT OPERATOR!!!
// (Or user-defined copy constructor!)
{
	if (a.Nc != b.Nr) {
		TRACE("Matrix-Matrix Mult Error!\n");
		//cerr << "Matrix-Vector Mult Error!\n";
		exit(1);
	}
	matrix c(a.Nr, b.Nc);
	for (int i = 0; i < a.Nr; i++) {
		for (int k = 0; k < b.Nc; k++) {
			c[i][k] = 0;
			for (int j = 0; j < a.Nc; j++) {
				c[i][k] += a.get(i,j) * b.get(j,k);
				///////////////////////////////////////////////////////////
				//c[i][k] += a[i][j] * b[j][k]; // this doesn't work, why?
			}
		}
	}

    return c; // there is no problem with this code!!!
}

#define fabs(x) ((x) >= 0 ? (x) : -(x))

void ludcmp(matrix& a, int_vector& indx, double &d)
{
	const double TINY = 1.0e-20; // a very small number
	int i, imax, j, k;
	double big, dum, sum, temp;
	
	int n = a.Nr; // number of rows
	vector vv = vector(n); // vector of n elements;
	//vector vv(n); // vector of n elements;
	d = 1.0;
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++)
			if ((temp = fabs(a[i][j])) > big)
				big = temp;
		if (big == 0.0) {
			cout << "Singular matrix in routine ludcmp!\n";
			TRACE("Singular matrix in routine ludcmp!\n");
			exit(1);
		}
		vv[i] = 1.0 / big;
	}
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			d = -d;
			vv[imax] = vv[j];
		}
		///////////////////////////
		// double vector -> int vector
		indx[j] = imax; // here is a problem 
		if (a[j][j] == 0.0)
			a[j][j] = TINY;
		if (j != n-1) {
			dum = 1.0 / (a[j][j]);
			for (i = j + 1; i < n; i++)
			a[i][j] *= dum;
		}
	}
	//free_vector(vv, 1, n);
}

void lubksb(matrix& a, int_vector& indx, vector &b)
{
  int i, ii = 0, ip, j;
  double sum;

  int n = a.Nr;
  for (i = 0; i < n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii != 0)
      for (j = ii-1; j < i; j++)
        sum -= a[i][j] * b[j];
    else if (sum != 0.0)
      ii = i+1;
    b[i] = sum;
  }
  for (i = n-1; i >= 0; i--) {
    sum = b[i];
    for (j = i + 1; j < n; j++)
      sum -= a[i][j] * b[j];
    b[i] = sum / a[i][i];
  }
}

void inverse(matrix &a, matrix &inverse) // Destroys matrix a!!!
{
	double d;
	int	i, j, n;

	n = a.Nr;
	int_vector indx(n);
	vector col(n);
	inverse.init(n, n);
	ludcmp(a, indx, d); // This way, matrix 'a' will be destroyed!
	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) col[i] = 0.0;
		col[j] = 1.0;
		lubksb(a, indx, col);
		for (i = 0; i < n; i++) inverse[i][j] = col[i];
	}
	
}

void inverse2(matrix &a, matrix &inverse) // does not destroy original matrix a!!!
{
	double d;
	int	i, j, n;
	matrix a_copy;

	a_copy = a.copy();

	n = a.Nr;
	int_vector indx(n);
	vector col(n);
	inverse.init(n, n);
	ludcmp(a_copy, indx, d); // This way, matrix 'a' will be destroyed!
	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) col[i] = 0.0;
		col[j] = 1.0;
		lubksb(a_copy, indx, col);
		for (i = 0; i < n; i++) inverse[i][j] = col[i];
	}
	
}

void transpose(matrix &a, matrix &transpose)
{
	int	i, j;

	//int_vector indx(n);
	//vector col(n);
	transpose.init(a.Nc, a.Nr);
	for (j = 0; j < a.Nc; j++) {
		for (i = 0; i < a.Nr; i++) {
			transpose[j][i] = a[i][j];
		}
	}
	
}


void merge_h(matrix &a, matrix &b, matrix &merged)
{
	int	i, j;

	if (a.Nr != b.Nr) {
		TRACE("matrix merge error!\n");
		exit(1);
	}
	merged.init(a.Nr, a.Nc + b.Nc);
	for (i = 0; i < a.Nr; i++) {
		for (j = 0; j < a.Nc; j++) 
			merged[i][j] = a[i][j];
		for (j = 0; j < b.Nc; j++) 
			merged[i][j+a.Nc] = b[i][j];
	}
}

void split_v(matrix &merged, int a_row, matrix &a, matrix &b)
{
	int	i, j;

	a.init(a_row, merged.Nc);
	b.init(merged.Nr - a_row, merged.Nc);
	for (j = 0; j < merged.Nc; j++) {
		for (i = 0; i < a.Nr; i++) 
			a[i][j] = merged[i][j]; 
		for (i = 0; i < b.Nr; i++) 
			b[i][j] = merged[i+a.Nr][j];
	}
}

matrix matrix::inverse()
{
	double d;
	int	i, j, n;

	if (Nr != Nc) {
		TRACE("matrix inverse error!\n");
		exit(1);
	}
	// First, copy this matrix to a
	// To avoid destroying original matrix
	matrix a(Nr, Nc);
	for (i = 0; i < Nr; i++) {
		for (j = 0; j < Nc; j++) {
			a[i][j] = p[i][j];
		}
	}

	n = a.Nr;
	matrix inverse(n, n);
	int_vector indx(n);
	vector col(n);
	ludcmp(a, indx, d);
	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) col[i] = 0.0;
		col[j] = 1.0;
		lubksb(a, indx, col);
		for (i = 0; i < n; i++) inverse[i][j] = col[i];
	}
	
    return inverse;
}

void matrix::print()
{
	for (int i = 0; i < Nr; i++) {
		for (int j = 0; j < Nc; j++) {
			TRACE("  [%d][%d] = %.2f", i, j, p[i][j]);
		}
		TRACE("\n");
	}
	TRACE("\n");
}

class test {
private:
	int N;
public:
	test(int n) {
		N = n;
	}
};

#define BS_1	1/6.0
#define BS_2	2/3.0

matrix P[10], Q[10], A[10], B[10]; // global variables

MRBspline::MRBspline()
// default constructor
{
	//stroke_mode = MIDDLE_THICK;
	stroke_mode = UNIFORM_THICK;
	//image_x = IMAGE_X, image_y = IMAGE_Y;
	//st_pt = end_pt = NULL;
	//total_pts = 0;
	//matrix A(2, 2); // Create 5x5 matrix
	//A[0][0] = 0.0; A[0][1] = 1.0; 
	//A[1][0] = 2.0; A[1][1] = 3.0; 
	//TRACE("A[0][0] = %f\n", A[0][0]);
	//TRACE("A[0][1] = %f\n", A[0][1]);
	//TRACE("A[1][0] = %f\n", A[1][0]);
	//TRACE("A[1][1] = %f\n", A[1][1]);
	/*
	int	N = 5;
	matrix a(N, N);
	matrix inv(N, N);
	vector b(N);
	vector indx(N);
	double d;
	
	a[0][0] = -1/2.; a[0][1] = 0; a[0][2] = 1/2.; a[0][3] = 0; a[0][4] = 0;
	a[1][0] = BS_1; a[1][1] = BS_2; a[1][2] = BS_1; a[1][3] = 0; a[1][4] = 0;
	a[2][0] = 0; a[2][1] = BS_1; a[2][2] = BS_2; a[2][3] = BS_1; a[2][4] = 0;
	a[3][0] = 0; a[3][1] = 0; a[3][2] = BS_1; a[3][3] = BS_2; a[3][4] = BS_1;
	a[4][0] = 0; a[4][1] = 0; a[4][2] = 1/2.; a[4][3] = 0; a[4][4] = -1/2.;
	matrix a2 = a.copy();
	//a2.print();
	inverse(a, inv);
	inv.print();
	matrix a3 = a2 * inv;
	a3.print();
	matrix d1(N, 2);
	d1[0][0] = 0;	d1[0][1] = 0;
	d1[1][0] = 1;	d1[1][1] = 1;
	d1[2][0] = 2;	d1[2][1] = 1;
	d1[3][0] = 3;	d1[3][1] = 2;
	d1[4][0] = 0;	d1[4][1] = 0;
	matrix c1 = inv * d1;
	//c1.print();
	*/
	/*
	int	N = 4;
	matrix a(N, N);
	matrix inv(N, N);
	vector b(N);
	vector indx(N);
	double d;
	
	a[0][0] = 0; a[0][1] = 0; a[0][2] = -1.; a[0][3] = 0;
	a[1][0] = 0; a[1][1] = 1; a[1][2] = 0; a[1][3] = 0;
	a[2][0] = 1; a[2][1] = 0; a[2][2] = 0; a[2][3] = -1;
	a[3][0] = 0; a[3][1] = 0; a[3][2] = 0; a[3][3] = 1; 
	matrix a2 = a.copy();
	//a2.print();
	inverse(a, inv);
	//inv.print();
	matrix a3 = a2 * inv;
	//a3.print();
	matrix d1(N, 1);
	d1[0][0] = 4;	
	d1[1][0] = 8;	
	d1[2][0] = -4;	
	d1[3][0] = 1;	
	matrix c1 = inv * d1;
	//matrix c1 = a2 * d1;
	c1.print();
	*/
	/*
	int	N2 = 6;
	matrix b(N2, N2);
	matrix inv2(N2, N2);
	b[0][0] = -1/2.; b[0][1] = 0; b[0][2] = 1/2.; b[0][3] = 0; b[0][4] = 0; b[0][5] = 0;
	b[1][0] = BS_1; b[1][1] = BS_2; b[1][2] = BS_1; b[1][3] = 0; b[1][4] = 0; b[1][5] = 0;
	b[2][0] = 0; b[2][1] = BS_1; b[2][2] = BS_2; b[2][3] = BS_1; b[2][4] = 0; b[2][5] = 0;
	b[3][0] = 0; b[3][1] = 0; b[3][2] = BS_1; b[3][3] = BS_2; b[3][4] = BS_1; b[3][5] = 0;
	b[4][0] = 0; b[4][1] = 0; b[4][2] = 0; b[4][3] = BS_1; b[4][4] = BS_2; b[4][5] = BS_1;
	b[5][0] = 0; b[5][1] = 0; b[5][2] = 0; b[5][3] = 1/2.; b[5][4] = 0; b[5][5] = -1/2.;
	matrix b2 = b.copy();
	//b2.print();
	inverse(b, inv2);
	inv2.print();
	matrix d2(N2, 2);
	d2[0][0] = 0;	d2[0][1] = 0;
	d2[1][0] = 1;	d2[1][1] = 1;
	d2[2][0] = 2;	d2[2][1] = 1;
	d2[3][0] = 3;	d2[3][1] = 2;
	d2[4][0] = 4;	d2[4][1] = 4;
	d2[5][0] = 0;	d2[5][1] = 0;
	matrix c2 = inv2 * d2;
	//c2.print();
	matrix b3 = b2 * inv2;
	b3.print();
	*/

	//ludcmp(a, indx, d);
	//lubksb(a, indx, b);

	/*
	// # curve segments = control points - 3 = 6 - 3 = 3
	knot.init(10); // # control points + k (k = 4 for cubic B-spline)
	knot[0] = knot[1] = knot[2] = knot[3] = 0.0;
	knot[4] = 1.0, knot[5] = 2.0; 
	knot[4] = knot[5] = knot[6] = knot[7] = 1.0;
	//cnt_pnts.init(4, 2);
	cnt_x.init(6);
	cnt_y.init(6);
	cnt_x[0] = 10.0;		cnt_y[0] = 10.0;
	cnt_x[1] = 40.0;		cnt_y[1] = 50.0;
	cnt_x[2] = 70.0;		cnt_y[2] = 50.0;
	cnt_x[3] = 100.0;		cnt_y[3] = 10.0;
	cnt_x[4] = 130.0;		cnt_y[4] = 10.0;
	cnt_x[5] = 160.0;		cnt_y[5] = 50.0;
	*/
	/*
	knot.init(11); // # control points + k (k = 4 for cubic B-spline)
	knot[0] = knot[1] = knot[2] = knot[3] = 0.0;
	knot[4] = 1.0, knot[5] = 2.0; 
	knot[6] = 3.0;
	//knot[6] = knot[7] = knot[8] = knot[9] = 3.0;
	knot[7] = knot[8] = knot[9] = knot[10] = 4.0;
	cnt_x.init(7);
	cnt_y.init(7);
	cnt_x[0] = 10.0;		cnt_y[0] = 10.0;
	cnt_x[1] = 40.0;		cnt_y[1] = 50.0;
	cnt_x[2] = 70.0;		cnt_y[2] = 50.0;
	cnt_x[3] = 100.0;		cnt_y[3] = 10.0;
	cnt_x[4] = 130.0;		cnt_y[4] = 10.0;
	cnt_x[5] = 130.0;		cnt_y[5] = 10.0;
	cnt_x[6] = 160.0;		cnt_y[6] = 50.0;
	*/
	//knot.print();
	//cnt_pnts.print();
	//DrawCurve();
	
	//GenMRMatrices();

	h.init(2);
	t.init(2);
}

double MRBspline::N(int k, int m, double t, vector& knot)
{
	double denom1, denom2, sum = 0.0;

	if (m == 1) {
		return (t >= knot[k] && t < knot[k+1]); // 1 or 0
	}
	// m exceeds 1 -> use recursion!
	denom1 = knot[k+m-1] - knot[k];
	if (denom1 != 0.0)
		sum = (t - knot[k]) * N(k, m-1, t, knot) / denom1;
	denom2 = knot[k+m] - knot[k+1];
	if (denom2 != 0.0)
		sum += (knot[k+m] - t) * N(k+1, m-1, t, knot) / denom2;
	//if (t <= 0.1)
	//	TRACE("N(%d, %d, %.3f) = %.3f\n", k, m, t, sum);
	return sum;
}

//#define floor(x)   ( (x) >=0.0 ? ((int)(x)) : (-((int)(1-(x)))) )
#define round(x) ((int) ((x) + 0.5))


void MRBspline::DrawCurve(CDC *dc, int r, int g, int b)
{
	double x, y;
	int xi, yi, i;
	int	old_x, old_y;
	double st_t, end_t;
	double N_val;
	int	pen_width;
	double KNOT_INC; 

	CPen pen;
		
	if (cnt_x.N < 4) { // can't generate cubic B-spline. Just display the polyline
		pen.CreatePen(PS_SOLID, 2, RGB(r, g, b)); // Thickened for PRESENCE paper! (from 2 to 3)
		CPen *pOldPen = (CPen *)dc->SelectObject(&pen);
		dc->MoveTo(round(cnt_x[0]), round(cnt_y[0]));
		for (i = 1; i < cnt_x.N; i++) {
			dc->LineTo(round(cnt_x[i]), round(cnt_y[i]));
		}
		dc->SelectObject(pOldPen);
		return;
	}
	pen_width = 2;
	pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	st_t = knot[0];
	end_t = knot[knot.getMax()-1];

	int polyline_count = GetPolylineCount();
	//int polyline_count = cnt_x.N;
	//KNOT_INC = (end_t - st_t) / 100.0;
	//KNOT_INC = (end_t - st_t) / (double)GetPolylineCount();
	//KNOT_INC = (end_t - st_t) / polyline_count / 5.0;
	KNOT_INC = (end_t - st_t) / polyline_count / 3.0;
	//KNOT_INC = (end_t - st_t) / polyline_count * 3.0;
	//TRACE("polyline_count = %d\n", GetPolylineCount());
	int polyline_count2 = (polyline_count > 100) ? 100: polyline_count;

	double tmp;
		
	for (double t = st_t; t < end_t; t += KNOT_INC) { // Final end_t is NOT included in the knot intervals!!!
	//for (double t = st_t; t < end_t; t += 1.0) { // Final end_t is NOT included in the knot intervals!!!
		x = y = 0.0;
		//for (int i = 0; i < cnt_x.getMax(); i++) {
		for (int i = (int)t; i <= (int)t+3; i++) { // Only 4 control points are involved
			N_val = N(i, 4, t, knot);
			x += cnt_x[i] * N_val;	
			y += cnt_y[i] * N_val;	
		}
		xi = round(x);	
		yi = round(y);	// find the nearest integer values
		if (t > st_t) {
			dc->MoveTo(xi, yi);
			dc->LineTo(old_x, old_y);
			/////////////////////////////////////////////////
			// change the line width !!!
			switch (stroke_mode) {
				case UNIFORM_THICK:
					pen_width = 3;
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case MIDDLE_THICK:
					tmp = ((t-st_t) < (end_t-t)) ? (t-st_t) : (end_t-t);
					//pen_width = (int)(3 + 10 * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/10. * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/15. * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/7. * tmp / (end_t-st_t));
					//////////////////////////////////////////////////////////////////
					pen_width = (int)(1 + polyline_count2/5. * tmp / (end_t-st_t));
                    if (pen_width > 3) pen_width = 3;
					/////////////////////////////////////////////////////////////////
					//pen_width = (int)(3 + polyline_count2/3. * tmp / (end_t-st_t));
					////////////////////////////////////////////////////////////////
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case BEGIN_THICK:
					//pen_width = (int)(3 + 10*(t-st_t)/(end_t-st_t)); // from free point to seed
					pen_width = (int)(3 + polyline_count2/17.*(t-st_t)/(end_t-st_t)); // from free point to seed
					//pen_width = (int)(3 + polyline_count2/20.*(t-st_t)/(end_t-st_t)); // from free point to seed
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case END_THICK:
					//pen_width = (int)(3 + 10*(end_t-t)/(end_t-st_t)); // from free point to seed
					//pen_width = (int)(3 + polyline_count2/20.*(end_t-t)/(end_t-st_t)); // from free point to seed
					pen_width = (int)(3 + polyline_count2/17.*(end_t-t)/(end_t-st_t)); // from free point to seed
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
			}
		}
		old_x = xi, old_y = yi;
		//TRACE("t = %.3f, [x, y] = [%.3f, %.3f]\n", t, x, y);
		//TRACE("xi = %d, yi = %d\n", xi, yi);
		//break;
	}
	// Connect to the final control point!!! (because this point is NOT included in the knot intervals!)
	xi = round(cnt_x[cnt_x.getMax()-1]);
	yi = round(cnt_y[cnt_y.getMax()-1]);
	dc->MoveTo(xi, yi);
	dc->LineTo(old_x, old_y); 
	//TRACE("# of cnt_pnts = %d\n", cnt_x.getMax());
	//TRACE("t = %.3f, [x, y] = [%.3f, %.3f] (I'm out!)\n", t, x, y);

	dc->SelectObject(pOldPen);
}

void MRBspline::DrawCurve2(CDC *dc, int r, int g, int b)
// Should use IMAGE_Y-1-cnt_y
{
	double x, y;
	int xi, yi, i;
	int	old_x, old_y;
	double st_t, end_t;
	double N_val;
	int	pen_width;
	double KNOT_INC; 

	CPen pen;
		
	if (cnt_x.N < 4) { // can't generate cubic B-spline. Just display the polyline
		pen.CreatePen(PS_SOLID, 2, RGB(r, g, b)); // Thickened for PRESENCE paper! (from 2 to 3)
		CPen *pOldPen = (CPen *)dc->SelectObject(&pen);
		dc->MoveTo(round(cnt_x[0]), round(cnt_y[0]));
		for (i = 1; i < cnt_x.N; i++) {
			dc->LineTo(round(cnt_x[i]), round(cnt_y[i]));
		}
		dc->SelectObject(pOldPen);
		return;
	}
	pen_width = 2;
	pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	st_t = knot[0];
	end_t = knot[knot.getMax()-1];

	//int polyline_count = GetPolylineCount();
	int polyline_count = cnt_x.N;
	//KNOT_INC = (end_t - st_t) / 100.0;
	//KNOT_INC = (end_t - st_t) / (double)GetPolylineCount();
	//KNOT_INC = (end_t - st_t) / polyline_count / 5.0;
	KNOT_INC = (end_t - st_t) / polyline_count / 3.0;
	//KNOT_INC = (end_t - st_t) / polyline_count * 3.0;
	//TRACE("polyline_count = %d\n", GetPolylineCount());
	int polyline_count2 = (polyline_count > 100) ? 100: polyline_count;

	double tmp;
		
	for (double t = st_t; t < end_t; t += KNOT_INC) { // Final end_t is NOT included in the knot intervals!!!
	//for (double t = st_t; t < end_t; t += 1.0) { // Final end_t is NOT included in the knot intervals!!!
		x = y = 0.0;
		//for (int i = 0; i < cnt_x.getMax(); i++) {
		for (int i = (int)t; i <= (int)t+3; i++) { // Only 4 control points are involved
			N_val = N(i, 4, t, knot);
			x += cnt_x[i] * N_val;	
			y += cnt_y[i] * N_val;	
		}
		xi = round(x);	
		yi = round(y);	// find the nearest integer values
		if (t > st_t) {
			dc->MoveTo(xi, IMAGE_Y-1-yi);
			dc->LineTo(old_x, IMAGE_Y-1-old_y);
			/////////////////////////////////////////////////
			// change the line width !!!
			switch (stroke_mode) {
				case UNIFORM_THICK:
					pen_width = 3;
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case MIDDLE_THICK:
					tmp = ((t-st_t) < (end_t-t)) ? (t-st_t) : (end_t-t);
					//pen_width = (int)(3 + 10 * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/10. * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/15. * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/7. * tmp / (end_t-st_t));
					//////////////////////////////////////////////////////////////////
					pen_width = (int)(1 + polyline_count2/5. * tmp / (end_t-st_t));
                    if (pen_width > 3) pen_width = 3;
					/////////////////////////////////////////////////////////////////
					//pen_width = (int)(3 + polyline_count2/3. * tmp / (end_t-st_t));
					////////////////////////////////////////////////////////////////
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case BEGIN_THICK:
					//pen_width = (int)(3 + 10*(t-st_t)/(end_t-st_t)); // from free point to seed
					pen_width = (int)(3 + polyline_count2/17.*(t-st_t)/(end_t-st_t)); // from free point to seed
					//pen_width = (int)(3 + polyline_count2/20.*(t-st_t)/(end_t-st_t)); // from free point to seed
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case END_THICK:
					//pen_width = (int)(3 + 10*(end_t-t)/(end_t-st_t)); // from free point to seed
					//pen_width = (int)(3 + polyline_count2/20.*(end_t-t)/(end_t-st_t)); // from free point to seed
					pen_width = (int)(3 + polyline_count2/17.*(end_t-t)/(end_t-st_t)); // from free point to seed
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
			}
		}
		old_x = xi, old_y = yi;
		//TRACE("t = %.3f, [x, y] = [%.3f, %.3f]\n", t, x, y);
		//TRACE("xi = %d, yi = %d\n", xi, yi);
		//break;
	}
	// Connect to the final control point!!! (because this point is NOT included in the knot intervals!)
	xi = round(cnt_x[cnt_x.getMax()-1]);
	yi = round(cnt_y[cnt_y.getMax()-1]);
	dc->MoveTo(xi, IMAGE_Y-1-yi);
	dc->LineTo(old_x, IMAGE_Y-1-old_y); 
	//TRACE("# of cnt_pnts = %d\n", cnt_x.getMax());
	//TRACE("t = %.3f, [x, y] = [%.3f, %.3f] (I'm out!)\n", t, x, y);

	dc->SelectObject(pOldPen);
}

void MRBspline::DrawInterpCurve(CDC *dc, int r, int g, int b)
{
	double x, y;
	int xi, yi, i;
	int	old_x, old_y;
	double st_t, end_t;
	double N_val;
	int	pen_width;
	double KNOT_INC; 

	CPen pen;
	
	pen_width = 3;
	if (cnt_x.N < 4) { // can't generate cubic B-spline. Just display the polyline
		pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); // Thickened for PRESENCE paper! (from 2 to 3)
		CPen *pOldPen = (CPen *)dc->SelectObject(&pen);
		dc->MoveTo(round(cnt_x[0]), round(cnt_y[0]));
		for (i = 1; i < cnt_x.N; i++) {
			dc->LineTo(round(cnt_x[i]), round(cnt_y[i]));
		}
		dc->SelectObject(pOldPen);
		return;
	}
	pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	st_t = knot[0];
	end_t = knot[knot.getMax()-1];

	int polyline_count = GetPolylineCount();
	//KNOT_INC = (end_t - st_t) / 100.0;
	//KNOT_INC = (end_t - st_t) / (double)GetPolylineCount();
	//KNOT_INC = (end_t - st_t) / polyline_count / 5.0;
	KNOT_INC = (end_t - st_t) / polyline_count / 3.0;
	//KNOT_INC = (end_t - st_t) / polyline_count * 3.0;
	//TRACE("polyline_count = %d\n", GetPolylineCount());
	int polyline_count2 = (polyline_count > 100) ? 100: polyline_count;

	double tmp;
	/*
	int t_int = 0; // integer value of current t;
	arc_length = dist2(cnt_x[0], cnt_y[0], cnt_x[1], cnt_y[1]); 
	arc_length *= 4; // approximation to the arc length between successive data points
	KNOT_INC = 1/arc_length;		
	*/
		
	for (double t = st_t; t < end_t; t += KNOT_INC) { // Final end_t is NOT included in the knot intervals!!!
	//for (double t = st_t; t < end_t; t += 1.0) { // Final end_t is NOT included in the knot intervals!!!
		/*
		if ( ((int)t) > t_int ) {
			t_int = (int)t; // udpated
			// new segment: recompute KNOT_INC
			arc_length = dist2(cnt_x[(int)t], cnt_y[(int)t], cnt_x[(int)t+1], cnt_y[(int)t+1]); 
			KNOT_INC = 1 / arc_length;		
		}
		*/
		x = y = 0.0;
		//for (int i = 0; i < cnt_x.getMax(); i++) {
		for (int i = (int)t; i <= (int)t+3; i++) { // Only 4 control points are involved
			N_val = N(i, 4, t, knot);
			x += cnt_x[i] * N_val;	
			y += cnt_y[i] * N_val;	
		}
		xi = round(x);	
		yi = round(y);	// find the nearest integer values
		if (t > st_t) {
			dc->MoveTo(xi, yi);
			dc->LineTo(old_x, old_y);
			/////////////////////////////////////////////////
			// change the line width !!!
			switch (stroke_mode) {
				case UNIFORM_THICK:
					pen_width = 3;
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case MIDDLE_THICK:
					tmp = ((t-st_t) < (end_t-t)) ? (t-st_t) : (end_t-t);
					//pen_width = (int)(3 + 10 * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/10. * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/15. * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/7. * tmp / (end_t-st_t));
					//////////////////////////////////////////////////////////////////
					pen_width = (int)(1 + polyline_count2/5. * tmp / (end_t-st_t));
                    if (pen_width > 3) pen_width = 3;
					/////////////////////////////////////////////////////////////////
					//pen_width = (int)(3 + polyline_count2/3. * tmp / (end_t-st_t));
					////////////////////////////////////////////////////////////////
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case BEGIN_THICK:
					//pen_width = (int)(3 + 10*(t-st_t)/(end_t-st_t)); // from free point to seed
					pen_width = (int)(3 + polyline_count2/17.*(t-st_t)/(end_t-st_t)); // from free point to seed
					//pen_width = (int)(3 + polyline_count2/20.*(t-st_t)/(end_t-st_t)); // from free point to seed
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case END_THICK:
					//pen_width = (int)(3 + 10*(end_t-t)/(end_t-st_t)); // from free point to seed
					//pen_width = (int)(3 + polyline_count2/20.*(end_t-t)/(end_t-st_t)); // from free point to seed
					pen_width = (int)(3 + polyline_count2/17.*(end_t-t)/(end_t-st_t)); // from free point to seed
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
			}
		}
		old_x = xi, old_y = yi;
		//TRACE("t = %.3f, [x, y] = [%.3f, %.3f]\n", t, x, y);
		//TRACE("xi = %d, yi = %d\n", xi, yi);
		//break;
	}
	// Connect to the final control point!!! (because this point is NOT included in the knot intervals!)
	xi = round(cnt_x[cnt_x.getMax()-1]);
	yi = round(cnt_y[cnt_y.getMax()-1]);
	dc->MoveTo(xi, yi);
	dc->LineTo(old_x, old_y); 
	//TRACE("# of cnt_pnts = %d\n", cnt_x.getMax());
	//TRACE("t = %.3f, [x, y] = [%.3f, %.3f] (I'm out!)\n", t, x, y);

	dc->SelectObject(pOldPen);
}

void MRBspline::DrawBezierSegments(CDC& dc, int r, int g, int b, int pen_size)
// Draw an interpolating B-spline as a sequence of Bezier curve segments
// Note that we need to convert y coord to GDI mode!
{
	double x, y;
	int xi, yi, i, k, j;
	int	old_x, old_y;
	//double st_t, end_t;
	int	pen_width;
	double KNOT_INC; 

	CPen pen;
	
	//////////////////
	//pen_width = 3; // pen width
	pen_width = pen_size; // pen width
	///////////////////
	//if (cnt_x.getMax() < 4) { // can't generate cubic B-spline. Just display the polyline
	if (cnt_x.getMax() < 4) { // can't generate cubic B-spline. Just display the polyline
		pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); // Thickened for PRESENCE paper! (from 2 to 3)
		CPen *pOldPen = (CPen *)dc.SelectObject(&pen);
		dc.MoveTo(round(cnt_x[0]), IMAGE_Y-1-round(cnt_y[0]));
		if (cnt_x.getMax() == 1) { // only a single data point!
			dc.LineTo(round(cnt_x[0]), IMAGE_Y-1-round(cnt_y[0]));
		}
		else { // at least 2 data points!
			for (i = 1; i < cnt_x.getMax(); i++) {
				dc.LineTo(round(cnt_x[i]), IMAGE_Y-1-round(cnt_y[i]));
			}
		}
		dc.SelectObject(pOldPen);
		return;
	}
	pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
	CPen *pOldPen = (CPen *)dc.SelectObject(&pen);

	//st_t = knot[0];
	//end_t = knot[knot.getMax()-1];
    
	int polyline_count = GetPolylineCount();
	//KNOT_INC = (end_t - st_t) / 100.0;
	//KNOT_INC = (end_t - st_t) / (double)GetPolylineCount();
	//KNOT_INC = (end_t - st_t) / polyline_count / 5.0;
	//KNOT_INC = (end_t - st_t) / polyline_count / 3.0;
	//KNOT_INC = (end_t - st_t) / polyline_count * 3.0;
	//TRACE("polyline_count = %d\n", GetPolylineCount());
	int polyline_count2 = (polyline_count > 100) ? 100: polyline_count;

	matrix B(4, 4);
	B[0][0] = -1.0;	B[0][1] = 3.0;	B[0][2] = -3.0;	B[0][3] = 1.0;
	B[1][0] = 3.0;	B[1][1] = -6.0;	B[1][2] = 3.0;	B[1][3] = 0.0;
	B[2][0] = -3.0;	B[2][1] = 3.0;	B[2][2] = 0.0;	B[2][3] = 0.0;
	B[3][0] = 1.0;	B[3][1] = 0.0;	B[3][2] = 0.0;	B[3][3] = 0.0;
	matrix T(1, 4);
	matrix Px(4, 1);
	matrix Py(4, 1);
	
	matrix tmp1_x(4, 1); // B * P
	matrix tmp1_y(4, 1); // B * P
	matrix tmp2(1, 1); // T * (B * P)

	//double tmp;
	double arc_length;
	double factor = 0.05;
	int segment_num = (cnt_x.getMax() - 1)/3; // total number of Bezier segments
	
	//st_t = 0.0;
	//end_t = cnt_x.getMax();
		
	for (j = 0; j < segment_num; j++) {
		k = j * 3; // starting cp index for this segment
		arc_length = dist2(cnt_x[k], cnt_y[k], cnt_x[k+1], cnt_y[k+1]); 
		arc_length *= 4; // approximation of actual arc length of this segment
		KNOT_INC = 1 / arc_length;
		//KNOT_INC = 0.1;
		Px[0][0] = cnt_x[k]; Px[1][0] = cnt_x[k+1]; Px[2][0] = cnt_x[k+2]; Px[3][0] = cnt_x[k+3];
		Py[0][0] = IMAGE_Y-1-cnt_y[k]; Py[1][0] = IMAGE_Y-1-cnt_y[k+1]; Py[2][0] = IMAGE_Y-1-cnt_y[k+2]; Py[3][0] = IMAGE_Y-1-cnt_y[k+3];
		MatMatMult(B, Px, tmp1_x);
		MatMatMult(B, Py, tmp1_y);
		old_x = round(cnt_x[k]); old_y = round(IMAGE_Y-1-cnt_y[k]);
		for (double t = 0.0; t < 1.0; t += KNOT_INC) { // Final end_t is NOT included in the knot intervals!!!
			//x = y = 0.0;
			T[0][0] = t * t * t; T[0][1] = t * t; T[0][2] = t; T[0][3] = 1;
			MatMatMult(T, tmp1_x, tmp2);
			x = tmp2[0][0]; 
			//tmp2.print();
			MatMatMult(T, tmp1_y, tmp2);
			y = tmp2[0][0];
			//tmp2.print();
			///////////////////////////
			xi = round(x);	
			yi = round(y);	// find the nearest integer values
			if (t > 0.0) {
				dc.MoveTo(xi, yi);
				dc.LineTo(old_x, old_y);
				/////////////////////////////////////////////////
				// change the line width !!!
				switch (stroke_mode) {
					case UNIFORM_THICK:
						/*
						pen_width = 3;
						pen.DeleteObject();
						pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
						dc.SelectObject(&pen);
						*/
						break;
					case MIDDLE_THICK:
						/*
						tmp = ((t-st_t) < (end_t-t)) ? (t-st_t) : (end_t-t);
						//pen_width = (int)(3 + 10 * tmp / (end_t-st_t));
						pen_width = (int)(3 + polyline_count2/10. * tmp / (end_t-st_t));
						//pen_width = (int)(3 + polyline_count2/15. * tmp / (end_t-st_t));
						//pen_width = (int)(3 + polyline_count2/7. * tmp / (end_t-st_t));
						//////////////////////////////////////////////////////////////////
						//pen_width = (int)(1 + polyline_count2/5. * tmp / (end_t-st_t));
						if (pen_width > 3) pen_width = 3;
						/////////////////////////////////////////////////////////////////
						//pen_width = (int)(3 + polyline_count2/3. * tmp / (end_t-st_t));
						////////////////////////////////////////////////////////////////
						pen.DeleteObject();
						pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
						dc.SelectObject(&pen);
						*/
						break;
					case BEGIN_THICK:
						/*
						//pen_width = (int)(3 + 10*(t-st_t)/(end_t-st_t)); // from free point to seed
						pen_width = (int)(3 + polyline_count2/17.*(t-st_t)/(end_t-st_t)); // from free point to seed
						//pen_width = (int)(3 + polyline_count2/20.*(t-st_t)/(end_t-st_t)); // from free point to seed
						pen.DeleteObject();
						pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
						dc.SelectObject(&pen);
						*/
						break;
					case END_THICK:
						/*
						//pen_width = (int)(3 + 10*(end_t-t)/(end_t-st_t)); // from free point to seed
						//pen_width = (int)(3 + polyline_count2/20.*(end_t-t)/(end_t-st_t)); // from free point to seed
						pen_width = (int)(3 + polyline_count2/17.*(end_t-t)/(end_t-st_t)); // from free point to seed
						pen.DeleteObject();
						pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
						dc.SelectObject(&pen);
						*/
						break;
				}
			}
			////////////////////////////////
			old_x = xi, old_y = yi;
			//TRACE("t = %.3f, [x, y] = [%.3f, %.3f]\n", t, x, y);
			//TRACE("xi = %d, yi = %d\n", xi, yi);
			// Connect to the final control point!!! (because this point is NOT included in the knot intervals!)
			//break;
		}
		xi = round(cnt_x[k+3]);
		yi = round(IMAGE_Y-1-cnt_y[k+3]);
		dc.MoveTo(xi, yi);
		dc.LineTo(old_x, old_y); 
	}
	//TRACE("# of cnt_pnts = %d\n", cnt_x.getMax());
	//TRACE("t = %.3f, [x, y] = [%.3f, %.3f] (I'm out!)\n", t, x, y);

	dc.SelectObject(pOldPen);
}


void MRBspline::DrawHermiteCurve(CDC& dc, int r, int g, int b)
{
	double x, y;
	int xi, yi, k;
	int	old_x, old_y;
	double st_t, end_t;
	int	pen_width;
	double KNOT_INC; 

	CPen pen;
		
	if (cnt_x.N <= 1) { // can't generate Hermite-spline. Just display the point
		pen.CreatePen(PS_SOLID, 2, RGB(r, g, b)); // Thickened for PRESENCE paper! (from 2 to 3)
		CPen *pOldPen = (CPen *)dc.SelectObject(&pen);
		dc.MoveTo(round(cnt_x[0]), round(cnt_y[0]));
		dc.LineTo(round(cnt_x[0]), round(cnt_y[0]));
		dc.SelectObject(pOldPen);
		return;
	}
	pen_width = 3;
	pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
	CPen *pOldPen = (CPen *)dc.SelectObject(&pen);

	//st_t = knot[0];
	//end_t = knot[knot.getMax()-1];
    
	int polyline_count = GetPolylineCount();
	//KNOT_INC = (end_t - st_t) / 100.0;
	//KNOT_INC = (end_t - st_t) / (double)GetPolylineCount();
	//KNOT_INC = (end_t - st_t) / polyline_count / 5.0;
	//KNOT_INC = (end_t - st_t) / polyline_count / 3.0;
	//KNOT_INC = (end_t - st_t) / polyline_count * 3.0;
	//TRACE("polyline_count = %d\n", GetPolylineCount());
	int polyline_count2 = (polyline_count > 100) ? 100: polyline_count;

	matrix H(4, 4);
	H[0][0] = 2.0;	H[0][1] = -2.0;	H[0][2] = 1.0;	H[0][3] = -1.0;
	H[1][0] = -3.0; H[1][1] = 3.0;	H[1][2] = -2.0; H[1][3] = -1.0;
	H[2][0] = 0.0;	H[2][1] = 0.0;	H[2][2] = 1.0;	H[2][3] = 0.0;
	H[3][0] = 1.0;	H[3][1] = 0.0;	H[3][2] = 0.0;	H[3][3] = 0.0;
	matrix T(1, 4);
	matrix Px(4, 1);
	matrix Py(4, 1);
	
	matrix tmp1_x(4, 1); // H * P
	matrix tmp1_y(4, 1); // H * P
	matrix tmp2(1, 1); // T * (H * P)

	double tmp;
	double arc_length;
	double factor = 0.05;

	st_t = 0.0;
	end_t = cnt_x.getMax();
		
	for (k = 0; k < cnt_x.getMax()-1; k++) {
		arc_length = dist2(cnt_x[k], cnt_y[k], cnt_x[k+1], cnt_y[k+1]); // approximation!
		KNOT_INC = 1 / arc_length;
		//KNOT_INC = 0.1;
		//Px[0][0] = cnt_x[k]; Px[1][0] = cnt_x[k+1]; Px[2][0] = tan_x[k]/arc_length; Px[3][0] = tan_x[k+1]/arc_length;
		//Py[0][0] = cnt_y[k]; Py[1][0] = cnt_y[k+1]; Py[2][0] = tan_y[k]/arc_length; Py[3][0] = tan_y[k+1]/arc_length;
		Px[0][0] = cnt_x[k]; Px[1][0] = cnt_x[k+1]; Px[2][0] = tan_x[k]*arc_length*factor; Px[3][0] = tan_x[k+1]*arc_length*factor;
		Py[0][0] = cnt_y[k]; Py[1][0] = cnt_y[k+1]; Py[2][0] = tan_y[k]*arc_length*factor; Py[3][0] = tan_y[k+1]*arc_length*factor;
		MatMatMult(H, Px, tmp1_x);
		MatMatMult(H, Py, tmp1_y);
		old_x = round(cnt_x[k]); old_y = round(cnt_y[k]);
		for (double t = 0.0; t < 1.0; t += KNOT_INC) { // Final end_t is NOT included in the knot intervals!!!
			//x = y = 0.0;
			T[0][0] = t * t * t; T[0][1] = t * t; T[0][2] = t; T[0][3] = 1;
			MatMatMult(T, tmp1_x, tmp2);
			x = tmp2[0][0]; 
			//tmp2.print();
			MatMatMult(T, tmp1_y, tmp2);
			y = tmp2[0][0];
			//tmp2.print();
			///////////////////////////
			xi = round(x);	
			yi = round(y);	// find the nearest integer values
			if (t > 0.0) {
				dc.MoveTo(xi, yi);
				dc.LineTo(old_x, old_y);
				/////////////////////////////////////////////////
				// change the line width !!!
				switch (stroke_mode) {
					case UNIFORM_THICK:
						pen_width = 3;
						pen.DeleteObject();
						pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
						dc.SelectObject(&pen);
						break;
					case MIDDLE_THICK:
						tmp = ((t-st_t) < (end_t-t)) ? (t-st_t) : (end_t-t);
						//pen_width = (int)(3 + 10 * tmp / (end_t-st_t));
						pen_width = (int)(3 + polyline_count2/10. * tmp / (end_t-st_t));
						//pen_width = (int)(3 + polyline_count2/15. * tmp / (end_t-st_t));
						//pen_width = (int)(3 + polyline_count2/7. * tmp / (end_t-st_t));
						//////////////////////////////////////////////////////////////////
						//pen_width = (int)(1 + polyline_count2/5. * tmp / (end_t-st_t));
						if (pen_width > 3) pen_width = 3;
						/////////////////////////////////////////////////////////////////
						//pen_width = (int)(3 + polyline_count2/3. * tmp / (end_t-st_t));
						////////////////////////////////////////////////////////////////
						pen.DeleteObject();
						pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
						dc.SelectObject(&pen);
						break;
					case BEGIN_THICK:
						//pen_width = (int)(3 + 10*(t-st_t)/(end_t-st_t)); // from free point to seed
						pen_width = (int)(3 + polyline_count2/17.*(t-st_t)/(end_t-st_t)); // from free point to seed
						//pen_width = (int)(3 + polyline_count2/20.*(t-st_t)/(end_t-st_t)); // from free point to seed
						pen.DeleteObject();
						pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
						dc.SelectObject(&pen);
						break;
					case END_THICK:
						//pen_width = (int)(3 + 10*(end_t-t)/(end_t-st_t)); // from free point to seed
						//pen_width = (int)(3 + polyline_count2/20.*(end_t-t)/(end_t-st_t)); // from free point to seed
						pen_width = (int)(3 + polyline_count2/17.*(end_t-t)/(end_t-st_t)); // from free point to seed
						pen.DeleteObject();
						pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
						dc.SelectObject(&pen);
						break;
				}
			}
			////////////////////////////////
			old_x = xi, old_y = yi;
			//TRACE("t = %.3f, [x, y] = [%.3f, %.3f]\n", t, x, y);
			//TRACE("xi = %d, yi = %d\n", xi, yi);
			//break;
		}
		// Connect to the final control point!!! (because this point is NOT included in the knot intervals!)
		xi = round(cnt_x[k+1]);
		yi = round(cnt_y[k+1]);
		dc.MoveTo(xi, yi);
		dc.LineTo(old_x, old_y); 
	}
	//TRACE("# of cnt_pnts = %d\n", cnt_x.getMax());
	//TRACE("t = %.3f, [x, y] = [%.3f, %.3f] (I'm out!)\n", t, x, y);

	dc.SelectObject(pOldPen);
}


void MRBspline::DrawMRCurve(CDC *dc, int level, int r, int g, int b)
// Draw the specific level's curve using MRCntPnts[level]
{
	double x, y;
	int xi, yi;
	int	old_x, old_y;
	double st_t, end_t;
	double N_val;
	int	pen_width;
	double KNOT_INC; 

	vector cnt_x;
	vector cnt_y;

	
	cnt_x.copy(mr_cnt_x[level]);
	cnt_y.copy(mr_cnt_y[level]);

	//mr_cnt_x[level].print();
	//cnt_x.print();
	//mr_cnt_y[level].print();
	//cnt_y.print();

	////////////////////////////////////////////
	// Knot definition
	int M = (int)pow((double)2, level) + 3; // number of control points
	//N = (int)pow(2, 0) + 3; // control points: 4

	vector knot(M+4); // n + 1 + K
	
	knot[0] = knot[1] = knot[2] = knot[3] = 0.0;
	for (int i = 4; i < M; i++)
		knot[i] = (double)i - 3;
	knot[M] = knot[M+1] = knot[M+2] = knot[M+3] = (double)M - 3;
	/////////////////////////////////////////////////////////
    
	CPen pen;
		
	pen_width = 2;
	pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	st_t = knot[0];
	end_t = knot[knot.getMax()-1];

	//knot.print();

	int polyline_count = GetPolylineCount();
	//KNOT_INC = (end_t - st_t) / 100.0;
	//KNOT_INC = (end_t - st_t) / (double)GetPolylineCount();
	KNOT_INC = (end_t - st_t) / polyline_count / 5.0;
	//TRACE("polyline_count = %d\n", GetPolylineCount());
	int polyline_count2 = (polyline_count > 100) ? 100: polyline_count;

	//TRACE("KNOT_INC = %f\n", KNOT_INC);
	//TRACE("st_t = %f\n", st_t);
	//TRACE("end_t = %f\n", end_t);
	//TRACE("polyline_count = %d\n", polyline_count);
	//TRACE("polyline_count2 = %d\n", polyline_count2);

	double tmp;

	old_x = round(cnt_x[0]);
	old_y = round(cnt_y[0]);
		
	for (double t = st_t; t < end_t; t += KNOT_INC) { // Final end_t is NOT included in the knot intervals!!!
	//for (double t = st_t; t < end_t; t += 1.0) { // Final end_t is NOT included in the knot intervals!!!
		x = y = 0.0;
		//for (int i = 0; i < cnt_x.getMax(); i++) {
		for (int i = (int)t; i <= (int)t+3; i++) { // Only 4 control points are involved
			N_val = N(i, 4, t, knot);
			
			x += cnt_x[i] * N_val;	
			y += cnt_y[i] * N_val;	
			/*
			if (t > 0.19 && t < 0.195) {
				//TRACE("[%f] i = %d, N_val = %f\n", t, i, N_val);
				//TRACE("[%f] cnt_x[%d] = %f\, cnt_y[%d] = %f\n", t, i, cnt_x[i], i, cnt_y[i]);
				//TRACE("[%f] x = %f, y = %f\n", t, x, y);
			}
			*/
		}
		xi = round(x);	
		yi = round(y);	// find the nearest integer values
		//TRACE("[%f] xi = %d\, yi = %d\n", t, xi, yi);
		//TRACE("[%f] xi = %d\, yi = %d\n", t, xi, yi);
		if (t > st_t) {
			dc->MoveTo(xi, yi);
			dc->LineTo(old_x, old_y);
			/////////////////////////////////////////////////
			// change the line width !!!
			switch (stroke_mode) {
				case UNIFORM_THICK:
					pen_width = 4;
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case MIDDLE_THICK:
					tmp = ((t-st_t) < (end_t-t)) ? (t-st_t) : (end_t-t);
					//pen_width = (int)(3 + 10 * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/10. * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/15. * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/7. * tmp / (end_t-st_t));
					pen_width = (int)(3 + polyline_count2/5. * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/3. * tmp / (end_t-st_t));
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case BEGIN_THICK:
					//pen_width = (int)(3 + 10*(t-st_t)/(end_t-st_t)); // from free point to seed
					pen_width = (int)(3 + polyline_count2/17.*(t-st_t)/(end_t-st_t)); // from free point to seed
					//pen_width = (int)(3 + polyline_count2/20.*(t-st_t)/(end_t-st_t)); // from free point to seed
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case END_THICK:
					//pen_width = (int)(3 + 10*(end_t-t)/(end_t-st_t)); // from free point to seed
					//pen_width = (int)(3 + polyline_count2/20.*(end_t-t)/(end_t-st_t)); // from free point to seed
					pen_width = (int)(3 + polyline_count2/17.*(end_t-t)/(end_t-st_t)); // from free point to seed
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
			}
		}
		old_x = xi, old_y = yi;
		//TRACE("t = %.3f, [x, y] = [%.3f, %.3f]\n", t, x, y);
		//TRACE("xi = %d, yi = %d\n", xi, yi);
		//break;
	}
	// Connect to the final control point!!! (because this point is NOT included in the knot intervals!)
	xi = round(cnt_x[cnt_x.getMax()-1]);
	yi = round(cnt_y[cnt_y.getMax()-1]);
	dc->MoveTo(xi, yi);
	dc->LineTo(old_x, old_y); 
	//TRACE("# of cnt_pnts = %d\n", cnt_x.getMax());
	//TRACE("t = %.3f, [x, y] = [%.3f, %.3f] (I'm out!)\n", t, x, y);

	dc->SelectObject(pOldPen);
}

void MRBspline::DrawBlendedCurve(CDC *dc, int level1, int level2, double u, int r, int g, int b)
// Draw the blended curve between two level curves
{
	double x1, y1, x2, y2;
	int xi, yi;
	int	old_x, old_y;
	double st_t, end_t, st_t2, end_t2, t, t2;
	double N_val;
	int	pen_width;
	double KNOT_INC; 
	int M;

	vector cnt_x1;
	vector cnt_y1;
	vector cnt_x2;
	vector cnt_y2;

	cnt_x1.copy(mr_cnt_x[level1]);
	cnt_y1.copy(mr_cnt_y[level1]);
	cnt_x2.copy(mr_cnt_x[level2]);
	cnt_y2.copy(mr_cnt_y[level2]);

	//mr_cnt_x[level].print();
	//cnt_x.print();
	//mr_cnt_y[level].print();
	//cnt_y.print();

	////////////////////////////////////////////
	// Knot1 definition
	M = (int)pow((double)2, level1) + 3; // number of control points
	//N = (int)pow(2, 0) + 3; // control points: 4

	vector knot1(M+4); // n + 1 + K
	
	knot1[0] = knot1[1] = knot1[2] = knot1[3] = 0.0;
	for (int i = 4; i < M; i++)
		knot1[i] = (double)i - 3;
	knot1[M] = knot1[M+1] = knot1[M+2] = knot1[M+3] = (double)M - 3;
	/////////////////////////////////////////////////////////

	////////////////////////////////////////////
	// Knot2 definition
	M = (int)pow((double)2, level2) + 3; // number of control points
	//N = (int)pow(2, 0) + 3; // control points: 4

	vector knot2(M+4); // n + 1 + K
	
	knot2[0] = knot2[1] = knot2[2] = knot2[3] = 0.0;
	for (int i = 4; i < M; i++)
		knot2[i] = (double)i - 3;
	knot2[M] = knot2[M+1] = knot2[M+2] = knot2[M+3] = (double)M - 3;
	/////////////////////////////////////////////////////////
    
	CPen pen;
		
	pen_width = 2;
	pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	st_t = knot1[0];
	end_t = knot1[knot1.getMax()-1];
	st_t2 = knot2[0];
	end_t2 = knot2[knot2.getMax()-1];

	//knot.print();

	int polyline_count = GetPolylineCount();
	//KNOT_INC = (end_t - st_t) / 100.0;
	//KNOT_INC = (end_t - st_t) / (double)GetPolylineCount();
	KNOT_INC = (end_t - st_t) / polyline_count / 5.0;
	//TRACE("polyline_count = %d\n", GetPolylineCount());
	int polyline_count2 = (polyline_count > 100) ? 100: polyline_count;

	//TRACE("KNOT_INC = %f\n", KNOT_INC);
	//TRACE("st_t = %f\n", st_t);
	//TRACE("end_t = %f\n", end_t);
	//TRACE("polyline_count = %d\n", polyline_count);
	//TRACE("polyline_count2 = %d\n", polyline_count2);

	double tmp;

	old_x = round(cnt_x1[0]);
	old_y = round(cnt_y1[0]);
		
	for (t = st_t; t < end_t; t += KNOT_INC) { // Final end_t is NOT included in the knot intervals!!!
	//for (double t = st_t; t < end_t; t += 1.0) { // Final end_t is NOT included in the knot intervals!!!
		x1 = y1 = 0.0;
		for (int i = (int)t; i <= (int)t+3; i++) { // Only 4 control points are involved
			N_val = N(i, 4, t, knot1);
			
			x1 += cnt_x1[i] * N_val;	
			y1 += cnt_y1[i] * N_val;	
		}
		//////////////////////////////////////////
        t2 = st_t2 + end_t2 * ( t / end_t );
		x2 = y2 = 0.0;
		for (int i = (int)t2; i <= (int)t2+3; i++) { // Only 4 control points are involved
			N_val = N(i, 4, t2, knot2);
			
			x2 += cnt_x2[i] * N_val;	
			y2 += cnt_y2[i] * N_val;	
		}
		////////////////////////////////////////////
		xi = round( (1-u) * x1 + u * x2 );	
		yi = round( (1-u) * y1 + u * y2 );	// find the nearest integer values
		//TRACE("[%f] xi = %d\, yi = %d\n", t, xi, yi);
		//TRACE("[%f] xi = %d\, yi = %d\n", t, xi, yi);
		if (t > st_t) {
			dc->MoveTo(xi, yi);
			dc->LineTo(old_x, old_y);
			/////////////////////////////////////////////////
			// change the line width !!!
			switch (stroke_mode) {
				case UNIFORM_THICK:
					pen_width = 4;
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case MIDDLE_THICK:
					tmp = ((t-st_t) < (end_t-t)) ? (t-st_t) : (end_t-t);
					//pen_width = (int)(3 + 10 * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/10. * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/15. * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/7. * tmp / (end_t-st_t));
					pen_width = (int)(3 + polyline_count2/5. * tmp / (end_t-st_t));
					//pen_width = (int)(3 + polyline_count2/3. * tmp / (end_t-st_t));
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case BEGIN_THICK:
					//pen_width = (int)(3 + 10*(t-st_t)/(end_t-st_t)); // from free point to seed
					pen_width = (int)(3 + polyline_count2/17.*(t-st_t)/(end_t-st_t)); // from free point to seed
					//pen_width = (int)(3 + polyline_count2/20.*(t-st_t)/(end_t-st_t)); // from free point to seed
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
				case END_THICK:
					//pen_width = (int)(3 + 10*(end_t-t)/(end_t-st_t)); // from free point to seed
					//pen_width = (int)(3 + polyline_count2/20.*(end_t-t)/(end_t-st_t)); // from free point to seed
					pen_width = (int)(3 + polyline_count2/17.*(end_t-t)/(end_t-st_t)); // from free point to seed
					pen.DeleteObject();
					pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
					dc->SelectObject(&pen);
					break;
			}
		}
		old_x = xi, old_y = yi;
		//TRACE("t = %.3f, [x, y] = [%.3f, %.3f]\n", t, x, y);
		//TRACE("xi = %d, yi = %d\n", xi, yi);
		//break;
	}
	// Connect to the final control point!!! (because this point is NOT included in the knot intervals!)
	xi = round(cnt_x1[cnt_x1.getMax()-1]);
	yi = round(cnt_y1[cnt_y1.getMax()-1]);
	dc->MoveTo(xi, yi);
	dc->LineTo(old_x, old_y); 
	//TRACE("# of cnt_pnts = %d\n", cnt_x1.getMax());
	//TRACE("t = %.3f, [x, y] = [%.3f, %.3f] (I'm out!)\n", t, x, y);

	dc->SelectObject(pOldPen);
}


void MRBspline::CopyLevel(int k, MRBspline& s)
// Copy level i control points and knots to target MRBspline
{
	int	i, j, N;

	j = s.GetMaxLevel(); // Maximum level of the source curve
	
	if (k < 0) k = 0;
	else if (k > j) k = j;
	
	//N = (int)pow(2, j) + 3;

	//old_cnt_x.copy(curve.cnt_x);

	vector cnt_x;
	vector cnt_y;
	vector t_cnt_x;
	vector t_cnt_y;
	vector D_x;
	vector D_y;
	vector t_D_x;
	vector t_D_y;

	cnt_x.copy(s.GetCnt0Pnts_x());
	cnt_y.copy(s.GetCnt0Pnts_y());

	N = (int)pow((double)2, 0) + 3; // Assume the level is 0

	for (i = 0; i < k; i++) {
		MatVecMult(P[i+1], cnt_x, t_cnt_x);
		MatVecMult(P[i+1], cnt_y, t_cnt_y);
		D_x.copy(s.GetDiffPnts_x(i));
		D_y.copy(s.GetDiffPnts_y(i));
		MatVecMult(Q[i+1], D_x, t_D_x);
		MatVecMult(Q[i+1], D_y, t_D_y);
		VecVecAdd(t_cnt_x, t_D_x, cnt_x);
		VecVecAdd(t_cnt_y, t_D_y, cnt_y);
		N = (int)pow((double)2, i+1) + 3;
	}
	
	//else return; // do not need to draw anything
	
	vector knot(N+4); // n + 1 + K

	// Knot definition
	knot[0] = knot[1] = knot[2] = knot[3] = 0.0;
	for (i = 4; i < N; i++)
		knot[i] = (double)i - 3;
	knot[N] = knot[N+1] = knot[N+2] = knot[N+3] = (double)N - 3;
			
	SetKnotVector(knot);
	SetCntPnts(cnt_x, cnt_y);
	//curve.SetDiffPnts(j, D_x, D_y);
	//DrawCurve(dc, r, g, b);
	//DrawCntPnts(dc, 0, 0, 0);
}

void MRBspline::DrawCntPnts(CDC *dc, int r, int g, int b)
{
	int x, y;
	int	old_x, old_y;
	int	pen_width;
	int	i;

	CPen pen, pen2;
		
	//pen_width = 3;
	pen_width = 1;
	pen.CreatePen(PS_SOLID, pen_width, RGB(0, 0, 255)); 
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	old_x = round(cnt_x[0]);
	old_y = round(cnt_y[0]);
	for (i = 0; i < cnt_x.getMax(); i++) {
		x = round(cnt_x[i]);
		y = round(cnt_y[i]);
		dc->MoveTo(x, y);
		dc->LineTo(old_x, old_y); 
		old_x = x;
		old_y = y;
	}
	dc->SelectObject(pOldPen);

	pen_width = 5;
	pen2.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
	pOldPen = (CPen *)dc->SelectObject(&pen2);

	for (i = 0; i < cnt_x.getMax(); i++) {
		x = round(cnt_x[i]);
		y = round(cnt_y[i]);
		dc->MoveTo(x, y);
		dc->LineTo(x, y); 
	}

	dc->SelectObject(pOldPen);
}

void MRBspline::DrawMRCntPnts(CDC *dc, int level, int r, int g, int b)
{
	int x, y;
	int	old_x, old_y;
	int	pen_width;
	int	i;
	vector cnt_x, cnt_y;

	CPen pen, pen2;
		
	//pen_width = 3;
	pen_width = 1;
	pen.CreatePen(PS_SOLID, pen_width, RGB(0, 0, 255)); 
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	cnt_x.copy(mr_cnt_x[level]);
	cnt_y.copy(mr_cnt_y[level]);

	old_x = round(cnt_x[0]);
	old_y = round(cnt_y[0]);
	for (i = 0; i < cnt_x.getMax(); i++) {
		x = round(cnt_x[i]);
		y = round(cnt_y[i]);
		dc->MoveTo(x, y);
		dc->LineTo(old_x, old_y); 
		old_x = x;
		old_y = y;
	}
	dc->SelectObject(pOldPen);

	pen_width = 5;
	pen2.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
	pOldPen = (CPen *)dc->SelectObject(&pen2);

	for (i = 0; i < cnt_x.getMax(); i++) {
		x = round(cnt_x[i]);
		y = round(cnt_y[i]);
		dc->MoveTo(x, y);
		dc->LineTo(x, y); 
	}

	dc->SelectObject(pOldPen);
}


void MRBspline::DrawCntPnt(CDC *dc, int i, int r, int g, int b)
{
	int x, y;
	//int	old_x, old_y;
	//double st_t, end_t;
	//double N_val;
	int	pen_width;
	//int	i;

	CPen pen;
		
	//pen_width = 3;
	pen_width = 5;
	pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	x = round(cnt_x[i]);
	y = round(cnt_y[i]);
	dc->MoveTo(x, y);
	dc->LineTo(x, y); 

	dc->SelectObject(pOldPen);
}

void MRBspline::DrawMRCntPnt(CDC *dc, int level, int i, int r, int g, int b)
{
	int x, y;
	//int	old_x, old_y;
	//double st_t, end_t;
	//double N_val;
	int	pen_width;
	//int	i;
	vector cnt_x, cnt_y;

	cnt_x.copy(mr_cnt_x[level]);
	cnt_y.copy(mr_cnt_y[level]);

	CPen pen;
		
	//pen_width = 3;
	pen_width = 5;
	pen.CreatePen(PS_SOLID, pen_width, RGB(r, g, b)); 
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	x = round(cnt_x[i]);
	y = round(cnt_y[i]);
	dc->MoveTo(x, y);
	dc->LineTo(x, y); 

	dc->SelectObject(pOldPen);
}


//#define dist2(x1, y1, x2, y2) sqrt( (((double)x1)-((double)x2))*(((double)x1)-((double)x2)) + (((double)y1)-((double)y2))*(((double)y1)-((double)y2)) )

int MRBspline::FindCntPnt(int x, int y)
{
	int mask = 50;
	int half;
	int i;
	int min_index = -1;
	double tmp_dist;
	double min_dist = 10000000;

	half = mask / 2;

	for (i = 0; i < cnt_x.getMax(); i++) {
		if (cnt_x[i] >= x - half && cnt_x[i] <= x + half 
			&& cnt_y[i] >= y - half && cnt_y[i] <= y + half) {
				tmp_dist = dist2(x, y, cnt_x[i], cnt_y[i]);
				if ( tmp_dist < min_dist ) {
					min_index = i;
					min_dist = tmp_dist;
				}
		}
	}

	return min_index;

}

int MRBspline::FindMRCntPnt(int level, int x, int y)
{
	int mask = 50;
	int half;
	int i;
	int min_index = -1;
	double tmp_dist;
	double min_dist = 10000000;
	vector cnt_x, cnt_y;

	half = mask / 2;

	cnt_x.copy(mr_cnt_x[level]);
	cnt_y.copy(mr_cnt_y[level]);

	for (i = 0; i < cnt_x.getMax(); i++) {
		if (cnt_x[i] >= x - half && cnt_x[i] <= x + half 
			&& cnt_y[i] >= y - half && cnt_y[i] <= y + half) {
				tmp_dist = dist2(x, y, cnt_x[i], cnt_y[i]);
				if ( tmp_dist < min_dist ) {
					min_index = i;
					min_dist = tmp_dist;
				}
		}
	}

	return min_index;

}

void MRBspline::SetKnotVector(vector& in)
{
	knot.copy(in);
}

void MRBspline::SetCntPnts(vector& in_x, vector& in_y)
{
	cnt_x.copy(in_x);
	cnt_y.copy(in_y);
}

void MRBspline::SetCntPnts() 
// Set normal B-spline Control Points from the input data points
{
	int i;
	int N = pnts.size(); // number of control points (data points)

	vector in_x(N);
	vector in_y(N);

	for (i = 0; i < N; i++) {
		in_x[i] = pnts[i].x;
		in_y[i] = pnts[i].y;
	}

	vector tmp_knot(N+4); // n + 1 + K

	// Knot definition
	tmp_knot[0] = tmp_knot[1] = tmp_knot[2] = tmp_knot[3] = 0.0;
	for (int i = 4; i < N; i++)
		tmp_knot[i] = (double)i - 3;
	tmp_knot[N] = tmp_knot[N+1] = tmp_knot[N+2] = tmp_knot[N+3] = (double)N - 3;
			
	SetKnotVector(tmp_knot);
	cnt_x.copy(in_x);
	cnt_y.copy(in_y);
}

void MRBspline::SetCntPnts(deque<double2D>& pnts)
{
	int i, N;
	
	N = pnts.size();

	cnt_x.init(N);
	cnt_y.init(N);
	//cnt_x.copy(in_x);
	//cnt_y.copy(in_y);
	for (i = 0; i < N; i++) {
		cnt_x[i] = pnts[i].x;
		cnt_y[i] = pnts[i].y;
	}
}

void MRBspline::SetDataPnts(deque<double2D>& input_pnts)
{
	int i, N;
	
	N = input_pnts.size();
	
	pnts.clear();

	for (i = 0; i < N; i++) {
		pnts.push_back(input_pnts[i]);
	}
}

void MRBspline::SetSamplePnts(deque<double2D>& pnts)
{
	int i, N;
	
	N = pnts.size();
	
	sample_pnts.clear();

	for (i = 0; i < N; i++) {
		//sample_pnts.insert(sample_pnts.end(), pnts[i]);
		sample_pnts.push_back(pnts[i]);
	}
}

void MRBspline::SetHeadTailDir() 
{
	int count, i;
	
	h[0] = h[1] = 0.0;
	count = 0;
	for (i = 0; i < 3; i++) {
		if (i+1 > (signed)pnts.size()-1) break;
		h[0] += pnts[i].x - pnts[i+1].x;
		h[1] += pnts[i].y - pnts[i+1].y;
		count++;
	}
	h[0] /= (double)count;
	h[1] /= (double)count;
	h.make_unit();

	t[0] = t[1] = 0.0;
	count = 0;
	//TRACE("curve->pnts.size() = %d\n", curve->pnts.size());
	for (i = 0; i < 3; i++) {
		if ((signed)pnts.size()-2-i < 0) break;
		//TRACE("i = %d\n", i);
		//TRACE("curve->pnts.size()-2-%d = %d\n", i, curve->pnts.size()-2-i);
		t[0] += pnts[pnts.size()-1-i].x - pnts[pnts.size()-2-i].x;
		t[1] += pnts[pnts.size()-1-i].y - pnts[pnts.size()-2-i].y;
		count++;
	}
	t[0] /= (double)count;
	t[1] /= (double)count;
	t.make_unit();
}

void MRBspline::SetHeadDir() 
{
	int count, i;
	
	h[0] = h[1] = 0.0;
	count = 0;
	for (i = 0; i < 3; i++) {
		if (i+1 > (signed)pnts.size()-1) break;
		h[0] += pnts[i].x - pnts[i+1].x;
		h[1] += pnts[i].y - pnts[i+1].y;
		count++;
	}
	h[0] /= (double)count;
	h[1] /= (double)count;
	h.make_unit();

	
}

void MRBspline::SetTailDir() 
{
	int count, i;
	
	t[0] = t[1] = 0.0;
	count = 0;
	//TRACE("curve->pnts.size() = %d\n", curve->pnts.size());
	for (i = 0; i < 3; i++) {
		if ((signed)pnts.size()-2-i < 0) break;
		//TRACE("i = %d\n", i);
		//TRACE("curve->pnts.size()-2-%d = %d\n", i, curve->pnts.size()-2-i);
		t[0] += pnts[pnts.size()-1-i].x - pnts[pnts.size()-2-i].x;
		t[1] += pnts[pnts.size()-1-i].y - pnts[pnts.size()-2-i].y;
		count++;
	}
	t[0] /= (double)count;
	t[1] /= (double)count;
	t.make_unit();
}

//void MRBspline::SetInterpCntPnts(deque<double2D>& pnts)
void MRBspline::SetInterpCntPnts()
{
	int i, N, M, k;
	
	N = pnts.size();
	if (N < 3) { // we need at least 3 data points to compute the control points!
		cnt_x.init(N); // use data points as control points
		cnt_y.init(N);
		for (i = 0; i < N; i++) {
			cnt_x[i] = pnts[i].x;
			cnt_y[i] = pnts[i].y;
		}
		return; 
	}

	//vector knot;
	vector Ax(N), Ay(N), B(N), dx(N), dy(N);
	//Ax.zero(); Ay.zero(); B.zero();
	dx[0] = (pnts[1].x - pnts[0].x) / 3;
	dy[0] = (pnts[1].y - pnts[0].y) / 3;
	dx[N-1] = (pnts[N-1].x - pnts[N-2].x) / 3;
	dy[N-1] = (pnts[N-1].y - pnts[N-2].y) / 3;

	B[1] = -0.25;
	Ax[1] = (pnts[2].x - pnts[0].x - dx[0]) / 4;
	Ay[1] = (pnts[2].y - pnts[0].y - dy[0]) / 4;
	for (i = 2; i < N-1; i++) {
		B[i] = -1/(4 + B[i-1]);
		Ax[i] = -(pnts[i+1].x - pnts[i-1].x - Ax[i-1])*B[i];
		Ay[i] = -(pnts[i+1].y - pnts[i-1].y - Ay[i-1])*B[i];
	}
	for (i = N-2; i > 0; i--) {
		dx[i] = Ax[i] + dx[i+1] * B[i];
		dy[i] = Ay[i] + dy[i+1] * B[i];
	}
	////////////////////////////////////////////////
	// M: the number of control points
	// N: number of data points
	// N-1: number of short curve segments between data points
	M = N + (N-1) * 2; 
	cnt_x.init(M);
	cnt_y.init(M);
	
	k = 0; // control point index
	cnt_x[0] = pnts[0].x; cnt_y[0] = pnts[0].y; // starting point
	for (i = 0; i < N-1; i++) { // for each segment between data points
		k++;
		cnt_x[k] = pnts[i].x + dx[i]; 
		cnt_y[k] = pnts[i].y + dy[i];
		k++;
		cnt_x[k] = pnts[i+1].x - dx[i+1]; 
		cnt_y[k] = pnts[i+1].y - dy[i+1];
		k++;
		cnt_x[k] = pnts[i+1].x; cnt_y[k] = pnts[i+1].y;
	}
	////////////////////////////////////////////////////////////////
	/////////////////////////////////
	/*
	knot.init(M+4); // n + 1 + K
	// Knot definition
	knot[0] = knot[1] = knot[2] = knot[3] = 0.0;
	for (i = 4; i < M; i++)
		knot[i] = (double)i - 3;
	knot[M] = knot[M+1] = knot[M+2] = knot[M+3] = (double)M - 3;
    			
	SetKnotVector(knot);
	*/
}


void MRBspline::SetTangents(deque<double2D>& pnts)
{
	int i, N;
	double dist;
	
	N = pnts.size();
	//N = cnt_x.getMax();
	if (N < 2) return; // just a point!

	tan_x.init(N);
	tan_y.init(N);
	//cnt_x.copy(in_x);
	//cnt_y.copy(in_y);
	///////////////////////////////
	// Set the unit tangent vector: between [0, 1]
	dist = dist2(pnts[1].x, pnts[1].y, pnts[0].x, pnts[0].y);
    tan_x[0] = (pnts[1].x - pnts[0].x) / dist;
	tan_y[0] = (pnts[1].y - pnts[0].y) / dist;
	dist = dist2(pnts[N-1].x, pnts[N-1].y, pnts[N-2].x, pnts[N-2].y);
	tan_x[N-1] = (pnts[N-1].x - pnts[N-2].x) / dist;
	tan_y[N-1] = (pnts[N-1].y - pnts[N-2].y) / dist;
	if (N >= 3) { // we have at least 3 data points
		for (i = 1; i <= N-2; i++) {
			dist = dist2(pnts[i+1].x, pnts[i+1].y, pnts[i-1].x, pnts[i-1].y);
			tan_x[i] = (pnts[i+1].x - pnts[i-1].x) / dist;
			tan_y[i] = (pnts[i+1].y - pnts[i-1].y) / dist;
		}
	}
	////////////////////////////////////////////////////////////////
	/*
	/////////////////////////////////////
	// follow Catmull-Rom convention
	tan_x[0] = (pnts[1].x - pnts[0].x) / 3.0;
	tan_y[0] = (pnts[1].y - pnts[0].y) / 3.0;
	tan_x[N-1] = (pnts[N-1].x - pnts[N-2].x) / 3.0;
	tan_y[N-1] = (pnts[N-1].y - pnts[N-2].y) / 3.0;
	if (N >= 3) { // we have at least 3 data points
		for (i = 1; i <= N-2; i++) {
			tan_x[i] = (pnts[i+1].x - pnts[i-1].x) / 2.0;
			tan_y[i] = (pnts[i+1].y - pnts[i-1].y) / 2.0;
		}
	}
	*/
}

void MRBspline::SetTangents(deque<double2D>& pnts, Field& gfield)
// Follow the tangent vector computed in gfield
{
	int i, N;
	int xi, yi;
	
	N = pnts.size();
	//N = cnt_x.getMax();
	if (N < 2) return; // just a point!

	tan_x.init(N);
	tan_y.init(N);
	//cnt_x.copy(in_x);
	//cnt_y.copy(in_y);
	///////////////////////////////
	// Set the unit tangent vector: between [0, 1]
	for (i = 0; i < N; i++) {
		xi = round(pnts[i].x);
		yi = round(pnts[i].y);
		tan_x[i] = -gfield[xi][yi].gy;
		tan_y[i] = gfield[xi][yi].gx;
		//tan_x[i] = gfield[xi][yi].gx * 50;
		//tan_y[i] = gfield[xi][yi].gy * 50;
	}
	////////////////////////////////////////////////////////////////
	/*
	/////////////////////////////////////
	// follow Catmull-Rom convention
	tan_x[0] = (pnts[1].x - pnts[0].x) / 3.0;
	tan_y[0] = (pnts[1].y - pnts[0].y) / 3.0;
	tan_x[N-1] = (pnts[N-1].x - pnts[N-2].x) / 3.0;
	tan_y[N-1] = (pnts[N-1].y - pnts[N-2].y) / 3.0;
	if (N >= 3) { // we have at least 3 data points
		for (i = 1; i <= N-2; i++) {
			tan_x[i] = (pnts[i+1].x - pnts[i-1].x) / 2.0;
			tan_y[i] = (pnts[i+1].y - pnts[i-1].y) / 2.0;
		}
	}
	*/
}


void MRBspline::SetMRCntPnts(int level, vector& in_x, vector& in_y)
{
	mr_cnt_x[level].copy(in_x);
	mr_cnt_y[level].copy(in_y);
}

void MRBspline::SetCnt0Pnts(vector& in_x, vector& in_y)
{
	C0_x.copy(in_x);
	C0_y.copy(in_y);
}

void MRBspline::SetDiffPnts(int j, vector* in_x, vector* in_y)
{
	D_x = new vector[j];
	D_y = new vector[j];

	for (int i = 0; i < j; i++) {
		D_x[i].copy(in_x[i]);
		D_y[i].copy(in_y[i]);
	}
}

void MRBspline::SetEditedDiffPnts(int j, vector* in_x, vector* in_y)
{
	//D_x = new vector[j];
	//D_y = new vector[j];

	for (int i = 0; i < j; i++) {
		D_x[i].copy(in_x[i]);
		D_y[i].copy(in_y[i]);
	}
}

void MRBspline::AddDataPointRear(double x, double y) 
{ 
	double2D p; 
	p.x = x; p.y = y; 
	pnts.insert(pnts.end(), p); // add at the end
}

void MRBspline::AddDataPointHead(double x, double y) 
{ 
	double2D p; 
	p.x = x; p.y = y; 
	pnts.insert(pnts.begin(), p); // add at the head
}

/*
#define NRANSI
#include "nrutil.h"
#define TINY 1.0e-20;

void ludcmp(Mat a, int n, int *indx, double *d)
{
  int i, imax, j, k;
  double big, dum, sum, temp;
  double *vv;

  vv = vector(1, n);
  *d = 1.0;
  for (i = 1; i <= n; i++) {
    big = 0.0;
    for (j = 1; j <= n; j++)
      if ((temp = fabs(a[i-1][j-1])) > big)
        big = temp;
    if (big == 0.0)
      nrerror("Singular matrix in routine ludcmp");
    vv[i] = 1.0 / big;
  }
  for (j = 1; j <= n; j++) {
    for (i = 1; i < j; i++) {
      sum = a[i-1][j-1];
      for (k = 1; k < i; k++)
        sum -= a[i-1][k-1] * a[k-1][j-1];
      a[i-1][j-1] = sum;
    }
    big = 0.0;
    for (i = j; i <= n; i++) {
      sum = a[i-1][j-1];
      for (k = 1; k < j; k++)
        sum -= a[i-1][k-1] * a[k-1][j-1];
      a[i-1][j-1] = sum;
      if ((dum = vv[i] * fabs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 1; k <= n; k++) {
        dum = a[imax-1][k-1];
        a[imax-1][k-1] = a[j-1][k-1];
        a[j-1][k-1] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax; // here is a problem 
    if (a[j-1][j-1] == 0.0)
      a[j-1][j-1] = TINY;
    if (j != n) {
      dum = 1.0 / (a[j-1][j-1]);
      for (i = j + 1; i <= n; i++)
        a[i-1][j-1] *= dum;
    }
  }
  free_vector(vv, 1, n);
}
#undef TINY
#undef NRANSI

void lubksb(Mat a, int n, int *indx, double b[])
{
  int i, ii = -1, ip, j;
  double sum;

  for (i = 1; i <= n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii > -1)
      for (j = ii; j <= i - 1; j++)
        sum -= a[i-1][j-1] * b[j];
    else if (sum)
      ii = i;
    b[i] = sum;
  }
  for (i = n; i >= 1; i--) {
    sum = b[i];
    for (j = i + 1; j <= n; j++)
      sum -= a[i-1][j-1] * b[j];
    b[i] = sum / a[i-1][i-1];
  }
}

void inverse_mat(Mat a, int n, Mat inverse)
{
	int	indx[MAX_CNT_PNTS];
	double		d, col[MAX_CNT_PNTS];
	int	i, j;

	ludcmp(a, n, indx, &d);
	for (j = 1; j <= n; j++) {
		for (i = 1; i <= n; i++) col[i] = 0.0;
		col[j] = 1.0;
		lubksb(a, n, indx, col);
		for (i = 1; i <= n; i++) inverse[i-1][j-1] = col[i];
	}
	
}


*/