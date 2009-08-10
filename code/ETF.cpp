#include "stdafx.h"

#include <cmath>
#include "globals.h"
#include "ETF.h"

void ETF::set(imatrix& image) 
{
	int i, j;
	double MAX_VAL = 1020.; // 255 + 2 * 255 + 255 (Maximum possible Sobel value)
	double v[2];

	max_grad = -1.;

	for (i = 1; i < Nr - 1; i++) { 
		for (j = 1; j < Nc - 1; j++) {
			////////////////////////////////////////////////////////////////
			p[i][j].tx = (image[i+1][j-1] + 2*(double)image[i+1][j] + image[i+1][j+1] 
				- image[i-1][j-1] - 2*(double)image[i-1][j] - image[i-1][j+1]) / MAX_VAL;
			p[i][j].ty = (image[i-1][j+1] + 2*(double)image[i][j+1] + image[i+1][j+1]
				- image[i-1][j-1] - 2*(double)image[i][j-1] - image[i+1][j-1]) / MAX_VAL;
			/////////////////////////////////////////////
			v[0] = p[i][j].tx;
			v[1] = p[i][j].ty;
			p[i][j].tx = -v[1];
			p[i][j].ty = v[0];
			//////////////////////////////////////////////
			p[i][j].mag = sqrt(p[i][j].tx * p[i][j].tx + p[i][j].ty * p[i][j].ty);

			if (p[i][j].mag > max_grad) {
				max_grad = p[i][j].mag;
			}
		}
	}

	for (i = 1; i <= Nr - 2; i++) {
		p[i][0].tx = p[i][1].tx;
		p[i][0].ty = p[i][1].ty;
		p[i][0].mag = p[i][1].mag;
		p[i][Nc - 1].tx = p[i][Nc - 2].tx;
		p[i][Nc - 1].ty = p[i][Nc - 2].ty;
		p[i][Nc - 1].mag = p[i][Nc - 2].mag;
	}
	
	for (j = 1; j <= Nc - 2; j++) {
		p[0][j].tx = p[1][j].tx;
		p[0][j].ty = p[1][j].ty;
		p[0][j].mag = p[1][j].mag;
		p[Nr - 1][j].tx = p[Nr - 2][j].tx;
		p[Nr - 1][j].ty = p[Nr - 2][j].ty;
		p[Nr - 1][j].mag = p[Nr - 2][j].mag;
	}
	
	p[0][0].tx = ( p[0][1].tx + p[1][0].tx ) / 2;
	p[0][0].ty = ( p[0][1].ty + p[1][0].ty ) / 2;
	p[0][0].mag = ( p[0][1].mag + p[1][0].mag ) / 2;
	p[0][Nc-1].tx = ( p[0][Nc-2].tx + p[1][Nc-1].tx ) / 2;
	p[0][Nc-1].ty = ( p[0][Nc-2].ty + p[1][Nc-1].ty ) / 2;
	p[0][Nc-1].mag = ( p[0][Nc-2].mag + p[1][Nc-1].mag ) / 2;
	p[Nr-1][0].tx = ( p[Nr-1][1].tx + p[Nr-2][0].tx ) / 2;
	p[Nr-1][0].ty = ( p[Nr-1][1].ty + p[Nr-2][0].ty ) / 2;
	p[Nr-1][0].mag = ( p[Nr-1][1].mag + p[Nr-2][0].mag ) / 2;
	p[Nr - 1][Nc - 1].tx = ( p[Nr - 1][Nc - 2].tx + p[Nr - 2][Nc - 1].tx ) / 2;
	p[Nr - 1][Nc - 1].ty = ( p[Nr - 1][Nc - 2].ty + p[Nr - 2][Nc - 1].ty ) / 2;
	p[Nr - 1][Nc - 1].mag = ( p[Nr - 1][Nc - 2].mag + p[Nr - 2][Nc - 1].mag ) / 2;

	normalize();
}

void ETF::set2(imatrix& image) 
// get gradients from gradient map!
{
	int i, j;
	double MAX_VAL = 1020.; // 255 + 2 * 255 + 255 (Maximum possible Sobel value)
	double v[2];

	max_grad = -1.;

	matrix tmp(Nr, Nc);

	for (i = 1; i < Nr - 1; i++) { 
		for (j = 1; j < Nc - 1; j++) {
			////////////////////////////////////////////////////////////////
			p[i][j].tx = (image[i+1][j-1] + 2*(double)image[i+1][j] + image[i+1][j+1] 
				- image[i-1][j-1] - 2*(double)image[i-1][j] - image[i-1][j+1]) / MAX_VAL;
			p[i][j].ty = (image[i-1][j+1] + 2*(double)image[i][j+1] + image[i+1][j+1]
				- image[i-1][j-1] - 2*(double)image[i][j-1] - image[i+1][j-1]) / MAX_VAL;
			/////////////////////////////////////////////
			v[0] = p[i][j].tx;
			v[1] = p[i][j].ty;
			//////////////////////////////////////////////
			tmp[i][j] = sqrt(p[i][j].tx * p[i][j].tx + p[i][j].ty * p[i][j].ty);

			if (tmp[i][j] > max_grad) {
				max_grad = tmp[i][j];
			}
		}
	}

	for (i = 1; i <= Nr - 2; i++) {
		tmp[i][0] = tmp[i][1];
		tmp[i][Nc - 1] = tmp[i][Nc - 2];
	}
	
	for (j = 1; j <= Nc - 2; j++) {
		tmp[0][j] = tmp[1][j];
		tmp[Nr - 1][j] = tmp[Nr - 2][j];
	}
	
	tmp[0][0] = ( tmp[0][1] + tmp[1][0] ) / 2;
	tmp[0][Nc-1] = ( tmp[0][Nc-2] + tmp[1][Nc-1] ) / 2;
	tmp[Nr-1][0] = ( tmp[Nr-1][1] + tmp[Nr-2][0] ) / 2;
	tmp[Nr - 1][Nc - 1] = ( tmp[Nr - 1][Nc - 2] + tmp[Nr - 2][Nc - 1] ) / 2;

	//TRACE("max_grad = %f\n", max_grad);

	imatrix gmag(Nr, Nc);

	// normalize the magnitude
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			tmp[i][j] /= max_grad;
			gmag[i][j] = round(tmp[i][j] * 255.0);
		}
	}

	for (i = 1; i < Nr - 1; i++) { 
		for (j = 1; j < Nc - 1; j++) {
			////////////////////////////////////////////////////////////////
			p[i][j].tx = (gmag[i+1][j-1] + 2*(double)gmag[i+1][j] + gmag[i+1][j+1] 
				- gmag[i-1][j-1] - 2*(double)gmag[i-1][j] - gmag[i-1][j+1]) / MAX_VAL;
			p[i][j].ty = (gmag[i-1][j+1] + 2*(double)gmag[i][j+1] + gmag[i+1][j+1]
				- gmag[i-1][j-1] - 2*(double)gmag[i][j-1] - gmag[i+1][j-1]) / MAX_VAL;
			/////////////////////////////////////////////
			v[0] = p[i][j].tx;
			v[1] = p[i][j].ty;
			p[i][j].tx = -v[1];
			p[i][j].ty = v[0];
			//////////////////////////////////////////////
			p[i][j].mag = sqrt(p[i][j].tx * p[i][j].tx + p[i][j].ty * p[i][j].ty);

			if (p[i][j].mag > max_grad) {
				max_grad = p[i][j].mag;
			}
		}
	}

	for (i = 1; i <= Nr - 2; i++) {
		p[i][0].tx = p[i][1].tx;
		p[i][0].ty = p[i][1].ty;
		p[i][0].mag = p[i][1].mag;
		p[i][Nc - 1].tx = p[i][Nc - 2].tx;
		p[i][Nc - 1].ty = p[i][Nc - 2].ty;
		p[i][Nc - 1].mag = p[i][Nc - 2].mag;
	}
	
	for (j = 1; j <= Nc - 2; j++) {
		p[0][j].tx = p[1][j].tx;
		p[0][j].ty = p[1][j].ty;
		p[0][j].mag = p[1][j].mag;
		p[Nr - 1][j].tx = p[Nr - 2][j].tx;
		p[Nr - 1][j].ty = p[Nr - 2][j].ty;
		p[Nr - 1][j].mag = p[Nr - 2][j].mag;
	}
	
	p[0][0].tx = ( p[0][1].tx + p[1][0].tx ) / 2;
	p[0][0].ty = ( p[0][1].ty + p[1][0].ty ) / 2;
	p[0][0].mag = ( p[0][1].mag + p[1][0].mag ) / 2;
	p[0][Nc-1].tx = ( p[0][Nc-2].tx + p[1][Nc-1].tx ) / 2;
	p[0][Nc-1].ty = ( p[0][Nc-2].ty + p[1][Nc-1].ty ) / 2;
	p[0][Nc-1].mag = ( p[0][Nc-2].mag + p[1][Nc-1].mag ) / 2;
	p[Nr-1][0].tx = ( p[Nr-1][1].tx + p[Nr-2][0].tx ) / 2;
	p[Nr-1][0].ty = ( p[Nr-1][1].ty + p[Nr-2][0].ty ) / 2;
	p[Nr-1][0].mag = ( p[Nr-1][1].mag + p[Nr-2][0].mag ) / 2;
	p[Nr - 1][Nc - 1].tx = ( p[Nr - 1][Nc - 2].tx + p[Nr - 2][Nc - 1].tx ) / 2;
	p[Nr - 1][Nc - 1].ty = ( p[Nr - 1][Nc - 2].ty + p[Nr - 2][Nc - 1].ty ) / 2;
	p[Nr - 1][Nc - 1].mag = ( p[Nr - 1][Nc - 2].mag + p[Nr - 2][Nc - 1].mag ) / 2;

	//TRACE("max_grad = %f\n", max_grad);

	normalize();

	//TRACE("max_grad = %f\n", max_grad);
}


inline void solve_eig(double A, double B, double C, double D, 
    double& lambda1, double& v1x, double& v1y, 
    double& lambda2, double& v2x, double& v2y ) 
// A (j11), B (j12), C (j21), D (j22)
// lambda1, (v1x, v1y)
// lambda2, (v2x, v2y)
{
    if(B*C <= 0.1e-20 ) {
        lambda1 = A; v1x = 1; v1y = 0;
        lambda2 = D; v2x = 0; v2y = 1;
        return;
    }

    double tr = A + D;
    double det = A * D - B * C;
    double S = sqrt( (tr/2)*(tr/2) - det );
    lambda1 = tr/2 + S;
    lambda2 = tr/2 - S;

    double SS = sqrt( MAX( ((A-D)/2) * ((A-D)/2) + B*C, 0.0 ) );
    if( A - D < 0 ) {
        v1x = C;
        v1y = - (A-D)/2 + SS;
        v2x = + (A-D)/2 - SS;
        v2y = B;
    } else {
        v2x = C;
        v2y = - (A-D)/2 - SS;
        v1x = + (A-D)/2 + SS;
        v1y = B;
    }

    double n1 = sqrt( v1x*v1x + v1y*v1y );
    v1x /= n1; v1y /= n1;
    double n2 = sqrt( v2x*v2x + v2y*v2y );
    v2x /= n2; v2y /= n2;
}

void ETF::tensor(imatrix& image) 
// create structure tensor
{
	int i, j;
	double MAX_VAL = 1020.; // 255 + 2 * 255 + 255 (Maximum possible Sobel value)
	double j11,j12,j21,j22, Ix, Iy;
	vector v1(2);
	//vector v2(2);
	//double e1, e2;

	max_grad = -1.;

	for (i = 1; i < Nr - 1; i++) { 
		for (j = 1; j < Nc - 1; j++) {
			////////////////////////////////////////////////////////////////
			Ix = (image[i+1][j-1] + 2*(double)image[i+1][j] + image[i+1][j+1] 
				- image[i-1][j-1] - 2*(double)image[i-1][j] - image[i-1][j+1]) / MAX_VAL;
			Iy = (image[i-1][j+1] + 2*(double)image[i][j+1] + image[i+1][j+1]
				- image[i-1][j-1] - 2*(double)image[i][j-1] - image[i+1][j-1]) / MAX_VAL;
			/////////////////////////////////////////////
			//Ix = p[i][j].tx;
			//Iy = p[i][j].ty;
			//p[i][j].tx = -v[1];
			//p[i][j].ty = v[0];
			//////////////////////////////////////////////
			j11 = Ix * Ix;
			j12 = Ix * Iy;
			j21 = Iy * Ix;
			j22 = Iy * Iy;

			//solve_eig(j11, j12, j21, j22, e1, v1[0], v1[1],	e2, v2[0], v2[1]);
			
			v1[0] = 2 * j12;
			v1[1] = j22 - j11 + sqrt( (j11-j22) * (j11-j22) + 4 * j12 * j12 );

			//v2[0] = 2 * j12;
			//v2[1] = j22 - j11 - sqrt( (j11-j22) * (j11-j22) + 4 * j12 * j12 );

			//e1 = 0.5 * ( j11 + j22 + sqrt( (j11 - j22) * (j11 - j22) + 4 * j12 * j12 ) );
			//e2 = 0.5 * ( j11 + j22 - sqrt( (j11 - j22) * (j11 - j22) + 4 * j12 * j12 ) );

			//TRACE("e1 = %f, [%f, %f]\n", e1, v1[0], v1[1]);
			//TRACE("e2 = %f, [%f, %f]\n\n", e2, v2[0], v2[1]);

			//////////////////////////////////////////////////////
			// v1 corresponds to the bigger eigen value!
			p[i][j].tx = -v1[1];
			p[i][j].ty = v1[0];
			/*
			if (e1 > e2) {
				p[i][j].tx = -v1[1];
				p[i][j].ty = v1[0];
			}
			else {
				p[i][j].tx = -v2[1];
				p[i][j].ty = v2[0];
			}
			*/

			/////////////////////////////////////
			p[i][j].mag = sqrt(Ix * Ix + Iy * Iy);

			if (p[i][j].mag > max_grad) {
				max_grad = p[i][j].mag;
			}
			///////////////////////////////////
		}
	}

	for (i = 1; i <= Nr - 2; i++) {
		p[i][0].tx = p[i][1].tx;
		p[i][0].ty = p[i][1].ty;
		p[i][0].mag = p[i][1].mag;
		p[i][Nc - 1].tx = p[i][Nc - 2].tx;
		p[i][Nc - 1].ty = p[i][Nc - 2].ty;
		p[i][Nc - 1].mag = p[i][Nc - 2].mag;
	}
	
	for (j = 1; j <= Nc - 2; j++) {
		p[0][j].tx = p[1][j].tx;
		p[0][j].ty = p[1][j].ty;
		p[0][j].mag = p[1][j].mag;
		p[Nr - 1][j].tx = p[Nr - 2][j].tx;
		p[Nr - 1][j].ty = p[Nr - 2][j].ty;
		p[Nr - 1][j].mag = p[Nr - 2][j].mag;
	}
	
	p[0][0].tx = ( p[0][1].tx + p[1][0].tx ) / 2;
	p[0][0].ty = ( p[0][1].ty + p[1][0].ty ) / 2;
	p[0][0].mag = ( p[0][1].mag + p[1][0].mag ) / 2;
	p[0][Nc-1].tx = ( p[0][Nc-2].tx + p[1][Nc-1].tx ) / 2;
	p[0][Nc-1].ty = ( p[0][Nc-2].ty + p[1][Nc-1].ty ) / 2;
	p[0][Nc-1].mag = ( p[0][Nc-2].mag + p[1][Nc-1].mag ) / 2;
	p[Nr-1][0].tx = ( p[Nr-1][1].tx + p[Nr-2][0].tx ) / 2;
	p[Nr-1][0].ty = ( p[Nr-1][1].ty + p[Nr-2][0].ty ) / 2;
	p[Nr-1][0].mag = ( p[Nr-1][1].mag + p[Nr-2][0].mag ) / 2;
	p[Nr - 1][Nc - 1].tx = ( p[Nr - 1][Nc - 2].tx + p[Nr - 2][Nc - 1].tx ) / 2;
	p[Nr - 1][Nc - 1].ty = ( p[Nr - 1][Nc - 2].ty + p[Nr - 2][Nc - 1].ty ) / 2;
	p[Nr - 1][Nc - 1].mag = ( p[Nr - 1][Nc - 2].mag + p[Nr - 2][Nc - 1].mag ) / 2;

	normalize();
}

void ETF::SmoothTensor(imatrix& image, double sigma)
// Smooth the 4 components of the tensor matrix
{
	int i, j;
	double MAX_VAL = 1020.; // 255 + 2 * 255 + 255 (Maximum possible Sobel value)
	double j11,j12,j21,j22, Ix, Iy;
	vector v1(2);
	//vector v2(2);
	//double e1, e2;
	int left, right, up, down;

	max_grad = -1.;
	
	mymatrix tensor[3];
	tensor[0].init(Nr, Nc);
	tensor[1].init(Nr, Nc);
	tensor[2].init(Nr, Nc);
	mymatrix tmp[3];
	tmp[0].init(Nr, Nc);
	tmp[1].init(Nr, Nc);
	tmp[2].init(Nr, Nc);

	for (j = 0; j < Nc; j++) {
		for (i = 0; i < Nr; i++) {		
			if (i == 0) left = 0; else left = i-1;
			if (i == Nr-1) right = Nr-1; else right = i+1;
			if (j == 0) down = 0; else down = j-1;
			if (j == Nc-1) up = Nc-1; else up = j+1;
			/////////////////////////////////////////////
			Ix = (image[right][down] + 2*(double)image[right][j] + image[right][up] 
				- image[left][down] - 2*(double)image[left][j] - image[left][up]) / MAX_VAL;
			Iy = (image[left][up] + 2*(double)image[i][up] + image[right][up]
				- image[left][down] - 2*(double)image[i][down] - image[right][down]) / MAX_VAL;
			/////////////////////////////////////////////
			p[i][j].mag = sqrt(Ix * Ix + Iy * Iy);

			if (p[i][j].mag > max_grad) {
				max_grad = p[i][j].mag;
			}
			//////////////////////////////////////////////
			j11 = Ix * Ix;
			j12 = Ix * Iy;
			j21 = Iy * Ix;
			j22 = Iy * Iy;

			tensor[0][i][j] = j11;
			tensor[1][i][j] = j12;
			tensor[2][i][j] = j22;
		}
	}
	
	//int MAX_GRADIENT = -1;
	double g0, g1, g2;
	int s, t;
	int x, y;
	double weight, w_sum;

	//int image_x = image.getRow();
	//int image_y = image.getCol();

	myvec GAU1;
	MakeGaussianVector(sigma, GAU1); 
	int half = GAU1.getMax()-1;

	//max_g = -1;
	//min_g = 10000000;
	for (j = 0; j < Nc; j++) {
		for (i = 0; i < Nr; i++) {
			g0 = g1 = g2 = 0.0;
			weight = w_sum = 0.0;
			for (s = -half; s <= half; s++) {
				x = i+s; y = j;
				if (x > Nr-1) x = Nr-1;
				else if (x < 0) x = 0;
				if (y > Nc-1) y = Nc-1;
				else if (y < 0) y = 0;
				weight = GAU1[ABS(s)];
				g0 += weight * tensor[0][x][y];
				g1 += weight * tensor[1][x][y];
				g2 += weight * tensor[2][x][y];
				w_sum += weight;
			}
			g0 /= w_sum;
			g1 /= w_sum;
			g2 /= w_sum;
			tmp[0][i][j] = g0;
			tmp[1][i][j] = g1;
			tmp[2][i][j] = g2;
		}
	}
	for (j = 0; j < Nc; j++) {
		for (i = 0; i < Nr; i++) {
			g0 = g1 = g2 = 0.0;
			weight = w_sum = 0.0;
			for (t = -half; t <= half; t++) {
					x = i; y = j+t;
					if (x > IMAGE_X-1) x = IMAGE_X-1;
					else if (x < 0) x = 0;
					if (y > IMAGE_Y-1) y = IMAGE_Y-1;
					else if (y < 0) y = 0;
					weight = GAU1[ABS(t)];
					g0 += weight * tmp[0][x][y];
					g1 += weight * tmp[1][x][y];
					g2 += weight * tmp[2][x][y];
					w_sum += weight;
			}
			g0 /= w_sum;
			g1 /= w_sum;
			g2 /= w_sum;
			tensor[0][i][j] = g0;
			tensor[1][i][j] = g1;
			tensor[2][i][j] = g2;
		}
	}
	for (j = 0; j < Nc; j++) {
		for (i = 0; i < Nr; i++) {		if (i == 0) left = 0; else left = i-1;
			//////////////////////////////////////////////
			j11 = tensor[0][i][j];
			j12 = tensor[1][i][j];
			j21 = j12;
			j22 = tensor[2][i][j];

			v1[0] = 2 * j12;
			v1[1] = j22 - j11 + sqrt( (j11-j22) * (j11-j22) + 4 * j12 * j12 );

			//v2[0] = 2 * j12;
			//v2[1] = j22 - j11 - sqrt( (j11-j22) * (j11-j22) + 4 * j12 * j12 );

			//e1 = 0.5 * ( j11 + j22 + sqrt( (j11 - j22) * (j11 - j22) + 4 * j12 * j12 ) );
			//e2 = 0.5 * ( j11 + j22 - sqrt( (j11 - j22) * (j11 - j22) + 4 * j12 * j12 ) );

			//TRACE("e1 = %f, [%f, %f]\n", e1, v1[0], v1[1]);
			//TRACE("e2 = %f, [%f, %f]\n\n", e2, v2[0], v2[1]);

			//////////////////////////////////////////////////////
			// v1 corresponds to the bigger eigen value!
			p[i][j].tx = -v1[1];
			p[i][j].ty = v1[0];
			/*
			if (e1 > e2) {
				p[i][j].tx = -v1[1];
				p[i][j].ty = v1[0];
			}
			else {
				p[i][j].tx = -v2[1];
				p[i][j].ty = v2[0];
			}
			*/

			
		}
	}

	normalize();

	
//	TRACE("max_g = %f\n", max_g);
//	TRACE("min_g = %f\n", min_g);
}

inline void make_unit(double& vx, double& vy)
{
	double mag = sqrt( vx*vx + vy*vy );
	if (mag != 0.0) { 
		vx /= mag; 
		vy /= mag;
	}
}

void ETF::normalize() 
{
	int i, j;
	
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			make_unit(p[i][j].tx, p[i][j].ty);
		}
	}
}

void ETF::orient180() 
// re-orient each vector so that it ranges in [0, 180]
{
	int i, j;
	
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			if (p[i][j].ty < 0.0) { 
				p[i][j].tx = -p[i][j].tx;
				p[i][j].ty = -p[i][j].ty;
			}
		}
	}
}

void ETF::Smooth(int half_w, int M)
{
	int	i, j, k;
	int MAX_GRADIENT = -1;
	double weight;
	int s, t;
	int x, y;
	double mag_diff;

	int image_x = getRow();
	int image_y = getCol();
    
	ETF e2; 

	e2.init(image_x, image_y); 
	e2.copy(*this); 

	//TRACE("image_x = %d, image_y = %d\n", image_x, image_y);

	double v[2], w[2], g[2];
	double angle;
	double factor;

	for (k = 0; k < M; k++) {
		////////////////////////
		// horizontal
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				g[0] = g[1] = 0.0;
				v[0] = p[i][j].tx;
				v[1] = p[i][j].ty;
				for (s = -half_w; s <= half_w; s++) {
					////////////////////////////////////////
					x = i+s; y = j;
					if (x > image_x-1) x = image_x-1;
					else if (x < 0) x = 0;
					if (y > image_y-1) y = image_y-1;
					else if (y < 0) y = 0;
					////////////////////////////////////////
					mag_diff = p[x][y].mag - p[i][j].mag; 
					//////////////////////////////////////////////////////
					w[0] = p[x][y].tx;
					w[1] = p[x][y].ty;
					////////////////////////////////
					factor = 1.0;
					angle = v[0] * w[0] + v[1] * w[1];
					if (angle < 0.0) {
						factor = -1.0; // reverse the direction
					}
					weight = mag_diff + 1; // [0, 2] // this part takes up 5 seconds more!
					//////////////////////////////////////////////////////
					g[0] += weight * p[x][y].tx * factor;
					g[1] += weight * p[x][y].ty * factor;
				}
				make_unit(g[0], g[1]);
				e2[i][j].tx = g[0];
				e2[i][j].ty = g[1];
			}
		}
		this->copy(e2);
		/////////////////////////////////
		// vertical
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				g[0] = g[1] = 0.0;
				v[0] = p[i][j].tx;
				v[1] = p[i][j].ty;
				for (t = -half_w; t <= half_w; t++) {
					////////////////////////////////////////
					x = i; y = j+t;
					if (x > image_x-1) x = image_x-1;
					else if (x < 0) x = 0;
					if (y > image_y-1) y = image_y-1;
					else if (y < 0) y = 0;
					////////////////////////////////////////
					mag_diff = p[x][y].mag - p[i][j].mag; 
					//////////////////////////////////////////////////////
					w[0] = p[x][y].tx;
					w[1] = p[x][y].ty;
					////////////////////////////////
					factor = 1.0;
					///////////////////////////////
					angle = v[0] * w[0] + v[1] * w[1];
					if (angle < 0.0) factor = -1.0; 
					/////////////////////////////////////////////////////////
					weight = mag_diff + 1; 
					//////////////////////////////////////////////////////
					g[0] += weight * p[x][y].tx * factor;
					g[1] += weight * p[x][y].ty * factor;
				}
				make_unit(g[0], g[1]);
				e2[i][j].tx = g[0];
				e2[i][j].ty = g[1];
			}
		}
		this->copy(e2);
	}
	////////////////////////////////////////////
}


#define dist2(x1, y1, x2, y2) sqrt( (((double)x1)-((double)x2))*(((double)x1)-((double)x2)) + (((double)y1)-((double)y2))*(((double)y1)-((double)y2)) )

void ETF::SmoothFull(int half_w, int M)
//////////////////////////////////////////////////////////////////////
// HERE we smooth with lower magnitude vectors too! (This is the BEST!)
// non-adaptive version
{
	int	i, j, k;
	int MAX_GRADIENT = -1;
	double weight;
	int s, t;
	int x, y;
	double mag_diff;
	//int d1;

	int image_x = getRow();
	int image_y = getCol();
    
	ETF e2; 

	e2.init(image_x, image_y); 
	e2.copy(*this); 

	double v[2], w[2], g[2];
	double angle;
	double factor;

	for (k = 0; k < M; k++) {
		////////////////////////
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				g[0] = g[1] = 0.0;
				v[0] = p[i][j].tx;
				v[1] = p[i][j].ty;
				for (t = -half_w; t <= half_w; t++) {
					for (s = -half_w; s <= half_w; s++) {
						////////////////////////////////////////
						x = i+s; y = j+t;
						if (x > IMAGE_X-1) x = IMAGE_X-1;
						else if (x < 0) x = 0;
						if (y > IMAGE_Y-1) y = IMAGE_Y-1;
						else if (y < 0) y = 0;
						/////////////////////////////////////////////////////////
						// circular kernel
						//d1 = (int)dist2(x, y, i, j);
						//if ( d1 > half_w ) continue; // this doesn't make any difference in speed
						/////////////////////////////////////
						////////////////////////////////////////
						mag_diff = p[x][y].mag - p[i][j].mag; 
						//////////////////////////////////////////////////////
						w[0] = p[x][y].tx;
						w[1] = p[x][y].ty;
						////////////////////////////////
						factor = 1.0;
						angle = v[0] * w[0] + v[1] * w[1];
						if (angle < 0.0) {
							factor = -1.0; // reverse the direction
						}
						weight = mag_diff + 1; // [0, 2] // this part takes up 5 seconds more!
						//////////////////////////////////////////////////////
						g[0] += weight * p[x][y].tx * factor;
						g[1] += weight * p[x][y].ty * factor;
					}
				}
				make_unit(g[0], g[1]);
				e2[i][j].tx = g[0];
				e2[i][j].ty = g[1];
			}
		}
		this->copy(e2);
	}
	////////////////////////////////////////////
}

void ETF::Smooth2(int half_w, int M)
// Smoothing without reversing the vector!
{
	int	i, j, k;
	int MAX_GRADIENT = -1;
	double weight;
	int s, t;
	int x, y;
	double mag_diff;

	int image_x = getRow();
	int image_y = getCol();
    
	ETF e2; 

	e2.init(image_x, image_y); 
	e2.copy(*this); 

	double v[2], w[2], g[2];

	for (k = 0; k < M; k++) {
		////////////////////////
		// horizontal
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				g[0] = g[1] = 0.0;
				v[0] = p[i][j].tx;
				v[1] = p[i][j].ty;
				for (s = -half_w; s <= half_w; s++) {
					////////////////////////////////////////
					x = i+s; y = j;
					if (x > IMAGE_X-1) x = IMAGE_X-1;
					else if (x < 0) x = 0;
					if (y > IMAGE_Y-1) y = IMAGE_Y-1;
					else if (y < 0) y = 0;
					////////////////////////////////////////
					mag_diff = p[x][y].mag - p[i][j].mag; 
					//////////////////////////////////////////////////////
					w[0] = p[x][y].tx;
					w[1] = p[x][y].ty;
					////////////////////////////////
					weight = mag_diff + 1; // [0, 2] // this part takes up 5 seconds more!
					//////////////////////////////////////////////////////
					g[0] += weight * p[x][y].tx;
					g[1] += weight * p[x][y].ty;
				}
				make_unit(g[0], g[1]);
				e2[i][j].tx = g[0];
				e2[i][j].ty = g[1];
			}
		}
		this->copy(e2);
		/////////////////////////////////
		// vertical
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				g[0] = g[1] = 0.0;
				v[0] = p[i][j].tx;
				v[1] = p[i][j].ty;
				for (t = -half_w; t <= half_w; t++) {
					////////////////////////////////////////
					x = i; y = j+t;
					if (x > IMAGE_X-1) x = IMAGE_X-1;
					else if (x < 0) x = 0;
					if (y > IMAGE_Y-1) y = IMAGE_Y-1;
					else if (y < 0) y = 0;
					////////////////////////////////////////
					mag_diff = p[x][y].mag - p[i][j].mag; 
					//////////////////////////////////////////////////////
					w[0] = p[x][y].tx;
					w[1] = p[x][y].ty;
					/////////////////////////////////////////////////////////
					weight = mag_diff + 1; 
					//weight = tanh(mag_diff) + 1; // [0, 2] // this part takes up 5 seconds more!
					//////////////////////////////////////////////////////
					g[0] += weight * p[x][y].tx;
					g[1] += weight * p[x][y].ty;
				}
				make_unit(g[0], g[1]);
				e2[i][j].tx = g[0];
				e2[i][j].ty = g[1];
			}
		}
		this->copy(e2);
	}
	////////////////////////////////////////////
}

void ETF::Smooth3(int half_w, int M)
// revserse the vector direction, but make sure it doesn't get out of the tensor range!
{
	int	i, j, k;
	int MAX_GRADIENT = -1;
	double weight;
	int s, t;
	int x, y;
	double mag_diff;

	int image_x = getRow();
	int image_y = getCol();
    
	ETF e2; 

	e2.init(image_x, image_y); 
	e2.copy(*this); 

	double v[2], w[2], g[2];
	double angle;
	double factor;

	for (k = 0; k < M; k++) {
		////////////////////////
		// horizontal
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				g[0] = g[1] = 0.0;
				v[0] = p[i][j].tx;
				v[1] = p[i][j].ty;
				for (s = -half_w; s <= half_w; s++) {
					////////////////////////////////////////
					x = i+s; y = j;
					if (x > IMAGE_X-1) x = IMAGE_X-1;
					else if (x < 0) x = 0;
					if (y > IMAGE_Y-1) y = IMAGE_Y-1;
					else if (y < 0) y = 0;
					////////////////////////////////////////
					mag_diff = p[x][y].mag - p[i][j].mag; 
					//////////////////////////////////////////////////////
					w[0] = p[x][y].tx;
					w[1] = p[x][y].ty;
					////////////////////////////////
					factor = 1.0;
					angle = v[0] * w[0] + v[1] * w[1];
					if (angle < 0.0) {
						factor = -1.0; // reverse the direction
					}
					weight = mag_diff + 1; // [0, 2] // this part takes up 5 seconds more!
					//////////////////////////////////////////////////////
					g[0] += weight * p[x][y].tx * factor;
					g[1] += weight * p[x][y].ty * factor;
				}
				make_unit(g[0], g[1]);
				/////////////////////////////////////////////////
				// Reverse if it is out of the tensor range!
				if ( g[0] > 0 ) { g[0] = -g[0], g[1] = -g[1]; }
				///////////////////////////////////////////////
				e2[i][j].tx = g[0];
				e2[i][j].ty = g[1];
			}
		}
		this->copy(e2);
		/////////////////////////////////
		// vertical
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				g[0] = g[1] = 0.0;
				v[0] = p[i][j].tx;
				v[1] = p[i][j].ty;
				for (t = -half_w; t <= half_w; t++) {
					////////////////////////////////////////
					x = i; y = j+t;
					if (x > IMAGE_X-1) x = IMAGE_X-1;
					else if (x < 0) x = 0;
					if (y > IMAGE_Y-1) y = IMAGE_Y-1;
					else if (y < 0) y = 0;
					////////////////////////////////////////
					mag_diff = p[x][y].mag - p[i][j].mag; 
					//////////////////////////////////////////////////////
					w[0] = p[x][y].tx;
					w[1] = p[x][y].ty;
					////////////////////////////////
					factor = 1.0;
					///////////////////////////////
					angle = v[0] * w[0] + v[1] * w[1];
					if (angle < 0.0) factor = -1.0; 
					/////////////////////////////////////////////////////////
					weight = mag_diff + 1; 
					//////////////////////////////////////////////////////
					g[0] += weight * p[x][y].tx * factor;
					g[1] += weight * p[x][y].ty * factor;
				}
				make_unit(g[0], g[1]);
				/////////////////////////////////////////////////
				// Reverse if it is out of the tensor range!
				if ( g[0] > 0 ) { g[0] = -g[0], g[1] = -g[1]; }
				///////////////////////////////////////////////
				e2[i][j].tx = g[0];
				e2[i][j].ty = g[1];
			}
		}
		this->copy(e2);
	}
	////////////////////////////////////////////
}
