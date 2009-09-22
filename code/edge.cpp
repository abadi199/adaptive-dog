//#include <stdio.h>
//#include <stdlib.h>
#include "stdafx.h"
#include <cmath>
#include <deque>
#include <iomanip>
#include <fstream>
using namespace std;

#include "Cube.h"
#include "CubeDoc.h"
#include "CubeView.h"
#include "NameDialog.h"
#include "globals.h"
#include "Field.h"

#define TRUE_ 1
#define FALSE_ 0

#define SQRT2	1.414

#define EPS		0.0001
#define EPS2	0.000001

//#define MAX_ELEMENTS 10000 // for EL
#define MAX_ELEMENTS 500000 // for NPR

#define HEAP_FULL(n) (n == MAX_ELEMENTS-1)
#define HEAP_EMPTY(n) (!n)

#define RGB_GETRED(rgb)    ((rgb) & 0xff) 
#define RGB_GETGREEN(rgb)    (((rgb) >> 8) & 0xff) 
#define RGB_GETBLUE(rgb)    (((rgb) >> 16) & 0xff) 

#define ABS(x)	( ((x)>0) ? (x) : (-(x)) )
#define dist2(x1, y1, x2, y2) sqrt( (((double)x1)-((double)x2))*(((double)x1)-((double)x2)) + (((double)y1)-((double)y2))*(((double)y1)-((double)y2)) )

int GLOBAL_CANNY_DONE = 0; // when global canny is called, it is set to 1

imatrix gray, gray2, emap, gray3, dog;
mymatrix output1, output2;
double lambda;
cimatrix cmap;
imatrix off_mark, on_mark, pixel_mark;

inline double gauss(double x, double mean, double sigma)
{
	return ( exp( (-(x-mean)*(x-mean)) / (2*sigma*sigma) ) / sqrt(PI * 2.0 * sigma * sigma) );
	//sqrt(PI * 2.0 * sigma * sigma) is the sum of integral over all x
	//return ( exp( (-x*x) / (2*sigma*sigma) ) / (PI * 2.0 * sigma * sigma) );
}

inline double gauss2(double x, double mean, double sigma)
{
	
	return ( exp( (-(x-mean)*(x-mean)) / (2*sigma*sigma) ) );
	//sqrt(PI * 2.0 * sigma * sigma) is the sum of integral over all x
	//return ( exp( (-x*x) / (2*sigma*sigma) ) / (PI * 2.0 * sigma * sigma) );
}

inline double gauss1D_norm(double x, double mean, double sigma)
// function value is now between [0, 1]	
// incorrect version
{
	
	//return ( exp( (-(x-mean)*(x-mean)) / (2*sigma*sigma) ) );
	//sqrt(PI * 2.0 * sigma * sigma) is the sum of integral over all x
	return ( exp( (-x*x) / (2*sigma*sigma) ) / (PI * 2.0 * sigma * sigma) );
}

inline double gauss1D_norm2(double x, double mean, double sigma)
// function value is now between [0, 1]	
// correct version
{
	
	//return ( exp( (-(x-mean)*(x-mean)) / (2*sigma*sigma) ) );
	//sqrt(PI * 2.0 * sigma * sigma) is the sum of integral over all x
	return ( exp( (-x*x) / (2*sigma*sigma) ) / sqrt(PI * 2.0 * sigma * sigma) );
}

inline double gauss1D_bl(double x, double mean, double sigma)
// used in bilateral filter (no normalization)
{
	
	//return ( exp( (-(x-mean)*(x-mean)) / (2*sigma*sigma) ) );
	//sqrt(PI * 2.0 * sigma * sigma) is the sum of integral over all x
	return ( exp( (-x*x) / (2*sigma*sigma) ) );
}

//double gau[30];

double gau[GAU_MAX];
double gau2[30];
int gau_w = 10;
int gau_w2 = 10;
double max_grad2;
extern NameDialog* dlg;

int MakeGaussMask(double sigma, double* gau)
{
	int i;
	double sum = 0.0;

	gau[0] = gauss((double)0.0, 0.0, sigma);
	sum += gau[0];
	//TRACE("\n%.5f ", gau[0]);
	for (i=1; i < GAU_MAX; i++) {
		gau[i] = gauss((double)i, 0.0, sigma);
		sum += 2 * gau[i];
		//TRACE("%.5f ", gau[i]);
		if (gau[i] < 0.005) {
			gau_w = i+1;
			break;
		}
	}
	//TRACE("\nsum = %.5f ", sum);
	//TRACE("sigma = %f ", sigma);
	//TRACE("gau_w = %d\n", gau_w);
	return gau_w;
}

matrix GAU;
int GAU_W[255];

void MakeGaussMatrix(double factor, matrix& GAU)
// factor is used when converting the index to sigma
{
	int i, j;
	double sum = 0.0;
	double sigma;
	//double factor;
	double max_sigma;
	double threshold;

	//factor = 3.0;

	threshold = 0.001;
	max_sigma = exp(factor) - 0.5;

	///////////////////
	// determine the maximum column number required
	i = 0;
	while(1) {
		i++;
		if ( gauss((double)i, 0.0, max_sigma) < threshold )
			break;
	}
	GAU.init(256, i+1); // size of GAU matrix
	GAU.zero(); 
	TRACE("max_col = %d\n", GAU.getCol());

	for (i = 0; i < 256; i++) {
		sigma = exp(factor * (double)i / 255.0) - 0.5;
		TRACE("sigma[%d] = %f\n", i, sigma);
		//sigma = max_sigma * (double)i / 255.0;
		GAU[i][0] = gauss((double)0.0, 0.0, sigma);
		//TRACE("\n%.5f ", gau[0]);
		for (j = 1; j < GAU.getCol(); j++) {
			GAU[i][j] = gauss((double)j, 0.0, sigma);
			//TRACE("%.5f ", gau[i]);
			if (GAU[i][j] < threshold) {
				GAU_W[i] = j+1;
				break;
			}
		}
	}
	//TRACE("\nsum = %.5f ", sum);
	//TRACE("sigma = %f ", sigma);
	//TRACE("gau_w = %d\n", gau_w);
}

void MakeGaussVector(double sigma, vector& GAU)
// factor is used when converting the index to sigma
{
	int i, j;
	double sum = 0.0;
	//double sigma;
	//double factor;
	//double max_sigma;
	double threshold;

	//factor = 3.0;

	threshold = 0.001;

	///////////////////
	// determine the maximum number of vector elements required
	i = 0;
	while(1) {
		i++;
		if ( gauss((double)i, 0.0, sigma) < threshold )
			break;
	}
	GAU.init(i+1); // size of GAU vector
	GAU.zero(); 
	//TRACE("max_elm = %d\n", GAU.getMax());

	GAU[0] = gauss((double)0.0, 0.0, sigma);
	for (j = 1; j < GAU.getMax(); j++) {
		GAU[j] = gauss((double)j, 0.0, sigma);
	}
}

inline void MakeGaussVectorBL(double sigma, vector& GAU, int max_index)
// used for speeding up bilateral filtering
// use index (instead of threshold value) to cut
{
	int i;
	//double sum = 0.0;
	//double sigma;
	//double factor;
	//double max_sigma;
	//double threshold;

	//factor = 3.0;

	//threshold = 0.001;

	///////////////////
	// determine the maximum number of vector elements required
	i = 0;
	//while(1) {
	//	i++;
	//	if ( gauss((double)i, 0.0, sigma) < threshold )
	//		break;
	//}
	GAU.init(max_index); // size of GAU vector
	GAU.zero(); 
	//TRACE("max_elm = %d\n", GAU.getMax());

	GAU[0] = gauss1D_bl((double)0.0, 0.0, sigma);
	for (i = 1; i < GAU.getMax(); i++) {
		GAU[i] = gauss1D_bl((double)i, 0.0, sigma);
	}
}

//matrix tmp_x, tmp_y;
matrix tmp_x, tmp_y;
matrix G_x, G_y, G_mag;
//matrix thin_edge;
imatrix thin_edge;
matrix scale_map;

void GaussSmooth(int image_x, int image_y, imatrix& image, int gau_w)
{
	int	i, j, k, i1, i2, j1, j2;
	int MAX_GRADIENT = -1;

	imatrix tmp(IMAGE_X, IMAGE_Y);
	double	x, y;
	
	/// copy image to tmp
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			tmp[i][j] = image[i][j];
		}
	}
		
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			x = gau[0] * tmp[i][j];	
			y = gau[0] * tmp[i][j];
			//TRACE("x = %f, y = %f\n", x, y);
			for (k = 1; k < gau_w; k++) {
				i1 = (i+k)%image_x;
				i2 = (i-k+image_x)%image_x;
				x += gau[k] * tmp[i1][j] + gau[k] * tmp[i2][j];
				j1 = (j+k)%image_y;
				j2 = (j-k+image_y)%image_y;
				y += gau[k] * tmp[i][j1] + gau[k] * tmp[i][j2];
			}
			tmp_x[i][j] = x;
			tmp_y[i][j] = y;
			//TRACE("x = %f, y = %f\n", x, y);
			//if (x > 255) x = 255;
			//if (y > 255) y = 255;
			image[i][j] = (int)((x + y)/2);
			if (image[i][j] > 255) image[i][j] = 255;
			//image[i][j] = (int)y;
		}
	}

}



void GaussSmooth5(imatrix& image, double sigma)
// circular kernel (equivalent to using 2D Gaussian mask)
// with normalization!
{
	int	i, j, k;
	int MAX_GRADIENT = -1;
	double g, max_g, min_g;
	int s, t;
	int x, y;
	//int half = gau_w-1;
	double weight, w_sum;

	int image_x = image.getRow();
	int image_y = image.getCol();

	vector GAU1;
	MakeGaussVector(sigma, GAU1); 
	int half = GAU1.getMax()-1;

	imatrix tmp(IMAGE_X, IMAGE_Y);
	
	/// copy image to tmp
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			tmp[i][j] = image[i][j];
		}
	}
		
	max_g = -1;
	min_g = 10000000;
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			//count = 0;
			g = 0.0;
			weight = w_sum = 0.0;
			for (s = -half; s <= half; s++) {
				for (t = -half; t <= half; t++) {
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					if ( k > half ) continue; 
					//////////////////////////////////////////////////
					if (x > IMAGE_X-1) x = IMAGE_X-1;
					else if (x < 0) x = 0;
					if (y > IMAGE_Y-1) y = IMAGE_Y-1;
					else if (y < 0) y = 0;
					//count++;					
					weight = GAU1[k];
					g += weight * tmp[x][y];
					w_sum += weight;
					//TRACE("k = %d\n", k);
					//TRACE("gau[%d] = %f\n", k, gau[k]);
				}
			}
			
			//g /= count;
			g /= w_sum;
			//TRACE("g = %f\n", g);
			//TRACE("x = %f, y = %f\n", x, y);
			//if (x > 255) x = 255;
			//if (y > 255) y = 255;
			if (g > max_g) max_g = g;
			if (g < min_g) min_g = g;
			//TRACE("max_g = %f\n", max_g);
			image[i][j] = round(g);
			//if (image[i][j] > 255) image[i][j] = 255;
			//image[i][j] = (int)y;
		}
	}
	
	TRACE("max_g = %f\n", max_g);
	TRACE("min_g = %f\n", min_g);
}



void GaussColorSmooth5(cimatrix& image, double sigma)
// circular kernel (equivalent to using 2D Gaussian mask)
// with normalization!
{
	int	i, j, k;
	int MAX_GRADIENT = -1;
	double r, g, b, max_g, min_g;
	int s, t;
	int x, y;
	//int half = gau_w-1;
	double weight, w_sum;

	int image_x = image.getRow();
	int image_y = image.getCol();

	vector GAU1;
	MakeGaussVector(sigma, GAU1); 
	int half = GAU1.getMax()-1;

	cimatrix tmp(IMAGE_X, IMAGE_Y);
	
	/// copy image to tmp
	tmp.copy(image);
		
	max_g = -1;
	min_g = 10000000;
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			//count = 0;
			r = g = b = 0.0;
			weight = w_sum = 0.0;
			for (s = -half; s <= half; s++) {
				for (t = -half; t <= half; t++) {
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					if ( k > half ) continue; 
					//////////////////////////////////////////////////
					if (x > IMAGE_X-1) x = IMAGE_X-1;
					else if (x < 0) x = 0;
					if (y > IMAGE_Y-1) y = IMAGE_Y-1;
					else if (y < 0) y = 0;
					//count++;					
					weight = GAU1[k];
					r += weight * tmp[x][y].r;
					g += weight * tmp[x][y].g;
					b += weight * tmp[x][y].b;
					w_sum += weight;
					//TRACE("k = %d\n", k);
					//TRACE("gau[%d] = %f\n", k, gau[k]);
				}
			}
			
			//g /= count;
			r /= w_sum;
			g /= w_sum;
			b /= w_sum;
			//TRACE("g = %f\n", g);
			//TRACE("x = %f, y = %f\n", x, y);
			//if (x > 255) x = 255;
			//if (y > 255) y = 255;
			//if (g > max_g) max_g = g;
			//if (g < min_g) min_g = g;
			//TRACE("max_g = %f\n", max_g);
			image[i][j].r = round(r);
			image[i][j].g = round(g);
			image[i][j].b = round(b);
			//if (image[i][j] > 255) image[i][j] = 255;
			//image[i][j] = (int)y;
		}
	}
}



void FindNearestPoint(imatrix& image, int i, int j, int &x, int &y)
{
	int dir_x[4] = {1, 0, -1, 0};
	int dir_y[4] = {0, 1, 0, -1};
	int cur_x, cur_y;
	int max_step, k;
	int cur_dir_idx;
	double dist, min_dist;
	int min_x, min_y;
	int flag, flag2;

	cur_x = i;
	cur_y = j;
	max_step = 1;
	cur_dir_idx = 0;
	flag = 0;
	flag2 = 0;
	min_dist = 1000000000000000;
	if (image[cur_x][cur_y] == 0) {
		x = cur_x;
		y = cur_y;
		return;
	}
	while(1) {
		for (k = 0; k < max_step; k++) {
			cur_x += dir_x[cur_dir_idx];
			cur_y += dir_y[cur_dir_idx];
			if (cur_x < 0 || cur_x > IMAGE_X-1 || cur_y < 0 || cur_y > IMAGE_Y-1)
				continue;
			if (image[cur_x][cur_y] == 0) {
				flag = 1;
				dist = dist2(i, j, cur_x, cur_y);
				//TRACE("dist = %f\n", dist);
				if ( dist < min_dist ) {
					min_dist = dist;
					min_x = cur_x;
					min_y = cur_y;
				}
			}
		}
		cur_dir_idx = (cur_dir_idx + 1) % 4;
		for (k = 0; k < max_step; k++) {
			cur_x += dir_x[cur_dir_idx];
			cur_y += dir_y[cur_dir_idx];
			if (cur_x < 0 || cur_x > IMAGE_X-1 || cur_y < 0 || cur_y > IMAGE_Y-1)
				continue;
			if (image[cur_x][cur_y] == 0) {
				flag = 1;
				dist = dist2(i, j, cur_x, cur_y);
				if ( dist < min_dist ) {
					min_dist = dist;
					min_x = cur_x;
					min_y = cur_y;
				}
			}
		}
		cur_dir_idx = (cur_dir_idx + 1) % 4;
		if (flag == 1 && flag2 == 0)
			flag2 = 1; // go one more round
		else if (flag == 1 && flag2 == 1) {
			x = min_x;
			y = min_y;
			return;
		}
		max_step++;
	}
}

void FindNearestPoint2(imatrix& image, int i, int j, int &x, int &y)
// Around curved edges, you have to go some more steps to make sure there're no closer points
{
	int dir_x[4] = {1, 0, -1, 0};
	int dir_y[4] = {0, 1, 0, -1};
	int cur_x, cur_y;
	int max_step, k;
	int cur_dir_idx;
	double dist, min_dist, first_dist;
	int min_x, min_y;
	int target_count, count;
	double target_dist;

	cur_x = i;
	cur_y = j;
	max_step = 1;
	cur_dir_idx = 0;
	count = -1;
	min_dist = 1000000000000000;
	first_dist = -1;
	if (image[cur_x][cur_y] == 0) {
		x = cur_x;
		y = cur_y;
		return;
	}
	while(1) {
		for (k = 0; k < max_step; k++) {
			cur_x += dir_x[cur_dir_idx];
			cur_y += dir_y[cur_dir_idx];
			if (cur_x < 0 || cur_x > IMAGE_X-1 || cur_y < 0 || cur_y > IMAGE_Y-1)
				continue;
			if (image[cur_x][cur_y] == 0) {
				dist = dist2(i, j, cur_x, cur_y);
				//TRACE("dist = %f\n", dist);
				if (count < 0) { // this is the first point found
					count = 0;
					target_dist = dist - dist/1.414; 
					target_count = round(2 * (target_dist+1.0) ); // steps to go more
				}
				if ( dist < min_dist ) {
					min_dist = dist;
					min_x = cur_x;
					min_y = cur_y;
				}
			}
		}
		cur_dir_idx = (cur_dir_idx + 1) % 4;
		for (k = 0; k < max_step; k++) {
			cur_x += dir_x[cur_dir_idx];
			cur_y += dir_y[cur_dir_idx];
			if (cur_x < 0 || cur_x > IMAGE_X-1 || cur_y < 0 || cur_y > IMAGE_Y-1)
				continue;
			if (image[cur_x][cur_y] == 0) {
				dist = dist2(i, j, cur_x, cur_y);
				if (count < 0) { // this is the first point found
					count = 0;
					target_dist = dist - dist/1.414; 
					target_count = round(2 * (target_dist+1.0) ); // steps to go more
				}
				if (first_dist < 0)
					first_dist = dist; // this is the first point found
				if ( dist < min_dist ) {
					min_dist = dist;
					min_x = cur_x;
					min_y = cur_y;
				}
			}
		}
		cur_dir_idx = (cur_dir_idx + 1) % 4;
		////////////////////////////////////////
		if (count > target_count) {
			x = min_x;
			y = min_y;
			return;
		}
		///////////////////////////////////////
		if (count >= 0) count++;
		max_step++;
	}
}


void DistanceField(int image_x, int image_y, imatrix& image, double factor)
{
	int	i, j;
	int MAX_GRADIENT = -1;
	double g, max_g, min_g;
	int count;
	int x, y;
	int half = gau_w-1;
	double dist;


	imatrix tmp(IMAGE_X, IMAGE_Y);
	matrix tmp_dist(IMAGE_X, IMAGE_Y);
	
	/// copy image to tmp
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			tmp[i][j] = image[i][j];
		}
	}
		
	max_g = -1;
	min_g = 10000000;
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			count = 0;
			g = 0.0;
			FindNearestPoint2(tmp, i, j, x, y);

			dist = dist2(x, y, i, j);
			if (dist > max_g) max_g = dist;
			if (dist < min_g) min_g = dist;
					
			//image[i][j] = (int)dist;
			tmp_dist[i][j] = dist;
			//if (image[i][j] > 255) image[i][j] = 255;
			//image[i][j] = (int)y;
		}
	}
	
	TRACE("max_g = %f\n", max_g);
	TRACE("min_g = %f\n", min_g);
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			//image[i][j] = (int)(255 * (image[i][j] / max_g));
			//////////////////////////////////////////////////
			tmp_dist[i][j] = 255.0 * ((tmp_dist[i][j]-min_g) / (max_g-min_g));
			//tmp_dist[i][j] += 255.0 * ( sin( PI / 2.0 * (tmp_dist[i][j]-min_g) / (max_g-min_g) ) );
			//if (tmp_dist[i][j] > 255) tmp_dist[i][j] = 255;
			//tmp_dist[i][j] = 55 + 200.0 * ((tmp_dist[i][j]-min_g) / (max_g-min_g));
			////////////////////////////////////////////////////////////////////
			// level adjustment
			image[i][j] = (int)(255.0 * (1.0-pow((double)(1.0-tmp_dist[i][j]/255.0), factor)) );
			//image[i][j] = 255.0 * sin( PI / 2.0 * image[i][j] / 255.0 );
			//if (image[i][j] > 255) image[i][j] = 255;
            //image[i][j] = (int)(255.0 * (1.0-pow((double)(1.0-tmp_dist[i][j]/255.0), factor)) );
			//image[i][j] = (int)(255 * (pow((double)2, (double)image[i][j]/255.0)-1));
			//image[i][j] = (int)(255 * pow((double)(image[i][j]/255.0), 5.0));
			//TRACE("image[%d][%d] = %d\n", i, j, image[i][j]);
			//if (image[i][j] > 0) image[i][j] = 255;
			///////////////////////////////////////////////////////////
		}
	}
}

inline double gauss2D(double x, double y, double mean, double sigma)
// function value is now between [0, 1] 
{
	double g;
	g = exp( -( (x-mean)*(x-mean) + (y-mean)*(y-mean) ) / (2*sigma*sigma) );
	g /= (2 * PI * sigma * sigma); 
	return g;
	//sqrt(PI * 2.0 * sigma * sigma) is the sum of integral over all x
	//return ( exp( (-x*x) / (2*sigma*sigma) ) / (PI * 2.0 * sigma * sigma) );
}

matrix gau_mask;

void Make2DGaussMask(double sigma)
// 2D Gauss mask
{
	int i, j;
	double sum = 0.0;
	int N;
	int half;
	double x, y;
	int size;
	double threshold = 0.0001;

	matrix tmp;

	i = 0;
	while (1) { // compute N (Gaussian mask size)
		//TRACE("[%d] Gauss value = %f\n", i, gauss2D((double)i, 0.0, 0.0, sigma));
		if (gauss2D((double)i, 0.0, 0.0, sigma) < threshold)
			break;
		i++;
	}

	N = i * 2 + 5;
	//TRACE("N = %d\n", N);
	//N = 31;
	half = N / 2;

	tmp.init(N, N);
	tmp.zero();
	
	//gau[0] = gauss((double)0.0, 0.0, sigma);
	//sum += gau[0];
	//TRACE("\n%.5f ", gau[0]);
	for (j = 0; j < N; j++) {
		for (i = 0; i < N; i++) {
			x = i - half; 
			y = j - half; 
			tmp[i][j] = gauss2D(x, y, 0.0, sigma);
			//sum += 2 * gau[i];
			//TRACE("%.5f ", gau[i]);
		}
	}
	for (i = half; i < N; i++) {
		if (tmp[i][half] < threshold) {
			size = 2 * (i - half) + 1;
			break;
		}
	}
	gau_mask.init(size, size);
	gau_mask.zero();

	for (j = 0; j < size; j++) {
		for (i = 0; i < size; i++) {
			gau_mask[i][j] = tmp[half-size/2+i][half-size/2+j];
		}
	}

	TRACE("%dx%d GAUSSIAN MASK\n", size, size);
	
	//TRACE("\nsum = %.5f ", sum);
	//TRACE("sigma = %f ", sigma);
	//TRACE("gau_w = %d\n", gau_w);
}

void GaussBlur(imatrix& image, double sigma)
// using 2D Gaussian mask
{
	int	i, j;
	int MAX_GRADIENT = -1;
	double g, max_g, min_g, w_sum;
	int count;
	int s, t;
	int x, y;
	int half;
	int N;

	int image_x = image.getRow();
	int image_y = image.getCol();

	Make2DGaussMask(sigma);

	N = gau_mask.getRow();
	half = N / 2;

	imatrix tmp(IMAGE_X, IMAGE_Y);
	matrix tmp2(IMAGE_X, IMAGE_Y);
	tmp2.zero();
	
	/// copy image to tmp
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			tmp[i][j] = image[i][j];
		}
	}
		
	max_g = -1;
	min_g = 10000000;
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			count = 0;
			g = 0.0;
			for (s = -half; s <= half; s++) {
				for (t = -half; t <= half; t++) {
					x = i+s; y = j+t;
                    if (x > IMAGE_X-1) x = IMAGE_X-1;
					else if (x < 0) x = 0;
					if (y > IMAGE_Y-1) y = IMAGE_Y-1;
					else if (y < 0) y = 0;
					//count++;					
					g += gau_mask[s+half][t+half] * tmp[x][y];
					//TRACE("k = %d\n", k);
					//TRACE("gau[%d] = %f\n", k, gau[k]);
				}
			}
			//TRACE("g = %f\n", g);
			//g /= count;
			//TRACE("x = %f, y = %f\n", x, y);
			//if (x > 255) x = 255;
			//if (y > 255) y = 255;
			if (g > max_g) max_g = g;
			if (g < min_g) min_g = g;
			//TRACE("max_g = %f\n", max_g);
			tmp2[i][j] = g;
			//if (image[i][j] > 255) image[i][j] = 255;
			//image[i][j] = (int)y;
		}
	}
	
	TRACE("max_g = %f\n", max_g);
	TRACE("min_g = %f\n", min_g);

	////////////////////////////////////
	// Normalization
	w_sum = 0.0;
	for (j = 0; j < N; j++) {
		for (i = 0; i < N; i++) {
			w_sum += gau_mask[i][j];
		}
	}
	//TRACE("N = %d\n", N);
	TRACE("w_sum = %f\n", w_sum);
	//gau_mask.print();

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			tmp2[i][j] /= w_sum;
			//image[i][j] = (int)(255 * (image[i][j] / max_g));
			image[i][j] = (int)tmp2[i][j];
		}
	}
	////////////////////////////////////////////
}

void GaussBlurMemDC(CDC& dc, int image_x, int image_y, double sigma)
{
	int	i, j;
	int MAX_GRADIENT = -1;
	double g, max_g, min_g, w_sum;
	int count;
	int s, t;
	int x, y;

	Make2DGaussMask(sigma);

	int half = gau_mask.getRow() / 2;

	imatrix tmp(image_x, image_y);
	//imatrix image(image_x, image_y);

	//int gau_w = MakeGaussMask(sigma, gau);
	//int half = gau_w-1;

	GLubyte r;
	
	/// copy image to tmp
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			//tmp[i][j] = image[i][j];
			tmp[i][j] = RGB_GETRED(dc.GetPixel(i, image_y-1-j));
		}
	}
		
	max_g = -1;
	min_g = 10000000;
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			count = 0;
			g = 0.0;
			w_sum = 0.0;
			for (s = -half; s <= half; s++) {
				for (t = -half; t <= half; t++) {
					x = i+s; y = j+t;
                    if (x > image_x-1) x = image_x-1;
					else if (x < 0) x = 0;
					if (y > image_y-1) y = image_y-1;
					else if (y < 0) y = 0;
					//count++;					
					g += gau_mask[s+half][t+half] * tmp[x][y];
					//count++;					
					//g += gau[k] * tmp[x][y];
					w_sum += gau_mask[s+half][t+half];
					//TRACE("k = %d\n", k);
					//TRACE("gau[%d] = %f\n", k, gau[k]);
				}
			}
			//TRACE("g = %f\n", g);
			g /= w_sum; // weight normalize
			//TRACE("x = %f, y = %f\n", x, y);
			//if (x > 255) x = 255;
			//if (y > 255) y = 255;
			if (g > max_g) max_g = g;
			if (g < min_g) min_g = g;
			//TRACE("max_g = %f\n", max_g);
			//image[i][j] = (int)g;
			r = (GLubyte)g;
			dc.SetPixelV(x, IMAGE_Y-1-y, RGB(r, r, r));
			//if (image[i][j] > 255) image[i][j] = 255;
			//image[i][j] = (int)y;
		}
	}
	
	TRACE("max_g = %f\n", max_g);
	TRACE("min_g = %f\n", min_g);

	/*
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			//image[i][j] = (int)(255 * (image[i][j] / max_g));
			//image[i][j] = (int)(255 * ((image[i][j]-min_g) / (max_g-min_g)));
			////////////////////////////////////////////////////////////////////
			// Level adjustment
			//image[i][j] = (int)(255 * pow((double)(image[i][j]/255.0), 5.0));
			//TRACE("image[%d][%d] = %d\n", i, j, image[i][j]);
			//if (image[i][j] > 0) image[i][j] = 255;
			///////////////////////////////////////////////////////////
		}
	}
	*/
}

void GaussBlurBilateral(int image_x, int image_y, imatrix& image, double sigma, double sigma2, int max_itr)
// using 2D Gaussian mask
{
	int	i, j, u, v;
	int MAX_GRADIENT = -1;
	double g, max_g, min_g, w_sum;
	int s, t, k;
	int x, y;
	int half;
	int N;
	double w_sim;

	Make2DGaussMask(sigma);

	N = gau_mask.getRow();
	half = N / 2;

	matrix tmp_mask(N, N);
	tmp_mask.zero();

	imatrix tmp(IMAGE_X, IMAGE_Y);
	tmp.copy(image);

	matrix tmp2(IMAGE_X, IMAGE_Y);
	tmp2.zero();
	
	max_g = -1;
	min_g = 10000000;
	for (k = 0; k < max_itr; k++) {
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				g = 0.0;
				for (s = 0; s < N; s++) {
					for (t = 0; t < N; t++) {
						x = i-half+s; y = j-half+t;
						if (x > IMAGE_X-1) x = IMAGE_X-1;
						else if (x < 0) x = 0;
						if (y > IMAGE_Y-1) y = IMAGE_Y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						w_sim = gauss1D_norm((double)(tmp[i][j] - tmp[x][y])/255.0, 0.0, sigma2);
						//count++;					
						tmp_mask[s][t] = w_sim * gau_mask[s][t];
						g += tmp_mask[s][t] * tmp[x][y];
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
					}
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				tmp2[i][j] = g;
				////////////////////////////////////
				// Normalization
				w_sum = 0.0;
				for (v = 0; v < N; v++) {
					for (u = 0; u < N; u++) {
						w_sum += tmp_mask[u][v];
					}
				}
				tmp2[i][j] /= w_sum;
				image[i][j] = round(tmp2[i][j]);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		tmp.copy(image);
	}
	//TRACE("max_g = %f\n", max_g);
	//TRACE("min_g = %f\n", min_g);
}

#define dist3(x, y, z) sqrt(((double)x)*((double)x) + ((double)y)*((double)y) + ((double)z)*((double)z))

void BilateralColor(int image_x, int image_y, cimatrix& image, double sigma, double sigma2, int max_itr)
{
	int	i, j;
	int MAX_GRADIENT = -1;
	double max_g, min_g;
	int s, t, k;
	int x, y;
	int half;
	int N;
	double weight;
	double val_r, val_g, val_b;
	double c_val_r, c_val_g, c_val_b;
	double w_sum, sum_r, sum_g, sum_b;
	double c_dist;

	Make2DGaussMask(sigma);

	N = gau_mask.getRow();
	half = N / 2;

	cimatrix tmp(IMAGE_X, IMAGE_Y);
	tmp.copy(image);
	
	max_g = -1;
	min_g = 10000000;
	for (k = 0; k < max_itr; k++) {
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum_r = sum_g = sum_b = 0.0;
				weight = w_sum = 0.0;
				c_val_r = tmp[i][j].r;
				c_val_g = tmp[i][j].g;
				c_val_b = tmp[i][j].b;
				//c_val = (c_val_r + c_val_g + c_val_b) / 3.0;
				for (s = 0; s < N; s++) {
					for (t = 0; t < N; t++) {
						x = i-half+s; y = j-half+t;
						if (x > IMAGE_X-1) x = IMAGE_X-1;
						else if (x < 0) x = 0;
						if (y > IMAGE_Y-1) y = IMAGE_Y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val_r = tmp[x][y].r;
						val_g = tmp[x][y].g;
						val_b = tmp[x][y].b;
						//val = (val_r + val_g + val_b) / 3.0;
						//weight = gauss1D_norm((c_val - val)/255.0, 0.0, sigma2);
						c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
						//weight = gauss1D_norm( c_dist/255.0, 0.0, sigma2);
						weight = gauss1D_norm( c_dist/255.0/1.732, 0.0, sigma2);
						//count++;					
						weight *= gau_mask[s][t];
						sum_r += weight * val_r;
						sum_g += weight * val_g;
						sum_b += weight * val_b;
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
					}
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum_r /= w_sum;
				sum_g /= w_sum;
				sum_b /= w_sum;
				//image[i][j] = round(tmp2);
				image[i][j].r = round(sum_r);
				image[i][j].g = round(sum_g);
				image[i][j].b = round(sum_b);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		tmp.copy(image);
	}
	//TRACE("max_g = %f\n", max_g);
	//TRACE("min_g = %f\n", min_g);
}

void BilateralColorTomasi(int image_x, int image_y, cimatrix& image, double sigma, double sigma2, int max_itr)
{
	int	i, j;
	int MAX_GRADIENT = -1;
	double max_g, min_g;
	int s, t, k;
	int x, y;
	int half;
	int N;
	double weight;
	double val_r, val_g, val_b;
	double c_val_r, c_val_g, c_val_b;
	double w_sum, sum_r, sum_g, sum_b;
	double c_dist;

	Make2DGaussMask(sigma);

	N = gau_mask.getRow();
	half = N / 2;

	vector GAU2;

	MakeGaussVectorBL(sigma2, GAU2, 300); // 300 array elements created

	cimatrix tmp(IMAGE_X, IMAGE_Y);
	tmp.copy(image);
	
	max_g = -1;
	min_g = 10000000;
	for (k = 0; k < max_itr; k++) {
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum_r = sum_g = sum_b = 0.0;
				weight = w_sum = 0.0;
				c_val_r = tmp[i][j].r;
				c_val_g = tmp[i][j].g;
				c_val_b = tmp[i][j].b;
				//c_val = (c_val_r + c_val_g + c_val_b) / 3.0;
				for (s = 0; s < N; s++) {
					for (t = 0; t < N; t++) {
						x = i-half+s; y = j-half+t;
						if (x > IMAGE_X-1) x = IMAGE_X-1;
						else if (x < 0) x = 0;
						if (y > IMAGE_Y-1) y = IMAGE_Y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val_r = tmp[x][y].r;
						val_g = tmp[x][y].g;
						val_b = tmp[x][y].b;
						//val = (val_r + val_g + val_b) / 3.0;
						//weight = gauss1D_norm((c_val - val)/255.0, 0.0, sigma2);
						c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
						//weight = gauss1D_norm( c_dist/255.0, 0.0, sigma2);
						//weight = gauss1D_norm( c_dist/255.0/1.732, 0.0, sigma2);
						//weight = gauss1D_bl( c_dist/1.732, 0.0, sigma2); // sigma2 = 10.0
						weight = GAU2[ (int)(c_dist/1.732) ]; // like Gaussian!
						//count++;					
						weight *= gau_mask[s][t];
						sum_r += weight * val_r;
						sum_g += weight * val_g;
						sum_b += weight * val_b;
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
					}
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum_r /= w_sum;
				sum_g /= w_sum;
				sum_b /= w_sum;
				//image[i][j] = round(tmp2);
				image[i][j].r = round(sum_r);
				image[i][j].g = round(sum_g);
				image[i][j].b = round(sum_b);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		tmp.copy(image);
	}
	//TRACE("max_g = %f\n", max_g);
	//TRACE("min_g = %f\n", min_g);
}

void BilateralColorChui(cimatrix& image, double sigma, double sigma2, int max_itr)
{
	int	i, j;
	int MAX_GRADIENT = -1;
	double max_g, min_g;
	int s, t, k;
	int x, y;
	int half;
	int N;
	double weight;
	double val_r, val_g, val_b;
	double c_val_r, c_val_g, c_val_b;
	double w_sum, sum_r, sum_g, sum_b;
	double c_dist, s_dist;

	int image_x = image.getRow();
	int image_y = image.getCol();

	Make2DGaussMask(sigma);

	N = gau_mask.getRow();
	half = N / 2;

	//vector GAU2;

	//MakeGaussVectorBL(sigma2, GAU2, 300); // 300 array elements created

	cimatrix tmp(IMAGE_X, IMAGE_Y);
	tmp.copy(image);
	
	max_g = -1;
	min_g = 10000000;
	for (k = 0; k < max_itr; k++) {
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum_r = sum_g = sum_b = 0.0;
				weight = w_sum = 0.0;
				c_val_r = tmp[i][j].r;
				c_val_g = tmp[i][j].g;
				c_val_b = tmp[i][j].b;
				//c_val = (c_val_r + c_val_g + c_val_b) / 3.0;
				for (s = 0; s < N; s++) {
					for (t = 0; t < N; t++) {
						x = i-half+s; y = j-half+t;
						if (x > IMAGE_X-1) x = IMAGE_X-1;
						else if (x < 0) x = 0;
						if (y > IMAGE_Y-1) y = IMAGE_Y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val_r = tmp[x][y].r;
						val_g = tmp[x][y].g;
						val_b = tmp[x][y].b;
						//val = (val_r + val_g + val_b) / 3.0;
						//weight = gauss1D_norm((c_val - val)/255.0, 0.0, sigma2);
						c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
						//weight = gauss1D_norm( c_dist/255.0, 0.0, sigma2);
						//weight = gauss1D_norm( c_dist/255.0/1.732, 0.0, sigma2);
						//weight = gauss1D_bl( c_dist/1.732, 0.0, sigma2); // sigma2 = 10.0
						c_dist /= 1.732;
						s_dist = dist2((double)i, (double)j, (double)x, (double)y);
						c_dist *= s_dist;
						weight = gauss1D_bl( c_dist, 0.0, sigma2 ); // sigma2 = 10.0
						//weight = GAU2[ (int)(c_dist/1.732) ]; // like Gaussian!
						//count++;					
						//weight *= gau_mask[s][t];
						sum_r += weight * val_r;
						sum_g += weight * val_g;
						sum_b += weight * val_b;
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
					}
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum_r /= w_sum;
				sum_g /= w_sum;
				sum_b /= w_sum;
				//image[i][j] = round(tmp2);
				image[i][j].r = round(sum_r);
				image[i][j].g = round(sum_g);
				image[i][j].b = round(sum_b);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		tmp.copy(image);
	}
	//TRACE("max_g = %f\n", max_g);
	//TRACE("min_g = %f\n", min_g);
}

void BilateralToon(cimatrix& cmap, double sigma, double sigma2, int max_itr)
// use bilateral weight function from smoothed image 
{
	int	i, j;
	int MAX_GRADIENT = -1;
	double max_g, min_g;
	int s, t, k;
	int x, y;
	int half;
	int N;
	double weight;
	double val_r, val_g, val_b;
	double c_val_r, c_val_g, c_val_b;
	double w_sum, sum_r, sum_g, sum_b;
	double c_dist;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	Make2DGaussMask(sigma);

	N = gau_mask.getRow();
	half = N / 2;

	cimatrix tmp(IMAGE_X, IMAGE_Y);
	tmp.copy(cmap);
	
	max_g = -1;
	min_g = 10000000;
	for (k = 0; k < max_itr; k++) {
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum_r = sum_g = sum_b = 0.0;
				weight = w_sum = 0.0;
				c_val_r = tmp[i][j].r;
				c_val_g = tmp[i][j].g;
				c_val_b = tmp[i][j].b;
				//c_val = (c_val_r + c_val_g + c_val_b) / 3.0;
				for (s = 0; s < N; s++) {
					for (t = 0; t < N; t++) {
						x = i-half+s; y = j-half+t;
						if (x > IMAGE_X-1) x = IMAGE_X-1;
						else if (x < 0) x = 0;
						if (y > IMAGE_Y-1) y = IMAGE_Y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val_r = tmp[x][y].r;
						val_g = tmp[x][y].g;
						val_b = tmp[x][y].b;
						//val = (val_r + val_g + val_b) / 3.0;
						//weight = gauss1D_norm((c_val - val)/255.0, 0.0, sigma2);
						c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
						//weight = gauss1D_norm( c_dist/255.0, 0.0, sigma2);
						//weight = gauss1D_norm( c_dist/255.0/1.732, 0.0, sigma2);
						weight = gauss1D_bl( c_dist/1.732, 0.0, sigma2); // sigma2 = 10.0
						//count++;					
						weight *= gau_mask[s][t];
						sum_r += weight * val_r;
						sum_g += weight * val_g;
						sum_b += weight * val_b;
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
					}
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum_r /= w_sum;
				sum_g /= w_sum;
				sum_b /= w_sum;
				//image[i][j] = round(tmp2);
				cmap[i][j].r = round(sum_r);
				cmap[i][j].g = round(sum_g);
				cmap[i][j].b = round(sum_b);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		tmp.copy(cmap);
	}
	//TRACE("max_g = %f\n", max_g);
	//TRACE("min_g = %f\n", min_g);
}

void BilateralColorSep(int image_x, int image_y, cimatrix& cmap, double sigma, double sigma2, int max_itr)
// Separable Bilateral filter
{
	int	i, j;
	int MAX_GRADIENT = -1;
	double max_g, min_g;
	int s, t, k;
	int x, y;
	double weight;
	double val_r, val_g, val_b;
	double c_val_r, c_val_g, c_val_b;
	double w_sum, sum_r, sum_g, sum_b;
	double c_dist;

	vector GAU1;
	MakeGaussVector(sigma, GAU1); // length of the kernel
	//MakeGaussVector(sigma2, GAU2); // width of the kernel
	//MakeGaussVector(sigma3, GAU3); // for similarity function

	//half_l = GAU1.getMax()-1;
	int half_w = GAU1.getMax()-1;

	//Make2DGaussMask(sigma);
	vector GAU2;
	MakeGaussVectorBL(sigma2, GAU2, 300); // 300 array elements created

	//N = gau_mask.getRow();
	//half = N / 2;

	cimatrix tmp(IMAGE_X, IMAGE_Y);
	tmp.copy(cmap);

	max_g = -1;
	min_g = 10000000;
	for (k = 0; k < max_itr; k++) {
		///////////////////////////////
		// x-direction bilateral filter
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum_r = sum_g = sum_b = 0.0;
				weight = w_sum = 0.0;
				c_val_r = tmp[i][j].r;
				c_val_g = tmp[i][j].g;
				c_val_b = tmp[i][j].b;
				//c_val = (c_val_r + c_val_g + c_val_b) / 3.0;
				for (s = -half_w; s <= half_w; s++) {
						x = i+s; 
						y = j;
						if (x > IMAGE_X-1) x = IMAGE_X-1;
						else if (x < 0) x = 0;
						if (y > IMAGE_Y-1) y = IMAGE_Y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val_r = tmp[x][y].r;
						val_g = tmp[x][y].g;
						val_b = tmp[x][y].b;
						//val = (val_r + val_g + val_b) / 3.0;
						weight = GAU1[ABS(s)];
						c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
						//weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma2);
						weight *= GAU2[ (int)(c_dist/1.732) ]; // like Gaussian!
						//count++;					
						//weight *= gau_mask[s][t];
						sum_r += weight * val_r;
						sum_g += weight * val_g;
						sum_b += weight * val_b;
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum_r /= w_sum;
				sum_g /= w_sum;
				sum_b /= w_sum;
				//image[i][j] = round(tmp2);
				cmap[i][j].r = round(sum_r);
				cmap[i][j].g = round(sum_g);
				cmap[i][j].b = round(sum_b);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		tmp.copy(cmap);

		///////////////////////////////
		// y-direction bilateral filter
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum_r = sum_g = sum_b = 0.0;
				weight = w_sum = 0.0;
				c_val_r = tmp[i][j].r;
				c_val_g = tmp[i][j].g;
				c_val_b = tmp[i][j].b;
				//c_val = (c_val_r + c_val_g + c_val_b) / 3.0;
				for (t = -half_w; t <= half_w; t++) {
						x = i; 
						y = j+t;
						if (x > IMAGE_X-1) x = IMAGE_X-1;
						else if (x < 0) x = 0;
						if (y > IMAGE_Y-1) y = IMAGE_Y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val_r = tmp[x][y].r;
						val_g = tmp[x][y].g;
						val_b = tmp[x][y].b;
						//val = (val_r + val_g + val_b) / 3.0;
						weight = GAU1[ABS(t)];
						//weight = gauss1D_norm((c_val - val)/255.0, 0.0, sigma2);
						c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
						//weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma2);
						weight *= GAU2[ (int)(c_dist/1.732) ]; // like Gaussian!
						//count++;					
						//weight *= gau_mask[s][t];
						sum_r += weight * val_r;
						sum_g += weight * val_g;
						sum_b += weight * val_b;
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum_r /= w_sum;
				sum_g /= w_sum;
				sum_b /= w_sum;
				//image[i][j] = round(tmp2);
				cmap[i][j].r = round(sum_r);
				cmap[i][j].g = round(sum_g);
				cmap[i][j].b = round(sum_b);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		tmp.copy(cmap);
	}
	//TRACE("max_g = %f\n", max_g);
	//TRACE("min_g = %f\n", min_g);
}

void BilateralColorFull(cimatrix& image, double sigma, double sigma2, int max_itr)
// full kernel
// use table for bilateral weight
{
	int	i, j;
	int MAX_GRADIENT = -1;
	double max_g, min_g;
	int s, t, k, z;
	int x, y;
	int half;
	double weight;
	double val_r, val_g, val_b;
	double c_val_r, c_val_g, c_val_b;
	double w_sum, sum_r, sum_g, sum_b;
	double c_dist;

	int image_x = image.getRow();
	int image_y = image.getCol();

	//Make2DGaussMask(sigma);

	//N = gau_mask.getRow();

	//half = N / 2;

	vector GAU1;
	MakeGaussVector(sigma, GAU1);
	half = GAU1.getMax()-1;

	TRACE("half = %d\n", half);

	vector GAU2;
	MakeGaussVectorBL(sigma2, GAU2, 300); // 300 array elements created

	cimatrix tmp(IMAGE_X, IMAGE_Y);
	tmp.copy(image);
	
	max_g = -1;
	min_g = 10000000;
	for (k = 0; k < max_itr; k++) {
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum_r = sum_g = sum_b = 0.0;
				weight = w_sum = 0.0;
				c_val_r = tmp[i][j].r;
				c_val_g = tmp[i][j].g;
				c_val_b = tmp[i][j].b;
				//c_val = (c_val_r + c_val_g + c_val_b) / 3.0;
				for (s = -half; s <= half; s++) {
					for (t = -half; t <= half; t++) {
						x = i+s; y = j+t;
						/////////////////////////////////////////////////////////
						// circular kernel
						z = (int)dist2(x, y, i, j);
						//k = round( dist2(x, y, i, j) );
						if ( z > half ) continue; 
						//////////////////////////////////////////////
						if (x > IMAGE_X-1) x = IMAGE_X-1;
						else if (x < 0) x = 0;
						if (y > IMAGE_Y-1) y = IMAGE_Y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val_r = tmp[x][y].r;
						val_g = tmp[x][y].g;
						val_b = tmp[x][y].b;
						//val = (val_r + val_g + val_b) / 3.0;
						//weight = gauss1D_norm((c_val - val)/255.0, 0.0, sigma2);
						//c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
						//weight = gauss1D_norm( c_dist/255.0, 0.0, sigma2);
						//weight = gauss1D_norm( c_dist/255.0/1.732, 0.0, sigma2);

						//weight = GAU1[ABS(s)];
						c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
						//weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma2);
						weight = GAU2[ (int)(c_dist/1.732) ]; // like Gaussian!

						//count++;				
						weight *= GAU1[z];
						//weight *= gau_mask[s][t];
						sum_r += weight * val_r;
						sum_g += weight * val_g;
						sum_b += weight * val_b;
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
					}
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum_r /= w_sum;
				sum_g /= w_sum;
				sum_b /= w_sum;
				//image[i][j] = round(tmp2);
				image[i][j].r = round(sum_r);
				image[i][j].g = round(sum_g);
				image[i][j].b = round(sum_b);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		tmp.copy(image);
	}
	//TRACE("max_g = %f\n", max_g);
	//TRACE("min_g = %f\n", min_g);
}



void BilateralGraySep(imatrix& gray, double sigma, double sigma2, int max_itr)
// Separable Bilateral filter
{
	int	i, j;
	int MAX_GRADIENT = -1;
	double max_g, min_g;
	int s, t, k;
	int x, y;
	double weight;
	double val;
	double c_val;
	double w_sum, sum;
	double c_dist;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	vector GAU1;
	MakeGaussVector(sigma, GAU1); // length of the kernel
	//MakeGaussVector(sigma2, GAU2); // width of the kernel
	//MakeGaussVector(sigma3, GAU3); // for similarity function

	//half_l = GAU1.getMax()-1;
	int half_w = GAU1.getMax()-1;

	//Make2DGaussMask(sigma);
	vector GAU2;
	MakeGaussVectorBL(sigma2, GAU2, 300); // 300 array elements created

	//N = gau_mask.getRow();
	//half = N / 2;

	imatrix tmp(image_x, image_y);
	tmp.copy(gray);

	max_g = -1;
	min_g = 10000000;
	for (k = 0; k < max_itr; k++) {
		///////////////////////////////
		// x-direction bilateral filter
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum = 0.0;
				weight = w_sum = 0.0;
				c_val = tmp[i][j];
				//c_val = (c_val_r + c_val_g + c_val_b) / 3.0;
				for (s = -half_w; s <= half_w; s++) {
						x = i+s; 
						y = j;
						if (x > image_x-1) x = image_x-1;
						else if (x < 0) x = 0;
						if (y > image_y-1) y = image_y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val = tmp[x][y];
						//val = (val_r + val_g + val_b) / 3.0;
						weight = GAU1[ABS(s)];
						c_dist = fabs(c_val-val);
						//c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
						//weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma2);
						weight *= GAU2[ (int)c_dist ]; // like Gaussian!
						//weight *= GAU2[ (int)(c_dist/1.732) ]; // like Gaussian!
						//count++;					
						//weight *= gau_mask[s][t];
						sum += weight * val;
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum /= w_sum;
				//image[i][j] = round(tmp2);
				gray[i][j] = round(sum);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		tmp.copy(gray);

		///////////////////////////////
		// y-direction bilateral filter
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum = 0.0;
				weight = w_sum = 0.0;
				c_val = tmp[i][j];
				for (t = -half_w; t <= half_w; t++) {
						x = i; 
						y = j+t;
						if (x > image_x-1) x = image_x-1;
						else if (x < 0) x = 0;
						if (y > image_y-1) y = image_y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val = tmp[x][y];
						//val = (val_r + val_g + val_b) / 3.0;
						weight = GAU1[ABS(t)];
						//weight = gauss1D_norm((c_val - val)/255.0, 0.0, sigma2);
						c_dist = fabs(c_val-val);
						//weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma2);
						weight *= GAU2[ (int)(c_dist) ]; // like Gaussian!
						//count++;					
						//weight *= gau_mask[s][t];
						sum += weight * val;
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum /= w_sum;
				//image[i][j] = round(tmp2);
				gray[i][j] = round(sum);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		tmp.copy(gray);
	}
	//TRACE("max_g = %f\n", max_g);
	//TRACE("min_g = %f\n", min_g);
}

void LinearInterp(cimatrix& cmap1, cimatrix& cmap2, double t, cimatrix& cmap3)
{
	int i, j;
	GLubyte r, g, b;
	double rd, gd, bd;
	double r1, g1, b1, r2, g2, b2;
	double h1, s1, v1, h2, s2, v2; 
	double v3;

	int image_x = cmap1.getRow();
	int image_y = cmap1.getCol();

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			r1 = (double)cmap1[i][j].r / 255.0;
			g1 = (double)cmap1[i][j].g / 255.0;
			b1 = (double)cmap1[i][j].b / 255.0;
			RGB2HSV(r1, g1, b1, h1, s1, v1);
			
			r2 = (double)cmap2[i][j].r / 255.0;
			g2 = (double)cmap2[i][j].g / 255.0;
			b2 = (double)cmap2[i][j].b / 255.0;
			RGB2HSV(r2, g2, b2, h2, s2, v2);
			
			//s3 = (1-t) * s1 + t * s2;
			
			//v3 = (1-t) * v1 + t * v2; 

			v3 = v2 + t * (v2 - v1);

			if (v3 > 1.0) v3 = 1.0; if (v3 < 0.0) v3 = 0.0;

			HSV2RGB(h2, s2, v3, rd, gd, bd);

			rd *= 255;
			gd *= 255;
			bd *= 255;

			r = (GLubyte)rd;
			g = (GLubyte)gd;
			b = (GLubyte)bd;
			
			if (r < 0) r = 0; if (r > 255) r = 255;
			if (g < 0) g = 0; if (g > 255) g = 255;
			if (b < 0) b = 0; if (b > 255) b = 255;
			
			cmap3[i][j].r = r;
			cmap3[i][j].g = g;
			cmap3[i][j].b = b;
		}
	}
}

void LinearInterp2(cimatrix& cmap1, cimatrix& cmap2, double t, cimatrix& cmap3)
{
	int i, j;
	GLubyte r, g, b;
	double rd, gd, bd;
	double r1, g1, b1, r2, g2, b2;

	int image_x = cmap1.getRow();
	int image_y = cmap1.getCol();

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			
			r1 = (double)cmap1[i][j].r;
			g1 = (double)cmap1[i][j].g;
			b1 = (double)cmap1[i][j].b;

			r2 = (double)cmap2[i][j].r;
			g2 = (double)cmap2[i][j].g;
			b2 = (double)cmap2[i][j].b;
			
			//s3 = (1-t) * s1 + t * s2;
			//v3 = (1-t) * v1 + t * v2;
			rd = (1-t) * r1 + t * r2;
			gd = (1-t) * g1 + t * g2;
			bd = (1-t) * b1 + t * b2;

			//if (v3 > 1.0) v3 = 1.0; if (v3 < 0.0) v3 = 0.0;

			//HSV2RGB(h2, s2, v3, rd, gd, bd);

			//rd *= 255;
			//gd *= 255;
			//bd *= 255;

			r = (GLubyte)rd;
			g = (GLubyte)gd;
			b = (GLubyte)bd;
			
			if (r < 0) r = 0; if (r > 255) r = 255;
			if (g < 0) g = 0; if (g > 255) g = 255;
			if (b < 0) b = 0; if (b > 255) b = 255;
			
			cmap3[i][j].r = r;
			cmap3[i][j].g = g;
			cmap3[i][j].b = b;
		}
	}
}

void LinearInterp3(cimatrix& cmap1, cimatrix& cmap2, double t, cimatrix& cmap3)
{
	int i, j;
	GLubyte r, g, b;
	double rd, gd, bd;
	double r1, g1, b1, r2, g2, b2;
	double h1, s1, v1, h2, s2, v2; 
	double v3;

	int image_x = cmap1.getRow();
	int image_y = cmap1.getCol();

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			r1 = (double)cmap1[i][j].r / 255.0;
			g1 = (double)cmap1[i][j].g / 255.0;
			b1 = (double)cmap1[i][j].b / 255.0;
			RGB2HSV(r1, g1, b1, h1, s1, v1);
			
			r2 = (double)cmap2[i][j].r / 255.0;
			g2 = (double)cmap2[i][j].g / 255.0;
			b2 = (double)cmap2[i][j].b / 255.0;
			RGB2HSV(r2, g2, b2, h2, s2, v2);
			
			//s3 = (1-t) * s1 + t * s2;
			
			//v3 = (1-t) * v1 + t * v2; 
			
			v3 = v2 + 0.5 * (v2 - v1);
			
			//if (v2 - v1 > 0.0) v3 = v2 + (1.0 - v2) * 0.1;	
			//else if (v2 - v1 < 0.0) v3 = v2 - (v2 - 0.0) * 0.1;	
			//v3 = max(v2, v1);

			if (v3 > 1.0) v3 = 1.0; if (v3 < 0.0) v3 = 0.0;

			HSV2RGB(h2, s2, v3, rd, gd, bd);

			rd *= 255;
			gd *= 255;
			bd *= 255;

			r = (GLubyte)rd;
			g = (GLubyte)gd;
			b = (GLubyte)bd;
			
			if (r < 0) r = 0; if (r > 255) r = 255;
			if (g < 0) g = 0; if (g > 255) g = 255;
			if (b < 0) b = 0; if (b > 255) b = 255;
			
			cmap3[i][j].r = r;
			cmap3[i][j].g = g;
			cmap3[i][j].b = b;
		}
	}
}

void BilateralToonSep(CDC& dc, cimatrix& cmap, double sigma, double sigma2, int max_itr)
// Separable Bilateral filter
// Toon version
{
	int	i, j;
	int MAX_GRADIENT = -1;
	double max_g, min_g;
	int s, t, k;
	int x, y;
	double weight;
	double val_r, val_g, val_b;
	double c_val_r, c_val_g, c_val_b;
	double w_sum, sum_r, sum_g, sum_b;
	double c_dist;
	//double z;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	vector GAU1;
	MakeGaussVector(sigma, GAU1); // length of the kernel
	//MakeGaussVector(sigma2, GAU2); // width of the kernel
	//MakeGaussVector(sigma3, GAU3); // for similarity function

	//half_l = GAU1.getMax()-1;
	int half_w = GAU1.getMax()-1;

	//Make2DGaussMask(sigma);

	//N = gau_mask.getRow();
	//half = N / 2;

	cimatrix tmp(image_x, image_y);
	tmp.copy(cmap);

	cimatrix cmap2(image_x, image_y);
	cmap2.copy(cmap);
	//GaussColSmoothSep(cmap2, 2.0);

	max_g = -1;
	min_g = 10000000;
	for (k = 0; k < max_itr; k++) {
		///////////////////////////////
		GaussColSmoothSep(cmap2, 2.0);
		///////////////////////////////
		// x-direction bilateral filter
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum_r = sum_g = sum_b = 0.0;
				weight = w_sum = 0.0;
				c_val_r = cmap2[i][j].r;
				c_val_g = cmap2[i][j].g;
				c_val_b = cmap2[i][j].b;
				//c_val = (c_val_r + c_val_g + c_val_b) / 3.0;
				for (s = -half_w; s <= half_w; s++) {
						x = i+s; 
						y = j;
						if (x > image_x-1) x = image_x-1;
						else if (x < 0) x = 0;
						if (y > image_y-1) y = image_y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val_r = cmap2[x][y].r;
						val_g = cmap2[x][y].g;
						val_b = cmap2[x][y].b;
						//val = (val_r + val_g + val_b) / 3.0;
						weight = GAU1[ABS(s)];
						//weight = gauss1D_norm((c_val - val)/255.0, 0.0, sigma2);
						c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
						//weight = gauss1D_norm( c_dist/255.0, 0.0, sigma2);
						//weight *= gauss1D_norm( c_dist/255.0/1.732, 0.0, sigma2);
						//weight *= gauss1D_norm( c_dist/255.0/1.732, 0.0, sigma2);
						//weight *= gauss1D_bl( c_dist/255.0/1.732, 0.0, sigma2);
						weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma2);
						//count++;					
						//weight *= gau_mask[s][t];
						sum_r += weight * tmp[x][y].r;
						sum_g += weight * tmp[x][y].g;
						sum_b += weight * tmp[x][y].b;
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum_r /= w_sum;
				sum_g /= w_sum;
				sum_b /= w_sum;
				//image[i][j] = round(tmp2);
				cmap[i][j].r = round(sum_r);
				cmap[i][j].g = round(sum_g);
				cmap[i][j].b = round(sum_b);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		tmp.copy(cmap);

		///////////////////////////////
		// y-direction bilateral filter
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum_r = sum_g = sum_b = 0.0;
				weight = w_sum = 0.0;
				c_val_r = cmap2[i][j].r;
				c_val_g = cmap2[i][j].g;
				c_val_b = cmap2[i][j].b;
				//c_val = (c_val_r + c_val_g + c_val_b) / 3.0;
				for (t = -half_w; t <= half_w; t++) {
						x = i; 
						y = j+t;
						if (x > image_x-1) x = image_x-1;
						else if (x < 0) x = 0;
						if (y > image_y-1) y = image_y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val_r = cmap2[x][y].r;
						val_g = cmap2[x][y].g;
						val_b = cmap2[x][y].b;
						//val = (val_r + val_g + val_b) / 3.0;
						weight = GAU1[ABS(t)];
						//weight = gauss1D_norm((c_val - val)/255.0, 0.0, sigma2);
						c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
						//weight = gauss1D_norm( c_dist/255.0, 0.0, sigma2);
						//weight *= gauss1D_norm( c_dist/255.0/1.732, 0.0, sigma2);
						//weight *= gauss1D_bl( c_dist/255.0/1.732, 0.0, sigma2);
						weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma2);
						//count++;					
						//weight *= gau_mask[s][t];
						sum_r += weight * tmp[x][y].r;
						sum_g += weight * tmp[x][y].g;
						sum_b += weight * tmp[x][y].b;
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum_r /= w_sum;
				sum_g /= w_sum;
				sum_b /= w_sum;
				//image[i][j] = round(tmp2);
				cmap[i][j].r = round(sum_r);
				cmap[i][j].g = round(sum_g);
				cmap[i][j].b = round(sum_b);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		
		//LinearInterp(cmap2, cmap, 1.5, cmap);
		//LinearInterp3(cmap2, cmap, 1.5, cmap);
		DrawColorImage(memDC, image_x, image_y, cmap);
		dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);

		tmp.copy(cmap);
		//LinearInterp2(cmap2, cmap, z, cmap);
	}
	//z = 2.0;
	//LinearInterp(cmap2, cmap, z, cmap);
	//TRACE("max_g = %f\n", max_g);
	//TRACE("min_g = %f\n", min_g);
}

void BilateralGrayToonSep(imatrix& gray, double sigma, double sigma2, int max_itr)
// Separable Bilateral filter
// Toon version: gray image
// Get bilateral weight function from smoothed image 
{
	int	i, j;
	int MAX_GRADIENT = -1;
	double max_g, min_g;
	int s, t, k;
	int x, y;
	double weight;
	double val;
	double c_val;
	double w_sum, sum;
	double c_dist;
	//double z;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	vector GAU1;
	MakeGaussVector(sigma, GAU1); // length of the kernel
	//MakeGaussVector(sigma2, GAU2); // width of the kernel
	//MakeGaussVector(sigma3, GAU3); // for similarity function

	//half_l = GAU1.getMax()-1;
	int half_w = GAU1.getMax()-1;

	//Make2DGaussMask(sigma);

	//N = gau_mask.getRow();
	//half = N / 2;

	imatrix tmp(image_x, image_y);
	tmp.copy(gray);

	imatrix gray2(image_x, image_y);
	gray2.copy(gray);
	//GaussColSmoothSep(gray2, 2.0);
	GaussSmoothSep(gray2, 2.0);

	max_g = -1;
	min_g = 10000000;
	for (k = 0; k < max_itr; k++) {
		///////////////////////////////
		//////////////////////////////////
		//GaussColSmoothSep(gray2, 1.0);
		/////////////////////////////////
		///////////////////////////////
		// x-direction bilateral filter
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum = 0.0;
				weight = w_sum = 0.0;
				c_val = gray2[i][j];
				//c_val = (c_val_r + c_val_g + c_val_b) / 3.0;
				for (s = -half_w; s <= half_w; s++) {
						x = i+s; 
						y = j;
						if (x > image_x-1) x = image_x-1;
						else if (x < 0) x = 0;
						if (y > image_y-1) y = image_y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val = gray2[x][y];
						//val = (val_r + val_g + val_b) / 3.0;
						weight = GAU1[ABS(s)];
						//weight = gauss1D_norm((c_val - val)/255.0, 0.0, sigma2);
						c_dist = sqrt( (c_val-val) * (c_val-val) );
						//weight = gauss1D_norm( c_dist/255.0, 0.0, sigma2);
						//weight *= gauss1D_norm( c_dist/255.0/1.732, 0.0, sigma2);
						//weight *= gauss1D_norm( c_dist/255.0/1.732, 0.0, sigma2);
						//weight *= gauss1D_bl( c_dist/255.0/1.732, 0.0, sigma2);
						weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma2);
						//count++;					
						//weight *= gau_mask[s][t];
						sum += weight * tmp[x][y];
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum /= w_sum;
				//image[i][j] = round(tmp2);
				gray[i][j] = round(sum);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		tmp.copy(gray);

		///////////////////////////////
		// y-direction bilateral filter
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				sum = 0.0;
				weight = w_sum = 0.0;
				c_val = gray2[i][j];
				//c_val = (c_val_r + c_val_g + c_val_b) / 3.0;
				for (t = -half_w; t <= half_w; t++) {
						x = i; 
						y = j+t;
						if (x > image_x-1) x = image_x-1;
						else if (x < 0) x = 0;
						if (y > image_y-1) y = image_y-1;
						else if (y < 0) y = 0;
						//w_sim = gauss2((double)tmp[i][j] - tmp[x][y], 0.0, sigma2);
						val = gray2[x][y];
						//val = (val_r + val_g + val_b) / 3.0;
						weight = GAU1[ABS(t)];
						//weight = gauss1D_norm((c_val - val)/255.0, 0.0, sigma2);
						c_dist = sqrt( (c_val-val) * (c_val-val) );
						//c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
						//weight = gauss1D_norm( c_dist/255.0, 0.0, sigma2);
						//weight *= gauss1D_norm( c_dist/255.0/1.732, 0.0, sigma2);
						//weight *= gauss1D_bl( c_dist/255.0/1.732, 0.0, sigma2);
						weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma2);
						//count++;					
						//weight *= gau_mask[s][t];
						sum += weight * tmp[x][y];
						//TRACE("k = %d\n", k);
						//TRACE("gau[%d] = %f\n", k, gau[k]);
						w_sum += weight;
				}
				//TRACE("i = %d, j = %d\n", i, j);
				//TRACE("g = %f\n", g);
				//tmp2 = g;
				////////////////////////////////////
				// Normalization
				sum /= w_sum;
				//image[i][j] = round(tmp2);
				gray[i][j] = round(sum);
				/////////////////////////////////////
				//TRACE("N = %d\n", N);
				//TRACE("w_sum = %f\n", w_sum);
				//gau_mask.print();
				//if (image[i][j] > 255) image[i][j] = 255;
				//image[i][j] = (int)y;
			}
		}
		
		//LinearInterp(cmap2, cmap, 1.5, cmap);
		//LinearInterp3(cmap2, cmap, 1.5, cmap);

		tmp.copy(gray);
		//LinearInterp2(cmap2, cmap, z, cmap);
	}
	//z = 2.0;
	//LinearInterp(cmap2, cmap, z, cmap);
	//TRACE("max_g = %f\n", max_g);
	//TRACE("min_g = %f\n", min_g);
}

void GetPSNR(CDC& dc, cimatrix& cmap, char* file1, char* file2)
{
	int	x, y;
	double r1, g1, b1, r2, g2, b2;
	double MSE, PSNR;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	LoadBMP(dc, (char *)LPCTSTR(file1), &memDC);

	MSE = 0.0;
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r1 = (double)cmap[x][y].r;
			g1 = (double)cmap[x][y].g;
			b1 = (double)cmap[x][y].b;
			r2 = (double)RGB_GETRED(memDC.GetPixel(x, image_y-1-y));
			g2 = (double)RGB_GETGREEN(memDC.GetPixel(x, image_y-1-y));
			b2 = (double)RGB_GETBLUE(memDC.GetPixel(x, image_y-1-y));
			MSE += (r1-r2) * (r1-r2);
			MSE += (g1-g2) * (g1-g2);
			MSE += (b1-b2) * (b1-b2);
		}
	}
	MSE = MSE / (double)( image_x * image_y ) / 3.0;
	PSNR = 20.0 * log10(255.0 / sqrt(MSE));
	
	TRACE("MSE1 = %f\n", MSE);
	TRACE("PSNR1 = %f\n", PSNR);
	
	LoadBMP(dc, (char *)LPCTSTR(file2), &memDC);

	MSE = 0.0;
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r1 = (double)cmap[x][y].r;
			g1 = (double)cmap[x][y].g;
			b1 = (double)cmap[x][y].b;
			r2 = (double)RGB_GETRED(memDC.GetPixel(x, image_y-1-y));
			g2 = (double)RGB_GETGREEN(memDC.GetPixel(x, image_y-1-y));
			b2 = (double)RGB_GETBLUE(memDC.GetPixel(x, image_y-1-y));
			MSE += (r1-r2) * (r1-r2);
			MSE += (g1-g2) * (g1-g2);
			MSE += (b1-b2) * (b1-b2);
		}
	}
	MSE = MSE / (double)( image_x * image_y ) / 3.0;
	PSNR = 20.0 * log10(255.0 / sqrt(MSE));
	
	TRACE("MSE2 = %f\n", MSE);
	TRACE("PSNR2 = %f\n\n", PSNR);
}

void GetDiff(CDC& dc, imatrix& gray, char* file1, char* file2)
// Get per pixel color difference value!
{
	int	x, y;
	double r1, g1, b1, r2, g2, b2;
	double MSE;
	int exist;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	cimatrix cmap(image_x, image_y);
	cimatrix cmap2(image_x, image_y);

	exist = LoadBMP(dc, (char *)LPCTSTR(file1), &memDC);
	if (exist == 0) { TRACE("file1 does not exist!\n"); exit(1); }

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			cmap[x][y].r = (GLubyte)RGB_GETRED(memDC.GetPixel(x, image_y-1-y));
			cmap[x][y].g = (GLubyte)RGB_GETGREEN(memDC.GetPixel(x, image_y-1-y));
			cmap[x][y].b = (GLubyte)RGB_GETBLUE(memDC.GetPixel(x, image_y-1-y));
		}
	}

	exist = LoadBMP(dc, (char *)LPCTSTR(file2), &memDC);
	if (exist == 0) { TRACE("file2 does not exist!\n"); exit(1); }

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			cmap2[x][y].r = (GLubyte)RGB_GETRED(memDC.GetPixel(x, image_y-1-y));
			cmap2[x][y].g = (GLubyte)RGB_GETGREEN(memDC.GetPixel(x, image_y-1-y));
			cmap2[x][y].b = (GLubyte)RGB_GETBLUE(memDC.GetPixel(x, image_y-1-y));
		}
	}

	TRACE("image_x = %d\n", image_x);
	TRACE("image_y = %d\n", image_y);

	MSE = 0.0;
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r1 = (double)cmap[x][y].r / 256.;
			g1 = (double)cmap[x][y].g / 256.;
			b1 = (double)cmap[x][y].b  / 256.;
			r2 = (double)cmap2[x][y].r / 256.;
			g2 = (double)cmap2[x][y].g / 256.;
			b2 = (double)cmap2[x][y].b / 256.;
			MSE += dist3(r1-r2, g1-g2, b1-b2) / 1.732;
			//TRACE("MSE = %f\n", MSE);
		}
	}
	
	MSE /= (image_x * image_y);
	TRACE("per pixel diff = %f\n", MSE);

}

void GetPSNR_partial(CDC& dc, cimatrix& cmap, char* file1, char* file2)
{
	int	x, y;
	double r1, g1, b1, r2, g2, b2;
	double MSE, PSNR;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	int count;

	LoadBMP(dc, (char *)LPCTSTR(file1), &memDC);

	MSE = 0.0;
	count = 0;
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r1 = (double)cmap[x][y].r;
			g1 = (double)cmap[x][y].g;
			b1 = (double)cmap[x][y].b;
			//////////////////////////////////////////
			// red area is ignored!
			if (r1 == 255 && g1 == 0 && b1 == 0) continue;
			////////////////////////////////////////
			count++;
			//////////
			r2 = (double)RGB_GETRED(memDC.GetPixel(x, image_y-1-y));
			g2 = (double)RGB_GETGREEN(memDC.GetPixel(x, image_y-1-y));
			b2 = (double)RGB_GETBLUE(memDC.GetPixel(x, image_y-1-y));
			MSE += (r1-r2) * (r1-r2);
			MSE += (g1-g2) * (g1-g2);
			MSE += (b1-b2) * (b1-b2);
		}
	}
	//MSE = MSE / (double)( image_x * image_y ) / 3.0;
	TRACE("count1 = %d\n", count);
	MSE = MSE / (double)( count ) / 3.0;
	PSNR = 20.0 * log10(255.0 / sqrt(MSE));
	
	TRACE("MSE1 = %f\n", MSE);
	TRACE("PSNR1 = %f\n", PSNR);
	
	LoadBMP(dc, (char *)LPCTSTR(file2), &memDC);

	MSE = 0.0;
	count = 0;
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r1 = (double)cmap[x][y].r;
			g1 = (double)cmap[x][y].g;
			b1 = (double)cmap[x][y].b;
			//////////////////////////////////////////
			// red area is ignored!
			if (r1 == 255 && g1 == 0 && b1 == 0) continue;
			////////////////////////////////////////
			count++;
			//////////
			r2 = (double)RGB_GETRED(memDC.GetPixel(x, image_y-1-y));
			g2 = (double)RGB_GETGREEN(memDC.GetPixel(x, image_y-1-y));
			b2 = (double)RGB_GETBLUE(memDC.GetPixel(x, image_y-1-y));
			MSE += (r1-r2) * (r1-r2);
			MSE += (g1-g2) * (g1-g2);
			MSE += (b1-b2) * (b1-b2);
		}
	}
	//MSE = MSE / (double)( image_x * image_y ) / 3.0;
	TRACE("count2 = %d\n", count);
	MSE = MSE / (double)( count ) / 3.0;
	PSNR = 20.0 * log10(255.0 / sqrt(MSE));
	
	TRACE("MSE2 = %f\n", MSE);
	TRACE("PSNR2 = %f\n", PSNR);
}

void MarkEdgePixels(int x, int y, imatrix& tmp, Field& gfield)
// following the flow, mark the pixels along the flow
{
	double d_i, d_j, tx, ty, nx, ny;
	int int_i, int_j;
	double t;

	//pnts.clear();
	
	//sum = 0; 
	//count = 1;
	//sum += (int)RGB_GETRED(double_buffer.GetPixel(x, IMAGE_Y-1-y));
	////////////////////////////////////
	t = 1.0;
	////////////////////////////////////////////////
	// One half
	d_i = (double)x; d_j = (double)y; 
	tx = -gfield[x][y].gy;
	ty = gfield[x][y].gx;
	int_i = x;	int_j = y;
	while (1) {
		if (tx == 0.0 && ty == 0.0) break;
		///////////////////////////////////////
		nx = tx / sqrt(tx*tx+ty*ty);  // x component of the unit direction vector 
		ny = ty / sqrt(tx*tx+ty*ty);  // y component of the unit direction vector
		////////////////////////////////
		//////////////////////////////
		d_i += nx * t; // accurate new location x
		d_j += ny * t; // accurate new location y
        ///////////////////////////
		if (round(d_i) == int_i && round(d_j) == int_j) // no change
			continue; // push some more
		/////////////////////////////////////////
		int_i = round(d_i); // integer version of new location x
		int_j = round(d_j); // integer version of new location y
		if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		//////////////////////////////////////////////////
		if (tmp[int_i][int_j] == 0) // already marked!
			break;
		tmp[int_i][int_j] = 0; // mark it!
		/////////////////////////////////
		//d_i = (double)int_i;
		//d_j = (double)int_j;
		//////////////////////////////////
		tx = -gfield[int_i][int_j].gy;
		ty = gfield[int_i][int_j].gx;
	}
	///////////////////////////////////////
	// other half
	///*
	d_i = (double)x; d_j = (double)y; 
	int_i = x;	int_j = y;
	tx = -gfield[x][y].gy;
	ty = gfield[x][y].gx;
	while (1) {
		if (tx == 0.0 && ty == 0.0) break;
		///////////////////////////////////////
		nx = tx / sqrt(tx*tx+ty*ty);  // x component of the direction vector 
		ny = ty / sqrt(tx*tx+ty*ty);  // y component of the direction vector
		////////////////////////////////
		//////////////////////////////
		d_i -= nx * t; // accurate new location x
		d_j -= ny * t; // accurate new location y
		///////////////////////////
		if (round(d_i) == int_i && round(d_j) == int_j) // no change
			continue; // push some more
		//////////////////////////////////////
		int_i = round(d_i); // integer version of new location x
		int_j = round(d_j); // integer version of new location y
		if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		////////////////////////////////////////////////
		if (tmp[int_i][int_j] == 0) // already marked!
			break;
		tmp[int_i][int_j] = 0; // mark it!
		/////////////////////////////////
		//d_i = (double)int_i;
		//d_j = (double)int_j;
		/////////////////////////////////
		tx = -gfield[int_i][int_j].gy;
		ty = gfield[int_i][int_j].gx;
	}
	//*/
		
	
}

		
void FlowConnectEdges(int image_x, int image_y, imatrix& image, Field& gfield)
// Connect edges after Canny's edge detection, based on Gfield
{
	//int	dx, dy;
	int	i, j;
	
	imatrix tmp(IMAGE_X, IMAGE_Y);
	tmp.copy(image);

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			//gx = gy = 0.0;
			//w_sum = 0.0;
			if (image[i][j] == 0) {
				MarkEdgePixels(i, j, tmp, gfield);
			}
		}
	}

	image.copy(tmp);
	////////////////////////////////////////////
}





void CopyMembuffer(int x, int y, int half, GLubyte* Membuffer, GLubyte* Doublebuffer)
{
	int	i, j;
	//int MAX_GRADIENT = -1;
    //int gau_w;
	int l, r, t, b;
	//GLubyte z;
	//Image tmp;
	//Image target;

	//Image tmp;
	l = x-half;
	if (l < 0) l = 0;
	r = x+half;
	if (r > IMAGE_X-1) r = IMAGE_X-1;
	t = y+half;
	if (t > IMAGE_Y-1) t = IMAGE_Y-1;
	b = y-half;
	if (b < 0) b = 0;

	/// copy memDC to tmp
	for (j = b; j <= t; j++) {
		for (i = l; i <= r; i++) {
			Doublebuffer[j * IMAGE_X + i] = Membuffer[j * IMAGE_X + i];
			//TRACE("tmp[%d][%d] = %d\n", i, j, tmp[i][j]);
		}
	}
	
}

void CopyCmap2Membuffer(cimatrix& cmap, GLubyte* Dbuffer)
{
	int	i, j;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			Dbuffer[(j * image_x + i) * 3 + 0] = cmap[i][j].r;
			Dbuffer[(j * image_x + i) * 3 + 1] = cmap[i][j].g;
			Dbuffer[(j * image_x + i) * 3 + 2] = cmap[i][j].b;
		}
	}
	
}

void CopyGray2Membuffer(imatrix& gray, GLubyte* Dbuffer)
{
	int	i, j;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			Dbuffer[(j * image_x + i) * 3 + 0] = gray[i][j];
			Dbuffer[(j * image_x + i) * 3 + 1] = gray[i][j];
			Dbuffer[(j * image_x + i) * 3 + 2] = gray[i][j];
		}
	}
	
}





void OnGlobalCanny(CDC& dc, int image_x, int image_y, imatrix& image, imatrix& image2)
{
	if (file_loaded) {
		//CClientDC dc(this);
		//GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);
		gau_w = MakeGaussMask(cur_sigma, gau);
		//cur_sigma = 1.0;
		//dlg->m_Scale_Scroll.SetScrollPos((int)(cur_sigma*10));
		//dlg->m_slider.SetScrollPos(SB_HORZ, (int)(cur_sigma*10), TRUE);
		//dlg->m_scale_val;
		//GaussSmooth(IMAGE_X, IMAGE_Y, image, gau_w);
		max_grad2 = GlobalCanny(IMAGE_X, IMAGE_Y, gray, gray2, gau_w);
		//NonmaxSuppressGray(IMAGE_X, IMAGE_Y);
		DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray2);
		//////////////////////////////
		emap.copy(gray2);
		////////////////////////////////
		/*
		//////////////////////////////////////////////////////
		///////////////////////////////////////////////////
		*/
		//DrawEdgeStrokes(memDC, IMAGE_X, IMAGE_Y, gray2);
        
		dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);

	}
}

double GetCannyGradient(int image_x, int image_y, imatrix& image, int gau_w)
{
	int	i, j, k, i1, i2, j1, j2;
	int MAX_GRADIENT = -1;

	//image2.init(image_x, image_y);
	
	imatrix tmp_x(image_x, image_y);
	imatrix tmp_y(image_x, image_y);
	G_x.init(image_x, image_y);
	G_y.init(image_x, image_y);
	G_mag.init(image_x, image_y);

	double	x, y;
	
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			x = gau[0] * image[i][j];	
			y = gau[0] * image[i][j];
			//TRACE("x = %f, y = %f\n", x, y);
			for (k = 1; k < gau_w; k++) {
				i1 = (i+k)%image_x;
				i2 = (i-k+image_x)%image_x;
				x += gau[k] * image[i1][j] + gau[k] * image[i2][j];
				j1 = (j+k)%image_y;
				j2 = (j-k+image_y)%image_y;
				y += gau[k] * image[i][j1] + gau[k] * image[i][j2];
			}
			tmp_x[i][j] = x;
			tmp_y[i][j] = y;
			//TRACE("x = %f, y = %f\n", x, y);
			if (x > 255) x = 255;
			if (y > 255) y = 255;
			//image[i][j] = (int)x;
			//image[i][j] = (int)y;
		}
	}

	for (j = 1; j < image_y - 1; j++) {
		for (i = 1; i < image_x - 1; i++) {
			G_x[i][j] = (tmp_x[i+1][j-1] + 2*tmp_x[i+1][j] + tmp_x[i+1][j+1] 
				- tmp_x[i-1][j-1] - 2*tmp_x[i-1][j] - tmp_x[i-1][j+1]);
			G_y[i][j] = (tmp_y[i-1][j+1] + 2*tmp_y[i][j+1] + tmp_y[i+1][j+1]
				- tmp_y[i-1][j-1] - 2*tmp_y[i][j-1] - tmp_y[i+1][j-1]);
			//G_mag[i][j] = sqrt(G_x[i][j] * G_x[i][j] + G_y[i][j] * G_y[i][j]);
			G_mag[i][j] = norm2(G_x[i][j], G_y[i][j]);

			if (G_mag[i][j] > MAX_GRADIENT) {
				MAX_GRADIENT = round(G_mag[i][j]);
				//TRACE("p[i][j].gx = %01.f\n", p[i][j].gx);
				//TRACE("p[i][j].gy = %01.f\n", p[i][j].gy);
				//TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);
			}
		}
	}

	return MAX_GRADIENT;
}


void hyster_visit(int image_x, int image_y, int x0, int y0, imatrix& image2)
{
	int	i, j, x, y;
	int	done_flag;
	double g;

	//image2[x0][y0] = 0; // now marked as edge
	done_flag = 0;
	for (j = -1; j <= 1; j++) {
		for (i = -1; i <= 1; i++) {
			x = x0 + i;
			y = y0 + j;
			if (x <= 0 || x >= image_x-1 || y <= 0 || y >= image_y-1)
				continue;
			if (image2[x][y] == 0) continue; // already marked as edge
			if (!thin_edge[x][y]) continue; // not a thin edge
			g = G_mag[x][y];
			if (g > lo_thres) {
				image2[x][y] = 0; // newly marked as edge pixel
				hyster_visit(image_x, image_y, x, y, image2);
				//done_flag = 1; // found an extending edge pixel
				//break;
			}
		}
		//if (done_flag) break; // found the right neighboring edge pixel. Get out!
	}
}



double prob(double r, double g, double b, double gr, edge_dist& s)
{
	double p_r, p_g, p_b, p_gr;
	double p_val1, p_val2, p_val;

	p_r = gauss2(r, s.r_mean, s.r_std);
	p_g = gauss2(g, s.g_mean, s.g_std);
	p_b = gauss2(b, s.b_mean, s.b_std);
	p_gr = gauss2(gr, s.gr_mean, s.gr_std);

	if (s.r_mean == 0.0 && s.g_mean == 0.0 && s.b_mean == 0.0)
		return 0.0;
	else {
		/*
		TRACE("\n");
		TRACE("r = %f\n", r);
		TRACE("g = %f\n", g);
		TRACE("b = %f\n", b);
		TRACE("gr = %f\n", gr);
		TRACE("r_mean = %f\n", s.r_mean);
		TRACE("g_mean = %f\n", s.g_mean);
		TRACE("b_mean = %f\n", s.b_mean);
		TRACE("gr_mean = %f\n", s.gr_mean);
		TRACE("r_std = %f\n", s.r_std);
		TRACE("g_std = %f\n", s.g_std);
		TRACE("b_std = %f\n", s.b_std);
		TRACE("gr_std = %f\n", s.gr_std);
		TRACE("p_r = %f\n", p_r);
		TRACE("p_g = %f\n", p_g);
		TRACE("p_b = %f\n", p_b);
		TRACE("p_gr = %f\n", p_gr);
		*/
		//p_val = (p_r + p_g + p_b + p_gr) / 4.0; 
		p_val1 = (p_r + p_g + p_b) / 3.0; 
		p_val2 = p_gr;
		p_val = 0.6 * p_val1 + 0.4 * p_val2;
		//TRACE("\ngr = %f\n", gr);
		//TRACE("gr_mean = %f\n", s.gr_mean);
		//TRACE("gr_std = %f\n", s.gr_std);
		//TRACE("p_val = %f\n", p_val);
        
		return p_val;
	}
}

bool prob2(double r, double g, double b, double gr, edge_dist& s, double factor)
{
	//double factor;

	//factor = 3.5;

	//p_r = gauss2(r, s.r_mean, s.r_std);
	//p_g = gauss2(g, s.g_mean, s.g_std);
	//p_b = gauss2(b, s.b_mean, s.b_std);
	//p_gr = gauss2(gr, s.gr_mean, s.gr_std);

	if (s.r_mean == 0.0 && s.g_mean == 0.0 && s.b_mean == 0.0)
		return false;
	else {
		if ( ABS(r-s.r_mean) <= factor * s.r_std 
			&& ABS(g-s.g_mean) <= factor * s.g_std 
			&& ABS(b-s.b_mean) <= factor * s.b_std 
			&& ABS(gr-s.gr_mean) <= factor * s.gr_std)
			return true;
		else return false;
		/*
		TRACE("\n");
		TRACE("r = %f\n", r);
		TRACE("g = %f\n", g);
		TRACE("b = %f\n", b);
		TRACE("gr = %f\n", gr);
		TRACE("r_mean = %f\n", s.r_mean);
		TRACE("g_mean = %f\n", s.g_mean);
		TRACE("b_mean = %f\n", s.b_mean);
		TRACE("gr_mean = %f\n", s.gr_mean);
		TRACE("r_std = %f\n", s.r_std);
		TRACE("g_std = %f\n", s.g_std);
		TRACE("b_std = %f\n", s.b_std);
		TRACE("gr_std = %f\n", s.gr_std);
		*/
	}
}

void off_edge_visit(int image_x, int image_y, int x0, int y0)
{
	int	i, j, x, y;
	double g;
	GLubyte rr, gg, bb;

	off_mark[x0][y0] = 255; // newly marked as edge pixel to be removed
	gray2[x0][y0] = 255; // now marked as non-edge
	for (j = -1; j <= 1; j++) {
		for (i = -1; i <= 1; i++) {
			x = x0 + i;
			y = y0 + j;
			if (x <= 0 || x >= image_x-1 || y <= 0 || y >= image_y-1)
				continue;
			if (off_mark[x][y] == 255) continue; // already marked as edge to be removed
			if (gray2[x][y] == 255) continue; // already marked as non-edge
			//if (!thin_edge[x][y]) continue; // not a thin edge
			//g = G_mag[x][y];
			//if (image2[x][y] == 0 && off_mark[x][y] == 0) { // marked as edge
			rr = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
			gg = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
			bb = Dbuffer[(y * IMAGE_X + x) * 3 + 2];
			g = G_mag[x][y];
			if ( prob2((double)rr, (double)gg, (double)bb, g, off_edge, factor2) ) {
				off_edge_visit(image_x, image_y, x, y);
			}
		}
		//if (done_flag) break; // found the right neighboring edge pixel. Get out!
	}
}

void ConvertThinEdgeMap2GrayImage(int image_x, int image_y, imatrix& thin_edge, imatrix& image)
{
	int i, j;

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			//TRACE("thin_edge[%d][%d] = %.3f\n", i, j, thin_edge[i][j]);
			//TRACE("thin_edge[%d][%d] = %d\n", i, j, thin_edge[i][j]);
			if (thin_edge[i][j]) {
				//TRACE("I'm thin edge!\n");
				image[i][j] = 0;
			}
			else image[i][j] = 255;
		}
	}
}

double GlobalCanny(int image_x, int image_y, imatrix& image, imatrix& image2, int gau_w)
{
	int	i, j, k, i1, i2, j1, j2;
	double MAX_GRADIENT = -1;
	GLubyte rr, gg, bb;

	//imatrix marked(image_x, image_y);
	//image2.init(image_x, image_y);
	
	//imatrix tmp(image_x, image_y);
	double	x, y;
	
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			x = gau[0] * image[i][j];	
			y = gau[0] * image[i][j];
			//TRACE("x = %f, y = %f\n", x, y);
			for (k = 1; k < gau_w; k++) {
				//i1 = (i+k)%image_x;
				//i2 = (i-k+image_x)%image_x;
				//j1 = (j+k)%image_y;
				//j2 = (j-k+image_y)%image_y;
				i1 = i+k;
				i2 = i-k;
				j1 = j+k;
				j2 = j-k;
				if (i1 > image_x-1) i1 = image_x-1;
				if (i2 < 0) i2 = 0;
				if (j1 > image_y-1) j1 = image_y-1;
				if (j2 < 0) j2 = 0;
				x += gau[k] * image[i1][j] + gau[k] * image[i2][j];
				y += gau[k] * image[i][j1] + gau[k] * image[i][j2];
			}
			tmp_x[i][j] = x;
			tmp_y[i][j] = y;
			//TRACE("x = %f, y = %f\n", x, y);
			if (x > 255) x = 255;
			if (y > 255) y = 255;
			//image[i][j] = (int)x;
			//image[i][j] = (int)y;
		}
	}

	for (j = 1; j < image_y - 1; j++) {
		for (i = 1; i < image_x - 1; i++) {
			G_x[i][j] = (tmp_x[i+1][j-1] + 2*tmp_x[i+1][j] + tmp_x[i+1][j+1] 
				- tmp_x[i-1][j-1] - 2*tmp_x[i-1][j] - tmp_x[i-1][j+1]);
			G_y[i][j] = (tmp_y[i-1][j+1] + 2*tmp_y[i][j+1] + tmp_y[i+1][j+1]
				- tmp_y[i-1][j-1] - 2*tmp_y[i][j-1] - tmp_y[i+1][j-1]);
			//G_mag[i][j] = sqrt(G_x[i][j] * G_x[i][j] + G_y[i][j] * G_y[i][j]);
			G_mag[i][j] = norm2(G_x[i][j], G_y[i][j]);

			/*
			if (G_mag[i][j] > MAX_GRADIENT) {
				MAX_GRADIENT = G_mag[i][j];
				//TRACE("p[i][j].gx = %01.f\n", p[i][j].gx);
				//TRACE("p[i][j].gy = %01.f\n", p[i][j].gy);
				//TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);
			}
			*/
		}
	}

	// Normalize each gradient value & init marked image
	//TRACE("MAX_GRADIENT = %f\n", MAX_GRADIENT);
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			if (i == 0 || i == image_x-1 || j == 0 || j == image_y-1) {
				image2[i][j] = 255;
				thin_edge[i][j] = 0; // init thin edge list
				continue;
			}
			G_mag[i][j] = (G_mag[i][j] / max_grad2); // G_mag between [0, 1]
			image2[i][j] = (int)(G_mag[i][j] * 255);

			//marked[i][j] = 0; // init marked image for hysteresis
			thin_edge[i][j] = 0; // init thin edge list
			scale_map[i][j] = cur_sigma; // for adaptive painterly rendering
		}
	}

	double	gx, gy;
	double	g, g1, g2, g3, g4;
	double	t; // interpolation parameter
	
	//////////////////////////////////////////
	// Nonmaxima suppression
	for (j = 1; j < image_y-1; j++) {
		for (i = 1; i < image_x-1; i++) {
			gx = G_x[i][j];
			gy = G_y[i][j];
			g = G_mag[i][j];
			//TRACE("gx = %f, gy = %f, g = %f\n", gx, gy, g);
			//if (gx < 0.01 && gy < 0.01) continue; // not an edge
			if (fabs(gx) >= fabs(gy)) { // close to horizontal (note: gy can be 0)
				t = fabs(gy) / fabs(gx); 
				//TRACE("t = %f\n", t);
				g1 = G_mag[i+1][j]; // right
				g2 = G_mag[i-1][j]; // left
				//TRACE("g1 = %f, g2 = %f\n", g1, g2);
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // right up
					g4 = G_mag[i-1][j-1]; // left down
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i+1][j-1]; // right down
					g4 = G_mag[i-1][j+1]; // left up
				}
			}
			else { // close to vertical (note: gx can be 0)
				t = fabs(gx) / fabs(gy);
				g1 = G_mag[i][j+1]; // up
				g2 = G_mag[i][j-1]; // down
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // up right
					g4 = G_mag[i-1][j-1]; // down left
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i-1][j+1]; // up left
					g4 = G_mag[i+1][j-1]; // down right
				}
			}
			if ( g > ((1-t)*g1 + t*g3) && g > ((1-t)*g2 + t*g4) ) {
				////////////////////////////////////////
				thin_edge[i][j] = 1; // it's a thin edge
				//////////////////////////////////////
				if (g > hi_thres) {
					
					rr = Dbuffer[(j * IMAGE_X + i) * 3 + 0];
					gg = Dbuffer[(j * IMAGE_X + i) * 3 + 1];
					bb = Dbuffer[(j * IMAGE_X + i) * 3 + 2];
					//p_val = prob((double)rr, (double)gg, (double)bb, g, off_edge);
					//TRACE("p_val = %f\n", p_val);
					/*
					if ( prob2((double)rr, (double)gg, (double)bb, g, off_edge) )
						image2[i][j] = 255; // off edge
					else 
						image2[i][j] = 0; // thin edge above hi_thres
					*/
					///*
					image2[i][j] = 0; // thin edge above hi_thres
					if ( prob2((double)rr, (double)gg, (double)bb, g, off_edge, factor1) )
						off_mark[i][j] = 255; // off edge
					//*/
					/*
					if ( p_val > 0.4 )
						image2[i][j] = 255; // off edge
					else 
						image2[i][j] = 0; // thin edge above hi_thres
					*/
				}
				else 
					image2[i][j] = 255; // thin edge below hi_thres
			}
			else { // non-maximum
				image2[i][j] = 255;
			}
			
		}
	}

	//////////////////////////////////////////
	// Hysteresis thresholding
	///*
	for (j = 1; j < image_y-1; j++) {
		for (i = 1; i < image_x-1; i++) {
			if (image2[i][j] == 0) { // computed thinned edges above hi_thres
				hyster_visit(image_x, image_y, i, j, image2); // visit neighboring pixels
			}
		}
	}
	//*/
	//////////////////////////////////////////
	// removing off_edges
	///*
	for (j = 1; j < image_y-1; j++) {
		for (i = 1; i < image_x-1; i++) {
			if (off_mark[i][j] == 255) { // marked off edge
				off_edge_visit(image_x, image_y, i, j); // visit neighboring pixels
			}
		}
	}
	//*/
	///////////////////////////////////////
	// clear the image border
	/*
	for (i = 0; i < image_x; i++) {
		image2[i][0] == 255;
		image2[i][image_y-1] == 255;
	}
	for (j = 0; j < image_y; j++) {
		image2[0][j] == 255;
		image2[image_x-1][j] == 255;
	}
	*/
	GLOBAL_CANNY_DONE = 1;
	
	return max_grad2;
}

void local_hyster_visit(int size, PixeL seed, int image_x, int image_y, int x0, int y0, imatrix& image2)
{
	int	i, j, x, y;
	int	done_flag;
	double g;
	int half;
	
	half = size/2;

	//image2[x0][y0] = 0; // now marked as edge
	done_flag = 0;
	for (j = -1; j <= 1; j++) {
		for (i = -1; i <= 1; i++) {
			x = x0 + i;
			y = y0 + j;
			if (x <= 0 || x >= image_x-1 || y <= 0 || y >= image_y-1)
				continue;
			if (image2[x][y] == 0) continue; // already marked as edge
			if (!thin_edge[x][y]) continue; // not a thin edge
			/////////////////////////////////////////////////////////
			// circular kernel
			if ( dist2(x, y, seed.x, seed.y) > half-1 ) continue; 
			//////////////////////////////////////////////////
			g = G_mag[x][y];
			if (g > lo_thres) {
				image2[x][y] = 0; // newly marked as edge pixel
				local_hyster_visit(size, seed, image_x, image_y, x, y, image2);
				//done_flag = 1; // found an extending edge pixel
				//break;
			}
		}
		//if (done_flag) break; // found the right neighboring edge pixel. Get out!
	}
}

edge_dist on_edge, off_edge, cur_edge;
pixel_dist pix_d[10];
int pix_d_count = -1;

void EdgeColor(CDC& dc, imatrix& gray2, int x0, int y0, GLubyte r, GLubyte g, GLubyte b)
{
	int i, j, x, y;
	GLubyte r1, g1, b1, r2, g2, b2;
	double a = 0.3; // learning rate
	double gr;

	dc.SetPixelV(x0, (IMAGE_Y-1)-y0, RGB(r, g, b));

	r1 = Dbuffer[(y0 * IMAGE_X + x0) * 3 + 0];
	g1 = Dbuffer[(y0 * IMAGE_X + x0) * 3 + 1];
	b1 = Dbuffer[(y0 * IMAGE_X + x0) * 3 + 2];
	gr = G_mag[x0][y0];

	//TRACE("gr = %.5f\n", gr);

	off_mark[x0][y0] = 1; // marked as sampled edge pixel
	//cur_edge.update((double)r1, (double)g1, (double)b1, gr);
	/*
	cur_edge.r_mean = (1-a) * cur_edge.r_mean + a * r1;
	cur_edge.g_mean = (1-a) * cur_edge.g_mean + a * g1;
	cur_edge.b_mean = (1-a) * cur_edge.b_mean + a * b1;
	cur_edge.gr_mean = (1-a) * cur_edge.gr_mean + a * gr;
	cur_edge.r_std = (1-a) * cur_edge.r_std + a * (r1 - cur_edge.r_mean);
	cur_edge.g_std = (1-a) * cur_edge.g_std + a * (g1 - cur_edge.g_mean);
	cur_edge.b_std = (1-a) * cur_edge.b_std + a * (b1 - cur_edge.b_mean);
	cur_edge.gr_std = (1-a) * cur_edge.gr_std + a * (gr - cur_edge.gr_mean);
	*/

	for (j = -1; j <= 1; j++) {
		for (i = -1; i <= 1; i++) {
			x = x0 + i;
			y = y0 + j;
			if (x <= 0 || x >= IMAGE_X-1 || y <= 0 || y >= IMAGE_Y-1)
				continue;
			if (gray2[x][y] == 255) continue; // not an edge pixel
			r2 = (GLubyte)RGB_GETRED(dc.GetPixel(x, (IMAGE_Y-1)-y));
			g2 = (GLubyte)RGB_GETGREEN(dc.GetPixel(x, (IMAGE_Y-1)-y));
			b2 = (GLubyte)RGB_GETBLUE(dc.GetPixel(x, (IMAGE_Y-1)-y));
			if (r2 == r && g2 == g && b2 == b) continue; // already marked as selected
			//if (!thin_edge[x][y]) continue; // not a thin edge
			/////////////////////////////////////////////////////////
			// circular kernel
			//if ( dist2(x, y, seed.x, seed.y) > half-1 ) continue; 
			//////////////////////////////////////////////////
			if (gray2[x][y] == 0) { // edge pixel
				EdgeColor(dc, gray2, x, y, r, g, b);
				//done_flag = 1; // found an extending edge pixel
				//break;
			}
		}
		//if (done_flag) break; // found the right neighboring edge pixel. Get out!
	}
}


	
double LocalCanny(int image_x, int image_y, int size, PixeL seed, imatrix& image, 
					   imatrix& image2, int gau_w)
{
	int	i, j, k, i1, i2, j1, j2, s, r;
	//int MAX_GRADIENT = -1;
	int	half = size/2;

	//tmp_x.init(size, size);
	//tmp_y.init(size, size);
	//G_x.init(size, size);
	//G_y.init(size, size);
	//G_mag.init(size, size);

	//imatrix tmp(image_x, image_y);
	double	x, y;

	for (s = 0; s < size; s++) {
		for (r = 0; r < size; r++) {
			
			i = seed.x-half+s;
			j = seed.y-half+r;
			if (i < 0 || i > image_x-1 || j < 0 || j > image_y-1)
				continue;
			x = gau[0] * image[i][j];	
			y = gau[0] * image[i][j];
			for (k = 1; k < gau_w; k++) {
				i1 = (i+k)%image_x;
				i2 = (i-k+image_x)%image_x;
				x += gau[k] * image[i1][j] + gau[k] * image[i2][j];
				j1 = (j+k)%image_y;
				j2 = (j-k+image_y)%image_y;
				y += gau[k] * image[i][j1] + gau[k] * image[i][j2];
			}
			tmp_x[i][j] = x;
			tmp_y[i][j] = y;
			//TRACE("x = %f, y = %f\n", x, y);
			if (x > 255) x = 255;
			if (y > 255) y = 255;
	
		}
	}

	for (s = 1; s < size-1; s++) {
		for (r = 1; r < size-1; r++) {

			i = seed.x-half+s;
			j = seed.y-half+r;
			if (i <= 0 || i >= image_x-1 || j <= 0 || j >= image_y-1)
				continue;
			G_x[i][j] = (tmp_x[i+1][j-1] + 2*tmp_x[i+1][j] + tmp_x[i+1][j+1] 
				- tmp_x[i-1][j-1] - 2*tmp_x[i-1][j] - tmp_x[i-1][j+1]);
			G_y[i][j] = (tmp_y[i-1][j+1] + 2*tmp_y[i][j+1] + tmp_y[i+1][j+1]
				- tmp_y[i-1][j-1] - 2*tmp_y[i][j-1] - tmp_y[i+1][j-1]);
			//G_mag[i][j] = sqrt(G_x[i][j] * G_x[i][j] + G_y[i][j] * G_y[i][j]);
			G_mag[i][j] = norm2(G_x[i][j], G_y[i][j]);

			/*
			if (G_mag[i][j] > max_grad2) {
				max_grad2 = G_mag[i][j];
				//TRACE("p[i][j].gx = %01.f\n", p[i][j].gx);
				//TRACE("p[i][j].gy = %01.f\n", p[i][j].gy);
				//TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);
			}
			*/
		}
	}

	//// normalize each gradient value
	TRACE("max_grad2 = %f\n", max_grad2);
	TRACE("hi_thres = %f\n", hi_thres);
	for (s = 1; s < size-1; s++) {
		for (r = 1; r < size-1; r++) {
			i = seed.x-half+s;
			j = seed.y-half+r;
			if (i <= 0 || i >= image_x-1 || j <= 0 || j >= image_y-1)
				continue;
			G_mag[i][j] = (G_mag[i][j] / max_grad2); // G_mag between [0, 1]
			//image2[i][j] = (int)(G_mag[i][j] * 255);
			thin_edge[i][j] = 0; // init thin edge list
		}
	}

	double	gx, gy;
	double	g, g1, g2, g3, g4;
	double	t; // interpolation parameter

	// nonmaxima suppress
	for (s = 1; s < size-1; s++) {
		for (r = 1; r < size-1; r++) {
			/////////////////////////////////////////////////////////
			// circular kernel
			if ( dist2(s, r, half, half) > half-1 ) continue; 
			
			i = seed.x-half+s;
			j = seed.y-half+r;
			if (i <= 0 || i >= image_x-1 || j <= 0 || j >= image_y-1)
				continue;
			gx = G_x[i][j];
			gy = G_y[i][j];
			g = G_mag[i][j];
			//////////////////////////////////////////////////
			// mark pixel image
			pixel_mark[i][j] = 1; // marked
			scale_map[i][j] = cur_sigma;
			/////////////////////////////////////////////
			//TRACE("gx = %f, gy = %f, g = %f\n", gx, gy, g);
			//if (gx < 0.01 && gy < 0.01) continue; // not an edge
			if (fabs(gx) >= fabs(gy)) { // close to horizontal (note: gy can be 0)
				t = fabs(gy) / fabs(gx); 
				//TRACE("t = %f\n", t);
				g1 = G_mag[i+1][j]; // right
				g2 = G_mag[i-1][j]; // left
				//TRACE("g1 = %f, g2 = %f\n", g1, g2);
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // right up
					g4 = G_mag[i-1][j-1]; // left down
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i+1][j-1]; // right down
					g4 = G_mag[i-1][j+1]; // left up
				}
			}
			else { // close to vertical (note: gx can be 0)
				t = fabs(gx) / fabs(gy);
				g1 = G_mag[i][j+1]; // up
				g2 = G_mag[i][j-1]; // down
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // up right
					g4 = G_mag[i-1][j-1]; // down left
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i-1][j+1]; // up left
					g4 = G_mag[i+1][j-1]; // down right
				}
			}
			if ( g > ((1-t)*g1 + t*g3) && g > ((1-t)*g2 + t*g4) ) {
				thin_edge[i][j] = 1; // it's a thin edge
				if (g > hi_thres) 
					image2[i][j] = 0; // thin edge above hi_thres
				else 
					image2[i][j] = 255; // thin edge below hi_thres
			}
			else { // non-maximum
				//image[i][j] = 0;
				image2[i][j] = 255;
			}
			
		}
	}

	//////////////////////////////////////////
	// Hysteresis thresholding
	///*
	for (s = 0; s < size; s++) {
		for (r = 0; r < size; r++) {
			/////////////////////////////////////////////////////////
			// circular kernel
			//if ( dist2(s, r, half, half) > half-1 ) continue; 
			//////////////////////////////////////////////////
			i = seed.x-half+s;
			j = seed.y-half+r;
			if (i <= 0 || i >= image_x-1 || j <= 0 || j >= image_y-1)
				continue;
			if (image2[i][j] == 0) { // computed thinned edges above hi_thres
				local_hyster_visit(size, seed, image_x, image_y, i, j, image2); // visit neighboring pixels
			}
		}			
	}
	//*/

	return max_grad2;
		
}


void NonmaxSuppress(int image_x, int image_y)
{
	int	i, j;
	int MAX_GRADIENT = -1;

	double	gx, gy;
	double	g, g1, g2, g3, g4;
	double	t; // interpolation parameter
	
	for (j = 1; j < image_y-1; j++) {
		for (i = 1; i < image_x-1; i++) {
			gx = G_x[i][j];
			gy = G_y[i][j];
			g = G_mag[i][j];
			//TRACE("gx = %f, gy = %f, g = %f\n", gx, gy, g);
			//if (gx < 0.01 && gy < 0.01) continue; // not an edge
			if (fabs(gx) >= fabs(gy)) { // close to horizontal (note: gy can be 0)
				t = fabs(gy) / fabs(gx); 
				//TRACE("t = %f\n", t);
				g1 = G_mag[i+1][j]; // right
				g2 = G_mag[i-1][j]; // left
				//TRACE("g1 = %f, g2 = %f\n", g1, g2);
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // right up
					g4 = G_mag[i-1][j-1]; // left down
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i+1][j-1]; // right down
					g4 = G_mag[i-1][j+1]; // left up
				}
			}
			else { // close to vertical (note: gx can be 0)
				t = fabs(gx) / fabs(gy);
				g1 = G_mag[i][j+1]; // up
				g2 = G_mag[i][j-1]; // down
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // up right
					g4 = G_mag[i-1][j-1]; // down left
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i-1][j+1]; // up left
					g4 = G_mag[i+1][j-1]; // down right
				}
			}
			if ( g > ((1-t)*g1 + t*g3) && g > ((1-t)*g2 + t*g4) && g > 0.1 ) {
				//image[i][j] = G_mag[i][j]*255;
				//if (g > 0.5)
				//image[i][j] = 255;
				image[i][j] = 0;

			}
			else { // non-maximum
				//image[i][j] = 0;
				image[i][j] = 255;
			}
			
		}
	}
}

void NonmaxSuppressGray(int image_x, int image_y)
{
	int	i, j;
	int MAX_GRADIENT = -1;

	double	gx, gy;
	double	g, g1, g2, g3, g4;
	double	t; // interpolation parameter
	
	for (j = 1; j < image_y-1; j++) {
		for (i = 1; i < image_x-1; i++) {
			gx = G_x[i][j];
			gy = G_y[i][j];
			g = G_mag[i][j];
			//TRACE("gx = %f, gy = %f, g = %f\n", gx, gy, g);
			//if (gx < 0.01 && gy < 0.01) continue; // not an edge
			if (fabs(gx) >= fabs(gy)) { // close to horizontal (note: gy can be 0)
				t = fabs(gy) / fabs(gx); 
				//TRACE("t = %f\n", t);
				g1 = G_mag[i+1][j]; // right
				g2 = G_mag[i-1][j]; // left
				//TRACE("g1 = %f, g2 = %f\n", g1, g2);
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // right up
					g4 = G_mag[i-1][j-1]; // left down
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i+1][j-1]; // right down
					g4 = G_mag[i-1][j+1]; // left up
				}
			}
			else { // close to vertical (note: gx can be 0)
				t = fabs(gx) / fabs(gy);
				g1 = G_mag[i][j+1]; // up
				g2 = G_mag[i][j-1]; // down
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // up right
					g4 = G_mag[i-1][j-1]; // down left
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i-1][j+1]; // up left
					g4 = G_mag[i+1][j-1]; // down right
				}
			}
			if ( g > ((1-t)*g1 + t*g3) && g > ((1-t)*g2 + t*g4) && g > 0.1 ) {
				//image[i][j] = G_mag[i][j]*255;
				//if (g > 0.5)
				//image[i][j] = 255;
				gray[i][j] = 0;

			}
			else { // non-maximum
				//image[i][j] = 0;
				gray[i][j] = 255;
			}
			
			
		}
	}
		
		
}

void GetSobelGradient(int image_x, int image_y, imatrix& image, imatrix& image2, matrix& G_mag) 
// Sobel gradient version
// makes edge lines black (not white)
{
	int i, j;
	double MAX_GRADIENT = -1.;
	//double MAX_VAL = -1.;
	//double MAX_VAL = 255.;
	double MAX_VAL = 1020.; // 255 + 2 * 255 + 255 (Maximum possible Sobel value)

	Field tmp;

	tmp.init(image_x, image_y);

	/*
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			if (image[i][j] > MAX_VAL) MAX_VAL = (double)image[i][j];
		}
	}
	TRACE("MAX_VAL = %f\n", MAX_VAL);
	*/
	
	for (i = 1; i < image_x - 1; i++) { 
		for (j = 1; j < image_y - 1; j++) {
			////////////////////////////////////////////////////////////////
			// Important!: the value of image intensity should be normalized to [0,1]
			tmp[i][j].gx = (image[i+1][j-1] + 2*(double)image[i+1][j] + image[i+1][j+1] 
				- image[i-1][j-1] - 2*(double)image[i-1][j] - image[i-1][j+1]) / MAX_VAL;
			tmp[i][j].gy = (image[i-1][j+1] + 2*(double)image[i][j+1] + image[i+1][j+1]
				- image[i-1][j-1] - 2*(double)image[i][j-1] - image[i+1][j-1]) / MAX_VAL;
			/////////////////////////////////////////////
			// Amplify!!!
			//p[i][j].gx = pow(p[i][j].gx, 2); 
			//p[i][j].gy = pow(p[i][j].gy, 2); 
			//p[i][j].gx = pow(p[i][j].gx, 2); 
			//TRACE("p[i][j].gx = %.1f\n", p[i][j].gx);
			//TRACE("p[i][j].gy = %.1f\n", p[i][j].gy);
			tmp[i][j].mag = sqrt(tmp[i][j].gx * tmp[i][j].gx + tmp[i][j].gy * tmp[i][j].gy);

			if (tmp[i][j].mag > MAX_GRADIENT) {
				MAX_GRADIENT = tmp[i][j].mag;
				//TRACE("p[i][j].gx = %01.f\n", p[i][j].gx);
				//TRACE("p[i][j].gy = %01.f\n", p[i][j].gy);
				//TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);
			}
		}
	}

	for (i = 1; i <= image_x - 2; i++) {
		tmp[i][0].gx = tmp[i][1].gx;
		tmp[i][0].gy = tmp[i][1].gy;
		tmp[i][0].mag = tmp[i][1].mag;
		tmp[i][image_y - 1].gx = tmp[i][image_y - 2].gx;
		tmp[i][image_y - 1].gy = tmp[i][image_y - 2].gy;
		tmp[i][image_y - 1].mag = tmp[i][image_y - 2].mag;
	}
	
	for (j = 1; j <= image_y - 2; j++) {
		tmp[0][j].gx = tmp[1][j].gx;
		tmp[0][j].gy = tmp[1][j].gy;
		tmp[0][j].mag = tmp[1][j].mag;
		tmp[image_x - 1][j].gx = tmp[image_x - 2][j].gx;
		tmp[image_x - 1][j].gy = tmp[image_x - 2][j].gy;
		tmp[image_x - 1][j].mag = tmp[image_x - 2][j].mag;
	}
	
	tmp[0][0].gx = ( tmp[0][1].gx + tmp[1][0].gx ) / 2;
	tmp[0][0].gy = ( tmp[0][1].gy + tmp[1][0].gy ) / 2;
	tmp[0][0].mag = ( tmp[0][1].mag + tmp[1][0].mag ) / 2;
	tmp[0][image_y-1].gx = ( tmp[0][image_y-2].gx + tmp[1][image_y-1].gx ) / 2;
	tmp[0][image_y-1].gy = ( tmp[0][image_y-2].gy + tmp[1][image_y-1].gy ) / 2;
	tmp[0][image_y-1].mag = ( tmp[0][image_y-2].mag + tmp[1][image_y-1].mag ) / 2;
	tmp[image_x-1][0].gx = ( tmp[image_x-1][1].gx + tmp[image_x-2][0].gx ) / 2;
	tmp[image_x-1][0].gy = ( tmp[image_x-1][1].gy + tmp[image_x-2][0].gy ) / 2;
	tmp[image_x-1][0].mag = ( tmp[image_x-1][1].mag + tmp[image_x-2][0].mag ) / 2;
	tmp[image_x - 1][image_y - 1].gx = ( tmp[image_x - 1][image_y - 2].gx + tmp[image_x - 2][image_y - 1].gx ) / 2;
	tmp[image_x - 1][image_y - 1].gy = ( tmp[image_x - 1][image_y - 2].gy + tmp[image_x - 2][image_y - 1].gy ) / 2;
	tmp[image_x - 1][image_y - 1].mag = ( tmp[image_x - 1][image_y - 2].mag + tmp[image_x - 2][image_y - 1].mag ) / 2;

	TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);

	max_grad = round(MAX_GRADIENT);

	// Normalize gradient magnitude
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			tmp[i][j].mag = tmp[i][j].mag / (double)MAX_GRADIENT; // place it between [0, 1]
			G_mag[i][j] = tmp[i][j].mag;
			if (tmp[i][j].mag < 0)
				TRACE("tmp[%d][%d].mag = %0.2f\n", i, j, tmp[i][j].mag);

			image2[i][j] = 255 - round(tmp[i][j].mag * 255.);
		}
	}

	//////////////////////////////////////////////////////////////
	/// Amplify the gradients (strong grad -> stronger, weak grad -> weaker)
	/*
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			p[i][j].mag = pow(2.0, p[i][j].mag);
		}
	}
	*/
}

void GetSobelGradient2(int image_x, int image_y, imatrix& image, imatrix& image2, matrix& G_mag) 
// Sobel gradient version
// makes edge lines white (not black)
{
	int i, j;
	double MAX_GRADIENT = -1.;
	//double MAX_VAL = -1.;
	//double MAX_VAL = 255.;
	double MAX_VAL = 1020.; // 255 + 2 * 255 + 255 (Maximum possible Sobel value)

	Field tmp;

	tmp.init(image_x, image_y);

	/*
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			if (image[i][j] > MAX_VAL) MAX_VAL = (double)image[i][j];
		}
	}
	TRACE("MAX_VAL = %f\n", MAX_VAL);
	*/
	
	for (i = 1; i < image_x - 1; i++) { 
		for (j = 1; j < image_y - 1; j++) {
			////////////////////////////////////////////////////////////////
			// Important!: the value of image intensity should be normalized to [0,1]
			tmp[i][j].gx = (image[i+1][j-1] + 2*(double)image[i+1][j] + image[i+1][j+1] 
				- image[i-1][j-1] - 2*(double)image[i-1][j] - image[i-1][j+1]) / MAX_VAL;
			tmp[i][j].gy = (image[i-1][j+1] + 2*(double)image[i][j+1] + image[i+1][j+1]
				- image[i-1][j-1] - 2*(double)image[i][j-1] - image[i+1][j-1]) / MAX_VAL;
			/////////////////////////////////////////////
			// Amplify!!!
			//p[i][j].gx = pow(p[i][j].gx, 2); 
			//p[i][j].gy = pow(p[i][j].gy, 2); 
			//p[i][j].gx = pow(p[i][j].gx, 2); 
			//TRACE("p[i][j].gx = %.1f\n", p[i][j].gx);
			//TRACE("p[i][j].gy = %.1f\n", p[i][j].gy);
			tmp[i][j].mag = sqrt(tmp[i][j].gx * tmp[i][j].gx + tmp[i][j].gy * tmp[i][j].gy);

			if (tmp[i][j].mag > MAX_GRADIENT) {
				MAX_GRADIENT = tmp[i][j].mag;
				//TRACE("p[i][j].gx = %01.f\n", p[i][j].gx);
				//TRACE("p[i][j].gy = %01.f\n", p[i][j].gy);
				//TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);
			}
		}
	}

	for (i = 1; i <= image_x - 2; i++) {
		tmp[i][0].gx = tmp[i][1].gx;
		tmp[i][0].gy = tmp[i][1].gy;
		tmp[i][0].mag = tmp[i][1].mag;
		tmp[i][image_y - 1].gx = tmp[i][image_y - 2].gx;
		tmp[i][image_y - 1].gy = tmp[i][image_y - 2].gy;
		tmp[i][image_y - 1].mag = tmp[i][image_y - 2].mag;
	}
	
	for (j = 1; j <= image_y - 2; j++) {
		tmp[0][j].gx = tmp[1][j].gx;
		tmp[0][j].gy = tmp[1][j].gy;
		tmp[0][j].mag = tmp[1][j].mag;
		tmp[image_x - 1][j].gx = tmp[image_x - 2][j].gx;
		tmp[image_x - 1][j].gy = tmp[image_x - 2][j].gy;
		tmp[image_x - 1][j].mag = tmp[image_x - 2][j].mag;
	}
	
	tmp[0][0].gx = ( tmp[0][1].gx + tmp[1][0].gx ) / 2;
	tmp[0][0].gy = ( tmp[0][1].gy + tmp[1][0].gy ) / 2;
	tmp[0][0].mag = ( tmp[0][1].mag + tmp[1][0].mag ) / 2;
	tmp[0][image_y-1].gx = ( tmp[0][image_y-2].gx + tmp[1][image_y-1].gx ) / 2;
	tmp[0][image_y-1].gy = ( tmp[0][image_y-2].gy + tmp[1][image_y-1].gy ) / 2;
	tmp[0][image_y-1].mag = ( tmp[0][image_y-2].mag + tmp[1][image_y-1].mag ) / 2;
	tmp[image_x-1][0].gx = ( tmp[image_x-1][1].gx + tmp[image_x-2][0].gx ) / 2;
	tmp[image_x-1][0].gy = ( tmp[image_x-1][1].gy + tmp[image_x-2][0].gy ) / 2;
	tmp[image_x-1][0].mag = ( tmp[image_x-1][1].mag + tmp[image_x-2][0].mag ) / 2;
	tmp[image_x - 1][image_y - 1].gx = ( tmp[image_x - 1][image_y - 2].gx + tmp[image_x - 2][image_y - 1].gx ) / 2;
	tmp[image_x - 1][image_y - 1].gy = ( tmp[image_x - 1][image_y - 2].gy + tmp[image_x - 2][image_y - 1].gy ) / 2;
	tmp[image_x - 1][image_y - 1].mag = ( tmp[image_x - 1][image_y - 2].mag + tmp[image_x - 2][image_y - 1].mag ) / 2;

	TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);

	max_grad = round(MAX_GRADIENT);

	// Normalize gradient magnitude
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			tmp[i][j].mag = tmp[i][j].mag / (double)MAX_GRADIENT; // place it between [0, 1]
			G_mag[i][j] = tmp[i][j].mag;
			if (tmp[i][j].mag < 0)
				TRACE("tmp[%d][%d].mag = %0.2f\n", i, j, tmp[i][j].mag);

			//image2[i][j] = 255 - round(tmp[i][j].mag * 255.);
			image2[i][j] = round(tmp[i][j].mag * 255.);
		}
	}

	//////////////////////////////////////////////////////////////
	/// Amplify the gradients (strong grad -> stronger, weak grad -> weaker)
	/*
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			p[i][j].mag = pow(2.0, p[i][j].mag);
		}
	}
	*/
}


void GetDOG(int image_x, int image_y, imatrix& image, matrix& G_mag, imatrix& image2, int index1, int index2,
			double tau) 
// Difference of Gaussians (just awesome!)
// G_mag gets the actual real values in [0, 1]
// image2 gets the integer values in [0, 255]
{
	int	i, j, k;
	double MAX_DOG = -1.0;
	int s, t;
	int x, y;

	int half1, half2; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	double sum1, sum2, w_sum1, w_sum2;

	//////////////////////////////////
	//Make2DGaussMask(sigma); // mask for unilateral filtering

	//N = gau_mask.getRow();
	//half = N / 2;

	//ClearMemDC(&dc); // clear the canvas white
	//////////////////////////////////////////////
	double factor = 2.0;
	MakeGaussMatrix(factor, GAU); // factor controls the Gaussian window size for each index

	//index1 = 30;
	TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
	//index2 = 100;
	TRACE("GAU_W_2[%d] = %d\n", index2, GAU_W[index2]);

	////////////////////////////////////////////////////////
	// sensitivity of edge detector. between [0, 1] 
	// the smaller, less noise
	//double tau = 0.99; 

	matrix tmp(image_x, image_y);
	//tmp.zero();

	/*
	for (i = 0; i < 256; i++) {
		for (j = 0; j < GAU.getCol(); j++) {
			TRACE("GAU[%d][%d] = %.3f ", i, j, GAU[i][j]);
		}
		TRACE("\n");
	}
	*/
	//gau_w = MakeGaussMask(4.0, gau); // window 20x20
		
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			half1 = GAU_W[index1]-1;
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					sum1 += image[x][y] * GAU[index1][k];
					w_sum1 += GAU[index1][k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			sum1 /= w_sum1; 
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			sum2 = 0;
			w_sum2 = 0.0;
			half2 = GAU_W[index2]-1;
			//TRACE("GAU_W_2[%d] = %d\n", index2, GAU_W[index2]);
			for (s = -half2; s <= half2; s++) {
				for (t = -half2; t <= half2; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					sum2 += image[x][y] * GAU[index2][k];
					w_sum2 += GAU[index2][k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			//TRACE("w_sum = %.6f\n", w_sum);
			sum2 /= w_sum2; 
			//////////////////////////////
			if (sum1 - tau * sum2 > 0)
				tmp[i][j] = 1.0;
			else
				tmp[i][j] = 1.0 + tanh(sum1 - tau * sum2);
			//if (tmp[i][j] > MAX_DOG) MAX_DOG = tmp[i][j];
		}
	}
	////////////////////////////////////////////

	// Normalize DOG value (not necessary! - MAX_DOG is always 1)
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			//tmp[i][j] = tmp[i][j] / MAX_DOG; // place it between [0, 1]
			G_mag[i][j] = 1-tmp[i][j]; // used for nonmaxima suppression
			//TRACE("G_mag[%d][%d] = %.1f\n", i, j, G_mag[i][j]);
			//TRACE("tmp[%d][%d] = %0.3f\n", i, j, tmp[i][j]);
			//TRACE("G_mag[%d][%d] = %0.3f\n", i, j, G_mag[i][j]);
			if (tmp[i][j] < 0)
				TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);
			image2[i][j] = round(tmp[i][j] * 255.);
			//image2[i][j] = 255 - round(tmp[i][j] * 255.);
		}
	}
}

void GetDOG2(imatrix& image, matrix& G_mag, imatrix& image2, double sigma, double tau) 
// Difference of Gaussians (just awesome!)
// Find the two indexes from sigma!
// G_mag gets the actual real values in [0, 1]
// image2 gets the integer values in [0, 255]
{
	int	i, j, k;
	double MAX_DOG = -1.0;
	int s, t;
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half1, half2; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	double sum1, sum2, w_sum1, w_sum2;

	//////////////////////////////////
	//Make2DGaussMask(sigma); // mask for unilateral filtering

	//N = gau_mask.getRow();
	//half = N / 2;

	//ClearMemDC(&dc); // clear the canvas white
	//////////////////////////////////////////////
	//double factor = 2.0;
	//MakeGaussMatrix(factor, GAU); // factor controls the Gaussian window size for each index
	vector GAU1, GAU2;
	MakeGaussVector(sigma, GAU1);
	MakeGaussVector(sigma*1.6, GAU2);
	half1 = GAU1.getMax()-1;
	half2 = GAU2.getMax()-1;

	////////////////////////////////////////////////////////
	matrix tmp(image_x, image_y);
	//tmp.zero();
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					sum1 += image[x][y] * GAU1[k];
					w_sum1 += GAU1[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			sum1 /= w_sum1; 
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			sum2 = 0;
			w_sum2 = 0.0;
			//half2 = GAU_W[index2]-1;
			//half2 = GAU1.getMax()-1; // Now this was a big mistake!!!
			//TRACE("GAU_W_2[%d] = %d\n", index2, GAU_W[index2]);
			for (s = -half2; s <= half2; s++) {
				for (t = -half2; t <= half2; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					//if ( k > half1 ) continue; 
					if ( k > half2 ) continue; 
					/////////////////////////////////////////////////////
					sum2 += image[x][y] * GAU2[k];
					w_sum2 += GAU2[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			//TRACE("w_sum = %.6f\n", w_sum);
			sum2 /= w_sum2; 
			//////////////////////////////
			/*
			// to incorporate both sides of the edge!
			if (sum1 - sum2 > 0) // non edge
				tmp[i][j] = ( 1.0 - tanh(sum1 - sum2) ); // [0, 1]
			else
				tmp[i][j] = 1.0 + tanh(sum1 - sum2); // [ 0, 1 ]
			*/
			///* handle only one side of the edge!
			//TRACE("sum1 = %f\n", sum1);
			//TRACE("sum2 = %f\n", sum2);
			//TRACE("sum1 - tau * sum2 = %f\n", sum1 - tau * sum2);
			/////////////////////////////////////////
			// POSITIVE INNTER CIRCLE
			if (sum1 - tau * sum2 > 0) // non edge
				tmp[i][j] = 1.0;
			else
				tmp[i][j] = 1.0 + tanh(sum1 - tau * sum2); // [ 0, 1 ]
				//tmp[i][j] = 1.0 + tanh( (sum1 - tau * sum2)/10 ); // make it darker!
			///////////////////////////////////////////
			// POSITIVE INNTER CIRCLE (opposite)
			/*
			if (sum1 - tau * sum2 < 0) // non edge
				tmp[i][j] = 1.0;
			else
				tmp[i][j] = -1.0 + tanh(sum1 - tau * sum2); // [ 0, 1 ]
				//tmp[i][j] = 1.0 + tanh( (sum1 - tau * sum2)/10 ); // make it darker!
			*/	
			////////////////////////////////////
			// NEGATIVE INNTER CIRCLE
			/*
			if (-sum1 + tau * sum2 > 0) // non edge
				tmp[i][j] = 1.0;
			else
				tmp[i][j] = 1.0 + tanh(-sum1 + tau * sum2); // [ 0, 1 ]
				//tmp[i][j] = 1.0 + tanh( (sum1 - tau * sum2)/10 ); // make it darker!
			*/
			//if (tmp[i][j] > MAX_DOG) MAX_DOG = tmp[i][j];
		}
	}
	////////////////////////////////////////////

	// Normalize DOG value (not necessary! - MAX_DOG is always 1)
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			//tmp[i][j] = tmp[i][j] / MAX_DOG; // place it between [0, 1]
			//////////////////////////////////////////////////
			G_mag[i][j] = 1-tmp[i][j]; // [0, 1] used for nonmaxima suppression
			// the higher, the stronger the edgeness!
			///////////////////////////////////////////////////
			//TRACE("G_mag[%d][%d] = %.1f\n", i, j, G_mag[i][j]);
			//TRACE("tmp[%d][%d] = %0.3f\n", i, j, tmp[i][j]);
			//TRACE("G_mag[%d][%d] = %0.3f\n", i, j, G_mag[i][j]);
			if (tmp[i][j] < 0)
				TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);
			image2[i][j] = round(tmp[i][j] * 255.);
			//image2[i][j] = 255 - round(tmp[i][j] * 255.);
		}
	}
}

void GetShock(imatrix& image, matrix& G_mag, imatrix& image2, double sigma, double tau) 
// Use DOG to approximate Image Laplacian value!
// G_mag gets the actual real values in [0, 1]
// image2 gets the integer values in [0, 255]
{
	int	i, j, k;
	double MAX_DOG = -1.0;
	int s, t;
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half1, half2; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	double sum1, sum2, w_sum1, w_sum2;

	//////////////////////////////////
	//Make2DGaussMask(sigma); // mask for unilateral filtering

	//N = gau_mask.getRow();
	//half = N / 2;

	//ClearMemDC(&dc); // clear the canvas white
	//////////////////////////////////////////////
	//double factor = 2.0;
	//MakeGaussMatrix(factor, GAU); // factor controls the Gaussian window size for each index
	vector GAU1, GAU2;
	MakeGaussVector(sigma, GAU1);
	MakeGaussVector(sigma*1.6, GAU2);
	half1 = GAU1.getMax()-1;
	half2 = GAU2.getMax()-1;

	////////////////////////////////////////////////////////
	matrix tmp(image_x, image_y);
	//tmp.zero();
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					sum1 += image[x][y] * GAU1[k];
					w_sum1 += GAU1[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			sum1 /= w_sum1; 
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			sum2 = 0;
			w_sum2 = 0.0;
			//half2 = GAU_W[index2]-1;
			//half2 = GAU1.getMax()-1; // Now this was a big mistake!!!
			//TRACE("GAU_W_2[%d] = %d\n", index2, GAU_W[index2]);
			for (s = -half2; s <= half2; s++) {
				for (t = -half2; t <= half2; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					//if ( k > half1 ) continue; 
					if ( k > half2 ) continue; 
					/////////////////////////////////////////////////////
					sum2 += image[x][y] * GAU2[k];
					w_sum2 += GAU2[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			//TRACE("w_sum = %.6f\n", w_sum);
			sum2 /= w_sum2; 
			//////////////////////////////
			/*
			// to incorporate both sides of the edge!
			if (sum1 - sum2 > 0) // non edge
				tmp[i][j] = ( 1.0 - tanh(sum1 - sum2) ); // [0, 1]
			else
				tmp[i][j] = 1.0 + tanh(sum1 - sum2); // [ 0, 1 ]
			*/
			///* handle only one side of the edge!
			//TRACE("sum1 = %f\n", sum1);
			//TRACE("sum2 = %f\n", sum2);
			//TRACE("sum1 - tau * sum2 = %f\n", sum1 - tau * sum2);
			/////////////////////////////////////////
			// POSITIVE INNTER CIRCLE
			if (sum1 - tau * sum2 >= 0) // non edge
				tmp[i][j] = 0.0;
			else 
				tmp[i][j] = 1.0;
			/*
			if (sum1 - tau * sum2 > 0) // non edge
				tmp[i][j] = 1.0;
			else
				tmp[i][j] = 1.0 + tanh(sum1 - tau * sum2); // [ 0, 1 ]
				//tmp[i][j] = 1.0 + tanh( (sum1 - tau * sum2)/10 ); // make it darker!
			*/
		}
	}
	////////////////////////////////////////////
	double weight;
	int sign;
	/////////////////////////////////////////
	// Smoothing step!
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			if (tmp[i][j] == 1.0) sign = 1;
			else sign = 0;
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					if (sign) 
						if (tmp[x][y] == 1.0) weight = 1.0;
						else weight = 0.0;
					else 
						if (tmp[x][y] == 1.0) weight = 0.0;
						else weight = 1.0;
					sum1 += image[x][y] * GAU1[k] * weight;
					w_sum1 += GAU1[k] * weight;
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			sum1 /= w_sum1; 
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			///* handle only one side of the edge!
			//TRACE("sum1 = %f\n", sum1);
			//TRACE("sum2 = %f\n", sum2);
			//TRACE("sum1 - tau * sum2 = %f\n", sum1 - tau * sum2);
			/////////////////////////////////////////
			image2[i][j] = round(sum1);
		}
	}
	////////////////////////////////////////////

	// Normalize DOG value (not necessary! because MAX_DOG is always 1)
	/*
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			//tmp[i][j] = tmp[i][j] / MAX_DOG; // place it between [0, 1]
			//////////////////////////////////////////////////
			G_mag[i][j] = 1-tmp[i][j]; // [0, 1] used for nonmaxima suppression
			// the higher, the stronger the edgeness!
			///////////////////////////////////////////////////
			//TRACE("G_mag[%d][%d] = %.1f\n", i, j, G_mag[i][j]);
			//TRACE("tmp[%d][%d] = %0.3f\n", i, j, tmp[i][j]);
			//TRACE("G_mag[%d][%d] = %0.3f\n", i, j, G_mag[i][j]);
			if (tmp[i][j] < 0)
				TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);
			image2[i][j] = round(tmp[i][j] * 255.);
			//image2[i][j] = 255 - round(tmp[i][j] * 255.);
		}
	}
	*/
}

void GetShock2(imatrix& image, matrix& G_mag, imatrix& image2, double sigma, double tau, double z) 
{
	int	i, j, k;
	double MAX_DOG = -1.0;
	int s, t;
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half1, half2; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	double sum1, sum2, w_sum1, w_sum2;

	//////////////////////////////////
	//Make2DGaussMask(sigma); // mask for unilateral filtering

	//N = gau_mask.getRow();
	//half = N / 2;

	//ClearMemDC(&dc); // clear the canvas white
	//////////////////////////////////////////////
	//double factor = 2.0;
	//MakeGaussMatrix(factor, GAU); // factor controls the Gaussian window size for each index
	vector GAU1, GAU2;
	MakeGaussVector(sigma, GAU1);
	MakeGaussVector(sigma*1.6, GAU2);
	half1 = GAU1.getMax()-1;
	half2 = GAU2.getMax()-1;

	deque<int> vec;
	////////////////////////////////////////////////////////
	matrix tmp(image_x, image_y);
	imatrix tmp2(image_x, image_y);
	//tmp.zero();
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					sum1 += image[x][y] * GAU1[k];
					w_sum1 += GAU1[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			sum1 /= w_sum1; 
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			sum2 = 0;
			w_sum2 = 0.0;
			//half2 = GAU_W[index2]-1;
			//half2 = GAU1.getMax()-1; // Now this was a big mistake!!!
			//TRACE("GAU_W_2[%d] = %d\n", index2, GAU_W[index2]);
			for (s = -half2; s <= half2; s++) {
				for (t = -half2; t <= half2; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					//if ( k > half1 ) continue; 
					if ( k > half2 ) continue; 
					/////////////////////////////////////////////////////
					sum2 += image[x][y] * GAU2[k];
					w_sum2 += GAU2[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			//TRACE("w_sum = %.6f\n", w_sum);
			sum2 /= w_sum2; 
			//////////////////////////////
			/////////////////////////////////////////
			// POSITIVE INNTER CIRCLE
			if (sum1 - tau * sum2 >= 0) // non edge
				tmp[i][j] = 1.0; // non edge
			else 
				tmp[i][j] = 0.0; // edge
		}
	}
	////////////////////////////////////////////
	//double weight;
	int sign;
	/////////////////////////////////////////
	// Smoothing step!
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			if (tmp[i][j] == 1.0) sign = 1; // Non edge region
			else sign = 0; // edge region
			vec.clear();
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					if (sign) { // non edge region
                        if (tmp[x][y] == 1.0) vec.push_back(image[x][y]);
					}
					else { // edge region
						if (tmp[x][y] == 0.0) vec.push_back(image[x][y]);
					}
					//sum1 += image[x][y] * GAU1[k] * weight;
					//w_sum1 += GAU1[k] * weight;
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			//sum1 /= w_sum1; 
			sort(vec.begin(), vec.end()); 
			//flow_median = vec[vec.size()/2]; 
			//double z = 0.3;
			if (sign) // non edge region
				tmp2[i][j] = vec[vec.size()-1]; 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
			else // edge region
				//tmp2[i][j] = vec[0]; 
				tmp2[i][j] = round( (1-z) * vec[0] + z * 0 ); // darkest value
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			///* handle only one side of the edge!
			//TRACE("sum1 = %f\n", sum1);
			//TRACE("sum2 = %f\n", sum2);
			//TRACE("sum1 - tau * sum2 = %f\n", sum1 - tau * sum2);
			/////////////////////////////////////////
			//image2[i][j] = round(sum1);
		}
	}
	////////////////////////////////////////////

	// Normalize DOG value (not necessary! because MAX_DOG is always 1)
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			//tmp[i][j] = tmp[i][j] / MAX_DOG; // place it between [0, 1]
			//////////////////////////////////////////////////
			//G_mag[i][j] = 1-tmp[i][j]; // [0, 1] used for nonmaxima suppression
			// the higher, the stronger the edgeness!
			///////////////////////////////////////////////////
			//TRACE("G_mag[%d][%d] = %.1f\n", i, j, G_mag[i][j]);
			//TRACE("tmp[%d][%d] = %0.3f\n", i, j, tmp[i][j]);
			//TRACE("G_mag[%d][%d] = %0.3f\n", i, j, G_mag[i][j]);
			//if (tmp[i][j] < 0)
			//	TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);
			image2[i][j] = tmp2[i][j];
			//image2[i][j] = 255 - round(tmp[i][j] * 255.);
		}
	}
}

void GetColShock(cimatrix& cmap, imatrix& image, double sigma, double tau, double z) 
{
	int	i, j, k;
	double MAX_DOG = -1.0;
	int s, t;
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half1, half2; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	double sum1, sum2, w_sum1, w_sum2;

	//////////////////////////////////
	//Make2DGaussMask(sigma); // mask for unilateral filtering

	//N = gau_mask.getRow();
	//half = N / 2;

	//ClearMemDC(&dc); // clear the canvas white
	//////////////////////////////////////////////
	//double factor = 2.0;
	//MakeGaussMatrix(factor, GAU); // factor controls the Gaussian window size for each index
	vector GAU1, GAU2;
	MakeGaussVector(sigma, GAU1);
	MakeGaussVector(sigma*1.6, GAU2);
	half1 = GAU1.getMax()-1;
	half2 = GAU2.getMax()-1;

	deque<int> vec[3];
	////////////////////////////////////////////////////////
	matrix tmp(image_x, image_y);
	cimatrix tmp2(image_x, image_y);
	//tmp.zero();
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					sum1 += image[x][y] * GAU1[k];
					w_sum1 += GAU1[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			sum1 /= w_sum1; 
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			sum2 = 0;
			w_sum2 = 0.0;
			//half2 = GAU_W[index2]-1;
			//half2 = GAU1.getMax()-1; // Now this was a big mistake!!!
			//TRACE("GAU_W_2[%d] = %d\n", index2, GAU_W[index2]);
			for (s = -half2; s <= half2; s++) {
				for (t = -half2; t <= half2; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					//if ( k > half1 ) continue; 
					if ( k > half2 ) continue; 
					/////////////////////////////////////////////////////
					sum2 += image[x][y] * GAU2[k];
					w_sum2 += GAU2[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			//TRACE("w_sum = %.6f\n", w_sum);
			sum2 /= w_sum2; 
			//////////////////////////////
			/////////////////////////////////////////
			// POSITIVE INNTER CIRCLE
			if (sum1 - tau * sum2 >= 0) // non edge
				tmp[i][j] = 1.0; // non edge
			else 
				tmp[i][j] = 0.0; // edge
		}
	}

	////////////////////////////////////////////
	//double weight;
	int sign;
	/////////////////////////////////////////
	// Smoothing step!
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			if (tmp[i][j] == 1.0) sign = 1; // Non edge region
			else sign = 0; // edge region
			vec[0].clear();
			vec[1].clear();
			vec[2].clear();
			vec[0].push_back(cmap[i][j].r);
			vec[1].push_back(cmap[i][j].g);
			vec[2].push_back(cmap[i][j].b);
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					if (sign) { // non edge region
						if (tmp[x][y] == 1.0) {
							vec[0].push_back(cmap[x][y].r);
							vec[1].push_back(cmap[x][y].g);
							vec[2].push_back(cmap[x][y].b);
						}
					}
					else { // edge region
						if (tmp[x][y] == 0.0) {
							vec[0].push_back(cmap[x][y].r);
							vec[1].push_back(cmap[x][y].g);
							vec[2].push_back(cmap[x][y].b);
						}
					}
					//sum1 += image[x][y] * GAU1[k] * weight;
					//w_sum1 += GAU1[k] * weight;
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			//sum1 /= w_sum1; 
			sort(vec[0].begin(), vec[0].end()); 
			sort(vec[1].begin(), vec[1].end()); 
			sort(vec[2].begin(), vec[2].end()); 
			//flow_median = vec[vec.size()/2]; 
			//double z = 0.3;
			if (sign) { // non edge region
				//tmp2[i][j].r = vec[0][vec[0].size()-1]; 
				//tmp2[i][j].g = vec[1][vec[1].size()-1]; 
				//tmp2[i][j].b = vec[2][vec[2].size()-1]; 
				tmp2[i][j].r = vec[0][vec[0].size()/2]; 
				tmp2[i][j].g = vec[1][vec[1].size()/2]; 
				tmp2[i][j].b = vec[2][vec[2].size()/2]; 
				//tmp2[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
				//tmp2[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
				//tmp2[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
			}
			else { // edge region
				//tmp2[i][j] = vec[0]; 
				tmp2[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
				tmp2[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
				tmp2[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
			}
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			///* handle only one side of the edge!
			//TRACE("sum1 = %f\n", sum1);
			//TRACE("sum2 = %f\n", sum2);
			//TRACE("sum1 - tau * sum2 = %f\n", sum1 - tau * sum2);
			/////////////////////////////////////////
			//image2[i][j] = round(sum1);
		}
	}
	////////////////////////////////////////////

	// Normalize DOG value (not necessary! because MAX_DOG is always 1)
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			//tmp[i][j] = tmp[i][j] / MAX_DOG; // place it between [0, 1]
			//////////////////////////////////////////////////
			//G_mag[i][j] = 1-tmp[i][j]; // [0, 1] used for nonmaxima suppression
			// the higher, the stronger the edgeness!
			///////////////////////////////////////////////////
			//TRACE("G_mag[%d][%d] = %.1f\n", i, j, G_mag[i][j]);
			//TRACE("tmp[%d][%d] = %0.3f\n", i, j, tmp[i][j]);
			//TRACE("G_mag[%d][%d] = %0.3f\n", i, j, G_mag[i][j]);
			//if (tmp[i][j] < 0)
			//	TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);
			cmap[i][j].r = tmp2[i][j].r;
			cmap[i][j].g = tmp2[i][j].g;
			cmap[i][j].b = tmp2[i][j].b;
			//image2[i][j] = 255 - round(tmp[i][j] * 255.);
		}
	}
}

void GetColShock2(cimatrix& cmap, imatrix& image, double sigma, double tau, double z) 
// Multiply the Shock value with the DOG value!
{
	int	i, j, k;
	double MAX_DOG = -1.0;
	int s, t;
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half1, half2; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	double sum1, sum2, w_sum1, w_sum2;
	double z2;

	//////////////////////////////////
	//Make2DGaussMask(sigma); // mask for unilateral filtering

	//N = gau_mask.getRow();
	//half = N / 2;

	//ClearMemDC(&dc); // clear the canvas white
	//////////////////////////////////////////////
	//double factor = 2.0;
	//MakeGaussMatrix(factor, GAU); // factor controls the Gaussian window size for each index
	vector GAU1, GAU2;
	MakeGaussVector(sigma, GAU1);
	MakeGaussVector(sigma*1.6, GAU2);
	half1 = GAU1.getMax()-1;
	half2 = GAU2.getMax()-1;

	deque<int> vec[3];
	////////////////////////////////////////////////////////
	matrix tmp(image_x, image_y);
	cimatrix tmp2(image_x, image_y);
	//tmp.zero();
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					sum1 += image[x][y] * GAU1[k];
					w_sum1 += GAU1[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			sum1 /= w_sum1; 
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			sum2 = 0;
			w_sum2 = 0.0;
			//half2 = GAU_W[index2]-1;
			//half2 = GAU1.getMax()-1; // Now this was a big mistake!!!
			//TRACE("GAU_W_2[%d] = %d\n", index2, GAU_W[index2]);
			for (s = -half2; s <= half2; s++) {
				for (t = -half2; t <= half2; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					//if ( k > half1 ) continue; 
					if ( k > half2 ) continue; 
					/////////////////////////////////////////////////////
					sum2 += image[x][y] * GAU2[k];
					w_sum2 += GAU2[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			//TRACE("w_sum = %.6f\n", w_sum);
			sum2 /= w_sum2; 
			//////////////////////////////
			/////////////////////////////////////////
			tmp[i][j] = sum1 - tau * sum2;
			//TRACE("tmp[%d][%d] = %f\n", i, j, tmp[i][j]);
			//if (sum1 - tau * sum2 >= 0) // non edge
			//	tmp[i][j] = 1.0; // non edge
			//else 
			//	tmp[i][j] = 0.0; // edge
		}
	}
	////////////////////////////////////////////
	//double weight;
	int sign;
	/////////////////////////////////////////
	// Smoothing step!
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			///////////////////////////////////
			if (tmp[i][j] >= 0.0) sign = 1; // Non edge region
			else sign = 0; // edge region
			/////////////////////////////////////
			vec[0].clear();
			vec[1].clear();
			vec[2].clear();
			vec[0].push_back(cmap[i][j].r);
			vec[1].push_back(cmap[i][j].g);
			vec[2].push_back(cmap[i][j].b);
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					if (sign) { // non edge region
						if (tmp[x][y] >= 0.0) {
							vec[0].push_back(cmap[x][y].r);
							vec[1].push_back(cmap[x][y].g);
							vec[2].push_back(cmap[x][y].b);
						}
					}
					else { // edge region
						if (tmp[x][y] < 0.0) {
							vec[0].push_back(cmap[x][y].r);
							vec[1].push_back(cmap[x][y].g);
							vec[2].push_back(cmap[x][y].b);
						}
					}
					//sum1 += image[x][y] * GAU1[k] * weight;
					//w_sum1 += GAU1[k] * weight;
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			//sum1 /= w_sum1; 
			sort(vec[0].begin(), vec[0].end()); 
			sort(vec[1].begin(), vec[1].end()); 
			sort(vec[2].begin(), vec[2].end()); 
			//flow_median = vec[vec.size()/2]; 
			//double z = 0.3;
			if (sign) { // non edge region
				//tmp2[i][j].r = vec[0][vec[0].size()-1]; 
				//tmp2[i][j].g = vec[1][vec[1].size()-1]; 
				//tmp2[i][j].b = vec[2][vec[2].size()-1]; 
				tmp2[i][j].r = vec[0][vec[0].size()/2]; 
				tmp2[i][j].g = vec[1][vec[1].size()/2]; 
				tmp2[i][j].b = vec[2][vec[2].size()/2]; 
				//tmp2[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
				//tmp2[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
				//tmp2[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
				z2 = tmp[i][j]/10.0; // positive
				//TRACE("tmp[%d][%d] = %f\n", i, j, tmp[i][j]);
				tmp2[i][j].r = round( (1-z2) * tmp2[i][j].r + z2 * 255 ); // darkest value
				tmp2[i][j].g = round( (1-z2) * tmp2[i][j].g + z2 * 255 ); // darkest value
				tmp2[i][j].b = round( (1-z2) * tmp2[i][j].b + z2 * 255 ); // darkest value
			
			}
			else { // edge region
				//tmp2[i][j] = vec[0]; 
				tmp2[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
				tmp2[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
				tmp2[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
				z2 = tmp[i][j]/10; // negative
				//TRACE("tmp[%d][%d] = %f\n", i, j, tmp[i][j]);
				tmp2[i][j].r = round( (1-z2) * tmp2[i][j].r + z2 * 0 ); // darkest value
				tmp2[i][j].g = round( (1-z2) * tmp2[i][j].g + z2 * 0 ); // darkest value
				tmp2[i][j].b = round( (1-z2) * tmp2[i][j].b + z2 * 0 ); // darkest value
			}
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			///* handle only one side of the edge!
			//TRACE("sum1 = %f\n", sum1);
			//TRACE("sum2 = %f\n", sum2);
			//TRACE("sum1 - tau * sum2 = %f\n", sum1 - tau * sum2);
			/////////////////////////////////////////
			//image2[i][j] = round(sum1);
		}
	}
	////////////////////////////////////////////

	// Normalize DOG value (not necessary! because MAX_DOG is always 1)
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			//tmp[i][j] = tmp[i][j] / MAX_DOG; // place it between [0, 1]
			//////////////////////////////////////////////////
			//G_mag[i][j] = 1-tmp[i][j]; // [0, 1] used for nonmaxima suppression
			// the higher, the stronger the edgeness!
			///////////////////////////////////////////////////
			//TRACE("G_mag[%d][%d] = %.1f\n", i, j, G_mag[i][j]);
			//TRACE("tmp[%d][%d] = %0.3f\n", i, j, tmp[i][j]);
			//TRACE("G_mag[%d][%d] = %0.3f\n", i, j, G_mag[i][j]);
			//if (tmp[i][j] < 0)
			//	TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);
			cmap[i][j].r = tmp2[i][j].r;
			cmap[i][j].g = tmp2[i][j].g;
			cmap[i][j].b = tmp2[i][j].b;
			//image2[i][j] = 255 - round(tmp[i][j] * 255.);
		}
	}
}

void GetColShock3(cimatrix& cmap, imatrix& image, double sigma, double tau, int half, double z) 
// use Min-Max filtering!
{
	int	i, j, k;
	double MAX_DOG = -1.0;
	int s, t;
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half1, half2; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	double sum1, sum2, w_sum1, w_sum2;

	vector GAU1, GAU2;
	MakeGaussVector(sigma, GAU1);
	MakeGaussVector(sigma*1.6, GAU2);
	half1 = GAU1.getMax()-1;
	half2 = GAU2.getMax()-1;

	deque<int> vec[3];
	////////////////////////////////////////////////////////
	matrix tmp(image_x, image_y);
	cimatrix tmp2(image_x, image_y);
	//tmp.zero();
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					sum1 += image[x][y] * GAU1[k];
					w_sum1 += GAU1[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			sum1 /= w_sum1; 
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			sum2 = 0;
			w_sum2 = 0.0;
			//half2 = GAU_W[index2]-1;
			//half2 = GAU1.getMax()-1; // Now this was a big mistake!!!
			//TRACE("GAU_W_2[%d] = %d\n", index2, GAU_W[index2]);
			for (s = -half2; s <= half2; s++) {
				for (t = -half2; t <= half2; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					//if ( k > half1 ) continue; 
					if ( k > half2 ) continue; 
					/////////////////////////////////////////////////////
					sum2 += image[x][y] * GAU2[k];
					w_sum2 += GAU2[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			//TRACE("w_sum = %.6f\n", w_sum);
			sum2 /= w_sum2; 
			//////////////////////////////
			/////////////////////////////////////////
			// POSITIVE INNTER CIRCLE
			if (sum1 - tau * sum2 >= 0) // non edge
				tmp[i][j] = 1.0; // non edge
			else 
				tmp[i][j] = 0.0; // edge
		}
	}

	////////////////////////////////////////////
	//double weight;
	int sign;
	half1 = half;
	/////////////////////////////////////////
	// Smoothing step!
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			if (tmp[i][j] == 1.0) sign = 1; // Non edge region
			else sign = 0; // edge region
			vec[0].clear();
			vec[1].clear();
			vec[2].clear();
			vec[0].push_back(cmap[i][j].r);
			vec[1].push_back(cmap[i][j].g);
			vec[2].push_back(cmap[i][j].b);
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					if (sign) { // non edge region
						if (tmp[x][y] == 1.0) {
							vec[0].push_back(cmap[x][y].r);
							vec[1].push_back(cmap[x][y].g);
							vec[2].push_back(cmap[x][y].b);
						}
					}
					else { // edge region
						if (tmp[x][y] == 0.0) {
							vec[0].push_back(cmap[x][y].r);
							vec[1].push_back(cmap[x][y].g);
							vec[2].push_back(cmap[x][y].b);
						}
					}
					//sum1 += image[x][y] * GAU1[k] * weight;
					//w_sum1 += GAU1[k] * weight;
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			//sum1 /= w_sum1; 
			sort(vec[0].begin(), vec[0].end()); 
			sort(vec[1].begin(), vec[1].end()); 
			sort(vec[2].begin(), vec[2].end()); 
			//flow_median = vec[vec.size()/2]; 
			//double z = 0.3;
			if (sign) { // non edge region
				//tmp2[i][j].r = vec[0][vec[0].size()-1]; 
				//tmp2[i][j].g = vec[1][vec[1].size()-1]; 
				//tmp2[i][j].b = vec[2][vec[2].size()-1]; 
				tmp2[i][j].r = vec[0][vec[0].size()/2]; 
				tmp2[i][j].g = vec[1][vec[1].size()/2]; 
				tmp2[i][j].b = vec[2][vec[2].size()/2]; 
				//tmp2[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
				//tmp2[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
				//tmp2[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
			}
			else { // edge region
				//tmp2[i][j].r = vec[0][0];
				//tmp2[i][j].g = vec[1][0];
				//tmp2[i][j].b = vec[2][0];
				tmp2[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
				tmp2[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
				tmp2[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
			}
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			///* handle only one side of the edge!
			//TRACE("sum1 = %f\n", sum1);
			//TRACE("sum2 = %f\n", sum2);
			//TRACE("sum1 - tau * sum2 = %f\n", sum1 - tau * sum2);
			/////////////////////////////////////////
			//image2[i][j] = round(sum1);
		}
	}
	////////////////////////////////////////////

	// Normalize DOG value (not necessary! because MAX_DOG is always 1)
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			//tmp[i][j] = tmp[i][j] / MAX_DOG; // place it between [0, 1]
			//////////////////////////////////////////////////
			//G_mag[i][j] = 1-tmp[i][j]; // [0, 1] used for nonmaxima suppression
			// the higher, the stronger the edgeness!
			///////////////////////////////////////////////////
			//TRACE("G_mag[%d][%d] = %.1f\n", i, j, G_mag[i][j]);
			//TRACE("tmp[%d][%d] = %0.3f\n", i, j, tmp[i][j]);
			//TRACE("G_mag[%d][%d] = %0.3f\n", i, j, G_mag[i][j]);
			//if (tmp[i][j] < 0)
			//	TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);
			cmap[i][j].r = tmp2[i][j].r;
			cmap[i][j].g = tmp2[i][j].g;
			cmap[i][j].b = tmp2[i][j].b;
			//image2[i][j] = 255 - round(tmp[i][j] * 255.);
		}
	}
}

void GetColShockFast(cimatrix& cmap, imatrix& image, double sigma, double tau, double z) 
// faster version
{
	int	i, j;
	double MAX_DOG = -1.0;
	int s, t;
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half1, half2; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	double sum1, w_sum1;

	vector GAU1, GAU2;
	MakeGaussVector(sigma, GAU1);
	MakeGaussVector(sigma*1.6, GAU2);
	half1 = GAU1.getMax()-1;
	half2 = GAU2.getMax()-1;

	deque<int> vec[3];
	////////////////////////////////////////////////////////
	//matrix tmp(image_x, image_y);
	cimatrix tmp2(image_x, image_y);
	
	matrix inner(image_x, image_y);
	GaussSmoothSepDouble(image, sigma, inner);

	matrix outer(image_x, image_y);
	GaussSmoothSepDouble(image, sigma*1.6, outer);

	imatrix tmp(image_x, image_y);

	double dog;
	//deque<int> vec[3];
	//////////////////////////////
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			dog = inner[i][j] - tau * outer[i][j];
			if (dog >= 0) tmp[i][j] = 1; 
			//if (dog >= thres) lap[i][j] = 255; 
			//else if (tanh(dog) + 1 > thres) lap[i][j] = 255; 
			//else if (tanh(dog) + 1 > 0.1) lap[i][j] = 255; 
			else tmp[i][j] = 0; 
		}
	}

	////////////////////////////////////////////
	//double weight;
	int sign;
	//half1 = half;
	/////////////////////////////////////////
	// Smoothing step!
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			if (tmp[i][j] == 1) sign = 1; // Non edge region
			else sign = 0; // edge region
			vec[0].clear();
			vec[1].clear();
			vec[2].clear();
			vec[0].push_back(cmap[i][j].r);
			vec[1].push_back(cmap[i][j].g);
			vec[2].push_back(cmap[i][j].b);
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (s = -half1; s <= half1; s++) {
				////////////////////////
				x = i+s; y = j;
				/////////////////////////////////////////////////////
				if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) 
					continue;
				/////////////////////////////////////////////////////////
				// circular kernel
				//k = (int)dist2(x, y, i, j);
				//k = round( dist2(x, y, i, j) );
				//if ( k > half1 ) continue; 
				/////////////////////////////////////////////////////
				if (sign) { // non edge region
					if (tmp[x][y] == 1) {
						vec[0].push_back(cmap[x][y].r);
						vec[1].push_back(cmap[x][y].g);
						vec[2].push_back(cmap[x][y].b);
					}
				}
				else { // edge region
					if (tmp[x][y] == 0) {
						vec[0].push_back(cmap[x][y].r);
						vec[1].push_back(cmap[x][y].g);
						vec[2].push_back(cmap[x][y].b);
					}
				}
				//sum1 += image[x][y] * GAU1[k] * weight;
				//w_sum1 += GAU1[k] * weight;
				//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
				//TRACE("k = %d\n", k);
			}
			for (t = -half1; t <= half1; t++) {
				////////////////////////
				x = i; y = j+t;
				/////////////////////////////////////////////////////
				if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) 
					continue;
				/////////////////////////////////////////////////////////
				// circular kernel
				//k = (int)dist2(x, y, i, j);
				//k = round( dist2(x, y, i, j) );
				//if ( k > half1 ) continue; 
				/////////////////////////////////////////////////////
				if (sign) { // non edge region
					if (tmp[x][y] == 1) {
						vec[0].push_back(cmap[x][y].r);
						vec[1].push_back(cmap[x][y].g);
						vec[2].push_back(cmap[x][y].b);
					}
				}
				else { // edge region
					if (tmp[x][y] == 0) {
						vec[0].push_back(cmap[x][y].r);
						vec[1].push_back(cmap[x][y].g);
						vec[2].push_back(cmap[x][y].b);
					}
				}
				//sum1 += image[x][y] * GAU1[k] * weight;
				//w_sum1 += GAU1[k] * weight;
				//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
				//TRACE("k = %d\n", k);
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			//sum1 /= w_sum1; 
			sort(vec[0].begin(), vec[0].end()); 
			sort(vec[1].begin(), vec[1].end()); 
			sort(vec[2].begin(), vec[2].end()); 
			//flow_median = vec[vec.size()/2]; 
			//double z = 0.3;
			if (sign) { // non edge region
				tmp2[i][j].r = vec[0][vec[0].size()-1]; 
				tmp2[i][j].g = vec[1][vec[1].size()-1]; 
				tmp2[i][j].b = vec[2][vec[2].size()-1]; 
				//tmp2[i][j].r = vec[0][vec[0].size()/2]; 
				//tmp2[i][j].g = vec[1][vec[1].size()/2]; 
				//tmp2[i][j].b = vec[2][vec[2].size()/2]; 
				//tmp2[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
				//tmp2[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
				//tmp2[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
			}
			else { // edge region
				tmp2[i][j].r = vec[0][0];
				tmp2[i][j].g = vec[1][0];
				tmp2[i][j].b = vec[2][0];
				//tmp2[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
				//tmp2[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
				//tmp2[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
			}
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
			///* handle only one side of the edge!
			//TRACE("sum1 = %f\n", sum1);
			//TRACE("sum2 = %f\n", sum2);
			//TRACE("sum1 - tau * sum2 = %f\n", sum1 - tau * sum2);
			/////////////////////////////////////////
			//image2[i][j] = round(sum1);
		}
	}
	////////////////////////////////////////////

	// Normalize DOG value (not necessary! because MAX_DOG is always 1)
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			//tmp[i][j] = tmp[i][j] / MAX_DOG; // place it between [0, 1]
			//////////////////////////////////////////////////
			//G_mag[i][j] = 1-tmp[i][j]; // [0, 1] used for nonmaxima suppression
			// the higher, the stronger the edgeness!
			///////////////////////////////////////////////////
			//TRACE("G_mag[%d][%d] = %.1f\n", i, j, G_mag[i][j]);
			//TRACE("tmp[%d][%d] = %0.3f\n", i, j, tmp[i][j]);
			//TRACE("G_mag[%d][%d] = %0.3f\n", i, j, G_mag[i][j]);
			//if (tmp[i][j] < 0)
			//	TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);
			cmap[i][j].r = tmp2[i][j].r;
			cmap[i][j].g = tmp2[i][j].g;
			cmap[i][j].b = tmp2[i][j].b;
			//image2[i][j] = 255 - round(tmp[i][j] * 255.);
		}
	}
}


void GetColShockFastHSV(cimatrix& cmap, imatrix& image, double sigma, double tau) 
// faster version
{
	int	i, j;
	double MAX_DOG = -1.0;
	int ss, t;
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half1; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	double sum1, w_sum1;

	vector GAU1, GAU2;
	MakeGaussVector(sigma, GAU1);
	//MakeGaussVector(sigma*1.6, GAU2);
	half1 = GAU1.getMax()-1;
	//half2 = GAU2.getMax()-1;

	deque<double> vec;
	////////////////////////////////////////////////////////
	//matrix tmp(image_x, image_y);
	cimatrix tmp2(image_x, image_y);
	
	matrix inner(image_x, image_y);
	GaussSmoothSepDouble(image, sigma, inner);

	matrix outer(image_x, image_y);
	GaussSmoothSepDouble(image, sigma*1.6, outer);

	imatrix tmp(image_x, image_y);

	double dog;
	//deque<int> vec[3];
	//////////////////////////////
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			dog = inner[i][j] - tau * outer[i][j];
			if (dog >= 0) tmp[i][j] = 1; 
			//if (dog >= thres) lap[i][j] = 255; 
			//else if (tanh(dog) + 1 > thres) lap[i][j] = 255; 
			//else if (tanh(dog) + 1 > 0.1) lap[i][j] = 255; 
			else tmp[i][j] = 0; 
		}
	}

	////////////////////////////////////////////
	//double weight;
	int sign;
	double r, g, b, h, s, v, h1, s1, v1;
	//half1 = half;
	/////////////////////////////////////////
	// Smoothing step!
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			if (tmp[i][j] == 1) sign = 1; // Non edge region
			else sign = 0; // edge region
			
			r = (double)cmap[i][j].r/255.;
			g = (double)cmap[i][j].g/255.;
			b = (double)cmap[i][j].b/255.;
			RGB2HSV(r, g, b, h1, s1, v1);

			vec.clear();
			//vec.push_back(v);
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			for (ss = -half1; ss <= half1; ss++) {
				////////////////////////
				x = i+ss; y = j;
				/////////////////////////////////////////////////////
				if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) 
					continue;
				/////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////
				r = (double)cmap[x][y].r/255.;
				g = (double)cmap[x][y].g/255.;
				b = (double)cmap[x][y].b/255.;
				RGB2HSV(r, g, b, h, s, v);
				vec.push_back(v);
				//sum1 += image[x][y] * GAU1[k] * weight;
				//w_sum1 += GAU1[k] * weight;
				//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
				//TRACE("k = %d\n", k);
			}
			for (t = -half1; t <= half1; t++) {
				////////////////////////
				x = i; y = j+t;
				/////////////////////////////////////////////////////
				if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) 
					continue;
				/////////////////////////////////////////////////////
				r = (double)cmap[x][y].r/255.;
				g = (double)cmap[x][y].g/255.;
				b = (double)cmap[x][y].b/255.;
				RGB2HSV(r, g, b, h, s, v);
				vec.push_back(v);
			}
			////////////////////////////////////
			sort(vec.begin(), vec.end()); 
			//double z = 0.3;
			if (sign) { // non edge region
				v1 = vec[vec.size()-1];
			}
			else { // edge region
				v1 = vec[0];
			}
			HSV2RGB(h1, s1, v1, r, g, b);
			tmp2[i][j].r = (GLubyte)(r*255);
			tmp2[i][j].g = (GLubyte)(g*255);
			tmp2[i][j].b = (GLubyte)(b*255);
			//TRACE("w_sum = %.6f\n", w_sum);
			//////////////////////////////////////////////
		}
	}
	////////////////////////////////////////////

	// Normalize DOG value (not necessary! because MAX_DOG is always 1)
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			//tmp[i][j] = tmp[i][j] / MAX_DOG; // place it between [0, 1]
			//////////////////////////////////////////////////
			//G_mag[i][j] = 1-tmp[i][j]; // [0, 1] used for nonmaxima suppression
			// the higher, the stronger the edgeness!
			///////////////////////////////////////////////////
			//TRACE("G_mag[%d][%d] = %.1f\n", i, j, G_mag[i][j]);
			//TRACE("tmp[%d][%d] = %0.3f\n", i, j, tmp[i][j]);
			//TRACE("G_mag[%d][%d] = %0.3f\n", i, j, G_mag[i][j]);
			//if (tmp[i][j] < 0)
			//	TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);
			cmap[i][j].r = tmp2[i][j].r;
			cmap[i][j].g = tmp2[i][j].g;
			cmap[i][j].b = tmp2[i][j].b;
			//image2[i][j] = 255 - round(tmp[i][j] * 255.);
		}
	}
}

void GetColShockFastMinMax(cimatrix& cmap, imatrix& image, double sigma, double tau) 
// we find min and max in terms of R + G + B
{
	int	i, j;
	double MAX_DOG = -1.0;
	int ss, t;
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half1; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	double sum1, w_sum1;

	vector GAU1, GAU2;
	MakeGaussVector(sigma, GAU1);
	//MakeGaussVector(sigma*1.6, GAU2);
	half1 = GAU1.getMax()-1;
	//half2 = GAU2.getMax()-1;

	//deque<double> vec;
	////////////////////////////////////////////////////////
	//matrix tmp(image_x, image_y);
	cimatrix tmp2(image_x, image_y);
	
	matrix inner(image_x, image_y);
	GaussSmoothSepDouble(image, sigma, inner);

	matrix outer(image_x, image_y);
	GaussSmoothSepDouble(image, sigma*1.6, outer);

	imatrix tmp(image_x, image_y);

	double dog;
	//deque<int> vec[3];
	//double thres = 0.7;
	//////////////////////////////
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			dog = inner[i][j] - tau * outer[i][j];
			if (dog >= 0.0) tmp[i][j] = 1; 
			//if (dog >= thres) tmp[i][j] = 1; 
			//else if (tanh(dog) + 1 > thres) tmp[i][j] = 1; 
			//else if (tanh(dog) + 1 > 0.1) lap[i][j] = 255; 
			else tmp[i][j] = 0; 
		}
	}

	////////////////////////////////////////////
	//double weight;
	int sign;
	int r, g, b;
	int min, max, sum;
	GLubyte min_r, min_g, min_b, max_r, max_g, max_b;
	//half1 = half;
	half1 = 1;
	/////////////////////////////////////////
	// Smoothing step!
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			if (tmp[i][j] == 1) sign = 1; // Non edge region
			else sign = 0; // edge region

			max = -1;
			min = 1000000;
			
			for (ss = -half1; ss <= half1; ss++) {
				////////////////////////
				x = i+ss; y = j;
				/////////////////////////////////////////////////////
				if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) 
					continue;
				/////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////
				if (sign == tmp[x][y]) {
					r = cmap[x][y].r;
					g = cmap[x][y].g;
					b = cmap[x][y].b;
					sum = r + g + b;
					if (sum > max) { 
						max = sum; 
						max_r = r; max_g = g; max_b = b;
					}
					if (r + g + b < min) { 
						min = sum; 
						min_r = r; min_g = g; min_b = b;
					}
					//TRACE("k = %d\n", k)
				}
			}
			////////////////////////////////////
			if (sign) { // non edge region
				tmp2[i][j].r = max_r;
				tmp2[i][j].g = max_g;
				tmp2[i][j].b = max_b;
			}
			else { // edge region
				tmp2[i][j].r = min_r;
				tmp2[i][j].g = min_g;
				tmp2[i][j].b = min_b;
			}

		}
	}
	////////////////////////////////////////////
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			if (tmp[i][j] == 1) sign = 1; // Non edge region
			else sign = 0; // edge region

			max = -1;
			min = 1000000;
			
			for (t = -half1; t <= half1; t++) {
				////////////////////////
				x = i; y = j+t;
				/////////////////////////////////////////////////////
				if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) 
					continue;
				/////////////////////////////////////////////////////
				if (sign == tmp[x][y]) {
					r = tmp2[x][y].r;
					g = tmp2[x][y].g;
					b = tmp2[x][y].b;
					sum = r + g + b;
					if (sum > max) { 
						max = sum; 
						max_r = r; max_g = g; max_b = b;
					}
					if (r + g + b < min) { 
						min = sum; 
						min_r = r; min_g = g; min_b = b;
					}
				}
			}
			////////////////////////////////////
			if (sign) { // non edge region
				cmap[i][j].r = max_r;
				cmap[i][j].g = max_g;
				cmap[i][j].b = max_b;
			}
			else { // edge region
				cmap[i][j].r = min_r;
				cmap[i][j].g = min_g;
				cmap[i][j].b = min_b;
			}

		}
	}
	////////////////////////////////////////////

}


#define min3(x, y, z) ((x) < (y))? (((x) < (z))? (x) : (z)) : (((y) < (z))? (y) : (z))
#define max3(x, y, z) ((x) > (y))? (((x) > (z))? (x) : (z)) : (((y) > (z))? (y) : (z))

void GetWeickertUpwind(cimatrix& cmap, imatrix& sign, double h, int itr) 
// Use Inoue's upwind scheme to do MinMax
// using Tensor flow to get the Laplacian sign
// we find min and max in terms of R + G + B
// control both the half window size and iteration number
{
	int	i, j;
	double MAX_DOG = -1.0;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	cimatrix tmp2(image_x, image_y);
	
	////////////////////////////////////////////
	//double weight;
	int m;
	int x_r, x_l, y_t, y_b;
	int rx, gx, bx, ry, gy, by;
	//half1 = half;
	
	/////////////////////////////
	//int half1 = 1;
	//double h = 0.3; // step size (must be smaller than 0.5)
	/////////////////////////////////////////
	for (m = 0; m < itr; m++) {
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				////////////////////////////
				
				if (sign[i][j] == 255) {
					x_r = i+1; if (x_r > image_x-1) x_r = image_x-1;
					x_l = i-1; if (x_l < 0) x_l = 0;
					rx = max3(cmap[x_r][j].r - cmap[i][j].r, cmap[x_l][j].r - cmap[i][j].r, 0);
					gx = max3(cmap[x_r][j].g - cmap[i][j].g, cmap[x_l][j].g - cmap[i][j].g, 0);
					bx = max3(cmap[x_r][j].b - cmap[i][j].b, cmap[x_l][j].b - cmap[i][j].b, 0);

					y_t = j+1; if (y_t > image_y-1) y_t = image_y-1;
					y_b = j-1; if (y_b < 0) y_b = 0;
					ry = max3(cmap[i][y_t].r - cmap[i][j].r, cmap[i][y_b].r - cmap[i][j].r, 0);
					gy = max3(cmap[i][y_t].g - cmap[i][j].g, cmap[i][y_b].g - cmap[i][j].g, 0);
					by = max3(cmap[i][y_t].b - cmap[i][j].b, cmap[i][y_b].b - cmap[i][j].b, 0);
				
					tmp2[i][j].r = (GLubyte)( cmap[i][j].r + h * sqrt( (double) rx*rx + ry*ry ) );
					tmp2[i][j].g = (GLubyte)( cmap[i][j].g + h * sqrt( (double) gx*gx + gy*gy ) );
					tmp2[i][j].b = (GLubyte)( cmap[i][j].b + h * sqrt( (double) bx*bx + by*by ) );

				}
				else { // sign[i][j] == 0
					x_r = i+1; if (x_r > image_x-1) x_r = image_x-1;
					x_l = i-1; if (x_l < 0) x_l = 0;
					rx = min3(cmap[x_r][j].r - cmap[i][j].r, cmap[x_l][j].r - cmap[i][j].r, 0);
					gx = min3(cmap[x_r][j].g - cmap[i][j].g, cmap[x_l][j].g - cmap[i][j].g, 0);
					bx = min3(cmap[x_r][j].b - cmap[i][j].b, cmap[x_l][j].b - cmap[i][j].b, 0);

					y_t = j+1; if (y_t > image_y-1) y_t = image_y-1;
					y_b = j-1; if (y_b < 0) y_b = 0;
					ry = min3(cmap[i][y_t].r - cmap[i][j].r, cmap[i][y_b].r - cmap[i][j].r, 0);
					gy = min3(cmap[i][y_t].g - cmap[i][j].g, cmap[i][y_b].g - cmap[i][j].g, 0);
					by = min3(cmap[i][y_t].b - cmap[i][j].b, cmap[i][y_b].b - cmap[i][j].b, 0);
				
					tmp2[i][j].r = (GLubyte)( cmap[i][j].r - h * sqrt( (double) rx*rx + ry*ry ) );
					tmp2[i][j].g = (GLubyte)( cmap[i][j].g - h * sqrt( (double) gx*gx + gy*gy ) );
					tmp2[i][j].b = (GLubyte)( cmap[i][j].b - h * sqrt( (double) bx*bx + by*by ) );
				}
			}
				
		}
		cmap.copy(tmp2);
	}
	////////////////////////////////////////////
}

void GetColContrastEnhance(cimatrix& cmap, double pivot, double adjust) 
{
	int	i, j;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	double r, g, b;
	double h, s, v;
	double new_v;

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			r = cmap[i][j].r / 255.;
			g = cmap[i][j].g / 255.;
			b = cmap[i][j].b / 255.;
			//TRACE("r = %f\n", r);
			//TRACE("g = %f\n", g);
			//TRACE("b = %f\n", b);
			RGB2HSV(r, g, b, h, s, v);
			///////////////////////////////
			if (v < pivot) {
				new_v = (adjust / pivot) * v;
			}
			else {
				new_v = ( (1-adjust) / (1-pivot) ) * (v - 1.0) + 1.0;
			}
			if (new_v > 1.0) new_v = 1.0;
			if (new_v < 0.0) new_v = 0.0;
			///////////////////////////////
			HSV2RGB(h, s, new_v, r, g, b);
			cmap[i][j].r = (GLubyte)(r*255);
			cmap[i][j].g = (GLubyte)(g*255);
			cmap[i][j].b = (GLubyte)(b*255);
			//TRACE("r2 = %f\n", r);
			//TRACE("g2 = %f\n", g);
			//TRACE("b2 = %f\n", b);
		}			
	}
}

void GetGrayContrastEnhanceDouble(matrix& map, double pivot, double adjust) 
{
	int	i, j;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	double v;
	double new_v;

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			v = map[i][j];
			///////////////////////////////
			if (v < pivot) {
				new_v = (adjust / pivot) * v;
			}
			else {
				new_v = ( (1-adjust) / (1-pivot) ) * (v - 1.0) + 1.0;
			}
			if (new_v > 1.0) new_v = 1.0;
			if (new_v < 0.0) new_v = 0.0;
			///////////////////////////////
			map[i][j] = new_v;
		}			
	}
}

int GetPointFlowShock(int i, int j, imatrix& image, Field& gfield, 
							vector& GAU1, vector& GAU2, vector& GAU3, double tau)
// following the flow, compute the DOG
// but do not adjust GVF directions! Just follow it!
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s, i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	int dd;
	double val, c_val;

	int half_w1, half_w2, half_l;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;
	/////////////////////////////////////////////////
	// INSIDE and OUTSIDE FLOW DOG
	////////////////////////////////////
	//count = 0; // number of pixels traversed
	sum1 = sum2 = 0.0;
	w_sum1 = w_sum2 = 0.0;
	weight1 = weight2 = 0.0;
	////////////////////////////////////////////////
	// One half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	//////////////////////////
	c_val = image[i_x][i_y];
	weight1 = GAU1[0];
	weight1 *= GAU3[0];
	sum1 += c_val * weight1;
	w_sum1 += weight1;
	weight2 = GAU2[0];
	weight2 *= GAU3[0];
	sum2 += c_val * weight2;
	w_sum2 += weight2;
	////////////////////////////
	for (k = 0; k < half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = -gfield[i_x][i_y].gy;
		vt[1] = gfield[i_x][i_y].gx;
		vt.make_unit();
		////////////////////////////////////////
		for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s;
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			val = image[x1][y1];
			/////////////////////////////////////////////////////////
			dd = ABS(s);
			if (dd > half_w1) weight1 = 0.0;
			else weight1 = GAU1[dd];
			//////////////////////////////////
			// The following Gaussian smoothing along main axis is essential for good quality!
			weight1 *= GAU3[k]; 
			////////////////////
			sum1 += val * weight1;
			w_sum1 += weight1;
			//d3 = round( ABS(val - c_val) );
			/////////////////////////////////////////////////////
			weight2 = GAU2[dd];
			weight2 *= GAU3[k];
			sum2 += val * weight2;
			w_sum2 += weight2;
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
	}
	////////////////////////////////////////////////
	// Other half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	for (k = 0; k < half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = gfield[i_x][i_y].gy;
		vt[1] = -gfield[i_x][i_y].gx;
		vt.make_unit();
		////////////////////////////////////
		for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s;
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			val = image[x1][y1];
			/////////////////////////////////////////////////////////
			dd = ABS(s);
			if (dd > half_w1) weight1 = 0.0;
			else weight1 = GAU1[dd];
			weight1 *= GAU3[k];
			sum1 += val * weight1;
			w_sum1 += weight1;
			//d3 = round( ABS(val - c_val) );
			/////////////////////////////////////////////////////
			weight2 = GAU2[dd];
			weight2 *= GAU3[k];
			sum2 += val * weight2;
			w_sum2 += weight2;
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
	}
	sum1 /= w_sum1; // normalize
	sum2 /= w_sum2; // normalize
	////////////////////////////////////////////////////
	if (sum1 == 0.0 && sum2 == 0.0)
		sum1 = 1.0; // make it a non-edge
	//////////////////////////////////////
	if (sum1 - tau * sum2 >= 0) // non edge
		flow_DOG_sign = 1; // non edge
	else  // edge
		flow_DOG_sign = 0; // edge
	//if (sum1 - tau * sum2 > 0) // non edge
	//	flow_DOG = 1.0; 
	//else // edge!
	//	flow_DOG = 1.0 + tanh(sum1 - tau * sum2);

	return flow_DOG_sign;
	
}

double GetPointFlowShock2(int i, int j, imatrix& image, Field& gfield, 
							vector& GAU1, vector& GAU2, vector& GAU3, imatrix& sign, double z)
// following the flow, compute the DOG
// but do not adjust GVF directions! Just follow it!
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s, i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	double val, c_val;

	int half_w1, half_w2, half_l;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	double shock_val = 0.0; // flow-based gradient magnitude

	///////////////////////////
	deque<int> vec;
	vec.clear();

	int c_sign;

	if (sign[i][j] == 1) c_sign = 1;
	else c_sign = 0;
		
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;
	/////////////////////////////////////////////////
	// INSIDE and OUTSIDE FLOW DOG
	////////////////////////////////////
	//count = 0; // number of pixels traversed
	sum1 = sum2 = 0.0;
	w_sum1 = w_sum2 = 0.0;
	weight1 = weight2 = 0.0;
	////////////////////////////////////////////////
	// One half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	//////////////////////////
	c_val = image[i_x][i_y];
	weight1 = GAU1[0];
	weight1 *= GAU3[0];
	sum1 += c_val * weight1;
	w_sum1 += weight1;
	weight2 = GAU2[0];
	weight2 *= GAU3[0];
	sum2 += c_val * weight2;
	w_sum2 += weight2;
	////////////////////////////
	for (k = 0; k < half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = -gfield[i_x][i_y].gy;
		vt[1] = gfield[i_x][i_y].gx;
		vt.make_unit();
		////////////////////////////////////////
		for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s;
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			val = image[x1][y1];
			/////////////////////////////////////////////////////////
			if (c_sign) { // non edge region
                if (sign[x1][y1] == 1) vec.push_back(image[x1][y1]);
			}
			else { // edge region
				if (sign[x1][y1] == 0) vec.push_back(image[x1][y1]);
			}
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
	}
	////////////////////////////////////////////////
	// Other half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	for (k = 0; k < half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = gfield[i_x][i_y].gy;
		vt[1] = -gfield[i_x][i_y].gx;
		vt.make_unit();
		////////////////////////////////////
		for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s;
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			val = image[x1][y1];
			/////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////
			if (c_sign) { // non edge region
                if (sign[x1][y1] == 1) vec.push_back(image[x1][y1]);
			}
			else { // edge region
				if (sign[x1][y1] == 0) vec.push_back(image[x1][y1]);
			}
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
	}
	//sum1 /= w_sum1; // normalize
	//sum2 /= w_sum2; // normalize
	sort(vec.begin(), vec.end()); 
	//flow_median = vec[vec.size()/2]; 
	//double z = 0.3;
	if (vec.size() == 0)
		return round(c_val); 
	if (c_sign) // non edge region
		shock_val = vec[vec.size()-1]; 
		//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
	else // edge region
		//image2[i][j] = vec[0]; 
		shock_val = round( (1-z) * vec[0] + z * 0 ); // darkest value
	////////////////////////////////////////////////////
	return shock_val;
}

void GetFlowShockDoG(imatrix& image, Field& gfield, matrix& dog, vector& GAU1, vector& GAU2, double tau)
// For each pixel, compute the directional DOG
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s;
	int x1, y1;
	int i, j;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	int dd;
	double val;

	int half_w1, half_w2;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;

	int image_x, image_y;

	image_x = image.getRow();
	image_y = image.getCol();
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	/////////////////////////////////////////////////
	// INSIDE and OUTSIDE FLOW DOG
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			sum1 = sum2 = 0.0;
			w_sum1 = w_sum2 = 0.0;
			weight1 = weight2 = 0.0;
	
			vn[0] = gfield[i][j].gx;
			vn[1] = gfield[i][j].gy;
			if (vn[0] == 0.0 && vn[1] == 0.0) {
				sum1 = 1.0;
				sum2 = 1.0;
				dog[i][j] = sum1 - tau * sum2;
				continue;
			}
			vn.make_unit();
			d_x = i; d_y = j;
			////////////////////////////////////////
			for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
				////////////////////////
				x = d_x + vn[0] * s;
				y = d_y + vn[1] * s;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					continue;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				val = image[x1][y1];
				/////////////////////////////////////////////////////////
				dd = ABS(s);
				if (dd > half_w1) weight1 = 0.0;
				else weight1 = GAU1[dd];
				//////////////////////////////////
				sum1 += val * weight1;
				w_sum1 += weight1;
				//d3 = round( ABS(val - c_val) );
				/////////////////////////////////////////////////////
				weight2 = GAU2[dd];
				sum2 += val * weight2;
				w_sum2 += weight2;
			}
			/////////////////////////
			sum1 /= w_sum1; // normalize
			sum2 /= w_sum2; // normalize
			////////////////////////////////////////////////////
			//////////////////////////////////////
			dog[i][j] = sum1 - tau * sum2;
			//if (sum1 - tau * sum2 > 0) dog[i][j] = 1.0;
			//else dog[i][j] = 1.0 + tanh(sum1 - tau * sum2);
		}
	}

}

void GetFlowShockDoGETF(imatrix& image, ETF& e, matrix& dog, vector& GAU1, vector& GAU2, double tau)
// For each pixel, compute the directional DOG from ETF
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s;
	int x1, y1;
	int i, j;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	int dd;
	double val;

	int half_w1, half_w2;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;

	int image_x, image_y;

	image_x = image.getRow();
	image_y = image.getCol();
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	/////////////////////////////////////////////////
	// INSIDE and OUTSIDE FLOW DOG
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			sum1 = sum2 = 0.0;
			w_sum1 = w_sum2 = 0.0;
			weight1 = weight2 = 0.0;
	
			vn[0] = -e[i][j].ty;
			vn[1] = e[i][j].tx;
			if (vn[0] == 0.0 && vn[1] == 0.0) {
				sum1 = 1.0;
				sum2 = 1.0;
				dog[i][j] = sum1 - tau * sum2;
				continue;
			}
			//vn.make_unit();
			d_x = i; d_y = j;
			////////////////////////////////////////
			for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
				////////////////////////
				x = d_x + vn[0] * s;
				y = d_y + vn[1] * s;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					continue;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				val = image[x1][y1];
				/////////////////////////////////////////////////////////
				dd = ABS(s);
				if (dd > half_w1) weight1 = 0.0;
				else weight1 = GAU1[dd];
				//////////////////////////////////
				sum1 += val * weight1;
				w_sum1 += weight1;
				//d3 = round( ABS(val - c_val) );
				/////////////////////////////////////////////////////
				weight2 = GAU2[dd];
				sum2 += val * weight2;
				w_sum2 += weight2;
			}
			/////////////////////////
			sum1 /= w_sum1; // normalize
			sum2 /= w_sum2; // normalize
			////////////////////////////////////////////////////
			//////////////////////////////////////
			dog[i][j] = sum1 - tau * sum2;
			//TRACE("sum1[%d][%d] = %f\n", i, j, sum1);
			//TRACE("sum2[%d][%d] = %f\n", i, j, sum2);
			//TRACE("tau = %f\n", tau);
			//TRACE("sum1 - tau * sum2 = %f\n", sum1 - tau * sum2);
			//TRACE("dog[%d][%d] = %f\n", i, j, dog[i][j]);
			//if (sum1 - tau * sum2 > 0) dog[i][j] = 1.0;
			//else dog[i][j] = 1.0 + tanh(sum1 - tau * sum2);
		}
	}

}

void GetFlowShockDoGETF2(imatrix& image, ETF& e, imatrix& sign, double sigma, double tau)
// For each pixel, compute the directional DOG from ETF
// Use sigma, and sign
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s;
	int x1, y1;
	int i, j;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	int dd;
	double val;

	int half_w1, half_w2;

	vector GAU1, GAU2;
	MakeGaussVector(sigma, GAU1);
	MakeGaussVector(sigma*1.6, GAU2);

	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;

	int image_x, image_y;

	image_x = image.getRow();
	image_y = image.getCol();
	
	//int flow_DOG_sign = 0; 
	double dog;
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	/////////////////////////////////////////////////
	// INSIDE and OUTSIDE FLOW DOG
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			sum1 = sum2 = 0.0;
			w_sum1 = w_sum2 = 0.0;
			weight1 = weight2 = 0.0;
	
			vn[0] = -e[i][j].ty;
			vn[1] = e[i][j].tx;
			if (vn[0] == 0.0 && vn[1] == 0.0) {
				//sum1 = 1.0;
				//sum2 = 1.0;
				//dog = sum1 - tau * sum2;
				//continue;
				vn[0] = 1.0; vn[1] = 0.0; // assign horizontal unit vector
			}
			//vn.make_unit();
			d_x = i; d_y = j;
			////////////////////////////////////////
			for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
				////////////////////////
				x = d_x + vn[0] * s;
				y = d_y + vn[1] * s;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					continue;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				val = image[x1][y1];
				/////////////////////////////////////////////////////////
				dd = ABS(s);
				if (dd > half_w1) weight1 = 0.0;
				else weight1 = GAU1[dd];
				//////////////////////////////////
				sum1 += val * weight1;
				w_sum1 += weight1;
				//d3 = round( ABS(val - c_val) );
				/////////////////////////////////////////////////////
				weight2 = GAU2[dd];
				sum2 += val * weight2;
				w_sum2 += weight2;
			}
			/////////////////////////
			sum1 /= w_sum1; // normalize
			sum2 /= w_sum2; // normalize
			////////////////////////////////////////////////////
			//////////////////////////////////////
			dog = sum1 - tau * sum2;
			//TRACE("sum1[%d][%d] = %f\n", i, j, sum1);
			//TRACE("sum2[%d][%d] = %f\n", i, j, sum2);
			//TRACE("tau = %f\n", tau);
			//TRACE("sum1 - tau * sum2 = %f\n", sum1 - tau * sum2);
			//TRACE("dog[%d][%d] = %f\n", i, j, dog[i][j]);
			if (dog > 0) sign[i][j] = 255;
			else sign[i][j] = 0;
			//else dog[i][j] = 1.0 + tanh(sum1 - tau * sum2);
		}
	}

}


void GetFlowShockDoGETF3(imatrix& image, ETF& e, imatrix& sign)
// For each pixel, compute the directional second derivative from ETF
// Use sigma, and sign
{
	vector vn(2);
	double d_x, d_y;
	double x_r, y_r, x_l, y_l;
	double val_r, val_l, val;

	int x1, y1;
	int i, j;

	int image_x, image_y;

	image_x = image.getRow();
	image_y = image.getCol();
	
	//int flow_DOG_sign = 0; 
	double dog;
	
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			vn[0] = -e[i][j].ty;
			vn[1] = e[i][j].tx;
			if (vn[0] == 0.0 && vn[1] == 0.0) {
				vn[0] = 1.0; vn[1] = 0.0; // assign horizontal unit vector
			}
			//vn.make_unit();
			d_x = (double)i; d_y = (double)j;

			val = (double)image[i][j];

			x_r = d_x + vn[0];
			y_r = d_y + vn[1];
			x1 = round(x_r);	if (x1 < 0) x1 = 0; if (x1 > image_x-1) x1 = image_x-1;
			y1 = round(y_r);	if (y1 < 0) y1 = 0; if (y1 > image_y-1) y1 = image_y-1;
			val_r = (double)image[x1][y1];
			
			x_l = d_x - vn[0];
			y_l = d_y - vn[1];
			x1 = round(x_l);	if (x1 < 0) x1 = 0; if (x1 > image_x-1) x1 = image_x-1;
			y1 = round(y_l);	if (y1 < 0) y1 = 0; if (y1 > image_y-1) y1 = image_y-1;
			val_l = (double)image[x1][y1];
			
			//////////////////////////////////////
			dog = (val_r - val) - (val - val_l);
			if (dog >= 0) sign[i][j] = 255;
			else sign[i][j] = 0;
		}
	}

}



void GetLaplacian5(imatrix& image, double sigma, imatrix& lap, double tau, double thres)
// For each pixel, compute Laplacian value (using DoG)
// Faster version!
// use tau
{
	int	i, j;
	double dog;

	int image_x = image.getRow();
	int image_y = image.getCol();

	//vector GAU1, GAU2;
	//MakeGaussVector(sigma, GAU1);
	//MakeGaussVector(sigma*1.6, GAU2);
	//half1 = GAU1.getMax()-1;
	//half2 = GAU2.getMax()-1;

	matrix inner(image_x, image_y);
	//inner.zero();
	//GaussSmoothSep(inner, sigma);
	GaussSmoothSepDouble(image, sigma, inner);

	matrix outer(image_x, image_y);
	//outer.zero();
	GaussSmoothSepDouble(image, sigma*1.6, outer);

	//deque<int> vec[3];
	////////////////////////////////////////////////////////
	//matrix tmp(image_x, image_y);
	//cimatrix tmp2(image_x, image_y);
	//tmp.zero();
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			dog = inner[i][j] - tau * outer[i][j];
			if (dog >= 0) lap[i][j] = 255; 
			else if (tanh(dog) + 1 > thres) lap[i][j] = 255; 
			else lap[i][j] = 0; 
		}
	}

}

void GetLaplacian6(imatrix& image, double sigma, imatrix& lap, double tau)
// For each pixel, compute Laplacian value (using DoG)
// Faster version!
// use tau
{
	int	i, j;
	double dog;

	int image_x = image.getRow();
	int image_y = image.getCol();

	//vector GAU1, GAU2;
	//MakeGaussVector(sigma, GAU1);
	//MakeGaussVector(sigma*1.6, GAU2);
	//half1 = GAU1.getMax()-1;
	//half2 = GAU2.getMax()-1;

	matrix inner(image_x, image_y);
	//inner.zero();
	//GaussSmoothSep(inner, sigma);
	GaussSmoothSepDouble(image, sigma, inner);

	matrix outer(image_x, image_y);
	//outer.zero();
	GaussSmoothSepDouble(image, sigma*1.6, outer);

	//deque<int> vec[3];
	////////////////////////////////////////////////////////
	//matrix tmp(image_x, image_y);
	//cimatrix tmp2(image_x, image_y);
	//tmp.zero();
	//double thres = -0.000001;
	double thres;
	//thres = -0.0;
	//thres = -0.000001;
	thres = 0.5;
	//////////////////////////////
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			dog = inner[i][j] - tau * outer[i][j];
			if (dog >= 0) lap[i][j] = 255; 
			//if (dog >= thres) lap[i][j] = 255; 
			else if (tanh(dog) + 1 > thres) lap[i][j] = 255; 
			//else if (tanh(dog) + 1 > 0.1) lap[i][j] = 255; 
			else lap[i][j] = 0; 
		}
	}

}


void GetLaplacianDouble(imatrix& image, double sigma, matrix& G_map, double tau, double thres)
// For each pixel, compute Laplacian value (using DoG)
// Faster version!
// use tau
{
	int	i, j;
	double dog;

	int image_x = image.getRow();
	int image_y = image.getCol();

	//vector GAU1, GAU2;
	//MakeGaussVector(sigma, GAU1);
	//MakeGaussVector(sigma*1.6, GAU2);
	//half1 = GAU1.getMax()-1;
	//half2 = GAU2.getMax()-1;

	matrix inner(image_x, image_y);
	//inner.zero();
	//GaussSmoothSep(inner, sigma);
	GaussSmoothSepDouble(image, sigma, inner);

	matrix outer(image_x, image_y);
	//outer.zero();
	GaussSmoothSepDouble(image, sigma*1.6, outer);

	//////////////////////////////
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			dog = inner[i][j] - tau * outer[i][j];
			if (dog >= 0) G_map[i][j] = 1.0; 
			else G_map[i][j] = tanh(dog) + 1; 
		}
	}
}

void GetDogSep(imatrix& image, double sigma, imatrix& lap, double tau)
// For each pixel, compute DoG
// Faster version!
// use tau
{
	int	i, j;
	double dog;

	int image_x = image.getRow();
	int image_y = image.getCol();

	matrix inner(image_x, image_y);
	//inner.zero();
	//GaussSmoothSep(inner, sigma);
	GaussSmoothSepDouble(image, sigma, inner);

	matrix outer(image_x, image_y);
	//outer.zero();
	GaussSmoothSepDouble(image, sigma*1.6, outer);

	//deque<int> vec[3];
	////////////////////////////////////////////////////////
	//matrix tmp(image_x, image_y);
	//cimatrix tmp2(image_x, image_y);
	//tmp.zero();
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			dog = inner[i][j] - tau * outer[i][j];
			if (dog > 0) lap[i][j] = 255; 
			else lap[i][j] = round( (tanh(dog) + 1.0) * 255.0 ); 
		}
	}
}

void GetDogSepDouble(imatrix& image, double sigma, mymatrix& dog, double tau)
// used for DOGLIC!
// For each pixel, compute DoG
// Faster version!
// use tau
{
	int	i, j;
	//double dog;

	int image_x = image.getRow();
	int image_y = image.getCol();

	matrix inner(image_x, image_y);
	//inner.zero();
	//GaussSmoothSep(inner, sigma);
	GaussSmoothSepDouble(image, sigma, inner);

	matrix outer(image_x, image_y);
	//outer.zero();
	GaussSmoothSepDouble(image, sigma*1.6, outer);

	//deque<int> vec[3];
	////////////////////////////////////////////////////////
	//matrix tmp(image_x, image_y);
	//cimatrix tmp2(image_x, image_y);
	//tmp.zero();
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			dog[i][j] = inner[i][j] - tau * outer[i][j];
			//if (dog > 0) lap[i][j] = 255; 
			//else lap[i][j] = round( (tanh(dog) + 1.0) * 255.0 ); 
		}
	}
}

void GetDirectionalDogSepDouble(imatrix& image, double sigma, mymatrix& dog, double tau, ETF& e, myvec& GAU1, myvec& GAU2)
// used for DOGLIC!
// For each pixel, compute DoG
// Faster version!
// use tau
{
	int	i, j;
	int image_x = image.getRow();
	int image_y = image.getCol();
	matrix inner(image_x, image_y);
	GaussSmoothSepDouble(image, sigma, inner);
	matrix outer(image_x, image_y);
	GaussSmoothSepDouble(image, sigma*1.6, outer);


	myvec vn(2);
	double x, y, d_x, d_y;
	double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s;
	int x1, y1;
	int dd;
	double val;
	int half_w1, half_w2;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;
	int flow_DOG_sign = 0; 
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			/*dog[i][j] = 0;
			if(lambda2 > 0) {
				dog[i][j] = lambda2 * (inner[i][j] - tau * outer[i][j]);
			}*/
//			if(lambda > 0) {
				sum1 = sum2 = 0.0;
				w_sum1 = w_sum2 = 0.0;
				weight1 = weight2 = 0.0;
		
				vn[0] = -e[i][j].ty;
				vn[1] = e[i][j].tx;

				if (vn[0] == 0.0 && vn[1] == 0.0) {
					sum1 = 255.0;
					sum2 = 255.0;
					dog[i][j] = sum1 - tau * sum2;
					continue;
				}
				d_x = i; d_y = j;
				////////////////////////////////////////
				for (s = -half_w2; s <= half_w2; s++) { 
					////////////////////////
					x = d_x + vn[0] * s;
					y = d_y + vn[1] * s;
					/////////////////////////////////////////////////////
					if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
						continue;
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					val = image[x1][y1];
					/////////////////////////////////////////////////////////
					dd = ABS(s);
					if (dd > half_w1) weight1 = 0.0;
					else weight1 = GAU1[dd];
					//////////////////////////////////
					sum1 += val * weight1;
					w_sum1 += weight1;
					/////////////////////////////////////////////////////
					weight2 = GAU2[dd];
					sum2 += val * weight2;
					w_sum2 += weight2;
				}
				/////////////////////////
				sum1 /= w_sum1;
				sum2 /= w_sum2;
				//////////////////////////////////////
				dog[i][j] = (1 - e[i][j].mag) * (sum1 - tau * sum2);
//			}
		}
	}
}

void GetDirectionalDogSepDouble2(imatrix& image, double sigma, mymatrix& dog, 
								double tau, ETF& e, myvec& GAU1, myvec& GAU2, imatrix& edge)
// used for DOGLIC!
// For each pixel, compute DoG
// Faster version!
// use tau
{
	int	i, j;
	int image_x = image.getRow();
	int image_y = image.getCol();
	matrix inner(image_x, image_y);
	GaussSmoothSepDouble(image, sigma, inner);
	matrix outer(image_x, image_y);
	GaussSmoothSepDouble(image, sigma*1.6, outer);


	myvec vn(2);
	double x, y, d_x, d_y;
	double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s;
	int x1, y1;
	int dd;
	double val;
	int half_w1, half_w2;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;
	int flow_DOG_sign = 0; 
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			/*dog[i][j] = 0;
			if(lambda2 > 0) {
				dog[i][j] = lambda2 * (inner[i][j] - tau * outer[i][j]);
			}*/
			//double lambda = ((double) edge[i][j] / 255.0);
			//double lambda2 = 1.0 - lambda;
			//if(edge[i][j] >= 254) {
				sum1 = sum2 = 0.0;
				w_sum1 = w_sum2 = 0.0;
				weight1 = weight2 = 0.0;
		
				vn[0] = -e[i][j].ty;
				vn[1] = e[i][j].tx;

				if (vn[0] == 0.0 && vn[1] == 0.0) {
					sum1 = 255.0;
					sum2 = 255.0;
					dog[i][j] = sum1 - tau * sum2;
					continue;
				}
				d_x = i; d_y = j;
				////////////////////////////////////////
				for (s = -half_w2; s <= half_w2; s++) { 
					////////////////////////
					x = d_x + vn[0] * s;
					y = d_y + vn[1] * s;
					/////////////////////////////////////////////////////
					if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
						continue;
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					val = image[x1][y1];
					/////////////////////////////////////////////////////////
					dd = ABS(s);
					if (dd > half_w1) weight1 = 0.0;
					else weight1 = GAU1[dd];
					//////////////////////////////////
					sum1 += val * weight1;
					w_sum1 += weight1;
					/////////////////////////////////////////////////////
					weight2 = GAU2[dd];
					sum2 += val * weight2;
					w_sum2 += weight2;
				}
				/////////////////////////
				sum1 /= w_sum1;
				sum2 /= w_sum2;
				//////////////////////////////////////
				dog[i][j] = (sum1 - tau * sum2);
			//}
		}
	}
}

void GetDogSepDouble2(imatrix& image, double sigma, matrix& dog, double tau)
// Faster version!
// normalized to [0,1]
{
	int	i, j;
	//double dog;

	int image_x = image.getRow();
	int image_y = image.getCol();

	matrix inner(image_x, image_y);
	//inner.zero();
	//GaussSmoothSep(inner, sigma);
	GaussSmoothSepDouble(image, sigma, inner);

	matrix outer(image_x, image_y);
	//outer.zero();
	GaussSmoothSepDouble(image, sigma*1.6, outer);

	//deque<int> vec[3];
	////////////////////////////////////////////////////////
	//matrix tmp(image_x, image_y);
	//cimatrix tmp2(image_x, image_y);
	//tmp.zero();
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			dog[i][j] = inner[i][j] - tau * outer[i][j];
			if (dog[i][j] >= 0) dog[i][j] = 1.0; 
			else dog[i][j] = tanh(dog[i][j]) + 1.0; 
		}
	}
}


void GetFlowShockInfluence(Field& gfield, matrix& dog, imatrix& sign, vector& GAU3)
// following the flow, compute the DOG
// FASTER version!
{
	vector vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, w_sum1, sum1;

	int i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	double val;

	int image_x = sign.getRow();
	int image_y = sign.getCol();

	int half_l;
	//half_w1 = GAU1.getMax()-1;
	//half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;

	int i, j;

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
	
			/////////////////////////////////////////////////
			// INSIDE and OUTSIDE FLOW DOG
			////////////////////////////////////
			//count = 0; // number of pixels traversed
			sum1 = 0.0;
			w_sum1 = 0.0;
			weight1 = 0.0;
			/////////////////////////////////
			val = dog[i][j];
			weight1 = GAU3[0]; 
			//weight1 = 1.0; 
			sum1 = val * weight1;
			w_sum1 += weight1;
			////////////////////////////////////////////////
			// One half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			////////////////////////////
			for (k = 0; k < half_l; k++) {
				vt[0] = -gfield[i_x][i_y].gy;
				vt[1] = gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				val = dog[x1][y1];
				//////////////////////////////
				weight1 = GAU3[k]; 
				//weight1 = 1.0; // uniform weight
				////////////////////
				sum1 += val * weight1;
				w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////////////
			// Other half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			for (k = 0; k < half_l; k++) {
				vt[0] = gfield[i_x][i_y].gy;
				vt[1] = -gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				val = dog[x1][y1];
				//////////////////////////////
				weight1 = GAU3[k]; 
				//weight1 = 1.0; // uniform weight
				////////////////////
				sum1 += val * weight1;
				w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////
			sum1 /= w_sum1; // normalize
			//////////////////////////////////////
			if (sum1 >= 0) // non edge
				sign[i][j] = 1; // non edge
			else  // edge
				sign[i][j] = 0; // edge
			//if (sum1 - tau * sum2 > 0) // non edge
			//	flow_DOG = 1.0; 
			//else // edge!
			//	flow_DOG = 1.0 + tanh(sum1 - tau * sum2);
		}
	}
	
}


void GetFlowShockInfluenceETF(ETF& e, matrix& dog, imatrix& sign, vector& GAU3)
// following the flow, compute the DOG
// FASTER version!
{
	vector vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, w_sum1, sum1;

	int i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	double val;

	int image_x = sign.getRow();
	int image_y = sign.getCol();

	int half_l;
	//half_w1 = GAU1.getMax()-1;
	//half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;

	int i, j;

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
	
			/////////////////////////////////////////////////
			// INSIDE and OUTSIDE FLOW DOG
			////////////////////////////////////
			//count = 0; // number of pixels traversed
			sum1 = 0.0;
			w_sum1 = 0.0;
			weight1 = 0.0;
			/////////////////////////////////
			val = dog[i][j];
			weight1 = GAU3[0]; 
			//weight1 = 1.0; 
			sum1 = val * weight1;
			w_sum1 += weight1;
			////////////////////////////////////////////////
			// One half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			////////////////////////////
			for (k = 0; k < half_l; k++) {
				vt[0] = e[i_x][i_y].tx;
				vt[1] = e[i_x][i_y].ty;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				val = dog[x1][y1];
				//////////////////////////////
				weight1 = GAU3[k]; 
				//weight1 = 1.0; // uniform weight
				////////////////////
				sum1 += val * weight1;
				w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////////////
			// Other half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			for (k = 0; k < half_l; k++) {
				vt[0] = -e[i_x][i_y].tx;
				vt[1] = -e[i_x][i_y].ty;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				val = dog[x1][y1];
				//////////////////////////////
				weight1 = GAU3[k]; 
				//weight1 = 1.0; // uniform weight
				////////////////////
				sum1 += val * weight1;
				w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////
			sum1 /= w_sum1; // normalize
			//////////////////////////////////////
			if (sum1 >= 0) // non edge
				sign[i][j] = 1; // non edge
			else  // edge
				sign[i][j] = 0; // edge
			//if (sum1 - tau * sum2 > 0) // non edge
			//	flow_DOG = 1.0; 
			//else // edge!
			//	flow_DOG = 1.0 + tanh(sum1 - tau * sum2);
		}
	}
	
}

void GetFlowDoGMainAxis(Field& gfield, matrix& dog, matrix& tmp, vector& GAU3)
// following the flow, get the Gaussian of DOG values
{
	vector vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, w_sum1, sum1;

	int i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	double val;

	int image_x = dog.getRow();
	int image_y = dog.getCol();

	int half_l;
	//half_w1 = GAU1.getMax()-1;
	//half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;

	int i, j;

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
	
			/////////////////////////////////////////////////
			// INSIDE and OUTSIDE FLOW DOG
			////////////////////////////////////
			//count = 0; // number of pixels traversed
			sum1 = 0.0;
			w_sum1 = 0.0;
			weight1 = 0.0;
			/////////////////////////////////
			val = dog[i][j];
			weight1 = GAU3[0]; 
			//weight1 = 1.0; 
			sum1 = val * weight1;
			w_sum1 += weight1;
			////////////////////////////////////////////////
			// One half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			////////////////////////////
			for (k = 0; k < half_l; k++) {
				vt[0] = -gfield[i_x][i_y].gy;
				vt[1] = gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				val = dog[x1][y1];
				//////////////////////////////
				weight1 = GAU3[k]; 
				//weight1 = 1.0; // uniform weight
				////////////////////
				sum1 += val * weight1;
				w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////////////
			// Other half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			for (k = 0; k < half_l; k++) {
				vt[0] = gfield[i_x][i_y].gy;
				vt[1] = -gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				val = dog[x1][y1];
				//////////////////////////////
				weight1 = GAU3[k]; 
				//weight1 = 1.0; // uniform weight
				////////////////////
				sum1 += val * weight1;
				w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////
			sum1 /= w_sum1; // normalize
			//////////////////////////////////////
			//if (sum1 == 0.0 && sum2 == 0.0)
			//	sum1 = 1.0; // make it a non-edge
			//////////////////////////////////////
			if (sum1 > 0) // non edge
				tmp[i][j] = 1.0; 
			else // edge!
				tmp[i][j] = 1.0 + tanh(sum1);
			
		}
	}
	
}

void GetFlowShockDoG2(imatrix& image, Field& gfield, matrix& dog, vector& GAU1, vector& GAU2, double tau)
// For each pixel, compute the directional DOG
// Optimized!
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s;
	int x1, y1;
	int i, j;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	int dd;
	double val;

	int half_w1, half_w2;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;

	int image_x, image_y;

	image_x = image.getRow();
	image_y = image.getCol();
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	//half_w1 = half_w2 = 25; // for TVCG

	/////////////////////////////////////////////////
	// INSIDE and OUTSIDE FLOW DOG
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			sum1 = sum2 = 0.0;
			w_sum1 = w_sum2 = 0.0;
			weight1 = weight2 = 0.0;
	
			vn[0] = gfield[i][j].gx;
			vn[1] = gfield[i][j].gy;
			if (vn[0] == 0.0 && vn[1] == 0.0) {
				sum1 = 1.0;
				sum2 = 1.0;
				dog[i][j] = sum1 - tau * sum2;
				continue;
			}
			//vn.make_unit();
			d_x = i; d_y = j;
			////////////////////////////////////////
			for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
				////////////////////////
				x = d_x + vn[0] * s;
				y = d_y + vn[1] * s;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					continue;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				val = image[x1][y1];
				/////////////////////////////////////////////////////////
				dd = ABS(s);
				if (dd > half_w1) weight1 = 0.0;
				else weight1 = GAU1[dd];
				//////////////////////////////////
				sum1 += val * weight1;
				w_sum1 += weight1;
				//d3 = round( ABS(val - c_val) );
				/////////////////////////////////////////////////////
				weight2 = GAU2[dd];
				sum2 += val * weight2;
				w_sum2 += weight2;
			}
			/////////////////////////
			sum1 /= w_sum1; // normalize
			sum2 /= w_sum2; // normalize
			////////////////////////////////////////////////////
			//////////////////////////////////////
			dog[i][j] = sum1 - tau * sum2;
			//if (sum1 - tau * sum2 > 0) dog[i][j] = 1.0;
			//else dog[i][j] = 1.0 + tanh(sum1 - tau * sum2);
		}
	}

}

void GetFlowDoGMainAxis2(Field& gfield, matrix& dog, matrix& tmp, vector& GAU3)
// following the flow, get the Gaussian of DOG values
// optimized!
{
	vector vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, w_sum1, sum1;

	int i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	double val;
	int i, j;

	int image_x = dog.getRow();
	int image_y = dog.getCol();

	int half_l;
	//half_w1 = GAU1.getMax()-1;
	//half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;
	
	//half_l = 25; // for TVCG

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
	
			/////////////////////////////////////////////////
			// INSIDE and OUTSIDE FLOW DOG
			////////////////////////////////////
			//count = 0; // number of pixels traversed
			sum1 = 0.0;
			w_sum1 = 0.0;
			weight1 = 0.0;
			/////////////////////////////////
			val = dog[i][j];
			weight1 = GAU3[0]; 
			//weight1 = 1.0; 
			sum1 = val * weight1;
			w_sum1 += weight1;
			////////////////////////////////////////////////
			// One half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			////////////////////////////
			for (k = 0; k < half_l; k++) {
				vt[0] = -gfield[i_x][i_y].gy;
				vt[1] = gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				//vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				val = dog[x1][y1];
				//////////////////////////////
				weight1 = GAU3[k]; 
				//weight1 = 1.0; // uniform weight
				////////////////////
				sum1 += val * weight1;
				w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////////////
			// Other half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			for (k = 0; k < half_l; k++) {
				vt[0] = gfield[i_x][i_y].gy;
				vt[1] = -gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				//vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				val = dog[x1][y1];
				//////////////////////////////
				weight1 = GAU3[k]; 
				//weight1 = 1.0; // uniform weight
				////////////////////
				sum1 += val * weight1;
				w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////
			sum1 /= w_sum1; // normalize
			//////////////////////////////////////
			//if (sum1 == 0.0 && sum2 == 0.0)
			//	sum1 = 1.0; // make it a non-edge
			//////////////////////////////////////
			if (sum1 > 0) // non edge
				tmp[i][j] = 1.0; 
			else // edge!
				tmp[i][j] = 1.0 + tanh(sum1);
			
		}
	}
	
}



void Get1DFlowColShockMinMax(cimatrix& cmap, Field& gfield, cimatrix& tmp, imatrix& sign, vector& GAU3, double z)
// following the flow, compute the MinMax
// Faster version!
{
	vector vt(2), vt_old(2);
	double x, y, d_x, d_y;
	//double weight1, w_sum1, sum1;
	//double r, g, b, h, s, v;

	int i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	//double val;
	int c_sign;
	deque<int> vec[3];
	int R, G, B;
	int inc_col, dec_col;

	int image_x = sign.getRow();
	int image_y = sign.getCol();

	int half_l;
	//half_w1 = GAU1.getMax()-1;
	//half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;

	int i, j;

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			//sum1 = 0.0;
			//w_sum1 = 0.0;
			//weight1 = 0.0;
			if (sign[i][j] == 1) c_sign = 1;
			else c_sign = 0;
			vec[0].clear();
			vec[1].clear();
			vec[2].clear();
			vec[0].push_back(cmap[i][j].r);
			vec[1].push_back(cmap[i][j].g);
			vec[2].push_back(cmap[i][j].b);
			/////////////////////////////////
			//val = dog[i][j];
			//weight1 = GAU3[0]; 
			//sum1 = val * weight1;
			//w_sum1 += weight1;
			////////////////////////////////////////////////
			// One half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			////////////////////////////
			for (k = 0; k < half_l; k++) {
				vt[0] = -gfield[i_x][i_y].gy;
				vt[1] = gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				//val = dog[x1][y1];
				//////////////////////////////
				//weight1 *= GAU3[k]; 
				////////////////////
				//sum1 += val * weight1;
				//w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////////////
			// Other half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			for (k = 0; k < half_l; k++) {
				vt[0] = gfield[i_x][i_y].gy;
				vt[1] = -gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				//val = dog[x1][y1];
				//////////////////////////////
				//weight1 *= GAU3[k]; 
				////////////////////
				//sum1 += val * weight1;
				//w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			//sum1 /= w_sum1; // normalize
			//////////////////////////////////////
			//if (sum1 >= 0) // non edge
			//	sign[i][j] = 1; // non edge
			//else  // edge
			//	sign[i][j] = 0; // edge
			sort(vec[0].begin(), vec[0].end()); 
			sort(vec[1].begin(), vec[1].end()); 
			sort(vec[2].begin(), vec[2].end()); 
			//flow_median = vec[vec.size()/2]; 
			//double z = 0.3;
			//if (vec[0].size() == 0 || vec[1].size() == 0 || vec[2].size() == 0) {
			//	shock_val.r = c_val_r;
			//	shock_val.g = c_val_g;
			//	shock_val.b = c_val_b;
			//	return shock_val; 
			//}
			if (c_sign) { // non edge region
				//tmp[i][j].r = vec[0][vec[0].size()-1]; 
				//tmp[i][j].g = vec[1][vec[1].size()-1]; 
				//tmp[i][j].b = vec[2][vec[2].size()-1]; 
				tmp[i][j].r = vec[0][vec[0].size()/2]; 
				tmp[i][j].g = vec[1][vec[1].size()/2]; 
				tmp[i][j].b = vec[2][vec[2].size()/2]; 
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
				//tmp[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
				//tmp[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
				//////////////////////////////////////////////////
				inc_col = 7;
				R = tmp[i][j].r; G = tmp[i][j].g; B = tmp[i][j].b;
				R += inc_col; if (R > 255) R = 255;
				G += inc_col; if (G > 255) G = 255;
				B += inc_col; if (B > 255) B = 255;
				tmp[i][j].r = (GLubyte)R;
				tmp[i][j].g = (GLubyte)G;
				tmp[i][j].b = (GLubyte)B;
				//////////////////////////////////////////
				/*
				r = tmp[i][j].r / 255.;
				g = tmp[i][j].g / 255.;
				b = tmp[i][j].b / 255.;
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 1.0;
				v += 0.05; if (v > 1.0) v = 1.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
				*/
			}
			else {// edge region
				tmp[i][j].r = vec[0][0]; // darkest value
				tmp[i][j].g = vec[1][0];
				tmp[i][j].b = vec[2][0];
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
				//////////////////////////////////////
				dec_col = 2;
				R = tmp[i][j].r; G = tmp[i][j].g; B = tmp[i][j].b;
				R -= dec_col; if (R < 0) R = 0;
				G -= dec_col; if (G < 0) G = 0;
				B -= dec_col; if (B < 0) B = 0;
				tmp[i][j].r = (GLubyte)R;
				tmp[i][j].g = (GLubyte)G;
				tmp[i][j].b = (GLubyte)B;
				//////////////////////////////////////////
				/*
				r = tmp[i][j].r / 255.;
				g = tmp[i][j].g / 255.;
				b = tmp[i][j].b / 255.;
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 0.0;
				v -= 0.001; if (v < 0.0) v = 0.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
				*/
			}
		}
	}
	
}


void Get1DPerpFlowShockMinmax(cimatrix& cmap, Field& gfield, cimatrix& tmp, imatrix& sign, vector& GAU2, double z)
// For each pixel, compute the directional MinMax in the perpendicular direction to the flow
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	//double r, g, b, h, s, v;

	int ss;
	int x1, y1;
	int i, j;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	//int dd;
	//double val;

	int half_w2;
	//half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;

	int image_x, image_y;

	image_x = sign.getRow();
	image_y = sign.getCol();
	
	int c_sign;
	deque<int> vec[3];
	int R, G, B;
	int inc_col, dec_col;
	
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			if (sign[i][j] == 1) c_sign = 1;
			else c_sign = 0;
			vec[0].clear();
			vec[1].clear();
			vec[2].clear();
			vec[0].push_back(cmap[i][j].r);
			vec[1].push_back(cmap[i][j].g);
			vec[2].push_back(cmap[i][j].b);
			
			vn[0] = gfield[i][j].gx;
			vn[1] = gfield[i][j].gy;
			//if (vn[0] == 0.0 && vn[1] == 0.0) {
			//	continue;
			//}
			vn.make_unit();
			d_x = i; d_y = j;
			////////////////////////////////////////
			for (ss = -half_w2; ss <= half_w2; ss++) { // width of Gaussian kernel
				////////////////////////
				x = d_x + vn[0] * ss;
				y = d_y + vn[1] * ss;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					continue;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				//val = image[x1][y1];
				/////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
			}
			sort(vec[0].begin(), vec[0].end()); 
			sort(vec[1].begin(), vec[1].end()); 
			sort(vec[2].begin(), vec[2].end()); 
			if (c_sign) { // non edge region
				//tmp[i][j].r = vec[0][vec[0].size()-1]; 
				//tmp[i][j].g = vec[1][vec[1].size()-1]; 
				//tmp[i][j].b = vec[2][vec[2].size()-1]; 
				tmp[i][j].r = vec[0][vec[0].size()/2]; 
				tmp[i][j].g = vec[1][vec[1].size()/2]; 
				tmp[i][j].b = vec[2][vec[2].size()/2]; 
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
				//tmp[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
				//tmp[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
				//////////////////////////////////////////////////
				inc_col = 7;
				R = tmp[i][j].r; G = tmp[i][j].g; B = tmp[i][j].b;
				R += inc_col; if (R > 255) R = 255;
				G += inc_col; if (G > 255) G = 255;
				B += inc_col; if (B > 255) B = 255;
				tmp[i][j].r = (GLubyte)R;
				tmp[i][j].g = (GLubyte)G;
				tmp[i][j].b = (GLubyte)B;
				//////////////////////////////////////////
				/*
				//////////////////////////////////////////
				r = tmp[i][j].r / 255.;
				g = tmp[i][j].g / 255.;
				b = tmp[i][j].b / 255.;
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 1.0;
				v += 0.05; if (v > 1.0) v = 1.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
				*/
			}
			else {// edge region
				tmp[i][j].r = vec[0][0]; // darkest value
				tmp[i][j].g = vec[1][0];
				tmp[i][j].b = vec[2][0];
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
				//////////////////////////////////////
				dec_col = 2;
				R = tmp[i][j].r; G = tmp[i][j].g; B = tmp[i][j].b;
				R -= dec_col; if (R < 0) R = 0;
				G -= dec_col; if (G < 0) G = 0;
				B -= dec_col; if (B < 0) B = 0;
				tmp[i][j].r = (GLubyte)R;
				tmp[i][j].g = (GLubyte)G;
				tmp[i][j].b = (GLubyte)B;
				//////////////////////////////////////////
				/*
				//////////////////////////////////////////
				r = tmp[i][j].r / 255.;
				g = tmp[i][j].g / 255.;
				b = tmp[i][j].b / 255.;
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 0.0;
				v -= 0.001; if (v < 0.0) v = 0.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
				*/
			}
		}
	}

}


void Get1DFlowColShockAverage(cimatrix& cmap, Field& gfield, cimatrix& tmp, imatrix& sign, vector& GAU3, double z)
// following the flow, compute the average
// Faster version!
{
	vector vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight, w_sum, sum_r, sum_g, sum_b;
	double val_r, val_g, val_b;

	int i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	//double val;
	int c_sign;
	deque<int> vec[3];

	int image_x = sign.getRow();
	int image_y = sign.getCol();

	int half_l;
	//half_w1 = GAU1.getMax()-1;
	//half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;

	int i, j;

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			sum_r = 0.0; sum_g = 0.0; sum_b = 0.0;
			w_sum = 0.0; 
			weight = 0.0; 
			///////////////////////////////
			if (sign[i][j] == 1) c_sign = 1;
			else c_sign = 0;
			/////////////////////////////////
			val_r = cmap[i][j].r;
			val_g = cmap[i][j].g;
			val_b = cmap[i][j].b;
			weight = GAU3[0]; 
			sum_r += val_r * weight;
			sum_g += val_g * weight;
			sum_b += val_b * weight;
			w_sum += weight;
			//w_sum1 += weight1;
			////////////////////////////////////////////////
			// One half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			////////////////////////////
			for (k = 0; k < half_l; k++) {
				vt[0] = -gfield[i_x][i_y].gy;
				vt[1] = gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						val_r = cmap[x1][y1].r;
						val_g = cmap[x1][y1].g;
						val_b = cmap[x1][y1].b;
						weight = GAU3[k]; 
						sum_r += val_r * weight;
						sum_g += val_g * weight;
						sum_b += val_b * weight;
						w_sum += weight;
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						val_r = cmap[x1][y1].r;
						val_g = cmap[x1][y1].g;
						val_b = cmap[x1][y1].b;
						weight = GAU3[k]; 
						sum_r += val_r * weight;
						sum_g += val_g * weight;
						sum_b += val_b * weight;
						w_sum += weight;
					}
				}
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////////////
			// Other half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			for (k = 0; k < half_l; k++) {
				vt[0] = gfield[i_x][i_y].gy;
				vt[1] = -gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						val_r = cmap[x1][y1].r;
						val_g = cmap[x1][y1].g;
						val_b = cmap[x1][y1].b;
						weight = GAU3[k]; 
						sum_r += val_r * weight;
						sum_g += val_g * weight;
						sum_b += val_b * weight;
						w_sum += weight;
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						val_r = cmap[x1][y1].r;
						val_g = cmap[x1][y1].g;
						val_b = cmap[x1][y1].b;
						weight = GAU3[k]; 
						sum_r += val_r * weight;
						sum_g += val_g * weight;
						sum_b += val_b * weight;
						w_sum += weight;
					}
				}
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			sum_r /= w_sum; // normalize
			sum_g /= w_sum; // normalize
			sum_b /= w_sum; // normalize
			//////////////////////////////////////
			if (c_sign) { // non edge region
				tmp[i][j].r = round(sum_r);
				tmp[i][j].g = round(sum_g); 
				tmp[i][j].b = round(sum_b); 
				//tmp[i][j].r = round( (1-z) * sum_r + z * 255 );
				//tmp[i][j].g = round( (1-z) * sum_g + z * 255 ); 
				//tmp[i][j].b = round( (1-z) * sum_b + z * 255 ); 
			}
			else {// edge region
				//tmp[i][j].r = round( (1-z) * sum_r + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * sum_g + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * sum_b + z * 0 ); // darkest value
				tmp[i][j].r = round( (1-z) * sum_r + z * 0 ); // darkest value
				tmp[i][j].g = round( (1-z) * sum_g + z * 0 ); // darkest value
				tmp[i][j].b = round( (1-z) * sum_b + z * 0 ); // darkest value
			}
		}
	}
}

void Get1DFlowColShockAverageHSV(cimatrix& cmap, Field& gfield, cimatrix& tmp, imatrix& sign, vector& GAU3, 
								 double z1, double z2)
// Faster version!
{
	vector vt(2), vt_old(2);
	double x, y, d_x, d_y;
	//double weight1, w_sum1, sum1;
	double r, g, b, h, s, v;
	double ave_r, ave_g, ave_b, min_r, min_g, min_b;
	int count;

	int i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	//double val;
	int c_sign;
	deque<int> vec[3];

	int image_x = sign.getRow();
	int image_y = sign.getCol();

	int half_l;
	//half_w1 = GAU1.getMax()-1;
	//half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;

	int i, j;

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			//sum1 = 0.0;
			//w_sum1 = 0.0;
			//weight1 = 0.0;
			if (sign[i][j] == 1) c_sign = 1;
			else c_sign = 0;
			//vec[0].clear();
			//vec[1].clear();
			//vec[2].clear();
			//vec[0].push_back(cmap[i][j].r);
			//vec[1].push_back(cmap[i][j].g);
			//vec[2].push_back(cmap[i][j].b);
			ave_r = ave_g = ave_b = 0.0;
			min_r = min_g = min_b = 1000.0;
			count = 0;

			ave_r += cmap[i][j].r;
			ave_g += cmap[i][j].g;
			ave_b += cmap[i][j].b;
			count++;
			/////////////////////////////////
			//val = dog[i][j];
			//weight1 = GAU3[0]; 
			//sum1 = val * weight1;
			//w_sum1 += weight1;
			////////////////////////////////////////////////
			// One half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			////////////////////////////
			for (k = 0; k < half_l; k++) {
				vt[0] = -gfield[i_x][i_y].gy;
				vt[1] = gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						//vec[0].push_back(cmap[x1][y1].r);
						//vec[1].push_back(cmap[x1][y1].g);
						//vec[2].push_back(cmap[x1][y1].b);
						ave_r += cmap[x1][y1].r;
						ave_g += cmap[x1][y1].g;
						ave_b += cmap[x1][y1].b;
						count++;
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						if (cmap[x1][y1].r < min_r) min_r = cmap[x1][y1].r;
						if (cmap[x1][y1].g < min_g) min_g = cmap[x1][y1].g;
						if (cmap[x1][y1].b < min_b) min_b = cmap[x1][y1].b;
					}
				}
				//val = dog[x1][y1];
				//////////////////////////////
				//weight1 *= GAU3[k]; 
				////////////////////
				//sum1 += val * weight1;
				//w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////////////
			// Other half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			for (k = 0; k < half_l; k++) {
				vt[0] = gfield[i_x][i_y].gy;
				vt[1] = -gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						ave_r += cmap[x1][y1].r;
						ave_g += cmap[x1][y1].g;
						ave_b += cmap[x1][y1].b;
						count++;
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						if (cmap[x1][y1].r < min_r) min_r = cmap[x1][y1].r;
						if (cmap[x1][y1].g < min_g) min_g = cmap[x1][y1].g;
						if (cmap[x1][y1].b < min_b) min_b = cmap[x1][y1].b;
					}
				}
				//val = dog[x1][y1];
				//////////////////////////////
				//weight1 *= GAU3[k]; 
				////////////////////
				//sum1 += val * weight1;
				//w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			//sort(vec[0].begin(), vec[0].end()); 
			//sort(vec[1].begin(), vec[1].end()); 
			//sort(vec[2].begin(), vec[2].end()); 
			if (c_sign) { // non edge region
				//tmp[i][j].r = vec[0][vec[0].size()/2]; 
				//tmp[i][j].g = vec[1][vec[1].size()/2]; 
				//tmp[i][j].b = vec[2][vec[2].size()/2]; 
				tmp[i][j].r = (GLubyte)(ave_r / count); 
				tmp[i][j].g = (GLubyte)(ave_g / count);  
				tmp[i][j].b = (GLubyte)(ave_b / count);  
				//////////////////////////////////////////
				r = tmp[i][j].r / 255.;
				g = tmp[i][j].g / 255.;
				b = tmp[i][j].b / 255.;
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 1.0;
				//v += 0.05; if (v > 1.0) v = 1.0;
				v += z1; if (v > 1.0) v = 1.0;
				//v += 0.5; if (v > 1.0) v = 1.0;
				//s += 0.5; if (s > 1.0) s = 1.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
			}
			else {// edge region
				tmp[i][j].r = (GLubyte)min_r; // darkest value
				tmp[i][j].g = (GLubyte)min_g;
				tmp[i][j].b = (GLubyte)min_b;
				//////////////////////////////////////////
				r = tmp[i][j].r / 255.;
				g = tmp[i][j].g / 255.;
				b = tmp[i][j].b / 255.;
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 0.0;
				v -= z2; if (v < 0.0) v = 0.0;
				//v -= 0.001; if (v < 0.0) v = 0.0;
				//v -= 0.1; if (v < 0.0) v = 0.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
			}
		}
	}
	
}

void Get1DFlowColShockAverageHSVETF(cimatrix& cmap, ETF& e, cimatrix& tmp, imatrix& sign, vector& GAU3, 
								 double z1, double z2)
// Faster version!
{
	vector vt(2), vt_old(2);
	double x, y, d_x, d_y;
	//double weight1, w_sum1, sum1;
	double r, g, b, h, s, v;
	double ave_r, ave_g, ave_b, min_r, min_g, min_b;
	int count;

	int i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	//double val;
	int c_sign;
	deque<int> vec[3];

	int image_x = sign.getRow();
	int image_y = sign.getCol();

	int half_l;
	//half_w1 = GAU1.getMax()-1;
	//half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;

	int i, j;

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			//sum1 = 0.0;
			//w_sum1 = 0.0;
			//weight1 = 0.0;
			if (sign[i][j] == 1) c_sign = 1;
			else c_sign = 0;
			//vec[0].clear();
			//vec[1].clear();
			//vec[2].clear();
			//vec[0].push_back(cmap[i][j].r);
			//vec[1].push_back(cmap[i][j].g);
			//vec[2].push_back(cmap[i][j].b);
			ave_r = ave_g = ave_b = 0.0;
			min_r = min_g = min_b = 1000.0;
			count = 0;

			ave_r += cmap[i][j].r;
			ave_g += cmap[i][j].g;
			ave_b += cmap[i][j].b;
			count++;
			/////////////////////////////////
			//val = dog[i][j];
			//weight1 = GAU3[0]; 
			//sum1 = val * weight1;
			//w_sum1 += weight1;
			////////////////////////////////////////////////
			// One half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			////////////////////////////
			for (k = 0; k < half_l; k++) {
				vt[0] = e[i_x][i_y].tx;
				vt[1] = e[i_x][i_y].ty;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						//vec[0].push_back(cmap[x1][y1].r);
						//vec[1].push_back(cmap[x1][y1].g);
						//vec[2].push_back(cmap[x1][y1].b);
						ave_r += cmap[x1][y1].r;
						ave_g += cmap[x1][y1].g;
						ave_b += cmap[x1][y1].b;
						count++;
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						if (cmap[x1][y1].r < min_r) min_r = cmap[x1][y1].r;
						if (cmap[x1][y1].g < min_g) min_g = cmap[x1][y1].g;
						if (cmap[x1][y1].b < min_b) min_b = cmap[x1][y1].b;
					}
				}
				//val = dog[x1][y1];
				//////////////////////////////
				//weight1 *= GAU3[k]; 
				////////////////////
				//sum1 += val * weight1;
				//w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////////////
			// Other half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			for (k = 0; k < half_l; k++) {
				vt[0] = -e[i_x][i_y].tx;
				vt[1] = -e[i_x][i_y].ty;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						ave_r += cmap[x1][y1].r;
						ave_g += cmap[x1][y1].g;
						ave_b += cmap[x1][y1].b;
						count++;
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						if (cmap[x1][y1].r < min_r) min_r = cmap[x1][y1].r;
						if (cmap[x1][y1].g < min_g) min_g = cmap[x1][y1].g;
						if (cmap[x1][y1].b < min_b) min_b = cmap[x1][y1].b;
					}
				}
				//val = dog[x1][y1];
				//////////////////////////////
				//weight1 *= GAU3[k]; 
				////////////////////
				//sum1 += val * weight1;
				//w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			//sort(vec[0].begin(), vec[0].end()); 
			//sort(vec[1].begin(), vec[1].end()); 
			//sort(vec[2].begin(), vec[2].end()); 
			if (c_sign) { // non edge region
				//tmp[i][j].r = vec[0][vec[0].size()/2]; 
				//tmp[i][j].g = vec[1][vec[1].size()/2]; 
				//tmp[i][j].b = vec[2][vec[2].size()/2]; 
				tmp[i][j].r = (GLubyte)(ave_r / count); 
				tmp[i][j].g = (GLubyte)(ave_g / count);  
				tmp[i][j].b = (GLubyte)(ave_b / count);  
				//////////////////////////////////////////
				r = tmp[i][j].r / 255.;
				g = tmp[i][j].g / 255.;
				b = tmp[i][j].b / 255.;
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 1.0;
				//v += 0.05; if (v > 1.0) v = 1.0;
				v += z1; if (v > 1.0) v = 1.0;
				//v += 0.5; if (v > 1.0) v = 1.0;
				//s += 0.5; if (s > 1.0) s = 1.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
			}
			else {// edge region
				tmp[i][j].r = (GLubyte)min_r; // darkest value
				tmp[i][j].g = (GLubyte)min_g;
				tmp[i][j].b = (GLubyte)min_b;
				//////////////////////////////////////////
				r = tmp[i][j].r / 255.;
				g = tmp[i][j].g / 255.;
				b = tmp[i][j].b / 255.;
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 0.0;
				v -= z2; if (v < 0.0) v = 0.0;
				//v -= 0.001; if (v < 0.0) v = 0.0;
				//v -= 0.1; if (v < 0.0) v = 0.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
			}
		}
	}
	
}


void Get1DPerpFlowShockAverageHSV(cimatrix& cmap, Field& gfield, cimatrix& tmp, imatrix& sign, 
								  vector& GAU2, double z1, double z2)
// For each pixel, compute the directional MinMax in the perpendicular direction to the flow
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double r, g, b, h, s, v;
	double ave_r, ave_g, ave_b, min_r, min_g, min_b;
	int count;

	int ss;
	int x1, y1;
	int i, j;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	//int dd;
	//double val;

	int half_w2;
	//half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;

	int image_x, image_y;

	image_x = sign.getRow();
	image_y = sign.getCol();
	
	int c_sign;
	deque<int> vec[3];
	
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			if (sign[i][j] == 1) c_sign = 1;
			else c_sign = 0;

			ave_r = ave_g = ave_b = 0.0;
			min_r = min_g = min_b = 1000.0;
			count = 0;

			ave_r += cmap[i][j].r;
			ave_g += cmap[i][j].g;
			ave_b += cmap[i][j].b;
			count++;
			
			vn[0] = gfield[i][j].gx;
			vn[1] = gfield[i][j].gy;
			//if (vn[0] == 0.0 && vn[1] == 0.0) {
			//	continue;
			//}
			vn.make_unit();
			d_x = i; d_y = j;
			////////////////////////////////////////
			for (ss = -half_w2; ss <= half_w2; ss++) { // width of Gaussian kernel
				////////////////////////
				x = d_x + vn[0] * ss;
				y = d_y + vn[1] * ss;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					continue;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				//val = image[x1][y1];
				/////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						ave_r += cmap[x1][y1].r;
						ave_g += cmap[x1][y1].g;
						ave_b += cmap[x1][y1].b;
						count++;
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						if (cmap[x1][y1].r < min_r) min_r = cmap[x1][y1].r;
						if (cmap[x1][y1].g < min_g) min_g = cmap[x1][y1].g;
						if (cmap[x1][y1].b < min_b) min_b = cmap[x1][y1].b;
					}
				}
			}
			sort(vec[0].begin(), vec[0].end()); 
			sort(vec[1].begin(), vec[1].end()); 
			sort(vec[2].begin(), vec[2].end()); 
			if (c_sign) { // non edge region
				r = (ave_r / count)  / 255.; 
				g = (ave_g / count) / 255.; 
				b = (ave_b / count) / 255.; 
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
				//tmp[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
				//tmp[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
				//////////////////////////////////////////
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 1.0;
				v += z1; if (v > 1.0) v = 1.0;
				//v += 0.5; if (v > 1.0) v = 1.0;
				//s += 0.5; if (s > 1.0) s = 1.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
			}
			else {// edge region
				r = min_r / 255.; // darkest value
				g = min_g / 255.;
				b = min_b / 255.;
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
				//////////////////////////////////////
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 0.0;
				v -= z2; if (v < 0.0) v = 0.0;
				//v -= 0.001; if (v < 0.0) v = 0.0;
				//v -= 0.1; if (v < 0.0) v = 0.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
			}
		}
	}
}

void Get1DPerpFlowShockAverageHSVETF(cimatrix& cmap, ETF& e, cimatrix& tmp, imatrix& sign, 
								  vector& GAU2, double z1, double z2)
// For each pixel, compute the directional MinMax in the perpendicular direction to the flow
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double r, g, b, h, s, v;
	double ave_r, ave_g, ave_b, min_r, min_g, min_b;
	int count;

	int ss;
	int x1, y1;
	int i, j;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	//int dd;
	//double val;

	int half_w2;
	//half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;

	int image_x, image_y;

	image_x = sign.getRow();
	image_y = sign.getCol();
	
	int c_sign;
	deque<int> vec[3];
	
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			if (sign[i][j] == 1) c_sign = 1;
			else c_sign = 0;

			ave_r = ave_g = ave_b = 0.0;
			min_r = min_g = min_b = 1000.0;
			count = 0;

			ave_r += cmap[i][j].r;
			ave_g += cmap[i][j].g;
			ave_b += cmap[i][j].b;
			count++;
			
			vn[0] = -e[i][j].ty;
			vn[1] = e[i][j].tx;
			//if (vn[0] == 0.0 && vn[1] == 0.0) {
			//	continue;
			//}
			vn.make_unit();
			d_x = i; d_y = j;
			////////////////////////////////////////
			for (ss = -half_w2; ss <= half_w2; ss++) { // width of Gaussian kernel
				////////////////////////
				x = d_x + vn[0] * ss;
				y = d_y + vn[1] * ss;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					continue;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				//val = image[x1][y1];
				/////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						ave_r += cmap[x1][y1].r;
						ave_g += cmap[x1][y1].g;
						ave_b += cmap[x1][y1].b;
						count++;
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						if (cmap[x1][y1].r < min_r) min_r = cmap[x1][y1].r;
						if (cmap[x1][y1].g < min_g) min_g = cmap[x1][y1].g;
						if (cmap[x1][y1].b < min_b) min_b = cmap[x1][y1].b;
					}
				}
			}
			sort(vec[0].begin(), vec[0].end()); 
			sort(vec[1].begin(), vec[1].end()); 
			sort(vec[2].begin(), vec[2].end()); 
			if (c_sign) { // non edge region
				r = (ave_r / count)  / 255.; 
				g = (ave_g / count) / 255.; 
				b = (ave_b / count) / 255.; 
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
				//tmp[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
				//tmp[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
				//////////////////////////////////////////
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 1.0;
				v += z1; if (v > 1.0) v = 1.0;
				//v += 0.5; if (v > 1.0) v = 1.0;
				//s += 0.5; if (s > 1.0) s = 1.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
			}
			else {// edge region
				r = min_r / 255.; // darkest value
				g = min_g / 255.;
				b = min_b / 255.;
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
				//////////////////////////////////////
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 0.0;
				v -= z2; if (v < 0.0) v = 0.0;
				//v -= 0.001; if (v < 0.0) v = 0.0;
				//v -= 0.1; if (v < 0.0) v = 0.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
			}
		}
	}
}


void Get1DFlowColShockMinMaxHSV(cimatrix& cmap, Field& gfield, cimatrix& tmp, imatrix& sign, vector& GAU3, double z)
// following the flow, compute the MinMax
// Faster version!
{
	vector vt(2), vt_old(2);
	double x, y, d_x, d_y;
	//double weight1, w_sum1, sum1;
	double r, g, b, h, s, v;

	int i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	//double val;
	int c_sign;
	deque<int> vec[3];

	int image_x = sign.getRow();
	int image_y = sign.getCol();

	int half_l;
	//half_w1 = GAU1.getMax()-1;
	//half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;

	int i, j;

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			//sum1 = 0.0;
			//w_sum1 = 0.0;
			//weight1 = 0.0;
			if (sign[i][j] == 1) c_sign = 1;
			else c_sign = 0;
			vec[0].clear();
			vec[1].clear();
			vec[2].clear();
			vec[0].push_back(cmap[i][j].r);
			vec[1].push_back(cmap[i][j].g);
			vec[2].push_back(cmap[i][j].b);
			/////////////////////////////////
			//val = dog[i][j];
			//weight1 = GAU3[0]; 
			//sum1 = val * weight1;
			//w_sum1 += weight1;
			////////////////////////////////////////////////
			// One half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			////////////////////////////
			for (k = 0; k < half_l; k++) {
				vt[0] = -gfield[i_x][i_y].gy;
				vt[1] = gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				//val = dog[x1][y1];
				//////////////////////////////
				//weight1 *= GAU3[k]; 
				////////////////////
				//sum1 += val * weight1;
				//w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////////////
			// Other half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			for (k = 0; k < half_l; k++) {
				vt[0] = gfield[i_x][i_y].gy;
				vt[1] = -gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				//val = dog[x1][y1];
				//////////////////////////////
				//weight1 *= GAU3[k]; 
				////////////////////
				//sum1 += val * weight1;
				//w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			sort(vec[0].begin(), vec[0].end()); 
			sort(vec[1].begin(), vec[1].end()); 
			sort(vec[2].begin(), vec[2].end()); 
			if (c_sign) { // non edge region
				//tmp[i][j].r = vec[0][vec[0].size()-1]; 
				//tmp[i][j].g = vec[1][vec[1].size()-1]; 
				//tmp[i][j].b = vec[2][vec[2].size()-1]; 
				tmp[i][j].r = vec[0][vec[0].size()/2]; 
				tmp[i][j].g = vec[1][vec[1].size()/2]; 
				tmp[i][j].b = vec[2][vec[2].size()/2]; 
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
				//tmp[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
				//tmp[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
				//////////////////////////////////////////
				r = tmp[i][j].r / 255.;
				g = tmp[i][j].g / 255.;
				b = tmp[i][j].b / 255.;
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 1.0;
				v += 0.05; if (v > 1.0) v = 1.0;
				//v += 0.5; if (v > 1.0) v = 1.0;
				//s += 0.5; if (s > 1.0) s = 1.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
			}
			else {// edge region
				tmp[i][j].r = vec[0][0]; // darkest value
				tmp[i][j].g = vec[1][0];
				tmp[i][j].b = vec[2][0];
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
				//////////////////////////////////////////
				r = tmp[i][j].r / 255.;
				g = tmp[i][j].g / 255.;
				b = tmp[i][j].b / 255.;
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 0.0;
				v -= 0.001; if (v < 0.0) v = 0.0;
				//v -= 0.1; if (v < 0.0) v = 0.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
			}
		}
	}
	
}


void Get1DPerpFlowShockMinmaxHSV(cimatrix& cmap, Field& gfield, cimatrix& tmp, imatrix& sign, vector& GAU2, double z)
// For each pixel, compute the directional MinMax in the perpendicular direction to the flow
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double r, g, b, h, s, v;

	int ss;
	int x1, y1;
	int i, j;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	//int dd;
	//double val;

	int half_w2;
	//half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;

	int image_x, image_y;

	image_x = sign.getRow();
	image_y = sign.getCol();
	
	int c_sign;
	deque<int> vec[3];
	
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			if (sign[i][j] == 1) c_sign = 1;
			else c_sign = 0;
			vec[0].clear();
			vec[1].clear();
			vec[2].clear();
			vec[0].push_back(cmap[i][j].r);
			vec[1].push_back(cmap[i][j].g);
			vec[2].push_back(cmap[i][j].b);
			
			vn[0] = gfield[i][j].gx;
			vn[1] = gfield[i][j].gy;
			//if (vn[0] == 0.0 && vn[1] == 0.0) {
			//	continue;
			//}
			vn.make_unit();
			d_x = i; d_y = j;
			////////////////////////////////////////
			for (ss = -half_w2; ss <= half_w2; ss++) { // width of Gaussian kernel
				////////////////////////
				x = d_x + vn[0] * ss;
				y = d_y + vn[1] * ss;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					continue;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				//val = image[x1][y1];
				/////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
			}
			sort(vec[0].begin(), vec[0].end()); 
			sort(vec[1].begin(), vec[1].end()); 
			sort(vec[2].begin(), vec[2].end()); 
			if (c_sign) { // non edge region
				//tmp[i][j].r = vec[0][vec[0].size()-1]; 
				//tmp[i][j].g = vec[1][vec[1].size()-1]; 
				//tmp[i][j].b = vec[2][vec[2].size()-1]; 
				r = vec[0][vec[0].size()/2] / 255.; 
				g = vec[1][vec[1].size()/2] / 255.; 
				b = vec[2][vec[2].size()/2] / 255.; 
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
				//tmp[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
				//tmp[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
				//////////////////////////////////////////
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 1.0;
				v += 0.05; if (v > 1.0) v = 1.0;
				//v += 0.5; if (v > 1.0) v = 1.0;
				//s += 0.5; if (s > 1.0) s = 1.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
			}
			else {// edge region
				r = vec[0][0] / 255.; // darkest value
				g = vec[1][0] / 255.;
				b = vec[2][0] / 255.;
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
				//////////////////////////////////////
				RGB2HSV(r, g, b, h, s, v);
				///////////////////////////////
				//v = (1 - z) * v + z * 0.0;
				v -= 0.001; if (v < 0.0) v = 0.0;
				//v -= 0.1; if (v < 0.0) v = 0.0;
				///////////////////////////////
				HSV2RGB(h, s, v, r, g, b);
				tmp[i][j].r = (GLubyte)(r*255);
				tmp[i][j].g = (GLubyte)(g*255);
				tmp[i][j].b = (GLubyte)(b*255);
				////////////////////////////////////
			}
		}
	}
}


void GetFlowColShockETF(cimatrix& cmap, ETF& e, imatrix& gray, 
			  double sigma, double sigma3, double tau, double z1, double z2) 
// Faster version!
{
	//int	i, j;
	double MAX_DOG = -1.0;
	//int s, t;
	//int x, y;

    int image_x = gray.getRow();
	int image_y = gray.getCol();

	//int half1, half2; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	//double sum1, sum2, w_sum1, w_sum2;

	//////////////////////////////////////////////
	//double factor = 2.0;
	//MakeGaussMatrix(factor, GAU); // factor controls the Gaussian window size for each index
	vector GAU1, GAU2, GAU3;
	MakeGaussVector(sigma, GAU1); // entire inner circle
	MakeGaussVector(sigma*1.6, GAU2); // entire outer circle

	int half_w1, half_w2, half_l;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;
	
	//double factor = 3.0;
	//MakeGaussVector(sigma*factor, GAU3); // flow length
	MakeGaussVector(sigma3, GAU3); // flow length
	half_l = GAU3.getMax()-1;
	
	//int half_l = round( half_w1 * factor );

	TRACE("half_w1 = %d\n", half_w1);
	TRACE("half_w2 = %d\n", half_w2);
	TRACE("half_l = %d\n", half_l);

	imatrix sign(image_x, image_y);
	cimatrix tmp(image_x, image_y);
	matrix dog(image_x, image_y);
	//tmp.zero();

	//StartTimer();

	GetFlowShockDoGETF(gray, e, dog, GAU1, GAU2, tau);
	//TRACE("DoG computation ended\n");
	//GetFlowShockInfluence(gfield, dog, sign, GAU3);
	GetFlowShockInfluenceETF(e, dog, sign, GAU3);
	//TRACE("Influence computation ended\n");

	//Get1DFlowColShockMinMaxHSV(cmap, gfield, tmp, sign, GAU3, 0.05);
	//Get1DPerpFlowShockMinmaxHSV(tmp, gfield, cmap, sign, GAU2, 0.05);
	Get1DFlowColShockAverageHSVETF(cmap, e, tmp, sign, GAU3, z1, z2);
	//Get1DPerpFlowShockAverageHSV(tmp, gfield, cmap, sign, GAU2, 0.05);
	Get1DPerpFlowShockAverageHSVETF(tmp, e, cmap, sign, GAU1, z1, z2);
	
	/*
	outfile.open("_FShock.txt", ios::app);
	outfile << fixed << setprecision(1);
	outfile << "Elapsed Time2 = " << ElapsedTime() << endl;
	outfile.close();
	
	exit(1);
	*/
}

void GetColMinMaxHSV(cimatrix& cmap, imatrix& lap)
// inspect 5 pixel neighborhood
{
	double r, g, b, h, s, v;
	double max_r, max_g, max_b;
	double min_r, min_g, min_b;
	double min_v, max_v;
	int i, j, x, y;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	cimatrix tmp(image_x, image_y);
	tmp.copy(cmap);

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			//////////////////////////
			if (lap[i][j] > 0) {
				r = tmp[i][j].r / 255.0;
				g = tmp[i][j].g / 255.0;
				b = tmp[i][j].b / 255.0;
				RGB2HSV(r, g, b, h, s, v);
				max_v = v;
				max_r = r; max_g = g; max_b = b;
				/////////////////////////////////////////////////////
				x = i + 1;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] > 0) {
					//////////////////////////////////////////////
					r = tmp[x][y].r / 255.0;
					g = tmp[x][y].g / 255.0;
					b = tmp[x][y].b / 255.0;
					RGB2HSV(r, g, b, h, s, v);
					if (v > max_v) { max_v = v; max_r = r; max_g = g; max_b = b; }
				}
                /////////////////////////////////////////////////////
				x = i - 1;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] > 0) {
					r = tmp[x][y].r / 255.0;
					g = tmp[x][y].g / 255.0;
					b = tmp[x][y].b / 255.0;
					RGB2HSV(r, g, b, h, s, v);
					if (v > max_v) { max_v = v; max_r = r; max_g = g; max_b = b; }
				}
				/////////////////////////////////////////////////////
				x = i;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j + 1;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] > 0) {
					r = tmp[x][y].r / 255.0;
					g = tmp[x][y].g / 255.0;
					b = tmp[x][y].b / 255.0;
					RGB2HSV(r, g, b, h, s, v);
					if (v > max_v) { max_v = v; max_r = r; max_g = g; max_b = b; }
				}
				/////////////////////////////////////////////////////
				x = i;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j - 1;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] > 0) {
					r = tmp[x][y].r / 255.0;
					g = tmp[x][y].g / 255.0;
					b = tmp[x][y].b / 255.0;
					RGB2HSV(r, g, b, h, s, v);
					if (v > max_v) { max_v = v; max_r = r; max_g = g; max_b = b; }
				}
				/////////////////////////////////////////////////////
				cmap[i][j].r = round(max_r * 255.0);
				cmap[i][j].g = round(max_g * 255.0);
				cmap[i][j].b = round(max_b * 255.0);
			}
			///////////////////////////
			/////////////////////////////
			else { // lap[i][j] == 0
				r = tmp[i][j].r / 255.0;
				g = tmp[i][j].g / 255.0;
				b = tmp[i][j].b / 255.0;
				RGB2HSV(r, g, b, h, s, v);
				min_v = v;
				min_r = r; min_g = g; min_b = b;
				/////////////////////////////////////////////////////
				x = i + 1;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] == 0) {
					//////////////////////////////////////////////
					r = tmp[x][y].r / 255.0;
					g = tmp[x][y].g / 255.0;
					b = tmp[x][y].b / 255.0;
					RGB2HSV(r, g, b, h, s, v);
					if (v < min_v) { min_v = v; min_r = r; min_g = g; min_b = b; }
				}
                /////////////////////////////////////////////////////
				x = i - 1;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] == 0) {
					r = tmp[x][y].r / 255.0;
					g = tmp[x][y].g / 255.0;
					b = tmp[x][y].b / 255.0;
					RGB2HSV(r, g, b, h, s, v);
					if (v < min_v) { min_v = v; min_r = r; min_g = g; min_b = b; }
				}
				/////////////////////////////////////////////////////
				x = i;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j + 1;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] == 0) {
					r = tmp[x][y].r / 255.0;
					g = tmp[x][y].g / 255.0;
					b = tmp[x][y].b / 255.0;
					RGB2HSV(r, g, b, h, s, v);
					if (v < min_v) { min_v = v; min_r = r; min_g = g; min_b = b; }
				}
				/////////////////////////////////////////////////////
				x = i;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j - 1;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] == 0) {
					r = tmp[x][y].r / 255.0;
					g = tmp[x][y].g / 255.0;
					b = tmp[x][y].b / 255.0;
					RGB2HSV(r, g, b, h, s, v);
					if (v < min_v) { min_v = v; min_r = r; min_g = g; min_b = b; }
				}
				/////////////////////////////////////////////////////
				cmap[i][j].r = round(min_r * 255.0);
				cmap[i][j].g = round(min_g * 255.0);
				cmap[i][j].b = round(min_b * 255.0);
			}

			
		}
	}
	
}

void GetMinMax(imatrix& gray, imatrix& lap)
// gray version
// inspect 5 pixel neighborhood
{
	double min_v, max_v, v;
	int i, j, x, y;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	imatrix tmp(image_x, image_y);
	tmp.copy(gray);

	max_v = -0.1;
	min_v = 1.1;
	////////////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			//////////////////////////
			if (lap[i][j] > 0) {
				v = tmp[i][j] / 255.0;
				max_v = v;
				/////////////////////////////////////////////////////
				x = i + 1;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] > 0) {
					//////////////////////////////////////////////
					v = tmp[x][y] / 255.0;
					if (v > max_v) { max_v = v; }
				}
                /////////////////////////////////////////////////////
				x = i - 1;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] > 0) {
					//////////////////////////////////////////////
					v = tmp[x][y] / 255.0;
					if (v > max_v) { max_v = v; }
				}
				/////////////////////////////////////////////////////
				x = i;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j + 1;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] > 0) {
					//////////////////////////////////////////////
					v = tmp[x][y] / 255.0;
					if (v > max_v) { max_v = v; }
				}
				/////////////////////////////////////////////////////
				x = i;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j - 1;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] > 0) {
					//////////////////////////////////////////////
					v = tmp[x][y] / 255.0;
					if (v > max_v) { max_v = v; }
				}
				/////////////////////////////////////////////////////
				gray[i][j] = round(max_v * 255.0);
			}
			///////////////////////////
			/////////////////////////////
			else { // lap[i][j] == 0
				v = tmp[i][j] / 255.0;
				min_v = v;
				/////////////////////////////////////////////////////
				x = i + 1;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] == 0) {
					//////////////////////////////////////////////
					v = tmp[x][y] / 255.0;
					if (v < min_v) { min_v = v; }
				}
                /////////////////////////////////////////////////////
				x = i - 1;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] == 0) {
					//////////////////////////////////////////////
					v = tmp[x][y] / 255.0;
					if (v < min_v) { min_v = v; }
				}
				/////////////////////////////////////////////////////
				x = i;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j + 1;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] == 0) {
					//////////////////////////////////////////////
					v = tmp[x][y] / 255.0;
					if (v < min_v) { min_v = v; }
				}
				/////////////////////////////////////////////////////
				x = i;	if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
				y = j - 1;		if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
				if (lap[x][y] == 0) {
					//////////////////////////////////////////////
					v = tmp[x][y] / 255.0;
					if (v < min_v) { min_v = v; }
				}
				/////////////////////////////////////////////////////
				gray[i][j] = round(min_v * 255.0);
			}
		}
	}
	
}


void GetColMinMaxHSV2(cimatrix& cmap, imatrix& lap)
// inspect 3x3 pixel neighborhood
{
	double r, g, b, h, s, v;
	double max_r, max_g, max_b;
	double min_r, min_g, min_b;
	double min_v, max_v;
	int i, j, x, y, t, u;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	cimatrix tmp(image_x, image_y);
	tmp.copy(cmap);

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			////////////////////////////////
			max_v = -0.1;
			min_v = 1.1;
			//////////////////////////
			for (u = -1; u <= 1; u++) {
				for (t = -1; t <= 1; t++) {
					/////////////////////////
					x = i + u;
					y = j + t;
					if (x < 0) x = 0; if (x > image_x-1) x = image_x-1;
					if (y < 0) y = 0; if (y > image_y-1) y = image_y-1;
					/////////////////////////////////////////
					if (lap[i][j] > 0 && lap[x][y] == 0) continue;
					else if (lap[i][j] == 0 && lap[x][y] > 0) continue;
					////////////////////////////////////////
					r = tmp[x][y].r / 255.0;
					g = tmp[x][y].g / 255.0;
					b = tmp[x][y].b / 255.0;
					RGB2HSV(r, g, b, h, s, v);
					/////////////////////////////////
					if (lap[x][y] > 0) {
						if (v > max_v) { max_v = v; max_r = r; max_g = g; max_b = b; }
					}
					else {
						if (v < min_v) { min_v = v; min_r = r; min_g = g; min_b = b; }
					}
					
				}
			}
			/////////////////////////////////////////////////////
			if (lap[i][j] > 0) {
				cmap[i][j].r = round(max_r * 255.0);
				cmap[i][j].g = round(max_g * 255.0);
				cmap[i][j].b = round(max_b * 255.0);
			}
			else {
				cmap[i][j].r = round(min_r * 255.0);
				cmap[i][j].g = round(min_g * 255.0);
				cmap[i][j].b = round(min_b * 255.0);
			}
			
		}
	}
}

void GetColMaxMin(CDC& dc, cimatrix& cmap, double sigma, int itr) 
// Faster version!
{
	//int	i, j;
	double MAX_DOG = -1.0;

    int image_x = gray.getRow();
	int image_y = gray.getCol();

	//////////////////////////////////////////////
	//imatrix sign(image_x, image_y);
	//cimatrix tmp(image_x, image_y);
	imatrix lap(image_x, image_y);
	imatrix gray(image_x, image_y);
	
	CopyCol2GrayImage(image_x, image_y, cmap, gray); 

	//GetLaplacian2(gray, lap);
	//GetLaplacian3(gray, sigma, lap);
	//GetLaplacian4(gray, sigma, lap);
	double tau = 0.99;
	double thres = 0.5;
	GetLaplacian5(gray, sigma, lap, tau, thres);

	for (int k = 0; k < itr; k++) {
		/////////////////////////////
		/*
		for (int i = 0; i < image_x; i++) {
			for (int j = 0; j < image_y; j++) {
				if (lap[i][j] >= 0) lap[i][j] = 255;
				else lap[i][j] = 0;
			}
		}
		*/
		//GetFlowShockDoG(gray, gfield, dog, GAU1, GAU2, tau);
		//GetFlowShockInfluence(gfield, dog, sign, GAU3);
		GetColMinMaxHSV(cmap, lap);
		//GetColMinMaxHSV2(cmap, lap);
	}

	DrawGrayImage(memDC, image_x, image_y, lap);
	dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
	//CopyGray2Membuffer(lap, Dbuffer);
	//Get1DFlowColShockAverageHSV(cmap, gfield, tmp, sign, GAU3, z1, z2);
	//Get1DPerpFlowShockAverageHSV(tmp, gfield, cmap, sign, GAU1, z1, z2);
}

void GetColMaxMin2(CDC& dc, cimatrix& cmap, double sigma, double tau, int itr) 
// Faster version!
// use tau
{
	//int	i, j;
	double MAX_DOG = -1.0;

    int image_x = gray.getRow();
	int image_y = gray.getCol();

	//////////////////////////////////////////////
	//imatrix sign(image_x, image_y);
	//cimatrix tmp(image_x, image_y);
	imatrix lap(image_x, image_y);
	imatrix gray(image_x, image_y);
	
	CopyCol2GrayImage(image_x, image_y, cmap, gray); 

	//GetLaplacian2(gray, lap);
	//GetLaplacian3(gray, sigma, lap);
	//GetLaplacian4(gray, sigma, lap);
	//double tau = 0.99;
	double thres = 0.2;
	GetLaplacian5(gray, sigma, lap, tau, thres);

	for (int k = 0; k < itr; k++) {
		/////////////////////////////
		/*
		for (int i = 0; i < image_x; i++) {
			for (int j = 0; j < image_y; j++) {
				if (lap[i][j] >= 0) lap[i][j] = 255;
				else lap[i][j] = 0;
			}
		}
		*/
		//GetFlowShockDoG(gray, gfield, dog, GAU1, GAU2, tau);
		//GetFlowShockInfluence(gfield, dog, sign, GAU3);
		GetColMinMaxHSV(cmap, lap);
		//GetColMinMaxHSV2(cmap, lap);
	}

	DrawGrayImage(memDC, image_x, image_y, lap);
	dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
	//CopyGray2Membuffer(lap, Dbuffer);
	//Get1DFlowColShockAverageHSV(cmap, gfield, tmp, sign, GAU3, z1, z2);
	//Get1DPerpFlowShockAverageHSV(tmp, gfield, cmap, sign, GAU1, z1, z2);
}




void Get1DFlowColShockMinMaxLAB(cimatrix& cmap, Field& gfield, cimatrix& tmp, imatrix& sign, vector& GAU3, double z)
// following the flow, compute the MinMax
// Faster version!
{
	vector vt(2), vt_old(2);
	double x, y, d_x, d_y;
	//double weight1, w_sum1, sum1;
	//double r, g, b, h, s, v;
	int L;

	int i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	//double val;
	int c_sign;
	deque<int> vec[3];

	int image_x = sign.getRow();
	int image_y = sign.getCol();

	int half_l;
	//half_w1 = GAU1.getMax()-1;
	//half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	int flow_DOG_sign = 0; 
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;

	int i, j;

	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			//sum1 = 0.0;
			//w_sum1 = 0.0;
			//weight1 = 0.0;
			if (sign[i][j] == 1) c_sign = 1;
			else c_sign = 0;
			vec[0].clear();
			vec[1].clear();
			vec[2].clear();
			vec[0].push_back(cmap[i][j].r);
			vec[1].push_back(cmap[i][j].g);
			vec[2].push_back(cmap[i][j].b);
			/////////////////////////////////
			//val = dog[i][j];
			//weight1 = GAU3[0]; 
			//sum1 = val * weight1;
			//w_sum1 += weight1;
			////////////////////////////////////////////////
			// One half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			////////////////////////////
			for (k = 0; k < half_l; k++) {
				vt[0] = -gfield[i_x][i_y].gy;
				vt[1] = gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				//val = dog[x1][y1];
				//////////////////////////////
				//weight1 *= GAU3[k]; 
				////////////////////
				//sum1 += val * weight1;
				//w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			////////////////////////////////////////////////
			// Other half
			d_x = (double)i; d_y = (double)j; 
			i_x = i; i_y = j;
			for (k = 0; k < half_l; k++) {
				vt[0] = gfield[i_x][i_y].gy;
				vt[1] = -gfield[i_x][i_y].gx;
				if (vt[0] == 0.0 && vt[1] == 0.0) {
					break;
				}
				vt.make_unit();
				x = d_x;
				y = d_y;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					break;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				//val = dog[x1][y1];
				//////////////////////////////
				//weight1 *= GAU3[k]; 
				////////////////////
				//sum1 += val * weight1;
				//w_sum1 += weight1;
				/////////////////////////////////////////
				d_x += vt[0] * step_size1; // accurate new location x
				d_y += vt[1] * step_size1; // accurate new location y
				/////////////////////////////////////////
				i_x = round(d_x); // integer version of new location x
				i_y = round(d_y); // integer version of new location y
				//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
				if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				/////////////////////////
			}
			//sum1 /= w_sum1; // normalize
			//////////////////////////////////////
			//if (sum1 >= 0) // non edge
			//	sign[i][j] = 1; // non edge
			//else  // edge
			//	sign[i][j] = 0; // edge
			sort(vec[0].begin(), vec[0].end()); 
			sort(vec[1].begin(), vec[1].end()); 
			sort(vec[2].begin(), vec[2].end()); 
			//flow_median = vec[vec.size()/2]; 
			//double z = 0.3;
			//if (vec[0].size() == 0 || vec[1].size() == 0 || vec[2].size() == 0) {
			//	shock_val.r = c_val_r;
			//	shock_val.g = c_val_g;
			//	shock_val.b = c_val_b;
			//	return shock_val; 
			//}
			if (c_sign) { // non edge region
				//tmp[i][j].r = vec[0][vec[0].size()-1]; 
				//tmp[i][j].g = vec[1][vec[1].size()-1]; 
				//tmp[i][j].b = vec[2][vec[2].size()-1]; 
				tmp[i][j].r = vec[0][vec[0].size()/2]; 
				tmp[i][j].g = vec[1][vec[1].size()/2]; 
				tmp[i][j].b = vec[2][vec[2].size()/2]; 
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
				//tmp[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
				//tmp[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
				//tmp[i][j].r += 15; if (tmp[i][j].r > 255) tmp[i][j].r = 255;
				L = (int)tmp[i][j].r;
				L += 15; if (L > 255) L = 255;
				tmp[i][j].r = (GLubyte)L;
			}
			else {// edge region
				tmp[i][j].r = vec[0][0]; // darkest value
				tmp[i][j].g = vec[1][0];
				tmp[i][j].b = vec[2][0];
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
				//////////////////////////////////////////
				//tmp[i][j].r -= 3; if (tmp[i][j].r < 0) tmp[i][j].r = 0;
				L = (int)tmp[i][j].r;
				L -= 3; if (L < 0) L = 0;
				tmp[i][j].r = (GLubyte)L;
				////////////////////////////////////
			}
		}
	}
	
}

void Get1DPerpFlowShockMinmaxLAB(cimatrix& cmap, Field& gfield, cimatrix& tmp, imatrix& sign, vector& GAU2, double z)
// For each pixel, compute the directional MinMax in the perpendicular direction to the flow
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	//double r, g, b, h, s, v;

	int ss;
	int x1, y1;
	int i, j;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	//int dd;
	//double val;
	int L;

	int half_w2;
	//half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;

	int image_x, image_y;

	image_x = sign.getRow();
	image_y = sign.getCol();
	
	int c_sign;
	deque<int> vec[3];
	
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			if (sign[i][j] == 1) c_sign = 1;
			else c_sign = 0;
			vec[0].clear();
			vec[1].clear();
			vec[2].clear();
			vec[0].push_back(cmap[i][j].r);
			vec[1].push_back(cmap[i][j].g);
			vec[2].push_back(cmap[i][j].b);
			
			vn[0] = gfield[i][j].gx;
			vn[1] = gfield[i][j].gy;
			//if (vn[0] == 0.0 && vn[1] == 0.0) {
			//	continue;
			//}
			vn.make_unit();
			d_x = i; d_y = j;
			////////////////////////////////////////
			for (ss = -half_w2; ss <= half_w2; ss++) { // width of Gaussian kernel
				////////////////////////
				x = d_x + vn[0] * ss;
				y = d_y + vn[1] * ss;
				/////////////////////////////////////////////////////
				if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					continue;
				x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
				y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
				//val = image[x1][y1];
				/////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////
				if (c_sign) { // non edge region
					if (sign[x1][y1] == 1) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
				else { // edge region
					if (sign[x1][y1] == 0) {
						vec[0].push_back(cmap[x1][y1].r);
						vec[1].push_back(cmap[x1][y1].g);
						vec[2].push_back(cmap[x1][y1].b);
					}
				}
			}
			sort(vec[0].begin(), vec[0].end()); 
			sort(vec[1].begin(), vec[1].end()); 
			sort(vec[2].begin(), vec[2].end()); 
			if (c_sign) { // non edge region
				//tmp[i][j].r = vec[0][vec[0].size()-1]; 
				//tmp[i][j].g = vec[1][vec[1].size()-1]; 
				//tmp[i][j].b = vec[2][vec[2].size()-1]; 
				tmp[i][j].r = vec[0][vec[0].size()/2]; 
				tmp[i][j].g = vec[1][vec[1].size()/2]; 
				tmp[i][j].b = vec[2][vec[2].size()/2]; 
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
				//tmp[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
				//tmp[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
				//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
				//////////////////////////////////////////
				L = (int)tmp[i][j].r;
				L += 15; if (L > 255) L = 255;
				tmp[i][j].r = (GLubyte)L;
				////////////////////////////////////
			}
			else {// edge region
				tmp[i][j].r = vec[0][0]; // darkest value
				tmp[i][j].g = vec[1][0];
				tmp[i][j].b = vec[2][0];
				//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][vec[0].size()/2] + z * 0 ); // darkest value
				//tmp[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
				//tmp[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
				//tmp[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
				//////////////////////////////////////////
				L = (int)tmp[i][j].r;
				L -= 3; if (L < 0) L = 0;
				tmp[i][j].r = (GLubyte)L;
				////////////////////////////////////
			}
		}
	}

}


void GetPointFlowColShock(int i, int j, cimatrix& cmap, Field& gfield, cimatrix& tmp, 
						  vector& GAU1, vector& GAU2, vector& GAU3, imatrix& sign, double z)
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	//double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s, i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	//double val_r, val_g, val_b;
	GLubyte c_val_r, c_val_g, c_val_b;

	int half_w1, half_w2, half_l;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	//double shock_val = 0.0; // flow-based gradient magnitude

	///////////////////////////
	deque<int> vec[3];
	vec[0].clear();
	vec[1].clear();
	vec[2].clear();

	int c_sign;

	if (sign[i][j] == 1) c_sign = 1;
	else c_sign = 0;

	//iRGB shock_val;
		
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;
	/////////////////////////////////////////////////
	// INSIDE and OUTSIDE FLOW DOG
	////////////////////////////////////
	//count = 0; // number of pixels traversed
	////////////////////////////////////////////////
	// One half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	//////////////////////////
	c_val_r = cmap[i_x][i_y].r;
	c_val_g = cmap[i_x][i_y].g;
	c_val_b = cmap[i_x][i_y].b;
	vec[0].push_back(c_val_r);
	vec[1].push_back(c_val_g);
	vec[2].push_back(c_val_b);
	////////////////////////////
	for (k = 0; k < half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = -gfield[i_x][i_y].gy;
		vt[1] = gfield[i_x][i_y].gx;
		vt.make_unit();
		////////////////////////////////////////
		for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s;
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			//val = image[x1][y1];
			/////////////////////////////////////////////////////////
			if (c_sign) { // non edge region
				if (sign[x1][y1] == 1) {
					vec[0].push_back(cmap[x1][y1].r);
					vec[1].push_back(cmap[x1][y1].g);
					vec[2].push_back(cmap[x1][y1].b);
				}
			}
			else { // edge region
				if (sign[x1][y1] == 0) {
					vec[0].push_back(cmap[x1][y1].r);
					vec[1].push_back(cmap[x1][y1].g);
					vec[2].push_back(cmap[x1][y1].b);
				}
			}
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
	}
	////////////////////////////////////////////////
	// Other half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	for (k = 0; k < half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = gfield[i_x][i_y].gy;
		vt[1] = -gfield[i_x][i_y].gx;
		vt.make_unit();
		////////////////////////////////////
		for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s;
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			//val = image[x1][y1];
			/////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////
			if (c_sign) { // non edge region
				if (sign[x1][y1] == 1) {
					vec[0].push_back(cmap[x1][y1].r);
					vec[1].push_back(cmap[x1][y1].g);
					vec[2].push_back(cmap[x1][y1].b);
				}
			}
			else { // edge region
				if (sign[x1][y1] == 0) {
					vec[0].push_back(cmap[x1][y1].r);
					vec[1].push_back(cmap[x1][y1].g);
					vec[2].push_back(cmap[x1][y1].b);
				}
			}
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
	}
	//sum1 /= w_sum1; // normalize
	//sum2 /= w_sum2; // normalize
	sort(vec[0].begin(), vec[0].end()); 
	sort(vec[1].begin(), vec[1].end()); 
	sort(vec[2].begin(), vec[2].end()); 
	//flow_median = vec[vec.size()/2]; 
	//double z = 0.3;
	//if (vec[0].size() == 0 || vec[1].size() == 0 || vec[2].size() == 0) {
	//	shock_val.r = c_val_r;
	//	shock_val.g = c_val_g;
	//	shock_val.b = c_val_b;
	//	return shock_val; 
	//}
	if (c_sign) { // non edge region
		//tmp[i][j].r = vec[0][vec[0].size()-1]; 
		//tmp[i][j].g = vec[1][vec[1].size()-1]; 
		//tmp[i][j].b = vec[2][vec[2].size()-1]; 
		tmp[i][j].r = vec[0][vec[0].size()/2]; 
		tmp[i][j].g = vec[1][vec[1].size()/2]; 
		tmp[i][j].b = vec[2][vec[2].size()/2]; 
		//tmp[i][j].r = round( (1-z) * vec[0][vec[0].size()/2] + z * 255); 
		//tmp[i][j].g = round( (1-z) * vec[1][vec[1].size()/2] + z * 255); 
		//tmp[i][j].b = round( (1-z) * vec[2][vec[2].size()/2] + z * 255); 
		//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
	}
	else {// edge region
		//image2[i][j] = vec[0]; 
		tmp[i][j].r = round( (1-z) * vec[0][0] + z * 0 ); // darkest value
		tmp[i][j].g = round( (1-z) * vec[1][0] + z * 0 ); // darkest value
		tmp[i][j].b = round( (1-z) * vec[2][0] + z * 0 ); // darkest value
	}
	////////////////////////////////////////////////////
	//return shock_val;
}

double GetPointFlowShock3(int i, int j, imatrix& image, Field& gfield, 
							vector& GAU1, vector& GAU2, vector& GAU3, imatrix& sign)
// following the flow, Gaussian sum the signed values
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s, i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	double val, c_val;
	int dd;

	int half_w1, half_w2, half_l;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	double shock_val = 0.0; // flow-based gradient magnitude

	///////////////////////////
	//deque<int> vec;
	//vec.clear();

	int c_sign;

	if (sign[i][j] == 1) c_sign = 1;
	else c_sign = 0;
		
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;
	/////////////////////////////////////////////////
	// INSIDE and OUTSIDE FLOW DOG
	////////////////////////////////////
	//count = 0; // number of pixels traversed
	sum1 = sum2 = 0.0;
	w_sum1 = w_sum2 = 0.0;
	weight1 = weight2 = 0.0;
	////////////////////////////////////////////////
	// One half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	//////////////////////////
	c_val = image[i_x][i_y];
	weight1 = GAU1[0];
	weight1 *= GAU3[0];
	sum1 += c_val * weight1;
	w_sum1 += weight1;
	////////////////////////////
	for (k = 0; k < half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = -gfield[i_x][i_y].gy;
		vt[1] = gfield[i_x][i_y].gx;
		vt.make_unit();
		////////////////////////////////////////
		for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s;
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			val = image[x1][y1];
			/////////////////////////////////////////////////////////
			dd = ABS(s);
			if (dd > half_w1) weight1 = 0.0;
			else if (c_sign) { // non edge region
                if (sign[x1][y1] == 1) weight1 = GAU1[dd];
				else weight1 = 0.0;
			}
			else { // edge region
				if (sign[x1][y1] == 0) weight1 = GAU1[dd];
				else weight1 = 0.0;
			}
			//////////////////////////////////
			// The following Gaussian smoothing along main axis is essential for good quality!
			weight1 *= GAU3[k]; 
			////////////////////
			sum1 += val * weight1;
			w_sum1 += weight1;
			//d3 = round( ABS(val - c_val) );
			/////////////////////////////////////////////////////
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
	}
	////////////////////////////////////////////////
	// Other half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	for (k = 0; k < half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = gfield[i_x][i_y].gy;
		vt[1] = -gfield[i_x][i_y].gx;
		vt.make_unit();
		////////////////////////////////////
		for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s;
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			val = image[x1][y1];
			/////////////////////////////////////////////////////////
			dd = ABS(s);
			if (dd > half_w1) weight1 = 0.0;
			else if (c_sign) { // non edge region
                if (sign[x1][y1] == 1) weight1 = GAU1[dd];
				else weight1 = 0.0;
			}
			else { // edge region
				if (sign[x1][y1] == 0) weight1 = GAU1[dd];
				else weight1 = 0.0;
			}
			//////////////////////////////////
			// The following Gaussian smoothing along main axis is essential for good quality!
			weight1 *= GAU3[k]; 
			////////////////////
			sum1 += val * weight1;
			w_sum1 += weight1;
			//d3 = round( ABS(val - c_val) );
			/////////////////////////////////////////////////////
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
	}
	sum1 /= w_sum1; // normalize
	//sum2 /= w_sum2; // normalize
	//sort(vec.begin(), vec.end()); 
	//flow_median = vec[vec.size()/2]; 
	double z = 0.3;
	if (c_sign) // non edge region
		//shock_val = vec[vec.size()-1]; 
		//shock_val = sum1;
		shock_val = sum1;
		//image2[i][j] = round( (1-z) * vec[vec.size()-1] + z * 255 ); // brightest value
	else // edge region
		//image2[i][j] = vec[0]; 
		shock_val = round( (1-z) * sum1 + z * 0 ); // darkest value
		//shock_val = sum1;
	////////////////////////////////////////////////////
	return shock_val;
}

void GetFlowShock(imatrix& image, Field& gfield, imatrix& image2, matrix& G_mag, 
			  double sigma, double sigma3, double tau, double z) 
// Flow DOG
{
	int	i, j;
	double MAX_DOG = -1.0;
	//int s, t;
	//int x, y;

    int image_x = image.getRow();
	int image_y = image.getCol();

	//int half1, half2; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	//double sum1, sum2, w_sum1, w_sum2;

	//////////////////////////////////////////////
	//double factor = 2.0;
	//MakeGaussMatrix(factor, GAU); // factor controls the Gaussian window size for each index
	vector GAU1, GAU2, GAU3;
	MakeGaussVector(sigma, GAU1); // entire inner circle
	MakeGaussVector(sigma*1.6, GAU2); // entire outer circle

	int half_w1, half_w2, half_l;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;
	
	//double factor = 3.0;
	//MakeGaussVector(sigma*factor, GAU3); // flow length
	MakeGaussVector(sigma3, GAU3); // flow length
	half_l = GAU3.getMax()-1;
	
	//int half_l = round( half_w1 * factor );

	TRACE("half_w1 = %d\n", half_w1);
	TRACE("half_w2 = %d\n", half_w2);
	TRACE("half_l = %d\n", half_l);

	imatrix sign(image_x, image_y);
	matrix tmp(image_x, image_y);
	//tmp.zero();

	//StartTimer();
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			//tmp[i][j] = GetPointFlowDOG(i, j, image, gfield, GAU1, GAU2, tau);
			//tmp[i][j] = GetPointFlowDOG3(i, j, image, gfield, GAU1, GAU2, GAU3, tau);
			sign[i][j] = GetPointFlowShock(i, j, gray, gfield, GAU1, GAU2, GAU3, tau);
			tmp[i][j] = GetPointFlowShock2(i, j, image, gfield, GAU1, GAU2, GAU3, sign, z);
			//tmp[i][j] = GetPointFlowShock3(i, j, image, gfield, GAU1, GAU2, GAU3, sign);
		}
	}
	//TRACE("Elapsed Time = %f\n", ElapsedTime());
	////////////////////////////////////////////

	// Normalize DOG value (not necessary! - MAX_DOG is always 1)
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			//tmp[i][j] = tmp[i][j] / MAX_DOG; // place it between [0, 1]
			//G_mag[i][j] = 1-tmp[i][j]; // used for nonmaxima suppression
			//TRACE("G_mag[%d][%d] = %.1f\n", i, j, G_mag[i][j]);
			//TRACE("tmp[%d][%d] = %0.3f\n", i, j, tmp[i][j]);
			//TRACE("G_mag[%d][%d] = %0.3f\n", i, j, G_mag[i][j]);
			if (tmp[i][j] < 0)
				TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);
			//image2[i][j] = round(tmp[i][j] * 255.);
			image2[i][j] = round(tmp[i][j]);
			//image2[i][j] = 255 - round(tmp[i][j] * 255.);
		}
	}
}

void GetFlowColShock(cimatrix& cmap, Field& gfield, imatrix& gray, 
			  double sigma, double sigma3, double tau, double z) 
// Flow DOG
{
	int	i, j;
	double MAX_DOG = -1.0;
	//int s, t;
	//int x, y;

    int image_x = gray.getRow();
	int image_y = gray.getCol();

	//int half1, half2; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	//double sum1, sum2, w_sum1, w_sum2;

	//////////////////////////////////////////////
	//double factor = 2.0;
	//MakeGaussMatrix(factor, GAU); // factor controls the Gaussian window size for each index
	vector GAU1, GAU2, GAU3;
	MakeGaussVector(sigma, GAU1); // entire inner circle
	MakeGaussVector(sigma*1.6, GAU2); // entire outer circle

	int half_w1, half_w2, half_l;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;
	
	//double factor = 3.0;
	//MakeGaussVector(sigma*factor, GAU3); // flow length
	MakeGaussVector(sigma3, GAU3); // flow length
	half_l = GAU3.getMax()-1;
	
	//int half_l = round( half_w1 * factor );

	TRACE("half_w1 = %d\n", half_w1);
	TRACE("half_w2 = %d\n", half_w2);
	TRACE("half_l = %d\n", half_l);

	imatrix sign(image_x, image_y);
	cimatrix tmp(image_x, image_y);
	//tmp.zero();

	//StartTimer();
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			//tmp[i][j] = GetPointFlowDOG(i, j, image, gfield, GAU1, GAU2, tau);
			//tmp[i][j] = GetPointFlowDOG3(i, j, image, gfield, GAU1, GAU2, GAU3, tau);
			sign[i][j] = GetPointFlowShock(i, j, gray, gfield, GAU1, GAU2, GAU3, tau);
			//tmp[i][j] = GetPointFlowShock2(i, j, image, gfield, GAU1, GAU2, GAU3, sign, z);
			GetPointFlowColShock(i, j, cmap, gfield, tmp, GAU1, GAU2, GAU3, sign, z);
			//tmp[i][j] = GetPointFlowShock3(i, j, image, gfield, GAU1, GAU2, GAU3, sign);
		}
	}
	//TRACE("Elapsed Time = %f\n", ElapsedTime());
	////////////////////////////////////////////

	// Normalize DOG value (not necessary! - MAX_DOG is always 1)
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//if (tmp[i][j] < 0)
			//	TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);
			//image2[i][j] = round(tmp[i][j] * 255.);
			cmap[i][j].r = tmp[i][j].r;
			cmap[i][j].g = tmp[i][j].g;
			cmap[i][j].b = tmp[i][j].b;
			//image2[i][j] = 255 - round(tmp[i][j] * 255.);
		}
	}
}


void GetContrast(imatrix& image, matrix& G_mag, imatrix& image2, double sigma, double tao) 
// Using a dodge-and-burn in tone mapping paper
{
	int	i, j, k;
	double MAX_DOG = -1.0;
	int s, t;
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half1; // small window and big window
	//int index1, index2; // index for Gaussian matrix
	double sum1, w_sum1;
	double c_val;

	//////////////////////////////////
	//Make2DGaussMask(sigma); // mask for unilateral filtering

	//N = gau_mask.getRow();
	//half = N / 2;

	//ClearMemDC(&dc); // clear the canvas white
	//////////////////////////////////////////////
	//double factor = 2.0;
	//MakeGaussMatrix(factor, GAU); // factor controls the Gaussian window size for each index
	vector GAU1;
	MakeGaussVector(sigma, GAU1);
	//MakeGaussVector(sigma*1.6, GAU2);

	TRACE("half_w = %d\n", GAU1.getMax());

	////////////////////////////////////////////////////////
	matrix tmp(image_x, image_y);
	//tmp.zero();

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			////////////////////////////
			//index = gray2[i][j];
			sum1 = 0;
			w_sum1 = 0.0;
			//half1 = GAU_W[index1]-1;
			half1 = GAU1.getMax()-1;
			//half = 8;
			//TRACE("GAU_W_1[%d] = %d\n", index1, GAU_W[index1]);
			c_val = (double)image[i][j] / 255.;
			for (s = -half1; s <= half1; s++) {
				for (t = -half1; t <= half1; t++) {
					////////////////////////
					x = i+s; y = j+t;
					/////////////////////////////////////////////////////
					if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
						continue;
					/////////////////////////////////////////////////////////
					// circular kernel
					k = (int)dist2(x, y, i, j);
					//k = round( dist2(x, y, i, j) );
					if ( k > half1 ) continue; 
					/////////////////////////////////////////////////////
					sum1 += image[x][y] * GAU1[k];
					w_sum1 += GAU1[k];
					//TRACE("GAU[%d][%d] = %f\n", index, k, GAU[index][k]);
					//TRACE("k = %d\n", k);
				}
			}
			////////////////////////////////////
			// Normalization
			//if (w_sum == 0.0) TRACE("w_sum = 0.0\n");
			sum1 /= w_sum1; 
			sum1 /= 255.; // now [0, 1]
			//TRACE("sum1 %.6f\n", sum1);
			//////////////////////////////
			
			/////////////////////////////////////////
			// POSITIVE INNTER CIRCLE
			//tmp[i][j] = c_val / (1 + pow(sum1, 30)); // dodge-and-burn, doesn't work!
			//c_val += 20.0 * (c_val - sum1); // works very well!!!
			//c_val += tanh( 5.0 * (c_val - sum1) ); // good!
			c_val += tanh( tao * (c_val - sum1) ); // good!
			//c_val += 0.001 / (c_val - sum1); // doesn't work
			//TRACE("c_val = %.6f\n", c_val);
			if (c_val > 1.0) c_val = 1.0;
			else if (c_val < 0.0) c_val = 0.0;
			tmp[i][j] = c_val;
		}
	}
	////////////////////////////////////////////

	// Normalize DOG value (not necessary! - MAX_DOG is always 1)
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			//tmp[i][j] = tmp[i][j] / MAX_DOG; // place it between [0, 1]
			//////////////////////////////////////////////////
			G_mag[i][j] = 1-tmp[i][j]; // [0, 1] used for nonmaxima suppression
			// the higher, the stronger the edgeness!
			///////////////////////////////////////////////////
			//TRACE("G_mag[%d][%d] = %.1f\n", i, j, G_mag[i][j]);
			//TRACE("tmp[%d][%d] = %0.3f\n", i, j, tmp[i][j]);
			//TRACE("G_mag[%d][%d] = %0.3f\n", i, j, G_mag[i][j]);
			if (tmp[i][j] < 0)
				TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);
			image2[i][j] = round(tmp[i][j] * 255.);
			//image2[i][j] = round( tmp[i][j] );
			//image2[i][j] = 255 - round(tmp[i][j] * 255.);
		}
	}
}




double GetPointFlowDOG3(int i, int j, imatrix& image, Field& gfield, 
							vector& GAU1, vector& GAU2, vector& GAU3, double tau)
// following the flow, compute the DOG
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s, i_x, i_y, k;
	int x1, y1, x2, y2;
	double lb, lt, rb, rt;
	double uu, vv, bb, tt, val;
	int dd;

	//int half_l = ker_l / 2 + 1; // half length of the curved kernel
	//double half_w = 1.414; // half width of the curved kernel
	//double half_w = 1.414; // half width of the curved kernel
	//double half_w = 1; // half width of the curved kernel
	//double half_w = ker_w / 2; // half width of the curved kernel
	int half_w1, half_w2, half_l;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	double flow_DOG = 0.0; // flow-based gradient magnitude
	
	double angle, angle_thres1, angle_thres2;
	angle_thres1 = 90;
	angle_thres2 = 180 - angle_thres1;
	//angle_thres2 = 150;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	step_size1 = step_size;
	/////////////////////////////////////////////////
	// INSIDE and OUTSIDE FLOW DOG
	////////////////////////////////////
	//count = 0; // number of pixels traversed
	sum1 = sum2 = 0.0;
	w_sum1 = w_sum2 = 0.0;
	weight1 = weight2 = 0.0;
	////////////////////////////////////////////////
	// One half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	//v[0] = z[0] = -gfield[int_i][int_j].gy;
	//v[1] = z[1] = gfield[int_i][int_j].gx;
	vt_old[0] = -gfield[i_x][i_y].gy;
	vt_old[1] = gfield[i_x][i_y].gx;
	for (k = 0; k < half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = -gfield[i_x][i_y].gy;
		vt[1] = gfield[i_x][i_y].gx;
		vt.make_unit();
		//////////////////////////////////////////////////////
		angle = vt * vt_old;
		if (angle > 1.0) angle = 1.0;
		if (angle < -1.0) angle = -1.0;
		angle = acos(angle) / PI * 180; // angle in [0, 180]
		//TRACE("angle = %.3f\n", angle);
		if ( angle <= angle_thres1  ) { // the angle between two vectors is between [0, 90]
			step_size1 = step_size; 
		}
		else if ( angle > angle_thres2  ) { // the angle between two vectors is between (90, 180]
		//if ( angle > angle_thres2  ) { // the angle between two vectors is between (150, 180]
			step_size1 = -step_size; // reverse the moving direction!
		}
		////////////////////////////////////////
		for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s;
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			/////////////////////////////////////////////////
			// BILINEAR INTERPOLATION!
			x1 = (int)x;	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			x2 = (int)x+1;	if (x2 < 0) x2 = 0; if (x2 > IMAGE_X-1) x2 = IMAGE_X-1;
			y1 = (int)y;	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			y2 = (int)y+1;	if (y2 < 0) y2 = 0; if (y2 > IMAGE_Y-1) y2 = IMAGE_Y-1;
			lb = (double)image[x1][y1];
			rb = (double)image[x2][y1];
			lt = (double)image[x1][y2];
			rt = (double)image[x2][y2];
			uu = x - x1;
			vv = y - y1;
			bb = (1-uu) * lb + uu * rb;
			tt = (1-uu) * lt + uu * rt;
			val = (1-vv) * bb + vv * tt;
			/////////////////////////////////////////////////////////
			dd = ABS(s);
			if (dd > half_w1) weight1 = 0.0;
			else weight1 = GAU1[dd];
			weight1 *= GAU3[k];
			sum1 += val * weight1;
			w_sum1 += weight1;
			//d3 = round( ABS(val - c_val) );
			//weight *= GAU3[d3];
			/////////////////////////////////////////////////////
			weight2 = GAU2[dd];
			weight2 *= GAU3[k];
			sum2 += val * weight2;
			w_sum2 += weight2;
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
		vt_old.copy(vt);
	}
	////////////////////////////////////////////////
	// Other half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	vt_old[0] = gfield[i_x][i_y].gy;
	vt_old[1] = -gfield[i_x][i_y].gx;
	for (k = 0; k < half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = gfield[i_x][i_y].gy;
		vt[1] = -gfield[i_x][i_y].gx;
		vt.make_unit();
		//////////////////////////////////////////////////////
		angle = vt * vt_old;
		if (angle > 1.0) angle = 1.0;
		if (angle < -1.0) angle = -1.0;
		angle = acos(angle) / PI * 180; // angle in [0, 180]
		//TRACE("angle = %.3f\n", angle);
		if ( angle <= angle_thres1  ) { // the angle between two vectors is between [0, 90]
			step_size1 = step_size; 
		}
		else if ( angle > angle_thres2  ) { // the angle between two vectors is between (90, 180]
			step_size1 = -step_size; // reverse the moving direction!
		}
		//if ( angle > angle_thres2  ) { // the angle between two vectors is between (150, 180]
		//	step_size1 = -step_size; // reverse the moving direction!
		//}
		////////////////////////////////////
		for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s;
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			/////////////////////////////////////////////////
			// BILINEAR INTERPOLATION!
			x1 = (int)x;	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			x2 = (int)x+1;	if (x2 < 0) x2 = 0; if (x2 > IMAGE_X-1) x2 = IMAGE_X-1;
			y1 = (int)y;	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			y2 = (int)y+1;	if (y2 < 0) y2 = 0; if (y2 > IMAGE_Y-1) y2 = IMAGE_Y-1;
			lb = (double)image[x1][y1];
			rb = (double)image[x2][y1];
			lt = (double)image[x1][y2];
			rt = (double)image[x2][y2];
			uu = x - x1;
			vv = y - y1;
			bb = (1-uu) * lb + uu * rb;
			tt = (1-uu) * lt + uu * rt;
			val = (1-vv) * bb + vv * tt;
			/////////////////////////////////////////////////////////
			dd = ABS(s);
			if (dd > half_w1) weight1 = 0.0;
			else weight1 = GAU1[dd];
			sum1 += val * weight1;
			w_sum1 += weight1;
			//d3 = round( ABS(val - c_val) );
			//weight *= GAU3[d3];
			/////////////////////////////////////////////////////
			weight2 = GAU2[dd];
			sum2 += val * weight2;
			w_sum2 += weight2;
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
		vt_old.copy(vt);
	}
	sum1 /= w_sum1; // normalize
	sum2 /= w_sum2; // normalize
	////////////////////////////////////////////////////
	if (sum1 == 0.0 && sum2 == 0.0)
		sum1 = 1.0; // make it a non-edge
	//////////////////////////////////////
	if (sum1 - tau * sum2 > 0) // non edge
		flow_DOG = 1.0; 
	else // edge!
		flow_DOG = 1.0 + tanh(sum1 - tau * sum2);

	return flow_DOG;
	
}

double GetPointFlowDOG4(int i, int j, imatrix& image, Field& gfield, 
							vector& GAU1, vector& GAU2, vector& GAU3, double tau)
// following the flow, compute the DOG
// but do not adjust GVF directions! Just follow it!
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s, i_x, i_y, k;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt;
	int dd;
	double val, c_val;

	int half_w1, half_w2, half_l;
	half_w1 = GAU1.getMax()-1;
	half_w2 = GAU2.getMax()-1;
	half_l = GAU3.getMax()-1;
	
	double flow_DOG = 0.0; // flow-based gradient magnitude
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	//step_size = 2.0;
	step_size1 = step_size;
	/////////////////////////////////////////////////
	// INSIDE and OUTSIDE FLOW DOG
	////////////////////////////////////
	//count = 0; // number of pixels traversed
	sum1 = sum2 = 0.0;
	w_sum1 = w_sum2 = 0.0;
	weight1 = weight2 = 0.0;
	////////////////////////////////////////////////
	// One half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	//////////////////////////
	c_val = image[i_x][i_y];
	weight1 = GAU1[0];
	weight1 *= GAU3[0];
	sum1 += c_val * weight1;
	w_sum1 += weight1;
	weight2 = GAU2[0];
	weight2 *= GAU3[0];
	sum2 += c_val * weight2;
	w_sum2 += weight2;
	////////////////////////////
	for (k = 0; k < half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = -gfield[i_x][i_y].gy;
		vt[1] = gfield[i_x][i_y].gx;
		vt.make_unit();
		////////////////////////////////////////
		for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s;
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			val = image[x1][y1];
			/////////////////////////////////////////////////////////
			dd = ABS(s);
			if (dd > half_w1) weight1 = 0.0;
			else weight1 = GAU1[dd];
			//////////////////////////////////
			// The following Gaussian smoothing along main axis is essential for good quality!
			weight1 *= GAU3[k]; 
			////////////////////
			sum1 += val * weight1;
			w_sum1 += weight1;
			//d3 = round( ABS(val - c_val) );
			/////////////////////////////////////////////////////
			weight2 = GAU2[dd];
			weight2 *= GAU3[k];
			sum2 += val * weight2;
			w_sum2 += weight2;
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
	}
	////////////////////////////////////////////////
	// Other half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	for (k = 0; k < half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = gfield[i_x][i_y].gy;
		vt[1] = -gfield[i_x][i_y].gx;
		vt.make_unit();
		////////////////////////////////////
		for (s = -half_w2; s <= half_w2; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s;
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			val = image[x1][y1];
			/////////////////////////////////////////////////////////
			dd = ABS(s);
			if (dd > half_w1) weight1 = 0.0;
			else weight1 = GAU1[dd];
			weight1 *= GAU3[k];
			sum1 += val * weight1;
			w_sum1 += weight1;
			//d3 = round( ABS(val - c_val) );
			/////////////////////////////////////////////////////
			weight2 = GAU2[dd];
			weight2 *= GAU3[k];
			sum2 += val * weight2;
			w_sum2 += weight2;
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
	}
	sum1 /= w_sum1; // normalize
	sum2 /= w_sum2; // normalize
	////////////////////////////////////////////////////
	if (sum1 == 0.0 && sum2 == 0.0)
		sum1 = 1.0; // make it a non-edge
	//////////////////////////////////////
	if (sum1 - tau * sum2 > 0) // non edge
		flow_DOG = 1.0; 
	else // edge!
		flow_DOG = 1.0 + tanh(sum1 - tau * sum2);

	return flow_DOG;
	
}

inline void EndTimer(char *filename)
{
	ofstream outfile;

	outfile.open(filename, ios::app);
	outfile << fixed << setprecision(1);
	outfile << "Elapsed Time = " << ElapsedTime() << endl;
	outfile.close();
}





int GetPointFlowMedian3(int i, int j, imatrix& image, Field& gfield, int half_w, int half_l)
// following the perpendicular flow
// give weight to pixels (min max filter?)
{
	vector vn(2), vt(2), vn_old(2);
	double x, y, d_x, d_y;
	//double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s, i_x, i_y, k;
	
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt, val;
	//int dd;

	int x1, y1, c_val, val;

	deque<int> vec;
	
	int flow_median; // return value
	
	double angle, angle_thres1, angle_thres2;
	angle_thres1 = 90;
	angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	step_size1 = step_size;
	/////////////////////////////////////////////////
	// FLOW MEDIAN
	////////////////////////////////////
	c_val = image[i][j];
	////////////////////////////////////////////////
	// One half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	vn_old[0] = gfield[i_x][i_y].gx;
	vn_old[1] = gfield[i_x][i_y].gy;
	for (k = 0; k <= half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		//////////////////////////////////////////////////////
		angle = vn * vn_old;
		if (angle > 1.0) angle = 1.0;
		if (angle < -1.0) angle = -1.0;
		angle = acos(angle) / PI * 180; // angle in [0, 180]
		//TRACE("angle = %.3f\n", angle);
		if ( angle <= angle_thres1  ) { // the angle between two vectors is between [0, 90]
			step_size1 = step_size; 
		}
		else if ( angle > angle_thres2  ) { // the angle between two vectors is between (90, 180]
			step_size1 = -step_size; // reverse the moving direction!
		}
		////////////////////////////////////////
		for (s = -half_w; s <= half_w; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vt[0] * s; // arm
			y = d_y + vt[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			/////////////////////////////////////////////////
			// BILINEAR INTERPOLATION!
			/*
			x1 = (int)x;	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			x2 = (int)x+1;	if (x2 < 0) x2 = 0; if (x2 > IMAGE_X-1) x2 = IMAGE_X-1;
			y1 = (int)y;	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			y2 = (int)y+1;	if (y2 < 0) y2 = 0; if (y2 > IMAGE_Y-1) y2 = IMAGE_Y-1;
			lb = (double)image[x1][y1];
			rb = (double)image[x2][y1];
			lt = (double)image[x1][y2];
			rt = (double)image[x2][y2];
			uu = x - x1;
			vv = y - y1;
			bb = (1-uu) * lb + uu * rb;
			tt = (1-uu) * lt + uu * rt;
			val = (1-vv) * bb + vv * tt;
			*/
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			val = image[x1][y1];
			/////////////////////////////////////////////////////////
			vec.push_back(val);
			//if (val < c_val) {
			//	vec.push_back(val);
			//}
			/////////////////////////////////////////////////////
		}
		/////////////////////////////////////////
		d_x += vn[0] * step_size1; // accurate new location x
		d_y += vn[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
		vn_old.copy(vn);
	}
	////////////////////////////////////////////////
	// Other half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	vn_old[0] = gfield[i_x][i_y].gx;
	vn_old[1] = gfield[i_x][i_y].gy;
	for (k = 0; k <= half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		//////////////////////////////////////////////////////
		angle = vn * vn_old;
		if (angle > 1.0) angle = 1.0;
		if (angle < -1.0) angle = -1.0;
		angle = acos(angle) / PI * 180; // angle in [0, 180]
		//TRACE("angle = %.3f\n", angle);
		if ( angle <= angle_thres1  ) { // the angle between two vectors is between [0, 90]
			step_size1 = step_size; 
		}
		else if ( angle > angle_thres2  ) { // the angle between two vectors is between (90, 180]
			step_size1 = -step_size; // reverse the moving direction!
		}
		////////////////////////////////////////
		for (s = -half_w; s <= half_w; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vt[0] * s; // arm
			y = d_y + vt[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			/////////////////////////////////////////////////
			// BILINEAR INTERPOLATION!
			/*
			x1 = (int)x;	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			x2 = (int)x+1;	if (x2 < 0) x2 = 0; if (x2 > IMAGE_X-1) x2 = IMAGE_X-1;
			y1 = (int)y;	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			y2 = (int)y+1;	if (y2 < 0) y2 = 0; if (y2 > IMAGE_Y-1) y2 = IMAGE_Y-1;
			lb = (double)image[x1][y1];
			rb = (double)image[x2][y1];
			lt = (double)image[x1][y2];
			rt = (double)image[x2][y2];
			uu = x - x1;
			vv = y - y1;
			bb = (1-uu) * lb + uu * rb;
			tt = (1-uu) * lt + uu * rt;
			val = (1-vv) * bb + vv * tt;
			*/
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			val = image[x1][y1];
			/////////////////////////////////////////////////////////
			vec.push_back(val);
			//if (val < c_val) {
			//	vec.push_back(val);
			//}
			/////////////////////////////////////////////////////
		}
		/////////////////////////////////////////
		d_x += vn[0] * step_size1; // accurate new location x
		d_y += vn[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
		vn_old.copy(vn);
	}
	////////////////////////////////////////////
	sort(vec.begin(), vec.end()); 
	//flow_median = vec[vec.size()/2]; 
	flow_median = vec[vec.size()-1]; 
	//for (k = count/2; k > 0; k--) {
	//	if (vec[k] != 255) break;
	//}
	//flow_median = vec[k]; 

	return flow_median;
}

int GetPointFlowMedian(int i, int j, imatrix& image, Field& gfield, int half_w, int half_l)
// following the flow, compute the DOG
// give weight to pixels
{
	vector vn(2), vt(2), vt_old(2);
	double x, y, d_x, d_y;
	//double weight1, weight2, w_sum1, sum1, sum2, w_sum2;

	int s, i_x, i_y, k;
	
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double uu, vv, bb, tt, val;
	//int dd;

	int x1, y1, c_val, val;

	deque<int> vec;
	
	int flow_median; // return value
	
	//double angle, angle_thres1, angle_thres2;
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;

	double step_size, step_size1; 
	//double t;

	step_size = 1.0;
	step_size1 = step_size;
	/////////////////////////////////////////////////
	// FLOW MEDIAN
	////////////////////////////////////
	c_val = image[i][j];
	////////////////////////////////////////////////
	// One half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	vt_old[0] = -gfield[i_x][i_y].gy;
	vt_old[1] = gfield[i_x][i_y].gx;
	vt_old.make_unit();
	for (k = 0; k <= half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = -gfield[i_x][i_y].gy;
		vt[1] = gfield[i_x][i_y].gx;
		vt.make_unit();
		//////////////////////////////////////////////////////
		/*
		angle = vt * vt_old;
		if (angle > 1.0) angle = 1.0;
		if (angle < -1.0) angle = -1.0;
		angle = acos(angle) / PI * 180; // angle in [0, 180]
		//TRACE("angle = %.3f\n", angle);
		if ( angle <= angle_thres1  ) { // the angle between two vectors is between [0, 90]
			step_size1 = step_size; 
		}
		else if ( angle > angle_thres2  ) { // the angle between two vectors is between (90, 180]
			step_size1 = -step_size; // reverse the moving direction!
		}
		//else break;
		*/
		////////////////////////////////////////
		for (s = -half_w; s <= half_w; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s; // arm
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			val = image[x1][y1];
			/////////////////////////////////////////////////////////
			vec.push_back(val);
			/////////////////////////////////////////////////////
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
		vt_old.copy(vt);
	}
	////////////////////////////////////////////////
	// Other half
	d_x = (double)i; d_y = (double)j; 
	i_x = i; i_y = j;
	vt_old[0] = gfield[i_x][i_y].gy;
	vt_old[1] = -gfield[i_x][i_y].gx;
	vt_old.make_unit();
	for (k = 0; k <= half_l; k++) {
		vn[0] = gfield[i_x][i_y].gx;
		vn[1] = gfield[i_x][i_y].gy;
		if (vn[0] == 0.0 && vn[1] == 0.0) break;
		vn.make_unit();
		vt[0] = gfield[i_x][i_y].gy;
		vt[1] = -gfield[i_x][i_y].gx;
		vt.make_unit();
		//////////////////////////////////////////////////////
		/*
		angle = vt * vt_old;
		if (angle > 1.0) angle = 1.0;
		if (angle < -1.0) angle = -1.0;
		angle = acos(angle) / PI * 180; // angle in [0, 180]
		//TRACE("angle = %.3f\n", angle);
		if ( angle <= angle_thres1  ) { // the angle between two vectors is between [0, 90]
			step_size1 = step_size; 
		}
		else if ( angle > angle_thres2  ) { // the angle between two vectors is between (90, 180]
			step_size1 = -step_size; // reverse the moving direction!
		}
		//else break;
		*/
		////////////////////////////////////////
		for (s = -half_w; s <= half_w; s++) { // width of Gaussian kernel
			////////////////////////
			x = d_x + vn[0] * s; // arm
			y = d_y + vn[1] * s;
			/////////////////////////////////////////////////////
			if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
				continue;
			x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
			y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
			val = image[x1][y1];
			/////////////////////////////////////////////////////////
			vec.push_back(val);
			/////////////////////////////////////////////////////
		}
		/////////////////////////////////////////
		d_x += vt[0] * step_size1; // accurate new location x
		d_y += vt[1] * step_size1; // accurate new location y
		/////////////////////////////////////////
		i_x = round(d_x); // integer version of new location x
		i_y = round(d_y); // integer version of new location y
		//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
		/////////////////////////
		vt_old.copy(vt);
	}
	////////////////////////////////////////////
	sort(vec.begin(), vec.end()); 
	flow_median = vec[vec.size()/2]; 
	//for (k = count/2; k > 0; k--) {
	//	if (vec[k] != 255) break;
	//}
	//flow_median = vec[k]; 

	return flow_median;
}





void Thresholding(int image_x, int image_y, imatrix& image, double thres) 
{
	int	i, j;
	double val;

	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			val = image[i][j] / 255.0; // [0, 1]
			if (val < thres)
				image[i][j] = 0; // black
			else image[i][j] = 255; // white
		}
	}
}


void EdgeThinning(CDC& dc, imatrix& image) 
// Perform Morphological Edge Thinning!
{
	int	i, j, k;
	int s, t;
	int x, y;
	int image_x, image_y;

	image_x = image.getRow();
	image_y = image.getCol();

	//imatrix tmp(image_x, image_y);
	//tmp.copy(image);

	imatrix p[8];
	for (k = 0; k < 8; k++) { p[k].init(3, 3); }
	/*
	///////////////////////////////////////////
	// weaker constraint for diagonal edges (remove less)
	p[0][0][0] = 1;	p[0][0][1] = 1;	p[0][0][2] = 2; 
	p[0][1][0] = 1;	p[0][1][1] = 2;	p[0][1][2] = 0; 
	p[0][2][0] = 2;	p[0][2][1] = 0;	p[0][2][2] = 0; 

	p[1][0][0] = 2;	p[1][0][1] = 0;	p[1][0][2] = 0; 
	p[1][1][0] = 1;	p[1][1][1] = 2;	p[1][1][2] = 0; 
	p[1][2][0] = 1;	p[1][2][1] = 1;	p[1][2][2] = 2; 
	
	p[2][0][0] = 0;	p[2][0][1] = 0;	p[2][0][2] = 2; 
	p[2][1][0] = 0;	p[2][1][1] = 2;	p[2][1][2] = 1; 
	p[2][2][0] = 2;	p[2][2][1] = 1;	p[2][2][2] = 1;

	p[3][0][0] = 2;	p[3][0][1] = 1;	p[3][0][2] = 1; 
	p[3][1][0] = 0;	p[3][1][1] = 2;	p[3][1][2] = 1; 
	p[3][2][0] = 0;	p[3][2][1] = 0;	p[3][2][2] = 2; 
	*/
	////////////////////////////////////////////
	// stronger constraint for diagonal edges (remove more)
	p[0][0][0] = 2;	p[0][0][1] = 1;	p[0][0][2] = 2; 
	p[0][1][0] = 1;	p[0][1][1] = 2;	p[0][1][2] = 0; 
	p[0][2][0] = 2;	p[0][2][1] = 0;	p[0][2][2] = 0; 

	p[1][0][0] = 2;	p[1][0][1] = 0;	p[1][0][2] = 0; 
	p[1][1][0] = 1;	p[1][1][1] = 2;	p[1][1][2] = 0; 
	p[1][2][0] = 2;	p[1][2][1] = 1;	p[1][2][2] = 2; 
	
	p[2][0][0] = 0;	p[2][0][1] = 0;	p[2][0][2] = 2; 
	p[2][1][0] = 0;	p[2][1][1] = 2;	p[2][1][2] = 1; 
	p[2][2][0] = 2;	p[2][2][1] = 1;	p[2][2][2] = 2;

	p[3][0][0] = 2;	p[3][0][1] = 1;	p[3][0][2] = 2; 
	p[3][1][0] = 0;	p[3][1][1] = 2;	p[3][1][2] = 1; 
	p[3][2][0] = 0;	p[3][2][1] = 0;	p[3][2][2] = 2; 

	p[4][0][0] = 1;	p[4][0][1] = 2;	p[4][0][2] = 0; 
	p[4][1][0] = 1;	p[4][1][1] = 2;	p[4][1][2] = 0; 
	p[4][2][0] = 1;	p[4][2][1] = 2;	p[4][2][2] = 0; 

	p[5][0][0] = 0;	p[5][0][1] = 0;	p[5][0][2] = 0; 
	p[5][1][0] = 2;	p[5][1][1] = 2;	p[5][1][2] = 2; 
	p[5][2][0] = 1;	p[5][2][1] = 1;	p[5][2][2] = 1; 

	p[6][0][0] = 0;	p[6][0][1] = 2;	p[6][0][2] = 1; 
	p[6][1][0] = 0;	p[6][1][1] = 2;	p[6][1][2] = 1; 
	p[6][2][0] = 0;	p[6][2][1] = 2;	p[6][2][2] = 1; 

	p[7][0][0] = 1;	p[7][0][1] = 1;	p[7][0][2] = 1; 
	p[7][1][0] = 2;	p[7][1][1] = 2;	p[7][1][2] = 2; 
	p[7][2][0] = 0;	p[7][2][1] = 0;	p[7][2][2] = 0; 

	bool removable =false;
	bool match_fail = false;
	bool remain = true;
	int round = 1;
	
	while (remain) { // there's still a possibility of remaining removable pixels
		remain = false; // assume that there's no more removable pixels
		for (j = 1; j < image_y-1; j++) {
			for (i = 1; i < image_x-1; i++) {
				//TRACE("image[%d][%d] == %d\n", i, j, image[i][j]);
				if (image[i][j] == 255) continue; // already white pixel
				////////////////////////////
				if (image[i][j] == 0) { // black pixel. let's check if it's removable
					removable = false; // assume this pixel is unremovable
					for (k = 0; k < 8; k++) {
						match_fail = false; // assume this pattern matches
						for (s = 0; s <= 2; s++) {
							for (t = 0; t <= 2; t++) {
								x = i+s-1; y = j+t-1;
								if (p[k][s][t] == 1) { // edge pixel
									if (image[x][y] == 0) continue; // match
									else match_fail = true; // match failed
								}
								else if (p[k][s][t] == 0) { // non-edge pixel
									if (image[x][y] == 255) continue; // match
									else match_fail = true; // match failed
								}
								if (match_fail) break; // get out of for t
							}
							if (match_fail) break; // get out of for s
						}
						if (match_fail == false) { // match success! let's get out!
							removable = true;
							break; // get out of (for k)
						}
					} // for k 
					if (removable) {
						//tmp[i][j] = 255; // incorrect: the connectivity is not preserved!
						// You have to change it and then evaluate the next pixel right away!
						image[i][j] = 255; // make it white 
						//TRACE("Round %d: [%d][%d] removable!\n", round, x, y);
						remain = true; // we have to run at least one more round!
					}
				} // if (tmp[i][j] == 0)
			} // for i
		} // for j
		//tmp.copy(image); // update tmp
		round++;
		//if (round > 100) break;
		//DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, image);
		//dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	} // while (remain)
	////////////////////////////////////////////
}

double GetPointFlowDOG(int x, int y, matrix& G_mag, Field& gfield, 
							int ker_l, double step_size)
// following the flow, compute the average DOG at (x, y) over the flow
// Here we simply add the DOG[i][j] cumulatively (for repetitive application!)
{
	double d_i, d_j, tx, ty, nx, ny;
	int k, int_i, int_j;
	//int half2;
	//double t;
	int count;

	int half_l = ker_l / 2 + 1; // half length of the curved kernel

	//double half_w = 1.414; // half width of the curved kernel
	//double half_w = 1.414; // half width of the curved kernel
	//double half_w = 1; // half width of the curved kernel
	//double half_w = ker_w / 2; // half width of the curved kernel

	double flow_dog = 0.0; // flow-based gradient magnitude

	////////////////////////////////////
	//t = 1.0;
	//t = 1.414; // step_size
	//t = 0.5; // step_size
	count = 0; // number of pixels traversed
	////////////////////////////////////////////////
	// One half
	d_i = (double)x; d_j = (double)y; 
	tx = -gfield[x][y].gy;
	ty = gfield[x][y].gx;
	int_i = x;	int_j = y;
	for (k = 0; k < half_l; k++) {
		if (tx == 0.0 && ty == 0.0) break;
		///////////////////////////////////////
		nx = tx / sqrt(tx*tx+ty*ty);  // x component of the unit direction vector 
		ny = ty / sqrt(tx*tx+ty*ty);  // y component of the unit direction vector
		////////////////////////////////
		/////////////////////////////////////////
		//flow_mag += ABS( image[right_i][right_j] - image[left_i][left_j] );
		//flow_mag += gfield[int_i][int_j].mag;
		flow_dog += G_mag[int_i][int_j];
		//TRACE("G_mag[%d][%d] = %.1f\n", int_i, int_j, G_mag[int_i][int_j]);
		count++;
		//////////////////////////////////
		//dc.SetPixelV(int_i, IMAGE_Y-1-int_j, RGB(r, g, b));
		//////////////////////////////
		d_i += nx * step_size; // accurate new location x
		d_j += ny * step_size; // accurate new location y
        ///////////////////////////
		//if (round(d_i) == int_i && round(d_j) == int_j) // no change
		//	continue; // push some more
		/////////////////////////////////////////
		int_i = round(d_i); // integer version of new location x
		int_j = round(d_j); // integer version of new location y
		if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		
		///////////////////////////////////
		//d_i = (double)int_i;
		//d_j = (double)int_j;
		//////////////////////////////////
		tx = -gfield[int_i][int_j].gy;
		ty = gfield[int_i][int_j].gx;
	}
	///////////////////////////////////////
	// other half
	// Note that center pixel gradient magnitude is computed twice (consistent with Sobel grad)
	///*
	d_i = (double)x; d_j = (double)y; 
	int_i = x;	int_j = y;
	tx = -gfield[x][y].gy;
	ty = gfield[x][y].gx;
	/////////////////////////////////////
	for (k = 0; k < half_l; k++) {
		if (tx == 0.0 && ty == 0.0) break;
		///////////////////////////////////////
		nx = tx / sqrt(tx*tx+ty*ty);  // x component of the unit direction vector 
		ny = ty / sqrt(tx*tx+ty*ty);  // y component of the unit direction vector
		////////////////////////////////
		/////////////////////////////////////////
		//flow_mag += gfield[int_i][int_j].mag;
		flow_dog += G_mag[int_i][int_j];
		count++;
		//////////////////////////////////
		//dc.SetPixelV(int_i, IMAGE_Y-1-int_j, RGB(r, g, b));
		//////////////////////////////
		d_i -= nx * step_size; // accurate new location x
		d_j -= ny * step_size; // accurate new location y
        ///////////////////////////
		//if (round(d_i) == int_i && round(d_j) == int_j) // no change
		//	continue; // push some more
		/////////////////////////////////////////
		int_i = round(d_i); // integer version of new location x
		int_j = round(d_j); // integer version of new location y
		if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		
		//////////////////////////////////
		tx = -gfield[int_i][int_j].gy;
		ty = gfield[int_i][int_j].gx;
	}
	//*/

	if (count > 0)
		flow_dog /= (double)count;

	return flow_dog;
	
}

void GetFlowDOG(int image_x, int image_y, matrix& G_mag, imatrix& dog, Field& gfield, 
				int ker_l, double step_size) 
// Flow DOG: simply get the average DOG value over the flow
{
	int i, j;
	double MAX_CUM_DOG = -1.;
	//double MAX_VAL = 1020.; // 255 + 2 * 255 + 255 (Maximum possible Sobel value)

	matrix tmp(image_x, image_y);

	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//tmp[i][j].mag = GetPointFlowGradient(dc, i, j, image, gfield, ker_l, ker_w, step_size);
			//TRACE("G_mag[%d][%d] = %.1f\n", i, j, G_mag[i][j]);
			tmp[i][j] = GetPointFlowDOG(i, j, G_mag, gfield, ker_l, step_size);
			if (tmp[i][j] > MAX_CUM_DOG) MAX_CUM_DOG = tmp[i][j];
			//if (tmp[i][j] > MAX_DOG) {
			//	MAX_GRADIENT = tmp[i][j].mag;
			//}
		}
	}
	
	TRACE("MAX_CUM_DOG = %0.1f\n", MAX_CUM_DOG);

	// Normalize DOG
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			tmp[i][j] = tmp[i][j] / MAX_CUM_DOG; // place it between [0, 1]
			//G_mag[i][j] = tmp[i][j]; // used for nonmaxima suppresion
			if (tmp[i][j] < 0)
				TRACE("tmp[%d][%d] = %0.2f\n", i, j, tmp[i][j]);

			dog[i][j] = 255 - round(tmp[i][j] * 255.);
			G_mag[i][j] = tmp[i][j];
		}
	}

}

double GetPointFlowGradient2(CDC& dc, int x, int y, imatrix& image, Field& gfield, 
							int ker_l, double step_size)
// following the flow, compute the gradient at (x, y)
// Here we simply add the gfield[i][j].mag cumulatively (for repetitive application!)
{
	double d_i, d_j, tx, ty, nx, ny;
	int k, int_i, int_j;
	//int half2;
	int count;

	int half_l = ker_l / 2 + 1; // half length of the curved kernel

	//double half_w = 1.414; // half width of the curved kernel
	//double half_w = 1.414; // half width of the curved kernel
	double half_w = 1; // half width of the curved kernel

	double flow_mag = 0.0; // flow-based gradient magnitude

	////////////////////////////////////
	//t = 1.0;
	//t = 1.414; // step_size
	//t = 0.5; // step_size
	count = 0; // number of pixels traversed
	////////////////////////////////////////////////
	// One half
	d_i = (double)x; d_j = (double)y; 
	tx = -gfield[x][y].gy;
	ty = gfield[x][y].gx;
	int_i = x;	int_j = y;
	for (k = 0; k < half_l; k++) {
		if (tx == 0.0 && ty == 0.0) break;
		///////////////////////////////////////
		nx = tx / sqrt(tx*tx+ty*ty);  // x component of the unit direction vector 
		ny = ty / sqrt(tx*tx+ty*ty);  // y component of the unit direction vector
		////////////////////////////////
		/////////////////////////////////////////
		//flow_mag += ABS( image[right_i][right_j] - image[left_i][left_j] );
		flow_mag += gfield[int_i][int_j].mag;
		count++;
		//////////////////////////////////
		//dc.SetPixelV(int_i, IMAGE_Y-1-int_j, RGB(r, g, b));
		//////////////////////////////
		d_i += nx * step_size; // accurate new location x
		d_j += ny * step_size; // accurate new location y
        ///////////////////////////
		//if (round(d_i) == int_i && round(d_j) == int_j) // no change
		//	continue; // push some more
		/////////////////////////////////////////
		int_i = round(d_i); // integer version of new location x
		int_j = round(d_j); // integer version of new location y
		if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		
		///////////////////////////////////
		//d_i = (double)int_i;
		//d_j = (double)int_j;
		//////////////////////////////////
		tx = -gfield[int_i][int_j].gy;
		ty = gfield[int_i][int_j].gx;
	}
	///////////////////////////////////////
	// other half
	// Note that center pixel gradient magnitude is computed twice (consistent with Sobel grad)
	///*
	d_i = (double)x; d_j = (double)y; 
	int_i = x;	int_j = y;
	tx = -gfield[x][y].gy;
	ty = gfield[x][y].gx;
	/////////////////////////////////////
	for (k = 0; k < half_l; k++) {
		if (tx == 0.0 && ty == 0.0) break;
		///////////////////////////////////////
		nx = tx / sqrt(tx*tx+ty*ty);  // x component of the unit direction vector 
		ny = ty / sqrt(tx*tx+ty*ty);  // y component of the unit direction vector
		////////////////////////////////
		/////////////////////////////////////////
		//flow_mag += ABS( image[right_i][right_j] - image[left_i][left_j] );
		flow_mag += gfield[int_i][int_j].mag;
		count++;
		//////////////////////////////////
		//dc.SetPixelV(int_i, IMAGE_Y-1-int_j, RGB(r, g, b));
		//////////////////////////////
		d_i -= nx * step_size; // accurate new location x
		d_j -= ny * step_size; // accurate new location y
        ///////////////////////////
		//if (round(d_i) == int_i && round(d_j) == int_j) // no change
		//	continue; // push some more
		/////////////////////////////////////////
		int_i = round(d_i); // integer version of new location x
		int_j = round(d_j); // integer version of new location y
		if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		
		//////////////////////////////////
		tx = -gfield[int_i][int_j].gy;
		ty = gfield[int_i][int_j].gx;
	}
	//*/

	if (count > 0)
		flow_mag /= (double)count;

	return flow_mag;
	
}

double GetPointFlowGradient(CDC& dc, int x, int y, imatrix& image, Field& gfield, 
							int ker_l, double ker_w, double step_size)
// This function is the key of my flow-based image processing algorithm!
// following the flow, compute the gradient at (x, y)
// this one keeps the current pixel location as float coordinates, so works correctly!
{
	double d_i, d_j, tx, ty, nx, ny;
	int k, int_i, int_j;
	//int half2;
	int left_i, left_j, right_i, right_j;
	int count;

	int half_l = ker_l / 2 + 1; // half length of the curved kernel

	//double half_w = 1.414; // half width of the curved kernel
	//double half_w = 1.414; // half width of the curved kernel
	//double half_w = 1; // half width of the curved kernel
	double half_w = ker_w / 2; // half width of the curved kernel

	double flow_mag = 0.0; // flow-based gradient magnitude

	/*
	r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
	g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
	b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];
	PerturbRGB(r, g, b, 70);
	*/

	////////////////////////////////////
	//t = 1.0;
	//t = 1.414; // step_size
	//t = 0.5; // step_size
	count = 0; // number of pixels traversed
	////////////////////////////////////////////////
	// One half
	d_i = (double)x; d_j = (double)y; 
	tx = -gfield[x][y].gy;
	ty = gfield[x][y].gx;
	int_i = x;	int_j = y;
	for (k = 0; k < half_l; k++) {
		if (tx == 0.0 && ty == 0.0) break;
		///////////////////////////////////////
		nx = tx / sqrt(tx*tx+ty*ty);  // x component of the unit direction vector 
		ny = ty / sqrt(tx*tx+ty*ty);  // y component of the unit direction vector
		////////////////////////////////
		//////////////////////////////////////////////////
		// Angle computation
		left_i = round(d_i + -ny * half_w); // left perpendicular
		left_j = round(d_j + nx * half_w); // left perpendicular
		if (left_i < 0) left_i = 0;
		if (left_i > IMAGE_X-1) left_i = IMAGE_X-1;
		if (left_j < 0) left_j = 0;
		if (left_j > IMAGE_Y-1) left_j = IMAGE_Y-1;
		right_i = round(d_i + ny * half_w); // right perpendicular
		right_j = round(d_j + -nx * half_w); // right perpendicular
		if (right_i < 0) right_i = 0;
		if (right_i > IMAGE_X-1) right_i = IMAGE_X-1;
		if (right_j < 0) right_j = 0;
		if (right_j > IMAGE_Y-1) right_j = IMAGE_Y-1;
		/////////////////////////////////////////
		flow_mag += ABS( image[right_i][right_j] - image[left_i][left_j] );
		count++;
		//////////////////////////////////
		//dc.SetPixelV(int_i, IMAGE_Y-1-int_j, RGB(r, g, b));
		//////////////////////////////
		d_i += nx * step_size; // accurate new location x
		d_j += ny * step_size; // accurate new location y
        ///////////////////////////
		//if (round(d_i) == int_i && round(d_j) == int_j) // no change
		//	continue; // push some more
		/////////////////////////////////////////
		int_i = round(d_i); // integer version of new location x
		int_j = round(d_j); // integer version of new location y
		if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		
		///////////////////////////////////
		//d_i = (double)int_i;
		//d_j = (double)int_j;
		//////////////////////////////////
		tx = -gfield[int_i][int_j].gy;
		ty = gfield[int_i][int_j].gx;
	}
	///////////////////////////////////////
	// other half
	// Note that center pixel gradient magnitude is computed twice (consistent with Sobel grad)
	///*
	d_i = (double)x; d_j = (double)y; 
	int_i = x;	int_j = y;
	tx = -gfield[x][y].gy;
	ty = gfield[x][y].gx;
	/////////////////////////////////////
	for (k = 0; k < half_l; k++) {
		if (tx == 0.0 && ty == 0.0) break;
		///////////////////////////////////////
		nx = tx / sqrt(tx*tx+ty*ty);  // x component of the unit direction vector 
		ny = ty / sqrt(tx*tx+ty*ty);  // y component of the unit direction vector
		////////////////////////////////
		//////////////////////////////////////////////////
		// Angle computation
		left_i = round(d_i + -ny * half_w); // left perpendicular
		left_j = round(d_j + nx * half_w); // left perpendicular
		if (left_i < 0) left_i = 0;
		if (left_i > IMAGE_X-1) left_i = IMAGE_X-1;
		if (left_j < 0) left_j = 0;
		if (left_j > IMAGE_Y-1) left_j = IMAGE_Y-1;
		right_i = round(d_i + ny * half_w); // right perpendicular
		right_j = round(d_j + -nx * half_w); // right perpendicular
		if (right_i < 0) right_i = 0;
		if (right_i > IMAGE_X-1) right_i = IMAGE_X-1;
		if (right_j < 0) right_j = 0;
		if (right_j > IMAGE_Y-1) right_j = IMAGE_Y-1;
		/////////////////////////////////////////
		flow_mag += ABS( image[right_i][right_j] - image[left_i][left_j] );
		count++;
		//////////////////////////////////
		//dc.SetPixelV(int_i, IMAGE_Y-1-int_j, RGB(r, g, b));
		//////////////////////////////
		d_i -= nx * step_size; // accurate new location x
		d_j -= ny * step_size; // accurate new location y
        ///////////////////////////
		//if (round(d_i) == int_i && round(d_j) == int_j) // no change
		//	continue; // push some more
		/////////////////////////////////////////
		int_i = round(d_i); // integer version of new location x
		int_j = round(d_j); // integer version of new location y
		if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		
		///////////////////////////////////
		//d_i = (double)int_i;
		//d_j = (double)int_j;
		//////////////////////////////////
		tx = -gfield[int_i][int_j].gy;
		ty = gfield[int_i][int_j].gx;
	}
	//*/

	if (count > 0)
		flow_mag /= (double)count;

	return flow_mag;
}

double GetPointFlowGradient3(CDC& dc, int x, int y, imatrix& image, Field& gfield, 
							int ker_l, double ker_w, double step_size)
// following the flow, compute the gradient at (x, y)
// this one keeps the current pixel location as float coordinates
// and this one uses weight, giving higher smoothing weight to higher gradient magnitudes!!!
{
	double d_i, d_j, tx, ty, nx, ny;
	int k, int_i, int_j;
	//int half2;
	int left_i, left_j, right_i, right_j;
	int count;

	int half_l = ker_l / 2 + 1; // half length of the curved kernel

	//double half_w = 1.414; // half width of the curved kernel
	//double half_w = 1.414; // half width of the curved kernel
	//double half_w = 1; // half width of the curved kernel
	double half_w = ker_w / 2; // half width of the curved kernel

	double flow_mag = 0.0; // flow-based gradient magnitude
	double weight, total_weight = 0.0;

	/*
	r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
	g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
	b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];
	PerturbRGB(r, g, b, 70);
	*/

	////////////////////////////////////
	//t = 1.0;
	//t = 1.414; // step_size
	//t = 0.5; // step_size
	count = 0; // number of pixels traversed
	////////////////////////////////////////////////
	// One half
	d_i = (double)x; d_j = (double)y; 
	tx = -gfield[x][y].gy;
	ty = gfield[x][y].gx;
	int_i = x;	int_j = y;
	for (k = 0; k < half_l; k++) {
		if (tx == 0.0 && ty == 0.0) break;
		///////////////////////////////////////
		nx = tx / sqrt(tx*tx+ty*ty);  // x component of the unit direction vector 
		ny = ty / sqrt(tx*tx+ty*ty);  // y component of the unit direction vector
		////////////////////////////////
		//////////////////////////////////////////////////
		// Angle computation
		left_i = round(d_i + -ny * half_w); // left perpendicular
		left_j = round(d_j + nx * half_w); // left perpendicular
		if (left_i < 0) left_i = 0;
		if (left_i > IMAGE_X-1) left_i = IMAGE_X-1;
		if (left_j < 0) left_j = 0;
		if (left_j > IMAGE_Y-1) left_j = IMAGE_Y-1;
		right_i = round(d_i + ny * half_w); // right perpendicular
		right_j = round(d_j + -nx * half_w); // right perpendicular
		if (right_i < 0) right_i = 0;
		if (right_i > IMAGE_X-1) right_i = IMAGE_X-1;
		if (right_j < 0) right_j = 0;
		if (right_j > IMAGE_Y-1) right_j = IMAGE_Y-1;
		/////////////////////////////////////////
		weight = ABS( image[right_i][right_j] - image[left_i][left_j] );
		total_weight += weight;
		flow_mag += weight * ABS( image[right_i][right_j] - image[left_i][left_j] );
		count++;
		//////////////////////////////////
		//dc.SetPixelV(int_i, IMAGE_Y-1-int_j, RGB(r, g, b));
		//////////////////////////////
		d_i += nx * step_size; // accurate new location x
		d_j += ny * step_size; // accurate new location y
        ///////////////////////////
		//if (round(d_i) == int_i && round(d_j) == int_j) // no change
		//	continue; // push some more
		/////////////////////////////////////////
		int_i = round(d_i); // integer version of new location x
		int_j = round(d_j); // integer version of new location y
		if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		
		///////////////////////////////////
		//d_i = (double)int_i;
		//d_j = (double)int_j;
		//////////////////////////////////
		tx = -gfield[int_i][int_j].gy;
		ty = gfield[int_i][int_j].gx;
	}
	///////////////////////////////////////
	// other half
	// Note that center pixel gradient magnitude is computed twice (consistent with Sobel grad)
	///*
	d_i = (double)x; d_j = (double)y; 
	int_i = x;	int_j = y;
	tx = -gfield[x][y].gy;
	ty = gfield[x][y].gx;
	/////////////////////////////////////
	for (k = 0; k < half_l; k++) {
		if (tx == 0.0 && ty == 0.0) break;
		///////////////////////////////////////
		nx = tx / sqrt(tx*tx+ty*ty);  // x component of the unit direction vector 
		ny = ty / sqrt(tx*tx+ty*ty);  // y component of the unit direction vector
		////////////////////////////////
		//////////////////////////////////////////////////
		// Angle computation
		left_i = round(d_i + -ny * half_w); // left perpendicular
		left_j = round(d_j + nx * half_w); // left perpendicular
		if (left_i < 0) left_i = 0;
		if (left_i > IMAGE_X-1) left_i = IMAGE_X-1;
		if (left_j < 0) left_j = 0;
		if (left_j > IMAGE_Y-1) left_j = IMAGE_Y-1;
		right_i = round(d_i + ny * half_w); // right perpendicular
		right_j = round(d_j + -nx * half_w); // right perpendicular
		if (right_i < 0) right_i = 0;
		if (right_i > IMAGE_X-1) right_i = IMAGE_X-1;
		if (right_j < 0) right_j = 0;
		if (right_j > IMAGE_Y-1) right_j = IMAGE_Y-1;
		/////////////////////////////////////////
		weight = ABS( image[right_i][right_j] - image[left_i][left_j] );
		total_weight += weight;
		flow_mag += weight * ABS( image[right_i][right_j] - image[left_i][left_j] );
		count++;
		//////////////////////////////////
		//dc.SetPixelV(int_i, IMAGE_Y-1-int_j, RGB(r, g, b));
		//////////////////////////////
		d_i -= nx * step_size; // accurate new location x
		d_j -= ny * step_size; // accurate new location y
        ///////////////////////////
		//if (round(d_i) == int_i && round(d_j) == int_j) // no change
		//	continue; // push some more
		/////////////////////////////////////////
		int_i = round(d_i); // integer version of new location x
		int_j = round(d_j); // integer version of new location y
		if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
		
		///////////////////////////////////
		//d_i = (double)int_i;
		//d_j = (double)int_j;
		//////////////////////////////////
		tx = -gfield[int_i][int_j].gy;
		ty = gfield[int_i][int_j].gx;
	}
	//*/

	if (count > 0) {
		//flow_mag /= (double)count;
		flow_mag /= total_weight;
	}

	return flow_mag;
	
}

void GetFlowGradient(CDC& dc, int image_x, int image_y, imatrix& image, Field& gfield, imatrix& image2, 
					 int ker_l, double ker_w, double step_size) 
// Flow gradient version
{
	int i, j;
	double MAX_GRADIENT = -1.;
	//double MAX_VAL = -1.;
	//double MAX_VAL = 255.;
	double MAX_VAL = 1020.; // 255 + 2 * 255 + 255 (Maximum possible Sobel value)

	Field tmp;

	tmp.init(image_x, image_y);

	/*
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			if (image[i][j] > MAX_VAL) MAX_VAL = (double)image[i][j];
		}
	}
	TRACE("MAX_VAL = %f\n", MAX_VAL);
	*/
	
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			////////////////////////////////////////////////////////////////
			// Important!: the value of image intensity should be normalized to [0,1]
			//tmp[i][j].gx = (image[i+1][j-1] + 2*(double)image[i+1][j] + image[i+1][j+1] 
			//	- image[i-1][j-1] - 2*(double)image[i-1][j] - image[i-1][j+1]) / MAX_VAL;
			//tmp[i][j].gy = (image[i-1][j+1] + 2*(double)image[i][j+1] + image[i+1][j+1]
			//	- image[i-1][j-1] - 2*(double)image[i][j-1] - image[i+1][j-1]) / MAX_VAL;
			/////////////////////////////////////////////
			// Amplify!!!
			//p[i][j].gx = pow(p[i][j].gx, 2); 
			//p[i][j].gy = pow(p[i][j].gy, 2); 
			//p[i][j].gx = pow(p[i][j].gx, 2); 
			//TRACE("p[i][j].gx = %.1f\n", p[i][j].gx);
			//TRACE("p[i][j].gy = %.1f\n", p[i][j].gy);
			//tmp[i][j].mag = sqrt(tmp[i][j].gx * tmp[i][j].gx + tmp[i][j].gy * tmp[i][j].gy);
			tmp[i][j].mag = GetPointFlowGradient(dc, i, j, image, gfield, ker_l, ker_w, step_size);
			//tmp[i][j].mag = GetPointFlowGradient3(dc, i, j, image, gfield, ker_l, ker_w, step_size);
			//tmp[i][j].mag = GetPointFlowGradient2(dc, i, j, image, gfield, ker_l, step_size);

			if (tmp[i][j].mag > MAX_GRADIENT) {
				MAX_GRADIENT = tmp[i][j].mag;
				//TRACE("p[i][j].gx = %01.f\n", p[i][j].gx);
				//TRACE("p[i][j].gy = %01.f\n", p[i][j].gy);
				//TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);
			}
		}
	}
	
	TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);

	max_grad = round(MAX_GRADIENT);

	// Normalize gradient magnitude
	for (i = 0; i < image_x; i++) { 
		for (j = 0; j < image_y; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			tmp[i][j].mag = tmp[i][j].mag / (double)MAX_GRADIENT; // place it between [0, 1]
			G_mag[i][j] = tmp[i][j].mag; // used for nonmaxima suppresion
			if (tmp[i][j].mag < 0)
				TRACE("tmp[%d][%d].mag = %0.2f\n", i, j, tmp[i][j].mag);

			image2[i][j] = 255 - round(tmp[i][j].mag * 255.);
			////////////////////////////////////////
			// update gfield.mag for repetitive application!
			//gfield[i][j].mag = tmp[i][j].mag; 
			////////////////////////////////////
		}
	}

	//////////////////////////////////////////////////////////////
	/// Amplify the gradients (strong grad -> stronger, weak grad -> weaker)
	/*
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			p[i][j].mag = pow(2.0, p[i][j].mag);
		}
	}
	*/
}


void NonmaximaSuppression(int image_x, int image_y, Field& gfield, matrix& G_mag, 
						  imatrix& image, double hi_thres)
{
	double	gx, gy;
	double	g, g1, g2, g3, g4;
	double	t; // interpolation parameter
	int i, j;

	///////////////////////////////////////
	for (j = 0; j < image_y; j++)
		for (i = 0; i < image_x; i++)
			image[i][j] = 255;
	//////////////////////////////////////////
	// Nonmaxima suppression
	for (j = 1; j < image_y-1; j++) {
		for (i = 1; i < image_x-1; i++) {
			gx = gfield[i][j].gx;
			gy = gfield[i][j].gy;
			g = G_mag[i][j];
			//TRACE("gx = %f, gy = %f, g = %f\n", gx, gy, g);
			//if (gx < 0.01 && gy < 0.01) continue; // not an edge
			if (fabs(gx) >= fabs(gy)) { // close to horizontal (note: gy can be 0)
				t = fabs(gy) / fabs(gx); 
				//TRACE("t = %f\n", t);
				g1 = G_mag[i+1][j]; // right
				g2 = G_mag[i-1][j]; // left
				//TRACE("g1 = %f, g2 = %f\n", g1, g2);
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // right up
					g4 = G_mag[i-1][j-1]; // left down
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i+1][j-1]; // right down
					g4 = G_mag[i-1][j+1]; // left up
				}
			}
			else { // close to vertical (note: gx can be 0)
				t = fabs(gx) / fabs(gy);
				g1 = G_mag[i][j+1]; // up
				g2 = G_mag[i][j-1]; // down
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // up right
					g4 = G_mag[i-1][j-1]; // down left
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i-1][j+1]; // up left
					g4 = G_mag[i+1][j-1]; // down right
				}
			}
			if ( g > ((1-t)*g1 + t*g3) && g > ((1-t)*g2 + t*g4) ) {
				////////////////////////////////////////
				//thin_edge[i][j] = 1; // it's a thin edge
				//image[i][j] = 0; // thin edge (maximum pixel)
				//////////////////////////////////////
				if (g > hi_thres) 
					image[i][j] = 0; // thin edge (maximum pixel)
				else 
					image[i][j] = 255; // thin edge below hi_thres
			}
			else { // non-maximum
				image[i][j] = 255;
			}
			
		}
	}
}

void NonmaximaSuppressionETF(matrix& gmag, ETF& e, imatrix& image, double thres)
// using ETF
{
	double	gx, gy;
	double	g, g1, g2, g3, g4;
	double	t; // interpolation parameter
	int i, j;

	int image_x = e.getRow();
	int image_y = e.getCol();

	///////////////////////////////////////
	for (j = 0; j < image_y; j++)
		for (i = 0; i < image_x; i++)
			image[i][j] = 255;
	//////////////////////////////////////////
	// Nonmaxima suppression
	for (j = 1; j < image_y-1; j++) {
		for (i = 1; i < image_x-1; i++) {
			//////////////////////////////////
			gx = -e[i][j].ty;
			gy = e[i][j].tx;
			g = gmag[i][j];
			//TRACE("gx = %f, gy = %f, g = %f\n", gx, gy, g);
			//if (gx < 0.01 && gy < 0.01) continue; 
			if (gx == 0.0 && gy == 0.0) { image[i][j] = 255; continue; } // not an edge 
			//////////////////////////////////////////////////////
			if (fabs(gx) >= fabs(gy)) { // close to horizontal (note: gy can be 0)
				t = fabs(gy) / fabs(gx); 
				//TRACE("t = %f\n", t);
				g1 = gmag[i+1][j]; // right
				g2 = gmag[i-1][j]; // left
				//TRACE("g1 = %f, g2 = %f\n", g1, g2);
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = gmag[i+1][j+1]; // right up
					g4 = gmag[i-1][j-1]; // left down
				}
				else { // 2 or 4 quadrant
					g3 = gmag[i+1][j-1]; // right down
					g4 = gmag[i-1][j+1]; // left up
				}
			}
			else { // close to vertical (note: gx can be 0)
				t = fabs(gx) / fabs(gy);
				g1 = gmag[i][j+1]; // up
				g2 = gmag[i][j-1]; // down
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = gmag[i+1][j+1]; // up right
					g4 = gmag[i-1][j-1]; // down left
				}
				else { // 2 or 4 quadrant
					g3 = gmag[i-1][j+1]; // up left
					g4 = gmag[i+1][j-1]; // down right
				}
			}
			if ( g < ((1-t)*g1 + t*g3) && g < ((1-t)*g2 + t*g4) ) {
				////////////////////////////////////////
				//thin_edge[i][j] = 1; // it's a thin edge
				//image[i][j] = 0; // thin edge (maximum pixel)
				//////////////////////////////////////
				if (g < thres) 
					image[i][j] = 0; // thin edge (maximum pixel)
				else 
					image[i][j] = 255; // thin edge below hi_thres
			}
			else { // non-maximum
				image[i][j] = 255;
			}
			
		}
	}
}

void Sobel(int image_x, int image_y, matrix& tmp_x, matrix& tmp_y, Image& gradient)
{
	double MAX_GRADIENT = -1.0;

	/*
	for (j = 1; j < image_y - 1; j++) {
		for (i = 1; i < image_x - 1; i++) {
			G_x[i][j] = (tmp_x[i+1][j-1] + 2*tmp_x[i+1][j] + tmp_x[i+1][j+1] 
				- tmp_x[i-1][j-1] - 2*tmp_x[i-1][j] - tmp_x[i-1][j+1]);
			G_y[i][j] = (tmp_y[i-1][j+1] + 2*tmp_y[i][j+1] + tmp_y[i+1][j+1]
				- tmp_y[i-1][j-1] - 2*tmp_y[i][j-1] - tmp_y[i+1][j-1]);
			G_mag[i][j] = sqrt(G_x[i][j] * G_x[i][j] + G_y[i][j] * g_y[i][j]);

			if (G_mag[i][j] > MAX_GRADIENT) {
				MAX_GRADIENT = G_mag[i][j];
				//TRACE("p[i][j].gx = %01.f\n", p[i][j].gx);
				//TRACE("p[i][j].gy = %01.f\n", p[i][j].gy);
				//TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);
			}
		}
	}
	*/
}

void InvertImage(int image_x, int image_y, imatrix& image) 
{
	int x, y;

	//for (y = IMAGE_Y - 1; y >= 0; y--) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			image[x][y] = 255 - image[x][y];
		}
	}
}


void DrawImage(CDC& dc, int image_x, int image_y, Image& image) 
{
	int x, y;
	GLubyte r;

	//for (y = IMAGE_Y - 1; y >= 0; y--) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			r = (GLubyte)image[x][y];
			/// Set Pixel in MemDC
			memDC2.SetPixelV(x, (IMAGE_Y-1)-y, RGB(r, r, r));
		}
	}
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC2, 0, 0, SRCCOPY);
}

void DrawGrayImage(CDC& dc, int image_x, int image_y, imatrix& image) 
{
	int x, y;
	GLubyte r;

	//for (y = IMAGE_Y - 1; y >= 0; y--) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			r = (GLubyte)image[x][y];
			dc.SetPixelV(x, (image_y-1)-y, RGB(r, r, r));
			//////////////////////////////////////////////////////////
			//if (r > 0 && r < 255) dc.SetPixelV(x, (IMAGE_Y-1)-y, RGB(255, 0, 0));
			/////////////////////////////////////////////////////////////
		}
	}
}

void CreateRedZone(imatrix& redzone, double sigma) 
// Find the boundary between white and black
{
	int x, y;
	GLubyte r;

	GaussSmoothSep(redzone, sigma); 

	int image_x = redzone.getRow();
	int image_y = redzone.getCol();

	//for (y = IMAGE_Y - 1; y >= 0; y--) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			////////////////////////
			r = (GLubyte)redzone[x][y];
			if (r > 0 && r < 255) redzone[x][y] = 255;
			else redzone[x][y] = 0;
			/////////////////////////////////////////////////////////////
		}
	}
}


void DrawCondInverseGrayImage(CDC& dc, imatrix& image, imatrix& bw) 
// Conditionally inverse the image
{
	int x, y;
	GLubyte r;

	int image_x = image.getRow();
	int image_y = image.getCol();

	//for (y = IMAGE_Y - 1; y >= 0; y--) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r = (GLubyte)image[x][y];
			if (bw[x][y] == 0) r = 255 - r;
			/// Set Pixel in MemDC
			dc.SetPixelV(x, (image_y-1)-y, RGB(r, r, r));
		}
	}
}

void DrawCondInverseGrayImageRedZone(CDC& dc, imatrix& image, imatrix& bw, imatrix& redzone) 
// Conditionally inverse the image
// If Redzone, do not inverse it!
{
	int x, y;
	GLubyte r;

	int image_x = image.getRow();
	int image_y = image.getCol();

	//for (y = IMAGE_Y - 1; y >= 0; y--) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r = (GLubyte)image[x][y];
			/////////////////////////////////////////////
			//if (redzone[x][y]) r = bw[x][y]; // redzone: keep original bw colors there!
			if (redzone[x][y] && bw[x][y] == 0) r = bw[x][y]; // redzone: only leave black area!
			//if (redzone[x][y]) r = 0; 
			// black area and no redzone, then inverse!
			else if (bw[x][y] == 0) r = 255 - r;
			///////////////////////////////////////////
			/// Set Pixel in MemDC
			dc.SetPixelV(x, (image_y-1)-y, RGB(r, r, r));
		}
	}
}

void DrawInverseGrayImage(CDC& dc, int image_x, int image_y, imatrix& image) 
{
	int x, y;
	GLubyte r;

	//for (y = IMAGE_Y - 1; y >= 0; y--) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			r = 255-(GLubyte)image[x][y];
			/// Set Pixel in MemDC
			dc.SetPixelV(x, (IMAGE_Y-1)-y, RGB(r, r, r));
		}
	}
}


void DrawColorImage(CDC& dc, int image_x, int image_y, cimatrix& image) 
{
	int x, y;
	GLubyte r, g, b;

	//for (y = IMAGE_Y - 1; y >= 0; y--) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			r = (GLubyte)image[x][y].r;
			g = (GLubyte)image[x][y].g;
			b = (GLubyte)image[x][y].b;
			/// Set Pixel in MemDC
			dc.SetPixelV(x, (image_y-1)-y, RGB(r, g, b));
		}
	}
}

void ColorQuantization(cimatrix& cmap, int level) 
// HSV color model
// Based on Gooch's method
{
	int x, y;
	GLubyte r1, g1, b1;
	double r, g, b;
	double h, s, v;
	double t;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	//imatrix tmp(image_x, image_y);

	//CopyCol2GrayImage(image_x, image_y, cmap, tmp);

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r = (double)cmap[x][y].r;
			g = (double)cmap[x][y].g;
			b = (double)cmap[x][y].b;
			r /= 255.;
			g /= 255.;
			b /= 255.;
			//TRACE("r = %f\n", r);
			//TRACE("g = %f\n", g);
			//TRACE("b = %f\n", b);
			RGB2HSV(r, g, b, h, s, v);
			///////////////////////////
			v = ( (int)(v * level) ) / (double)level + 1.0/(double)level/2.0;
			///////////////////////////////
			////////////////////////////////
			// extrapolate! s and v
			t = 1.1;
			s = (1-t)*0.0 + t*s;
			v = (1-t)*0.0 + t*v;
			if (s > 1.0) s = 1.0;
			if (v > 1.0) v = 1.0;
			///////////////////////////////
			HSV2RGB(h, s, v, r, g, b);
			r1 = (GLubyte)(r*255);
			g1 = (GLubyte)(g*255);
			b1 = (GLubyte)(b*255);
			cmap[x][y].r = r1;
			cmap[x][y].g = g1;
			cmap[x][y].b = b1;
		}
	}
}

void ColorQuantization2(cimatrix& cmap, int level) 
// HSV color model
// We quantize not just V but also H
{
	int x, y;
	GLubyte r1, g1, b1;
	double r, g, b;
	double h, s, v;
	//double t;

	int level1 = 10;
	//int level2 = 20; // bannana
	int level2 = 20; // eagle
	int level3 = 10;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	//imatrix tmp(image_x, image_y);

	//CopyCol2GrayImage(image_x, image_y, cmap, tmp);

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r = (double)cmap[x][y].r;
			g = (double)cmap[x][y].g;
			b = (double)cmap[x][y].b;
			r /= 255.;
			g /= 255.;
			b /= 255.;
			//TRACE("r = %f\n", r);
			//TRACE("g = %f\n", g);
			//TRACE("b = %f\n", b);
			RGB2HSV(r, g, b, h, s, v);
			///////////////////////////
			if (v == 0.0) v = 0.0;
			else if (v == 1.0) v = 1.0;
			v = ( (int)(v * level1) ) / (double)level1 + 1.0/(double)level1/2.0;
			h /= 360.0; // [0, 1]
			if (h == 0.0) h = 0.0;
			else if (h == 1.0) h = 1.0;
			else h = ( (int)(h * level2) ) / (double)level2 + 1.0/(double)level2/2.0;
			if (s == 0.0) s = 0.0;
			else if (s == 1.0) s = 1.0;
			else s = ( (int)(s * level3) ) / (double)level3 + 1.0/(double)level3/2.0;
			///////////////////////////////
			if (h > 1.0) h = 1.0;
			if (h < 0.0) h = 0.0;
			h *= 360.0; // [0, 360]
			////////////////////////////////
			// extrapolate! s and v
			double t = 1.1;
			//if (h > 360.0) h = 360.0;
			s = (1-t)*0.0 + t*s;
			if (s > 1.0) s = 1.0;
			if (s < 0.0) s = 0.0;
			//s = 1.0; // saturation always maximum!
			v = (1-t)*0.0 + t*v;
			if (v > 1.0) v = 1.0;
			if (v < 0.0) v = 0.0;
			///////////////////////////////
			HSV2RGB(h, s, v, r, g, b);
			r1 = (GLubyte)(r*255);
			g1 = (GLubyte)(g*255);
			b1 = (GLubyte)(b*255);
			cmap[x][y].r = r1;
			cmap[x][y].g = g1;
			cmap[x][y].b = b1;
		}
	}
}

void ColorQuantizationRGB(cimatrix& cmap, int level) 
// RGB color model
{
	int x, y;
	GLubyte r1, g1, b1;
	double r, g, b;
	//double t;

	int level1 = level;
	int level2 = level; // eagle
	int level3 = level;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	//imatrix tmp(image_x, image_y);

	//CopyCol2GrayImage(image_x, image_y, cmap, tmp);

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r = (double)cmap[x][y].r;
			g = (double)cmap[x][y].g;
			b = (double)cmap[x][y].b;
			r /= 255.;
			g /= 255.;
			b /= 255.;
			//TRACE("r = %f\n", r);
			//TRACE("g = %f\n", g);
			//TRACE("b = %f\n", b);
			//RGB2HSV(r, g, b, h, s, v);
			///////////////////////////
			if (r == 0.0) r = 0.0;
			else if (r == 1.0) r = 1.0;
			r = ( (int)(r * level1) ) / (double)level1 + 1.0/(double)level1/2.0;
			if (g == 0.0) g = 0.0;
			else if (g == 1.0) g = 1.0;
			else g = ( (int)(g * level2) ) / (double)level2 + 1.0/(double)level2/2.0;
			if (b == 0.0) b = 0.0;
			else if (b == 1.0) b = 1.0;
			else b = ( (int)(b * level3) ) / (double)level3 + 1.0/(double)level3/2.0;
			///////////////////////////////
			//if (h > 1.0) h = 1.0;
			//if (h < 0.0) h = 0.0;
			//h *= 360.0; // [0, 360]
			////////////////////////////////
			// extrapolate! s and v
			//double t = 1.1;
			//if (h > 360.0) h = 360.0;
			//s = (1-t)*0.0 + t*s;
			//if (s > 1.0) s = 1.0;
			//if (s < 0.0) s = 0.0;
			//s = 1.0; // saturation always maximum!
			//v = (1-t)*0.0 + t*v;
			//if (v > 1.0) v = 1.0;
			//if (v < 0.0) v = 0.0;
			///////////////////////////////
			r1 = (GLubyte)(r*255);
			g1 = (GLubyte)(g*255);
			b1 = (GLubyte)(b*255);
			cmap[x][y].r = r1;
			cmap[x][y].g = g1;
			cmap[x][y].b = b1;
		}
	}
}

void GrayQuantization(imatrix& image, int level) 
// Based on Gooch's method
{
	int x, y;
	double v;
	double t;

	int image_x = image.getRow();
	int image_y = image.getCol();

	//imatrix tmp(image_x, image_y);

	//CopyCol2GrayImage(image_x, image_y, cmap, tmp);

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			v = (double)image[x][y];
			v /= 255.;
			///////////////////////////
			v = ( (int)(v * level) ) / (double)level + 1.0/(double)level/2.0;
			///////////////////////////////
			// extrapolate! s and v
			t = 1.1;
			//s = (1-t)*0.0 + t*s;
			v = (1-t)*0.0 + t*v;
			//if (s > 1.0) s = 1.0;
			if (v > 1.0) v = 1.0;
			///////////////////////////////
			image[x][y] = round(v * 255.);
		}
	}
}


typedef struct node_item {
	int val;
	int count;
	struct node_item *next;
} node;

void Dilation(imatrix& image, int itr) 
// check neighboring pixels
{
	int i, j, val;

	int image_x = image.getRow();
	int image_y = image.getCol();

	node *list, *p;
	bool found;

	list = NULL;
	
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			found = false;
			for (p = list; p; p=p->next) {
				if (p->val == image[i][j]) {
					found = true; break;
				}
			}
			if (found) p->count++;
			else {
				p = (node *)malloc(sizeof(node));
				p->count = 0;
				p->val = image[i][j];
				p->next = list; // add to the list (at the head)
				list = p;
			}
		}
	}

	int max_count = 0;
	node *max_p = NULL;
	for (p = list; p; p=p->next) {
		if (p->val > 240 || p->val < 30) continue;
		if (p->count > max_count) {
			max_count = p->count; 
			max_p = p;
		}
	}

	//matrix Ix, Iy;
	//Ix.init(image_x, image_y);
	//Iy.init(image_x, image_y);

	imatrix tmp(image_x, image_y);
	tmp.copy(image);

	int x, y, m;

	int thres = 100, s, t;
	bool processed;

	for (m = 0; m < itr; m++) {
		//////////////////////////////
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				///////////////////////////////////
				//if (Ix[i][j] == 0.0 && Iy[i][j] == 0.0) continue;
				val = tmp[i][j];
				if (val <= 30) continue; // stop at lines
				if ( ABS(val - max_p->val) > thres ) continue; // too different colors
				//TRACE("val = %d\n", val);
				if (val != max_p->val) {
					processed = false;
					for (t = -1; t <= 1; t++) {
						for (s = -1; s <= 1; s++) {
							x = i + s; y = j + t;
							if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
							if (tmp[x][y] == max_p->val) {
								image[i][j] = max_p->val;
								processed = true;
							}
							if (processed) break;
						}
						if (processed) break;
					}
				}
			} // for i
		} // for j
		tmp.copy(image);
	} // for m

}

typedef struct nodeRGB_item {
	iRGB col;
	int count;
	struct nodeRGB_item *next;
} nodeRGB;

void DilationColor(CDC& dc, cimatrix& cmap) 
// check neighboring pixels
{
	int i, j;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	nodeRGB *list, *p, *q;
	iRGB col;

	bool found;

	list = NULL;
	
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			found = false;
			for (p = list; p; p=p->next) { // already in the list
				if (p->col.r == cmap[i][j].r && 
					p->col.g == cmap[i][j].g && 
					p->col.b == cmap[i][j].b) {
					found = true; break;
				}
			}
			if (found) p->count++;
			else { // newly insert
				p = (nodeRGB *)malloc(sizeof(nodeRGB));
				p->count = 0;
				p->col = cmap[i][j];
				p->next = list; // add to the list (at the head)
				list = p;
			}
		}
	}

	

	//TRACE("max_p->col = [%d, %d, %d]\n", max_p->col.r, max_p->col.g, max_p->col.b);
	//TRACE("max_p->count = %d\n", max_p->count);

	cimatrix tmp(image_x, image_y);
	tmp.copy(cmap);

	int x, y;

	int s, t;
	bool processed, remain;
	double dist;
	double thres = 100;

	int max_count;
	nodeRGB *max_p;

	//for (m = 0; m < itr; m++) {
	while (1) {
		/////////////////////////////////////
		// find max node
		max_count = 0;
		max_p = NULL;
		for (p = list; p; p=p->next) {
			//if (p->val > 240 || p->val < 30) continue;
			if (p->col.r < 10 && p->col.g < 10 && p->col.b < 10) continue; // line
			if (p->count > max_count) {
				max_count = p->count; 
				max_p = p;
			}
		}
		if (max_p == NULL) break; // exhausted
		/////////////////////////////
		thres *= 0.9;
		//////////////////////////////
		// Dilation!
		while (1) {
			remain = false;
			for (j = 0; j < image_y; j++) {
				for (i = 0; i < image_x; i++) {
					///////////////////////////////////
					//if (Ix[i][j] == 0.0 && Iy[i][j] == 0.0) continue;
					col = tmp[i][j];
					if (col.r < 10 && col.g < 10 && col.b < 10) continue; // stop at lines
					dist = dist3(col.r-max_p->col.r, col.g-max_p->col.g, col.b-max_p->col.b);
					if ( dist > thres ) continue; // colors are too different; stop
					//TRACE("val = %d\n", val);
					if (col.r != max_p->col.r || col.g != max_p->col.g || col.b != max_p->col.b) {
						processed = false;
						for (t = -1; t <= 1; t++) {
							for (s = -1; s <= 1; s++) {
								x = i + s; y = j + t;
								if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
								if (tmp[x][y].r == max_p->col.r && 
									tmp[x][y].g == max_p->col.g && 
									tmp[x][y].b == max_p->col.b) {
										////////////////////////////
										cmap[i][j] = max_p->col;
										processed = true;
										remain = true;
								}
								if (processed) break;
							}
							if (processed) break;
						}
					}
				} // for i
			} // for j
			if (!remain) break;
			tmp.copy(cmap);
			//////////////////////////////
			DrawColorImage(memDC, image_x, image_y, cmap);
			dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		}
		/////////////////////////////////////
		// delete max node
		for (p=list, q=NULL; p; q=p, p=p->next) {
			//if (p->val > 240 || p->val < 30) continue;
			if (p == max_p) {
				if (q == NULL) { // max_p is the first node
					list = p->next; break;
				}
				q->next = p->next;
				break;
			}
		}
		//////////////////////////////////
		free(max_p);
		////////////////////////////////////////
	} // for m

}

void DilationColor2(CDC& dc, cimatrix& cmap, imatrix& line, double thres, double factor) 
// use line maps as stopping criteria
{
	int i, j;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	nodeRGB *list, *p, *q;
	iRGB col;

	bool found;

    list = NULL;
	//list2 = NULL;
	
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			found = false;
			for (p = list; p; p=p->next) { // already in the list
				if (p->col.r == cmap[i][j].r && 
					p->col.g == cmap[i][j].g && 
					p->col.b == cmap[i][j].b) {
					found = true; break;
				}
			}
			if (found) p->count++;
			else { // newly insert
				p = (nodeRGB *)malloc(sizeof(nodeRGB));
				p->count = 0;
				p->col = cmap[i][j];
				p->next = list; // add to the list (at the head)
				list = p;
			}
		}
	}

	//TRACE("max_p->col = [%d, %d, %d]\n", max_p->col.r, max_p->col.g, max_p->col.b);
	//TRACE("max_p->count = %d\n", max_p->count);

	cimatrix tmp(image_x, image_y);
	tmp.copy(cmap);

	imatrix done(image_x, image_y);
	done.zero();

	int x, y;

	int s, t;
	bool processed, remain;
	double dist;
	//double thres = 100;

	int max_count;
	nodeRGB *max_p;

	//for (m = 0; m < itr; m++) {
	while (1) {
		/////////////////////////////////////
		// find max node
		max_count = 0;
		max_p = NULL;
		for (p = list; p; p=p->next) {
			//if (p->val > 240 || p->val < 30) continue;
			//if (p->col.r < 10 && p->col.g < 10 && p->col.b < 10) continue; // line
			if (p->count > max_count) {
				max_count = p->count; 
				max_p = p;
			}
		}
		if (max_p == NULL) break; // exhausted
		/////////////////////////////
		thres *= factor;
		//////////////////////////////
		// Dilation!
		while (1) {
			remain = false;
			for (j = 0; j < image_y; j++) {
				for (i = 0; i < image_x; i++) {
					///////////////////////////////////
					//if (Ix[i][j] == 0.0 && Iy[i][j] == 0.0) continue;
					if (done[i][j]) continue;
					col = tmp[i][j];
					if (col.r == max_p->col.r && col.g == max_p->col.g && col.b == max_p->col.b) {
						done[i][j] = 1; continue;
					}
					//if (col.r < 10 && col.g < 10 && col.b < 10) continue; // stop at lines
					if (line[i][j] == 0) continue; // stop at lines
					//////////////////////////////////////////////
					dist = dist3(col.r-max_p->col.r, col.g-max_p->col.g, col.b-max_p->col.b);
					if ( dist > thres ) continue; // colors are too different; stop
					//TRACE("val = %d\n", val);
					if (col.r != max_p->col.r || col.g != max_p->col.g || col.b != max_p->col.b) {
						processed = false;
						for (t = -1; t <= 1; t++) {
							for (s = -1; s <= 1; s++) {
								x = i + s; y = j + t;
								if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
								if (tmp[x][y].r == max_p->col.r && 
									tmp[x][y].g == max_p->col.g && 
									tmp[x][y].b == max_p->col.b) {
										////////////////////////////
										cmap[i][j] = max_p->col;
										done[i][j] = 1;
										processed = true;
										remain = true;
								}
								if (processed) break;
							}
							if (processed) break;
						}
					}
				} // for i
			} // for j
			if (!remain) break;
			tmp.copy(cmap);
			//////////////////////////////
			DrawColorImage(memDC, image_x, image_y, cmap);
			dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		}
		/////////////////////////////////////
		// delete max node
		for (p=list, q=NULL; p; q=p, p=p->next) {
			//if (p->val > 240 || p->val < 30) continue;
			if (p == max_p) {
				if (q == NULL) { // max_p is the first node
					list = p->next; break;
				}
				q->next = p->next;
				break;
			}
		}
		//free(max_p);
		// move to list2
		//if (list2 == NULL) max_p->next = NULL;
		//else max_p->next = list2;
		///////////////////////////////////////////////////
	} // for m

	/*
	for (p=list2; p; p=q) {
		q = p->next;
		free(p);
	}
	*/
}

void GetFlowDilation(imatrix& line, ETF& e, int half_l, int max_itr)
// flow-based Morphological dilation
// extend the black line along flow direction
{
	int k, m;
	int i, j, i_x, i_y;

	vector vt(2);
	double x, y, d_x, d_y;
	int x1, y1;

	int image_x = line.getRow();
	int image_y = line.getCol();

	imatrix tmp;
	tmp.copy(line);

	double step_size;
	bool found;
	int val;

	step_size = 1.0;

	for (m = 0; m < max_itr; m++) {
		/////////////////////////////////////////
		// FBL along main axis!
		//StartTimer();
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				/////////////////////////////
				if (line[i][j] == 0) continue; // already black
				////////////////////////////////
				found = false;
				////////////////////////////////////////////////
				// One half
				d_x = (double)i; d_y = (double)j; 
				i_x = i;	i_y = j;
				////////////////////////////////////////////////
				for (k = 0; k < half_l; k++) {
					vt[0] = e[i_x][i_y].tx;
					vt[1] = e[i_x][i_y].ty;
					if (vt[0] == 0.0 && vt[1] == 0.0) break;
					//vt.make_unit();
					////////////////////////////////////////////////
					x = d_x;
					y = d_y;
					/////////////////////////////////////////////////////
					if (x > (double)image_x-1 || x < 0.0 || y > (double)image_y-1 || y < 0.0) 
						continue;
					/////////////////////////////////////////////////
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > image_x-1) x1 = image_x-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > image_y-1) y1 = image_y-1;
					//////////////////////////////
					val = tmp[x1][y1];
					if (val == 0) { found = true; break; }
					/////////////////////////////////////////////////////////
					d_x += vt[0] * step_size; // accurate new location x
					d_y += vt[1] * step_size; // accurate new location y
					/////////////////////////////////////////
					i_x = round(d_x); // integer version of new location x
					i_y = round(d_y); // integer version of new location y
					//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
					if (d_x < 0 || d_x > image_x-1 || d_y < 0 || d_y > image_y-1) break;
				}
				/////////////////////////////
				if (found) line[i][j] = 0;
				////////////////////////////////////////////////
				// Other half
				d_x = (double)i; d_y = (double)j; 
				i_x = i;	i_y = j;
				//////////////////////////////////////
				for (k = 0; k < half_l; k++) {
					vt[0] = -e[i_x][i_y].tx;
					vt[1] = -e[i_x][i_y].ty;
					if (vt[0] == 0.0 && vt[1] == 0.0) break;
					//vt.make_unit();
					///////////////////////////////
					x = d_x;
					y = d_y;
					/////////////////////////////////////////////////////
					if (x > (double)image_x-1 || x < 0.0 || y > (double)image_y-1 || y < 0.0) 
						continue;
					/////////////////////////////////////////////////
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > image_x-1) x1 = image_x-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > image_y-1) y1 = image_y-1;
					val = tmp[x1][y1];
					if (val == 0) { found = true; break; }
					/////////////////////////////////////////
					d_x += vt[0] * step_size; // accurate new location x
					d_y += vt[1] * step_size; // accurate new location y
					/////////////////////////////////////////
					i_x = round(d_x); // integer version of new location x
					i_y = round(d_y); // integer version of new location y
					//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
					if (d_x < 0 || d_x > image_x-1 || d_y < 0 || d_y > image_y-1) break;
				}
				//////////////////////////////////
				if (found) line[i][j] = 0;
			}
		}
		tmp.copy(line); // update tmp
	}
}

iRGB *CreatePaletteFromImage(cimatrix& cmap, int& N) 
{
	int i, j;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	nodeRGB *list, *p, *q;
    
	bool found;

    list = NULL;
	
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			found = false;
			for (p = list; p; p=p->next) { // already in the list
				if (p->col.r == cmap[i][j].r && 
					p->col.g == cmap[i][j].g && 
					p->col.b == cmap[i][j].b) {
					found = true; break;
				}
			}
			if (found) p->count++;
			else { // newly insert
				p = (nodeRGB *)malloc(sizeof(nodeRGB));
				p->count = 0;
				p->col = cmap[i][j];
				p->next = list; // add to the list (at the head)
				list = p;
			}
		}
	}

	N = 0;
	for (p = list; p; p=p->next) N++;

	TRACE("N = %d\n", N);

	iRGB *pal;
	pal = new iRGB[N];

	for (p=list, i=0; p; p=p->next, i++) {
		pal[i].r = p->col.r;
		pal[i].g = p->col.g;
		pal[i].b = p->col.b;
		//TRACE("pal[%d] = [%d, %d, %d]\n", p->col.r, p->col.g, p->col.b);
		//TRACE("pal[%d] = [%d, %d, %d]\n", i, pal[i].r, pal[i].g, pal[i].b);
	}

	for (p=list; p; p=q) {
		q = p->next;
		free(p);
	}

	return pal;
}

void ColorQuantizationUsingPalette(cimatrix& cmap, iRGB* pal, int N) 
// use only the colors in the palette
{
	int i, j, k;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	iRGB col, p_col, min_col;
	double dist, min_dist;

	cimatrix tmp(image_x, image_y);
	tmp.copy(cmap);

	TRACE("N = %d\n", N);

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			col = tmp[i][j];
			min_dist = 1000000000;
			for (k = 0; k < N; k++) {
				p_col = pal[k];
				//TRACE("pal[%d] = [%d, %d, %d]\n", k, pal[i].r, pal[i].g, pal[i].b);
				dist = dist3(col.r-p_col.r, col.g-p_col.g, col.b-p_col.b);
				if (dist < min_dist) {
					min_dist = dist;
					min_col = p_col;
				}
			}
			cmap[i][j] = min_col;
		}
	}

}

void ClosingColor(CDC& dc, cimatrix& cmap, imatrix& line, int itr) 
// use line maps as stopping criteria
{
	int i, j;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	nodeRGB *list, *p, *q;
	iRGB col;

	bool found;

    list = NULL;
	//list2 = NULL;
	
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			found = false;
			for (p = list; p; p=p->next) { // already in the list
				if (p->col.r == cmap[i][j].r && 
					p->col.g == cmap[i][j].g && 
					p->col.b == cmap[i][j].b) {
					found = true; break;
				}
			}
			if (found) p->count++;
			else { // newly insert
				p = (nodeRGB *)malloc(sizeof(nodeRGB));
				p->count = 0;
				p->col = cmap[i][j];
				p->next = list; // add to the list (at the head)
				list = p;
			}
		}
	}

	//TRACE("max_p->col = [%d, %d, %d]\n", max_p->col.r, max_p->col.g, max_p->col.b);
	//TRACE("max_p->count = %d\n", max_p->count);

	cimatrix tmp(image_x, image_y);
	tmp.copy(cmap);

	imatrix done(image_x, image_y);
	done.zero();

	int x, y;

	int s, t, m;
	bool processed;
	//bool remain;
	//double dist;
	//double thres = 100;

	int max_count;
	nodeRGB *max_p;

	//for (m = 0; m < itr; m++) {
	while (1) {
		/////////////////////////////////////
		// find max node
		max_count = 0;
		max_p = NULL;
		for (p = list; p; p=p->next) {
			//if (p->val > 240 || p->val < 30) continue;
			//if (p->col.r < 10 && p->col.g < 10 && p->col.b < 10) continue; // line
			if (p->count > max_count) {
				max_count = p->count; 
				max_p = p;
			}
		}
		if (max_p == NULL) break; // exhausted
		/////////////////////////////
		//thres *= factor;
		//////////////////////////////
		// Dilation!
		//while (1) {
		for (m = 0; m < itr; m++) {
			//remain = false;
			for (j = 0; j < image_y; j++) {
				for (i = 0; i < image_x; i++) {
					///////////////////////////////////
					//if (Ix[i][j] == 0.0 && Iy[i][j] == 0.0) continue;
					//if (done[i][j]) continue;
					if (line[i][j] == 0) continue; // stop at lines
					col = tmp[i][j];
					if (col.r == max_p->col.r && col.g == max_p->col.g && col.b == max_p->col.b) {
						//done[i][j] = 1; 
						continue;
					}
					//////////////////////////////////////////////
					//dist = dist3(col.r-max_p->col.r, col.g-max_p->col.g, col.b-max_p->col.b);
					//if ( dist > thres ) continue; // colors are too different; stop
					processed = false;
					for (t = -1; t <= 1; t++) {
						for (s = -1; s <= 1; s++) {
							x = i + s; y = j + t;
							if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
							if (tmp[x][y].r == max_p->col.r && 
								tmp[x][y].g == max_p->col.g && 
								tmp[x][y].b == max_p->col.b) {
									////////////////////////////
									cmap[i][j] = max_p->col;
									processed = true;
									//remain = true;
							}
							if (processed) break;
						}
						if (processed) break;
					}
				} // for i
			} // for j
			//if (!remain) break;
			tmp.copy(cmap);
			//////////////////////////////////////////////
			// Erosion!
			for (j = 0; j < image_y; j++) {
				for (i = 0; i < image_x; i++) {
					///////////////////////////////////
					if (line[i][j] == 0) continue; // stop at lines
					col = tmp[i][j];
					if (col.r != max_p->col.r || col.g != max_p->col.g || col.b != max_p->col.b) {
						//done[i][j] = 1; 
						continue;
					}
					//////////////////////////////////////////////
					//dist = dist3(col.r-max_p->col.r, col.g-max_p->col.g, col.b-max_p->col.b);
					//if ( dist > thres ) continue; // colors are too different; stop
					processed = false;
					for (t = -1; t <= 1; t++) {
						for (s = -1; s <= 1; s++) {
							x = i + s; y = j + t;
							if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
							if (tmp[x][y].r != max_p->col.r || 
								tmp[x][y].g != max_p->col.g || 
								tmp[x][y].b != max_p->col.b) {
									////////////////////////////
									cmap[i][j] = tmp[x][y];
									processed = true;
									//remain = true;
							}
							if (processed) break;
						}
						if (processed) break;
					}
				} // for i
			} // for j
			//if (!remain) break;
			tmp.copy(cmap);
			//////////////////////////////
			DrawColorImage(memDC, image_x, image_y, cmap);
			dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
			TRACE("opening closing iteration %d\n", m);
		} // for m
		/////////////////////////////////////
		// delete max node
		for (p=list, q=NULL; p; q=p, p=p->next) {
			//if (p->val > 240 || p->val < 30) continue;
			if (p == max_p) {
				if (q == NULL) { // max_p is the first node
					list = p->next; break;
				}
				q->next = p->next;
				break;
			}
		}
		//free(max_p);
		// move to list2
		//if (list2 == NULL) max_p->next = NULL;
		//else max_p->next = list2;
		///////////////////////////////////////////////////
	} // for m

	/*
	for (p=list2; p; p=q) {
		q = p->next;
		free(p);
	}
	*/
}

void ErosionGray(CDC& dc, imatrix& gray, int itr) 
{
	int i, j;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	imatrix tmp(image_x, image_y);
	tmp.copy(gray);

	int x, y;

	int s, t, m;
	bool processed;

	for (m = 0; m < itr; m++) {
		//remain = false;
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				///////////////////////////////////
				if (tmp[i][j] != 0) continue;
				processed = false;
				for (t = -1; t <= 1; t++) {
					for (s = -1; s <= 1; s++) {
						x = i + s; y = j + t;
						if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
						if (tmp[x][y] == 255) {
							////////////////////////////
							gray[i][j] = 255;
							processed = true;
							//TRACE("processed!\n");
						}
						if (processed) break;
					}
					if (processed) break;
				}
			} // for i
		} // for j
		tmp.copy(gray);
		//////////////////////////////
		//DrawGrayImage(memDC, image_x, image_y, gray);
		//dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		//TRACE("Erosion iteration %d\n", m);
	} // for m
}

void ErosionGrayCircle(CDC& dc, imatrix& gray, int itr) 
// use circle structuring element
{
	int i, j;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	imatrix tmp(image_x, image_y);
	tmp.copy(gray);

	int x, y;

	int s, t, m;
	bool processed;

	int size = 2;

	int mask[5][5] = { 
		{0, 1, 1, 1, 0},
		{1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1},
		{0, 1, 1, 1, 0}
	};

	/*
	int mask[3][3] = { 
		{0, 1, 0},
		{1, 1, 1},
		{0, 1, 0}
	};
	*/

	for (m = 0; m < itr; m++) {
		//remain = false;
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				///////////////////////////////////
				if (tmp[i][j] != 0) continue;
				processed = false;
				for (t = -size; t <= size; t++) {
					for (s = -size; s <= size; s++) {
						x = i + s; y = j + t;
						if (!mask[s+size][t+size]) continue;
						if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
						if (tmp[x][y] == 255) {
							////////////////////////////
							gray[i][j] = 255;
							processed = true;
							//TRACE("processed!\n");
						}
						if (processed) break;
					}
					if (processed) break;
				}
			} // for i
		} // for j
		tmp.copy(gray);
		//////////////////////////////
		//DrawGrayImage(memDC, image_x, image_y, gray);
		//dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		//TRACE("Erosion iteration %d\n", m);
	} // for m
}


void RemoveIsolatedRegions(CDC& dc, imatrix& gray, int size_max) 
// using morphology, remove isolated regions
// start with the smalled possible isolated region, then test for bigger regions
{
	int i, j;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	imatrix tmp(image_x, image_y);
	tmp.copy(gray);

	int x, y;
	//int size = 1;

	int s, t;
	bool isolated;

	int size = 0;

	while (1) {
		//remain = false;
		size++;
		if (size > size_max) break;
		/////////////////////////////////
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				///////////////////////////////////
				if (tmp[i][j] != 0) continue; // ignore white pixels
				isolated = true;
				for (t = -size; t <= size; t++) { // left border
					x = i - size; y = j + t;
					if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
					if (tmp[x][y] == 0) { // we found neighboring black pixel
						isolated = false;
					}
					if (!isolated) break;
				}
				if (!isolated) continue;
				/////////////////////////////////
				for (t = -size; t <= size; t++) { // right border
					x = i + size; y = j + t;
					if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
					if (tmp[x][y] == 0) { // we found neighboring black pixel
						isolated = false;
					}
					if (!isolated) break;
				}
				if (!isolated) continue;
				for (s = -size; s <= size; s++) { // left border
					x = i + s; y = j - size;
					if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
					if (tmp[x][y] == 0) { // we found neighboring black pixel
						isolated = false;
					}
					if (!isolated) break;
				}
				if (!isolated) continue;
				/////////////////////////////////
				for (s = -size; s <= size; s++) { // right border
					x = i + s; y = j + size;
					if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
					if (tmp[x][y] == 0) { // we found neighboring black pixel
						isolated = false;
					}
					if (!isolated) break;
				}
				if (!isolated) continue;
				/////////////////////////////
				// isolated, confirmed
				gray[i][j] = 255;
			} // for i
		} // for j
		tmp.copy(gray);
		//////////////////////////////
		//DrawGrayImage(memDC, image_x, image_y, gray);
		//dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		//TRACE("Erosion iteration %d\n", m);
	} // for m
}

void DilationGray(CDC& dc, imatrix& gray, int itr) 
{
	int i, j;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	imatrix tmp(image_x, image_y);
	tmp.copy(gray);

	int x, y;

	int s, t, m;
	bool processed;

	for (m = 0; m < itr; m++) {
		//remain = false;
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				///////////////////////////////////
				if (tmp[i][j] == 0) continue;
				processed = false;
				for (t = -1; t <= 1; t++) {
					for (s = -1; s <= 1; s++) {
						x = i + s; y = j + t;
						if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
						if (tmp[x][y] == 0) {
							////////////////////////////
							gray[i][j] = 0;
							processed = true;
							//TRACE("processed!\n");
						}
						if (processed) break;
					}
					if (processed) break;
				}
			} // for i
		} // for j
		tmp.copy(gray);
		//////////////////////////////
		//DrawGrayImage(memDC, image_x, image_y, gray);
		//dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		//TRACE("Erosion iteration %d\n", m);
	} // for m
}

void GetOffsetLines(CDC& dc, imatrix& gray) 
{
	int i, j, k;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	imatrix tmp(image_x, image_y);
	tmp.copy(gray);

	deque<PixeL> pnts;
	pnts.clear();

	PixeL p;

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			if (tmp[i][j] == 0) {
				p.x = i; p.y = j;
				pnts.push_back(p);
			}
		}
	}

	TRACE("pnts.size() = %d\n", pnts.size());

	double max_dist = dist2(0, 0, image_x-1, image_y-1);
	double min_dist, dist;

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			min_dist = 1000000;
			if (tmp[i][j] != 0) {
				for (k = 0; k < (signed)pnts.size(); k++) {
					dist = dist2(i, j, pnts[k].x, pnts[k].y);
					if (dist < min_dist) min_dist = dist;
				}
				gray[i][j] = round( (min_dist / max_dist) * 255.0 );
			}
			DrawGrayImage(memDC, image_x, image_y, gray);
			dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		}
	}
	
}

void DilationGraySep(CDC& dc, imatrix& gray, int itr) 
// separable version
{
	int i, j;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	imatrix tmp(image_x, image_y);
	tmp.copy(gray);

	int x, y;

	int s, t, m;
	bool processed;

	for (m = 0; m < itr; m++) {
		//remain = false;
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				///////////////////////////////////
				if (tmp[i][j] == 0) continue;
				processed = false;
				for (s = -1; s <= 1; s++) {
					x = i + s; y = j;
					if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
					if (tmp[x][y] == 0) {
						////////////////////////////
						gray[i][j] = 0;
						processed = true;
						//TRACE("processed!\n");
					}
					if (processed) break;
				}
			} // for i
		} // for j
		tmp.copy(gray);

		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				///////////////////////////////////
				if (tmp[i][j] == 0) continue;
				processed = false;
				for (t = -1; t <= 1; t++) {
					x = i; y = j + t;
					if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
					if (tmp[x][y] == 0) {
						////////////////////////////
						gray[i][j] = 0;
						processed = true;
						//TRACE("processed!\n");
					}
					if (processed) break;
				}
			} // for i
		} // for j
		tmp.copy(gray);
		//////////////////////////////
		//DrawGrayImage(memDC, image_x, image_y, gray);
		//dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		//TRACE("Erosion iteration %d\n", m);
	} // for m
}


void ReplaceTextureWithDilation(CDC& dc, imatrix& gray, ETF& e, imatrix& tex, int itr) 
{
	int i, j;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	imatrix tmp(image_x, image_y);
	tmp.copy(gray);

	imatrix tmp2(image_x, image_y);
	tmp2.copy(gray);

	imatrix t_y_loc(image_x, image_y);

	t_y_loc.zero();
	gray.white();

	int x, y;

	int s, t, m;
	bool processed;

	double tx, ty, gx, gy;
	int t_x, t_y, x1, y1;

	int tex_x = tex.getRow();
	int tex_y = tex.getCol();

	int freq_x = 20;
	int freq_y = 10;

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			if (tmp[i][j] == 0) { // thin line
				////////////////////////////
				t_x = 0;
				t_x = t_x % tex_x;
				tx = e[i][j].tx;
				ty = e[i][j].ty;
				if ( fabs(tx) >= fabs(ty) ) t_y = i;
				else t_y = j;
				t_y = t_y * freq_y;
				t_y = t_y % tex_y;
				t_y_loc[i][j] = t_y;
				gray[i][j] = tex[t_x][t_y];
				processed = true;
				//TRACE("processed!\n");
			}
		}
	}
    
	for (m = 0; m < itr; m++) {
		//remain = false;
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				///////////////////////////////////
				if (tmp[i][j] == 0) continue;
				processed = false;
				for (t = -1; t <= 1; t++) {
					for (s = -1; s <= 1; s++) {
						x = i + s; y = j + t;
						if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
						if (tmp[x][y] == 0) {
							////////////////////////////
							tmp2[i][j] = 0;
							t_x = m * freq_x;
							t_x = t_x % tex_x;
							gx = -e[x][y].ty;
							gy = e[x][y].tx;
							x1 = round(x + gx);
							y1 = round(y + gy);
							t_y = t_y_loc[x1][y1];
							gray[i][j] = tex[t_x][t_y];
							processed = true;
							//TRACE("processed!\n");
						}
						if (processed) break;
					}
					if (processed) break;
				}
			} // for i
		} // for j
		tmp.copy(tmp2);
		//////////////////////////////
		//DrawGrayImage(memDC, image_x, image_y, gray);
		//dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		//TRACE("Erosion iteration %d\n", m);
	} // for m
}

void OpeningClosingColor(CDC& dc, cimatrix& cmap, imatrix& line, int itr) 
// this corresponds to mean curvature flow
// use line maps as stopping criteria
{
	int i, j;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	nodeRGB *list, *p, *q;
	iRGB col;

	bool found;

    list = NULL;
	//list2 = NULL;

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			found = false;
			for (p = list; p; p=p->next) { // already in the list
				if (p->col.r == cmap[i][j].r && 
					p->col.g == cmap[i][j].g && 
					p->col.b == cmap[i][j].b) {
					found = true; break;
				}
			}
			if (found) p->count++;
			else { // newly insert
				p = (nodeRGB *)malloc(sizeof(nodeRGB));
				p->count = 0;
				p->col = cmap[i][j];
				p->next = list; // add to the list (at the head)
				list = p;
			}
		}
	}

	//TRACE("max_p->col = [%d, %d, %d]\n", max_p->col.r, max_p->col.g, max_p->col.b);
	//TRACE("max_p->count = %d\n", max_p->count);

	cimatrix tmp(image_x, image_y);
	tmp.copy(cmap);

	imatrix done(image_x, image_y);
	done.zero();

	int x, y;

	int s, t, m;
	bool processed;
	//bool remain;
	//double dist;
	//double thres = 100;

	int size = 2; // size of the structuring element

	int mask[5][5] = { 
		{0, 1, 1, 1, 0},
		{1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1},
		{0, 1, 1, 1, 0}
	};

	int max_count;
	nodeRGB *max_p;

	//for (m = 0; m < itr; m++) {
	while (1) {
		/////////////////////////////////////
		// find max node
		max_count = 0;
		max_p = NULL;
		for (p = list; p; p=p->next) {
			//if (p->val > 240 || p->val < 30) continue;
			//if (p->col.r < 10 && p->col.g < 10 && p->col.b < 10) continue; // line
			if (p->count > max_count) {
				max_count = p->count; 
				max_p = p;
			}
		}
		if (max_p == NULL) break; // exhausted
		/////////////////////////////
		//thres *= factor;
		//while (1) {
		for (m = 0; m < itr; m++) {
			//////////////////////////////
			// Opening = Erosion + Dilation	
			//////////////////////////////////////////////
			// Erosion!
			for (j = 0; j < image_y; j++) {
				for (i = 0; i < image_x; i++) {
					///////////////////////////////////
					if (line[i][j] == 0) continue; // stop at lines
					col = tmp[i][j];
					if (col.r != max_p->col.r || col.g != max_p->col.g || col.b != max_p->col.b) {
						//done[i][j] = 1; 
						continue;
					}
					//////////////////////////////////////////////
					//dist = dist3(col.r-max_p->col.r, col.g-max_p->col.g, col.b-max_p->col.b);
					//if ( dist > thres ) continue; // colors are too different; stop
					processed = false;
					for (t = -size; t <= size; t++) {
						for (s = -size; s <= size; s++) {
							if (mask[s+size][t+size] == 0) continue;
							x = i + s; y = j + t;
							if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
							if (tmp[x][y].r != max_p->col.r || 
								tmp[x][y].g != max_p->col.g || 
								tmp[x][y].b != max_p->col.b) {
									////////////////////////////
									cmap[i][j] = tmp[x][y];
									processed = true;
									//remain = true;
							}
							if (processed) break;
						}
						if (processed) break;
					}
				} // for i
			} // for j
			//if (!remain) break;
			tmp.copy(cmap);
			/////////////////////////////////////
			// Dilation!	
			//remain = false;
			for (j = 0; j < image_y; j++) {
				for (i = 0; i < image_x; i++) {
					///////////////////////////////////
					//if (Ix[i][j] == 0.0 && Iy[i][j] == 0.0) continue;
					//if (done[i][j]) continue;
					if (line[i][j] == 0) continue; // stop at lines
					col = tmp[i][j];
					if (col.r == max_p->col.r && col.g == max_p->col.g && col.b == max_p->col.b) {
						//done[i][j] = 1; 
						continue;
					}
					//////////////////////////////////////////////
					//dist = dist3(col.r-max_p->col.r, col.g-max_p->col.g, col.b-max_p->col.b);
					//if ( dist > thres ) continue; // colors are too different; stop
					processed = false;
					for (t = -size; t <= size; t++) {
						for (s = -size; s <= size; s++) {
							if (mask[s+size][t+size] == 0) continue;
							x = i + s; y = j + t;
							if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
							if (tmp[x][y].r == max_p->col.r && 
								tmp[x][y].g == max_p->col.g && 
								tmp[x][y].b == max_p->col.b) {
									////////////////////////////
									cmap[i][j] = max_p->col;
									processed = true;
									//remain = true;
							}
							if (processed) break;
						}
						if (processed) break;
					}
				} // for i
			} // for j
			//if (!remain) break;
			tmp.copy(cmap);
			

			///////////////////////////////////
			// Closing = Dilation + Erosion	
			//////////////////////////////
			// Dilation!	
			//remain = false;
			for (j = 0; j < image_y; j++) {
				for (i = 0; i < image_x; i++) {
					///////////////////////////////////
					//if (Ix[i][j] == 0.0 && Iy[i][j] == 0.0) continue;
					//if (done[i][j]) continue;
					if (line[i][j] == 0) continue; // stop at lines
					col = tmp[i][j];
					if (col.r == max_p->col.r && col.g == max_p->col.g && col.b == max_p->col.b) {
						//done[i][j] = 1; 
						continue;
					}
					//////////////////////////////////////////////
					//dist = dist3(col.r-max_p->col.r, col.g-max_p->col.g, col.b-max_p->col.b);
					//if ( dist > thres ) continue; // colors are too different; stop
					processed = false;
					for (t = -size; t <= size; t++) {
						for (s = -size; s <= size; s++) {
							if (mask[s+size][t+size] == 0) continue;
							x = i + s; y = j + t;
							if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
							if (tmp[x][y].r == max_p->col.r && 
								tmp[x][y].g == max_p->col.g && 
								tmp[x][y].b == max_p->col.b) {
									////////////////////////////
									cmap[i][j] = max_p->col;
									processed = true;
									//remain = true;
							}
							if (processed) break;
						}
						if (processed) break;
					}
				} // for i
			} // for j
			//if (!remain) break;
			tmp.copy(cmap);
			//////////////////////////////////////////////
			// Erosion!
			for (j = 0; j < image_y; j++) {
				for (i = 0; i < image_x; i++) {
					///////////////////////////////////
					if (line[i][j] == 0) continue; // stop at lines
					col = tmp[i][j];
					if (col.r != max_p->col.r || col.g != max_p->col.g || col.b != max_p->col.b) {
						//done[i][j] = 1; 
						continue;
					}
					//////////////////////////////////////////////
					//dist = dist3(col.r-max_p->col.r, col.g-max_p->col.g, col.b-max_p->col.b);
					//if ( dist > thres ) continue; // colors are too different; stop
					processed = false;
					for (t = -size; t <= size; t++) {
						for (s = -size; s <= size; s++) {
							if (mask[s+size][t+size] == 0) continue;
							x = i + s; y = j + t;
							if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
							if (tmp[x][y].r != max_p->col.r || 
								tmp[x][y].g != max_p->col.g || 
								tmp[x][y].b != max_p->col.b) {
									////////////////////////////
									cmap[i][j] = tmp[x][y];
									processed = true;
									//remain = true;
							}
							if (processed) break;
						}
						if (processed) break;
					}
				} // for i
			} // for j
			//if (!remain) break;
			tmp.copy(cmap);
			//////////////////////////////
			DrawColorImage(memDC, image_x, image_y, cmap);
			dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		} // for m
		/////////////////////////////////////
		// delete max node
		for (p=list, q=NULL; p; q=p, p=p->next) {
			//if (p->val > 240 || p->val < 30) continue;
			if (p == max_p) {
				if (q == NULL) { // max_p is the first node
					list = p->next; break;
				}
				q->next = p->next;
				break;
			}
		}
		//free(max_p);
		// move to list2
		//if (list2 == NULL) max_p->next = NULL;
		//else max_p->next = list2;
		///////////////////////////////////////////////////
	} // for m

	/*
	for (p=list2; p; p=q) {
		q = p->next;
		free(p);
	}
	*/
}

void OpeningClosingBlack(CDC& dc, imatrix& gray, int itr) 
// this corresponds to mean curvature flow
{
	int i, j;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	int val;

	imatrix tmp(image_x, image_y);
	tmp.copy(gray);

	int x, y;

	int s, t, m;
	bool processed;
	//bool remain;
	//double dist;
	//double thres = 100;

	int size = 3; // size of the structuring element

	int mask[7][7] = { 
		{0, 0, 1, 1, 1, 0, 0},
		{0, 1, 1, 1, 1, 1, 0},
		{1, 1, 1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1, 1, 1},
		{0, 1, 1, 1, 1, 1, 0},
		{0, 0, 1, 1, 1, 0, 0}
	};

	for (m = 0; m < itr; m++) {
		//////////////////////////////
		// Opening = Erosion + Dilation	
		//////////////////////////////////////////////
		// Erosion!
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				///////////////////////////////////
				val = tmp[i][j];
				if (val != 0) {
					continue;
				}
				processed = false;
				for (t = -size; t <= size; t++) {
					for (s = -size; s <= size; s++) {
						if (mask[s+size][t+size] == 0) continue;
						x = i + s; y = j + t;
						if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
						if (tmp[x][y] != 0) {
							////////////////////////////
							gray[i][j] = tmp[x][y];
							processed = true;
							//remain = true;
						}
						if (processed) break;
					}
					if (processed) break;
				}
			} // for i
		} // for j
		//if (!remain) break;
		tmp.copy(gray);
		/////////////////////////////////////
		// Dilation!	
		//remain = false;
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				val = tmp[i][j];
				if (val == 0) {
					continue;
				}
				//////////////////////////////////////////////
				//dist = dist3(col.r-max_p->col.r, col.g-max_p->col.g, col.b-max_p->col.b);
				//if ( dist > thres ) continue; // colors are too different; stop
				processed = false;
				for (t = -size; t <= size; t++) {
					for (s = -size; s <= size; s++) {
						if (mask[s+size][t+size] == 0) continue;
						x = i + s; y = j + t;
						if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
						if (tmp[x][y] == 0) {
							////////////////////////////
							gray[i][j] = 0;
							processed = true;
							//remain = true;
						}
						if (processed) break;
					}
					if (processed) break;
				}
			} // for i
		} // for j
		//if (!remain) break;
		tmp.copy(gray);
		

		///////////////////////////////////
		// Closing = Dilation + Erosion	
		//////////////////////////////
		// Dilation!	
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				val = tmp[i][j];
				if (val == 0) {
					continue;
				}
				//////////////////////////////////////////////
				processed = false;
				for (t = -size; t <= size; t++) {
					for (s = -size; s <= size; s++) {
						if (mask[s+size][t+size] == 0) continue;
						x = i + s; y = j + t;
						if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
						if (tmp[x][y] == 0) {
								////////////////////////////
								gray[i][j] = 0;
								processed = true;
								//remain = true;
						}
						if (processed) break;
					}
					if (processed) break;
				}
			} // for i
		} // for j
		//if (!remain) break;
		tmp.copy(gray);
		//////////////////////////////////////////////
		// Erosion!
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				///////////////////////////////////
				val = tmp[i][j];
				if (val != 0) {
					continue;
				}
				//////////////////////////////////////////////
				processed = false;
				for (t = -size; t <= size; t++) {
					for (s = -size; s <= size; s++) {
						if (mask[s+size][t+size] == 0) continue;
						x = i + s; y = j + t;
						if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
						if (tmp[x][y] != 0) {
								////////////////////////////
								gray[i][j] = tmp[x][y];
								processed = true;
								//remain = true;
						}
						if (processed) break;
					}
					if (processed) break;
				}
			} // for i
		} // for j
		//if (!remain) break;
		tmp.copy(gray);
		//////////////////////////////
		DrawGrayImage(memDC, image_x, image_y, gray);
		dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		TRACE("OpeningClosing %d\n", m);
		///////////////////////////////////////////////////////////
	} // for m
}

void Closing(imatrix& image, int size, int itr) 
// check neighboring pixels
{
	int i, j, val;

	int image_x = image.getRow();
	int image_y = image.getCol();

	node *list, *p;
	bool found;

	list = NULL;
	
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			found = false;
			for (p = list; p; p=p->next) {
				if (p->val == image[i][j]) {
					found = true; break;
				}
			}
			if (found) p->count++;
			else {
				p = (node *)malloc(sizeof(node));
				p->count = 0;
				p->val = image[i][j];
				p->next = list; // add to the list (at the head)
				list = p;
			}
		}
	}

	int max_count = 0;
	node *max_p = NULL;
	for (p = list; p; p=p->next) {
		if (p->val > 240 || p->val < 30) continue;
		if (p->count > max_count) {
			max_count = p->count; 
			max_p = p;
		}
	}

	//matrix Ix, Iy;
	//Ix.init(image_x, image_y);
	//Iy.init(image_x, image_y);

	imatrix tmp(image_x, image_y);
	tmp.copy(image);

	int x, y, m;

	int thres = 100, s, t;
	bool processed;

	for (m = 0; m < itr; m++) {
		////////////////////////////
		// Dilation
		//////////////////////////////
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				///////////////////////////////////
				//if (Ix[i][j] == 0.0 && Iy[i][j] == 0.0) continue;
				val = tmp[i][j];
				if (val <= 30) continue; // stop at lines
				if ( ABS(val - max_p->val) > thres ) continue; // too different colors
				//TRACE("val = %d\n", val);
				if (val != max_p->val) {
					processed = false;
					for (t = -size; t <= size; t++) {
						for (s = -size; s <= size; s++) { // size: size of the structuring element
							x = i + s; y = j + t;
							if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
							if (tmp[x][y] == max_p->val) {
								image[i][j] = max_p->val;
								processed = true;
							}
							if (processed) break;
						}
						if (processed) break;
					}
				}
			} // for i
		} // for j
		tmp.copy(image);
		/////////////////////////////
		// Erosion
		//////////////////////////////
		///*
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				///////////////////////////////////
				//if (Ix[i][j] == 0.0 && Iy[i][j] == 0.0) continue;
				val = tmp[i][j];
				if (val <= 30) continue; // stop at lines
				if ( ABS(val - max_p->val) > thres ) continue; // too different colors
				//TRACE("val = %d\n", val);
				if (val == max_p->val) {
					processed = false;
					for (t = -1; t <= 1; t++) {
						for (s = -1; s <= 1; s++) {
							x = i + s; y = j + t;
							if (x > image_x-1 || x < 0 || y > image_y-1 || y < 0) continue;
							if (tmp[x][y] != max_p->val) {
								image[i][j] = tmp[x][y];
								processed = true;
							}
							if (processed) break;
						}
						if (processed) break;
					}
				}
			} // for i
		} // for j
		tmp.copy(image);
		//*/
	} // for m

}

void RemoveSmallRegions(CDC& dc, imatrix& line, int size)
{
	int i, j;

	int image_x = line.getRow();
	int image_y = line.getCol();

	imatrix processed(image_x, image_y);
	processed.zero();

	deque<PixeL> pnts;
	deque<PixeL> pnts2;

	PixeL p, r, l, t, b;
	int count, k;

	DrawGrayImage(memDC, image_x, image_y, line);
	dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
	
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			if (line[i][j] == 255) continue;
			if (processed[i][j]) continue;
			pnts.clear();
			pnts2.clear();
			count = 0;
			p.x = i; p.y = j;
			pnts.push_back(p);
			pnts2.push_back(p);
			processed[p.x][p.y] = 255;
			while (1) {
				if ((signed)pnts.size() == 0) break;
				p = pnts[0];
				//TRACE("[before] pnts.size() = %d\n", (signed)pnts.size());
				//TRACE("p.x = %d, p.y = %d\n", p.x, p.y);
				//TRACE("processed[%d][%d] = %d\n", p.x, p.y, processed[p.x][p.y]);
				pnts.pop_front();
				count++;
				//if (count > size) break;
				//TRACE("[after] pnts.size() = %d\n", (signed)pnts.size());
				////////////////////////////////
				r.x = p.x+1; r.y = p.y; 
				if (r.x <= image_x-1)
					if (line[r.x][r.y] == 0 && processed[r.x][r.y] == 0) {
						pnts.push_back(r);
						pnts2.push_back(r);
						processed[r.x][r.y] = 255;
					}
				l.x = p.x-1; l.y = p.y; 
				if (l.x >= 0)
					if (line[l.x][l.y] == 0 && processed[l.x][l.y] == 0) {
						pnts.push_back(l);
						pnts2.push_back(l);
						processed[l.x][l.y] = 255;
					}
				t.x = p.x;   t.y = p.y+1; 
				if (t.y <= image_y-1)
					if (line[t.x][t.y] == 0 && processed[t.x][t.y] == 0) {
						pnts.push_back(t);
						pnts2.push_back(t);
						processed[t.x][t.y] = 255;
					}
				b.x = p.x;   b.y = p.y-1; 
				if (b.y >= 0)
					if (line[b.x][b.y] == 0 && processed[b.x][b.y] == 0) {
						pnts.push_back(b);
						pnts2.push_back(b);
						processed[b.x][b.y] = 255;
					}
			}
			if (count <= size) { // remove this black region
				for (k = 0; k < (signed)pnts2.size(); k++) {
					p = pnts2[k];
					//TRACE("k = %d, p.x = %d, p.y = %d\n", k, p.x, p.y);
					line[p.x][p.y] = 255;
				}
			}
			DrawGrayImage(memDC, image_x, image_y, line);
			dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		}
	}
}


void ColorQuantizationLAB(cimatrix& cmap, int level) 
// LAB color model
// Based on Gooch's method
{
	int x, y;
	GLubyte r1, g1, b1;
	double r, g, b;
	double L, A, B;
	//double t;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r = (double)cmap[x][y].r;
			g = (double)cmap[x][y].g;
			b = (double)cmap[x][y].b;
			//r /= 255.;
			//g /= 255.;
			//b /= 255.;
			//TRACE("r = %f\n", r);
			//TRACE("g = %f\n", g);
			//TRACE("b = %f\n", b);
			//RGB2HSV(r, g, b, h, s, v);
			RGB2LAB(r, g, b, L, A, B);
			///////////////////////////
			L /= 100.; // [0, 1]
			L = ( (int)(L * level) ) / (double)level + 1.0/(double)level/2.0;
			L *= 100.;
			if (L > 100.) L = 100.;
			///////////////////////////////
			////////////////////////////////
			// extrapolate! s and v
			/*
			t = 1.1;
			s = (1-t)*0.0 + t*s;
			v = (1-t)*0.0 + t*v;
			if (s > 1.0) s = 1.0;
			if (v > 1.0) v = 1.0;
			*/
			///////////////////////////////
			LAB2RGB(L, A, B, r, g, b);
			//HSV2RGB(h, s, v, r, g, b);
			r1 = (GLubyte)r;
			g1 = (GLubyte)g;
			b1 = (GLubyte)b;
			cmap[x][y].r = r1;
			cmap[x][y].g = g1;
			cmap[x][y].b = b1;
		}
	}
}

void ColorQuantizationLAB2(cimatrix& cmap, int level1, int level2, int level3) 
// LAB color model
// L: [0, 100], a: [-86.184636, 98.254219], b: [-107.863681, 94.482485] 
{
	int x, y;
	GLubyte r1, g1, b1;
	double r, g, b;
	double L, A, B;
	//double t;

	//int level1 = 5;
	//int level2 = 10;
	//int level3 = 10;

	double A_min = -86.184636;
	double A_max = 98.254219;
	double B_min = -107.863681;
	double B_max = 94.482485;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r = (double)cmap[x][y].r;
			g = (double)cmap[x][y].g;
			b = (double)cmap[x][y].b;
			//r /= 255.;
			//g /= 255.;
			//b /= 255.;
			//TRACE("r = %f\n", r);
			//TRACE("g = %f\n", g);
			//TRACE("b = %f\n", b);
			//RGB2HSV(r, g, b, h, s, v);
			RGB2LAB(r, g, b, L, A, B);
			///////////////////////////
			L /= 100.; // [0, 1]
			if (L > 0.99) L = 1.0;
			else if (L < 0.01) L = 0.0;
			else L = ( (int)(L * level1) ) / (double)level1 + 1.0/(double)level1/2.0;
			//t = 1.1;
			//L = (1-t)*0.0 + t*L;
			if (L > 1.0) L = 1.0;
			if (L < 0.0) L = 0.0;
			L *= 100.;

			///////////////////////////////
			A += fabs(A_min); // make it positive
			A /= (A_max - A_min); // [0, 1]
			if (A == 1.0) A = 1.0;
			else if (A == 0.0) A = 0.0;
			else A = ( (int)(A * level2) ) / (double)level2 + 1.0/(double)level2/2.0;
			if (A > 1.0) A = 1.0;
			if (A < 0.0) A = 0.0;
			A *= (A_max - A_min);
			A -= fabs(A_min);
			////////////////////////////////
			B += fabs(B_min); // make it positive
			B /= (B_max - B_min); // [0, 1]
			if (B == 1.0) B = 1.0;
			else if (B == 0.0) B = 0.0;
			else B = ( (int)(B * level3) ) / (double)level3 + 1.0/(double)level3/2.0;
			if (B > 1.0) B = 1.0;
			if (B < 0.0) B = 0.0;
			B *= (B_max - B_min);
			B -= fabs(B_min);
			////////////////////////////////
			// extrapolate! s and v
			/*
			t = 1.1;
			s = (1-t)*0.0 + t*s;
			v = (1-t)*0.0 + t*v;
			if (s > 1.0) s = 1.0;
			if (v > 1.0) v = 1.0;
			*/
			///////////////////////////////
			LAB2RGB(L, A, B, r, g, b);
			//HSV2RGB(h, s, v, r, g, b);
			r1 = (GLubyte)r;
			g1 = (GLubyte)g;
			b1 = (GLubyte)b;
			cmap[x][y].r = r1;
			cmap[x][y].g = g1;
			cmap[x][y].b = b1;
		}
	}
}


void DrawCartoonImage(CDC& dc, int image_x, int image_y, cimatrix& image, imatrix& gray) 
// Merge color background and black outlines
{
	int x, y;
	GLubyte r, g, b;

	//for (y = IMAGE_Y - 1; y >= 0; y--) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			if (gray[x][y] == 0) {
				r = g = b = 0;
			}
			else {
				r = (GLubyte)image[x][y].r;
				g = (GLubyte)image[x][y].g;
				b = (GLubyte)image[x][y].b;
			}
			/// Set Pixel in MemDC
			dc.SetPixelV(x, (IMAGE_Y-1)-y, RGB(r, g, b));
		}
	}
}

void ConstructCartoonImage(cimatrix& image, imatrix& gray, cimatrix& merged) 
// Merge color background and black outlines
{
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	//for (y = IMAGE_Y - 1; y >= 0; y--) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			if (gray[x][y] == 0) {
				merged[x][y].r = 0;
				merged[x][y].g = 0;
				merged[x][y].b = 0;
			}
			else {
				merged[x][y].r = image[x][y].r;
				merged[x][y].g = image[x][y].g;
				merged[x][y].b = image[x][y].b;
			}
		}
	}
}

void ConstructCartoonImage2(cimatrix& image, imatrix& gray, cimatrix& merged, double threshold) 
// Merge color background and black outlines
/// with threshold [0, 1]
{
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	//for (y = IMAGE_Y - 1; y >= 0; y--) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			if (gray[x][y] < threshold * 255) {
				merged[x][y].r = (GLubyte)gray[x][y];
				merged[x][y].g = (GLubyte)gray[x][y];
				merged[x][y].b = (GLubyte)gray[x][y];
			}
			else {
				merged[x][y].r = image[x][y].r;
				merged[x][y].g = image[x][y].g;
				merged[x][y].b = image[x][y].b;
			}
		}
	}
}

void ConstructCartoonImage3(cimatrix& image, imatrix& gray, cimatrix& merged) 
// Merge color background and black outlines
// with smooth blending!
{
	int x, y;
	double t;

	int image_x = image.getRow();
	int image_y = image.getCol();

	//for (y = IMAGE_Y - 1; y >= 0; y--) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			t = gray[x][y] / 255.0;
			merged[x][y].r = (GLubyte)( (1-t) * gray[x][y] + t * image[x][y].r );
			merged[x][y].g = (GLubyte)( (1-t) * gray[x][y] + t * image[x][y].g );
			merged[x][y].b = (GLubyte)( (1-t) * gray[x][y] + t * image[x][y].b );
		}
	}
}

void DrawGrayImageMask(CDC& dc, int image_x, int image_y, int size, PixeL seed, imatrix& image) 
{
	GLubyte r2;
	int	s, r, half;
	int	i, j;

	half = size/2;

	for (s = 0; s < size; s++) {
		for (r = 0; r < size; r++) {
			/////////////////////////////////////////////////////////
			// circular kernel
			if ( dist2(s, r, half, half) > half-1 ) continue; 
			//////////////////////////////////////////////////
			i = seed.x-half+s;
			j = seed.y-half+r;
			if (i <= 0 || i >= image_x-1 || j <= 0 || j >= image_y-1)
				continue;
			r2 = (GLubyte)image[i][j];
			/// Set Pixel in MemDC
            dc.SetPixelV(i, (IMAGE_Y-1)-j, RGB(r2, r2, r2));
			/////////////////////////////////
		}
	}
}

void DrawEdgeStrokes(CDC& dc, int image_x, int image_y, imatrix& image, int width) 
{
	int x, y;
	double d_i, d_j, tx, ty, nx, ny;
	int	st_x, st_y, end_x, end_y;
	int int_i, int_j;
	double t;
	CPen pen, *pOldPen;
	imatrix visit;
	
	visit.init(image_x, image_y);
	visit.zero();
	ClearMemDC(&dc); // clear the canvas white

	//gx = G_x[i][j];
	//gy = G_y[i][j];
	//g = G_mag[i][j];
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			if (image[x][y] == 0 && visit[x][y] == 0) {
				visit[x][y] = 1;
				//////////////////////////////////////////////////////////////
				// draw strokes
				////////////////////////////////////////////////
				// One half
				pen.CreatePen(PS_SOLID, width, RGB(0, 0, 0));
				pOldPen = (CPen *)dc.SelectObject(&pen);
				d_i = (double)x; d_j = (double)y; 
				tx = -G_y[x][y];
				ty = G_x[x][y];
				st_x = x;
				st_y = y;
				while(1) {
					if (tx == 0.0 && ty == 0.0) break;
					nx = tx / sqrt(tx*tx+ty*ty);
					ny = ty / sqrt(tx*tx+ty*ty);
					if (ABS(nx) > ABS(ny)) t = 1/ABS(nx);
					else t = 1/ABS(ny);
					d_i += nx*t;
					d_j += ny*t;
					//(-ny, nx)
					if (d_i < 0 || d_i > IMAGE_X-1 || d_j < 0 || d_j > IMAGE_Y-1) break;
					int_i = round(d_i);
					int_j = round(d_j);
					////////////////////////////////////
					visit[int_i][int_j] = 1;
					if (image[int_i][int_j] != 0) break;
					///////////////////////////////////
					//TRACE("i = %.1f, j = %.1f\n", i, j);
					end_x = int_i; end_y = int_j;
					dc.MoveTo(st_x, IMAGE_Y-1-st_y);
					dc.LineTo(end_x, IMAGE_Y-1-end_y);
					st_x = end_x; st_y = end_y;
					////////////////////////////
					//sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, IMAGE_Y-1-int_j));
					//TRACE("sum = %d\n", sum);
					//count++;
					tx = -G_y[x][y];
					ty = G_x[x][y];
				}
				dc.SelectObject(pOldPen);
				pen.DeleteObject(); 
				////////////////////////////////////////
				// Other half
				pen.CreatePen(PS_SOLID, width, RGB(0, 0, 0));
				pOldPen = (CPen *)dc.SelectObject(&pen);
				d_i = (double)x; d_j = (double)y; 
				tx = -G_y[x][y];
				ty = G_x[x][y];
				st_x = x;
				st_y = y;
				while(1) {
					if (tx == 0.0 && ty == 0.0) break;
					nx = tx / sqrt(tx*tx+ty*ty);
					ny = ty / sqrt(tx*tx+ty*ty);
					if (ABS(nx) > ABS(ny)) t = 1/ABS(nx);
					else t = 1/ABS(ny);
					d_i -= nx*t;
					d_j -= ny*t;
					//(-ny, nx)
					if (d_i < 0 || d_i > IMAGE_X-1 || d_j < 0 || d_j > IMAGE_Y-1) break;
					int_i = round(d_i);
					int_j = round(d_j);
					////////////////////////////////////
					visit[int_i][int_j] = 1;
					if (image[int_i][int_j] != 0) break;
					///////////////////////////////////
					//TRACE("i = %.1f, j = %.1f\n", i, j);
					end_x = int_i; end_y = int_j;
					dc.MoveTo(st_x, IMAGE_Y-1-st_y);
					dc.LineTo(end_x, IMAGE_Y-1-end_y);
					st_x = end_x; st_y = end_y;
					////////////////////////////
					//sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, IMAGE_Y-1-int_j));
					//TRACE("sum = %d\n", sum);
					//count++;
					tx = -G_y[x][y];
					ty = G_x[x][y];
				}
				dc.SelectObject(pOldPen);
				pen.DeleteObject(); 
			}
			/// Set Pixel in MemDC
			//dc.SetPixelV(x, (IMAGE_Y-1)-y, RGB(r, r, r));
		}
	}
}

void DrawEdgeStrokes2(CDC& dc, int image_x, int image_y, imatrix& image, int width) 
// try to maximize the connectivity
{
	int x, y;
	//double d_i, d_j, tx, ty, nx, ny;
	//int	st_x, st_y, end_x, end_y;
	int int_i, int_j;
	CPen pen, *pOldPen;
	imatrix visit;
	int i, j, xx, yy, found;
	
	visit.init(image_x, image_y); // is this pixel already visited?
	visit.zero(); // initially no pixel is visited
	ClearMemDC(&dc); // clear the canvas white

	//gx = G_x[i][j];
	//gy = G_y[i][j];
	//g = G_mag[i][j];
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			if (image[x][y] == 0 && visit[x][y] == 0) {
				visit[x][y] = 1;
				//////////////////////////////////////////////////////////////
				// draw strokes
				////////////////////////////////////////////////
				// One half
				pen.CreatePen(PS_SOLID, width, RGB(0, 0, 0));
				pOldPen = (CPen *)dc.SelectObject(&pen);
				int_i = x;
				int_j = y;
				dc.MoveTo(int_i, IMAGE_Y-1-int_j);
				while(1) {
					found = 0;
					/////////////////////////////////
					// search the 8-neighbor to find the connecting edge pixel
					for (i = -1; i <= 1; i++) {
						for (j = -1; j <= 1; j++) {
							xx = int_i + i;
							yy = int_j + j;
							if (xx < 0 || xx > IMAGE_X-1 || yy < 0 || yy > IMAGE_Y-1) continue;
							if (image[xx][yy] == 0 && visit[xx][yy] == 0) {
								found = 1;
								break; // found
							}
						}
						if (found) break;
					}
					if (!found) break; // dead end; no need to search any more
					///////////////////////////////
					int_i = xx;
					int_j = yy;
					visit[int_i][int_j] = 1;
					///////////////////////////////////
					dc.LineTo(int_i, IMAGE_Y-1-int_j);
					////////////////////////////
				}
				dc.SelectObject(pOldPen);
				pen.DeleteObject(); 
				////////////////////////////////////////
				// Other half
				pen.CreatePen(PS_SOLID, width, RGB(0, 0, 0));
				pOldPen = (CPen *)dc.SelectObject(&pen);
				int_i = x;
				int_j = y;
				dc.MoveTo(int_i, IMAGE_Y-1-int_j);
				while(1) {
					found = 0;
					for (i = -1; i <= 1; i++) {
						for (j = -1; j <= 1; j++) {
							xx = int_i + i;
							yy = int_j + j;
							if (xx < 0 || xx > IMAGE_X-1 || yy < 0 || yy > IMAGE_Y-1) continue;
							if (image[xx][yy] == 0 && visit[xx][yy] == 0) {
								found = 1;
								break; // found
							}
						}
						if (found) break;
					}
					if (!found) break; // dead end
					///////////////////////////////
					int_i = xx;
					int_j = yy;
					visit[int_i][int_j] = 1;
					///////////////////////////////////
					dc.LineTo(int_i, IMAGE_Y-1-int_j);
					////////////////////////////
				}
				dc.SelectObject(pOldPen);
				pen.DeleteObject(); 
				
			}
			/// Set Pixel in MemDC
			//dc.SetPixelV(x, (IMAGE_Y-1)-y, RGB(r, r, r));
		}
	}
}

//deque<PntsDeque> pnt_array;
/*
struct double2D {
	double x, y;
};
*/

//typedef deque<double2D> PntsDeque;
//deque<PntsDeque*> big_array;
deque<deque<double2D>*> big_array; // a one-dimensional array of strokes (each stroke as point set)
//deque<double2D> *pnts_array;

void StoreEdgeStrokes(int image_x, int image_y, imatrix& image) 
{
	int x, y;
	double d_i, d_j, tx, ty, nx, ny;
	int	st_x, st_y, end_x, end_y;
	int int_i, int_j;
	double t;
	imatrix visit;
	double2D p;
	//PntsDeque* dq;
	deque<double2D>* dq;
	
	visit.init(image_x, image_y);
	visit.zero();
	//ClearMemDC(&dc); // clear the canvas white

	//gx = G_x[i][j];
	//gy = G_y[i][j];
	//g = G_mag[i][j];
	//deque<PntsDeque>::iterator it;
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//TRACE("i'm here!\n");
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			if (image[x][y] == 0 && visit[x][y] == 0) {
				//dq = new PntsDeque;
				dq = new deque<double2D>;
				big_array.insert(big_array.end(), dq);
				//TRACE("big_array.size() = %d\n", big_array.size());
				visit[x][y] = 1;
				///////////////////////////////////
				p.x = (double)x; 
				p.y = (double)y; 
				//pnt_array.insert(pnt_array.end(), p);
				//it = pnts_array.end();
				//it->insert(it->end(), p);
				dq->insert(dq->end(), p);
				//////////////////////////////////////////////////////////////
				// draw strokes
				////////////////////////////////////////////////
				// One half
				d_i = (double)x; d_j = (double)y; 
				tx = -G_y[x][y];
				ty = G_x[x][y];
				st_x = x;
				st_y = y;
				while(1) {
					if (tx == 0.0 && ty == 0.0) break;
					nx = tx / sqrt(tx*tx+ty*ty);
					ny = ty / sqrt(tx*tx+ty*ty);
					if (ABS(nx) > ABS(ny)) t = 1/ABS(nx);
					else t = 1/ABS(ny);
					d_i += nx*t;
					d_j += ny*t;
					//(-ny, nx)
					if (d_i < 0 || d_i > IMAGE_X-1 || d_j < 0 || d_j > IMAGE_Y-1) break;
					int_i = round(d_i);
					int_j = round(d_j);
					////////////////////////////////////
					visit[int_i][int_j] = 1;
					if (image[int_i][int_j] != 0) break;
					///////////////////////////////////
					//TRACE("i = %.1f, j = %.1f\n", i, j);
					end_x = int_i; end_y = int_j;
					st_x = end_x; st_y = end_y;
					///////////////////////////////////
					p.x = (double)int_i; 
					p.y = (double)int_j; 
					//pnt_array.insert(pnt_array.end(), p);
					//it = pnts_array.end();
					//it->insert(it->end(), p);
					dq->insert(dq->end(), p);
					////////////////////////////
					//sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, IMAGE_Y-1-int_j));
					//TRACE("sum = %d\n", sum);
					//count++;
					tx = -G_y[x][y];
					ty = G_x[x][y];
				}
				////////////////////////////////////////
				// Other half
				d_i = (double)x; d_j = (double)y; 
				tx = -G_y[x][y];
				ty = G_x[x][y];
				st_x = x;
				st_y = y;
				while(1) {
					if (tx == 0.0 && ty == 0.0) break;
					nx = tx / sqrt(tx*tx+ty*ty);
					ny = ty / sqrt(tx*tx+ty*ty);
					if (ABS(nx) > ABS(ny)) t = 1/ABS(nx);
					else t = 1/ABS(ny);
					d_i -= nx*t;
					d_j -= ny*t;
					//(-ny, nx)
					if (d_i < 0 || d_i > IMAGE_X-1 || d_j < 0 || d_j > IMAGE_Y-1) break;
					int_i = round(d_i);
					int_j = round(d_j);
					////////////////////////////////////
					visit[int_i][int_j] = 1;
					if (image[int_i][int_j] != 0) break;
					///////////////////////////////////
					//TRACE("i = %.1f, j = %.1f\n", i, j);
					end_x = int_i; end_y = int_j;
					st_x = end_x; st_y = end_y;
					///////////////////////////////////
					p.x = (double)int_i; 
					p.y = (double)int_j; 
					//pnt_array.insert(pnt_array.begin(), p);
					//it = pnts_array.end();
					//it->insert(it->begin(), p);
					dq->insert(dq->begin(), p);
					////////////////////////////
					//sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, IMAGE_Y-1-int_j));
					//TRACE("sum = %d\n", sum);
					//count++;
					tx = -G_y[x][y];
					ty = G_x[x][y];
				} // while
			} // if
			//big_array.insert(big_array.end(), dq);
			/// Set Pixel in MemDC
			//dc.SetPixelV(x, (IMAGE_Y-1)-y, RGB(r, r, r));
		}
	}

}


void StoreEdgeStrokes2(int image_x, int image_y, imatrix& image) 
// try to maximize the length of connected points
{
	int x, y;
	//double d_i, d_j, tx, ty, nx, ny;
	//int	st_x, st_y, end_x, end_y;
	int int_i, int_j;
	imatrix visit;
	int i, j, xx, yy, found;
	double2D p;
	//PntsDeque* dq;
	deque<double2D>* dq;
	
	visit.init(image_x, image_y);
	visit.zero();
	//ClearMemDC(&dc); // clear the canvas white

	//gx = G_x[i][j];
	//gy = G_y[i][j];
	//g = G_mag[i][j];
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			if (image[x][y] == 0 && visit[x][y] == 0) {
				visit[x][y] = 1;
				dq = new deque<double2D>;
				big_array.insert(big_array.end(), dq);
				p.x = (double)x; 
				p.y = (double)y; 
				dq->insert(dq->end(), p);
				//////////////////////////////////////////////////////////////
				// draw strokes
				////////////////////////////////////////////////
				// One half
				//pen.CreatePen(PS_SOLID, 3, RGB(255, 255, 0));
				//pOldPen = (CPen *)dc.SelectObject(&pen);
				int_i = x;
				int_j = y;
				//dc.MoveTo(int_i, IMAGE_Y-1-int_j);
				while(1) {
					found = 0;
					for (i = -1; i <= 1; i++) {
						for (j = -1; j <= 1; j++) {
							xx = int_i + i;
							yy = int_j + j;
							if (xx < 0 || xx > IMAGE_X-1 || yy < 0 || yy > IMAGE_Y-1) continue;
							if (image[xx][yy] == 0 && visit[xx][yy] == 0) {
								found = 1;
								break; // found
							}
						}
						if (found) break;
					}
					///////////////////////////////
					if (!found) break; // dead end
					///////////////////////////////
					int_i = xx;
					int_j = yy;
					visit[int_i][int_j] = 1;
					///////////////////////////////////
					p.x = (double)int_i; 
					p.y = (double)int_j; 
					dq->insert(dq->end(), p);
					///////////////////////////////////
					//dc.LineTo(int_i, IMAGE_Y-1-int_j);
					////////////////////////////
				}
				//dc.SelectObject(pOldPen);
				//pen.DeleteObject(); 
				////////////////////////////////////////
				// Other half
				//pen.CreatePen(PS_SOLID, 3, RGB(255, 0, 0));
				//pOldPen = (CPen *)dc.SelectObject(&pen);
				int_i = x;
				int_j = y;
				//dc.MoveTo(int_i, IMAGE_Y-1-int_j);
				while(1) {
					found = 0;
					for (i = -1; i <= 1; i++) {
						for (j = -1; j <= 1; j++) {
							xx = int_i + i;
							yy = int_j + j;
							if (xx < 0 || xx > IMAGE_X-1 || yy < 0 || yy > IMAGE_Y-1) continue;
							if (image[xx][yy] == 0 && visit[xx][yy] == 0) {
								found = 1;
								break; // found
							}
						}
						if (found) break;
					}
					if (!found) break; // dead end
					///////////////////////////////
					int_i = xx;
					int_j = yy;
					visit[int_i][int_j] = 1;
					///////////////////////////////////
					p.x = (double)int_i; 
					p.y = (double)int_j; 
					dq->insert(dq->begin(), p);
					///////////////////////////////////
					//dc.LineTo(int_i, IMAGE_Y-1-int_j);
					////////////////////////////
				}
				//dc.SelectObject(pOldPen);
				//pen.DeleteObject(); 
				
			}
			/// Set Pixel in MemDC
			//dc.SetPixelV(x, (IMAGE_Y-1)-y, RGB(r, r, r));
		}
	}
}

//////////////////////////
MyList<MRBspline> bstrokes; // a list of B-spline strokes!
//MyList<MRBspline> bstrokes; // a list of Hermite spline strokes!

qmatrix<MRBspline*> b_map; 


void adaptive_hyster_visit(int x0, int y0, imatrix& image2, double lo_thres, double lo_thres2)
{
	int	i, j, x, y;
	int	done_flag;
	double g;

	//image2[x0][y0] = 0; // now marked as edge
	done_flag = 0;
	for (j = -1; j <= 1; j++) {
		for (i = -1; i <= 1; i++) {
			x = x0 + i;
			y = y0 + j;
			if (x <= 0 || x >= IMAGE_X-1 || y <= 0 || y >= IMAGE_Y-1)
				continue;
			if (image2[x][y] == 0) continue; // already marked as edge
			if (!thin_edge[x][y]) continue; // not a thin edge
			g = G_mag[x][y];
			if (pixel_mark[x][y] <= 102) { // foreground pixels
				if (g > lo_thres) {
					image2[x][y] = 0; // newly marked as edge pixel
					adaptive_hyster_visit(x, y, image2, lo_thres, lo_thres2);
					//done_flag = 1; // found an extending edge pixel
					//break;
				}
			}
			else {
				if (g > lo_thres2) {
					image2[x][y] = 0; // newly marked as edge pixel
					adaptive_hyster_visit(x, y, image2, lo_thres, lo_thres2);
					//done_flag = 1; // found an extending edge pixel
					//break;
				}
			}
		}
		//if (done_flag) break; // found the right neighboring edge pixel. Get out!
	}
}

double GlobalAdaptiveCanny(imatrix& image, imatrix& image2)
{
	int	i, j, k, i1, i2, j1, j2;
	double MAX_GRADIENT = -1;
	
	int gau_w, gau_w2;
	double hi_thres, hi_thres2, lo_thres, lo_thres2;

	gau_w = MakeGaussMask(1.0, gau);
	gau_w2 = MakeGaussMask(0.4, gau2);
	hi_thres = 0.1;
	lo_thres = 0.1;
	hi_thres2 = 0.17;
	lo_thres2 = 0.07;

	//imatrix marked(IMAGE_X, IMAGE_Y);
	//image2.init(IMAGE_X, IMAGE_Y);
	
	//imatrix tmp(IMAGE_X, IMAGE_Y);
	double	x, y;
	
	//////////////////////////////////////////////
	// Gaussian smoothing (adaptive)
	for (j = 0; j < IMAGE_Y; j++) {
		for (i = 0; i < IMAGE_X; i++) {
			if (pixel_mark[i][j]) { // selected pixels
				x = gau[0] * image[i][j];	
				y = gau[0] * image[i][j];
				//TRACE("x = %f, y = %f\n", x, y);
				for (k = 1; k < gau_w; k++) {
					i1 = (i+k)%IMAGE_X;
					i2 = (i-k+IMAGE_X)%IMAGE_X;
					x += gau[k] * image[i1][j] + gau[k] * image[i2][j];
					j1 = (j+k)%IMAGE_Y;
					j2 = (j-k+IMAGE_Y)%IMAGE_Y;
					y += gau[k] * image[i][j1] + gau[k] * image[i][j2];
				}
				tmp_x[i][j] = x;
				tmp_y[i][j] = y;
			}
			else {
				x = gau2[0] * image[i][j];	
				y = gau2[0] * image[i][j];
				//TRACE("x = %f, y = %f\n", x, y);
				for (k = 1; k < gau_w2; k++) {
					i1 = (i+k)%IMAGE_X;
					i2 = (i-k+IMAGE_X)%IMAGE_X;
					x += gau2[k] * image[i1][j] + gau2[k] * image[i2][j];
					j1 = (j+k)%IMAGE_Y;
					j2 = (j-k+IMAGE_Y)%IMAGE_Y;
					y += gau2[k] * image[i][j1] + gau2[k] * image[i][j2];
				}
				tmp_x[i][j] = x;
				tmp_y[i][j] = y;
			}
			//TRACE("x = %f, y = %f\n", x, y);
			if (x > 255) x = 255;
			if (y > 255) y = 255;
			//image[i][j] = (int)x;
			//image[i][j] = (int)y;
		}
	}

	for (j = 1; j < IMAGE_Y - 1; j++) {
		for (i = 1; i < IMAGE_X - 1; i++) {
			G_x[i][j] = (tmp_x[i+1][j-1] + 2*tmp_x[i+1][j] + tmp_x[i+1][j+1] 
				- tmp_x[i-1][j-1] - 2*tmp_x[i-1][j] - tmp_x[i-1][j+1]);
			G_y[i][j] = (tmp_y[i-1][j+1] + 2*tmp_y[i][j+1] + tmp_y[i+1][j+1]
				- tmp_y[i-1][j-1] - 2*tmp_y[i][j-1] - tmp_y[i+1][j-1]);
			//G_mag[i][j] = sqrt(G_x[i][j] * G_x[i][j] + G_y[i][j] * G_y[i][j]);
			G_mag[i][j] = norm2(G_x[i][j], G_y[i][j]);
			
			/*
			if (G_mag[i][j] > MAX_GRADIENT) {
				MAX_GRADIENT = G_mag[i][j];
				//TRACE("p[i][j].gx = %01.f\n", p[i][j].gx);
				//TRACE("p[i][j].gy = %01.f\n", p[i][j].gy);
				//TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);
			}
			*/
		}
	}

	// Normalize each gradient value & init marked image
	//TRACE("MAX_GRADIENT = %f\n", MAX_GRADIENT);
	//TRACE("max_grad2 = %f\n", max_grad2);
	//TRACE("hi_thres = %f\n", hi_thres);
	for (j = 0; j < IMAGE_Y; j++) {
		for (i = 0; i < IMAGE_X; i++) {
			if (i == 0 || i == IMAGE_X-1 || j == 0 || j == IMAGE_Y-1) {
				image2[i][j] = 255;
				thin_edge[i][j] = 0; // init thin edge list
				continue;
			}
			G_mag[i][j] = (G_mag[i][j] / max_grad2); // G_mag between [0, 1]
			image2[i][j] = (int)(G_mag[i][j] * 255);

			//marked[i][j] = 0; // init marked image for hysteresis
			thin_edge[i][j] = 0; // init thin edge list
		}
	}

	double	gx, gy;
	double	g, g1, g2, g3, g4;
	double	t; // interpolation parameter
	
	//////////////////////////////////////////
	// Nonmaxima suppression
	for (j = 1; j < IMAGE_Y-1; j++) {
		for (i = 1; i < IMAGE_X-1; i++) {
			gx = G_x[i][j];
			gy = G_y[i][j];
			g = G_mag[i][j];
			//TRACE("gx = %f, gy = %f, g = %f\n", gx, gy, g);
			//if (gx < 0.01 && gy < 0.01) continue; // not an edge
			if (fabs(gx) >= fabs(gy)) { // close to horizontal (note: gy can be 0)
				t = fabs(gy) / fabs(gx); 
				//TRACE("t = %f\n", t);
				g1 = G_mag[i+1][j]; // right
				g2 = G_mag[i-1][j]; // left
				//TRACE("g1 = %f, g2 = %f\n", g1, g2);
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // right up
					g4 = G_mag[i-1][j-1]; // left down
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i+1][j-1]; // right down
					g4 = G_mag[i-1][j+1]; // left up
				}
			}
			else { // close to vertical (note: gx can be 0)
				t = fabs(gx) / fabs(gy);
				g1 = G_mag[i][j+1]; // up
				g2 = G_mag[i][j-1]; // down
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // up right
					g4 = G_mag[i-1][j-1]; // down left
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i-1][j+1]; // up left
					g4 = G_mag[i+1][j-1]; // down right
				}
			}
			if ( g > ((1-t)*g1 + t*g3) && g > ((1-t)*g2 + t*g4) ) {
                thin_edge[i][j] = 1; // it's a thin edge
				if (pixel_mark[i][j]) { // selected pixels
					//TRACE("g = %f\n", g);
					if (g > hi_thres) {
					//if (g > 0.0) {
						//rr = Dbuffer[(j * IMAGE_X + i) * 3 + 0];
						//gg = Dbuffer[(j * IMAGE_X + i) * 3 + 1];
						//bb = Dbuffer[(j * IMAGE_X + i) * 3 + 2];
						/*
						if ( prob2((double)rr, (double)gg, (double)bb, g, off_edge, factor1) )
							off_mark[i][j] = 255; // off edge
						*/
						image2[i][j] = 0; // thin edge above hi_thres
					}
					else 
						image2[i][j] = 255; // thin edge below hi_thres
				}
				else {
					if (g > hi_thres2) {
						image2[i][j] = 0; // thin edge above hi_thres
					}
					else 
						image2[i][j] = 255; // thin edge below hi_thres
				}
			}
			else { // non-maximum
				image2[i][j] = 255;
			}
			
		}
	}

	//////////////////////////////////////////
	// Hysteresis thresholding
	/*
	for (j = 1; j < IMAGE_Y-1; j++) {
		for (i = 1; i < IMAGE_X-1; i++) {
			if (image2[i][j] == 0) { // computed thinned edges above hi_thres
				adaptive_hyster_visit(i, j, image2, lo_thres, lo_thres2); // visit neighboring pixels
			}
		}
	}
	*/
	//////////////////////////////////////////
	// removing off_edges
	/*
	for (j = 1; j < IMAGE_Y-1; j++) {
		for (i = 1; i < IMAGE_X-1; i++) {
			if (off_mark[i][j] == 255) { // marked off edge
				off_edge_visit(IMAGE_X, IMAGE_Y, i, j); // visit neighboring pixels
			}
		}
	}
	*/
	
	return MAX_GRADIENT;
}

double GlobalAdaptiveCanny2(imatrix& image, imatrix& image2)
{
	int	i, j, k, i1, i2, j1, j2;
	double MAX_GRADIENT = -1;
	
	int gau_w, gau_w2;
	double hi_thres, hi_thres2, lo_thres, lo_thres2;

	gau_w = MakeGaussMask(0.8, gau);
	gau_w2 = MakeGaussMask(0.4, gau2);
	hi_thres = 0.21;
	lo_thres = 0.16;
	hi_thres2 = 0.17;
	lo_thres2 = 0.07;

	//imatrix marked(IMAGE_X, IMAGE_Y);
	//image2.init(IMAGE_X, IMAGE_Y);
	
	//imatrix tmp(IMAGE_X, IMAGE_Y);
	double	x, y;
	
	//////////////////////////////////////////////
	// Gaussian smoothing (adaptive)
	for (j = 0; j < IMAGE_Y; j++) {
		for (i = 0; i < IMAGE_X; i++) {
			if (pixel_mark[i][j] <= 102) { // foreground pixels
				x = gau[0] * image[i][j];	
				y = gau[0] * image[i][j];
				//TRACE("x = %f, y = %f\n", x, y);
				for (k = 1; k < gau_w; k++) {
					i1 = (i+k)%IMAGE_X;
					i2 = (i-k+IMAGE_X)%IMAGE_X;
					x += gau[k] * image[i1][j] + gau[k] * image[i2][j];
					j1 = (j+k)%IMAGE_Y;
					j2 = (j-k+IMAGE_Y)%IMAGE_Y;
					y += gau[k] * image[i][j1] + gau[k] * image[i][j2];
				}
				tmp_x[i][j] = x;
				tmp_y[i][j] = y;
			}
			else { // other area
				x = gau2[0] * image[i][j];	
				y = gau2[0] * image[i][j];
				//TRACE("x = %f, y = %f\n", x, y);
				for (k = 1; k < gau_w2; k++) {
					i1 = (i+k)%IMAGE_X;
					i2 = (i-k+IMAGE_X)%IMAGE_X;
					x += gau2[k] * image[i1][j] + gau2[k] * image[i2][j];
					j1 = (j+k)%IMAGE_Y;
					j2 = (j-k+IMAGE_Y)%IMAGE_Y;
					y += gau2[k] * image[i][j1] + gau2[k] * image[i][j2];
				}
				tmp_x[i][j] = x;
				tmp_y[i][j] = y;
			}
			//TRACE("x = %f, y = %f\n", x, y);
			if (x > 255) x = 255;
			if (y > 255) y = 255;
			//image[i][j] = (int)x;
			//image[i][j] = (int)y;
		}
	}

	for (j = 1; j < IMAGE_Y - 1; j++) {
		for (i = 1; i < IMAGE_X - 1; i++) {
			G_x[i][j] = (tmp_x[i+1][j-1] + 2*tmp_x[i+1][j] + tmp_x[i+1][j+1] 
				- tmp_x[i-1][j-1] - 2*tmp_x[i-1][j] - tmp_x[i-1][j+1]);
			G_y[i][j] = (tmp_y[i-1][j+1] + 2*tmp_y[i][j+1] + tmp_y[i+1][j+1]
				- tmp_y[i-1][j-1] - 2*tmp_y[i][j-1] - tmp_y[i+1][j-1]);
			//G_mag[i][j] = sqrt(G_x[i][j] * G_x[i][j] + G_y[i][j] * G_y[i][j]);
			G_mag[i][j] = norm2(G_x[i][j], G_y[i][j]);
			
			/*
			if (G_mag[i][j] > MAX_GRADIENT) {
				MAX_GRADIENT = G_mag[i][j];
				//TRACE("p[i][j].gx = %01.f\n", p[i][j].gx);
				//TRACE("p[i][j].gy = %01.f\n", p[i][j].gy);
				//TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);
			}
			*/
		}
	}

	// Normalize each gradient value & init marked image
	//TRACE("MAX_GRADIENT = %f\n", MAX_GRADIENT);
	//TRACE("max_grad2 = %f\n", max_grad2);
	//TRACE("hi_thres = %f\n", hi_thres);
	for (j = 0; j < IMAGE_Y; j++) {
		for (i = 0; i < IMAGE_X; i++) {
			if (i == 0 || i == IMAGE_X-1 || j == 0 || j == IMAGE_Y-1) {
				image2[i][j] = 255;
				thin_edge[i][j] = 0; // init thin edge list
				continue;
			}
			G_mag[i][j] = (G_mag[i][j] / max_grad2); // G_mag between [0, 1]
			image2[i][j] = (int)(G_mag[i][j] * 255);

			//marked[i][j] = 0; // init marked image for hysteresis
			thin_edge[i][j] = 0; // init thin edge list
		}
	}

	double	gx, gy;
	double	g, g1, g2, g3, g4;
	double	t; // interpolation parameter
	
	//////////////////////////////////////////
	// Nonmaxima suppression
	for (j = 1; j < IMAGE_Y-1; j++) {
		for (i = 1; i < IMAGE_X-1; i++) {
			gx = G_x[i][j];
			gy = G_y[i][j];
			g = G_mag[i][j];
			//TRACE("gx = %f, gy = %f, g = %f\n", gx, gy, g);
			//if (gx < 0.01 && gy < 0.01) continue; // not an edge
			if (fabs(gx) >= fabs(gy)) { // close to horizontal (note: gy can be 0)
				t = fabs(gy) / fabs(gx); 
				//TRACE("t = %f\n", t);
				g1 = G_mag[i+1][j]; // right
				g2 = G_mag[i-1][j]; // left
				//TRACE("g1 = %f, g2 = %f\n", g1, g2);
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // right up
					g4 = G_mag[i-1][j-1]; // left down
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i+1][j-1]; // right down
					g4 = G_mag[i-1][j+1]; // left up
				}
			}
			else { // close to vertical (note: gx can be 0)
				t = fabs(gx) / fabs(gy);
				g1 = G_mag[i][j+1]; // up
				g2 = G_mag[i][j-1]; // down
				if (gx*gy >= 0.0) { // 1 or 3 quadrant
					g3 = G_mag[i+1][j+1]; // up right
					g4 = G_mag[i-1][j-1]; // down left
				}
				else { // 2 or 4 quadrant
					g3 = G_mag[i-1][j+1]; // up left
					g4 = G_mag[i+1][j-1]; // down right
				}
			}
			if ( g > ((1-t)*g1 + t*g3) && g > ((1-t)*g2 + t*g4) ) {
                thin_edge[i][j] = 1; // it's a thin edge
				if (pixel_mark[i][j] <= 102) { // foreground pixels
					//TRACE("g = %f\n", g);
					if (g > hi_thres) {
					//if (g > 0.0) {
						//rr = Dbuffer[(j * IMAGE_X + i) * 3 + 0];
						//gg = Dbuffer[(j * IMAGE_X + i) * 3 + 1];
						//bb = Dbuffer[(j * IMAGE_X + i) * 3 + 2];
						/*
						if ( prob2((double)rr, (double)gg, (double)bb, g, off_edge, factor1) )
							off_mark[i][j] = 255; // off edge
						*/
						image2[i][j] = 0; // thin edge above hi_thres
					}
					else 
						image2[i][j] = 255; // thin edge below hi_thres
				}
				else {
					if (g > hi_thres2) {
						image2[i][j] = 0; // thin edge above hi_thres
					}
					else 
						image2[i][j] = 255; // thin edge below hi_thres
				}
			}
			else { // non-maximum
				image2[i][j] = 255;
			}
			
		}
	}

	//////////////////////////////////////////
	// Hysteresis thresholding
	///*
	for (j = 1; j < IMAGE_Y-1; j++) {
		for (i = 1; i < IMAGE_X-1; i++) {
			if (image2[i][j] == 0) { // computed thinned edges above hi_thres
				adaptive_hyster_visit(i, j, image2, lo_thres, lo_thres2); // visit neighboring pixels
			}
		}
	}
	//*/
	//////////////////////////////////////////
	// removing off_edges
	/*
	for (j = 1; j < IMAGE_Y-1; j++) {
		for (i = 1; i < IMAGE_X-1; i++) {
			if (off_mark[i][j] == 255) { // marked off edge
				off_edge_visit(IMAGE_X, IMAGE_Y, i, j); // visit neighboring pixels
			}
		}
	}
	*/
	
	return MAX_GRADIENT;
}

////////////////////////////////////////////////////////////
// time measuring
static time_t StartingTime;   // starting time of clock 

// Function to initialize static global variable StartingTime
// to the number of ticks of the clock
// Pre:  none
// Post: StartingTime has been initialized to clock()
//

void StartTimer(void)
{
   StartingTime = clock();
}

//
// Function to return the elapsed time in seconds since 
// static global variable StartingTime was initialized
// Pre:  StartingTime has been initialized to a value of clock.
// Post: The elapsed time since the initialization of StartingTime 
//       has been returned.
//

double ElapsedTime(void)
{
   //return (double(clock() - StartingTime));
   return (double(clock() - StartingTime)/CLOCKS_PER_SEC);
}

void GetGrayImageFromMemDC(int image_x, int image_y, imatrix& image, CDC& dc)
{
	GLubyte	r, g, b;
	int 	min, max;
	int 	x, y;

	//image.init(image_x, image_y);

	for (y=0; y < image_y; y++) {
		for (x=0; x < image_x; x++) {
			r = (GLubyte)RGB_GETRED(dc.GetPixel(x, image_y-1-y));
			g = (GLubyte)RGB_GETGREEN(dc.GetPixel(x, image_y-1-y));
			b = (GLubyte)RGB_GETBLUE(dc.GetPixel(x, image_y-1-y));
			if (r > g) {
				max = r;
				min = g;
			}
			else {	
				max = g;
				min = r;	
			}
			if (b > max)
				max = b;
			else if (b < min)
				min = b;
			image[x][y] = (int) ((min + max) / 2.0);

		}
	}
}

void GetColorImageFromMemDC(int image_x, int image_y, cimatrix& cmap, CDC& dc)
{
	//GLubyte	r, g, b;
	int 	x, y;

	//image.init(image_x, image_y);

	for (y=0; y < image_y; y++) {
		for (x=0; x < image_x; x++) {
			cmap[x][y].r = (GLubyte)RGB_GETRED(dc.GetPixel(x, image_y-1-y));
			cmap[x][y].g = (GLubyte)RGB_GETGREEN(dc.GetPixel(x, image_y-1-y));
			cmap[x][y].b = (GLubyte)RGB_GETBLUE(dc.GetPixel(x, image_y-1-y));
		}
	}
}
