#include "stdafx.h"
//#include <cstdio>
#include <cstdlib>
#include <cmath>
//#include <iostream>
#include <iomanip>
#include <fstream>
#include <deque>
using namespace std;

#include "Cube.h"
#include "CubeDoc.h"
#include "CubeView.h"
#include "globals.h"
#include "minheap.h"
#include "Stroke.h"
#include "SimpleList.h"
//#include "Heap.h"
#include "Field.h"

#define ABS(x)	( ((x)>0) ? (x) : (-(x)) )

//int max_grad;

extern void ClearMemDC(CDC *dc);

#define RGB_GETRED(rgb)    ((rgb) & 0xff) 
#define RGB_GETGREEN(rgb)    (((rgb) >> 8) & 0xff) 
#define RGB_GETBLUE(rgb)    (((rgb) >> 16) & 0xff)  
#define dist3(x, y, z) sqrt(((double)x)*((double)x) + ((double)y)*((double)y) + ((double)z)*((double)z))
#define dist2(x1, y1, x2, y2) sqrt( (((double)x1)-((double)x2))*(((double)x1)-((double)x2)) + (((double)y1)-((double)y2))*(((double)y1)-((double)y2)) )

#define MAX_OBJ_VALUE  (sqrt((double)3)*255.0*IMAGE_X*IMAGE_Y)

//#define STROKE_CUTLINE 1400
//#define STROKE_CUTLINE 800 // Peppers
//#define STROKE_CUTLINE 700 // Rose
//#define STROKE_CUTLINE 250 // Rose
#define STROKE_CUTLINE 1000 // Rose
//#define STROKE_CUTLINE 900 // Eagle

inline double drand48() 
{
	return (rand()/(double)RAND_MAX); 
}

int intrand(int lo, int hi)
{
	return (int)(0.5 + lo + (hi - lo) * drand48()); 
}      


// uday_data is defined in types.h
uday_data data[500];


void CopyDbuffertoMemDC(CDC& target)
{
	int x, y;
	unsigned char r, g, b;

	for (y = IMAGE_Y - 1; y >= 0; y--) {
		for (x = 0; x < IMAGE_X; x++) {
			r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
			g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
			b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];
			/// Set Pixel in MemDC
			target.SetPixelV(x, (IMAGE_Y-1)-y, RGB(r, g, b));
		}
	}
}

void CopyMemDCtoDbuffer(CDC& dc)
{
	int x, y;
	GLubyte r, g, b;

	for (y = 0; y < IMAGE_Y; y++) {
		for (x = 0; x < IMAGE_X; x++) {
			COLORREF rgbcol = dc.GetPixel(x, (IMAGE_Y-1)-y);
			r = (GLubyte)RGB_GETRED(rgbcol);
			g = (GLubyte)RGB_GETGREEN(rgbcol);
			b = (GLubyte)RGB_GETBLUE(rgbcol);
			Dbuffer[(y * IMAGE_X + x) * 3 + 0] = r;
			Dbuffer[(y * IMAGE_X + x) * 3 + 1] = g;
			Dbuffer[(y * IMAGE_X + x) * 3 + 2] = b;
		}
	}
}

void CopyMemDCtoImatrix(CDC& dc, imatrix& image)
// color image to gray image
{
	int x, y;
	int r;

	image.init(IMAGE_X, IMAGE_Y);

	int image_x = image.getRow();
	int image_y = image.getCol();

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			COLORREF rgbcol = dc.GetPixel(x, (image_y-1)-y);
			r = (int)RGB_GETRED(rgbcol);
			//g = (int)RGB_GETGREEN(rgbcol);
			//b = (int)RGB_GETBLUE(rgbcol);
			image[x][y] = r;
		}
	}
}

void InvertImatrix(imatrix& image)
// white color to black color
{
	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			image[x][y] = 255 - image[x][y];
			if (image[x][y] > 255) image[x][y] = 255;
			if (image[x][y] < 0) image[x][y] = 0;
		}
	}
}





void ClearDC(CDC& dc)
{
	// Create a solid brush object...
	//HBRUSH brush;
	CBrush brush, *OldBrush;
	CPen	pen, *pOldPen;

	brush.CreateSolidBrush(RGB(255, 255, 255)); // color to fill the rectangle
	// ...and attach it into the device context.
	OldBrush = (CBrush *)dc.SelectObject(brush);
	pen.CreatePen(PS_SOLID, 3, RGB(255, 255, 255)); // color to draw the rectangle borders
	pOldPen = (CPen *)dc.SelectObject(&pen);

	dc.Rectangle(0, 0, IMAGE_X-1, IMAGE_Y-1);
	dc.SelectObject(pOldPen);
	dc.SelectObject(OldBrush);
	pen.DeleteObject(); 
	brush.DeleteObject();
}


void PerturbRGB(GLubyte& r, GLubyte& g, GLubyte& b)
{
	int	perturb_r, perturb_g, perturb_b;
	const int MAX = 60; // size of perturbation

	////////////////////////////////////////////////////////
	// Perturb the R, G, B color
	perturb_r = (int) ( (MAX) * (float)rand() / RAND_MAX );
	perturb_r -= MAX/2;
	if (r + perturb_r > 255) { r = 255; } // 0 <= GLubyte <= 255
	else if (r + perturb_r < 0) { r = 0; } // 0 <= GLubyte <= 255
	else r = r + perturb_r;
	perturb_g = (int) ( (MAX) * (float)rand() / RAND_MAX );
	perturb_g -= MAX/2;
	if (g + perturb_g > 255) { g = 255; } // 0 <= GLubyte <= 255
	else if (g + perturb_g < 0) { g = 0; } // 0 <= GLubyte <= 255
	else g = g + perturb_g;
	perturb_b = (int) ( (MAX) * (float)rand() / RAND_MAX );
	perturb_b -= MAX/2;
	if (b + perturb_b > 255) { b = 255; } // 0 <= GLubyte <= 255
	else if (b + perturb_b < 0) { b = 0; } // 0 <= GLubyte <= 255
	else b = b + perturb_b;
	/////////////////////////////////////////////////////////
}

void PerturbRGB(GLubyte& r, GLubyte& g, GLubyte& b, int MAX)
{
	int	perturb_r, perturb_g, perturb_b;
	//const int MAX = 60; // size of perturbation

	////////////////////////////////////////////////////////
	// Perturb the R, G, B color
	perturb_r = (int) ( (MAX) * (float)rand() / RAND_MAX );
	perturb_r -= MAX/2;
	if (r + perturb_r > 255) { r = 255; } // 0 <= GLubyte <= 255
	else if (r + perturb_r < 0) { r = 0; } // 0 <= GLubyte <= 255
	else r = r + perturb_r;
	perturb_g = (int) ( (MAX) * (float)rand() / RAND_MAX );
	perturb_g -= MAX/2;
	if (g + perturb_g > 255) { g = 255; } // 0 <= GLubyte <= 255
	else if (g + perturb_g < 0) { g = 0; } // 0 <= GLubyte <= 255
	else g = g + perturb_g;
	perturb_b = (int) ( (MAX) * (float)rand() / RAND_MAX );
	perturb_b -= MAX/2;
	if (b + perturb_b > 255) { b = 255; } // 0 <= GLubyte <= 255
	else if (b + perturb_b < 0) { b = 0; } // 0 <= GLubyte <= 255
	else b = b + perturb_b;
	/////////////////////////////////////////////////////////
}


CImage m_Image;

///*
int	LoadBMP(CDC& dc, char* filename, CDC *memDC)
// Open BMP File
{
	int		m_nMul, m_nDiv;

	//CClientDC dc(this);

	m_nMul = m_nDiv = 1;
	if ( m_Image.Load(filename) == FALSE ) // the file does not exist
		return 0;

	//IMAGE_X = m_Image.GetWidth();
	bitmap.DeleteObject();
	bitmap.CreateCompatibleBitmap(&dc, m_Image.GetWidth(), m_Image.GetHeight());
	//memDC.CreateCompatibleDC(&dc);
	memDC->SelectObject(&bitmap);

	CRect rcDIB;
	rcDIB.top = rcDIB.left = 0;
	rcDIB.right = m_Image.GetWidth();
	rcDIB.bottom = m_Image.GetHeight();
	CRect rcDest;
	rcDest = rcDIB;
	rcDest.left = rcDest.left * m_nMul / m_nDiv;
	rcDest.top = rcDest.top * m_nMul / m_nDiv;
	rcDest.right= rcDest.right * m_nMul / m_nDiv;
	rcDest.bottom= rcDest.bottom * m_nMul / m_nDiv;

		//m_Image.Draw(pDC->m_hDC,&rcDIB,&rcDest);
		// Draw the BMP image into a memDC!!!
	
	// Draw the BMP image on memDC
	m_Image.Draw(memDC->m_hDC, &rcDIB, &rcDest);
//		dc.BitBlt(0, 0, m_Image.GetWidth(), m_Image.GetHeight(), 
//			memDC, 0, 0, SRCCOPY);

	return 1; // file does exist
}
//*/

int	LoadBMP3(CDC& dc, char* filename, CDC *memDC)
// Open BMP File
{
	int	m_nMul, m_nDiv;

	//CClientDC dc(this);

	m_nMul = m_nDiv = 1;
	if ( m_Image.Load(filename) == FALSE ) // the file does not exist
		return 0;

	//IMAGE_X = m_Image.GetWidth();
	bitmap.DeleteObject();
	bitmap.CreateCompatibleBitmap(&dc, m_Image.GetWidth(), m_Image.GetHeight());
	//memDC.CreateCompatibleDC(&dc);
	memDC->SelectObject(&bitmap);

	CRect rcDIB;
	rcDIB.top = rcDIB.left = 0;
	rcDIB.right = m_Image.GetWidth();
	rcDIB.bottom = m_Image.GetHeight();
	CRect rcDest;
	rcDest = rcDIB;
	rcDest.left = rcDest.left * m_nMul / m_nDiv;
	rcDest.top = rcDest.top * m_nMul / m_nDiv;
	rcDest.right= rcDest.right * m_nMul / m_nDiv;
	rcDest.bottom= rcDest.bottom * m_nMul / m_nDiv;

		//m_Image.Draw(pDC->m_hDC,&rcDIB,&rcDest);
		// Draw the BMP image into a memDC!!!
	
	// Draw the BMP image on memDC
	m_Image.Draw(memDC->m_hDC, &rcDIB, &rcDest);
//		dc.BitBlt(0, 0, m_Image.GetWidth(), m_Image.GetHeight(), 
//			memDC, 0, 0, SRCCOPY);

	return 1; // file does exist
}

int	LoadBMP4(CDC& dc, char* filename, CDC *memDC, int &image_x, int &image_y)
// Open BMP File
{
	int	m_nMul, m_nDiv;

	//CClientDC dc(this);

	m_nMul = m_nDiv = 1;
	if ( m_Image.Load(filename) == FALSE ) // the file does not exist
		return 0;

	image_x = m_Image.GetWidth(); 
	image_y = m_Image.GetHeight();

	//IMAGE_X = m_Image.GetWidth();
	bitmap.DeleteObject();
	bitmap.CreateCompatibleBitmap(&dc, m_Image.GetWidth(), m_Image.GetHeight());
	//memDC.CreateCompatibleDC(&dc);
	memDC->SelectObject(&bitmap);

	CRect rcDIB;
	rcDIB.top = rcDIB.left = 0;
	rcDIB.right = m_Image.GetWidth();
	rcDIB.bottom = m_Image.GetHeight();
	CRect rcDest;
	rcDest = rcDIB;
	rcDest.left = rcDest.left * m_nMul / m_nDiv;
	rcDest.top = rcDest.top * m_nMul / m_nDiv;
	rcDest.right= rcDest.right * m_nMul / m_nDiv;
	rcDest.bottom= rcDest.bottom * m_nMul / m_nDiv;

		//m_Image.Draw(pDC->m_hDC,&rcDIB,&rcDest);
		// Draw the BMP image into a memDC!!!
	
	// Draw the BMP image on memDC
	m_Image.Draw(memDC->m_hDC, &rcDIB, &rcDest);
//		dc.BitBlt(0, 0, m_Image.GetWidth(), m_Image.GetHeight(), 
//			memDC, 0, 0, SRCCOPY);

	return 1; // file does exist
}

int key_pressed = 0;

//////////////////////////////////////////////////////////////////////////////
RGB *pal, *pal_s;
XYRGB *pal2;

#define dist5(x, y, r, g, b) sqrt(((double)x)*((double)x) + ((double)y)*((double)y) + ((double)r)*((double)r) + ((double)g)*((double)g) + ((double)b)*((double)b))


void SelectColors(CDC& dc, int N)
// maximize the "average" distance between colors
{
	int x, y;
	double min_dist;
	double dist;
	int min_i, i, j;
	GLubyte r, g, b;
	GLubyte r2, g2, b2;
	GLubyte old_r, old_g, old_b;
	double fitness, old_fitness;
	int fail_count;
	int flag;

	//RGB *pal;
	int *number;
	RGB *sum;
    
	delete [] pal;
	pal = new RGB[N];
	number = new int[N]; // total number of pixels of the class
	sum = new RGB[N]; // sum of RGB values

	for (i = 0; i < N; i++) {
		while(1) {
			x = intrand(0, IMAGE_X-1); // 0 <= x <= IMAGE_X-1
			y = intrand(0, IMAGE_Y-1);

			flag = 0;
			r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
			g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
			b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];

			for (j = 0; j < i; j++) {
				if ((GLubyte)pal[j].r == r && (GLubyte)pal[j].g == g && (GLubyte)pal[j].b == b)
					flag = 1;
			}
			if (flag) {
				TRACE("[%d] identical_1\n", i);
				continue; // try again
			}
			break;
		}
		pal[i].r = (double)r;
		pal[i].g = (double)g;
		pal[i].b = (double)b;
	}

	dist = 0.0;
	for (i = 0; i < N; i++) {
		r = (GLubyte)pal[i].r;
		g = (GLubyte)pal[i].g;
		b = (GLubyte)pal[i].b;
		for (j = 0; j < N; j++) {
			r2 = (GLubyte)pal[j].r;
			g2 = (GLubyte)pal[j].g;
			b2 = (GLubyte)pal[j].b;
			dist += dist3(r-r2, g-g2, b-b2);
		}
	}
	old_fitness = dist / (N*N);
	//TRACE("fitness = %f\n", fitness);

	// iteratively try new color
	fail_count = 0;
	while (1) {
		while(1) {
			x = intrand(0, IMAGE_X-1); // 0 <= x <= IMAGE_X-1
			y = intrand(0, IMAGE_Y-1);

			flag = 0;
			r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
			g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
			b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];

			for (j = 0; j < N; j++) {
				if ((GLubyte)pal[j].r == r && (GLubyte)pal[j].g == g && (GLubyte)pal[j].b == b)
					flag = 1;
			}
			if (flag) {
				//TRACE("[%d] identical_2\n", i);
				continue; // try again
			}
			break;
		}
			
		min_i = -1;
		min_dist = 1000000000;
		for (i = 0; i < N; i++) {
			dist = dist3(pal[i].r-(double)r, pal[i].g-(double)g, pal[i].b-(double)b);
			//TRACE("dist = %f\n", dist);
			if ( dist < min_dist ) {
				min_dist = dist;
				min_i = i;
			}
		}
		// store the old values
		old_r = (GLubyte)pal[min_i].r;
		old_g = (GLubyte)pal[min_i].g;
		old_b = (GLubyte)pal[min_i].b;
		pal[min_i].r = (double)r;
		pal[min_i].g = (double)g;
		pal[min_i].b = (double)b;
		//TRACE("pal[%d].r = %f\n", col_count, pal[col_count].r);
		//TRACE("pal[%d].g = %f\n", col_count, pal[col_count].g);
		//TRACE("pal[%d].b = %f\n", col_count, pal[col_count].b);
		/////////////////////////////////////
		// compute new fitness
		dist = 0.0;
		for (i = 0; i < N; i++) {
			r = (GLubyte)pal[i].r;
			g = (GLubyte)pal[i].g;
			b = (GLubyte)pal[i].b;
			for (j = 0; j < N; j++) {
				r2 = (GLubyte)pal[j].r;
				g2 = (GLubyte)pal[j].g;
				b2 = (GLubyte)pal[j].b;
				dist += dist3(r-r2, g-g2, b-b2);
			}
		}
		fitness = dist / (N*N);
		if (fitness > old_fitness) { // update
			fail_count = 0; // reset
			//if (fitness - old_fitness < 0.05) break;
			TRACE("fitness = %f\n", fitness);
			old_fitness = fitness;
			///////////////////////////////////////////////
			// display result
			///*
			for (y = 0; y < IMAGE_Y; y++) {
				for (x = 0; x < IMAGE_X; x++) {
					r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
					g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
					b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];
					min_i = -1;
					min_dist = 1000000000;
					for (i = 0; i < N; i++) {
						dist = dist3(pal[i].r-(double)r, pal[i].g-(double)g, pal[i].b-(double)b);
						//TRACE("dist = %f\n", dist);
						if ( dist < min_dist ) {
							min_dist = dist;
							min_i = i;
						}
					}
					r = (GLubyte)pal[min_i].r;
					g = (GLubyte)pal[min_i].g;
					b = (GLubyte)pal[min_i].b;
					memDC2.SetPixelV(x, IMAGE_Y-1-y, RGB(r, g, b));
				}
			}
			dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC2, 0, 0, SRCCOPY); 
			//*/
		}
		else { // restore
			fail_count++;
			pal[min_i].r = (double)old_r;
			pal[min_i].g = (double)old_g;
			pal[min_i].b = (double)old_b;
			if (fail_count > 1000) break;
		}
	}

}

void SelectColors2(CDC& dc, int N)
// maximize the "minimum" distance between colors
{
	int x, y;
	double min_dist;
	double dist;
	int min_i, i, j, min_k;
	GLubyte r, g, b;
	GLubyte r2, g2, b2;
	GLubyte old_r, old_g, old_b;
	double fitness, old_fitness;
	int fail_count;
	int flag;

	//RGB *pal;
	int *number;
	RGB *sum;
    
	delete [] pal;
	pal = new RGB[N];
	number = new int[N]; // total number of pixels of the class
	sum = new RGB[N]; // sum of RGB values

	///////////////////////////////////////
	// select random N colors (but no identical ones)
	for (i = 0; i < N; i++) {
		while(1) {
			x = intrand(0, IMAGE_X-1); // 0 <= x <= IMAGE_X-1
			y = intrand(0, IMAGE_Y-1);

			flag = 0;
			r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
			g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
			b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];

			for (j = 0; j < i; j++) {
				if ((GLubyte)pal[j].r == r && (GLubyte)pal[j].g == g && (GLubyte)pal[j].b == b)
					flag = 1;
			}
			if (flag) {
				TRACE("[%d]identical\n", i);
				continue; // try again
			}
			break;
		}

		pal[i].r = (double)r;
		pal[i].g = (double)g;
		pal[i].b = (double)b;
	}
	///////////////////////////////////

	min_dist = 1000000.0;
	//min_i = -1;
	for (i = 0; i < N; i++) {
		r = (GLubyte)pal[i].r;
		g = (GLubyte)pal[i].g;
		b = (GLubyte)pal[i].b;
		for (j = 0; j < N; j++) {
			if (i == j) continue;
			r2 = (GLubyte)pal[j].r;
			g2 = (GLubyte)pal[j].g;
			b2 = (GLubyte)pal[j].b;
			dist = dist3(r-r2, g-g2, b-b2);
			if (dist < min_dist) {
				min_dist = dist;
				//min_i = j;
			}
		}
	}
	//old_fitness = dist / (N*N);
	old_fitness = min_dist;
	TRACE("old_fitness = %f\n", old_fitness);

	// iteratively try new color
	fail_count = 0;
	while (1) {
		while(1) {
			x = intrand(0, IMAGE_X-1); // 0 <= x <= IMAGE_X-1
			y = intrand(0, IMAGE_Y-1);

			flag = 0;
			r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
			g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
			b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];

			for (j = 0; j < N; j++) {
				if ((GLubyte)pal[j].r == r && (GLubyte)pal[j].g == g && (GLubyte)pal[j].b == b)
					flag = 1;
			}
			if (flag) {
				TRACE("[%d]identical\n", i);
				continue; // try again
			}
			break;
		}
	
		min_dist = 1000000000;
		min_i = -1;
		for (i = 0; i < N; i++) {
			dist = dist3(pal[i].r-(double)r, pal[i].g-(double)g, pal[i].b-(double)b);
			//TRACE("dist = %f\n", dist);
			if ( dist < min_dist ) {
				min_dist = dist;
				min_i = i;
			}
		}
		//min_i = intrand(0, N-1); // 0 <= x <= N-1
		// store the old values
		old_r = (GLubyte)pal[min_i].r;
		old_g = (GLubyte)pal[min_i].g;
		old_b = (GLubyte)pal[min_i].b;
		// update pal
		pal[min_i].r = (double)r;
		pal[min_i].g = (double)g;
		pal[min_i].b = (double)b;
		//TRACE("pal[%d].r = %f\n", col_count, pal[col_count].r);
		//TRACE("pal[%d].g = %f\n", col_count, pal[col_count].g);
		//TRACE("pal[%d].b = %f\n", col_count, pal[col_count].b);
		/////////////////////////////////////
		// compute new fitness
		min_dist = 1000000000;
		for (i = 0; i < N; i++) {
			r = (GLubyte)pal[i].r;
			g = (GLubyte)pal[i].g;
			b = (GLubyte)pal[i].b;
			for (j = 0; j < N; j++) {
				if (i == j) continue;
				r2 = (GLubyte)pal[j].r;
				g2 = (GLubyte)pal[j].g;
				b2 = (GLubyte)pal[j].b;
				dist = dist3(r-r2, g-g2, b-b2);
				if (dist < min_dist) {
					min_dist = dist;
					//TRACE("min_dist = %f\n", min_dist);
				}
			}
		}
		fitness = min_dist;
		if (fitness > old_fitness) { // update
			fail_count = 0; // reset
			//if (fitness - old_fitness < 0.05) break;
			TRACE("fitness = %f\n", fitness);
			old_fitness = fitness;
			///////////////////////////////////////////////
			// display result
			///*
			for (y = 0; y < IMAGE_Y; y++) {
				for (x = 0; x < IMAGE_X; x++) {
					r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
					g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
					b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];
					min_k = -1;
					min_dist = 1000000000;
					for (i = 0; i < N; i++) {
						dist = dist3(pal[i].r-(double)r, pal[i].g-(double)g, pal[i].b-(double)b);
						//TRACE("dist = %f\n", dist);
						if ( dist < min_dist ) {
							min_dist = dist;
							min_k = i;
						}
					}
					r = (GLubyte)pal[min_k].r;
					g = (GLubyte)pal[min_k].g;
					b = (GLubyte)pal[min_k].b;
					memDC2.SetPixelV(x, IMAGE_Y-1-y, RGB(r, g, b));
				}
			}
			dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC2, 0, 0, SRCCOPY); 
			//*/
		}
		else { // restore
			//TRACE("failure: fitness = %f\n", fitness);
			fail_count++;
			pal[min_i].r = (double)old_r;
			pal[min_i].g = (double)old_g;
			pal[min_i].b = (double)old_b;
			if (fail_count > 10000) break;
		}
	}

}

void SelectColors3(CDC& dc, int N)
// maximize the "average" and "minimum" distance between colors
{
	int x, y;
	double min_dist;
	double dist, tot_dist;
	int min_i, i, j;
	GLubyte r, g, b;
	GLubyte r2, g2, b2;
	GLubyte old_r, old_g, old_b;
	double fitness, old_fitness;
	int fail_count;
	int flag;
	double t = 0.1;

	//RGB *pal;
	int *number;
	RGB *sum;
    
	delete [] pal;
	pal = new RGB[N];
	number = new int[N]; // total number of pixels of the class
	sum = new RGB[N]; // sum of RGB values

	for (i = 0; i < N; i++) {
		while(1) {
			x = intrand(0, IMAGE_X-1); // 0 <= x <= IMAGE_X-1
			y = intrand(0, IMAGE_Y-1);

			flag = 0;
			r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
			g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
			b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];

			for (j = 0; j < i; j++) {
				if ((GLubyte)pal[j].r == r && (GLubyte)pal[j].g == g && (GLubyte)pal[j].b == b)
					flag = 1;
			}
			if (flag) {
				TRACE("[%d] identical_1\n", i);
				continue; // try again
			}
			break;
		}
		pal[i].r = (double)r;
		pal[i].g = (double)g;
		pal[i].b = (double)b;
	}

	tot_dist = 0.0;
	min_dist = 100000000;
	for (i = 0; i < N; i++) {
		r = (GLubyte)pal[i].r;
		g = (GLubyte)pal[i].g;
		b = (GLubyte)pal[i].b;
		for (j = 0; j < N; j++) {
			if (i == j) continue;
			r2 = (GLubyte)pal[j].r;
			g2 = (GLubyte)pal[j].g;
			b2 = (GLubyte)pal[j].b;
			dist = dist3(r-r2, g-g2, b-b2);
			tot_dist += dist;
			if (dist < min_dist) {
				min_dist = dist;
				//min_i = j;
			}
		}
	}
	old_fitness = t * tot_dist / (N*N) + (1-t) * min_dist;
	TRACE("old_fitness = %f\n", old_fitness);
	TRACE("min_dist = %f\n", min_dist);

	// iteratively try new color
	fail_count = 0;
	while (1) {
		while(1) {
			x = intrand(0, IMAGE_X-1); // 0 <= x <= IMAGE_X-1
			y = intrand(0, IMAGE_Y-1);

			flag = 0;
			r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
			g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
			b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];

			for (j = 0; j < N; j++) {
				if ((GLubyte)pal[j].r == r && (GLubyte)pal[j].g == g && (GLubyte)pal[j].b == b)
					flag = 1;
			}
			if (flag) {
				//TRACE("[%d] identical_2\n", i);
				continue; // try again
			}
			break;
		}
			
		min_i = -1;
		min_dist = 1000000000;
		for (i = 0; i < N; i++) {
			dist = dist3(pal[i].r-(double)r, pal[i].g-(double)g, pal[i].b-(double)b);
			//TRACE("dist = %f\n", dist);
			if ( dist < min_dist ) {
				min_dist = dist;
				min_i = i;
			}
		}
		// store the old values
		old_r = (GLubyte)pal[min_i].r;
		old_g = (GLubyte)pal[min_i].g;
		old_b = (GLubyte)pal[min_i].b;
		pal[min_i].r = (double)r;
		pal[min_i].g = (double)g;
		pal[min_i].b = (double)b;
		//TRACE("pal[%d].r = %f\n", col_count, pal[col_count].r);
		//TRACE("pal[%d].g = %f\n", col_count, pal[col_count].g);
		//TRACE("pal[%d].b = %f\n", col_count, pal[col_count].b);
		/////////////////////////////////////
		// compute new fitness
		tot_dist = 0.0;
		min_dist = 10000000;
		for (i = 0; i < N; i++) {
			r = (GLubyte)pal[i].r;
			g = (GLubyte)pal[i].g;
			b = (GLubyte)pal[i].b;
			for (j = 0; j < N; j++) {
				if (i == j) continue;
				r2 = (GLubyte)pal[j].r;
				g2 = (GLubyte)pal[j].g;
				b2 = (GLubyte)pal[j].b;
				dist = dist3(r-r2, g-g2, b-b2);
				tot_dist += dist;
				if (dist < min_dist) {
					min_dist = dist;
				}
			}
		}
		fitness = t * tot_dist / (N*N) + (1-t) * min_dist;
		//TRACE("min_dist = %f\n", min_dist);
		if (fitness > old_fitness) { // update
			fail_count = 0; // reset
			//if (fitness - old_fitness < 0.05) break;
			TRACE("fitness = %f\n", fitness);
			old_fitness = fitness;
			///////////////////////////////////////////////
			// display result
			///*
			for (y = 0; y < IMAGE_Y; y++) {
				for (x = 0; x < IMAGE_X; x++) {
					r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
					g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
					b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];
					min_i = -1;
					min_dist = 1000000000;
					for (i = 0; i < N; i++) {
						dist = dist3(pal[i].r-(double)r, pal[i].g-(double)g, pal[i].b-(double)b);
						//TRACE("dist = %f\n", dist);
						if ( dist < min_dist ) {
							min_dist = dist;
							min_i = i;
						}
					}
					r = (GLubyte)pal[min_i].r;
					g = (GLubyte)pal[min_i].g;
					b = (GLubyte)pal[min_i].b;
					memDC2.SetPixelV(x, IMAGE_Y-1-y, RGB(r, g, b));
				}
			}
			dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC2, 0, 0, SRCCOPY); 
			//*/
		}
		else { // restore
			fail_count++;
			pal[min_i].r = (double)old_r;
			pal[min_i].g = (double)old_g;
			pal[min_i].b = (double)old_b;
			if (fail_count > 1000) break;
		}
	}

}

void AdjustColor(int& col_count, int size, GLubyte& r, GLubyte& g, GLubyte& b)
{
	int i;
	double min_dist = 1000000000;
	double dist;
	int min_i;

	min_i = -1;
	for (i = 0; i < col_count; i++) {
		dist = dist3(pal[i].r-(double)r, pal[i].g-(double)g, pal[i].b-(double)b);
		//TRACE("dist = %f\n", dist);
		if ( dist < min_dist ) {
			min_dist = dist;
			min_i = i;
		}
	}
	if (min_dist < 50 || col_count >= size) { // use stored rgb
		r = (GLubyte)pal[min_i].r;
		g = (GLubyte)pal[min_i].g;
		b = (GLubyte)pal[min_i].b;
	} 
	else { // add new entry
		pal[col_count].r = (double)r;
		pal[col_count].g = (double)g;
		pal[col_count].b = (double)b;
		col_count++;
	}

}

void AdjustColor2(int N, GLubyte& r, GLubyte& g, GLubyte& b)
{
	int i;
	double min_dist = 1000000000;
	double dist;
	int min_i;

	min_i = -1;
	for (i = 0; i < N; i++) {
		//TRACE("pal[%d].r = %f\n", i, pal[i].r);
		//TRACE("pal[%d].g = %f\n", i, pal[i].g);
		//TRACE("pal[%d].b = %f\n", i, pal[i].b);
		dist = dist3(pal[i].r-(double)r, pal[i].g-(double)g, pal[i].b-(double)b);
		//TRACE("dist = %f\n", dist);
		if ( dist < min_dist ) {
			min_dist = dist;
			min_i = i;
		}
	}
	r = (GLubyte)pal[min_i].r;
	g = (GLubyte)pal[min_i].g;
	b = (GLubyte)pal[min_i].b;
}

void AdjustColor5(int N, int x, int y, GLubyte& r, GLubyte& g, GLubyte& b)
// x, y, r, g, b
{
	int i;
	double min_dist = 1000000000;
	double dist;
	int min_i;

	min_i = -1;
	for (i = 0; i < N; i++) {
		//TRACE("pal[%d].r = %f\n", i, pal[i].r);
		//TRACE("pal[%d].g = %f\n", i, pal[i].g);
		//TRACE("pal[%d].b = %f\n", i, pal[i].b);
		//dist = dist3(pal[i].r-(double)r, pal[i].g-(double)g, pal[i].b-(double)b);
		dist = dist5(pal2[i].x-(double)x, pal2[i].y-(double)y, pal2[i].r-(double)r, pal2[i].g-(double)g, pal2[i].b-(double)b);
		//TRACE("dist = %f\n", dist);
		if ( dist < min_dist ) {
			min_dist = dist;
			min_i = i;
		}
	}
	r = (GLubyte)pal2[min_i].r;
	g = (GLubyte)pal2[min_i].g;
	b = (GLubyte)pal2[min_i].b;
}

void AdjustColor3(int N, GLubyte& r, GLubyte& g, GLubyte& b)
{
	int i;
	double min_dist = 1000000000;
	double dist;
	int min_i;

	min_i = -1;
	for (i = 0; i < N; i++) {
		//TRACE("pal[%d].r = %f\n", i, pal[i].r);
		//TRACE("pal[%d].g = %f\n", i, pal[i].g);
		//TRACE("pal[%d].b = %f\n", i, pal[i].b);
		dist = dist3(pal[i].r-(double)r, pal[i].g-(double)g, pal[i].b-(double)b);
		//TRACE("dist = %f\n", dist);
		if ( dist < min_dist ) {
			min_dist = dist;
			min_i = i;
		}
	}
	r = (GLubyte)pal_s[min_i].r;
	g = (GLubyte)pal_s[min_i].g;
	b = (GLubyte)pal_s[min_i].b;
}

#define min3(x, y, z) ((x) < (y))? (((x) < (z))? (x) : (z)) : (((y) < (z))? (y) : (z))
#define max3(x, y, z) ((x) > (y))? (((x) > (z))? (x) : (z)) : (((y) > (z))? (y) : (z))


inline void RGB2HSV( double r, double g, double b, double& h, double& s, double& v )
// h = [0,360], s = [0,1], v = [0,1]
// r,g,b values are from 0 to 1
{
	double min, max, delta;

	// r,g,b values are from 0 to 1
	// h = [0,360], s = [0,1], v = [0,1]
	//		if s == 0, then h = -1 (undefined)
	min = min3( r, g, b );
	max = max3( r, g, b );
    //TRACE("[r, g, b] = [%f, %f, %f]\n", r, g, b);
	//TRACE("min = %f\n", min);
	//TRACE("max = %f\n", max);
	v = max;				// v
	delta = max - min;
	if( max != 0 )
		s = delta / max;		// s
	else {
		// r = g = b = 0		// s = 0, v is undefined
		s = 0;
		h = -1;
		return;
	}
	if( r == max )
		h = ( g - b ) / delta;		// between yellow & magenta
	else if( g == max )
		h = 2 + ( b - r ) / delta;	// between cyan & yellow
	else
		h = 4 + ( r - g ) / delta;	// between magenta & cyan
	h *= 60;				// degrees
	if( h < 0 )
		h += 360;
}

inline void HSV2RGB( double h, double s, double v, double& r, double& g, double& b )
// h = [0,360], s = [0,1], v = [0,1]
// r,g,b values are from 0 to 1
{
	int i;
	double f, p, q, t;
	if( s == 0 ) {
		// achromatic (grey)
		r = g = b = v;
		return;
	}
	h /= 60;			// sector 0 to 5
	i = (int)floor( h );
	f = h - i;			// factorial part of h
	p = v * ( 1 - s );
	q = v * ( 1 - s * f );
	t = v * ( 1 - s * ( 1 - f ) );
	switch( i ) {
		case 0:
			r = v;
			g = t;
			b = p;
			break;
		case 1:
			r = q;
			g = v;
			b = p;
			break;
		case 2:
			r = p;
			g = v;
			b = t;
			break;
		case 3:
			r = p;
			g = q;
			b = v;
			break;
		case 4:
			r = t;
			g = p;
			b = v;
			break;
		default:		// case 5:
			r = v;
			g = p;
			b = q;
			break;
	}
}

inline void RGB2LAB( double R, double G, double B, double& L, double& a, double& b )
// R: [0, 255], G: [0, 255], B: [0, 255] 
// L: [0, 100], a: [-86.184636, 98.254219], b: [-107.863681, 94.482485] 
{
	double XX, YY, ZZ;
	double var_R, var_G, var_B;
	double var_X, var_Y, var_Z;
	double ref_X, ref_Y, ref_Z;

	var_R = ( R / 255. );        //Where R = [0, 255]
	var_G = ( G / 255. );        //Where G = [0, 255]
	var_B = ( B / 255. );        //Where B = [0, 255]

	if ( var_R > 0.04045 ) var_R = pow( ( var_R + 0.055 ) / 1.055 , 2.4 );
	else                   var_R = var_R / 12.92;
	if ( var_G > 0.04045 ) var_G = pow( ( var_G + 0.055 ) / 1.055 , 2.4 );
	else                   var_G = var_G / 12.92;
	if ( var_B > 0.04045 ) var_B = pow( ( var_B + 0.055 ) / 1.055 , 2.4 );
	else                   var_B = var_B / 12.92;

	var_R = var_R * 100.;
	var_G = var_G * 100.;
	var_B = var_B * 100.;

	//Observer. = 2¡Æ, Illuminant = D65
	XX = var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805;
	YY = var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722;
	ZZ = var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505;

	///////////////////////////////////////
	// XYZ to CIELab

	ref_X = 95.047;
	ref_Y = 100.0;
	ref_Z = 108.883;

	var_X = XX / ref_X;          //ref_X =  95.047  Observer= 2¡Æ, Illuminant= D65
	var_Y = YY / ref_Y;          //ref_Y = 100.000
	var_Z = ZZ / ref_Z;          //ref_Z = 108.883

	if ( var_X > 0.008856 ) var_X = pow( var_X, 1/3. );
	else                    var_X = ( 7.787 * var_X ) + ( 16 / 116. );
	if ( var_Y > 0.008856 ) var_Y = pow( var_Y, 1/3. );
	else                    var_Y = ( 7.787 * var_Y ) + ( 16 / 116. );
	if ( var_Z > 0.008856 ) var_Z = pow( var_Z, 1/3. );
	else                    var_Z = ( 7.787 * var_Z ) + ( 16 / 116. );

	L = ( 116. * var_Y ) - 16;
	a = 500. * ( var_X - var_Y );
	b = 200. * ( var_Y - var_Z );

	//TRACE("R = %f G = %f B = %f\n", R, G, B);
	//TRACE("L = %f a = %f b = %f\n", L, a, b);
}

inline void LAB2RGB( double L, double a, double b, double& R, double& G, double& B )
// L: [0, 100], a: [-86.184636, 98.254219], b: [-107.863681, 94.482485] 
// R: [0, 255], G: [0, 255], B: [0, 255] 
{
	double XX, YY, ZZ;
	double var_R, var_G, var_B;
	double var_X, var_Y, var_Z;
	double ref_X, ref_Y, ref_Z;

	var_Y = ( L + 16. ) / 116.;
	var_X = a / 500. + var_Y;
	var_Z = var_Y - b / 200.;

	if ( pow(var_Y, 3.0) > 0.008856 ) var_Y = pow(var_Y, 3.0);
	else                      var_Y = ( var_Y - 16 / 116. ) / 7.787;
	if ( pow(var_X, 3.0) > 0.008856 ) var_X = pow(var_X, 3.0);
	else                      var_X = ( var_X - 16 / 116. ) / 7.787;
	if ( pow(var_Z, 3.0) > 0.008856 ) var_Z = pow(var_Z, 3.0);
	else                      var_Z = ( var_Z - 16 / 116. ) / 7.787;

	ref_X = 95.047;
	ref_Y = 100.0;
	ref_Z = 108.883;

	XX = ref_X * var_X;     //ref_X =  95.047  Observer= 2¡Æ, Illuminant= D65
	YY = ref_Y * var_Y;     //ref_Y = 100.000
	ZZ = ref_Z * var_Z;     //ref_Z = 108.883

	///////////////////////////////////////
	// XYZ to RGB
	
	var_X = XX / 100.;        //Where X = [0, 95.047]
	var_Y = YY / 100.;        //Where Y = [0, 100.000]
	var_Z = ZZ / 100.;        //Where Z = [0, 108.883]

	var_R = var_X *  3.2406 + var_Y * -1.5372 + var_Z * -0.4986;
	var_G = var_X * -0.9689 + var_Y *  1.8758 + var_Z *  0.0415;
	var_B = var_X *  0.0557 + var_Y * -0.2040 + var_Z *  1.0570;

	if ( var_R > 0.0031308 ) var_R = 1.055 * pow(var_R, 1/2.4) - 0.055;
	else                     var_R = 12.92 * var_R;
	if ( var_G > 0.0031308 ) var_G = 1.055 * pow(var_G, 1/2.4) - 0.055;
	else                     var_G = 12.92 * var_G;
	if ( var_B > 0.0031308 ) var_B = 1.055 * pow(var_B, 1/2.4) - 0.055;
	else                     var_B = 12.92 * var_B;

	R = var_R * 255;
	G = var_G * 255;
	B = var_B * 255; 

	if (R > 255.0) R = 255.0; if (R < 0.0) R = 0.0;
	if (G > 255.0) G = 255.0; if (G < 0.0) G = 0.0;
	if (B > 255.0) B = 255.0; if (B < 0.0) B = 0.0;

}

inline void RGB2LAB256( double R, double G, double B, double& L, double& a, double& b )
// R: [0, 255], G: [0, 255], B: [0, 255] 
// L: [0, 255], a: [0, 255], b: [0, 255] 
{
	double XX, YY, ZZ;
	double var_R, var_G, var_B;
	double var_X, var_Y, var_Z;
	double ref_X, ref_Y, ref_Z;
	double max_L, max_a, max_b, min_L, min_a, min_b;

	var_R = ( R / 255. );        //Where R = [0, 255]
	var_G = ( G / 255. );        //Where G = [0, 255]
	var_B = ( B / 255. );        //Where B = [0, 255]

	if ( var_R > 0.04045 ) var_R = pow( ( var_R + 0.055 ) / 1.055 , 2.4 );
	else                   var_R = var_R / 12.92;
	if ( var_G > 0.04045 ) var_G = pow( ( var_G + 0.055 ) / 1.055 , 2.4 );
	else                   var_G = var_G / 12.92;
	if ( var_B > 0.04045 ) var_B = pow( ( var_B + 0.055 ) / 1.055 , 2.4 );
	else                   var_B = var_B / 12.92;

	var_R = var_R * 100.;
	var_G = var_G * 100.;
	var_B = var_B * 100.;

	//Observer. = 2¡Æ, Illuminant = D65
	XX = var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805;
	YY = var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722;
	ZZ = var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505;

	///////////////////////////////////////
	// XYZ to CIELab

	ref_X = 95.047;
	ref_Y = 100.0;
	ref_Z = 108.883;

	var_X = XX / ref_X;          //ref_X =  95.047  Observer= 2¡Æ, Illuminant= D65
	var_Y = YY / ref_Y;          //ref_Y = 100.000
	var_Z = ZZ / ref_Z;          //ref_Z = 108.883

	if ( var_X > 0.008856 ) var_X = pow( var_X, 1/3. );
	else                    var_X = ( 7.787 * var_X ) + ( 16 / 116. );
	if ( var_Y > 0.008856 ) var_Y = pow( var_Y, 1/3. );
	else                    var_Y = ( 7.787 * var_Y ) + ( 16 / 116. );
	if ( var_Z > 0.008856 ) var_Z = pow( var_Z, 1/3. );
	else                    var_Z = ( 7.787 * var_Z ) + ( 16 / 116. );

	L = ( 116. * var_Y ) - 16;
	a = 500. * ( var_X - var_Y );
	b = 200. * ( var_Y - var_Z );

	// L: [0, 100], a: [-86.184636, 98.254219], b: [-107.863681, 94.482485] 
	min_L = 0.0;
	max_L = 100.0;
	min_a = -86.184636;
	max_a = 98.254219;
	min_b = -107.863681;
	max_b = 94.482485;

	L = 255. * ( (L + ABS(min_L)) / (max_L - min_L) ); // [0, 255];
	a = 255. * ( (a + ABS(min_a)) / (max_a - min_a) ); // [0, 255];
	b = 255. * ( (b + ABS(min_b)) / (max_b - min_b) ); // [0, 255];

	if (L > 255.) L = 255; if (L < 0.0) L = 0.0;
	if (a > 255.) a = 255; if (a < 0.0) a = 0.0;
	if (b > 255.) b = 255; if (b < 0.0) b = 0.0;

	//TRACE("R = %f G = %f B = %f\n", R, G, B);
	//TRACE("L = %f a = %f b = %f\n", L, a, b);
}

inline void LAB2RGB256( double L, double a, double b, double& R, double& G, double& B )
// L: [0, 255], a: [0, 255], b: [0, 255] 
// R: [0, 255], G: [0, 255], B: [0, 255] 
{
	double XX, YY, ZZ;
	double var_R, var_G, var_B;
	double var_X, var_Y, var_Z;
	double ref_X, ref_Y, ref_Z;
	double max_L, max_a, max_b, min_L, min_a, min_b;

	min_L = 0.0;
	max_L = 100.0;
	min_a = -86.184636;
	max_a = 98.254219;
	min_b = -107.863681;
	max_b = 94.482485;

	L = ( (L / 255.) * (max_L - min_L) ) - ABS(min_L); // [0, 100]
	a = ( (a / 255.) * (max_a - min_a) ) - ABS(min_a); // [-86.184636, 98.254219]
	b = ( (b / 255.) * (max_b - min_b) ) - ABS(min_b); // [-107.863681, 94.482485]

	var_Y = ( L + 16. ) / 116.;
	var_X = a / 500. + var_Y;
	var_Z = var_Y - b / 200.;

	if ( pow(var_Y, 3.0) > 0.008856 ) var_Y = pow(var_Y, 3.0);
	else                      var_Y = ( var_Y - 16 / 116. ) / 7.787;
	if ( pow(var_X, 3.0) > 0.008856 ) var_X = pow(var_X, 3.0);
	else                      var_X = ( var_X - 16 / 116. ) / 7.787;
	if ( pow(var_Z, 3.0) > 0.008856 ) var_Z = pow(var_Z, 3.0);
	else                      var_Z = ( var_Z - 16 / 116. ) / 7.787;

	ref_X = 95.047;
	ref_Y = 100.0;
	ref_Z = 108.883;

	XX = ref_X * var_X;     //ref_X =  95.047  Observer= 2¡Æ, Illuminant= D65
	YY = ref_Y * var_Y;     //ref_Y = 100.000
	ZZ = ref_Z * var_Z;     //ref_Z = 108.883

	///////////////////////////////////////
	// XYZ to RGB
	
	var_X = XX / 100.;        //Where X = [0, 95.047]
	var_Y = YY / 100.;        //Where Y = [0, 100.000]
	var_Z = ZZ / 100.;        //Where Z = [0, 108.883]

	var_R = var_X *  3.2406 + var_Y * -1.5372 + var_Z * -0.4986;
	var_G = var_X * -0.9689 + var_Y *  1.8758 + var_Z *  0.0415;
	var_B = var_X *  0.0557 + var_Y * -0.2040 + var_Z *  1.0570;

	if ( var_R > 0.0031308 ) var_R = 1.055 * pow(var_R, 1/2.4) - 0.055;
	else                     var_R = 12.92 * var_R;
	if ( var_G > 0.0031308 ) var_G = 1.055 * pow(var_G, 1/2.4) - 0.055;
	else                     var_G = 12.92 * var_G;
	if ( var_B > 0.0031308 ) var_B = 1.055 * pow(var_B, 1/2.4) - 0.055;
	else                     var_B = 12.92 * var_B;

	R = var_R * 255;
	G = var_G * 255;
	B = var_B * 255; 

	if (R > 255.) R = 255.; if (R < 0.0) R = 0.0;
	if (G > 255.) G = 255.; if (G < 0.0) G = 0.0;
	if (B > 255.) B = 255.; if (B < 0.0) B = 0.0;
}

void TestRGB2LAB()
// test it for all RGB combinations
{
	double R, G, B, L, AA, BB;
	int r, g, b;
	double max_L, max_a, max_b;
	double min_L, min_a, min_b;
	double max_R, max_G, max_B;
	double min_R, min_G, min_B;

	max_L = max_a = max_b = max_R = max_G = max_B = -1000;
	min_L = min_a = min_b = min_R = min_G = min_B = 1000;
	for (r = 0; r < 256; r++)
		for (g = 0; g < 256; g++)
			for (b = 0; b < 256; b++) {
				R = (double)r, G = (double)g, B = (double)b;
				RGB2LAB256(R, G, B, L, AA, BB);
				//RGB2LAB2(R, G, B, L, AA, BB);
				//TRACE("[%d, %d, %d]->[%f,%f,%f]\n", r, g, b, L, AA, BB);
				if (L > max_L) max_L = L;
				if (L < min_L) min_L = L;
				if (AA > max_a) max_a = AA;
				if (AA < min_a) min_a = AA;
				if (BB > max_b) max_b = BB;
				if (BB < min_b) min_b = BB;
				LAB2RGB256(L, AA, BB, R, G, B);
				if (R > max_R) max_R = R;
				if (R < min_R) min_R = R;
				if (G > max_G) max_G = G;
				if (G < min_G) min_G = G;
				if (B > max_B) max_B = B;
				if (B < min_B) min_B = B;

			}
	TRACE("max LAB: [%f,%f,%f]\n", max_L, max_a, max_b);
	TRACE("min LAB: [%f,%f,%f]\n", min_L, min_a, min_b);
	TRACE("max RGB: [%f,%f,%f]\n", max_R, max_G, max_G);
	TRACE("min RGB: [%f,%f,%f]\n", min_R, min_G, min_G);
}

void ConvertRGBtoLAB(cimatrix& cmap, cimatrix& target)
{
	int i, j;
	int image_x, image_y;
	double R, G, B;
	double L, a, b;

	double max_L, max_a, max_b;
	double min_L, min_a, min_b;
	
	max_L = max_a = max_b = -1000;
	min_L = min_a = min_b = 1000;

	image_x = cmap.getRow();
	image_y = cmap.getCol();
	
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			R = (double)cmap[i][j].r;
			G = (double)cmap[i][j].g;
			B = (double)cmap[i][j].b;
			RGB2LAB(R, G, B, L, a, b);
			//TRACE("L = %f, a = %f, b = %f\n", L, a, b);
			//TRACE("r = %d, g = %d, b = %d\n", (int)R, (int)G, (int)B);
			/*
			if (L > max_L) max_L = L;
			if (a > max_a) max_a = a;
			if (b > max_b) max_b = b;
			if (L < min_L) min_L = L;
			if (a < min_a) min_a = a;
			if (b < min_b) min_b = b;
			*/

			//LAB2RGB(L, a, b, R, G, B);
			//TRACE("R = %f, G = %f, B = %f\n", R, G, B);
			target[i][j].r = (GLubyte)(L * 255/100.);
			target[i][j].g = (GLubyte)((a + 100) * 255/200.);
			target[i][j].b = (GLubyte)((b + 100) * 255/200.);
		}
	}
	//TRACE("max_L = %f, max_a = %f, max_b = %f\n", max_L, max_a, max_b);
	//TRACE("min_L = %f, min_a = %f, min_b = %f\n", min_L, min_a, min_b);
}

void ConvertLABtoRGB(cimatrix& cmap, cimatrix& target)
{
	int i, j;
	int image_x, image_y;
	double R, G, B;
	double L, a, b;

	image_x = cmap.getRow();
	image_y = cmap.getCol();
	
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			L = (double)cmap[i][j].r;
			a = (double)cmap[i][j].g;
			b = (double)cmap[i][j].b;
			/////////////////////////////////
			L = L * 100/255.;
			a = (a * 200/255.) - 100.;
			b = (b * 200/255.) - 100.;
			/////////////////////////////
			//RGB2LAB(L, a, b, R, G, B);
			LAB2RGB(L, a, b, R, G, B);
			//TRACE("R = %f, G = %f, B = %f\n", R, G, B);
			target[i][j].r = (GLubyte)R;
			target[i][j].g = (GLubyte)G;
			target[i][j].b = (GLubyte)B;
		}
	}
}

void ConvertRGBtoLAB256(cimatrix& cmap, cimatrix& target)
{
	int i, j;
	int image_x, image_y;
	double R, G, B;
	double L, a, b;

	image_x = cmap.getRow();
	image_y = cmap.getCol();
	
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			R = (double)cmap[i][j].r;
			G = (double)cmap[i][j].g;
			B = (double)cmap[i][j].b;
			RGB2LAB256(R, G, B, L, a, b);
			//RGB2LAB(R, G, B, L, a, b);
			//TRACE("L = %f, a = %f, b = %f\n", L, a, b);
			//TRACE("r = %d, g = %d, b = %d\n", (int)R, (int)G, (int)B);
			//LAB2RGB(L, a, b, R, G, B);
			//TRACE("R = %f, G = %f, B = %f\n", R, G, B);
			target[i][j].r = (GLubyte)L;
			target[i][j].g = (GLubyte)a;
			target[i][j].b = (GLubyte)b;
			//TRACE("L = %d, a = %d, b = %d\n", target[i][j].r, target[i][j].g, target[i][j].b);
		}
	}
	//TRACE("max_L = %f, max_a = %f, max_b = %f\n", max_L, max_a, max_b);
	//TRACE("min_L = %f, min_a = %f, min_b = %f\n", min_L, min_a, min_b);
}

void ConvertLABtoRGB256(cimatrix& cmap, cimatrix& target)
{
	int i, j;
	int image_x, image_y;
	double R, G, B;
	double L, a, b;

	image_x = cmap.getRow();
	image_y = cmap.getCol();
	
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			L = (double)cmap[i][j].r;
			a = (double)cmap[i][j].g;
			b = (double)cmap[i][j].b;
			/////////////////////////////////
			/////////////////////////////
			LAB2RGB256(L, a, b, R, G, B);
			//TRACE("R = %f, G = %f, B = %f\n", R, G, B);
			target[i][j].r = (GLubyte)R;
			target[i][j].g = (GLubyte)G;
			target[i][j].b = (GLubyte)B;
			//TRACE("R = %d, G = %d, B = %d\n", target[i][j].r, target[i][j].g, target[i][j].b);
		}
	}
}

void ConvertRGBtoLAB_double(cimatrix& cmap, cmatrix& target)
// convert to double values
{
	int i, j;
	int image_x, image_y;
	double R, G, B;
	double L, a, b;

	image_x = cmap.getRow();
	image_y = cmap.getCol();
	
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			R = (double)cmap[i][j].r;
			G = (double)cmap[i][j].g;
			B = (double)cmap[i][j].b;
			//RGB2LAB256(R, G, B, L, a, b);
			RGB2LAB(R, G, B, L, a, b);
			//TRACE("L = %f, a = %f, b = %f\n", L, a, b);
			//TRACE("r = %d, g = %d, b = %d\n", (int)R, (int)G, (int)B);
			//LAB2RGB(L, a, b, R, G, B);
			//TRACE("R = %f, G = %f, B = %f\n", R, G, B);
			target[i][j].r = L;
			target[i][j].g = a;
			target[i][j].b = b;
			//TRACE("L = %f, a = %f, b = %f\n", target[i][j].r, target[i][j].g, target[i][j].b);
		}
	}
	//TRACE("max_L = %f, max_a = %f, max_b = %f\n", max_L, max_a, max_b);
	//TRACE("min_L = %f, min_a = %f, min_b = %f\n", min_L, min_a, min_b);
}

void ConvertLABtoRGB_double(cmatrix& cmap, cimatrix& target)
// convert from double to GLubyte
{
	int i, j;
	int image_x, image_y;
	double R, G, B;
	double L, a, b;

	image_x = cmap.getRow();
	image_y = cmap.getCol();
	
	////////////////////////////
	for (i = 0; i < image_x; i++) {
		for (j = 0; j < image_y; j++) {
			L = cmap[i][j].r;
			a = cmap[i][j].g;
			b = cmap[i][j].b;
			/////////////////////////////////
			/////////////////////////////
			LAB2RGB(L, a, b, R, G, B);
			//TRACE("R = %f, G = %f, B = %f\n", R, G, B);
			target[i][j].r = (GLubyte)R;
			target[i][j].g = (GLubyte)G;
			target[i][j].b = (GLubyte)B;
			//TRACE("R = %d, G = %d, B = %d\n", target[i][j].r, target[i][j].g, target[i][j].b);
		}
	}
}

void SaturatePalette(int N, double t)
{
	int i;
	double r, g, b;
	double h, s, v;
	//double t;

	delete [] pal_s;
	pal_s = new RGB[N];
	for (i = 0; i < N; i++) {
		r = pal[i].r / 255.;
		g = pal[i].g / 255.;
		b = pal[i].b / 255.;
		//TRACE("r = %f\n", r);
		//TRACE("g = %f\n", g);
		//TRACE("b = %f\n", b);
		RGB2HSV(r, g, b, h, s, v);
		///////////////////////////////
		// extrapolate! s and v
		//t = 1.5;
		s = (1-t)*0.0 + t*s;
		v = (1-t)*0.0 + t*v;
		if (s > 1.0) s = 1.0;
		if (v > 1.0) v = 1.0;
		///////////////////////////////
		HSV2RGB(h, s, v, r, g, b);
		pal_s[i].r = r*255;
		pal_s[i].g = g*255;
		pal_s[i].b = b*255;
		//TRACE("r2 = %f\n", r);
		//TRACE("g2 = %f\n", g);
		//TRACE("b2 = %f\n", b);
	}
}

void SaturateRGB(GLubyte& r1, GLubyte& g1, GLubyte& b1)
{
	double r, g, b;
	double h, s, v;
	double t;

	r = r1 / 255.;
	g = g1 / 255.;
	b = b1 / 255.;
	//TRACE("r = %f\n", r);
	//TRACE("g = %f\n", g);
	//TRACE("b = %f\n", b);
	RGB2HSV(r, g, b, h, s, v);
	///////////////////////////////
	// extrapolate! s and v
	t = 1.5;
	s = (1-t)*0.0 + t*s;
	v = (1-t)*0.0 + t*v;
	if (s > 1.0) s = 1.0;
	if (v > 1.0) v = 1.0;
	///////////////////////////////
	HSV2RGB(h, s, v, r, g, b);
	r1 = (GLubyte)(r*255);
	g1 = (GLubyte)(g*255);
	b1 = (GLubyte)(b*255);
	//TRACE("r2 = %f\n", r);
	//TRACE("g2 = %f\n", g);
	//TRACE("b2 = %f\n", b);
}



void ClearBuffer(GLubyte* buffer) 
{
	int x, y;

	for (y = 0; y < IMAGE_Y; y++) {
		for (x = 0; x < IMAGE_X; x++) {
			buffer[(y * IMAGE_X + x) * 3 + 0] = 255;
			buffer[(y * IMAGE_X + x) * 3 + 1] = 255;
			buffer[(y * IMAGE_X + x) * 3 + 2] = 255;
		}
	}
}

void ClearGrayBuffer(GLubyte* buffer) 
{
	int x, y;

	for (y = 0; y < IMAGE_Y; y++) {
		for (x = 0; x < IMAGE_X; x++) {
			buffer[y * IMAGE_X + x] = 255;
		}
	}
}

void BlackBuffer(GLubyte* buffer) 
{
	int x, y;

	for (y = 0; y < IMAGE_Y; y++) {
		for (x = 0; x < IMAGE_X; x++) {
			buffer[(y * IMAGE_X + x) * 3 + 0] = 0;
			buffer[(y * IMAGE_X + x) * 3 + 1] = 0;
			buffer[(y * IMAGE_X + x) * 3 + 2] = 0;
		}
	}
}


void MyDelay() 
{
	TRACE("MyDelay\n");
	for (int i = 0; i < DELAY; i++);
}


void CreateWhiteNoise(CDC& dc, int width, int height, double prob)
{
	int x, y;
	
	for (x = 0; x < width; x++)
		for (y = 0; y < height; y++) {
			if (drand48() > prob) 
				dc.SetPixelV(x, y, RGB(0, 0, 0));
		}
	
}

void DrawETF_LIC(CDC& dc, ETF& e, int length)
// this one is incorrect, because only integer locations are considered
// but result looks all right (strange!)
// I think this is because, this one is more faithful to the discrete nature of the image!
{
	//Node *q, *p;
	int k, x, y, half;
	//float tmp_cost;
	GLubyte	r, g, b;
	float theta, w, l;
	float w_ran, l_ran;
	int	fail_count;
	int	stroke_count = 0;
	//int	llc, rlc; // left, right length count
	double tx, ty;
	int	sum, count;
	double i, j;
	int	int_i, int_j;
	
	////////////////////////////////
	// Linked List of Strokes!
	//SimpleList<Stroke> list;

	theta = (float)(PI / 4.0); // 45 degree
	w_ran = 2; // random_width
	l_ran = 10; // random_length
	//w = 3;
	//w = 4; // 256x256
	w = 5; // 256x256
	//w = 8; // 128x128
	//w = 16; // 64x64
	l = 20; // 256x256
	//l = 30; // 128x128
	//l = 60; // 64x64
	//half = l; // actually, it should bigger than l/1.414
	half = length/2 + 1; // actually, it should bigger than l/1.414 / 2

	//Dbuffer: original
	//memDC: test and check compute difference
	//double_buffer: update and/or transfer back to memDC2
	///////////////////////////////////////
	/// for NSF proposal
	//LoadBMP(dc, (char *)LPCTSTR("Eagle_level2.bmp"), &memDC);
	///////////////////////////////////////////
	ClearMemDC(&memDC); // clear the canvas white
	//////////////////////////////////////////////////////
	ClearMemDC(&double_buffer); // clear the canvas white
	ClearMemDC(&memDC2); // clear the canvas white
	//ClearImage(IMAGE_X, IMAGE_Y, in_heap); // whether or not the pixel has been visited

	int image_x = e.getRow();
	int image_y = e.getCol();

	CreateWhiteNoise(double_buffer, image_x, image_y, 0.8);

	fail_count = 0;

	//CClientDC dc(this);

	
	srand( (unsigned)time( NULL ) );

	// first put enough strokes randomly
	//for (i = 0; i < IMAGE_X * IMAGE_Y / 20; i++) {
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//x = (int) ( (IMAGE_X-1) * (float)rand() / RAND_MAX );
			//y = (int) ( (IMAGE_Y-1) * (float)rand() / RAND_MAX );
			i = x; j = y; 
			sum = 0; 
			count = 1;
			sum += (int)RGB_GETRED(double_buffer.GetPixel(x, image_y-1-y));
			//tx = -gfield[x][y].gy;
			//ty = gfield[x][y].gx;
			tx = e[x][y].tx;
			ty = e[x][y].ty;
			//TRACE("sum = %d\n", sum);
			for (k = 0; k < half; k++) {
				if (tx == 0.0 && ty == 0.0) break;
				if (tx == 0.0) { // vertical
					j += 1;
					if (i < 0 || i > image_x-1 || j < 0 || j > image_y-1) break;
					int_i = round(i);
					//int_j = image_y-1-round(j);
					int_j = round(j);
					sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, image_y-1-int_j));
					//TRACE("sum = %d\n", sum);
					count++;
				}
				else if ( ABS(ty / tx) < 1 ) { 
					i += 1;
					j = j + ty/tx;
					if (i < 0 || i > image_x-1 || j < 0 || j > image_y-1) break;
					int_i = round(i);
					//int_j = image_y-1-round(j);
					//sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, int_j));
					int_j = round(j);
					sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, image_y-1-int_j));
					//TRACE("sum = %d\n", sum);
					count++;
				}
				else { // ABS(ty / tx) >= 1
					j += 1;
					i = i + tx/ty;
					if (i < 0 || i > image_x-1 || j < 0 || j > image_y-1) break;
					int_i = round(i);
					//int_j = image_y-1-round(j);
					//sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, int_j));
					int_j = round(j);
					sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, image_y-1-int_j));
					//TRACE("sum = %d\n", sum);
					count++;
				}
				tx = e[int_i][int_j].tx;
				ty = e[int_i][int_j].ty;
			}
			i = x; j = y; 
			tx = e[x][y].tx;
			ty = e[x][y].ty;
			//////////////////////////////////////////
			// Other half
			///*
			for (k = 0; k < half; k++) {
				if (tx == 0.0 && ty == 0.0) break;
				if (tx == 0.0) { // vertical
					j -= 1;
					if (i < 0 || i > image_x-1 || j < 0 || j > image_y-1) break;
					int_i = round(i);
					//int_j = image_y-1-round(j);
					int_j = round(j);
					sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, image_y-1-int_j));
					//TRACE("sum = %d\n", sum);
					count++;
				}
				else if ( ABS(ty / tx) < 1 ) { 
					i -= 1;
					j = j - ty/tx;
					if (i < 0 || i > image_x-1 || j < 0 || j > image_y-1) break;
					int_i = round(i);
					//int_j = image_y-1-round(j);
					//sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, int_j));
					int_j = round(j);
					sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, image_y-1-int_j));
					//TRACE("sum = %d\n", sum);
					count++;
				}
				else { // ABS(ty / tx) >= 1
					j -= 1;
					i = i - tx/ty;
					if (i < 0 || i > image_x-1 || j < 0 || j > image_y-1) break;
					int_i = round(i);
					//int_j = image_y-1-round(j);
					//sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, int_j));
					int_j = round(j);
					sum += (int)RGB_GETRED(double_buffer.GetPixel(int_i, image_y-1-int_j));
					//TRACE("sum = %d\n", sum);
					count++;
				}
				tx = e[int_i][int_j].tx;
				ty = e[int_i][int_j].ty;
			}
			//*/
			//TRACE("sum_total = %d\n", sum);
			sum = sum / count;
			//TRACE("sum = %d, count = %d\n", sum, count);

			r = g = b = (GLubyte)sum;
			//TRACE("r = %d, g = %d, b = %d\n", r, g, b);

			//PerturbRGB(r, g, b, 200);
				
			//pen.CreatePen(PS_SOLID, (int)w, RGB(0, 255, 255));
			//pen.CreatePen(PS_SOLID, (int)w, RGB(r+perturb_r, g+perturb_g, b+perturb_b));
			memDC.SetPixelV(x, image_y-1-y, RGB(r, g, b));
		
			//dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
			///////////////////////////////////////////////////////////////
			// Copy randomly created strokes in memDC2
			//memDC2.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
		}
	}
	dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
	//double_buffer.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY); // copy memDC to double_buffer
	
	//////////////////////////////////////////////////////////////
	// Copy double_buffer to Dbuffer
	/*
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r = (GLubyte)RGB_GETRED(double_buffer.GetPixel(x, y));
			g = (GLubyte)RGB_GETGREEN(double_buffer.GetPixel(x, y));
			b = (GLubyte)RGB_GETBLUE(double_buffer.GetPixel(x, y));
			Dbuffer[((image_y-1-y) * image_x + x) * 3 + 0] = r;
			Dbuffer[((image_y-1-y) * image_x + x) * 3 + 1] = g;
			Dbuffer[((image_y-1-y) * image_x + x) * 3 + 2] = b;
		}
	}
	*/
}

void DrawETFColor360(CDC& dc, ETF& e)
// Display vector as a color, to show how smooth the field is
{
	int x, y;
	GLubyte	r, g, b;
	double tx, ty;
	double h, s, v, r2, g2, b2;

	int image_x = e.getRow();
	int image_y = e.getCol();

	v = s = 1.0;
	
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {

			tx = e[x][y].tx;
			ty = e[x][y].ty;

			h = atan2(ty, tx); // [-PI, PI]
			if (h >= 0.0)
				h = (h / PI) * 180; // [0, 180]
			else 
				h = 360 + (h / PI) * 180; // [180, 360]

			if (h > 360) h = 360;
			if (h < 0) h = 0;

			///////////////////////////////
			HSV2RGB(h, s, v, r2, g2, b2);
			///////////////////////////////
			
			r = (GLubyte)round(r2 * 255);
			g = (GLubyte)round(g2 * 255);
			b = (GLubyte)round(b2 * 255);
			//TRACE("r = %d, g = %d, b = %d\n", r, g, b);

			dc.SetPixelV(x, image_y-1-y, RGB(r, g, b));
		
		}
	}
}

void DrawETFColorTensor(CDC& dc, ETF& e)
// Display vector as a color, to show how smooth the field is
// tensor version: 90 <= h <= 270
{
	int x, y;
	GLubyte	r, g, b;
	double tx, ty;
	double h, s, v, r2, g2, b2;

	int image_x = e.getRow();
	int image_y = e.getCol();

	v = s = 1.0;
	
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {

			tx = e[x][y].tx;
			ty = e[x][y].ty;

			h = atan2(ty, tx); // [-PI, PI]
			if (h >= 0.0)
				h = (h / PI) * 180; // [0, 180]
			else 
				h = 360 + (h / PI) * 180; // [180, 360]

			//if (h > 360) h = 360;
			//if (h < 0) h = 0;

			if (h > 270) h = 270;
			if (h < 90) h = 90;

			////////////////////////////////
			// Out of the tensor range!
			/*
            if (h > 270 || h < 90) {
				if (h != 0.0) // 0 is an exception
					TRACE("h = %f, [%f, %f]\n", h, tx, ty);
			}
			*/
			///////////////////////////////

			//////////////////////////////////
			// Readjust the color range!
			//if (h != 0.0) h -= 90; // [0, 180]
			h -= 90; // [0, 180]
			h *= 2; // [0, 360]
			if (h > 360) h = 360;
			if (h < 0) h = 0;
			///////////////////////////////
			HSV2RGB(h, s, v, r2, g2, b2);
			///////////////////////////////
			
			r = (GLubyte)round(r2 * 255);
			g = (GLubyte)round(g2 * 255);
			b = (GLubyte)round(b2 * 255);
			//TRACE("r = %d, g = %d, b = %d\n", r, g, b);

			dc.SetPixelV(x, image_y-1-y, RGB(r, g, b));
		
		}
	}
}

void DrawETFColor180(CDC& dc, ETF& e)
// Display vector as a color, to show how smooth the field is
// 180 degree version: 0 <= h <= 180
{
	int x, y;
	GLubyte	r, g, b;
	double tx, ty;
	double h, s, v, r2, g2, b2;

	int image_x = e.getRow();
	int image_y = e.getCol();

	v = s = 1.0;
	
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {

			tx = e[x][y].tx;
			ty = e[x][y].ty;

			h = atan2(ty, tx); // [-PI, PI]
			if (h >= 0.0)
				h = (h / PI) * 180; // [0, 180]
			else 
				h = 360 + (h / PI) * 180; // [180, 360]

			//if (h > 360) h = 360;
			//if (h < 0) h = 0;

			if (h > 180) h = 180;
			if (h < 0) h = 0;

			///////////////////////////////
			//////////////////////////////////
			// Readjust the color range!
			h *= 2; // [0, 360]
			///////////////////////////////
			HSV2RGB(h, s, v, r2, g2, b2);
			///////////////////////////////
			
			r = (GLubyte)round(r2 * 255);
			g = (GLubyte)round(g2 * 255);
			b = (GLubyte)round(b2 * 255);
			//TRACE("r = %d, g = %d, b = %d\n", r, g, b);

			dc.SetPixelV(x, image_y-1-y, RGB(r, g, b));
		
		}
	}
}



#define gauss1D_bl_macro(x, y) ( exp( (-(x)*(x)) / (2*(y)*(y)) ) )

void GetFBL(cimatrix& image, ETF& e, double sigma, double sigma2, 
						  double sigma3, double sigma4, int max_itr)
// flow-based color bilateral filter
// Bilateral + bilateral (alternating)
// Each bilateral filter has a different range sigma values! (sigma3 and sigma4)
// FAST & Separate version, horizontal then vertical, same number of iterations!
{
	int k, m;
	GLubyte	r, g, b;
	//double tx, ty;
	int i, j, i_x, i_y;

	vector vn(2), vt(2), w(2);
	double x, y, d_x, d_y;
	int s;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double bb, tt, val, uu, vv;
	//double sum, c_val;
	//double weight, w_sum;
	int d1;

	double sum_r, sum_g, sum_b;
	//double weight_r, weight_g, weight_b;
	//double w_sum_r, w_sum_g, w_sum_b;
	double weight, w_sum;
	double c_val_r, c_val_g, c_val_b;
	double val_r, val_g, val_b;
	//double angle, angle1, angle2;
	double c_dist;


	int image_x = image.getRow();
	int image_y = image.getCol();

	int half_l, half_w;

	//half_l = length/2 + 1; // actually, it should bigger than l/1.414 / 2
	//half_l = length/2;
		
	cimatrix tmp;
	tmp.copy(image);

	//double sigma;
	//double sigma2;
	
	//sigma = 1.0;
	//sigma2 = sigma * 1.0;
	//sigma2 = sigma * 10.0;
	//sigma2 = sigma * 2.0;
	//sigma3 = sigma;

	vector GAU1, GAU2, GAU3, GAU4;
	MakeGaussVector(sigma, GAU1); // length of the kernel
	MakeGaussVector(sigma2, GAU2); // width of the kernel
	//MakeGaussVector(sigma3, GAU3); // for similarity function

	half_l = GAU1.getMax()-1;
	half_w = GAU2.getMax()-1;

	MakeGaussVectorBL(sigma3, GAU3, 300); // 300 array elements created
	MakeGaussVectorBL(sigma4, GAU4, 300); // 300 array elements created
	
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;
	double step_size;
	double step_size1, step_size2;

	//step_size = 1.0;
	//step_size = 2.0;
	step_size = 1.0;
	step_size1 = step_size2 = 1.0;

	//half_l = half_w = 25; // for TVCG
	
	for (m = 0; m < max_itr; m++) {
		/////////////////////////////////////////
		// FBL along main axis!
		//StartTimer();
		for (j = 0; j < IMAGE_Y; j++) {
			for (i = 0; i < IMAGE_X; i++) {
				//sum = 0.0;
				sum_r = sum_g = sum_b = 0.0;
				weight = 0.0;
				w_sum = 0.0;
				///////////////////////////////////////////
				// First add the center intensity
				//weight = 1.0;
				//c_val = (double)tmp[i][j];
				//sum += weight * c_val;
				//w_sum += weight;
				////////////////////////////////////////////////
				// One half
				d_x = (double)i; d_y = (double)j; 
				i_x = i;	i_y = j;
				////////////////////////////////////////////////
				// Use the center color at each minor axis!
				c_val_r = (double)tmp[i_x][i_y].r;
				c_val_g = (double)tmp[i_x][i_y].g;
				c_val_b = (double)tmp[i_x][i_y].b;
				//////////////////////
				for (k = 0; k < half_l; k++) {
					vt[0] = e[i_x][i_y].tx;
					vt[1] = e[i_x][i_y].ty;
					if (vt[0] == 0.0 && vt[1] == 0.0) break;
					//vt.make_unit();
					///////////////////////////
					//d_x = c_x + vt[0] * t * step_size; // how narrow the height of rectangle kernel
					//d_y = c_y + vt[1] * t * step_size; // depends on step_size
					////////////////////////////////////////////////
					////////////////////////////////////////////////
					x = d_x;
					y = d_y;
					/////////////////////////////////////////////////////
					if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
						continue;
					/////////////////////////////////////////////////
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					//val = (double)tmp[x1][y1];
					val_r = (double)tmp[x1][y1].r;
					val_g = (double)tmp[x1][y1].g;
					val_b = (double)tmp[x1][y1].b;
					//val = (val_r + val_g + val_b) / 3.0;
					/////////////////////////////////////////////////////////
					//weight_r = weight_g = weight_b = GAU1[k];
					//weight = GAU2[d1];
					//weight = 0.0;
					//if (k == 0) weight = 0.0;
					//else 
					weight = GAU1[k]; // length
					//d3 = round( ABS(val - c_val) );
					c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
					//weight *= gauss1D_bl_macro( c_dist/1.732, sigma3); // like Gaussian!
					//weight *= gauss1D_norm2( c_dist/1.732, 0.0, sigma3); // like Gaussian!
					//weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma3); // like Gaussian!
					//TRACE("c_dist/1.732 = %f\n", c_dist/1.732);
					weight *= GAU3[ (int)(c_dist/1.732) ]; // like Gaussian!
					/////////////////////////////////////////////////////////
					/////////////////////////////////////////////////////
					//sum1 += image[xi][yi] * weight;
					sum_r += val_r * weight;
					sum_g += val_g * weight;
					sum_b += val_b * weight;
					//sum1 += image[xi][yi] * weight;
					//sum += val * weight;
					w_sum += weight;
					/////////////////////////////////////////
					d_x += vt[0] * step_size1; // accurate new location x
					d_y += vt[1] * step_size1; // accurate new location y
					/////////////////////////////////////////
					i_x = round(d_x); // integer version of new location x
					i_y = round(d_y); // integer version of new location y
					//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
					if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				}
				////////////////////////////////////////////////
				// Other half
				d_x = (double)i; d_y = (double)j; 
				i_x = i;	i_y = j;
				for (k = 0; k < half_l; k++) {
					vt[0] = -e[i_x][i_y].tx;
					vt[1] = -e[i_x][i_y].ty;
					if (vt[0] == 0.0 && vt[1] == 0.0) break;
					//vt.make_unit();
					///////////////////////////
					///////////////////////////////
					x = d_x;
					y = d_y;
					/////////////////////////////////////////////////////
					if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
						continue;
					/////////////////////////////////////////////////
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					//val = (double)tmp[x1][y1];
					val_r = (double)tmp[x1][y1].r;
					val_g = (double)tmp[x1][y1].g;
					val_b = (double)tmp[x1][y1].b;
					//val = (val_r + val_g + val_b) / 3.0;
					/////////////////////////////////////////////////////////
					//weight_r = weight_g = weight_b = GAU1[k];
					//if (k == 0) weight = 0.0;
					//else 
					weight = GAU1[k];
					//d3 = round( ABS(val - c_val) );
					c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
					//weight *= gauss1D_bl_macro( c_dist/1.732, sigma3); // like Gaussian!
					//weight *= gauss1D_norm2( c_dist/1.732, 0.0, sigma3); // like Gaussian!
					//weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma3); // like Gaussian!
					weight *= GAU3[ (int)(c_dist/1.732) ]; // like Gaussian!
					/////////////////////////////////////////////////////////
					/////////////////////////////////////////////////////
					//sum1 += image[xi][yi] * weight;
					sum_r += val_r * weight;
					sum_g += val_g * weight;
					sum_b += val_b * weight;
					//sum1 += image[xi][yi] * weight;
					//sum += val * weight;
					w_sum += weight;
					/////////////////////////////////////////
					d_x += vt[0] * step_size1; // accurate new location x
					d_y += vt[1] * step_size1; // accurate new location y
					/////////////////////////////////////////
					i_x = round(d_x); // integer version of new location x
					i_y = round(d_y); // integer version of new location y
					//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
					if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
				}
				//TRACE("sum_total = %d\n", sum);
				//sum = (int) ( sum / (double)count );
				//sum /= (double)count;
				sum_r /= (double)w_sum; // normalize!
				sum_g /= (double)w_sum; // normalize!
				sum_b /= (double)w_sum; // normalize!
				//TRACE("sum = %d, count = %d\n", sum, count);
				/////////////////////////////////////
				// Didn't add a single color!
				//if (sum_r == 0.0 && sum_g == 0.0 && sum_b == 0.0) {
				if (e[i][j].tx == 0.0 && e[i][j].ty == 0.0) {
					sum_r = tmp[i][j].r;
					sum_g = tmp[i][j].g;
					sum_b = tmp[i][j].b;
				}
				/////////////////////////////////////
				//if (sum_r > 255.0) sum_r = 255.0;
				//if (sum_g > 255.0) sum_g = 255.0;
				//if (sum_b > 255.0) sum_b = 255.0;
				//if (sum_r < 0.0) sum_r = 0.0;
				//if (sum_g < 0.0) sum_g = 0.0;
				//if (sum_b < 0.0) sum_b = 0.0;
				//sum_i = (GLubyte)sum;
				//r = g = b = (GLubyte)sum;
				r = (GLubyte)sum_r;
				g = (GLubyte)sum_g;
				b = (GLubyte)sum_b;
				//TRACE("r = %d, g = %d, b = %d\n", r, g, b);
				//memDC.SetPixelV(x, IMAGE_Y-1-y, RGB(r, g, b));
				image[i][j].r = r;
				image[i][j].g = g;
				image[i][j].b = b;
				//
				///////////////////////////////////////////////////////////////
				// Copy randomly created strokes in memDC2
				//memDC2.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
			}
		}
		tmp.copy(image); // update tmp

		/////////////////////////////////////////
		// FBL along minor axis!
		for (j = 0; j < IMAGE_Y; j++) {
			for (i = 0; i < IMAGE_X; i++) {
				//sum = 0.0;
				sum_r = sum_g = sum_b = 0.0;
				weight = 0.0;
				w_sum = 0.0;
				///////////////////////////////////////////
				// First add the center intensity
				//weight = 1.0;
				//c_val = (double)tmp[i][j];
				//sum += weight * c_val;
				//w_sum += weight;
				////////////////////////////////////////////////
				// One half
				d_x = (double)i; d_y = (double)j; 
				i_x = i;	i_y = j;
				//c_val = (double)tmp[i][j];
				//c_val_r = (double)tmp[i][j].r;
				//c_val_g = (double)tmp[i][j].g;
				//c_val_b = (double)tmp[i][j].b;
				vn[0] = -e[i_x][i_y].ty;
				vn[1] = e[i_x][i_y].tx;
				if (vn[0] == 0.0 && vn[1] == 0.0) break;
				//vn.make_unit();
				///////////////////////////
				//d_x = c_x + vt[0] * t * step_size; // how narrow the height of rectangle kernel
				//d_y = c_y + vt[1] * t * step_size; // depends on step_size
				////////////////////////////////////////////////
				// Use the center color at each minor axis!
				c_val_r = (double)tmp[i_x][i_y].r;
				c_val_g = (double)tmp[i_x][i_y].g;
				c_val_b = (double)tmp[i_x][i_y].b;
				////////////////////////////////////////////////
				for (s = -half_w; s <= half_w; s++) { // width of Gaussian kernel
					////////////////////////
					x = d_x + vn[0] * s;
					y = d_y + vn[1] * s;
					/////////////////////////////////////////////////////
					if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
						continue;
					/////////////////////////////////////////////////
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					//val = (double)tmp[x1][y1];
					val_r = (double)tmp[x1][y1].r;
					val_g = (double)tmp[x1][y1].g;
					val_b = (double)tmp[x1][y1].b;
					//val = (val_r + val_g + val_b) / 3.0;
					/////////////////////////////////////////////////////////
					//weight_r = weight_g = weight_b = GAU1[k];
					d1 = ABS(s);
					//weight = GAU2[d1];
					weight = GAU2[d1];
					//d3 = round( ABS(val - c_val) );
					c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
					//weight *= gauss1D_norm2( c_dist/1.732, 0.0, sigma4);
					//weight *= gauss1D_bl_macro( c_dist/1.732, sigma4);
					//weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma4);
					weight *= GAU4[ (int)(c_dist/1.732) ]; // like Gaussian!
					/////////////////////////////////////////////////////
					//sum1 += image[xi][yi] * weight;
					sum_r += val_r * weight;
					sum_g += val_g * weight;
					sum_b += val_b * weight;
					//sum1 += image[xi][yi] * weight;
					//sum += val * weight;
					w_sum += weight;
				}
				/////////////////////////////////////////
				//TRACE("sum_total = %d\n", sum);
				//sum = (int) ( sum / (double)count );
				//sum /= (double)count;
				sum_r /= (double)w_sum; // normalize!
				sum_g /= (double)w_sum; // normalize!
				sum_b /= (double)w_sum; // normalize!
				//TRACE("sum = %d, count = %d\n", sum, count);
				/////////////////////////////////////
				// Didn't add a single color!
				//if (sum_r == 0.0 && sum_g == 0.0 && sum_b == 0.0) {
				if (e[i][j].tx == 0.0 && e[i][j].ty == 0.0) {
					sum_r = tmp[i][j].r;
					sum_g = tmp[i][j].g;
					sum_b = tmp[i][j].b;
				}
				/////////////////////////////////////
				//if (sum_r > 255.0) sum_r = 255.0;
				//if (sum_g > 255.0) sum_g = 255.0;
				//if (sum_b > 255.0) sum_b = 255.0;
				//if (sum_r < 0.0) sum_r = 0.0;
				//if (sum_g < 0.0) sum_g = 0.0;
				//if (sum_b < 0.0) sum_b = 0.0;
				//sum_i = (GLubyte)sum;
				//r = g = b = (GLubyte)sum;
				r = (GLubyte)sum_r;
				g = (GLubyte)sum_g;
				b = (GLubyte)sum_b;
				//TRACE("r = %d, g = %d, b = %d\n", r, g, b);
				//memDC.SetPixelV(x, IMAGE_Y-1-y, RGB(r, g, b));
				image[i][j].r = r;
				image[i][j].g = g;
				image[i][j].b = b;
				//
				///////////////////////////////////////////////////////////////
				// Copy randomly created strokes in memDC2
				//memDC2.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
			}
		}
		//EndTimer("_Result_FBL.txt");
		tmp.copy(image); // update tmp
	}
}

void GetSABL(CDC& dc, cimatrix& image, ETF& e, double sigma, double sigma2, double sigma3, double sigma4, int max_itr)
// structure_adaptive bilateral filter (Tuan Pham)
// use parabolic arc!
{
	int k, m;
	GLubyte	rr, gg, bb;
	//double tx, ty;
	int i, j, i_x, i_y;

	vector vn(2), vt(2), w(2);
	double x, y, d_x, d_y, c_x, c_y;
	int s;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double bb, tt, val, uu, vv;
	//double sum, c_val;
	//double weight, w_sum;
	int d1;

	double sum_r, sum_g, sum_b;
	//double weight_r, weight_g, weight_b;
	//double w_sum_r, w_sum_g, w_sum_b;
	double weight, w_sum;
	double c_val_r, c_val_g, c_val_b;
	double val_r, val_g, val_b;
	//double angle, angle1, angle2;
	double c_dist;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half_l, half_w;

	//half_l = length/2 + 1; // actually, it should bigger than l/1.414 / 2
	//half_l = length/2;
		
	cimatrix tmp;
	tmp.copy(image);

	//double sigma;
	//double sigma2;
	
	//sigma = 1.0;
	//sigma2 = sigma * 1.0;
	//sigma2 = sigma * 10.0;
	//sigma2 = sigma * 2.0;
	//sigma3 = sigma;

	vector GAU1, GAU2, GAU3, GAU4;
	MakeGaussVector(sigma, GAU1); // length of the kernel
	MakeGaussVector(sigma2, GAU2); // width of the kernel
	//MakeGaussVector(sigma3, GAU3); // for similarity function

	half_l = GAU1.getMax()-1;
	TRACE("half_l = %d\n", half_l);
	half_w = GAU2.getMax()-1;

	MakeGaussVectorBL(sigma3, GAU3, 300); // 300 array elements created
	MakeGaussVectorBL(sigma4, GAU4, 300); // 300 array elements created
	
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;
	double step_size;
	double step_size1, step_size2;

	//step_size = 1.0;
	//step_size = 2.0;
	step_size = 1.0;
	step_size1 = step_size2 = 1.0;

	//half_l = half_w = 25; // for TVCG

	vector a(2), b(2), c(2);
	double cos_theta, theta, cos_angle, angle, cos_val, sin_val, alpha;

	bool random_ok = false;

	CPen pen, *pOldPen;

	pen.CreatePen(PS_SOLID, 3, RGB(255, 0, 0));
	pOldPen = (CPen *)dc.SelectObject(&pen);

	for (m = 0; m < max_itr; m++) {
		/////////////////////////////////////////
		// FBL along main axis!
		//StartTimer();
		for (j = 0; j < IMAGE_Y; j++) {
			for (i = 0; i < IMAGE_X; i++) {
				//sum = 0.0;
				
				if (e[i][j].tx == 0.0 && e[i][j].ty == 0.0) {
					sum_r = tmp[i][j].r;
					sum_g = tmp[i][j].g;
					sum_b = tmp[i][j].b;

					goto here; // jump to here!
				}

				sum_r = sum_g = sum_b = 0.0;
				weight = 0.0;
				w_sum = 0.0;

				///////////////////////////////////////////
				// compute Parabola parameters
				a[0] = e[i][j].tx;
				a[1] = e[i][j].ty;

				d_x = (double)i; d_y = (double)j; 
				d_x += a[0]; // accurate new location x
				d_y += a[1]; // accurate new location y
				if (d_x < 0) d_x = 0; if (d_x > IMAGE_X-1) d_x = IMAGE_X-1;
				if (d_y < 0) d_y = 0; if (d_y > IMAGE_Y-1) d_y = IMAGE_Y-1;
				a[0] = e[round(d_x)][round(d_y)].tx;
				a[1] = e[round(d_x)][round(d_y)].ty;
				if(a[0] == 0.0 && a[1] == 0.0) {
					sum_r = tmp[i][j].r;
					sum_g = tmp[i][j].g;
					sum_b = tmp[i][j].b;

					goto here; // jump to here!
				}
				d_x += a[0]; // accurate new location x
				d_y += a[1]; // accurate new location y
				if (d_x < 0) d_x = 0; if (d_x > IMAGE_X-1) d_x = IMAGE_X-1;
				if (d_y < 0) d_y = 0; if (d_y > IMAGE_Y-1) d_y = IMAGE_Y-1;
				a[0] = d_x - (double)i;
				a[1] = d_y - (double)j;
				/////////////////////////////////////////
				b[0] = -e[i][j].tx;
				b[1] = -e[i][j].ty;
				d_x = (double)i; d_y = (double)j; 
				d_x += b[0]; // accurate new location x
				d_y += b[1]; // accurate new location y
				if (d_x < 0) d_x = 0; if (d_x > IMAGE_X-1) d_x = IMAGE_X-1;
				if (d_y < 0) d_y = 0; if (d_y > IMAGE_Y-1) d_y = IMAGE_Y-1;
				b[0] = -e[round(d_x)][round(d_y)].tx;
				b[1] = -e[round(d_x)][round(d_y)].ty;
				if(b[0] == 0.0 && b[1] == 0.0) {
					sum_r = tmp[i][j].r;
					sum_g = tmp[i][j].g;
					sum_b = tmp[i][j].b;

					goto here; // jump to here!
				}
				d_x += b[0]; // accurate new location x
				d_y += b[1]; // accurate new location y
				if (d_x < 0) d_x = 0; if (d_x > IMAGE_X-1) d_x = IMAGE_X-1;
				if (d_y < 0) d_y = 0; if (d_y > IMAGE_Y-1) d_y = IMAGE_Y-1;
				b[0] = d_x - (double)i;
				b[1] = d_y - (double)j;
				/////////////////////////////////////////////////
				a.make_unit();
				b.make_unit();
				cos_theta = a[0] * b[0] + a[1] * b[1];
				theta = acos(cos_theta); // [0, pi]
				alpha = sin( (PI - theta) / 2 ); // parameter for parabola: y = alpha * x^2
				c[0] = a[0] + b[0];
				c[1] = a[1] + b[1];
				c.make_unit();
				cos_angle = c[0] * 0 + c[1] * 1; // angle formed with y axis vector (0, 1)
				angle = acos(cos_angle);
				if (c[0] >= 0.0) angle = -angle; // Clockwise rotation
				cos_val = cos(angle);
				sin_val = sin(angle);
				///////////////////////////////////////////////////////////
				
				////////////////////////////////////////////////
				c_x = (double)i; c_y = (double)j; 
				i_x = i;	i_y = j;
				////////////////////////////////////////////////
				// Use the center color at each minor axis!
				c_val_r = (double)tmp[i_x][i_y].r;
				c_val_g = (double)tmp[i_x][i_y].g;
				c_val_b = (double)tmp[i_x][i_y].b;
				//////////////////////

				//random_ok = (drand48() < 0.01) ? true : false;

				//////////////////////////////////////
				for (k = -half_l; k <= half_l; k++) {
					////////////////////////////////////////////////
					d_x = k;
					d_y = alpha * d_x * d_x;
					/////////////////////////////////
					x = c_x + ( cos_val * d_x - sin_val * d_y );
					y = c_y + ( sin_val * d_x + cos_val * d_y );
					/////////////////////////////////////////////////////
					//if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
					//	continue;
					/////////////////////////////////////////////////
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					/////////////////////////////////////////////////////////
					/*
					if (random_ok && k == -half_l) {
						//TRACE("MoveTo!\n");
						dc.MoveTo(x1, IMAGE_Y-1-y1);
					}
					//TRACE("k = %d\n", k);
					if (random_ok && k > -half_l) {
						//TRACE("LineTo!\n");
						dc.LineTo(x1, IMAGE_Y-1-y1);
					}
					*/
					//////////////////////////////////////////////////
					//val = (double)tmp[x1][y1];
					val_r = (double)tmp[x1][y1].r;
					val_g = (double)tmp[x1][y1].g;
					val_b = (double)tmp[x1][y1].b;
					//val = (val_r + val_g + val_b) / 3.0;
					/////////////////////////////////////////////////////////
					//weight_r = weight_g = weight_b = GAU1[k];
					//weight = GAU2[d1];
					//weight = 0.0;
					//if (k == 0) weight = 0.0;
					//else 
					weight = GAU1[ABS(k)]; // length
					//d3 = round( ABS(val - c_val) );
					c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
					//weight *= gauss1D_bl_macro( c_dist/1.732, sigma3); // like Gaussian!
					//weight *= gauss1D_norm2( c_dist/1.732, 0.0, sigma3); // like Gaussian!
					//weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma3); // like Gaussian!
					//TRACE("c_dist/1.732 = %f\n", c_dist/1.732);
					weight *= GAU3[ (int)(c_dist/1.732) ]; // like Gaussian!
					/////////////////////////////////////////////////////////
					/////////////////////////////////////////////////////
					//sum1 += image[xi][yi] * weight;
					sum_r += val_r * weight;
					sum_g += val_g * weight;
					sum_b += val_b * weight;
					//sum1 += image[xi][yi] * weight;
					//sum += val * weight;
					w_sum += weight;
					/////////////////////////////////////////
				}
				
			here: // jump to here!
				//TRACE("sum_total = %d\n", sum);
				//sum = (int) ( sum / (double)count );
				//sum /= (double)count;
				sum_r /= (double)w_sum; // normalize!
				sum_g /= (double)w_sum; // normalize!
				sum_b /= (double)w_sum; // normalize!
				//TRACE("sum = %d, count = %d\n", sum, count);
				/////////////////////////////////////
				// Didn't add a single color!
				/*
				if (e[i][j].tx == 0.0 && e[i][j].ty == 0.0) {
					sum_r = tmp[i][j].r;
					sum_g = tmp[i][j].g;
					sum_b = tmp[i][j].b;
				}
				*/
				rr = (GLubyte)sum_r;
				gg = (GLubyte)sum_g;
				bb = (GLubyte)sum_b;
				//TRACE("r = %d, g = %d, b = %d\n", r, g, b);
				//memDC.SetPixelV(x, IMAGE_Y-1-y, RGB(r, g, b));
				image[i][j].r = rr;
				image[i][j].g = gg;
				image[i][j].b = bb;
				//
				///////////////////////////////////////////////////////////////
				// Copy randomly created strokes in memDC2
				//memDC2.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
			}
		}
		tmp.copy(image); // update tmp

		/////////////////////////////////////////
		// FBL along minor axis!
		for (j = 0; j < IMAGE_Y; j++) {
			for (i = 0; i < IMAGE_X; i++) {
				//sum = 0.0;
				sum_r = sum_g = sum_b = 0.0;
				weight = 0.0;
				w_sum = 0.0;
				///////////////////////////////////////////
				// First add the center intensity
				//weight = 1.0;
				//c_val = (double)tmp[i][j];
				//sum += weight * c_val;
				//w_sum += weight;
				////////////////////////////////////////////////
				// One half
				d_x = (double)i; d_y = (double)j; 
				i_x = i;	i_y = j;
				//c_val = (double)tmp[i][j];
				//c_val_r = (double)tmp[i][j].r;
				//c_val_g = (double)tmp[i][j].g;
				//c_val_b = (double)tmp[i][j].b;
				vn[0] = -e[i_x][i_y].ty;
				vn[1] = e[i_x][i_y].tx;
				if (vn[0] == 0.0 && vn[1] == 0.0) break;
				//vn.make_unit();
				///////////////////////////
				//d_x = c_x + vt[0] * t * step_size; // how narrow the height of rectangle kernel
				//d_y = c_y + vt[1] * t * step_size; // depends on step_size
				////////////////////////////////////////////////
				// Use the center color at each minor axis!
				c_val_r = (double)tmp[i_x][i_y].r;
				c_val_g = (double)tmp[i_x][i_y].g;
				c_val_b = (double)tmp[i_x][i_y].b;
				////////////////////////////////////////////////
				for (s = -half_w; s <= half_w; s++) { // width of Gaussian kernel
					////////////////////////
					x = d_x + vn[0] * s;
					y = d_y + vn[1] * s;
					/////////////////////////////////////////////////////
					if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
						continue;
					/////////////////////////////////////////////////
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					//val = (double)tmp[x1][y1];
					val_r = (double)tmp[x1][y1].r;
					val_g = (double)tmp[x1][y1].g;
					val_b = (double)tmp[x1][y1].b;
					//val = (val_r + val_g + val_b) / 3.0;
					/////////////////////////////////////////////////////////
					//weight_r = weight_g = weight_b = GAU1[k];
					d1 = ABS(s);
					//weight = GAU2[d1];
					weight = GAU2[d1];
					//d3 = round( ABS(val - c_val) );
					c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
					//weight *= gauss1D_norm2( c_dist/1.732, 0.0, sigma4);
					//weight *= gauss1D_bl_macro( c_dist/1.732, sigma4);
					//weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma4);
					weight *= GAU4[ (int)(c_dist/1.732) ]; // like Gaussian!
					/////////////////////////////////////////////////////
					//sum1 += image[xi][yi] * weight;
					sum_r += val_r * weight;
					sum_g += val_g * weight;
					sum_b += val_b * weight;
					//sum1 += image[xi][yi] * weight;
					//sum += val * weight;
					w_sum += weight;
				}
				/////////////////////////////////////////
				//TRACE("sum_total = %d\n", sum);
				//sum = (int) ( sum / (double)count );
				//sum /= (double)count;
				sum_r /= (double)w_sum; // normalize!
				sum_g /= (double)w_sum; // normalize!
				sum_b /= (double)w_sum; // normalize!
				//TRACE("sum = %d, count = %d\n", sum, count);
				/////////////////////////////////////
				// Didn't add a single color!
				//if (sum_r == 0.0 && sum_g == 0.0 && sum_b == 0.0) {
				if (e[i][j].tx == 0.0 && e[i][j].ty == 0.0) {
					sum_r = tmp[i][j].r;
					sum_g = tmp[i][j].g;
					sum_b = tmp[i][j].b;
				}
				/////////////////////////////////////
				//if (sum_r > 255.0) sum_r = 255.0;
				//if (sum_g > 255.0) sum_g = 255.0;
				//if (sum_b > 255.0) sum_b = 255.0;
				//if (sum_r < 0.0) sum_r = 0.0;
				//if (sum_g < 0.0) sum_g = 0.0;
				//if (sum_b < 0.0) sum_b = 0.0;
				//sum_i = (GLubyte)sum;
				//r = g = b = (GLubyte)sum;
				rr = (GLubyte)sum_r;
				gg = (GLubyte)sum_g;
				bb = (GLubyte)sum_b;
				//TRACE("r = %d, g = %d, b = %d\n", r, g, b);
				//memDC.SetPixelV(x, IMAGE_Y-1-y, RGB(r, g, b));
				image[i][j].r = rr;
				image[i][j].g = gg;
				image[i][j].b = bb;
				//
				///////////////////////////////////////////////////////////////
				// Copy randomly created strokes in memDC2
				//memDC2.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
			}
		}
		//EndTimer("_Result_FBL.txt");
		tmp.copy(image); // update tmp
	}

	dc.SelectObject(pOldPen);
	pen.DeleteObject(); 

}

void GetAnisotropicBL(cimatrix& image, ETF& e, double sigma, double sigma2, double sigma3, double sigma4, int max_itr)
// Anisotropic Bilateral filter
// Uses straight line axis in both directions
{
	int m;
	GLubyte	r, g, b;
	//double tx, ty;
	int i, j, i_x, i_y;

	vector vn(2), vt(2), w(2);
	double x, y, d_x, d_y;
	int s;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double bb, tt, val, uu, vv;
	//double sum, c_val;
	//double weight, w_sum;
	int d1;

	double sum_r, sum_g, sum_b;
	//double weight_r, weight_g, weight_b;
	//double w_sum_r, w_sum_g, w_sum_b;
	double weight, w_sum;
	double c_val_r, c_val_g, c_val_b;
	double val_r, val_g, val_b;
	//double angle, angle1, angle2;
	double c_dist;


	int image_x = image.getRow();
	int image_y = image.getCol();

	int half_l, half_w;

	//half_l = length/2 + 1; // actually, it should bigger than l/1.414 / 2
	//half_l = length/2;
		
	cimatrix tmp;
	tmp.copy(image);

	//double sigma;
	//double sigma2;
	
	//sigma = 1.0;
	//sigma2 = sigma * 1.0;
	//sigma2 = sigma * 10.0;
	//sigma2 = sigma * 2.0;
	//sigma3 = sigma;

	vector GAU1, GAU2, GAU3, GAU4;
	MakeGaussVector(sigma, GAU1); // length of the kernel
	MakeGaussVector(sigma2, GAU2); // width of the kernel
	//MakeGaussVector(sigma3, GAU3); // for similarity function

	half_l = GAU1.getMax()-1;
	half_w = GAU2.getMax()-1;

	MakeGaussVectorBL(sigma3, GAU3, 300); // 300 array elements created
	MakeGaussVectorBL(sigma4, GAU4, 300); // 300 array elements created
	
	//angle_thres1 = 90;
	//angle_thres2 = 180 - angle_thres1;
	double step_size;
	double step_size1, step_size2;

	//step_size = 1.0;
	//step_size = 2.0;
	step_size = 1.0;
	step_size1 = step_size2 = 1.0;

	//half_l = half_w = 25; // for TVCG
	
	for (m = 0; m < max_itr; m++) {
		/////////////////////////////////////////
		// FBL along main axis!
		//StartTimer();
		/////////////////////////////////////////
		// FBL along minor axis!
		for (j = 0; j < IMAGE_Y; j++) {
			for (i = 0; i < IMAGE_X; i++) {
				//sum = 0.0;
				sum_r = sum_g = sum_b = 0.0;
				weight = 0.0;
				w_sum = 0.0;
				///////////////////////////////////////////
				// First add the center intensity
				//weight = 1.0;
				//c_val = (double)tmp[i][j];
				//sum += weight * c_val;
				//w_sum += weight;
				////////////////////////////////////////////////
				// One half
				d_x = (double)i; d_y = (double)j; 
				i_x = i;	i_y = j;
				//c_val = (double)tmp[i][j];
				//c_val_r = (double)tmp[i][j].r;
				//c_val_g = (double)tmp[i][j].g;
				//c_val_b = (double)tmp[i][j].b;
				vt[0] = e[i_x][i_y].tx;
				vt[1] = e[i_x][i_y].ty;
				if (vt[0] == 0.0 && vt[1] == 0.0) break;
				//vn.make_unit();
				///////////////////////////
				//d_x = c_x + vt[0] * t * step_size; // how narrow the height of rectangle kernel
				//d_y = c_y + vt[1] * t * step_size; // depends on step_size
				////////////////////////////////////////////////
				// Use the center color at each minor axis!
				c_val_r = (double)tmp[i_x][i_y].r;
				c_val_g = (double)tmp[i_x][i_y].g;
				c_val_b = (double)tmp[i_x][i_y].b;
				////////////////////////////////////////////////
				for (s = -half_l; s <= half_l; s++) { // width of Gaussian kernel
					////////////////////////
					x = d_x + vt[0] * s;
					y = d_y + vt[1] * s;
					/////////////////////////////////////////////////////
					if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
						continue;
					/////////////////////////////////////////////////
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					//val = (double)tmp[x1][y1];
					val_r = (double)tmp[x1][y1].r;
					val_g = (double)tmp[x1][y1].g;
					val_b = (double)tmp[x1][y1].b;
					//val = (val_r + val_g + val_b) / 3.0;
					/////////////////////////////////////////////////////////
					//weight_r = weight_g = weight_b = GAU1[k];
					d1 = ABS(s);
					//weight = GAU2[d1];
					weight = GAU2[d1];
					//d3 = round( ABS(val - c_val) );
					c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
					//weight *= gauss1D_norm2( c_dist/1.732, 0.0, sigma4);
					//weight *= gauss1D_bl_macro( c_dist/1.732, sigma4);
					//weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma4);
					weight *= GAU4[ (int)(c_dist/1.732) ]; // like Gaussian!
					/////////////////////////////////////////////////////
					//sum1 += image[xi][yi] * weight;
					sum_r += val_r * weight;
					sum_g += val_g * weight;
					sum_b += val_b * weight;
					//sum1 += image[xi][yi] * weight;
					//sum += val * weight;
					w_sum += weight;
				}
				/////////////////////////////////////////
				//TRACE("sum_total = %d\n", sum);
				//sum = (int) ( sum / (double)count );
				//sum /= (double)count;
				sum_r /= (double)w_sum; // normalize!
				sum_g /= (double)w_sum; // normalize!
				sum_b /= (double)w_sum; // normalize!
				//TRACE("sum = %d, count = %d\n", sum, count);
				/////////////////////////////////////
				// Didn't add a single color!
				//if (sum_r == 0.0 && sum_g == 0.0 && sum_b == 0.0) {
				if (e[i][j].tx == 0.0 && e[i][j].ty == 0.0) {
					sum_r = tmp[i][j].r;
					sum_g = tmp[i][j].g;
					sum_b = tmp[i][j].b;
				}
				/////////////////////////////////////
				//if (sum_r > 255.0) sum_r = 255.0;
				//if (sum_g > 255.0) sum_g = 255.0;
				//if (sum_b > 255.0) sum_b = 255.0;
				//if (sum_r < 0.0) sum_r = 0.0;
				//if (sum_g < 0.0) sum_g = 0.0;
				//if (sum_b < 0.0) sum_b = 0.0;
				//sum_i = (GLubyte)sum;
				//r = g = b = (GLubyte)sum;
				r = (GLubyte)sum_r;
				g = (GLubyte)sum_g;
				b = (GLubyte)sum_b;
				//TRACE("r = %d, g = %d, b = %d\n", r, g, b);
				//memDC.SetPixelV(x, IMAGE_Y-1-y, RGB(r, g, b));
				image[i][j].r = r;
				image[i][j].g = g;
				image[i][j].b = b;
				//
				///////////////////////////////////////////////////////////////
				// Copy randomly created strokes in memDC2
				//memDC2.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
			}
		}
		//EndTimer("_Result_FBL.txt");
		tmp.copy(image); // update tmp

		/////////////////////////////////////////
		// FBL along minor axis!
		for (j = 0; j < IMAGE_Y; j++) {
			for (i = 0; i < IMAGE_X; i++) {
				//sum = 0.0;
				sum_r = sum_g = sum_b = 0.0;
				weight = 0.0;
				w_sum = 0.0;
				///////////////////////////////////////////
				// First add the center intensity
				//weight = 1.0;
				//c_val = (double)tmp[i][j];
				//sum += weight * c_val;
				//w_sum += weight;
				////////////////////////////////////////////////
				// One half
				d_x = (double)i; d_y = (double)j; 
				i_x = i;	i_y = j;
				//c_val = (double)tmp[i][j];
				//c_val_r = (double)tmp[i][j].r;
				//c_val_g = (double)tmp[i][j].g;
				//c_val_b = (double)tmp[i][j].b;
				vn[0] = -e[i_x][i_y].ty;
				vn[1] = e[i_x][i_y].tx;
				if (vn[0] == 0.0 && vn[1] == 0.0) break;
				//vn.make_unit();
				///////////////////////////
				//d_x = c_x + vt[0] * t * step_size; // how narrow the height of rectangle kernel
				//d_y = c_y + vt[1] * t * step_size; // depends on step_size
				////////////////////////////////////////////////
				// Use the center color at each minor axis!
				c_val_r = (double)tmp[i_x][i_y].r;
				c_val_g = (double)tmp[i_x][i_y].g;
				c_val_b = (double)tmp[i_x][i_y].b;
				////////////////////////////////////////////////
				for (s = -half_w; s <= half_w; s++) { // width of Gaussian kernel
					////////////////////////
					x = d_x + vn[0] * s;
					y = d_y + vn[1] * s;
					/////////////////////////////////////////////////////
					if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
						continue;
					/////////////////////////////////////////////////
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					//val = (double)tmp[x1][y1];
					val_r = (double)tmp[x1][y1].r;
					val_g = (double)tmp[x1][y1].g;
					val_b = (double)tmp[x1][y1].b;
					//val = (val_r + val_g + val_b) / 3.0;
					/////////////////////////////////////////////////////////
					//weight_r = weight_g = weight_b = GAU1[k];
					d1 = ABS(s);
					//weight = GAU2[d1];
					weight = GAU2[d1];
					//d3 = round( ABS(val - c_val) );
					c_dist = dist3(c_val_r-val_r, c_val_g-val_g, c_val_b-val_b);
					//weight *= gauss1D_norm2( c_dist/1.732, 0.0, sigma4);
					//weight *= gauss1D_bl_macro( c_dist/1.732, sigma4);
					//weight *= gauss1D_bl( c_dist/1.732, 0.0, sigma4);
					weight *= GAU4[ (int)(c_dist/1.732) ]; // like Gaussian!
					/////////////////////////////////////////////////////
					//sum1 += image[xi][yi] * weight;
					sum_r += val_r * weight;
					sum_g += val_g * weight;
					sum_b += val_b * weight;
					//sum1 += image[xi][yi] * weight;
					//sum += val * weight;
					w_sum += weight;
				}
				/////////////////////////////////////////
				//TRACE("sum_total = %d\n", sum);
				//sum = (int) ( sum / (double)count );
				//sum /= (double)count;
				sum_r /= (double)w_sum; // normalize!
				sum_g /= (double)w_sum; // normalize!
				sum_b /= (double)w_sum; // normalize!
				//TRACE("sum = %d, count = %d\n", sum, count);
				/////////////////////////////////////
				// Didn't add a single color!
				//if (sum_r == 0.0 && sum_g == 0.0 && sum_b == 0.0) {
				if (e[i][j].tx == 0.0 && e[i][j].ty == 0.0) {
					sum_r = tmp[i][j].r;
					sum_g = tmp[i][j].g;
					sum_b = tmp[i][j].b;
				}
				/////////////////////////////////////
				//if (sum_r > 255.0) sum_r = 255.0;
				//if (sum_g > 255.0) sum_g = 255.0;
				//if (sum_b > 255.0) sum_b = 255.0;
				//if (sum_r < 0.0) sum_r = 0.0;
				//if (sum_g < 0.0) sum_g = 0.0;
				//if (sum_b < 0.0) sum_b = 0.0;
				//sum_i = (GLubyte)sum;
				//r = g = b = (GLubyte)sum;
				r = (GLubyte)sum_r;
				g = (GLubyte)sum_g;
				b = (GLubyte)sum_b;
				//TRACE("r = %d, g = %d, b = %d\n", r, g, b);
				//memDC.SetPixelV(x, IMAGE_Y-1-y, RGB(r, g, b));
				image[i][j].r = r;
				image[i][j].g = g;
				image[i][j].b = b;
				//
				///////////////////////////////////////////////////////////////
				// Copy randomly created strokes in memDC2
				//memDC2.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
			}
		}
		//EndTimer("_Result_FBL.txt");
		tmp.copy(image); // update tmp
	}
}










void GetMedian(imatrix& image, double sigma, int max_itr)
// median filter
{
	int m;
	//double tx, ty;
	int i, j;

	vector vn(2), vt(2), vt_old(2);
	int x, y;
	int s, t;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double bb, tt, uu, vv;
	int val;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half_w;
		
	imatrix tmp;
	tmp.copy(image);

	vector GAU1;
	MakeGaussVector(sigma, GAU1); // entire inner circle
	//MakeGaussVector(sigma2, GAU2); // entire outer circle
	//MakeGaussVector(sigma3, GAU3); // for similarity function

	//////////////////////////////
	//half_w = GAU1.getMax()-1;
	//if (half_w == 1) half_w = 0;
	/////////////////////////////////
	
	half_w = GAU1.getMax()-1;
	//half_l = 5;

	//TRACE("half_w = %d\n", half_w);
	TRACE("half_w = %d\n", half_w);

	deque<int> vec;
	int flow_median;
		
	double step_size1 = 1.0;
	
	for (m = 0; m < max_itr; m++) {
		for (j = 0; j < IMAGE_Y; j++) {
			for (i = 0; i < IMAGE_X; i++) {
				vec.clear();
				///////////////////////////////
				for (s = -half_w; s <= half_w; s++) {
					for (t = -half_w; t <= half_w; t++) {
						x = i+s;
						y = j+t;
						/////////////////////////////////////////////////////
						if (x > IMAGE_X-1 || x < 0 || y > IMAGE_Y-1 || y < 0) 
							continue;
						val = tmp[x][y];
						vec.push_back(val);
						/////////////////////////////////////////////////////////
					}
				}
				sort(vec.begin(), vec.end()); 
				flow_median = vec[vec.size()/2]; 
				//TRACE("r = %d, g = %d, b = %d\n", r, g, b);
				//memDC.SetPixelV(x, IMAGE_Y-1-y, RGB(r, g, b));
				image[i][j] = flow_median;
				//
				///////////////////////////////////////////////////////////////
				// Copy randomly created strokes in memDC2
				//memDC2.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
			}
		}
		tmp.copy(image); // update tmp
	} // end for m
}



void CopyCol2GrayImage(int image_x, int image_y, cimatrix& cmap, imatrix& image) 
{
	int x, y;
	GLubyte r;

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r = (GLubyte) ( ((float)cmap[x][y].r + (float)cmap[x][y].g + (float)cmap[x][y].b) / 3 );
			image[x][y] = r;
		}
	}
}

void CopyCol2GrayImage2(cimatrix& cmap, imatrix& image) 
{
	int x, y;
	GLubyte r;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			r = (GLubyte) ( ((float)cmap[x][y].r + (float)cmap[x][y].g + (float)cmap[x][y].b) / 3 );
			image[x][y] = r;
		}
	}
}

void DownsampleGrayImage(imatrix& gray, imatrix& gray_s) 
{
	int x, y;
	//GLubyte r;

	int image_x = gray.getRow();
	int image_y = gray.getCol();

	double lt, rt, lb, rb, avg;

	image_x = image_x / 2;
	image_y = image_y / 2;

	gray_s.init(image_x, image_y);

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//r = (GLubyte) ( ((float)cmap[x][y].r + (float)cmap[x][y].g + (float)cmap[x][y].b) / 3 );
			lt = (double)gray[x*2][y*2];
			rt = (double)gray[x*2+1][y*2];
			lb = (double)gray[x*2][y*2+1];
			rb = (double)gray[x*2+1][y*2+1];
			avg = (lt + rt + lb + rb) / 4.0;
			gray_s[x][y] = round(avg);
		}
	}
}

void UpsampleGrayImage(imatrix& gray_s, imatrix& gray) 
{
	int x, y;
	//GLubyte r;

	int image_x = gray_s.getRow();
	int image_y = gray_s.getCol();

	double l, r, b, t;

	image_x = image_x * 2;
	image_y = image_y * 2;

    gray.init(image_x, image_y);

	matrix tmp(image_x, image_y);

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			tmp[x][y] = (double)gray_s[x/2][y/2];
		}
	}

	for (y = 0; y < image_y-1; y++) {
		for (x = 0; x < image_x-1; x++) {
			if (x % 2 == 1 && y % 2 == 0) { // odd x, even y
				tmp[x][y] = (tmp[x-1][y] + tmp[x+1][y]) / 2.0;
			}
			else if (x % 2 == 0 && y % 2 == 1) { // even x, odd y
				tmp[x][y] = (tmp[x][y-1] + tmp[x][y+1]) / 2.0;
			}
		}
	}
	for (y = 1; y < image_y-1; y+=2) {
		for (x = 1; x < image_x-1; x+=2) {
			l = tmp[x-1][y];
			r = tmp[x+1][y];
			t = tmp[x][y+1];
			b = tmp[x][y-1];
			tmp[x][y] = (l + r + t + b) / 4.0;
		}
	}
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			gray[x][y] = round(tmp[x][y]);
		}
	}
}

void UpsampleColorImage(cimatrix& cmap_s, cimatrix& cmap) 
{
	int i, j;
	//GLubyte r;

	int image_x_s = cmap_s.getRow();
	int image_y_s = cmap_s.getCol();

	int image_x = image_x_s * 2;
	int image_y = image_y_s * 2;

    cmap.init(image_x, image_y);

	cmatrix tmp(image_x, image_y);
	
	double x, y, lb, rb, lt, rt, uu, vv, bb, tt, val;
	int x1, x2, y1, y2;

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
				x = 1.0/(double)(image_x-1) * (double)i;	
				y = 1.0/(double)(image_y-1) * (double)j;	
				x1 = (int)(x * (double)(image_x_s-1)); 
				x2 = x1 + 1; if (x2 > image_x_s-1) x2 = image_x_s-1;
				y1 = (int)(y * (double)(image_y_s-1)); 
				y2 = y1 + 1; if (y2 > image_y_s-1) y2 = image_y_s-1;
				uu = x * (double)(image_x_s-1) - (double)x1;
				vv = y * (double)(image_y_s-1) - (double)y1;
				////////////////////////////////////////
				lb = (double)cmap_s[x1][y1].r;
				rb = (double)cmap_s[x2][y1].r;
				lt = (double)cmap_s[x1][y2].r;
				rt = (double)cmap_s[x2][y2].r;
				bb = (1-uu) * lb + uu * rb;
				tt = (1-uu) * lt + uu * rt;
				val = (1-vv) * bb + vv * tt;
				tmp[i][j].r = val;
				///////////////////////////////////
				lb = (double)cmap_s[x1][y1].g;
				rb = (double)cmap_s[x2][y1].g;
				lt = (double)cmap_s[x1][y2].g;
				rt = (double)cmap_s[x2][y2].g;
				bb = (1-uu) * lb + uu * rb;
				tt = (1-uu) * lt + uu * rt;
				val = (1-vv) * bb + vv * tt;
				tmp[i][j].g = val;
				///////////////////////////////////
				lb = (double)cmap_s[x1][y1].b;
				rb = (double)cmap_s[x2][y1].b;
				lt = (double)cmap_s[x1][y2].b;
				rt = (double)cmap_s[x2][y2].b;
				bb = (1-uu) * lb + uu * rb;
				tt = (1-uu) * lt + uu * rt;
				val = (1-vv) * bb + vv * tt;
				tmp[i][j].b = val;
				///////////////////////////////////
		}
	}
	
	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			cmap[i][j].r = round(tmp[i][j].r);
			cmap[i][j].g = round(tmp[i][j].g);
			cmap[i][j].b = round(tmp[i][j].b);
		}
	}
}

void UpsampleColorImageNearest(cimatrix& cmap_s, cimatrix& cmap) 
// using nearest neighbor
{
	int x, y;

	int image_x_s = cmap_s.getRow();
	int image_y_s = cmap_s.getCol();

	int image_x = image_x_s * 2;
	int image_y = image_y_s * 2;

    cmap.init(image_x, image_y);

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			cmap[x][y].r = cmap_s[x/2][y/2].r;
			cmap[x][y].g = cmap_s[x/2][y/2].g;
			cmap[x][y].b = cmap_s[x/2][y/2].b;
		}
	}
}

void UpsampleGrayImageNearest(imatrix& gray_s, imatrix& gray) 
// using nearest neighbor
{
	int x, y;

	int image_x_s = gray_s.getRow();
	int image_y_s = gray_s.getCol();

	int image_x = image_x_s * 2;
	int image_y = image_y_s * 2;

    gray.init(image_x, image_y);

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			gray[x][y] = gray_s[x/2][y/2];
		}
	}
}


void DownsampleColorImage(cimatrix& cmap, cimatrix& cmap_s) 
{
	int x, y;
	//GLubyte r, g, b;

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();

	double lt, rt, lb, rb, avg;

	image_x = image_x / 2;
	image_y = image_y / 2;

	cmap_s.init(image_x, image_y);

	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//r = (GLubyte) ( ((float)cmap[x][y].r + (float)cmap[x][y].g + (float)cmap[x][y].b) / 3 );
			lt = (double)cmap[x*2][y*2].r;
			rt = (double)cmap[x*2+1][y*2].r;
			lb = (double)cmap[x*2][y*2+1].r;
			rb = (double)cmap[x*2+1][y*2+1].r;
			avg = (lt + rt + lb + rb) / 4.0;
			cmap_s[x][y].r = (GLubyte)round(avg);
			///////////////////////////////
			lt = (double)cmap[x*2][y*2].g;
			rt = (double)cmap[x*2+1][y*2].g;
			lb = (double)cmap[x*2][y*2+1].g;
			rb = (double)cmap[x*2+1][y*2+1].g;
			avg = (lt + rt + lb + rb) / 4.0;
			cmap_s[x][y].g = (GLubyte)round(avg);
			////////////////////////////////
			lt = (double)cmap[x*2][y*2].b;
			rt = (double)cmap[x*2+1][y*2].b;
			lb = (double)cmap[x*2][y*2+1].b;
			rb = (double)cmap[x*2+1][y*2+1].b;
			avg = (lt + rt + lb + rb) / 4.0;
			cmap_s[x][y].b = (GLubyte)round(avg);
		}
	}
}


void GetFlowLineContrast(imatrix& image, Field& gfield, double sigma, int max_itr)
// flow-based line contrast enhancing filter
{
	int k, m;
	//double tx, ty;
	int i, j, i_x, i_y;

	vector vn(2), vn_old(2);
	double x, y, d_x, d_y;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double bb, tt, uu, vv;
	double sum, c_val;
	double weight, w_sum;
	double val;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half_l, half_w;
		
	imatrix tmp;
	tmp.copy(image);

	vector GAU1;
	MakeGaussVector(sigma, GAU1); 

	half_l = GAU1.getMax()-1;
	//half_l = 5;

	//TRACE("half_w = %d\n", half_w);
	TRACE("half_l = %d\n", half_l);
	half_w = 0;

	double angle, angle_thres1, angle_thres2;
	angle_thres1 = 60;
	angle_thres2 = 180 - angle_thres1;
		
	double step_size, step_size1;
	
	step_size = 1.0;

	for (m = 0; m < max_itr; m++) {
		for (j = 0; j < IMAGE_Y; j++) {
			for (i = 0; i < IMAGE_X; i++) {
				sum = 0.0;
				weight = 0.0;
				w_sum = 0.0;
				////////////////////////////////////////////////
				// One half
				d_x = (double)i; d_y = (double)j; 
				i_x = i;	i_y = j;
				c_val = (double)tmp[i][j] / 255.0;
				vn_old[0] = gfield[i_x][i_y].gx;
				vn_old[1] = gfield[i_x][i_y].gy;
				vn_old.make_unit();
				///////////////////////////////
				for (k = 0; k <= half_l; k++) {
					vn[0] = gfield[i_x][i_y].gx;
					vn[1] = gfield[i_x][i_y].gy;
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
					else break;
					//////////////////////////////////////////////////////
					x = d_x;
					y = d_y;
					/////////////////////////////////////////////////////
					if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
						continue;
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					val = (double)tmp[x1][y1] / 255.0;
					/////////////////////////////////////////////////////////
					weight = GAU1[k];
					//d3 = round( ABS(val - c_val) );
					//weight *= GAU3[d3];
					//weight *= gauss2( val - c_val, 0.0, sigma3);
					//weight *= gauss1D_norm( (val - c_val)/255.0, 0.0, sigma3);
					/////////////////////////////////////////////////////
					//sum1 += image[xi][yi] * weight;
					sum += val * weight;
					w_sum += weight;
					/////////////////////////////////////////
					d_x += vn[0] * step_size1; // accurate new location x
					d_y += vn[1] * step_size1; // accurate new location y
					/////////////////////////////////////////
					i_x = round(d_x); // integer version of new location x
					i_y = round(d_y); // integer version of new location y
					//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
					if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
					//////////////////////////////////
					vn_old.copy(vn);
				}
				////////////////////////////////////////////////
				// One half
				d_x = (double)i; d_y = (double)j; 
				i_x = i;	i_y = j;
				vn_old[0] = -gfield[i_x][i_y].gx;
				vn_old[1] = -gfield[i_x][i_y].gy;
				vn_old.make_unit();
				///////////////////////////////
				for (k = 0; k <= half_l; k++) {
					vn[0] = -gfield[i_x][i_y].gx;
					vn[1] = -gfield[i_x][i_y].gy;
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
					else break;
					//////////////////////////////////////////////////////
					x = d_x;
					y = d_y;
					/////////////////////////////////////////////////////
					if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
						continue;
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					val = (double)tmp[x1][y1] / 255.0;
					/////////////////////////////////////////////////////////
					weight = GAU1[k];
					//d3 = round( ABS(val - c_val) );
					//weight *= GAU3[d3];
					//weight *= gauss2( val - c_val, 0.0, sigma3);
					//weight *= gauss1D_norm( (val - c_val)/255.0, 0.0, sigma3);
					/////////////////////////////////////////////////////
					//sum1 += image[xi][yi] * weight;
					sum += val * weight;
					w_sum += weight;
					/////////////////////////////////////////
					d_x += vn[0] * step_size1; // accurate new location x
					d_y += vn[1] * step_size1; // accurate new location y
					/////////////////////////////////////////
					i_x = round(d_x); // integer version of new location x
					i_y = round(d_y); // integer version of new location y
					//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
					if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
					//////////////////////////////////
					vn_old.copy(vn);
				}
				//TRACE("sum_total = %d\n", sum);
				//sum = (int) ( sum / (double)count );
				//sum /= (double)count;
				sum /= (double)w_sum; // normalize!
				//TRACE("sum = %d, count = %d\n", sum, count);
				//c_val += 5.0 * (c_val - sum); 
				c_val += tanh( 5.0 * (c_val - sum) );
				//TRACE("c_val = %.6f\n", c_val);
				if (c_val > 1.0) c_val = 1.0;
				else if (c_val < 0.0) c_val = 0.0;
				image[i][j] = round(c_val * 255.);
				//image[i][j] = r;
			}
		}
		tmp.copy(image); // update tmp
	} // end for m
}

void GetFlowLineContrast2(imatrix& image, Field& gfield, imatrix& dog, double sigma, int max_itr)
// flow-based line contrast enhancing filter
// use DOG black area
{
	int k, m;
	//double tx, ty;
	int i, j, i_x, i_y;

	vector vn(2), vn_old(2);
	double x, y, d_x, d_y;
	int x1, y1;
	//int x2, y2;
	//double lb, lt, rb, rt;
	//double bb, tt, uu, vv;
	double sum, c_val;
	double weight, w_sum;
	double val;

	int image_x = image.getRow();
	int image_y = image.getCol();

	int half_l, half_w;
		
	imatrix tmp;
	tmp.copy(image);

	vector GAU1;
	MakeGaussVector(sigma, GAU1); 

	half_l = GAU1.getMax()-1;
	//half_l = 5;

	//TRACE("half_w = %d\n", half_w);
	TRACE("half_l = %d\n", half_l);
	half_w = 0;

	double angle, angle_thres1, angle_thres2;
	angle_thres1 = 60;
	angle_thres2 = 180 - angle_thres1;
		
	double step_size, step_size1;
	
	step_size = 1.0;
		
	bool dog_flag = false;

	for (m = 0; m < max_itr; m++) {
		for (j = 0; j < IMAGE_Y; j++) {
			for (i = 0; i < IMAGE_X; i++) {
				sum = 0.0;
				weight = 0.0;
				w_sum = 0.0;
				dog_flag = false;
				////////////////////////////////////////////////
				// One half
				d_x = (double)i; d_y = (double)j; 
				i_x = i;	i_y = j;
				c_val = (double)tmp[i][j] / 255.0;
				vn_old[0] = gfield[i_x][i_y].gx;
				vn_old[1] = gfield[i_x][i_y].gy;
				vn_old.make_unit();
				///////////////////////////////
				for (k = 0; k <= half_l; k++) {
					vn[0] = gfield[i_x][i_y].gx;
					vn[1] = gfield[i_x][i_y].gy;
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
					else break;
					//////////////////////////////////////////////////////
					x = d_x;
					y = d_y;
					/////////////////////////////////////////////////////
					if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
						continue;
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					val = (double)tmp[x1][y1] / 255.0;
					//////////////////////////////////////
					if (dog[x1][y1] == 0) dog_flag = true;
					/////////////////////////////////////////////////////////
					weight = GAU1[k];
					//d3 = round( ABS(val - c_val) );
					//weight *= GAU3[d3];
					//weight *= gauss2( val - c_val, 0.0, sigma3);
					//weight *= gauss1D_norm( (val - c_val)/255.0, 0.0, sigma3);
					/////////////////////////////////////////////////////
					//sum1 += image[xi][yi] * weight;
					sum += val * weight;
					w_sum += weight;
					/////////////////////////////////////////
					d_x += vn[0] * step_size1; // accurate new location x
					d_y += vn[1] * step_size1; // accurate new location y
					/////////////////////////////////////////
					i_x = round(d_x); // integer version of new location x
					i_y = round(d_y); // integer version of new location y
					//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
					if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
					//////////////////////////////////
					vn_old.copy(vn);
				}
				////////////////////////////////////////////////
				// Other half
				d_x = (double)i; d_y = (double)j; 
				i_x = i;	i_y = j;
				vn_old[0] = -gfield[i_x][i_y].gx;
				vn_old[1] = -gfield[i_x][i_y].gy;
				vn_old.make_unit();
				///////////////////////////////
				for (k = 0; k <= half_l; k++) {
					vn[0] = -gfield[i_x][i_y].gx;
					vn[1] = -gfield[i_x][i_y].gy;
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
					else break;
					//////////////////////////////////////////////////////
					x = d_x;
					y = d_y;
					/////////////////////////////////////////////////////
					if (x > (double)IMAGE_X-1 || x < 0.0 || y > (double)IMAGE_Y-1 || y < 0.0) 
						continue;
					x1 = round(x);	if (x1 < 0) x1 = 0; if (x1 > IMAGE_X-1) x1 = IMAGE_X-1;
					y1 = round(y);	if (y1 < 0) y1 = 0; if (y1 > IMAGE_Y-1) y1 = IMAGE_Y-1;
					val = (double)tmp[x1][y1] / 255.0;
					//////////////////////////////////////
					if (dog[x1][y1] == 0) dog_flag = true;
					/////////////////////////////////////////////////////////
					weight = GAU1[k];
					/////////////////////////////////////////////////////
					//sum1 += image[xi][yi] * weight;
					sum += val * weight;
					w_sum += weight;
					/////////////////////////////////////////
					d_x += vn[0] * step_size1; // accurate new location x
					d_y += vn[1] * step_size1; // accurate new location y
					/////////////////////////////////////////
					i_x = round(d_x); // integer version of new location x
					i_y = round(d_y); // integer version of new location y
					//if (int_i < 0 || int_i > IMAGE_X-1 || int_j < 0 || int_j > IMAGE_Y-1) break;
					if (d_x < 0 || d_x > IMAGE_X-1 || d_y < 0 || d_y > IMAGE_Y-1) break;
					//////////////////////////////////
					vn_old.copy(vn);
				}
				sum /= (double)w_sum; // normalize!
				//TRACE("sum = %d, count = %d\n", sum, count);
				//c_val += 5.0 * (c_val - sum); 
				c_val += tanh( 2.0 * (c_val - sum) );
				//TRACE("c_val = %.6f\n", c_val);
				if (c_val > 1.0) c_val = 1.0;
				else if (c_val < 0.0) c_val = 0.0;
				////////////////////////////////////
				if (dog_flag) image[i][j] = round(c_val * 255.);
				//image[i][j] = r;
			}
		}
		tmp.copy(image); // update tmp
	} // end for m
}


void GetAnisotropic(imatrix& image, int max_itr)
// Perona-Malik anisotropic diffusion
// Use the gradient magnitude
{
	int i, j, k, m;

	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	matrix tmp, tmp2;
	tmp.init(image_x, image_y);
	tmp2.init(image_x, image_y);

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			tmp[i][j] = (double)image[i][j];
		}
	}

	double c_val, val;
	double diff, weight, lambda;
	double n_val;

	lambda = 0.25;

	for (m = 0; m < max_itr; m++) {
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				////////////////////////////////////////////////
				c_val = tmp[i][j];
				n_val = c_val;
				///////////////////////////////
				for (k = 0; k < 4; k++) {
					switch (k) {
						case 0: x = i+1; y = j; break;
						case 1: x = i-1; y = j; break;
						case 2: x = i; y = j+1; break;
						case 3: x = i; y = j-1; break;
					}
					if (x > IMAGE_X-1)	x = IMAGE_X-1;
					if (x < 0)			x = 0;
					if (y > IMAGE_Y-1)	y = IMAGE_Y-1;
					if (y < 0)			y = 0;
					val = tmp[x][y];
					diff = (val - c_val);
					weight = 1 / (1 + diff*diff);
					n_val += lambda * weight * diff;
					//TRACE("n_val = %f\n", n_val);
				}
				/////////////////////////////////////////////////////
				tmp2[i][j] = n_val;
				//TRACE("n_val = %f\n", n_val);
				//TRACE("tmp2[%d][%d] = %f\n", i, j, tmp2[i][j]);
			}
		}
		tmp.copy(tmp2); // update tmp
	} // end for m

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			//TRACE("tmp[%d][%d] = %f\n", i, j, tmp[i][j]);
			if (tmp[i][j] > 255.0) {
				TRACE("tmp[%d][%d] = %f\n", i, j, tmp[i][j]);
				tmp[i][j] = 255.0;
			}
			else if (tmp[i][j] < 0.0) {
				TRACE("tmp[%d][%d] = %f\n", i, j, tmp[i][j]);
				tmp[i][j] = 0.0;
			}
			image[i][j] = round(tmp[i][j]);
		}
	}
}

void GetAnisotropic2(imatrix& image, imatrix& line, int max_itr)
// Perona-Malik anisotropic diffusion
// Use constraint image!
{
	int i, j, k, m;

	int x, y;

	int image_x = image.getRow();
	int image_y = image.getCol();

	matrix tmp, tmp2;
	tmp.init(image_x, image_y);
	tmp2.init(image_x, image_y);

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			tmp[i][j] = (double)image[i][j];
		}
	}

	double c_val, val;
	double diff, weight, lambda;
	double n_val;

	//////////////////////////
	//lambda = 0.25; // suggested value by paper
	lambda = 0.25; // suggested value by paper
	//lambda = 1.0;
	//////////////////////////

	for (m = 0; m < max_itr; m++) {
		for (j = 0; j < image_y; j++) {
			for (i = 0; i < image_x; i++) {
				////////////////////////////////////////////////
				c_val = tmp[i][j];
				n_val = c_val;
				//weight = line[i][j] / 255.0; // [0, 1]
				//TRACE("weight = %f\n", weight);
				///////////////////////////////
				if (line[i][j] == 255) {
					for (k = 0; k < 4; k++) {
						switch (k) {
							case 0: x = i+1; y = j; break;
							case 1: x = i-1; y = j; break;
							case 2: x = i; y = j+1; break;
							case 3: x = i; y = j-1; break;
						}
						if (x > IMAGE_X-1)	x = IMAGE_X-1;
						if (x < 0)			x = 0;
						if (y > IMAGE_Y-1)	y = IMAGE_Y-1;
						if (y < 0)			y = 0;
						val = tmp[x][y];
						diff = (val - c_val);
						//////////////////////////////////
						if (line[x][y] == 0) weight = 0.0;
						else weight = 1.0;
						//weight = 1 / (1 + diff*diff);
						///////////////////////////////////
						n_val += lambda * weight * diff;
						//TRACE("n_val = %f\n", n_val);
					}
				}
				/////////////////////////////////////////////////////
				tmp2[i][j] = n_val;
				//TRACE("n_val = %f\n", n_val);
				//TRACE("tmp2[%d][%d] = %f\n", i, j, tmp2[i][j]);
			}
		}
		tmp.copy(tmp2); // update tmp
	} // end for m

	for (j = 0; j < image_y; j++) {
		for (i = 0; i < image_x; i++) {
			//TRACE("tmp[%d][%d] = %f\n", i, j, tmp[i][j]);
			if (tmp[i][j] > 255.0) {
				TRACE("tmp[%d][%d] = %f\n", i, j, tmp[i][j]);
				tmp[i][j] = 255.0;
			}
			else if (tmp[i][j] < 0.0) {
				TRACE("tmp[%d][%d] = %f\n", i, j, tmp[i][j]);
				tmp[i][j] = 0.0;
			}
			image[i][j] = round(tmp[i][j]);
		}
	}
}



void DrawGradientFieldArrow(CDC& dc, Field& gfield, GLubyte r, GLubyte g, GLubyte b)
// Show arrows to clear show the vector direction!
{
	int x, y;
	//GLubyte	r2, g2, b2;
	//int	grad;
	CPen pen, *pOldPen;
	int i, j;
	double gx, gy;

	int image_x = gfield.getRow();
	int image_y = gfield.getCol();
	
	///////////////////////////////////////////
	//ClearMemDC(&memDC); // clear the canvas white
	//////////////////////////////////////////////////////
	//ClearMemDC(&double_buffer); // clear the canvas white
	//////////////////////////////////////////////

	int box_size = 9; // 7x7 size
	int half = box_size / 2;

	int max_x = image_x / box_size; 
	TRACE("max_x = %d\n", max_x);
	int max_y = image_y / box_size;
	TRACE("max_y = %d\n", max_y);

	int st_x, st_y, d_x, d_y;
	double t;

	for (j = 0; j < max_y; j++) {
		for (i = 0; i < max_x; i++) {
			///////////////////////
			// draw the center point
			st_x = i * box_size; 
			st_y = j * box_size; 
			x = st_x + half;
			y = st_y + half;
			pen.DeleteObject();
            pen.CreatePen(PS_SOLID, 3, RGB(r, g, b));
			pOldPen = (CPen *)dc.SelectObject(&pen);
			dc.MoveTo(x, image_y-1-y);
			dc.LineTo(x, image_y-1-y);
			////////////////////////////////////
			// draw the arrow line
			gx = gfield[x][y].gx;
			gy = gfield[x][y].gy;
			if (gx == 0.0 && gy == 0.0) continue; // draw no arrow
			if (ABS(gx) >= ABS(gy)) { // meets right wall or left wall
				if (gx > 0.0) { // left wall (we draw arrow tail, not head, so reversed!)
					//d_x = x + 3;
					d_x = x - half; 
					t = (d_x - x) / gx; 
					d_y = round(y + gy * t);
				}
				else { // gx < 0.0, right wall
					//d_x = x - 3;
					d_x = x + half;
					t = (d_x - x) / gx; 
					d_y = round(y + gy * t);
				}
			}
			else if (ABS(gx) < ABS(gy)) { // meets upper wall or bottom wall
				if (gy > 0.0) { // bottom wall
					//d_y = y + 3;
					d_y = y - half;
					t = (d_y - y) / gy; 
					d_x = round(x + gx * t);
				}
				else { // gy < 0.0, up wall
					//d_y = y - 3;
					d_y = y + half;
					t = (d_y - y) / gy; 
					d_x = round(x + gx * t);
				}
			}
			pen.DeleteObject();
			pen.CreatePen(PS_SOLID, 1, RGB(r, g, b));
			pOldPen = (CPen *)dc.SelectObject(&pen);
			dc.MoveTo(x, image_y-1-y);
			dc.LineTo(d_x, image_y-1-d_y);
		}
	}
	dc.SelectObject(pOldPen);
}

void DrawTangentFieldArrow(CDC& dc, Field& gfield, GLubyte r, GLubyte g, GLubyte b)
// Show arrows to clear show the vector direction!
{
	int x, y;
	//GLubyte	r2, g2, b2;
	//int	grad;
	CPen pen, *pOldPen;
	int i, j;
	double gx, gy;

	int image_x = gfield.getRow();
	int image_y = gfield.getCol();
	
	///////////////////////////////////////////
	//ClearMemDC(&memDC); // clear the canvas white
	//////////////////////////////////////////////////////
	//ClearMemDC(&double_buffer); // clear the canvas white
	//////////////////////////////////////////////

	int box_size = 9; // 7x7 size
	int half = box_size / 2;

	int max_x = image_x / box_size; 
	TRACE("max_x = %d\n", max_x);
	int max_y = image_y / box_size;
	TRACE("max_y = %d\n", max_y);

	int st_x, st_y, d_x, d_y;
	double t;

	for (j = 0; j < max_y; j++) {
		for (i = 0; i < max_x; i++) {
			///////////////////////
			// draw the center point
			st_x = i * box_size; 
			st_y = j * box_size; 
			x = st_x + half;
			y = st_y + half;
			pen.DeleteObject();
            pen.CreatePen(PS_SOLID, 3, RGB(r, g, b));
			pOldPen = (CPen *)dc.SelectObject(&pen);
			dc.MoveTo(x, image_y-1-y);
			dc.LineTo(x, image_y-1-y);
			////////////////////////////////////
			// draw the arrow line
			gx = -gfield[x][y].gy; // tangent direction!
			gy = gfield[x][y].gx;
			if (gx == 0.0 && gy == 0.0) continue; // draw no arrow
			if (ABS(gx) >= ABS(gy)) { // meets right wall or left wall
				if (gx > 0.0) { // left wall (we draw arrow tail, not head, so reversed!)
					//d_x = x + 3;
					d_x = x - half; 
					t = (d_x - x) / gx; 
					d_y = round(y + gy * t);
				}
				else { // gx < 0.0, right wall
					//d_x = x - 3;
					d_x = x + half;
					t = (d_x - x) / gx; 
					d_y = round(y + gy * t);
				}
			}
			else if (ABS(gx) < ABS(gy)) { // meets upper wall or bottom wall
				if (gy > 0.0) { // bottom wall
					//d_y = y + 3;
					d_y = y - half;
					t = (d_y - y) / gy; 
					d_x = round(x + gx * t);
				}
				else { // gy < 0.0, up wall
					//d_y = y - 3;
					d_y = y + half;
					t = (d_y - y) / gy; 
					d_x = round(x + gx * t);
				}
			}
			pen.DeleteObject();
			pen.CreatePen(PS_SOLID, 1, RGB(r, g, b));
			pOldPen = (CPen *)dc.SelectObject(&pen);
			dc.MoveTo(x, image_y-1-y);
			dc.LineTo(d_x, image_y-1-d_y);
		}
	}
	dc.SelectObject(pOldPen);
}

void DrawETFArrow(CDC& dc, ETF& e, GLubyte r, GLubyte g, GLubyte b)
// Show arrows to clear show the vector direction!
{
	int x, y;
	//GLubyte	r2, g2, b2;
	//int	grad;
	CPen pen, *pOldPen;
	int i, j;
	double gx, gy;

	int image_x = e.getRow();
	int image_y = e.getCol();
	
	///////////////////////////////////////////
	//ClearMemDC(&memDC); // clear the canvas white
	//////////////////////////////////////////////////////
	//ClearMemDC(&double_buffer); // clear the canvas white
	//////////////////////////////////////////////

	int box_size = 9; // 7x7 size
	int half = box_size / 2;

	int max_x = image_x / box_size; 
	TRACE("max_x = %d\n", max_x);
	int max_y = image_y / box_size;
	TRACE("max_y = %d\n", max_y);

	int st_x, st_y, d_x, d_y;
	double t;

	for (j = 0; j < max_y; j++) {
		for (i = 0; i < max_x; i++) {
			///////////////////////
			// draw the center point
			st_x = i * box_size; 
			st_y = j * box_size; 
			x = st_x + half;
			y = st_y + half;
			pen.DeleteObject();
            pen.CreatePen(PS_SOLID, 3, RGB(r, g, b));
			pOldPen = (CPen *)dc.SelectObject(&pen);
			dc.MoveTo(x, image_y-1-y);
			dc.LineTo(x, image_y-1-y);
			////////////////////////////////////
			// draw the arrow line
			gx = e[x][y].tx; // tangent direction!
			gy = e[x][y].ty;
			if (gx == 0.0 && gy == 0.0) continue; // draw no arrow
			if (ABS(gx) >= ABS(gy)) { // meets right wall or left wall
				if (gx > 0.0) { // left wall (we draw arrow tail, not head, so reversed!)
					//d_x = x + 3;
					d_x = x - half; 
					t = (d_x - x) / gx; 
					d_y = round(y + gy * t);
				}
				else { // gx < 0.0, right wall
					//d_x = x - 3;
					d_x = x + half;
					t = (d_x - x) / gx; 
					d_y = round(y + gy * t);
				}
			}
			else if (ABS(gx) < ABS(gy)) { // meets upper wall or bottom wall
				if (gy > 0.0) { // bottom wall
					//d_y = y + 3;
					d_y = y - half;
					t = (d_y - y) / gy; 
					d_x = round(x + gx * t);
				}
				else { // gy < 0.0, up wall
					//d_y = y - 3;
					d_y = y + half;
					t = (d_y - y) / gy; 
					d_x = round(x + gx * t);
				}
			}
			pen.DeleteObject();
			pen.CreatePen(PS_SOLID, 1, RGB(r, g, b));
			pOldPen = (CPen *)dc.SelectObject(&pen);
			dc.MoveTo(x, image_y-1-y);
			dc.LineTo(d_x, image_y-1-d_y);
		}
	}
	dc.SelectObject(pOldPen);
}

	



