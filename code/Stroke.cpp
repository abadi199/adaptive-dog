#include "stdafx.h"
#include "Cube.h"
#include "CubeDoc.h"
#include "CubeView.h"

#include "globals.h"
#include "Stroke.h"

#include <cmath>
#include <climits>

#define dist3(x, y, z) sqrt((x)*(x) + (y)*(y) + (z)*(z))
//#define PI	3.1415

void Stroke::Draw(CDC& dc) 
{
	CPen pen, *pOldPen;

	pen.CreatePen(PS_SOLID, (int)w, RGB(r, g, b));
	pOldPen = (CPen *)dc.SelectObject(&pen);

	dc.MoveTo(st_x, IMAGE_Y-1-st_y);
	dc.LineTo(end_x, IMAGE_Y-1-end_y);

	dc.SelectObject(pOldPen);
	pen.DeleteObject(); 
}

void Stroke::DrawTileWhite(CDC& dc) 
// equivalent to "remove" a tile
{
	CPen pen, *pOldPen;
	CBrush brush, *pOldBrush;

	////////////////////////////////////////////
	//pen.CreatePen(PS_SOLID, 1, RGB(r, g, b));
	pen.CreatePen(PS_SOLID, 1, RGB(255, 255, 255)); // white boundary
	pOldPen = (CPen *)dc.SelectObject(&pen);

	brush.CreateSolidBrush(RGB(255, 255, 255)); // white region
	pOldBrush = dc.SelectObject(&brush);
	
	/////////////////////////////////////////
	dc.Polygon(P, 4);
	///////////////////////////////////////////////////////////////////////////////////////
		
	dc.SelectObject(pOldPen);
	dc.SelectObject(pOldBrush);
	pen.DeleteObject(); 
	brush.DeleteObject(); 
}

bool Stroke::DrawTileFast(CDC& dc, int overlap_test) 
// already know P[4]
{
	CPen pen, *pOldPen;
	CBrush brush, *pOldBrush;

	//////////////////////////////////////
	if (overlap_test) {
		if (SquareOverlap(dc, P)) return true; // there is an overlap!
	}

	////////////////////////////////////////////
	//pen.CreatePen(PS_SOLID, 1, RGB(r, g, b));
	pen.CreatePen(PS_SOLID, 1, RGB(255, 255, 255)); // white boundary
	pOldPen = (CPen *)dc.SelectObject(&pen);

	brush.CreateSolidBrush(RGB(r, g, b));
	pOldBrush = dc.SelectObject(&brush);
	
	/////////////////////////////////////////
	dc.Polygon(P, 4);
	///////////////////////////////////////////////////////////////////////////////////////
		
	dc.SelectObject(pOldPen);
	dc.SelectObject(pOldBrush);
	pen.DeleteObject(); 
	brush.DeleteObject(); 

	return false;

}

bool Stroke::DrawMosaic(CDC& dc, int overlap_test) 
{
	CPen pen, *pOldPen;
	CBrush brush, *pOldBrush;

	//pen.CreatePen(PS_SOLID, (int)w, RGB(r, g, b));
	//pOldPen = (CPen *)dc.SelectObject(&pen);

	//dc.MoveTo(st_x, IMAGE_Y-1-st_y);
	//dc.LineTo(end_x, IMAGE_Y-1-end_y);

	//dc.SelectObject(pOldPen);
	//pen.DeleteObject(); 

	double tx, ty, theta;
	double2D p1, p2, p3, p4;
	double2D q1, q2, q3, q4;
	POINT P[4];
	int k;

	//////////////////////////////////////////////
	tx = -gfield[xc][yc].gy;
	ty = gfield[xc][yc].gx;

	theta = atan(ty/tx);

	p1.x = w/2; p1.y = w/2;
	p2.x = -w/2; p2.y = w/2;
	p3.x = -w/2; p3.y = -w/2;
	p4.x = w/2; p4.y = -w/2;
	
	q1.x = cos(theta) * p1.x - sin(theta) * p1.y;  
	q1.y = sin(theta) * p1.x + cos(theta) * p1.y;  
	q2.x = cos(theta) * p2.x - sin(theta) * p2.y;  
	q2.y = sin(theta) * p2.x + cos(theta) * p2.y;  
	q3.x = cos(theta) * p3.x - sin(theta) * p3.y;  
	q3.y = sin(theta) * p3.x + cos(theta) * p3.y;  
	q4.x = cos(theta) * p4.x - sin(theta) * p4.y;  
	q4.y = sin(theta) * p4.x + cos(theta) * p4.y;  

	q1.x += xc;  q1.y += yc;  
	q2.x += xc;  q2.y += yc;  
	q3.x += xc;  q3.y += yc;  
	q4.x += xc;  q4.y += yc;  
	
	P[0].x = round(q1.x); P[0].y = IMAGE_Y-1-round(q1.y); 
	P[1].x = round(q2.x); P[1].y = IMAGE_Y-1-round(q2.y); 
	P[2].x = round(q3.x); P[2].y = IMAGE_Y-1-round(q3.y); 
	P[3].x = round(q4.x); P[3].y = IMAGE_Y-1-round(q4.y); 

	for (k = 0; k < 4; k++) {
		if (P[k].x < 0) P[k].x = 0;
		if (P[k].x > IMAGE_X-1) P[k].x = IMAGE_X-1;
		if (P[k].y < 0) P[k].y = 0;
		if (P[k].y > IMAGE_Y-1) P[k].y = IMAGE_Y-1;
	}

	//////////////////////////////////////
	if (overlap_test) {
		if (SquareOverlap(dc, P)) return true; // there is an overlap!
	}

	///////////////////////////////
	//r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
	//g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
	//b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];

	//PerturbRGB(r, g, b, 70);

	////////////////////////////////////////////
	//pen.CreatePen(PS_SOLID, 1, RGB(r, g, b));
	pen.CreatePen(PS_SOLID, 1, RGB(255, 255, 255)); // white boundary
	pOldPen = (CPen *)dc.SelectObject(&pen);

	brush.CreateSolidBrush(RGB(r, g, b));
	pOldBrush = dc.SelectObject(&brush);
	
	/////////////////////////////////////////
	dc.Polygon(P, 4);
	///////////////////////////////////////////////////////////////////////////////////////
		
	dc.SelectObject(pOldPen);
	dc.SelectObject(pOldBrush);
	pen.DeleteObject(); 
	brush.DeleteObject(); 

	return false;

}

void Stroke::SetValues(int xc, int yc, int w, int st_x, int st_y, int end_x, int end_y, unsigned char r, unsigned char g, unsigned char b)
{
	this->xc = xc;
	this->yc = yc;
	//this->theta = theta;
	this->w = w;
	this->st_x = st_x;
	this->st_y = st_y;
	this->end_x = end_x;
	this->end_y = end_y;
	this->r = r;
	this->g = g;
	this->b = b;
}

void Stroke::SetValues(int xc, int yc, int w, POINT *P, unsigned char r, unsigned char g, unsigned char b)
{
	this->xc = xc;
	this->yc = yc;
	//this->theta = theta;
	this->w = w;
	//this->st_x = st_x;
	//this->st_y = st_y;
	//this->end_x = end_x;
	//this->end_y = end_y;
	this->r = r;
	this->g = g;
	this->b = b;
	this->P[0].x = P[0].x;
	this->P[0].y = P[0].y;
	this->P[1].x = P[1].x;
	this->P[1].y = P[1].y;
	this->P[2].x = P[2].x;
	this->P[2].y = P[2].y;
	this->P[3].x = P[3].x;
	this->P[3].y = P[3].y;
}


/*
void CTIPView::ImageMoment(int size, int center_x, int center_y, Stroke *stroke_ptr)
{
	int i, j;
	int	z;
	int half_x, half_y;	
	//unsigned char r, g, b;
	//unsigned char center_r, center_g, center_b;
	double r, g, b, center_r, center_g, center_b;
	double *Tbuffer;
	double	M00, M01, M10, M11, M02, M20;
	double	A, B, C;
	double	xc, yc, theta, w, l;
	double	distance;

	Tbuffer = new double [size * size]; // Color Difference Image

	half_x = half_y = (int) size / 2;

	center_r = Dbuffer[(center_y * IMAGE_X + center_x) * 3 + 0];
	center_g = Dbuffer[(center_y * IMAGE_X + center_x) * 3 + 1];
	center_b = Dbuffer[(center_y * IMAGE_X + center_x) * 3 + 2];

	M00 = M10 = M01 = M11 = M20 = M02 = 0.0;

	CClientDC dc(this);

	bitmap3.DeleteObject();
	bitmap3.CreateCompatibleBitmap(&dc, size, size);
	//memDC.CreateCompatibleDC(&dc);
	memDC2.SelectObject(&bitmap3);
	
	TRACE("half_x = %d, half_y = %d\n", half_x, half_y);
	for (j = 0; j < size; j++) {
		z = ((center_y-half_y+j) * IMAGE_X + (center_x-half_x)) * 3;
		for (i = 0; i < size; i++) {
			//r = Dbuffer[z + i * 3 + 0] / 255.0;
			//g = Dbuffer[z + i * 3 + 1] / 255.0;
			//b = Dbuffer[z + i * 3 + 2] / 255.0;
			r = Dbuffer[z + i * 3 + 0];
			g = Dbuffer[z + i * 3 + 1];
			b = Dbuffer[z + i * 3 + 2];
			distance = dist3(r - center_r, g - center_g, b - center_b);
			if (distance <= 150) // bright if similar
				Tbuffer[j * size + i] = (1 - (distance/150) * (distance/150)) * (1 - (distance/150) * (distance/150));
			else // dark if different
				Tbuffer[j * size + i] = 0.0;
			//Tbuffer[j * size + i] = dist3(r - center_r, g - center_g, b - center_b);
			//memDC2.SetPixelV(center_x - half_x + i, (IMAGE_Y-1)-(center_y - half_y + j), RGB(r, g, b));
			memDC2.SetPixelV(i, size-1-j, RGB(r, g, b));
		}
	}
	dc.BitBlt(0, 0, size, size, &memDC2, 0, 0, SRCCOPY);
	//void GetImage(int width, int height, Image image, GLubyte *Dbuffer)
	//GetImage(size, size, image, Tbuffer);
	for (j = 0; j < size; j++) {
		for (i = 0; i < size; i++) {
			r = (GLubyte) (255 * Tbuffer[j * size + i]);
			g = b = r;
			memDC2.SetPixelV(i, size-1-j, RGB(r, g, b));
		}
	}
	dc.BitBlt(0, size, size, size, &memDC2, 0, 0, SRCCOPY);
	//for (j = IMAGE_Y - 1; j >= 0; j--) {
	for (j = 0; j < size; j++) {
		for (i = 0; i < size; i++) {
			M00 += Tbuffer[j * size + i];
			M10 += i * Tbuffer[j * size + i];
			M01 += j * Tbuffer[j * size + i];
			M11 += i*j * Tbuffer[j * size + i];
			M20 += i*i * Tbuffer[j * size + i];
			M02 += j*j * Tbuffer[j * size + i];
		}
	}

	TRACE("M00 = %f, M10 = %f, M01 = %f, M11 = %f, M20 = %f, M02 = %f\n", M00, M10, M01, M11, M20, M02);
	xc = M10 / M00;
	yc = M01 / M00;
	A = M20 / M00 - xc*xc;
	B = 2 * (M11 / M00 - xc * yc);
	C = M02 / M00 - yc*yc;
	TRACE("A = %f, B = %f, C = %f\n", A, B, C);
	theta = atan(B / (A - C)) / 2;
	w = sqrt( 6 * (A + C - sqrt( B*B + (A-C)*(A-C) ) ) );
	l = sqrt( 6 * (A + C + sqrt( B*B + (A-C)*(A-C) ) ) );
	xc += (center_x - half_x);
	yc += (center_y - half_y);
	//TRACE("xc = %f, yc = %f, theta = %f, w = %f, l = %f\n", xc, yc, theta, w, l);
	TRACE("xc = %f, yc = %f, theta = %f, w = %f, l = %f\n", xc, yc, theta/PI * 180, w, l);
	stroke_ptr->SetValues(xc, yc, theta, w, l, (unsigned char)center_r, (unsigned char)center_g, (unsigned char)center_b);
}

void BlockMatch(int center_x, int center_y, int *new_x, int *new_y, int MASK_SIZE, 
				int search_radius)
{
	int x, y;
	unsigned char r, g, b;
	int half_x, half_y;
	int i, j, z;
	GLubyte *Tbuffer;
	double corr;
	double max_corr = -10.0;
	double E_mask[3], E_image[3];
	
	Tbuffer = new GLubyte [MASK_SIZE * MASK_SIZE * 3];
	half_x = half_y = (int) MASK_SIZE / 2;
	TRACE("half_x = half_y = %d\n", half_x);

	E_mask[0] = E_mask[1] = E_mask[2] = 0.0;
	//z = ((center_y-half_y) * IMAGE_X + (center_x-half_x)) * 3;
	for (j = 0; j < MASK_SIZE; j++) {
		z = ((center_y-half_y+j) * IMAGE_X + (center_x-half_x)) * 3;
		for (i = 0; i < MASK_SIZE; i++) {
			r = Dbuffer[z + i * 3 + 0];
			g = Dbuffer[z + i * 3 + 1];
			b = Dbuffer[z + i * 3 + 2];
			Tbuffer[(j * MASK_SIZE + i) * 3 + 0] = r;
			Tbuffer[(j * MASK_SIZE + i) * 3 + 1] = g;
			Tbuffer[(j * MASK_SIZE + i) * 3 + 2] = b;
			E_mask[0] += (double)r;
			E_mask[1] += (double)g;
			E_mask[2] += (double)b;
		}
	}
	E_mask[0] /= (MASK_SIZE * MASK_SIZE);
	E_mask[1] /= (MASK_SIZE * MASK_SIZE);
	E_mask[2] /= (MASK_SIZE * MASK_SIZE);
		
	glRasterPos2i(0, 0);
	glDrawPixels(MASK_SIZE, MASK_SIZE, GL_RGB, GL_UNSIGNED_BYTE, Tbuffer);
	swap();

	//LoadImage2("fightchar *filename, char *maskname, char *backname);

	TRACE("center_x = %d, center_y = %d\n", center_x, center_y);
	// 414, 378
	//for (y = 0; y < IMAGE_Y; y++) {
	//for (y = 350; y < 400; y++) {
	for (y = center_y - search_radius; y < center_y + search_radius; y++) {
		//for (x = 0; x < IMAGE_X; x++) {
		//for (x = 400; x < 450; x++) {
		for (x = center_x - search_radius; x < center_x + search_radius; x++) {
			//r = Dbuffer[(y * IMAGE_X + x) * 3 + 0];
			//g = Dbuffer[(y * IMAGE_X + x) * 3 + 1];
			//b = Dbuffer[(y * IMAGE_X + x) * 3 + 2];
			if (y - half_y < 0 || y + half_y >= IMAGE_Y) continue;
			if (x - half_x < 0 || x + half_x >= IMAGE_X) continue;
			corr = 0;
			E_image[0] = E_image[1] = E_image[2] = 0.0;
			//TRACE("x = %d, y = %d\n", x, y);
			for (j = 0; j < MASK_SIZE; j++) {
				for (i = 0; i < MASK_SIZE; i++) {
					z = ((y-half_y+j) * IMAGE_X + (x-half_x+i)) * 3 + 0;
					E_image[0] += Maskbuffer[z+0]; 
					E_image[1] += Maskbuffer[z+1]; 
					E_image[2] += Maskbuffer[z+2]; 
				}
			}
			//glRasterPos2i(0, 0);
			//glDrawPixels(MASK_SIZE, MASK_SIZE, GL_RGB, GL_UNSIGNED_BYTE, Testbuffer);
			//swap();
			E_image[0] /= (MASK_SIZE * MASK_SIZE);
			E_image[1] /= (MASK_SIZE * MASK_SIZE);
			E_image[2] /= (MASK_SIZE * MASK_SIZE);
			//TRACE("E_mask[0] = %f\n", E_mask[0]);
			//TRACE("E_mask[1] = %f\n", E_mask[1]);
			//TRACE("E_mask[2] = %f\n", E_mask[2]);
			//TRACE("E_image[0] = %f\n", E_image[0]);
			//TRACE("E_image[1] = %f\n", E_image[1]);
			//TRACE("E_image[2] = %f\n", E_image[2]);
			for (j = 0; j < MASK_SIZE; j++) {
				//if (y - half_y + j < 0 || y - half_y + j >= IMAGE_Y) continue;
				for (i = 0; i < MASK_SIZE; i++) {
					//if (x - half_x + i < 0 || x - half_x + i >= IMAGE_X) continue;
					z = ((y-half_y+j) * IMAGE_X + (x-half_x+i)) * 3 + 0;
					//corr += 
					//  (Tbuffer[(j*MASK_SIZE+i)*3+0]-E_mask[0]) * (Dbuffer[z+0]-E_image[0])
					//+ (Tbuffer[(j*MASK_SIZE+i)*3+1]-E_mask[1]) * (Dbuffer[z+1]-E_image[1]) 
					//+ (Tbuffer[(j*MASK_SIZE+i)*3+2]-E_mask[2]) * (Dbuffer[z+2]-E_image[2]); 
					corr += 
					  (Tbuffer[(j*MASK_SIZE+i)*3+0]-E_mask[0]) * (Maskbuffer[z+0]-E_image[0])
					+ (Tbuffer[(j*MASK_SIZE+i)*3+1]-E_mask[1]) * (Maskbuffer[z+1]-E_image[1]) 
					+ (Tbuffer[(j*MASK_SIZE+i)*3+2]-E_mask[2]) * (Maskbuffer[z+2]-E_image[2]); 
				}
			}
			if (corr > max_corr) {
				max_corr = corr; 
				*new_x = x;
				*new_y = y;
			}
		}
	}

	delete [] Tbuffer;

	glRasterPos2i(0, 0);
	glDrawPixels(IMAGE_X, IMAGE_Y, GL_RGB, GL_UNSIGNED_BYTE, Maskbuffer);
	swap();


	TRACE("center_x = %d, center_y = %d\n", center_x, center_y);
	TRACE("new_x = %d, new_y = %d\n", *new_x, *new_y);
}
*/
