#include "stdafx.h"

#include <cmath>
#include "globals.h"
#include "Field.h"

/*
void Field::set(Image& image, Image& gradient) 
{
	int i, j;
	
	for (i = 0; i < Nr - 1; i++) { 
		for (j = 0; j < Nc - 1; j++) {
			p[i][j].gx = image[i+1][j] - (double)image[i][j];
			p[i][j].gy = image[i][j+1] - (double)image[i][j];
			p[i][j].mag = gradient[i][j] / (double)max_grad; // [0, 1]
		}
	}
	
	for (i = 0; i < Nr - 1; i++) {
		p[i][Nc - 1].gx = p[i][Nc - 2].gx;
		p[i][Nc - 1].gy = p[i][Nc - 2].gy;
		p[i][Nc - 1].mag = p[i][Nc - 2].mag;
	}
	
	for (j = 0; j < Nc - 1; j++) {
		p[Nr - 1][j].gx = p[Nr - 2][j].gx;
		p[Nr - 1][j].gy = p[Nr - 2][j].gy;
		p[Nr - 1][j].mag = p[Nr - 2][j].mag;
	}
	
	p[Nr - 1][Nc - 1].gx = ( p[Nr - 1][Nc - 2].gx + p[Nr - 2][Nc - 1].gx ) / 2;
	p[Nr - 1][Nc - 1].gy = ( p[Nr - 1][Nc - 2].gy + p[Nr - 2][Nc - 1].gy ) / 2;
	p[Nr - 1][Nc - 1].mag = ( p[Nr - 1][Nc - 2].mag + p[Nr - 2][Nc - 1].mag ) / 2;

}
*/

void Field::set(Image& image) 
// Sobel gradient version
{
	int i, j;
	double MAX_GRADIENT = -1.;
	//double MAX_VAL = -1.;
	//double MAX_VAL = 255.;
	double MAX_VAL = 1020.; // 255 + 2 * 255 + 255 (Maximum possible Sobel value)

	/*
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			if (image[i][j] > MAX_VAL) MAX_VAL = (double)image[i][j];
		}
	}
	TRACE("MAX_VAL = %f\n", MAX_VAL);
	*/
	
	for (i = 1; i < Nr - 1; i++) { 
		for (j = 1; j < Nc - 1; j++) {
			////////////////////////////////////////////////////////////////
			// Important!: the value of image intensity should be normalized to [0,1]
			p[i][j].gx = (image[i+1][j-1] + 2*(double)image[i+1][j] + image[i+1][j+1] 
				- image[i-1][j-1] - 2*(double)image[i-1][j] - image[i-1][j+1]) / MAX_VAL;
			p[i][j].gy = (image[i-1][j+1] + 2*(double)image[i][j+1] + image[i+1][j+1]
				- image[i-1][j-1] - 2*(double)image[i][j-1] - image[i+1][j-1]) / MAX_VAL;
			/////////////////////////////////////////////
			// Amplify!!!
			//p[i][j].gx = pow(p[i][j].gx, 2); 
			//p[i][j].gy = pow(p[i][j].gy, 2); 
			//p[i][j].gx = pow(p[i][j].gx, 2); 
			//TRACE("p[i][j].gx = %.1f\n", p[i][j].gx);
			//TRACE("p[i][j].gy = %.1f\n", p[i][j].gy);
			p[i][j].mag = sqrt(p[i][j].gx * p[i][j].gx + p[i][j].gy * p[i][j].gy);

			if (p[i][j].mag > MAX_GRADIENT) {
				MAX_GRADIENT = p[i][j].mag;
				//TRACE("p[i][j].gx = %01.f\n", p[i][j].gx);
				//TRACE("p[i][j].gy = %01.f\n", p[i][j].gy);
				//TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);
			}
		}
	}

	for (i = 1; i <= Nr - 2; i++) {
		p[i][0].gx = p[i][1].gx;
		p[i][0].gy = p[i][1].gy;
		p[i][0].mag = p[i][1].mag;
		p[i][Nc - 1].gx = p[i][Nc - 2].gx;
		p[i][Nc - 1].gy = p[i][Nc - 2].gy;
		p[i][Nc - 1].mag = p[i][Nc - 2].mag;
	}
	
	for (j = 1; j <= Nc - 2; j++) {
		p[0][j].gx = p[1][j].gx;
		p[0][j].gy = p[1][j].gy;
		p[0][j].mag = p[1][j].mag;
		p[Nr - 1][j].gx = p[Nr - 2][j].gx;
		p[Nr - 1][j].gy = p[Nr - 2][j].gy;
		p[Nr - 1][j].mag = p[Nr - 2][j].mag;
	}
	
	p[0][0].gx = ( p[0][1].gx + p[1][0].gx ) / 2;
	p[0][0].gy = ( p[0][1].gy + p[1][0].gy ) / 2;
	p[0][0].mag = ( p[0][1].mag + p[1][0].mag ) / 2;
	p[0][Nc-1].gx = ( p[0][Nc-2].gx + p[1][Nc-1].gx ) / 2;
	p[0][Nc-1].gy = ( p[0][Nc-2].gy + p[1][Nc-1].gy ) / 2;
	p[0][Nc-1].mag = ( p[0][Nc-2].mag + p[1][Nc-1].mag ) / 2;
	p[Nr-1][0].gx = ( p[Nr-1][1].gx + p[Nr-2][0].gx ) / 2;
	p[Nr-1][0].gy = ( p[Nr-1][1].gy + p[Nr-2][0].gy ) / 2;
	p[Nr-1][0].mag = ( p[Nr-1][1].mag + p[Nr-2][0].mag ) / 2;
	p[Nr - 1][Nc - 1].gx = ( p[Nr - 1][Nc - 2].gx + p[Nr - 2][Nc - 1].gx ) / 2;
	p[Nr - 1][Nc - 1].gy = ( p[Nr - 1][Nc - 2].gy + p[Nr - 2][Nc - 1].gy ) / 2;
	p[Nr - 1][Nc - 1].mag = ( p[Nr - 1][Nc - 2].mag + p[Nr - 2][Nc - 1].mag ) / 2;

	TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);

	max_grad = MAX_GRADIENT;

	vector v(2);

	// Normalize gradient magnitude
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			///////////////////////
			p[i][j].mag = p[i][j].mag / (double)MAX_GRADIENT; // place it between [0, 1]
			//TRACE("[%d][%d] p.gx = %0.3f, p.gy = %0.3f\n", i, j, p[i][j].gx, p[i][j].gy);
			if (p[i][j].mag < 0)
				TRACE("p[%d][%d].mag = %0.2f\n", i, j, p[i][j].mag);
			/////////////////////////////////////////
			v[0] = p[i][j].gx;
			v[1] = p[i][j].gy;
			//v.print();
			v.make_unit();
			///////////////////////////////////////
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

void Field::set(imatrix& image) 
// Sobel gradient version
// you can use this imatrix version also
{
	int i, j;
	double MAX_GRADIENT = -1.;
	//double MAX_VAL = -1.;
	//double MAX_VAL = 255.;
	double MAX_VAL = 1020.; // 255 + 2 * 255 + 255 (Maximum possible Sobel value)

	/*
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			if (image[i][j] > MAX_VAL) MAX_VAL = (double)image[i][j];
		}
	}
	TRACE("MAX_VAL = %f\n", MAX_VAL);
	*/
	
	for (i = 1; i < Nr - 1; i++) { 
		for (j = 1; j < Nc - 1; j++) {
			////////////////////////////////////////////////////////////////
			// Important!: the value of image intensity should be normalized to [0,1]
			p[i][j].gx = (image[i+1][j-1] + 2*(double)image[i+1][j] + image[i+1][j+1] 
				- image[i-1][j-1] - 2*(double)image[i-1][j] - image[i-1][j+1]) / MAX_VAL;
			p[i][j].gy = (image[i-1][j+1] + 2*(double)image[i][j+1] + image[i+1][j+1]
				- image[i-1][j-1] - 2*(double)image[i][j-1] - image[i+1][j-1]) / MAX_VAL;
			/////////////////////////////////////////////
			// Amplify!!!
			//p[i][j].gx = pow(p[i][j].gx, 2); 
			//p[i][j].gy = pow(p[i][j].gy, 2); 
			//p[i][j].gx = pow(p[i][j].gx, 2); 
			//TRACE("p[i][j].gx = %.1f\n", p[i][j].gx);
			//TRACE("p[i][j].gy = %.1f\n", p[i][j].gy);
			p[i][j].mag = sqrt(p[i][j].gx * p[i][j].gx + p[i][j].gy * p[i][j].gy);

			if (p[i][j].mag > MAX_GRADIENT) {
				MAX_GRADIENT = p[i][j].mag;
				//TRACE("p[i][j].gx = %01.f\n", p[i][j].gx);
				//TRACE("p[i][j].gy = %01.f\n", p[i][j].gy);
				//TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);
			}
		}
	}

	for (i = 1; i <= Nr - 2; i++) {
		p[i][0].gx = p[i][1].gx;
		p[i][0].gy = p[i][1].gy;
		p[i][0].mag = p[i][1].mag;
		p[i][Nc - 1].gx = p[i][Nc - 2].gx;
		p[i][Nc - 1].gy = p[i][Nc - 2].gy;
		p[i][Nc - 1].mag = p[i][Nc - 2].mag;
	}
	
	for (j = 1; j <= Nc - 2; j++) {
		p[0][j].gx = p[1][j].gx;
		p[0][j].gy = p[1][j].gy;
		p[0][j].mag = p[1][j].mag;
		p[Nr - 1][j].gx = p[Nr - 2][j].gx;
		p[Nr - 1][j].gy = p[Nr - 2][j].gy;
		p[Nr - 1][j].mag = p[Nr - 2][j].mag;
	}
	
	p[0][0].gx = ( p[0][1].gx + p[1][0].gx ) / 2;
	p[0][0].gy = ( p[0][1].gy + p[1][0].gy ) / 2;
	p[0][0].mag = ( p[0][1].mag + p[1][0].mag ) / 2;
	p[0][Nc-1].gx = ( p[0][Nc-2].gx + p[1][Nc-1].gx ) / 2;
	p[0][Nc-1].gy = ( p[0][Nc-2].gy + p[1][Nc-1].gy ) / 2;
	p[0][Nc-1].mag = ( p[0][Nc-2].mag + p[1][Nc-1].mag ) / 2;
	p[Nr-1][0].gx = ( p[Nr-1][1].gx + p[Nr-2][0].gx ) / 2;
	p[Nr-1][0].gy = ( p[Nr-1][1].gy + p[Nr-2][0].gy ) / 2;
	p[Nr-1][0].mag = ( p[Nr-1][1].mag + p[Nr-2][0].mag ) / 2;
	p[Nr - 1][Nc - 1].gx = ( p[Nr - 1][Nc - 2].gx + p[Nr - 2][Nc - 1].gx ) / 2;
	p[Nr - 1][Nc - 1].gy = ( p[Nr - 1][Nc - 2].gy + p[Nr - 2][Nc - 1].gy ) / 2;
	p[Nr - 1][Nc - 1].mag = ( p[Nr - 1][Nc - 2].mag + p[Nr - 2][Nc - 1].mag ) / 2;

	TRACE("MAX_GRADIENT = %0.1f\n", MAX_GRADIENT);

	max_grad = MAX_GRADIENT;

	vector v(2);

	// Normalize gradient magnitude
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			//TRACE("p[%d][%d].mag = %0.3f\n", i, j, p[i][j].mag);
			///////////////////////
			p[i][j].mag = p[i][j].mag / (double)MAX_GRADIENT; // place it between [0, 1]
			if (p[i][j].mag < 0)
				TRACE("p[%d][%d].mag = %0.2f\n", i, j, p[i][j].mag);
			/////////////////////////////////////////
			v[0] = p[i][j].gx;
			v[1] = p[i][j].gy;
			//v.print();
			v.make_unit();
			//TRACE("[%d][%d] v[0] = %f, v[1] = %f\n", i, j, v[0], v[1]);
			p[i][j].gx = v[0];
			p[i][j].gy = v[1];
			//*/
			///////////////////////////////////////
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

void Field::normalize() 
// Make (gx, gy) a normal vector
{
	int i, j;
	
	vector v(2);

	// Normalize gradient magnitude
	for (i = 0; i < Nr; i++) { 
		for (j = 0; j < Nc; j++) {
			/////////////////////////////////////////
			v[0] = p[i][j].gx;
			v[1] = p[i][j].gy;
			//v.print();
			v.make_unit();
			//TRACE("[%d][%d] v[0] = %.1f, v[1] = %.1f\n", i, j, v[0], v[1]);
			p[i][j].gx = v[0];
			p[i][j].gy = v[1];
			//*/
			///////////////////////////////////////
		}
	}
}