#include <stdio.h>
#include <stdlib.h>

#include "stdafx.h"
#include <cmath>
#include <deque>
using namespace std;

#include "Cube.h"
#include "CubeDoc.h"
#include "CubeView.h"

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

#define ABS(x)	( ((x)>0) ? (x) : (-(x)) )

#define wZ	0.43
#define wD 	0.43
#define wG	0.14

Node *minHeap[MAX_ELEMENTS];

int 	TotalNodes = 0;
int	max_nodes = 0;
Image	in_heap;

int		xw, yw;
int		tex_xw, tex_yw;

int max_grad;

PixeL seeds[MAX_SEEDS];
int seed_count;
PixeL new_seed;

PixeL 	seed, Free, st_point, OldFree;
int		st_x, st_y, motif_exist;
int		s_x, s_y, e_x, e_y, o_x, o_y;
int		px, py;

Nodes 	item;

PxlImage	next_pxl;
PxlImage	cumulative_next_pxl;
CostImage	pxl_cost;

PixeL seed2;
Nodes 	item2;
PxlImage	next_pxl2;
CostImage	pxl_cost2;

//Image		seed_map;
Image		boundary;
Image		tmp_img;
Image		boundary2;
//Image		motif;
//Image	 	laplacian;
Image	 	gradient;
Image	 	image;

/* VARIABLE MASK SIZE */

int OLD_MASK_SIZE = 49;



////////////////////////////////////////////////
///// Maintaining routines for Minheap of EL

void ClearMinHeap()
{
		int	i;

		for (i = 0; i < MAX_ELEMENTS; i++)	{
			minHeap[i] = NULL;	
		}
}

void insert_min_heap(Node *item)
{
	int i;
	if (HEAP_FULL(TotalNodes)) {
		//The heap is full
	}
	i = ++TotalNodes; // starting from index 1, not 0
	if (TotalNodes >= max_nodes)	max_nodes = TotalNodes;
	while ( (i!=1) && (item->cost < minHeap[i/2]->cost) ) {
		minHeap[i/2]->no_in_heap = i;
		minHeap[i] = minHeap[i/2];

		i /= 2;
	}
	item->no_in_heap = i;
	minHeap[i] = item;

	in_heap[item->loc.x][item->loc.y] = TRUE_;
}

Node *delete_root()
{
	int parent, child;
	Node *item, *temp;
	if(HEAP_EMPTY(TotalNodes)) {
		//The heap is empty
	}
	/* save the value of the element with the highest coord */
	item = minHeap[1];	/* item <- root */
	/* use last element in heap to adjust heap */
	temp = minHeap[TotalNodes--];
	parent = 1;
	child = 2;
	while (child <= TotalNodes) {
		/* find the larget child of the current parent */
		if ((child < TotalNodes) && (minHeap[child]->cost > minHeap[child+1]->cost))
			child++;
		if (temp->cost <= minHeap[child]->cost) break;
		/* move to the next lower level */
		minHeap[child]->no_in_heap = parent;
		minHeap[parent] = minHeap[child];
		parent = child;
		child *= 2;
	}
	temp->no_in_heap = parent;
	minHeap[parent] = temp;

	in_heap[item->loc.x][item->loc.y] = FALSE_;
	item->no_in_heap = 0;

	return item;
}

void raise_to_root(int no)
{
	int parent, child;
	Node *item;

	if (no==1) return; /* nothing to do */

	item = minHeap[no];

	child = no; /* child >= 2 */

	while (1) {
		parent = child/2;
		minHeap[parent]->no_in_heap = child;
		minHeap[child] = minHeap[parent];
		if (parent == 1) break;
		child = parent;
	}

	item->no_in_heap = 1;
	minHeap[1] = item; /* raise the target item to root */
}

void delete_item(int no)
{
	raise_to_root(no);
	delete_root();
}

Node *extract_min()
{
	return( minHeap[1] );
}

int In(Node *node)
{
	int i;

	for(i=1; i<=TotalNodes; ++i) {
		if( minHeap[i] == node )
			return i; /* location number */
	}
	return 0;
}

int No_(Node *node)
{
	int i;

	for(i=1; i<=TotalNodes; ++i) {
		/* YOU CAN CHANGE THIS PART */
		if( minHeap[i] == node )
			return i;
	}
	//not in
	exit(1);
}

void print_node(Node *item)
{
	printf("(%d, %d) : %f\n",item->coord.x, item->coord.y, item->cost);
}

void print_min_heap()
{
	int i;
	for(i=1; i<=TotalNodes; ++i) 
		printf("%5.1f ",minHeap[i]->cost);
	printf("\n");
}

/////////////////////////////////////////////////////////////
///// Image Handling routines

void ClearImage(int width, int height, Image image)
{
	int	x, y;

	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++)
			image[x][y] = 0;
}

void	ClearPxlImage(int width, int height, PxlImage& image)
{
	int	x, y;

	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++) {
			image[x][y][0] = x;
			image[x][y][1] = y;
		}
}

void	ClearCostImage(int width, int height, CostImage& image)
{
	int	x, y;

	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++) {
			image[x][y] = MAX_COST;
		}
}

void GetImage(int width, int height, Image image, GLubyte *Dbuffer)
{
	GLubyte	*glPixel;
	int 	r, g, b;
	int 	min, max;
	int 	i, j;

	glPixel = Dbuffer;
	for (j=0; j < height; j++)
		for (i=0; i < width; i++) {
			r = *glPixel;	glPixel++;
			g = *glPixel;	glPixel++;
			b = *glPixel;	glPixel++;
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
			image[i][j] = (int) ((min + max) / 2.0);

		}
}

void GetGrayImage(int width, int height, imatrix& image, GLubyte *Dbuffer)
{
	GLubyte	*glPixel;
	int 	r, g, b;
	int 	min, max;
	int 	i, j;

	image.init(width, height);

	glPixel = Dbuffer;
	for (j=0; j < height; j++)
		for (i=0; i < width; i++) {
			r = *glPixel;	glPixel++;
			g = *glPixel;	glPixel++;
			b = *glPixel;	glPixel++;
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
			image[i][j] = (int) ((min + max) / 2.0);

		}
}



void GetColorImage(int width, int height, cimatrix& cmap, GLubyte *Dbuffer)
{
	GLubyte	*glPixel;
	int 	r, g, b;
	int 	i, j;

	cmap.init(width, height);

	glPixel = Dbuffer;
	for (j=0; j < height; j++)
		for (i=0; i < width; i++) {
			r = *glPixel;	glPixel++;
			g = *glPixel;	glPixel++;
			b = *glPixel;	glPixel++;
			cmap[i][j].r = (GLubyte)r;
			cmap[i][j].g = (GLubyte)g;
			cmap[i][j].b = (GLubyte)b;
		}
}

int getGradient(int width, int height, Image gradient, Image image)
{
	int	dx, dy;
	int	i, j;
	int MAX_GRADIENT = -1;

	for (j = 0; j < height - 1; j++)
		for (i = 0; i < width - 1; i++) {
			dx = image[i+1][j] - image[i][j];
			dy = image[i][j+1] - image[i][j];
			gradient[i][j] = dx * dx + dy * dy;

			if (gradient[i][j] > MAX_GRADIENT) MAX_GRADIENT = gradient[i][j];
		}

	for (i = 0; i < width - 1; i++)
		gradient[i][height - 1] = gradient[i][height - 2];

	for (j = 0; j < height - 1; j++)
		gradient[width - 1][j] = gradient[width - 2][j];

	gradient[width - 1][height - 1]
		= (gradient[width - 1][height - 2] 
			+ gradient[width - 2][height - 1]) / 2;

	return (MAX_GRADIENT);
}

int getSobelGradient(int width, int height, Image gradient, Image image)
{
	int	dx, dy;
	int	i, j;
	int MAX_GRADIENT = -1;

	for (j = 1; j <= height - 2; j++)
		for (i = 1; i <= width - 2; i++) {
			dx = image[i+1][j] - image[i][j];
			dy = image[i][j+1] - image[i][j];
			gradient[i][j] = dx * dx + dy * dy;

			if (gradient[i][j] > MAX_GRADIENT) MAX_GRADIENT = gradient[i][j];
		}

	for (i = 0; i < width - 1; i++) {
		gradient[i][0] = gradient[i][1];
		gradient[i][height - 1] = gradient[i][height - 2];
	}

	for (j = 0; j < height - 1; j++) {
		gradient[0][j] = gradient[1][j];
		gradient[width - 1][j] = gradient[width - 2][j];
	}

	gradient[0][0]	= (gradient[0][1]+gradient[1][0])/2;
	gradient[0][height-1]	= (gradient[0][height-2]+gradient[1][height-1])/2;
	gradient[width-1][0]	= (gradient[width-1][1]+gradient[width-2][0])/2;
	gradient[width-1][height-1]	= (gradient[width-1][height-2]+gradient[width-2][height-1])/2;

	return (MAX_GRADIENT);
}

int getBlurredGradient(int width, int height, Image& gradient, Image& image, int N)
{
	int	dx, dy;
	int	i, j;
	int MAX_GRADIENT = -1;

	Image tmp;
	
	/// copy image to tmp
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			tmp[i][j] = image[i][j];
		}
	}

	for (int k = 0; k < N; k++) { // Do the blurring N times 
		
		for (j = 1; j < height-1; j++) {
			for (i = 1; i < width-1; i++) {
				gradient[i][j] = round(
					(tmp[i-1][j-1]+tmp[i][j-1]+tmp[i+1][j-1]
					+tmp[i-1][j]+tmp[i][j]+tmp[i+1][j]
					+tmp[i-1][j+1]+tmp[i][j+1]+tmp[i+1][j+1]) 
					/ 9.0 );
				
			}
		}

		// Boundaries
		for (i = 0; i < width - 1; i++) {
			gradient[i][0] = gradient[i][1];
			gradient[i][height - 1] = gradient[i][height - 2];
		}
		for (j = 0; j < height - 1; j++) {
			gradient[0][j] = gradient[1][j];
			gradient[width-1][j] = gradient[width-2][j];
		}

		// copy gradient to tmp
		if (k != N-1) { // last step can be skipped
			for (j = 0; j < height - 1; j++)
				for (i = 0; i < width - 1; i++) {
					tmp[i][j] = gradient[i][j];
				}
		}
	}

	//////////////////////////////////////////////////////////
	// Get the gradient values
	for (j = 0; j < height - 1; j++)
			for (i = 0; i < width - 1; i++) {
				dx = tmp[i+1][j] - tmp[i][j];
				dy = tmp[i][j+1] - tmp[i][j];
				//gradient[i][j] = sqrt(dx * dx + dy * dy);
				gradient[i][j] = dx * dx + dy * dy;
				if (gradient[i][j] > MAX_GRADIENT) MAX_GRADIENT = gradient[i][j];
			}

	// Assign values on Boundaries
	for (i = 0; i < width - 1; i++)
		gradient[i][height - 1] = gradient[i][height - 2];

	for (j = 0; j < height - 1; j++)
		gradient[width - 1][j] = gradient[width - 2][j];

	gradient[width - 1][height - 1]
		= (gradient[width - 1][height - 2] 
			+ gradient[width - 2][height - 1]) / 2;
	
	return (MAX_GRADIENT);
}


Field gfield, gfield2;

void InitGradField(int image_x, int image_y, Field& gfield, Image& image)
{
	gfield.init(image_x, image_y); // allocate gradient field
	//gfield.set(image, gradient);
	gfield.set(image);
	
}

void InitGradField(int image_x, int image_y, Field& gfield, imatrix& gray)
{
	gfield.init(image_x, image_y); // allocate gradient field
	//gfield.set(image, gradient);
	gfield.set(gray);
	
}

void AlignGradient(int image_x, int image_y, Field& gfield, int size)
{
	int	x, y;
	int	i, j;
	Grad	sum;
	int N;
	
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			sum.gx = sum.gy = sum.mag = 0.0;
			N = 0;
			for (j = -size/2; j <= size/2; j++) {
				for (i = -size/2; i <= size/2; i++) {
					if (x+i >= 0 && x+i < image_x && y+j >= 0 && y+j < image_y) {
						N++;
						sum.gx += gfield[x+i][y+j].gx * gfield[x+i][y+j].mag;
						sum.gy += gfield[x+i][y+j].gy * gfield[x+i][y+j].mag;
						sum.mag += gfield[x+i][y+j].mag;
						//if (gfield[x+i][y+j].mag > 0.5)
						//	TRACE("mag[%d][%d] = %.2f\n", x+i, y+j, gfield[x+i][y+j].mag);
					}
				}
			}
			sum.mag /= (double)N;
			//if (sum.mag > gfield[x][y].mag) { 
			// Gradient magnitude got bigger! Now update changes
				gfield[x][y].mag = sum.mag;
				gfield[x][y].gx = sum.gx;
				gfield[x][y].gy = sum.gy;
			//}
		}
	}

}

void AlignGradient2(int image_x, int image_y, Field& gfield, int size)
// Make changes to each pixel independently using the 'source' information
{
	int	x, y;
	int	i, j;
	Grad	sum;
	int N;
	Field tmp;

	tmp = gfield;
	
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			sum.gx = sum.gy = sum.mag = 0.0;
			N = 0;
			for (j = -size/2; j <= size/2; j++) {
				for (i = -size/2; i <= size/2; i++) {
					if (x+i >= 0 && x+i < image_x && y+j >= 0 && y+j < image_y) {
						N++;
						//TRACE("mag[%d][%d] = %.2f\n", x+i, y+j, tmp[x+i][y+j].mag);
						sum.gx += tmp[x+i][y+j].gx * tmp[x+i][y+j].mag;
						sum.gy += tmp[x+i][y+j].gy * tmp[x+i][y+j].mag;
						sum.mag += tmp[x+i][y+j].mag;
						//if (tmp[x+i][y+j].mag > 0.5)
						//TRACE("mag[%d][%d] = %.2f\n", x+i, y+j, gfield[x+i][y+j].mag);
					}
				}
			}
			sum.mag /= (double)N;
			//if (sum.mag > gfield[x][y].mag) { 
			// Gradient magnitude got bigger! Now update changes
				gfield[x][y].mag = sum.mag;
				gfield[x][y].gx = sum.gx;
				gfield[x][y].gy = sum.gy;
			//}
		}
	}

}

void GVF(int image_x, int image_y, Field& gfield, int N)
{
	int	x, y;
	int	i;
	//int N = 150;
	//double mu = 0.5;
	double mu = 0.2;
	Field gvf, gvf_tmp;
	double	fx, fy, u, v;
	double	sqr_mag_f, lap_u, lap_v;
	//double THRES = 0.8;

	
	gvf.init(gfield.getRow(), gfield.getCol()); // init gvf
	gvf_tmp.init(gfield.getRow(), gfield.getCol()); // init gvf
	gvf.copy(gfield); // copy gfield to gvf_tmp
	gvf_tmp.copy(gfield); // copy gfield to gvf_tmp
		
	for (i = 0; i < N; i++) {
		for (y = 0; y < image_y; y++) {
			for (x = 0; x < image_x; x++) {
				//sum.gx = sum.gy = sum.mag = 0.0;
				fx = gfield[x][y].gx;
				fy = gfield[x][y].gy;
				u = gvf[x][y].gx;
				v = gvf[x][y].gy;
				sqr_mag_f = fx*fx + fy*fy; // This is fixed throughout the iteration
				////////////////////////////////////////////
				// Compute the laplacian of the current u and v
				lap_u = lap_v = 0;
				// lap_u
				if (x-1 >= 0) lap_u += gvf[x-1][y].gx;
				else lap_u += gvf[x][y].gx;
				if (x+1 <= image_x-1) lap_u += gvf[x+1][y].gx;
				else lap_u += gvf[x][y].gx;
				if (y-1 >= 0) lap_u += gvf[x][y-1].gx;
				else lap_u += gvf[x][y].gx;
				if (y+1 <= image_y-1) lap_u += gvf[x][y+1].gx;
				else lap_u += gvf[x][y].gx;
				lap_u -= 4*gvf[x][y].gx;
				// lap_v
				if (x-1 >= 0) lap_v += gvf[x-1][y].gy;
				else lap_v += gvf[x][y].gy;
				if (x+1 <= image_x-1) lap_v += gvf[x+1][y].gy;
				else lap_v += gvf[x][y].gy;
				if (y-1 >= 0) lap_v += gvf[x][y-1].gy;
				else lap_v += gvf[x][y].gy;
				if (y+1 <= image_y-1) lap_v += gvf[x][y+1].gy;
				else lap_v += gvf[x][y].gy;
				lap_v -= 4*gvf[x][y].gy;
				
				u = u + mu*lap_u - sqr_mag_f * (u - fx);
				v = v + mu*lap_v - sqr_mag_f * (v - fy);
				//u = u + mu*lap_u - sqr_mag_f * pow((u - fx), 2);
				//v = v + mu*lap_v - sqr_mag_f * pow((v - fy), 2);
				//if (i > 10) {
				//	TRACE("[%d][%d]u = %.1f\n", x, y, u);
				//	TRACE("[%d][%d]v = %.1f\n", x, y, v);
				//}
				
				//if ( sqrt(u*u + v*v) / gvf_tmp.GetMaxGrad() > gvf_tmp[x][y].mag ) {
				//	gvf_tmp[x][y].mag = sqrt(u*u + v*v) / gvf_tmp.GetMaxGrad();
					gvf_tmp[x][y].gx = u;
					gvf_tmp[x][y].gy = v;
				//}
				//TRACE("u = %0.1f, v = %0.1f\n", u, v);
				//TRACE("sqrt(u*u + v*v) = %0.1f\n", sqrt(u*u + v*v));
				//TRACE("gvf_tmp.GetMaxGrad() = %0.1f\n", gvf_tmp.GetMaxGrad());
				//TRACE("gvf_tmp[x][y].mag = %0.2f\n", gvf_tmp[0][0].mag);
				
			}
			
		}
		gvf.copy(gvf_tmp);
		//TRACE("[%d]: gvf_tmp[19][0].mag = %0.2f\n", i, gvf_tmp[19][0].mag);
	}
	gfield.copy(gvf);

	/*
	double sum2 = 0.0;
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//sum.gx = sum.gy = sum.mag = 0.0;
			sum2 += gfield[x][y].mag;
			//TRACE("gfield[%d][%d].mag = %0.2f\n", x, y, gfield[x][y].mag);
		}
	}
	TRACE("sum2 = %f\n", sum2);
	sum2 /= ((double)image_x * image_y);
	TRACE("sum = %f\n", sum2);
	*/
}

void GVF2(int image_x, int image_y, Field& gfield, int N, double mu)
// included mu as a parameter
{
	int	x, y;
	int	i;
	//int N = 150;
	//double mu = 0.5;
	//double mu = 0.2;
	Field gvf, gvf_tmp;
	double	fx, fy, u, v;
	double	sqr_mag_f, lap_u, lap_v;
	//double THRES = 0.8;
	
	gvf.init(gfield.getRow(), gfield.getCol()); // init gvf
	gvf_tmp.init(gfield.getRow(), gfield.getCol()); // init gvf
	gvf.copy(gfield); // copy gfield to gvf_tmp
	gvf_tmp.copy(gfield); // copy gfield to gvf_tmp
		
	for (i = 0; i < N; i++) {
		for (y = 0; y < image_y; y++) {
			for (x = 0; x < image_x; x++) {
				//sum.gx = sum.gy = sum.mag = 0.0;
				fx = gfield[x][y].gx;
				fy = gfield[x][y].gy;
				u = gvf[x][y].gx;
				v = gvf[x][y].gy;
				sqr_mag_f = fx*fx + fy*fy; // This is fixed throughout the iteration
				////////////////////////////////////////////
				// Compute the laplacian of the current u and v
				lap_u = lap_v = 0;
				// lap_u
				if (x-1 >= 0) lap_u += gvf[x-1][y].gx;
				else lap_u += gvf[x][y].gx;
				if (x+1 <= image_x-1) lap_u += gvf[x+1][y].gx;
				else lap_u += gvf[x][y].gx;
				if (y-1 >= 0) lap_u += gvf[x][y-1].gx;
				else lap_u += gvf[x][y].gx;
				if (y+1 <= image_y-1) lap_u += gvf[x][y+1].gx;
				else lap_u += gvf[x][y].gx;
				lap_u -= 4*gvf[x][y].gx;
				// lap_v
				if (x-1 >= 0) lap_v += gvf[x-1][y].gy;
				else lap_v += gvf[x][y].gy;
				if (x+1 <= image_x-1) lap_v += gvf[x+1][y].gy;
				else lap_v += gvf[x][y].gy;
				if (y-1 >= 0) lap_v += gvf[x][y-1].gy;
				else lap_v += gvf[x][y].gy;
				if (y+1 <= image_y-1) lap_v += gvf[x][y+1].gy;
				else lap_v += gvf[x][y].gy;
				lap_v -= 4*gvf[x][y].gy;
				
				u = u + mu*lap_u - sqr_mag_f * (u - fx);
				v = v + mu*lap_v - sqr_mag_f * (v - fy);
				//u = u + mu*lap_u - sqr_mag_f * pow((u - fx), 2);
				//v = v + mu*lap_v - sqr_mag_f * pow((v - fy), 2);
				//if (i > 10) {
				//	TRACE("[%d][%d]u = %.1f\n", x, y, u);
				//	TRACE("[%d][%d]v = %.1f\n", x, y, v);
				//}
				
				//if ( sqrt(u*u + v*v) / gvf_tmp.GetMaxGrad() > gvf_tmp[x][y].mag ) {
				//	gvf_tmp[x][y].mag = sqrt(u*u + v*v) / gvf_tmp.GetMaxGrad();
					gvf_tmp[x][y].gx = u;
					gvf_tmp[x][y].gy = v;
				//}
				//TRACE("u = %0.1f, v = %0.1f\n", u, v);
				//TRACE("sqrt(u*u + v*v) = %0.1f\n", sqrt(u*u + v*v));
				//TRACE("gvf_tmp.GetMaxGrad() = %0.1f\n", gvf_tmp.GetMaxGrad());
				//TRACE("gvf_tmp[x][y].mag = %0.2f\n", gvf_tmp[0][0].mag);
				
			}
			
		}
		gvf.copy(gvf_tmp);
		//TRACE("[%d]: gvf_tmp[19][0].mag = %0.2f\n", i, gvf_tmp[19][0].mag);
	}
	gfield.copy(gvf);

	/*
	double sum2 = 0.0;
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//sum.gx = sum.gy = sum.mag = 0.0;
			sum2 += gfield[x][y].mag;
			//TRACE("gfield[%d][%d].mag = %0.2f\n", x, y, gfield[x][y].mag);
		}
	}
	TRACE("sum2 = %f\n", sum2);
	sum2 /= ((double)image_x * image_y);
	TRACE("sum = %f\n", sum2);
	*/
}

void GVF3(int image_x, int image_y, Field& gfield, int N, double factor)
// Adaptive GVF
{
	int	x, y;
	int	i;
	//int N = 150;
	//double mu = 0.5;
	//double mu = 0.2;
	Field gvf, gvf_tmp;
	double	fx, fy, u, v;
	double	sqr_mag_f, lap_u, lap_v;
	double mu;
	//double THRES = 0.8;

	gvf.init(gfield.getRow(), gfield.getCol()); // init gvf
	gvf_tmp.init(gfield.getRow(), gfield.getCol()); // init gvf
	gvf.copy(gfield); // copy gfield to gvf_tmp
	gvf_tmp.copy(gfield); // copy gfield to gvf_tmp
		
	for (i = 0; i < N; i++) {
		for (y = 0; y < image_y; y++) {
			for (x = 0; x < image_x; x++) {
				//sum.gx = sum.gy = sum.mag = 0.0;
				fx = gfield[x][y].gx;
				fy = gfield[x][y].gy;
				u = gvf[x][y].gx;
				v = gvf[x][y].gy;
				sqr_mag_f = fx*fx + fy*fy; // This is fixed throughout the iteration
				////////////////////////////////////////////
				// Compute the laplacian of the current u and v
				lap_u = lap_v = 0;
				// lap_u
				if (x-1 >= 0) lap_u += gvf[x-1][y].gx;
				else lap_u += gvf[x][y].gx;
				if (x+1 <= image_x-1) lap_u += gvf[x+1][y].gx;
				else lap_u += gvf[x][y].gx;
				if (y-1 >= 0) lap_u += gvf[x][y-1].gx;
				else lap_u += gvf[x][y].gx;
				if (y+1 <= image_y-1) lap_u += gvf[x][y+1].gx;
				else lap_u += gvf[x][y].gx;
				lap_u -= 4*gvf[x][y].gx;
				// lap_v
				if (x-1 >= 0) lap_v += gvf[x-1][y].gy;
				else lap_v += gvf[x][y].gy;
				if (x+1 <= image_x-1) lap_v += gvf[x+1][y].gy;
				else lap_v += gvf[x][y].gy;
				if (y-1 >= 0) lap_v += gvf[x][y-1].gy;
				else lap_v += gvf[x][y].gy;
				if (y+1 <= image_y-1) lap_v += gvf[x][y+1].gy;
				else lap_v += gvf[x][y].gy;
				lap_v -= 4*gvf[x][y].gy;
				
				//u = u + mu*lap_u - sqr_mag_f * (u - fx);
				//v = v + mu*lap_v - sqr_mag_f * (v - fy);
				/////////////////////////////////////////////
				// Adaptive!
				mu = ((double)gray2[x][y] / 255.0) * factor + 0.01;
				u = u + mu*lap_u - sqr_mag_f * (u - fx);
				v = v + mu*lap_v - sqr_mag_f * (v - fy);
				//////////////////////////////////////////
				//u = u + mu*lap_u - sqr_mag_f * pow((u - fx), 2);
				//v = v + mu*lap_v - sqr_mag_f * pow((v - fy), 2);
				//if (i > 10) {
				//	TRACE("[%d][%d]u = %.1f\n", x, y, u);
				//	TRACE("[%d][%d]v = %.1f\n", x, y, v);
				//}
				
				//if ( sqrt(u*u + v*v) / gvf_tmp.GetMaxGrad() > gvf_tmp[x][y].mag ) {
				//	gvf_tmp[x][y].mag = sqrt(u*u + v*v) / gvf_tmp.GetMaxGrad();
					gvf_tmp[x][y].gx = u;
					gvf_tmp[x][y].gy = v;
				//}
				//TRACE("u = %0.1f, v = %0.1f\n", u, v);
				//TRACE("sqrt(u*u + v*v) = %0.1f\n", sqrt(u*u + v*v));
				//TRACE("gvf_tmp.GetMaxGrad() = %0.1f\n", gvf_tmp.GetMaxGrad());
				//TRACE("gvf_tmp[x][y].mag = %0.2f\n", gvf_tmp[0][0].mag);
				
			}
			
		}
		gvf.copy(gvf_tmp);
		//TRACE("[%d]: gvf_tmp[19][0].mag = %0.2f\n", i, gvf_tmp[19][0].mag);
	}
	gfield.copy(gvf);

	/*
	double sum2 = 0.0;
	for (y = 0; y < image_y; y++) {
		for (x = 0; x < image_x; x++) {
			//sum.gx = sum.gy = sum.mag = 0.0;
			sum2 += gfield[x][y].mag;
			//TRACE("gfield[%d][%d].mag = %0.2f\n", x, y, gfield[x][y].mag);
		}
	}
	TRACE("sum2 = %f\n", sum2);
	sum2 /= ((double)image_x * image_y);
	TRACE("sum = %f\n", sum2);
	*/
}

void getLaplacian(int width, int height, Image laplacian, Image image)
{
	int	i, j;

	for (j = 1; j < height - 1; j++)
		for (i = 1; i < width - 1; i++) {
			laplacian[i][j] = 4 * image[i][j] -	image[i][j-1] - image[i-1][j] - image[i+1][j] - image[i][j+1];
		}

	for (i = 1; i < width - 1; i++) {
		laplacian[i][0] = laplacian[i][1];
		laplacian[i][height - 1] = laplacian[i][height - 2];
	}

	for (j = 1; j < height - 1; j++) {
		laplacian[0][j] = laplacian[1][j];
		laplacian[width - 1][j] = laplacian[width - 2][j];
	}

	laplacian[0][0] = (laplacian[0][1] + laplacian[1][0]) / 2;
	laplacian[0][height - 1] 
		= (laplacian[0][height - 2] + laplacian[1][height - 1]) / 2;
	laplacian[width - 1][0] 
		= (laplacian[width - 1][1] + laplacian[width - 2][0]) / 2;
	laplacian[width - 1][height - 1]
		= (laplacian[width - 1][height - 2] + laplacian[width - 2][height - 1]) / 2;
}


//////////////////////////////////////////////////////////////
/// Drawing routines

void DrawPoint(int cur_x, int cur_y, int size)
{
	long	vert[2];

	glPointSize((float)size);

	vert[0] = cur_x;	vert[1] = cur_y;
	glBegin(GL_POINTS);
		glVertex2i(vert[0], vert[1]);
	glEnd();
}

void DrawPointDC(CClientDC *dc, int x, int y, int size, int r, int g, int b)
{
	CPen pen;
		
	pen.CreatePen(PS_SOLID, size, RGB(r, g, b));
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);
	dc->MoveTo(x, y);
	dc->LineTo(x, y);
	dc->SelectObject(pOldPen);
}

void DrawGridLinesDC(CClientDC *dc, int y, int size, int r, int g, int b)
{
	int	i;
	
	CPen pen;
		
	pen.CreatePen(PS_SOLID, size, RGB(r, g, b));
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	if(fg_num)
	for(i=1; i<=9; i++){
		dc->MoveTo(i*(IMAGE_X/10), 0);
		dc->LineTo(i*(IMAGE_X/10), y);
		dc->LineTo((int)(i*(IMAGE_X*3/10)-IMAGE_X), (int)(IMAGE_Y-1));
		
	}
	for(i=1; i<=4; i++){
		dc->MoveTo(0, i*(y/5));
		dc->LineTo(IMAGE_X-1, i*(y/5));
		
	}
	
	for(i=1; i<=4; i++){
		
		dc->MoveTo(0, (int)(y + ((IMAGE_Y-y)/20.)*i*i) );
		dc->LineTo(IMAGE_X-1, (int)(y + ((IMAGE_Y-y)/20.)*i*i) );
		
	}
	dc->SelectObject(pOldPen);
}

void DrawPointMemDC(CDC *dc, int x, int y, int size, int r, int g, int b)
{
	CPen pen;
		
	pen.CreatePen(PS_SOLID, size, RGB(r, g, b));
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);
	dc->MoveTo(x, y);
	dc->LineTo(x, y);
	dc->SelectObject(pOldPen);
}

void DrawEdgeMemDC(CDC *dc, PixeL free, PixeL seed, int r, int g, int b)
{
	short redvec[3] = {255, 0, 0}; 
	int 	vert[2], old[2];

	CPen pen;
		
	pen.CreatePen(PS_SOLID, 2, RGB(r, g, b));
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	old[0] = free.x;
	old[1] = free.y;
	while(1) {
		vert[0] = next_pxl[old[0]][old[1]][0];
		vert[1] = next_pxl[old[0]][old[1]][1];
		dc->MoveTo(old[0], IMAGE_Y-1-old[1]);
		dc->LineTo(vert[0], IMAGE_Y-1-vert[1]);
		old[0] = vert[0];
		old[1] = vert[1];
		if (vert[0] == seed.x && vert[1] == seed.y) break;
	}

	dc->SelectObject(pOldPen);
}


void DrawCumulativeEdgeMemDC(CDC *dc, PixeL free, PixeL seed, int r, int g, int b)
{
	short redvec[3] = {255, 0, 0}; 
	long 	vert[2];
	long 	old[2];

	CPen pen;
		
	pen.CreatePen(PS_SOLID, 2, RGB(r, g, b));
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	old[0] = seed.x;
	old[1] = seed.y;
	while(1) {
		vert[0] = cumulative_next_pxl[old[0]][old[1]][0];
		vert[1] = cumulative_next_pxl[old[0]][old[1]][1];
		dc->MoveTo(old[0], IMAGE_Y-1-old[1]);
		dc->LineTo(vert[0], IMAGE_Y-1-vert[1]);
		old[0] = vert[0];
		old[1] = vert[1];
		if (vert[0] == free.x && vert[1] == free.y) break;
	}

	dc->SelectObject(pOldPen);
}

////////////////////////////////////////////////////////////////
// Draw Uniform Cubic B-spline from free-hand curve

inline bool even(int a)
{
  return ((a%2) == 0); // it returns 1 or 0
}





////////////////////////////////////////////////////////////////

void DrawSeedsMemDC(CDC *dc, int sc, int r, int g, int b)
{
	int k;

	for (k = 0; k < sc; k++) {
		DrawPointMemDC(dc, seeds[k].x, IMAGE_Y-1-seeds[k].y, FIS_POINT_SIZE, r, g, b);
	}
}

void DrawBoxMemDC(CDC *dc, int size, int x, int y, int r, int g, int b)
{
	long 	vert[4][2];
	int	half_x, half_y;

	half_x = half_y = size/2;

	CPen pen;
		
	pen.CreatePen(PS_SOLID, 2, RGB(r, g, b));
	CPen *pOldPen = (CPen *)dc->SelectObject(&pen);

	vert[0][0] = x - half_x;
	vert[0][1] = y - half_y;
	vert[1][0] = x + half_x;
	vert[1][1] = y - half_y;
	vert[2][0] = x + half_x;
	vert[2][1] = y + half_y;
	vert[3][0] = x - half_x;
	vert[3][1] = y + half_y;

	dc->MoveTo(vert[0][0], vert[0][1]);
	dc->LineTo(vert[1][0], vert[1][1]);
	dc->LineTo(vert[2][0], vert[2][1]);
	dc->LineTo(vert[3][0], vert[3][1]);
	dc->LineTo(vert[0][0], vert[0][1]);

	dc->SelectObject(pOldPen);


}

inline double norm2(double vx, double vy)
{
	return sqrt(vx*vx + vy*vy);
}

void ClearMemDC(CDC *dc)
{
	int	x, y;

	for (y = 0; y < IMAGE_Y; y++) {
		for (x = 0; x < IMAGE_X; x++) {
			/// Set Pixel in MemDC
			dc->SetPixelV(x, y, RGB(255, 255, 255));
		}
	}
}

