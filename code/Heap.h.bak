#ifndef HEAP_H
#define HEAP_H

#include "stdafx.h"
//#include <stdio.h>
#include <cstdlib>
#include <cmath>
//#include <iostream>
//using namespace std;

//#include "TIP.h"
//#include "TIPDoc.h"
//#include "TIPView.h"
//#include "globals.h"

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

/*

extern Node *minHeap[MAX_ELEMENTS];

extern int 	TotalNodes = 0;
extern int	max_nodes = 0;
extern Image	in_heap;

extern int		xw, yw;
extern int		tex_xw, tex_yw;

extern int max_grad;

extern PixeL seeds[MAX_SEEDS];
extern int seed_count;
extern PixeL new_seed;

extern PixeL 	seed, Free, st_point, OldFree;
extern int		st_x, st_y, motif_exist;
extern int		s_x, s_y, e_x, e_y, o_x, o_y;
extern int		px, py;

extern Nodes 	item;

extern PxlImage	next_pxl;
extern PxlImage	cumulative_next_pxl;
extern CostImage	pxl_cost;

extern PixeL seed2;
extern Nodes 	item2;
extern PxlImage	next_pxl2;
extern CostImage	pxl_cost2;

extern Image		seed_map;
extern Image		boundary;
extern Image		tmp_img;
extern Image		boundary2;
extern Image		motif;
extern Image	 	laplacian;
extern Image	 	gradient;
extern Image	 	image;


// VARIABLE MASK SIZE 
extern int	MASK_SIZE = 49;
extern int OLD_MASK_SIZE = 49;

*/

template <class T> class Heap {
public:
	int TotalNodes, max_nodes;
	T* p; // p is a main array of the heap!
public:
	// Constructor
	Heap() {
		TotalNodes = max_nodes = 0;
	}
	// Destructor
	~Heap() {
		delete[] p;
	}
	
};

#endif


////////////////////////////////////////////////
///// Maintaining routines for Minheap of EL

/*
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
	// save the value of the element with the highest coord 
	item = minHeap[1];	// item <- root 
	// use last element in heap to adjust heap 
	temp = minHeap[TotalNodes--];
	parent = 1;
	child = 2;
	while (child <= TotalNodes) {
		// find the larget child of the current parent 
		if ((child < TotalNodes) && (minHeap[child]->cost > minHeap[child+1]->cost))
			child++;
		if (temp->cost <= minHeap[child]->cost) break;
		// move to the next lower level 
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

	if (no==1) return; // nothing to do 

	item = minHeap[no];

	child = no; // child >= 2 

	while (1) {
		parent = child/2;
		minHeap[parent]->no_in_heap = child;
		minHeap[child] = minHeap[parent];
		if (parent == 1) break;
		child = parent;
	}

	item->no_in_heap = 1;
	minHeap[1] = item; // raise the target item to root 
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
			return i; // location number 
	}
	return 0;
}

int No_(Node *node)
{
	int i;

	for(i=1; i<=TotalNodes; ++i) {
		// YOU CAN CHANGE THIS PART 
		if( minHeap[i] == node )
			return i;
	}
	//not in
	exit(1);
}

void print_node(Node *item)
{
	//printf("(%d, %d) : %f\n",item->coord.x, item->coord.y, item->cost);
}

void print_min_heap()
{
	int i;
	for(i=1; i<=TotalNodes; ++i) {
		//printf("%5.1f ",minHeap[i]->cost);
	}
	//printf("\n");

}
*/

/////////////////////////////////////////////////////////////
///// Image Handling routines



//////////////////////////////////////////////////////////////
/// Drawing routines

////////////////////////////////////////////////////////
/////////// EL main routines

////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

void ImagetoNodesNPR(int width, int height, Image image, int max_grad, Image gradient)
// this is for TIPP NPR version
{
	int	x, y;

	ClearMinHeap();

	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++) {
			// item[i][j] = (Node *)malloc(sizeof(Node));

			item[x][y].no_in_heap = 0;	// not in the heap yet 

			//if (i == 0 || i == size-1 || j == 0 || j == size-1 || x <= 0 || x >= width-1 || y <= 0 || y >= height-1) 
			//item[i][j].expanded = TRUE_;
			//else 
			item[x][y].expanded = FALSE_;

			// loc is the relative position in the mask 
			item[x][y].loc.x = x;   item[x][y].loc.y = y;
			// coord is the absolute position in the image 
			item[x][y].coord.x = x;  item[x][y].coord.y = y;
			item[x][y].next.x = x;  item[x][y].next.y = y;

			//item[x][y].cost = MAX_COST;
			item[x][y].cost = 255 - 255; // initialized as zero: maximum difference
			item[x][y].gradient = 1 - sqrt( (double)gradient[x][y] / (double)max_grad );
			item[x][y].laplacian = 0;

			item[x][y].N[1].x = 0;
			item[x][y].N[1].y = 0;
			item[x][y].N[2].x = 0;
			item[x][y].N[2].y = 0;
			item[x][y].N[3].x = 0;
			item[x][y].N[3].y = 0;
			item[x][y].N[4].x = 0;
			item[x][y].N[4].y = 0;
			item[x][y].N[5].x = 0;
			item[x][y].N[5].y = 0;
			item[x][y].N[6].x = 0;
			item[x][y].N[6].y = 0;
			item[x][y].N[7].x = 0;
			item[x][y].N[7].y = 0;
			item[x][y].N[8].x = 0;
			item[x][y].N[8].y = 0;

			item[x][y].intensity = image[x][y];
		}      
}

void ClearNodesNPR(int width, int height, Nodes item)
// TIPP NPR 
{
   int   i, j;

   TotalNodes = 0; // the total number of nodes in the heap now
   max_nodes = 0; // the maximum number of nodes reached during the heap growing and shrinking

   for (i = 0; i < width; i++)
      for (j = 0; j < height; j++) {
         item[i][j].next.x = i;  item[i][j].next.y = j;
         if (i == 0 || i == width-1 || j == 0 || j == height-1)
            item[i][j].expanded = TRUE_;
         else item[i][j].expanded = FALSE_;

         item[i][j].cost = MAX_COST;
         item[i][j].no_in_heap = 0;

         in_heap[i][j] = FALSE_;
      }
}

void ClearNodesNPR2(int width, int height, Nodes item)
{
   TotalNodes = 0;
   max_nodes = 0;

}

#define RGB_GETRED(rgb)    ((rgb) & 0xff) 
#define RGB_GETGREEN(rgb)    (((rgb) >> 8) & 0xff) 
#define RGB_GETBLUE(rgb)    (((rgb) >> 16) & 0xff)  
#define dist3(x, y, z) sqrt((x)*(x) + (y)*(y) + (z)*(z))
#define dist2(x1, y1, x2, y2) sqrt( ((x1)-(x2))*((x1)-(x2)) + ((y1)-(y2))*((y1)-(y2)) )

void InitNodesNPR(int width, int height)
// TIPP NPR 
{
	int   x, y;

	TotalNodes = 0; // the total number of nodes in the heap now
	max_nodes = 0; // the maximum number of nodes reached during the heap growing and shrinking

	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++) {
			//item[loc.x][loc.y].cost = 0.0;
			//q = &item[loc.x][loc.y];
			//pxl_cost[q->coord.x][q->coord.y] = 0.0;
			insert_min_heap(&item[x][y]);
			//item[i][j].expanded = FALSE_;
			
			//item[i][j].cost = MAX_COST;
			//item[i][j].no_in_heap = 0; // no_in_heap starts from 1, not 0: 0 means not in heap

			//in_heap[i][j] = FALSE_;
		}
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

*/

