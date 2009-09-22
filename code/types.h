#ifndef _TYPES_H_
#define _TYPES_H_

typedef int		Image[MAXWIDTH][MAXHEIGHT];
typedef short 	Image2[MAXWIDTH][MAXHEIGHT]; /* <- region_map */

typedef char	ColorImage[MAXWIDTH][MAXHEIGHT][3];

typedef struct PixeL {
   short x, y;
} PixeL;

typedef short	PxlImage[MAXWIDTH][MAXHEIGHT][2];

typedef float	CostImage[MAXWIDTH][MAXHEIGHT];

typedef int 	Pxl[2];

typedef struct Node {
		int	no_in_heap;
		PixeL loc;
   	PixeL coord;
   	double intensity;
   	double gradient;
   	double laplacian;
   	double cost;
   	int expanded;
   	PixeL N[9];
   	PixeL next;
} Node; /* used in intelligent scissors */

//typedef Node Nodes[MAX_MASK_SIZE][MAX_MASK_SIZE]; // for EL
typedef Node Nodes[MAXWIDTH][MAXHEIGHT]; // for NPR

typedef struct FaceNode	{
	int	no_in_heap;
	double	cost;
	int	expanded;
} FaceNode;

typedef short     Point_[2];

struct Region_ {
         int   col;
         int   max_cd;
         Point_	start;
         struct Region_ *child[MAX_CHILD_REGIONS];
         struct CurveData  *cdArray[MAX_CD_IN_REGION];
};

struct CurveData {
         int region_no;
         int   max_index;
         Point_   data[MAX_CURVE_PNTS];
};

typedef struct Region_     Region_;
typedef struct CurveData   CurveData;

typedef struct{
	float x, y;
} point2D;

typedef struct{
	float x, y, z;
} point3D;

typedef struct{
	int obj, pnt;
} ObjPnt;

typedef struct{
	float x, y, z, w;
} hCoord;

typedef struct {
	int leftmost, rightmost, topmost, bottommost;
} Box;

////////////////////////////
// Genetic painter
struct uday_data {
	double	x, y;
};

//struct ONE_STROKE {
//	double	x, y;
//};

struct ONE_STROKE { 
	double x;
	double y;
	GLubyte red;
	GLubyte green;
	GLubyte blue;
};

/*
struct edge_dist {
	double r_sum, r_mean, r_std;
	double g_sum, g_mean, g_std;
	double b_sum, b_mean, b_std;
	double gr_sum, gr_mean, gr_std;
};
*/
struct RGB {
	double r, g, b;
};

struct iRGB {
	//unsigned char r, g, b;
	GLubyte r, g, b;
};


struct XYRGB {
	double x, y, r, g, b;
};

struct double2D {
	double x, y;
};

//typedef deque<double2D> PntsDeque;

#endif



