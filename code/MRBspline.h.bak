#ifndef _MRBSPLINE_H_
#define _MRBSPLINE_H_

#include "vector.h"
#include "matrix.h"

#include <deque>
using namespace std;

#include "gl\gl.h"
#include "gl\glu.h"
#include "defines2.h"
#include "defines3.h"
#include "types.h"
#include "Field.h"
#include "MyList.h"

/*
struct PixeL {
	short x, y;
};
*/


//#define MAX_MR_LEVEL 10
#define MAX_MR_LEVEL 1

class MRBspline {
private:
	//IBSNode* pts; // 1D array
	//IBSNode *st_pt, *end_pt;
	//int	image_x, image_y; // size of the image
	//int cur_min_index; // current min_node index
	//int total_pts;
	//matrix A(5, 5); // Create 5x5 matrix
	int polyline_count;
	
	///////////////////////
	
	int max_level;
	vector knot; // knot vector
	//matrix cnt_pnts; // control points
	vector cnt_x, cnt_y; // control points
	vector tan_x, tan_y; // tangent points
	vector C0_x, C0_y; // control points at level 0
	vector *D_x, *D_y; // Difference points
	vector mr_cnt_x[MAX_MR_LEVEL], mr_cnt_y[MAX_MR_LEVEL];

public:
	int cur_level; // current display level
	int stroke_mode;
	//deque<PixeL> pnts;
	deque<double2D> pnts; // original data points
	deque<double2D> sample_pnts; // selected points for display
	vector h, t; // head and tail direction (unit vectors)
	////////////////////////////////////
	MyListNode<MRBspline>* node; // used in line drawing project!
	///////////////////////////////////
	MRBspline(); 
	//~MRBspline() { }
	//IBSnake(int x, int y);
	void DrawHermiteCurve(CDC& dc, int r, int g, int b);
	void SetPolylineCount(int count) { polyline_count = count; }
	int GetPolylineCount() { return polyline_count; }
	void AddDataPointRear(double x, double y);
	void AddDataPointHead(double x, double y);
	void FindCurve(); // Compute fitting B-spline control points from a sequence of data points
	inline double	N(int k, int m, double t, vector& knot);
	//double	N(int k, int m, double t, vector& knot);
	void DrawCurve(CDC *dc, int r, int g, int b);
	void DrawCurve2(CDC *dc, int r, int g, int b);
	void DrawInterpCurve(CDC *dc, int r, int g, int b);
	//void DrawBezierSegments(CDC& dc, int r, int g, int b);
	void DrawBezierSegments(CDC& dc, int r, int g, int b, int pen_size);
	void DrawMRCurve(CDC *dc, int level, int r, int g, int b);
	void DrawBlendedCurve(CDC *dc, int level1, int level2, double t, int r, int g, int b);
	void DrawCntPnts(CDC *dc, int r, int g, int b);
	void DrawMRCntPnts(CDC *dc, int level, int r, int g, int b);
	void SetKnotVector(vector& in);
	void SetMRCntPnts(int level, vector& in_x, vector& in_y);
	void SetCntPnts(vector& in_x, vector& in_y);
	void SetCntPnts(deque<double2D>& pnts);
	void SetCntPnts(); 
	void SetDataPnts(deque<double2D>& pnts);
	void SetSamplePnts(deque<double2D>& pnts);
	void SetInterpCntPnts();
	void SetHeadTailDir();
	void SetHeadDir(); 
	void SetTailDir(); 
	void SetTangents(deque<double2D>& pnts);
	void SetTangents(deque<double2D>& pnts, Field& gfield);
	vector& GetCnt0Pnts_x() { return C0_x; }
	vector& GetCnt0Pnts_y() { return C0_y; }
	vector& GetCntPnts_x() { return cnt_x; }
	vector& GetCntPnts_y() { return cnt_y; }
	vector& GetDiffPnts_x(int i) { return D_x[i]; }
	vector& GetDiffPnts_y(int i) { return D_y[i]; }
	void SetCnt0Pnts(vector& in_x, vector& in_y);
	void SetDiffPnts(int j, vector* D_x, vector* D_y);
	void SetMaxLevel(int j) { max_level = j; }
	int GetMaxLevel() { return max_level; }
	void SetCurLevel(int j) { cur_level = j; }
	int GetCurLevel() { return cur_level; }
	int FindCntPnt(int x, int y); // find a control point that is closest to the (x, y) position
	int FindMRCntPnt(int level, int x, int y);
	void SetCntPnt(int i, double x, double y) { cnt_x[i] = x; cnt_y[i] = y; }
	void SetMRCntPnt(int level, int i, double x, double y) { mr_cnt_x[level][i] = x; mr_cnt_y[level][i] = y; }
	void CopyLevel(int i, MRBspline& s);
	void SetEditedDiffPnts(int j, vector* in_x, vector* in_y);
	void DrawCntPnt(CDC *dc, int i, int r, int g, int b);
	void DrawMRCntPnt(CDC *dc, int level, int i, int r, int g, int b);
	//void Display(CDC *dc, int r, int g, int b);
	//void Optimize();
	//doubleComputeEnergy();
	//void ClearAll();
	//void Insert(IBSNode*);
	//IBSNode* DeleteMin();
	//void DeleteItem(IBSNode*);
	//int IfEmpty();

	//int	C; // number of total elements - 1
};

#endif

