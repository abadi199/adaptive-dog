#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include "gl\gl.h"
#include "gl\glu.h"

#include "defines2.h"
#include "defines3.h"
#include "types.h"

#include "IPen.h"
#include "MRBspline.h"
#include "myvec.h"
#include "ivectormatrix.h"
#include "cimatrix.h"
#include "Field.h"
#include "ETF.h"
#include "fdog.h"
#include "edge_dist.h"
#include "pixel_dist.h"
#include "SimpleList.h"
#include "MyList.h"
#include "qmatrix.h"
#include <time.h>

#include <deque>
#include <algorithm>
using namespace std;

//#include "strokes.h"

#define ON_PLAY 0
#define OFF_PLAY 1
#define FORWARD 1
#define BACKWARD 2
#define ON_VIDEO 0
#define OFF_VIDEO 1

extern CDC	memDC;
extern HDC	HmemDC;
extern CDC			memDC2; // foreground texture
extern CDC			memDC3; // foreground mask
extern CDC			memDC4; // second image
extern CDC			double_buffer; // double buffer
extern CBitmap		bitmap;
extern CBitmap		bitmap2; // double buffer
extern CBitmap		bitmap3;

extern int	image_mode;	

extern int FOCAL;

extern int cur_fg;	// currently specifying foreground ojbect number

extern int IMAGE_X, IMAGE_Y;

extern PxlImage	next_pxl;
extern PxlImage	cumulative_next_pxl;
extern CostImage	pxl_cost;

extern Image		seed_map;
extern Image		boundary;
extern Image		boundary2;
extern Image		tmp_img;

extern Image		motif;
extern Image	 	image;
extern Image	 	laplacian;
extern Image	 	gradient;

extern int	OLD_MASK_SIZE;

extern int status, status2;
extern int max_grad;
extern PixeL 	seed, Free, st_point, OldFree;
extern int		st_x, st_y, motif_exist;
extern int		s_x, s_y, e_x, e_y, o_x, o_y;
extern int		px, py;

extern PixeL seeds[MAX_SEEDS];
extern int seed_count;
extern PixeL new_seed;

extern Nodes 	item;

extern PixeL seed2;
extern Nodes	item2;
extern PxlImage	next_pxl2;
extern CostImage	pxl_cost2;

// From file.cpp

extern int file_loaded;
extern int motif_exist;

extern float x_size, y_size;
extern char prefix[256];
extern char sys_msg[256];

extern int IMAGE_X, IMAGE_Y;

extern int snapped; 

extern GLuint texName[MAX_TIP_FG_NUM];

extern char fgText[MAX_TIP_FG_NUM][255];
extern char fg_maskText[MAX_TIP_FG_NUM][255];
extern char fg_maskTextTIV[MAX_TIP_FG_NUM][MAX_FRAMES][255];
extern point2D fgSize[MAX_TIP_FG_NUM];
extern GLubyte *fgImage[MAX_TIP_FG_NUM];

extern int fg_num;
extern int frame_num;
extern int fis_frame_num;
extern GLubyte BGtexture[BACK_TEXTURE_WIDTH][TEXTURE][4];
extern GLubyte FGtexture[MAX_TIP_FG_NUM][TEXTURE][TEXTURE][4];

extern GLubyte *Dbuffer;
extern GLubyte *Maskbuffer;
extern GLubyte *Backbuffer;
extern GLubyte *Membuffer;
extern GLubyte *Doublebuffer;

#define MIN(A, B) (((A) < (B))? (A) : (B))
#define MAX(A, B) (((A) > (B))? (A) : (B))

extern int cur_t;
extern int play_mode;
extern int video_mode;
extern int play_direction;
// 
extern int fg_num;

extern int scene_model_type;
extern int panoramic_view_created;

extern int		src_x, src_y; // the view-starting coordinates of the source bitmap
extern float		displace_x, displace_y; // the displacement in navigating
extern GLUquadricObj	*qobj;
extern GLuint texture_name[10];
extern GLuint scene_model_texture_name;
extern GLubyte *background_texture;
extern GLubyte *scene_model_texture;
extern float	model_rotate_x[MAX_TRANSFORM], model_rotate_y[MAX_TRANSFORM];
extern float	model_translate_z[MAX_TRANSFORM];
extern float	model_translate_z_single;
extern float	cumulative_model_rotate_y;
extern int		max_transform_num, transform_toggle, translate_direction;
extern float	aspect;
extern GLuint	startList;
extern float	focal_length, fovy;
extern float	sphere_radius;
extern float	vl_latitude;
extern int		vl_index, vl_array_index;
extern int		vanish_y;
extern double	vl_phi;
extern float	fg_x[10][5], fg_y[10][5], fg_z[10][5]; // foreground objects 3D coordinates
extern int longitude_size, latitude_size;

extern point2D fg[MAX_TIP_FG_NUM][MAX_FRAMES][11];
extern point2D TC[MAX_TIP_FG_NUM][MAX_FRAMES][4];

extern int stroke_mode;

//extern double min_obj_value;

//////////////////////////////////////////////////////////////////////////////

extern void OpenImage(char *filename, char *maskname, char *backname);

extern void swap(void);
extern void renderReshape(int width, int height);
extern void RGBtoHSV( float r, float g, float b, float *h, float *s, float *v );

extern void	ClearImage(int width, int height, Image image);
extern void	ClearPxlImage(int width, int height, PxlImage& image);
extern void	ClearCostImage(int width, int height, CostImage& image);
extern void ClearNodes(int width, int height, Nodes item);
extern void MasktoNodes(int width, int height, int size, PixeL seed, 
		Image& image, int max_grad, Image& gradient, Image& laplacian, 
		Nodes& item);
extern void Intelligent_Scissor(int size, int width, int height, Nodes& item);
extern void SaveCurrentRegionToPreviousAndUpdateCumulative(int xw, int yw, 
	int size, PixeL free, Image& previous, Image& cumulative);
extern void GetImage(int width, int height, Image image, GLubyte *Dbuffer);
extern int getGradient(int width, int height, Image gradient, Image image);

extern void swap(void);
extern void DrawPoint(int, int, int);
extern void DrawPointDC(CClientDC *, int, int, int, int, int, int);
extern void DrawBoxDC(CClientDC *dc, int size, int x, int y, int r, int g, int b);
extern void DrawEdgeDC(CClientDC *dc, PixeL, PixeL, int r, int g, int b);
extern void DrawEdgeDC2(CClientDC *dc, PixeL free, PixeL seed, int r, int g, int b);
extern void DrawCumulativeEdgeDC(CClientDC *dc, PixeL free, PixeL seed, int r, int g, int b);
extern void DrawGridLinesDC(CClientDC *dc, int y, int size, int r, int g, int b);
extern void DrawBoundingBoxDC(CClientDC *dc, int size, Box fg_box, int r, int g, int b);

extern void ClearRegion(int xw, int yw, int size, PixeL free, Image& image);
extern int FindIntersectingRegion3(int xw, int yw, int size, PixeL free, 
					  Image& cumulative, Image& interbound);
extern int FindIntersectingRegion4(int xw, int yw, int size, PixeL free, PixeL st_point);
extern float FindMinCostInIntersectingBoundary(int xw, int yw, int old_size, 
			   int size, PixeL free, PixeL old_free);
extern void ClearNodes2(int width, int height, Nodes& item);
extern void MasktoNodes3(int width, int height, int size, PixeL seed, 
				Image& image, Image& intersect, float min_cost, int max_grad, 
				Image& gradient, Image& laplacian, Nodes& item);
extern void Intelligent_Scissor2(int size, int width, int height, Nodes& item);
extern void SaveCumulativeNextPxl(PixeL free, PixeL seed, PxlImage& next_pxl,
			PxlImage& cumulative_next_pxl);
extern void SaveCumulativeNextPxl2(PixeL free, PixeL seed, PxlImage next_pxl,
			PxlImage cumulative_next_pxl);

extern void renderImage();
extern void makeTexture();
extern void LoadFGImage();

extern void DrawCumulativeMotifEdge(PixeL seed, PixeL free, Image image);

extern void MakeSeedMap();

extern void renderInit(void);

extern void TIPDrawScene(CDC *pDC);
extern void TIPPanoramicDrawScene(CDC *pDC);
extern void TIPDrawSceneModelFromCylinder(GLenum mode);
extern void TIPDrawSceneModelFromVL(GLenum mode);
extern void TIPPanoramicRenderInit() ;
extern void TIPMakeCorrectSceneModelFromCylinder(double radius, int h_size);
extern void TIPMakeSceneModelFromVL(double radius, int h_size, int v_size);
extern void TIPCalFGCoord(void);

extern void TIPMakeCylinder(double radius, int h_size);
extern void TIPMakeSphere(double radius, int h_size, int v_size);
extern void TIPDrawSphere(GLenum mode);
extern void TIPDrawCylinder(GLenum mode);
extern void TIPPanoramaInit();
extern void TIPCopyMemDCfromDbuffer();

extern void DrawPointMemDC(CDC *dc, int x, int y, int size, int r, int g, int b);
extern void DrawEdgeMemDC(CDC *dc, PixeL free, PixeL seed, int r, int g, int b);
extern void DrawCumulativeEdgeMemDC(CDC *dc, PixeL free, PixeL seed, int r, int g, int b);
extern void DrawSeedsMemDC(CDC *dc, int sc, int r, int g, int b);
extern void DrawBoxMemDC(CDC *dc, int size, int x, int y, int r, int g, int b);
extern void MasktoNodes2AvoidingIntersection(int width, int height, int size, PixeL seed, Image& image, int max_grad, 
				 Image& gradient, Image& laplacian, PxlImage& cumulative_next_pxl, Nodes& item);
extern void MasktoNodes3AvoidingIntersection(int width, int height, int size, PixeL seed, Image& image, Image& intersect, 
		float min_cost, int max_grad, Image& gradient, Image& laplacian, PxlImage& cumulative_next_pxl, Nodes& item);
extern int CheckIfInfiniteLoopExists(PixeL free, PixeL seed);

//extern void ImageMoment(int size, int center_x, int center_y, Stroke *stroke_ptr);
extern void ImagetoNodesNPR(int width, int height, Image image, int max_grad, Image gradient);
extern void InitNodesNPR(int width, int height);
extern void getLaplacian(int width, int height, Image laplacian, Image image);

//extern BOOL SaveBMP(BITMAPINFO& bi, const char* img, const char* lpszFileName);



//////////////////////////////////////////
// Genetic Painter
extern uday_data data[500];
extern double f_diff(CDC& dc, uday_data data[]);
extern int getBlurredGradient(int width, int height, Image& gradient, Image& image, int N);

extern void SaveFileBMP2( int width, int height, const char* fn );
extern int GeneticPainter(CDC& dc);
extern void CopyDbuffertoMemDC(CDC& target);
extern void CopyMemDCtoDbuffer(CDC& dc);
extern void PerturbRGB(GLubyte& r, GLubyte& g, GLubyte& b);
extern void PerturbRGB(GLubyte& r, GLubyte& g, GLubyte& b, int MAX);
extern inline double drand48();
extern int intrand(int lo, int hi);

extern Field gfield, gfield2;

extern void InitGradField(int image_x, int image_y, Field& gfield, Image& image);
extern void InitGradField(int image_x, int image_y, Field& gfield, imatrix& gray);
extern void AlignGradient(int image_x, int image_y, Field& gfield, int size);
extern void AlignGradient2(int image_x, int image_y, Field& gfield, int size);
extern void GVF(int image_x, int image_y, Field& gfield, int N);
extern int getSobelGradient(int width, int height, Image gradient, Image image);

extern HGLRC		hrc, hrc2; 
extern HGLRC       memRC;
extern CRect       m_oldRect;
extern float       m_fRadius;
extern void OpenGL_SwitchToMemDC();

extern matrix P[10], Q[10], A[10], B[10]; // global variables
extern int pc; // intelligent pen counter
extern MRBspline curve, tmp_curve;
extern int DrawMRcurveMemDC(CDC *dc, MRBspline& curve, int level_down, deque<PixeL>& pnts, int r, int g, int b);
extern int WaveletTransform(CDC *dc, MRBspline& curve, int level_down, deque<PixeL>& pnts, int r, int g, int b);
extern void WaveletTransform2(CDC *dc, MRBspline& curve, deque<PixeL>& pnts, int r, int g, int b);
extern int DrawFilteredMRCurveMemDC(CDC *dc, MRBspline& curve, int level_down, int r, int g, int b);
extern void CopyEditedCurve(int j, MRBspline& tmp_curve, MRBspline& curve);

extern void DrawGradientField(CDC& dc, int width, int height, Field& gfield);
extern void DrawGradientField_LIC(CDC& dc, int width, int height, Field& gfield, int length);
extern void DrawGradientField2(CDC& dc, int width, int height, Field& gfield);
extern void GVF2(int image_x, int image_y, Field& gfield, int N, double mu);

extern GLubyte al_texture2[TEXTURE][TEXTURE];
extern GLubyte al_texture[TEXTURE][TEXTURE][4];
extern GLubyte al_texture3[4][TEXTURE][TEXTURE];
extern GLuint tex_name;
extern GLuint tex_names[4];
extern GLubyte *brush_buffer;
extern void InitAlphaMap(CDC& dc);
extern int	LoadBMP(CDC& dc, char* filename, CDC *memDC);

extern deque<PixeL> pnts;

extern void FindCurve(IPen* ip, int pc, PixeL free, int& min_ip, int& min_curve);

extern void ClearMemDC(CDC *dc);
///////////////////////////////////////////////////////
// Intelligent Pen
extern void DrawIPenCumulativeMRSampledEdgeMemDC(CDC *dc, IPen* ip, int pc, int mode, int r, int g, int b);
extern void ClearAllImagesDepthFirst(int x, int y);
extern void DrawBWireMRSampledEdgeMemDC(CDC *dc, MRBspline& curve, int level_down, PixeL free, PixeL seed, int r, int g, int b);
extern void Intelligent_Scissor_Oriented(int size, int width, int height, Nodes& item);
extern void Intelligent_Scissor2_Oriented(int size, int width, int height, Nodes& item);

extern int MakeGaussMask(double sigma, double* gau);
extern void DrawImage(CDC& dc, int image_x, int image_y, Image& image);
//extern double GaussSmooth(int image_x, int image_y, Image& image, int gau_w);
extern void GaussSmooth(int image_x, int image_y, imatrix& image, int gau_w);
extern int gau_w;
extern double gau[GAU_MAX], gau2[30];
extern inline double norm2(double vx, double vy);
extern void NonmaxSuppress(int image_x, int image_y);
extern void GetGrayImage(int width, int height, imatrix& image, GLubyte *Dbuffer);
extern void DrawGrayImage(CDC& dc, int image_x, int image_y, imatrix& image);
extern double GlobalCanny(int image_x, int image_y, imatrix& image, imatrix& image2, int gau_w);
extern void NonmaxSuppressGray(int image_x, int image_y);
extern imatrix gray, gray2, emap, gray3, dog;
extern mymatrix output1, output2;
extern double lambda;
extern cimatrix cmap;
extern double max_grad2;
extern double LocalCanny(int image_x, int image_y, int size, PixeL seed, imatrix& image, 
					   imatrix& image2, int gau_w);
extern void DrawGrayImageMask(CDC& dc, int image_x, int image_y, int size, PixeL seed, imatrix& image);
extern double cur_sigma;
//extern NameDialog* dlg;
extern double hi_thres, lo_thres;
extern double factor1, factor2;
extern int	MASK_SIZE, MAX_MASK_SIZE, MIN_MASK_SIZE;
extern double max_cur_sigma, max_hi_thres, max_lo_thres;
extern matrix tmp_x, tmp_y, G_x, G_y, G_mag, scale_map;
extern imatrix thin_edge;
extern imatrix off_mark, on_mark;
extern imatrix pixel_mark;
extern int pix_d_count;
extern pixel_dist pix_d[10];

extern void OnGlobalCanny(CDC& dc, int image_x, int image_y, imatrix& image, imatrix& image2);
extern double GetCannyGradient(int image_x, int image_y, imatrix& image, int gau_w);
extern void EdgeColor(CDC& dc, imatrix& gray2, int x0, int y0, GLubyte r, GLubyte g, GLubyte b);
extern edge_dist on_edge, off_edge, cur_edge;

extern void SaveBMPfromMemDC( CDC& dc, int width, int height, const char* fn );

extern RGB *pal;
extern void InvertImage(int image_x, int image_y, imatrix& image);
extern void DrawStrokesRandomOrientedMR(CDC& dc, int width, int height, Field& gfield, int max_w);
extern void DrawStrokesAdaptive(CDC& dc, int width, int height, Field& gfield);
extern void GaussSmooth2(int image_x, int image_y, imatrix& image, int gau_w);
extern void PerturbRGB(GLubyte& r, GLubyte& g, GLubyte& b);
extern inline void RGB2HSV( double r, double g, double b, double& h, double& s, double& v );
extern inline void HSV2RGB( double h, double s, double v, double& r, double& g, double& b );
extern void StoreEdgeStrokes(int image_x, int image_y, imatrix& image);
extern void DrawBsplineStrokes(CDC& dc, GLubyte r, GLubyte g, GLubyte b);
extern void DrawBsplineStrokes2(CDC& dc, GLubyte r, GLubyte g, GLubyte b);
extern void DrawEdgeStrokes2(CDC& dc, int image_x, int image_y, imatrix& image);
extern void StoreEdgeStrokes2(int image_x, int image_y, imatrix& image);
extern void ClearBuffer(GLubyte* buffer); 
extern void ClearGrayBuffer(GLubyte* buffer); 
extern void CopyMembuffer(int x, int y, int half, GLubyte* Membuffer, GLubyte* Doublebuffer);
extern void DrawEdgeStrokes(CDC& dc, int image_x, int image_y, imatrix& image, int width);
extern void DrawEdgeStrokes2(CDC& dc, int image_x, int image_y, imatrix& image, int width); 

// bilateral filtering
extern void GaussBlur(imatrix& image, double sigma);
extern void GaussBlurBilateral(int image_x, int image_y, imatrix& image, double sigma, double sigma2, int max_itr);
extern void GVF3(int image_x, int image_y, Field& gfield, int N, double factor);
extern void DistanceField(int image_x, int image_y, imatrix& image, double factor);

extern int SquareOverlap(CDC& dc, POINT* P);
extern int GLOBAL_CANNY_DONE; 

extern void FlowConnectEdges(int image_x, int image_y, imatrix& image, Field& gfield);
extern void ConvertThinEdgeMap2GrayImage(int image_x, int image_y, imatrix& thin_edge, imatrix& image);
extern void GetSobelGradient(int image_x, int image_y, imatrix& image, imatrix& image2, matrix& G_mag); 
extern void GetFlowGradient(CDC& dc, int image_x, int image_y, imatrix& image, Field& gfield, imatrix& image2, 
							int ker_l, double ker_w, double step_size); 
extern void NonmaximaSuppression(int image_x, int image_y, Field& gfield, matrix& G_mag, imatrix& image, double hi_thres);
extern void GetDOG(int image_x, int image_y, imatrix& image, matrix& G_mag, imatrix& image2, int index1, int index2,
			double tau);
extern void GetFlowDOG(int image_x, int image_y, matrix& G_mag, imatrix& dog, Field& gfield, 
				int ker_l, double step_size); 
extern MyList<MRBspline> bstrokes; // a list of B-spline strokes!
extern void StoreEdgeStrokes3(int image_x, int image_y, imatrix& image, MyList<MRBspline>& bstrokes);
extern void DrawBStrokes(CDC& dc, MyList<MRBspline>& bstrokes, GLubyte r, GLubyte g, GLubyte b);
extern void DrawHermiteStrokes(CDC& dc, MyList<MRBspline>& hstrokes, Field& gfield, 
						GLubyte r, GLubyte g, GLubyte b);
extern void DrawInterpBStrokes(CDC& dc, MyList<MRBspline>& bstrokes, GLubyte r, GLubyte g, GLubyte b);
extern void Thresholding(int image_x, int image_y, imatrix& image, double thres); 
extern void DrawGradientField_LIC3(CDC& dc, int width, int height, Field& gfield, int length);
extern void DrawInterpBStrokesRandCol(CDC& dc, MyList<MRBspline>& bstrokes);
extern void StoreEdgeStrokes4(int image_x, int image_y, imatrix& image, MyList<MRBspline>& bstrokes,
					   int thres); 
extern void EdgeThinning(CDC& dc, imatrix& image); 
extern void ConstructBStrokesMap(MyList<MRBspline>& bstrokes, qmatrix<MRBspline*>& b_map);
extern qmatrix<MRBspline*> b_map; 
extern void CleanUpBStrokes(MyList<MRBspline>& bstrokes, qmatrix<MRBspline*>& b_map);
extern void DrawGradientFieldArrow(CDC& dc, Field& gfield, GLubyte r, GLubyte g, GLubyte b);
extern void DrawTangentFieldArrow(CDC& dc, Field& gfield, GLubyte r, GLubyte g, GLubyte b);
extern void GetSobelGradient2(int image_x, int image_y, imatrix& image, imatrix& image2, matrix& G_mag); 
extern void GetDOG2(imatrix& image, matrix& G_mag, imatrix& image2, double sigma, double tau); 
extern void SampleBStrokes(MyList<MRBspline>& bstrokes, int interval);
extern void ConstructBStrokesMap2(MyList<MRBspline>& bstrokes, qmatrix<MRBspline*>& b_map);
extern void MakeGaussVector(double sigma, vector& GAU);
extern void SplitBstrokesAll(CDC& dc, MyList<MRBspline>& bstrokes, double sigma, double angle_thres);
extern void ConstructBStrokesMap3(MyList<MRBspline>& bstrokes, qmatrix<MRBspline*>& b_map);
extern void SmoothBstrokesAll(CDC& dc, MyList<MRBspline>& bstrokes, matrix& G_mag);
void MergeSingleBStroke(MRBspline* c1, MyList<MRBspline>& bstrokes, qmatrix<MRBspline*>& b_map, 
						double sigma, int half_w, double merge_thres);
extern void MergeBStrokesAll(MyList<MRBspline>& bstrokes, qmatrix<MRBspline*>& b_map, 
					  double sigma, double merge_thres, int half_w);
extern void MergeBStrokesAll2(MyList<MRBspline>& bstrokes, double sigma, double merge_thres, 
							  int half_w);
extern void SimplifyBstrokesAll(CDC& dc, MyList<MRBspline>& bstrokes, matrix& G_mag, double remove_thres);
extern void DrawBStrokesPoints(CDC& dc, MyList<MRBspline>& bstrokes, GLubyte r, GLubyte g, GLubyte b);
extern void DrawInterpBStrokes2(CDC& dc, MyList<MRBspline>& bstrokes, GLubyte r, GLubyte g, GLubyte b,
						 int pen_width);
extern void DrawInterpBStrokesRandCol2(CDC& dc, MyList<MRBspline>& bstrokes, int pen_width);
extern void LengthenBstrokesAll(CDC& dc, CDC& dc2, MyList<MRBspline>& bstrokes, Field& gfield);
extern void ConnectBstrokesAll(CDC& dc, CDC& dc2, MyList<MRBspline>& bstrokes, Field& gfield);
extern bool MergeSingleBStroke2(MRBspline* c1, int join1, MyList<MRBspline>& bstrokes, 
						 double sigma, double merge_thres, int half_w);
extern void GaussSmooth4(CDC& dc, int image_x, int image_y, double sigma);
extern void GaussBlurMemDC(CDC& dc, int image_x, int image_y, double sigma);
extern bool MergeSingleBStroke3(CDC& dc, MRBspline* c1, int join1, MyList<MRBspline>& bstrokes, 
						 double sigma, double merge_thres, int half_w);
extern void MergeBStrokesAll3(CDC& dc, MyList<MRBspline>& bstrokes, double merge_thres, int half_w);
extern void DrawNormalBStrokes(CDC& dc, MyList<MRBspline>& bstrokes, GLubyte r, GLubyte g, GLubyte b);
extern void InvertImatrix(imatrix& image);
extern void SnapBstrokesAll(CDC& dc, CDC& dc2, MyList<MRBspline>& bstrokes, Field& gfield);
extern void CopyMemDCtoImatrix(CDC& dc, imatrix& image);
extern void DrawGradientField_LIC4(CDC& dc, int width, int height, Field& gfield, int length);
extern void DrawGradientField_LIC5(CDC& dc, int width, int height, Field& gfield, int length);
extern void DrawGradientField_LIC6(CDC& dc, int width, int height, Field& gfield, int length);
extern void GetFlowBlur(CDC& dc, Field& gfield, imatrix& image, int length);
extern void GetFlowDiffusion(CDC& dc, Field& gfield, imatrix& image, matrix& G_mag, int half_l, int max_itr);
extern void GetFlowBilateral(CDC& dc, Field& gfield, imatrix& image, double sigma, 
							 double sigma3, int max_itr);
extern inline double gauss2(double x, double mean, double sigma);
extern inline double gauss1D_norm(double x, double mean, double sigma);
extern void GetFlowMedian(imatrix& image, Field& gfield, imatrix& image2, matrix& G_mag, 
				   int width, int length, int max_itr);
extern void GetFlowMedianDiffusion(imatrix& image, Field& gfield, imatrix& image2, 
								   matrix& G_mag, int max_itr);
extern void GetFlowBilateral2(CDC& dc, Field& gfield, imatrix& image, double sigma, 
					   double sigma3, int max_itr);
extern void GetFlowMedian2(CDC& dc, Field& gfield, imatrix& image, int width, int length, int max_itr);
extern void GetFlowBilateral3(Field& gfield, imatrix& image, imatrix& dog, double sigma, double sigma3, int max_itr);
extern void GetAnisotropic(imatrix& image, int max_itr);
extern void GetFlowLineBilateral(imatrix& image, Field& gfield, double sigma, double sigma3, int max_itr);
extern void GetFlowLineContrast(imatrix& image, Field& gfield, double sigma, int max_itr);
extern void GetFlowLineContrast2(imatrix& image, Field& gfield, imatrix& dog, double sigma, int max_itr);
extern void GetFlowLineBilateral2(imatrix& image, Field& gfield, double sigma, double sigma3, int max_itr);
extern void GetFlowLineBilateral3(imatrix& image, Field& gfield, double sigma, double sigma3, int max_itr);
extern void GetColorImage(int width, int height, cimatrix& cmap, GLubyte *Dbuffer);
extern void GetFlowLineColBilateral3(cimatrix& image, Field& gfield, double sigma, double sigma3, int max_itr);
extern void CopyCol2GrayImage(int image_x, int image_y, cimatrix& cmap, imatrix& image);
extern void DrawColorImage(CDC& dc, int image_x, int image_y, cimatrix& image); 
extern void BilateralColor(int image_x, int image_y, cimatrix& image, double sigma, double sigma2, int max_itr);
extern void DrawCartoonImage(CDC& dc, int image_x, int image_y, cimatrix& image, imatrix& gray); 
extern void ConstructCartoonImage(cimatrix& image, imatrix& gray, cimatrix& merged); 
extern void ConstructCartoonImage2(cimatrix& image, imatrix& gray, cimatrix& merged, double threshold); 
extern void ConstructCartoonImage3(cimatrix& image, imatrix& gray, cimatrix& merged);
extern void GetFlowLineMedian(imatrix& image, Field& gfield, double sigma, int max_itr);
extern void GetMedian(imatrix& image, double sigma, int max_itr);
extern void GaussSmooth5(imatrix& image, double sigma);
extern void GaussSmooth3(int image_x, int image_y, imatrix& image, int gau_w);
extern void ConstructMergedImage(imatrix& image, imatrix& gray, imatrix& merged); 
extern void GaussSmooth6(imatrix& image, double sigma);
extern void GetContrast(imatrix& image, matrix& G_mag, imatrix& image2, double sigma, double tao); 
extern void GetDOG3(imatrix& image, matrix& G_mag, imatrix& image2, double sigma, double tau); 
extern void GetDOG4(imatrix& image, matrix& G_mag, imatrix& image2, double sigma, double tau, double phi); 
extern void ColorQuantization(cimatrix& cmap, int level); 
extern void GetShock(imatrix& image, matrix& G_mag, imatrix& image2, double sigma, double tau); 
extern void GetShock2(imatrix& image, matrix& G_mag, imatrix& image2, double sigma, double tau, double z);
extern void GetFlowShock(imatrix& image, Field& gfield, imatrix& image2, matrix& G_mag, 
			  double sigma, double sigma3, double tau, double z); 
extern void GetColShock(cimatrix& cmap, imatrix& image, double sigma, double tau, double z); 
extern void GetFlowColShock(cimatrix& cmap, Field& gfield, imatrix& image, 
			  double sigma, double sigma3, double tau, double z); 
extern void GetColContrastEnhance(cimatrix& cmap, double pivot, double adjust); 
extern void GetColShock2(cimatrix& cmap, imatrix& image, double sigma, double tau, double z); 
extern void GetColShock3(cimatrix& cmap, imatrix& image, double sigma, double tau, int half, double z); 
extern void ConvertRGBtoLAB(cimatrix& cmap, cimatrix& target);
extern void ConvertLABtoRGB(cimatrix& cmap, cimatrix& target);
extern void GetFlowColBilateral3(cimatrix& image, Field& gfield, double sigma, double sigma2, double sigma3, int max_itr);
extern void GetFlowColBilateral4(cimatrix& image, Field& gfield, double sigma, double sigma2, double sigma3, 
						  int max_itr1, int max_itr2);
extern void GetFlowColBilateral5(cimatrix& image, Field& gfield, double sigma, double sigma2, double sigma3, 
						  int max_itr);
extern void CopyCmap2Membuffer(cimatrix& cmap, GLubyte* Dbuffer);
extern void GetFlowColBilateral6(cimatrix& image, Field& gfield, double sigma, double sigma2, double sigma3, 
						  int max_itr1, int max_itr2);
extern void GetFlowColBilateral7(cimatrix& image, Field& gfield, double sigma, double sigma2, double sigma3, 
						  int max_itr);
extern void GetPSNR(CDC& dc, cimatrix& cmap, char* file1, char* file2);
extern void GetFlowColBilateral9(cimatrix& image, Field& gfield, double sigma, double sigma2, 
								 double sigma3, double sigma4, int max_itr);
extern void GetPSNR_partial(CDC& dc, cimatrix& cmap, char* file1, char* file2);
extern void GetFlowColBilateralFlowAxis(cimatrix& image, Field& gfield, double sigma, double sigma3, int max_itr);
extern void GetFlowColBilateralGradAxis(cimatrix& image, Field& gfield, double sigma2, double sigma4, int max_itr);
extern void BilateralColorSep(int image_x, int image_y, cimatrix& image, double sigma, double sigma2, int max_itr);
extern inline void MakeGaussVectorBL(double sigma, vector& GAU, int max_index);
extern inline void RGB2LAB2( double R, double G, double B, double& L, double& a, double& b );
extern inline void LAB2RGB( double L, double a, double b, double& R, double& G, double& B );
extern inline void RGB2LAB( double R, double G, double B, double& L, double& a, double& b );
extern inline void RGB2LAB256( double R, double G, double B, double& L, double& a, double& b );
extern inline void LAB2RGB256( double L, double a, double b, double& R, double& G, double& B );
extern void ConvertRGBtoLAB256(cimatrix& cmap, cimatrix& target);
extern void ConvertLABtoRGB256(cimatrix& cmap, cimatrix& target);
extern void ConvertRGBtoLAB_double(cimatrix& cmap, cmatrix& target);
extern void ConvertLABtoRGB_double(cmatrix& cmap, cimatrix& target);
extern void GetFlowColBilateralLAB(cimatrix& cmap, Field& gfield, double sigma, double sigma2, 
						  double sigma3, double sigma4, int max_itr);
extern void ColorQuantizationLAB(cimatrix& cmap, int level); 
extern void DrawETF_LIC(CDC& dc, ETF& e, int length);
//extern void GetFlowDOG5(imatrix& image, ETF& e, double sigma, double sigma3, double tau); 
extern void GetFDoG(imatrix& image, ETF& e, double sigma, double sigma3, double tau); 
extern void GaussSmoothSep(imatrix& image, double sigma);
extern void DrawInverseGrayImage(CDC& dc, int image_x, int image_y, imatrix& image); 
extern void GetLDoG(imatrix& image, ETF& e, double sigma, double sigma3, double tau); 
extern void DrawCondInverseGrayImage(CDC& dc, imatrix& image, imatrix& bw);
extern void CreateRedZone(imatrix& redzone, double sigma); 
extern void GetFDoGRedZone(imatrix& image, ETF& e, imatrix& redzone, imatrix& bw, double sigma, double sigma3, double tau); 
extern void DrawCondInverseGrayImageRedZone(CDC& dc, imatrix& image, imatrix& bw, imatrix& redzone); 
extern void DrawETFColor360(CDC& dc, ETF& e);
extern void MakeGaussianVector(double sigma, myvec& GAU);
extern void DrawETFColorTensor(CDC& dc, ETF& e);
extern void DrawETFColor180(CDC& dc, ETF& e);
extern void BilateralToon(cimatrix& cmap, double sigma, double sigma2, int max_itr);
extern void BilateralToonSep(CDC& dc, cimatrix& cmap, double sigma, double sigma2, int max_itr);
extern void GaussColSmoothSep(cimatrix& image, double sigma);
extern void GetFMedian(imatrix& image, ETF& e, int itr);
extern void BilateralGrayToonSep(imatrix& gray, double sigma, double sigma2, int max_itr);
extern void GetFMedian2(imatrix& image, ETF& e, int itr); 
extern void GetFMedian3(imatrix& image, ETF& e, int itr);
extern void GetFMedian4(imatrix& image, int itr); 
extern void GetFMedian5(CDC& dc, imatrix& image, ETF& e, int itr); 
extern void CopyGray2Membuffer(imatrix& gray, GLubyte* Dbuffer);
extern void GaussSmoothSepDouble(imatrix& image, double sigma, matrix& output);
extern void GetColMaxMin(CDC& dc, cimatrix& cmap, double sigma, int itr); 
extern void GetColMaxMin2(CDC& dc, cimatrix& cmap, double sigma, double tau, int itr); 
extern void GetDirectionalDoG(imatrix& image, ETF& e, mymatrix& dog, myvec& GAU1, myvec& GAU2, double tau);
extern void GetFlowColBilateral10(cimatrix& image, ETF& e, double sigma, double sigma2, 
						  double sigma3, double sigma4, int max_itr);
extern void GetFlowInfluence(ETF& e, mymatrix& dog, mymatrix& tmp, myvec& GAU3);
extern void MeanCurvatureFlow(imatrix& image, int max_itr);
extern void MeanCurvatureFlow2(imatrix& image, imatrix& line, int max_itr);
extern void GetDogSep(imatrix& image, double sigma, imatrix& lap, double tau);
extern void GrayQuantization(imatrix& image, int level);
extern void GetLaplacian6(imatrix& image, double sigma, imatrix& lap, double tau);
extern void GetAnisotropic2(imatrix& image, imatrix& line, int max_itr);
extern void BilateralColorTomasi(int image_x, int image_y, cimatrix& image, double sigma, double sigma2, int max_itr);
extern void BilateralColorChui(cimatrix& image, double sigma, double sigma2, int max_itr);
extern void MinMaxFlow(imatrix& image, double sigma, int max_itr);
extern void MinMaxFlow2(imatrix& image, int radius, double thres, int max_itr);
extern void MinMaxFlow3(imatrix& image, int radius, int max_itr);
extern void MeanCurvatureFlow3(imatrix& image, ETF& e, int max_itr);
extern void MeanCurvatureFlowSapiro(imatrix& image, int max_itr);
extern void MinMaxFlow4(imatrix& image, int radius, int max_itr);
extern void Dilation(imatrix& image, int itr); 
extern void Closing(imatrix& image, int size, int itr); 
extern void ConstructMergedImageColor(cimatrix& image, imatrix& gray, cimatrix& merged); 
extern void MeanCurvatureFlowColor2(cimatrix& cmap, imatrix& line, int max_itr);
extern void MeanCurvatureFlowColor(cimatrix& cmap, int max_itr);
extern void MeanCurvatureFlowColorLAB2(cimatrix& cmap, imatrix& line, int max_itr);
extern void GrayThresholding(imatrix& image, double thres); 
extern void GetFlowDoG(ETF& e, mymatrix& dog, mymatrix& tmp, myvec& GAU3);
extern void GetFlowDoGLIC(ETF& e, mymatrix& dog, mymatrix& tmp1, mymatrix& tmp2, myvec& GAU3, matrix& inner, matrix& outer, double tau);
extern void GetFlowDoGLIC2(ETF& e, mymatrix& dog, mymatrix& tmp, myvec& GAU3, matrix& inner, matrix& outer, double tau, imatrix& lambda);
extern void GetDogSepDouble(imatrix& image, double sigma, mymatrix& dog, double tau);
extern void GetDirectionalDogSepDouble(imatrix& image, double sigma, mymatrix& dog, double tau, ETF& e, myvec& GAU1, myvec& GAU2);
extern void GetDirectionalDogSepDouble2(imatrix& image, double sigma, mymatrix& dog, double tau, ETF& e, myvec& GAU1, myvec& GAU2, imatrix& edge);
extern void GetEdgeStructure(ETF& e, imatrix& image, imatrix& edge, double sigma, double tau);
extern void GetDoGLIC(imatrix& image, ETF& e, double sigma, double sigma3, double tau);
extern void GetDoGLIC2(imatrix& image, ETF& e, imatrix& edge, double sigma, double sigma3, double tau);
extern void GetDoGLICAdjust(imatrix& image, double lambda);
extern int LoadBMP3(CDC& dc, char* filename, CDC *memDC);
extern void GetGrayImageFromMemDC(int image_x, int image_y, imatrix& image, CDC& dc);
extern void ConstructMergedImageMult(imatrix& image, imatrix& gray, imatrix& merged);
extern void GetFDoGDouble(imatrix& image, ETF& e, matrix& gmag, double sigma, double sigma3, double tau);
extern void NonmaximaSuppressionETF(matrix& gmag, ETF& e, imatrix& image, double thres);
extern void DrawInterpBStrokes3(CDC& dc, MyList<MRBspline>& bstrokes, GLubyte r, GLubyte g, GLubyte b,
						 int pen_width);
extern void ReplaceStrokeTexture(matrix& gmag, ETF& e, imatrix& image, double thres);
extern void ReplaceStrokeTexture2(matrix& fdog, ETF& e, imatrix& image, double thres);
extern void GaussSmoothSepDouble2(matrix& input, imatrix& gray, double sigma, matrix& output);
extern void ReplaceStrokeTexture3(matrix& fdog, ETF& e, imatrix& image);
extern void ReplaceStrokeTexture4(matrix& fdog, ETF& e, imatrix& image, imatrix& tex);
extern int	LoadBMP4(CDC& dc, char* filename, CDC *memDC, int &image_x, int &image_y);
extern void GetGrayContrastEnhanceDouble(matrix& map, double pivot, double adjust);
extern void ReplaceStrokeTexture5(matrix& fdog, ETF& e, imatrix& image, imatrix& tex);
extern void ReplaceStrokeTexture6(matrix& fdog, ETF& e, imatrix& image, imatrix& tex);
extern void ReplaceStrokeTexture7(matrix& fdog, ETF& e, imatrix& image, imatrix& tex);
extern void ReplaceStrokeTexture9(matrix& fdog, ETF& e, imatrix& image, imatrix& tex, bool sumie);
extern void GetFDoGDoubleSubpixel(imatrix& image, ETF& e, matrix& gmag, double sigma, double sigma3, double tau); 
extern void ReplaceStrokeTexture10(matrix& fdog, ETF& e, imatrix& image, imatrix& tex, bool sumie);
extern void DownsampleGrayImage(imatrix& gray, imatrix& gray_s);
extern void DownsampleColorImage(cimatrix& cmap, cimatrix& cmap_s);
extern void CopyCol2GrayImage2(cimatrix& cmap, imatrix& image); 
extern void UpsampleGrayImage(imatrix& gray_s, imatrix& gray); 
extern void UpsampleColorImage(cimatrix& cmap_s, cimatrix& cmap);
extern void ColorQuantization2(cimatrix& cmap, int level);
extern void ColorQuantizationRGB(cimatrix& cmap, int level);
extern void ColorQuantizationLAB2(cimatrix& cmap, int level1, int level2, int level3);
extern void DilationColor(CDC& dc, cimatrix& cmap); 
extern void DilationColor2(CDC& dc, cimatrix& cmap, imatrix& line, double thres, double factor); 
extern void ClosingColor(CDC& dc, cimatrix& cmap, imatrix& line, int itr);
extern void OpeningClosingColor(CDC& dc, cimatrix& cmap, imatrix& line, int itr);
extern iRGB *CreatePaletteFromImage(cimatrix& cmap, int& N);
extern void ColorQuantizationUsingPalette(cimatrix& cmap, iRGB* pal, int N);
extern void GetColorImageFromMemDC(int image_x, int image_y, cimatrix& cmap, CDC& dc);
extern void OpeningClosingBlack(CDC& dc, imatrix& gray, int itr);
extern void UpsampleColorImageNearest(cimatrix& cmap_s, cimatrix& cmap);
extern void GetFlowDilation(imatrix& line, ETF& e, int half_l, int max_itr);
extern void ReplaceStrokeTexture11(matrix& fdog, ETF& e, imatrix& image, imatrix& thin, imatrix& tex, bool sumie);
extern void GaussSmoothSepDouble3(imatrix& image, double sigma, matrix& output);
extern void ErosionGray(CDC& dc, imatrix& gray, int itr);
extern void DilationGray(CDC& dc, imatrix& gray, int itr);
extern void RemoveIsolatedRegions(CDC& dc, imatrix& gray, int size_max); 
extern void UpsampleGrayImageNearest(imatrix& gray_s, imatrix& gray);
extern void GetLaplacianDouble(imatrix& image, double sigma, matrix& G_map, double tau, double thres);
extern void ErosionGrayCircle(CDC& dc, imatrix& gray, int itr);
extern void ReplaceStrokeTexture12(matrix& fdog, ETF& e, imatrix& image, imatrix& thin, imatrix& tex, bool sumie);
extern void ReplaceTextureWithDilation(CDC& dc, imatrix& gray, ETF& e, imatrix& tex, int itr);
extern void DilationGraySep(CDC& dc, imatrix& gray, int itr); 
extern void GetDogSepDouble2(imatrix& image, double sigma, matrix& dog, double tau);
extern void GetOffsetLines(CDC& dc, imatrix& gray);
extern void ReplaceStrokeTexture8(matrix& fdog, ETF& e, imatrix& image, imatrix& thin, imatrix& tex, bool sumie);
extern void ReplaceStrokeTexture13(matrix& fdog, ETF& e, imatrix& image, imatrix& thin, imatrix& tex, bool sumie);
extern void RemoveSmallRegions(CDC& dc, imatrix& line, int size);
extern void MeanCurvatureFlowColor3(cimatrix& cmap, ETF& e, int max_itr);
extern void GetColShockFast(cimatrix& cmap, imatrix& image, double sigma, double tau, double z); 
extern void GetColShockFastHSV(cimatrix& cmap, imatrix& image, double sigma, double tau);
extern void GetColShockFastMinMax(cimatrix& cmap, imatrix& image, double sigma, double tau);
extern void GetFlowShockDoGETF(imatrix& image, ETF& e, matrix& dog, vector& GAU1, vector& GAU2, double tau);
extern void GetFlowShockDoGETF2(imatrix& image, ETF& e, imatrix& sign, double sigma, double tau);
extern void GetFlowShockDoGETF3(imatrix& image, ETF& e, imatrix& sign);
extern void GetWeickertUpwind(cimatrix& cmap, imatrix& sign, double h, int itr);
extern void BilateralGraySep(imatrix& gray, double sigma, double sigma2, int max_itr);
extern void DrawETFArrow(CDC& dc, ETF& e, GLubyte r, GLubyte g, GLubyte b);
extern void GetDiff(CDC& dc, imatrix& gray, char* file1, char* file2);
extern void GetFBL(cimatrix& image, ETF& e, double sigma, double sigma2, double sigma3, double sigma4, int max_itr);
extern void BilateralColorFull(cimatrix& image, double sigma, double sigma2, int max_itr);
extern void GetSABL(CDC& dc, cimatrix& image, ETF& e, double sigma, double sigma2, double sigma3, double sigma4, int max_itr);
extern void GetAnisotropicBL(cimatrix& image, ETF& e, double sigma, double sigma2, double sigma3, double sigma4, int max_itr);
extern void GetFlowShockInfluenceETF(ETF& e, matrix& dog, imatrix& sign, vector& GAU3);
extern void Get1DFlowColShockAverageHSVETF(cimatrix& cmap, ETF& e, cimatrix& tmp, imatrix& sign, vector& GAU3, 
								 double z1, double z2);
extern void Get1DPerpFlowShockAverageHSVETF(cimatrix& cmap, ETF& e, cimatrix& tmp, imatrix& sign, 
								  vector& GAU2, double z1, double z2);
extern void GetFlowColShockETF(cimatrix& cmap, ETF& e, imatrix& gray, 
			  double sigma, double sigma3, double tau, double z1, double z2); 

extern inline double gauss1D_norm2(double x, double mean, double sigma);
extern inline double gauss1D_bl(double x, double mean, double sigma);
extern void StartTimer(void);
extern double ElapsedTime(void);
inline extern void EndTimer(char *file);

#endif
