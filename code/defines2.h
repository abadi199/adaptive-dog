#define BLACK_		0
#define YELLOW_	1
#define CYAN_		2
#define RED_		3
#define BLUE_		4
#define GREEN_		5
#define MAGENTA_	6
#define WHITE_		7

#define PLANE		0
#define CYLINDER	1
#define SPHERE		2

/// status list

#define TWOD	0
#define RENDER	2
#define PANORAMA_RENDER_MODE	8
#define FIS		3
#define INTEL_IS	5
#define LW_ 		7
#define INTEL_ 10
#define INTEL2_ 11
#define LL_ 	  12
#define MR_FIS_	15 
#define FIS2_ 17
#define FIS_BG 18
#define FIS_BG2 19
#define FORWARD_ 20
#define BACKWARD_ 30
#define ACTIVATED_MR_ 33
#define ACTIVATED_LL_ 35
#define ACTIVATED_LW 37
#define ACTIVATED_IS 39
#define ACTIVATED_ 40
#define ACTIVATED2_ 45
#define ACTIVATED3_ 47
#define MRC_   50
#define MRC2_   55
#define MRC_AUTO_ 	57
#define PICK_REGION_ 60
#define PICK_MOTIF_ 65
#define PICK_QUANT_ 67
#define BLEED_POINT_	70
#define TEXTURE_MANUAL_	71
#define ABSTRACT_POINT_	75
#define MOSAIC_POINT_	77
#define LINE_START_	80
#define LINE_MOTION_ 90
#define ERODE_REGION_ 100
#define DILATE_REGION_ 110
#define STROKE_AUTO_ 120
#define STROKE_MANUAL_ 130
#define CROSSHATCH_MANUAL_ 132
#define STROKE_INREGION_ 135
#define BOUNDARY_PIXMAP_ 140
#define NO_STATUS	170
#define DIFF	200
#define FIS_MOVE	210
#define FIS_CREATE	220
#define FIS_DELETE	230
#define SEED_MOVE	240
#define SEED_MOVE_ACTIVATED 250
#define SEED_DELETE 260
#define STROKES 270
#define MY_STROKES 280
#define AUTO_STROKES 290
#define FREEHAND_CURVE 300
#define FREEHAND_ACTIVATED 310
#define CURVE_EDIT 320
#define CURVE_EDIT_ACTIVATED 330
#define IPEN_EXTERNAL	340
#define ACTIVATED_IPEN_EXTERNAL	350
#define IPEN_INTERNAL	360
#define ACTIVATED_IPEN_INTERNAL	370
#define STROKE_EDIT	380
#define STROKE_EDIT_ACTIVATED	390
#define STROKE_EDIT_ACTIVATED2	400
#define STROKE_EDIT_INDIV	410
#define STROKE_EDIT_INDIV_ACTIVATED	420
#define STROKE_EDIT_INDIV_ACTIVATED2	430
#define EDGE_PAINT 440
#define EDGE_PAINT_ACTIVATED 450
#define EDGE_TUNE 460
#define EDGE_SAMPLE 470
#define EDGE_PICK 480
#define EDGE_SAMPLE_ACTIVATED 490
#define PIXEL_SAMPLE 500
#define PIXEL_SAMPLE_ACTIVATED 510

////////////////////////
// status2
#define SCALE_UP 1000
#define SCALE_DOWN 1010
#define HI_THRES_UP 1020
#define HI_THRES_DOWN 1030
#define LO_THRES_UP 1040
#define LO_THRES_DOWN 1050
#define SIZE_UP 1060
#define SIZE_DOWN 1070

///////////////////////////
// stroke_mode
#define UNIFORM_THICK 10
#define MIDDLE_THICK 20
#define BEGIN_THICK 30
#define END_THICK 40


///////////////////////// 
// IMAGE RENDER MODE
#define OPENGL_MODE 10
#define BITMAP_MODE 20

//#define MAXWIDTH 2048
//#define MAXHEIGHT 1024
//#define MAXWIDTH 1000
//#define MAXHEIGHT 1000
#define MAXWIDTH 700
#define MAXHEIGHT 700

#define MAX_SCROLL_VIEW_X 1250
#define MAX_SCROLL_VIEW_Y 750

#define MAX_CURVE_PNTS     300
#define MAX_REGIONS        200
#define MAX_CHILD_REGIONS  5  
#define MAX_CD_IN_REGION   20

#define MAX_COST	10000

#define MAX_NOISE_SIZE		40

#define MAX_QUANT_COLOR		10

/////////////////////////////////////////////////////////////////////
// Image sizes for TIP

#define DELAY 1000000000

//#define TEXTURE	512
//#define TEXTURE_WIDTH	512
//#define TEXTURE_HEIGHT	256
//#define TEXTURE	32
#define TEXTURE	64
#define TEXTURE_WIDTH	8
#define TEXTURE_HEIGHT	8

//#define MAX_TIP_FG_NUM 18
#define MAX_TIP_FG_NUM 1

//#define BACK_TEXTURE_WIDTH 256
#define BACK_TEXTURE_WIDTH 8

//#define MAX_FRAMES 51
//#define MAX_FIS_FRAMES 51
//#define MAX_TEX_NAMES 200

#define MAX_FRAMES 1
#define MAX_FIS_FRAMES 1
#define MAX_TEX_NAMES 1

#define FIS_POINT_SIZE	5

#define GAU_MAX 100

/////////////////////////////////////////////////////////

#define MAX_SEEDS 50

#define MAX_IS_MASK_SIZE			10
//#define MAX_MASK_SIZE				30

#define MAX_RAND		pow(2., 15.)

#define MAX_USR_COL_NUM    8

#define MAX_RGB_COLORS_    10

#define COLORSHIFT 8
#define COLORSHIFT2 16
#define COLORSHIFT3 24

#define TkGetColorFromRGB(r, g, b) \
            ((r) | ((g) << COLORSHIFT) | ((b) << COLORSHIFT2))

#define TkGetColorFromRGBA(r, g, b, a) \
    ((r) | ((g) << COLORSHIFT) | ((b) << COLORSHIFT2) | ((a) << COLORSHIFT3))

#define GLRED(c) ((c) & 0xFF)
#define GLGREEN(c) (((c) >> COLORSHIFT) & 0xFF)
#define GLBLUE(c) (((c) >> COLORSHIFT2) & 0xFF)

#define TRUE_		1
#define FALSE_ 	0

#define ABS(x) ( ((x)>0) ? (x) : (-(x)) )

//#define X	0
//#define Y	1
//#define Z	2

#define Sign(x)   ((x) >= 0) ? (((x) > 0) ? 1 : 0) : (-1)

#define round(x) ((int) ((x) + 0.5))

#define HEAD 0 // head of stroke
#define TAIL 1 // tail of stroke


