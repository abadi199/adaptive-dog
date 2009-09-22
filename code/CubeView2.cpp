// TIPView2.cpp : implementation of the CCubeView class
//
#include "stdafx.h"

#include "Cube.h"
#include "CubeDoc.h"
#include "CubeView.h"
#include "NameDialog.h"

#include "defines.h"
#include "globals.h"

extern NameDialog* dlg;

#include <cmath>
#include <deque>
using namespace std;

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#include <limits.h>

#define PI 3.1415926535
#define INFI 	10000

#define STEP	10

extern int cal3D;
extern point3D camera, Prevcamera, vanVec;
extern float camRot;
extern float camRot2;
extern float viewAngle;
extern char name[256];
extern char prefix[256];

extern GLuint texName[MAX_TIP_FG_NUM];

extern int isVanish, fg_num;

static GLfloat winWidth, winHeight;
extern void swap(void);
extern int file_loaded;
extern GLubyte *Dbuffer;
extern GLubyte *Maskbuffer;
extern GLubyte *Backbuffer;

extern float x_size, y_size;
extern point2D ImageCoord[7];

extern hCoord WC[7];
extern hCoord FG[MAX_TIP_FG_NUM][MAX_FRAMES][MAX_TIP_FG_NUM];

extern point3D Wc, ray;
extern point2D Vc;
extern point3D v[MAX_TIP_FG_NUM];
extern float half_x, half_y;

extern void CalImageCoord();
extern void CalWorldCoord();
extern void renderNone(void);

extern point2D vanish;
extern point2D movement;
extern int isVanish, isPnt;
extern ObjPnt which;
extern int pcount;

extern int status;
extern int fg_num;

extern int IMAGE_X, IMAGE_Y;
extern int view_x, view_y;

int play_mode = OFF_PLAY;
int video_mode = OFF_VIDEO;
int play_direction = FORWARD;

//#define LEVEL_DOWN 5  // for multiresolution curve
int LEVEL_DOWN = 0;  // for multiresolution curve

//MRBspline curve[20]; // IB-snakes

IPen ip[500]; // intelligent pen strokes
int pc; // pen counter

deque<PixeL> pnts;
MRBspline curve, tmp_curve;

//int stroke_mode = UNIFORM_THICK;
int stroke_mode = MIDDLE_THICK;
int min_ip, min_curve;

/////////////////////////////////////////////////////////////////////////////
// CCubeView construction/destruction

void CCubeView::OnLButtonDown(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	button_pressed(nFlags, 1, point.x, point.y);
	
	CView::OnLButtonDown(nFlags, point);
}

void CCubeView::OnRButtonDown(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	button_pressed(nFlags, 2, point.x, point.y);

	CView::OnRButtonDown(nFlags, point);
}

void CCubeView::OnLButtonUp(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	button_released(1, point.x, point.y);

	CView::OnLButtonUp(nFlags, point);
}

void CCubeView::OnRButtonUp(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	button_released(2, point.x, point.y);

	CView::OnRButtonUp(nFlags, point);
}

void CCubeView::OnMouseMove(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	mouse_motion(nFlags, point.x, point.y);	

	CView::OnMouseMove(nFlags, point);
}

extern int key_pressed;

void CCubeView::OnKeyUp(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	// TODO: Add your message handler code here and/or call default
	//if (nChar != VK_PROCESSKEY) {
		key_input(nChar);
	//}
	
	CView::OnKeyUp(nChar, nRepCnt, nFlags);
}

int cur_level = 0;

void CCubeView::key_input(UINT nChar)
{
	CClientDC dc(this);

	switch (nChar) {
		case 68: // 'd' (demonstration of continuous curve blending)
			break;
		case 88: // Stop the program and save the current memDC
			key_pressed = 1;
			break;
		case 83: // 's'
			TRACE("save BMP\n");
			SaveBMPfromMemDC(dc, IMAGE_X, IMAGE_Y, "_result.bmp");
			break;
		case 84: // 't'
			break;
		case 67: // 'c'
			break;
		case VK_DOWN:
			break;
		case VK_UP:
			break;
			
		case 'w':
		case '8':
			camera.z-=STEP;
    		break;
		case VK_NUMPAD8:
			camRot2 += 1.0;
			break;
		case VK_NUMPAD2:
			camRot2 -= 1.0;
			break;
		case '2':
			camera.z+=STEP;
    		break;
		case VK_LEFT:
			camRot += 1.0; // M_PI/180
			break;
		case '4':
			camera.x-=STEP;
    		break;
		case VK_RIGHT:
			camRot -= 1.0; // M_PI/180
			break;
		case '6':
			camera.x+=STEP;
    		break;
		case '1' : 
			camera.x-=STEP/2;
			camera.z+=STEP;
			break;
		case '3' :
			camera.x+=STEP/2;
			camera.z+=STEP;
			break;
		case '7' :
			camera.x-=STEP/2;
			camera.z-=STEP;
			break;
		case '9' :
			camera.x+=STEP/2;
			camera.z-=STEP;
			break;
		case VK_NUMPAD6:
			viewAngle+=2;
			break;
		case VK_NUMPAD4:
			viewAngle-=2;
			break;
	}

	if (status == RENDER)
		Invalidate(FALSE);
}

void CCubeView::button_pressed(UINT flag, int button, int x, int y)
{
	switch (button) {
		case 1: LeftButtonPressed(x, y);
				break;
		case 2: RightButtonPressed(flag, x, y);
				break;
	}
}

//Stroke *stroke_ptr = new Stroke;
int idx = -1; // curve control point index

void CCubeView::LeftButtonPressed(int x, int y)
{
	//key_pressed = 1;

	CClientDC dc(this);

}

void CCubeView::RightButtonPressed(UINT flag, int x, int y)
{
	CClientDC dc(this);

}

void CCubeView::button_released(int button, int x, int y)
{
	CClientDC dc(this);
}

void CCubeView::mouse_motion(UINT flag, int x, int y)
{
	static int old_snapped;
			
	CClientDC dc(this);


}


