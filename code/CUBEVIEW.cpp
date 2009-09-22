#include "stdafx.h"
#include "cube.h"
#include "cubedoc.h"
//#include "cubeview.h"
#include "NameDialog.h"

#include ".\cubeview.h"

#include "gl\gl.h"
#include "gl\glu.h"

#include <string>
#include <cmath>
#include <iostream>
#include <string>
using namespace std;

#include "defines.h"
#include "globals.h"

NameDialog* dlg;
CRect dlg_rect;


#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

unsigned char threeto8[8] =
{
	0, 0111>>1, 0222>>1, 0333>>1, 0444>>1, 0555>>1, 0666>>1, 0377
};

unsigned char twoto8[4] =
{
	0, 0x55, 0xaa, 0xff
};

unsigned char oneto8[2] =
{
	0, 255
};

static int defaultOverride[13] =
{
	0, 3, 24, 27, 64, 67, 88, 173, 181, 236, 247, 164, 91
};

static PALETTEENTRY defaultPalEntry[20] =
{
	{ 0,   0,   0,    0 },
	{ 0x80,0,   0,    0 },
	{ 0,   0x80,0,    0 },
	{ 0x80,0x80,0,    0 },
	{ 0,   0,   0x80, 0 },
	{ 0x80,0,   0x80, 0 },
	{ 0,   0x80,0x80, 0 },
	{ 0xC0,0xC0,0xC0, 0 },

	{ 192, 220, 192,  0 },
	{ 166, 202, 240,  0 },
	{ 255, 251, 240,  0 },
	{ 160, 160, 164,  0 },

	{ 0x80,0x80,0x80, 0 },
	{ 0xFF,0,   0,    0 },
	{ 0,   0xFF,0,    0 },
	{ 0xFF,0xFF,0,    0 },
	{ 0,   0,   0xFF, 0 },
	{ 0xFF,0,   0xFF, 0 },
	{ 0,   0xFF,0xFF, 0 },
	{ 0xFF,0xFF,0xFF, 0 }
};

///////////////////////////////////////////////////////////////////////////////
// Global variables

#define STEP	10
#define STEP2	10
//#define MAX_FRAMES 30

int file_loaded;
float x_size, y_size;
char prefix[256];
char sys_msg[256];

extern int cal3D;
extern point3D camera, Prevcamera, vanVec;
extern float camRot;
extern float camRot2;
extern float viewAngle;
extern int YZflag[MAX_TIP_FG_NUM];
extern char name[256];
extern char prefix[256];

point2D vanish;
point2D fg[MAX_TIP_FG_NUM][MAX_FRAMES][11];
point2D TC[MAX_TIP_FG_NUM][MAX_FRAMES][4];

point2D movement;
int isVanish = 0;
int isPnt;
ObjPnt which;
int pcount=1;

int status = NO_STATUS;
int status2 = NO_STATUS;
int fg_num = 0;
int frame_num = 1;
int fis_frame_num = 0;
int cur_fg = 1;

int	snapped; 

int scene_model_type = 0; // default scene_model_type is plane
int panoramic_view_created = 0; // if a parnoramic view has been created

// This is for test
CDC	memDC;

CDC	memDC2; // foreground texture
CDC	memDC3; // foreground mask
CDC	memDC4; // second image
CDC	double_buffer; // double buffer

CBitmap		bitmap;
CBitmap		bitmap2;
CBitmap		bitmap3;

extern float x_size, y_size;

// image size
extern int IMAGE_X, IMAGE_Y;
// view size
extern int view_x, view_y;


BITMAPINFO	m_bmi;
LPVOID		m_pBitmapBits;
HBITMAP	m_hDib;
CSize		m_szPage;

HDC	HmemDC;
HGLRC       memRC;
HDC   m_hOldDC;
HGLRC   m_hOldRC;

CRect       m_oldRect;
float       m_fRadius;

CPalette    m_cPalette2;
CPalette    *m_pOldPalette2;
CRect		m_oldRect2;

string cur_file_name;

/////////////////////////////////////////////////////////////////////////////
// CCubeView

IMPLEMENT_DYNCREATE(CCubeView, CView)

BEGIN_MESSAGE_MAP(CCubeView, CView)
	//{{AFX_MSG_MAP(CCubeView)
	ON_COMMAND(ID_FILE_PLAY, OnFilePlay)
	ON_UPDATE_COMMAND_UI(ID_FILE_PLAY, OnUpdateFilePlay)
	ON_WM_CREATE()
	ON_WM_DESTROY()
	ON_WM_SIZE()
	ON_WM_TIMER()
	ON_WM_ERASEBKGND()
	//}}AFX_MSG_MAP
	ON_COMMAND(ID_FILE_OPEN, OnFileOpen)
	ON_COMMAND(ID_CURVE_FREEHAND, OnCurveFreehand)
	ON_WM_LBUTTONDOWN()
	ON_WM_RBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_RBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_KEYUP()
	//	ON_WM_MBUTTONDOWN()
	ON_WM_KEYDOWN()
	ON_COMMAND(ID_EDGE_PAINT, OnEdgePaint)
	ON_COMMAND(ID_EDGE_BOOST, OnEdgeBoost)
	ON_COMMAND(ID_EDGE_THINNING, OnEdgeThinning)
	ON_COMMAND(ID_EDGE_GOOCH, OnEdgeGooch)
	ON_COMMAND(ID_EDGE_BILATERAL32816, OnEdgeBilateral)
	ON_COMMAND(ID_EDGE_ETF, OnEdgeEtf)
	ON_COMMAND(ID_EDGE_CLD, OnEdgeCld)
	ON_COMMAND(ID_TOON_SHOCK, OnToonShock)
	ON_COMMAND(ID_TOON_WOODCUT2, OnToonWoodcut)
	ON_COMMAND(ID_EDGE_CLDITR, OnEdgeClditr)
	ON_COMMAND(ID_EDGE_WEICKERT, OnEdgeWeickert)
	ON_COMMAND(ID_EDGE_TENSOR, OnEdgeTensor)
	ON_COMMAND(ID_EDGE_FBL, OnEdgeFbl)
	ON_COMMAND(ID_EDGE_FABSTRACT, OnEdgeFabstract)
	ON_COMMAND(ID_EDGE_DOG, OnEdgeDog)
	ON_COMMAND(ID_EDGE_HYBRID, OnEdgeHybrid)
	ON_COMMAND(ID_LEFT, OnKeyLeft)
	ON_COMMAND(ID_RIGHT, OnKeyRight)
	ON_COMMAND(ID_SAVE, OnSave)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CCubeView construction/destruction

CCubeView::CCubeView()
{
	m_pDC = NULL;
	m_pOldPalette = NULL;
	m_play = FALSE;
	test = 12345;
}

CCubeView::~CCubeView()
{
}

/////////////////////////////////////////////////////////////////////////////
// CCubeView drawing

void CCubeView::OnDraw(CDC* pDC)
{
	CCubeDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	//DrawScene();
	TIPDrawScene(pDC);
	CClientDC dc(this);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
}

/////////////////////////////////////////////////////////////////////////////
// CCubeView diagnostics

#ifdef _DEBUG
void CCubeView::AssertValid() const
{
	CView::AssertValid();
}

void CCubeView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CCubeDoc* CCubeView::GetDocument() // non-debug version is inline
{
	return STATIC_DOWNCAST(CCubeDoc, m_pDocument);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CCubeView message handlers

void CCubeView::OnFilePlay()
{
	m_play = m_play ? FALSE : TRUE;
	if (m_play)
		SetTimer(1, 15, NULL);
	else
		KillTimer(1);
}

void CCubeView::OnUpdateFilePlay(CCmdUI* pCmdUI)
{
	pCmdUI->SetCheck(m_play);
}

BOOL CCubeView::PreCreateWindow(CREATESTRUCT& cs)
{
	// An OpenGL window must be created with the following flags and must not
	// include CS_PARENTDC for the class style. Refer to SetPixelFormat
	// documentation in the "Comments" section for further information.
	cs.style |= WS_CLIPSIBLINGS | WS_CLIPCHILDREN;

	return CView::PreCreateWindow(cs);
}

void ComputeStdDev()
{
	double x[5] = { 
1366120.9,
1356020.1,
1368124.6,
1351647.1,
1361258.3
	};
	//double mean = 3606536.5;
	double mean = 0.0;
	for (int i = 0; i < 5; i++) {
		mean += x[i];
	}
	mean /= 5.0;

	double std;
	double var = 0.0;
	for (int i = 0; i < 5; i++) {
		var += (x[i] - mean) * (x[i] - mean);
	}
	var /= 4.0;
	std = sqrt(var);
	TRACE("mean = %f, var = %f, std = %f\n", mean, var, std);
}

int CCubeView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CView::OnCreate(lpCreateStruct) == -1)
		return -1;

	Init(); // initialize OpenGL

	image_mode = OPENGL_MODE; // image render mode initialization
	//image_mode = BITMAP_MODE; // image render mode initialization

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	swap();
	InitMemDCandBitmap();
	CClientDC dc(this);
	InitAlphaMap(dc);

	GetParentFrame()->SetWindowPos(NULL, 0, 0, IMAGE_X+CHILDFRM_RIGHT_BAR_WIDTH, 
		IMAGE_Y+CHILDFRM_TITLE_BAR_WIDTH, SWP_NOMOVE);
	GetParentFrame()->SetWindowText(LPCTSTR("Paint"));

	cur_sigma = 1.0;
	hi_thres = 0.1;
	lo_thres = 0.01;
	MASK_SIZE = 71;


		
	Invalidate(FALSE);
	
	return 0;
}

void CCubeView::OnDestroy()
{
	HGLRC   hrc;

	KillTimer(1);

	hrc = ::wglGetCurrentContext();

	::wglMakeCurrent(NULL,  NULL);

	if (hrc)
		::wglDeleteContext(hrc);

	if (m_pOldPalette)
		m_pDC->SelectPalette(m_pOldPalette, FALSE);

	if (m_pDC)
		delete m_pDC;

	CView::OnDestroy();
}

void CCubeView::OnSize(UINT nType, int cx, int cy)
{
	CView::OnSize(nType, cx, cy);

	if (cy > 0)
	{
		glViewport(0, 0, cx, cy);

		if((m_oldRect.right > cx) || (m_oldRect.bottom > cy))
			RedrawWindow();

		m_oldRect.right = cx;
		m_oldRect.bottom = cy;

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(45.0f, (GLdouble)cx/cy, 3.0f, 7.0f);
		glMatrixMode(GL_MODELVIEW);
	}
	
	CClientDC dc(this);
	TIPDrawScene(&dc);
	//dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
}

void CCubeView::OnTimer(UINT nIDEvent)
{
	CString str_val;

	if (status == EDGE_TUNE) {
		switch (status2) {
			case SCALE_DOWN:
				cur_sigma -= 0.05;
				if (cur_sigma < 0.1) cur_sigma = 0.1;
				str_val.Format("%.1f", cur_sigma);
				dlg->SetDlgItemText(IDC_STATIC4, str_val);
				dlg->m_Scale_Scroll.SetScrollPos((int)(cur_sigma*10));
				break;
			case SCALE_UP:
				cur_sigma += 0.05;
				if (cur_sigma > max_cur_sigma) cur_sigma = max_cur_sigma;
				str_val.Format("%.1f", cur_sigma);
				dlg->SetDlgItemText(IDC_STATIC4, str_val);
				dlg->m_Scale_Scroll.SetScrollPos((int)(cur_sigma*10));
				break;
			case HI_THRES_DOWN:
				hi_thres -= 0.005;
				if (hi_thres < 0.0) hi_thres = 0.0;
				str_val.Format("%.2f", hi_thres);
				dlg->SetDlgItemText(IDC_STATIC_HTHRES, str_val);
				dlg->m_scrollbar_hthres.SetScrollPos((int)(hi_thres*100));
				break;
			case HI_THRES_UP:
				hi_thres += 0.005;
				if (hi_thres > max_hi_thres) hi_thres = max_hi_thres;
				str_val.Format("%.2f", hi_thres);
				dlg->SetDlgItemText(IDC_STATIC_HTHRES, str_val);
				dlg->m_scrollbar_hthres.SetScrollPos((int)(hi_thres*100));
				break;
			case LO_THRES_DOWN:
				lo_thres -= 0.005;
				if (lo_thres < 0.0) lo_thres = 0.0;
				str_val.Format("%.2f", lo_thres);
				dlg->SetDlgItemText(IDC_STATIC_LTHRES, str_val);
				dlg->m_scrollbar_lthres.SetScrollPos((int)(lo_thres*100));
				break;
			case LO_THRES_UP:
				lo_thres += 0.005;
				if (lo_thres > max_lo_thres) lo_thres = max_lo_thres;
				str_val.Format("%.2f", lo_thres);
				dlg->SetDlgItemText(IDC_STATIC_LTHRES, str_val);
				dlg->m_scrollbar_lthres.SetScrollPos((int)(lo_thres*100));
				break;
			case SIZE_DOWN:
				MASK_SIZE -= 1;
				if (MASK_SIZE < MIN_MASK_SIZE) MASK_SIZE = MIN_MASK_SIZE;
				str_val.Format("%d", MASK_SIZE);
				dlg->SetDlgItemText(IDC_STATIC_BRUSH, str_val);
				dlg->m_scrollbar_brush.SetScrollPos((int)(MASK_SIZE));
				break;
			case SIZE_UP:
				MASK_SIZE += 1;
				if (MASK_SIZE > MAX_MASK_SIZE) MASK_SIZE = MAX_MASK_SIZE;
				str_val.Format("%d", MASK_SIZE);
				dlg->SetDlgItemText(IDC_STATIC_BRUSH, str_val);
				dlg->m_scrollbar_brush.SetScrollPos((int)(MASK_SIZE));
				break;
		}
		/////////////////////////////////////////////
		CClientDC dc(this);
		double_buffer.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
		gau_w = MakeGaussMask(cur_sigma, gau);
		max_grad2 = LocalCanny(IMAGE_X, IMAGE_Y, MASK_SIZE, seed, gray, gray2, gau_w);
		//DrawGrayImage(double_buffer, IMAGE_X, IMAGE_Y, gray2);
		DrawGrayImageMask(double_buffer, IMAGE_X, IMAGE_Y, MASK_SIZE, seed, gray2);
		memDC.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &double_buffer, 0, 0, SRCCOPY);
		DrawBoxMemDC(&double_buffer, MASK_SIZE, seed.x, IMAGE_Y-1-seed.y, 255, 0, 255);
		//NonmaxSuppressMask(IMAGE_X, IMAGE_Y);
		dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &double_buffer, 0, 0, SRCCOPY);
	}
	//DrawScene();

	CView::OnTimer(nIDEvent);

	// Eat spurious WM_TIMER messages
	MSG msg;
	while(::PeekMessage(&msg, m_hWnd, WM_TIMER, WM_TIMER, PM_REMOVE));
}

///////////////////////////////////////////////////////////////
void CCubeView::InitMemDCandBitmap()
{
	CClientDC dc(this);

	bitmap.DeleteObject();
	bitmap.CreateCompatibleBitmap(&dc, IMAGE_X, IMAGE_Y);
	memDC.CreateCompatibleDC(&dc);
	memDC.SelectObject(&bitmap);
	
	// init double_buffer
    bitmap2.DeleteObject();
	bitmap2.CreateCompatibleBitmap(&dc, IMAGE_X, IMAGE_Y);
	double_buffer.CreateCompatibleDC(&dc);
	double_buffer.SelectObject(&bitmap2);

	memDC2.CreateCompatibleDC(&dc);
	
	bitmap3.DeleteObject();
	bitmap3.CreateCompatibleBitmap(&dc, TEXTURE, TEXTURE);
	memDC3.CreateCompatibleDC(&dc);
	memDC3.SelectObject(&bitmap3);

	memDC4.CreateCompatibleDC(&dc);
}


/////////////////////////////////////////////////////////////////////////////
// GL helper functions

void CCubeView::Init()
{
	PIXELFORMATDESCRIPTOR pfd;
	int         n;
	HGLRC       hrc;
	GLfloat     fMaxObjSize, fAspect;
	GLfloat     fNearPlane, fFarPlane;

	m_pDC = new CClientDC(this); 

	ASSERT(m_pDC != NULL);

	if (!bSetupPixelFormat())
		return;

	n = ::GetPixelFormat(m_pDC->GetSafeHdc());
	::DescribePixelFormat(m_pDC->GetSafeHdc(), n, sizeof(pfd), &pfd);

	CreateRGBPalette();

	hrc = wglCreateContext(m_pDC->GetSafeHdc());
	wglMakeCurrent(m_pDC->GetSafeHdc(), hrc);

	GetClientRect(&m_oldRect);
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);

	if (m_oldRect.bottom)
		fAspect = (GLfloat)m_oldRect.right/m_oldRect.bottom;
	else    // don't divide by zero, not that we should ever run into that...
		fAspect = 1.0f;
	fNearPlane = 3.0f;
	fFarPlane = 7.0f;
	fMaxObjSize = 3.0f;
	m_fRadius = fNearPlane + fMaxObjSize / 2.0f;

	renderInit();

}

BOOL CCubeView::bSetupPixelFormat()
{
	static PIXELFORMATDESCRIPTOR pfd =
	{
		sizeof(PIXELFORMATDESCRIPTOR),  // size of this pfd
		1,                              // version number
		PFD_DRAW_TO_WINDOW |            // support window
		  //PFD_SUPPORT_OPENGL,          // support OpenGL
		  PFD_SUPPORT_OPENGL |          // support OpenGL
		  PFD_DOUBLEBUFFER,             // double buffered
		PFD_TYPE_RGBA,                  // RGBA type
		24,                             // 24-bit color depth
		0, 0, 0, 0, 0, 0,               // color bits ignored
		0,                              // no alpha buffer
		0,                              // shift bit ignored
		0,                              // no accumulation buffer
		0, 0, 0, 0,                     // accum bits ignored
		32,                             // 32-bit z-buffer
		0,                              // no stencil buffer
		0,                              // no auxiliary buffer
		PFD_MAIN_PLANE,                 // main layer
		0,                              // reserved
		0, 0, 0                         // layer masks ignored
	};
	int pixelformat;

	if ( (pixelformat = ChoosePixelFormat(m_pDC->GetSafeHdc(), &pfd)) == 0 )
	{
		MessageBox("ChoosePixelFormat failed");
		return FALSE;
	}

	if (SetPixelFormat(m_pDC->GetSafeHdc(), pixelformat, &pfd) == FALSE)
	{
		MessageBox("SetPixelFormat failed");
		return FALSE;
	}

	return TRUE;
}

unsigned char CCubeView::ComponentFromIndex(int i, UINT nbits, UINT shift)
{
	unsigned char val;

	val = (unsigned char) (i >> shift);
	switch (nbits)
	{

	case 1:
		val &= 0x1;
		return oneto8[val];
	case 2:
		val &= 0x3;
		return twoto8[val];
	case 3:
		val &= 0x7;
		return threeto8[val];

	default:
		return 0;
	}
}


void CCubeView::CreateRGBPalette()
{
	PIXELFORMATDESCRIPTOR pfd;
	LOGPALETTE *pPal;
	int n, i;

	n = ::GetPixelFormat(m_pDC->GetSafeHdc());
	::DescribePixelFormat(m_pDC->GetSafeHdc(), n, sizeof(pfd), &pfd);

	if (pfd.dwFlags & PFD_NEED_PALETTE)
	{
		n = 1 << pfd.cColorBits;
		pPal = (PLOGPALETTE) new char[sizeof(LOGPALETTE) + n * sizeof(PALETTEENTRY)];

		ASSERT(pPal != NULL);

		pPal->palVersion = 0x300;
		pPal->palNumEntries = n;
		for (i=0; i<n; i++)
		{
			pPal->palPalEntry[i].peRed =
					ComponentFromIndex(i, pfd.cRedBits, pfd.cRedShift);
			pPal->palPalEntry[i].peGreen =
					ComponentFromIndex(i, pfd.cGreenBits, pfd.cGreenShift);
			pPal->palPalEntry[i].peBlue =
					ComponentFromIndex(i, pfd.cBlueBits, pfd.cBlueShift);
			pPal->palPalEntry[i].peFlags = 0;
		}

		/* fix up the palette to include the default GDI palette */
		if ((pfd.cColorBits == 8)                           &&
			(pfd.cRedBits   == 3) && (pfd.cRedShift   == 0) &&
			(pfd.cGreenBits == 3) && (pfd.cGreenShift == 3) &&
			(pfd.cBlueBits  == 2) && (pfd.cBlueShift  == 6)
		   )
		{
			for (i = 1 ; i <= 12 ; i++)
				pPal->palPalEntry[defaultOverride[i]] = defaultPalEntry[i];
		}

		m_cPalette.CreatePalette(pPal);
		delete pPal;

		m_pOldPalette = m_pDC->SelectPalette(&m_cPalette, FALSE);
		m_pDC->RealizePalette();
	}
}

void CCubeView::DrawScene(void)
{
	static BOOL     bBusy = FALSE;
	static GLfloat  wAngleY = 10.0f;
	static GLfloat  wAngleX = 1.0f;
	static GLfloat  wAngleZ = 5.0f;

	if(bBusy)
		return;
	bBusy = TRUE;

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix();

		glTranslatef(0.0f, 0.0f, -m_fRadius);
		glRotatef(wAngleX, 1.0f, 0.0f, 0.0f);
		glRotatef(wAngleY, 0.0f, 1.0f, 0.0f);
		glRotatef(wAngleZ, 0.0f, 0.0f, 1.0f);

		wAngleX += 1.0f;
		wAngleY += 10.0f;
		wAngleZ += 5.0f;


		glBegin(GL_QUAD_STRIP);
			glColor3f(1.0f, 0.0f, 1.0f);
			glVertex3f(-0.5f, 0.5f, 0.5f);

			glColor3f(1.0f, 0.0f, 0.0f);
			glVertex3f(-0.5f, -0.5f, 0.5f);

			glColor3f(1.0f, 1.0f, 1.0f);
			glVertex3f(0.5f, 0.5f, 0.5f);

			glColor3f(1.0f, 1.0f, 0.0f);
			glVertex3f(0.5f, -0.5f, 0.5f);

			glColor3f(0.0f, 1.0f, 1.0f);
			glVertex3f(0.5f, 0.5f, -0.5f);

			glColor3f(0.0f, 1.0f, 0.0f);
			glVertex3f(0.5f, -0.5f, -0.5f);

			glColor3f(0.0f, 0.0f, 1.0f);
			glVertex3f(-0.5f, 0.5f, -0.5f);

			glColor3f(0.0f, 0.0f, 0.0f);
			glVertex3f(-0.5f, -0.5f,  -0.5f);

			glColor3f(1.0f, 0.0f, 1.0f);
			glVertex3f(-0.5f, 0.5f, 0.5f);

			glColor3f(1.0f, 0.0f, 0.0f);
			glVertex3f(-0.5f, -0.5f, 0.5f);

		glEnd();

		glBegin(GL_QUADS);
			glColor3f(1.0f, 0.0f, 1.0f);
			glVertex3f(-0.5f, 0.5f, 0.5f);

			glColor3f(1.0f, 1.0f, 1.0f);
			glVertex3f(0.5f, 0.5f, 0.5f);

			glColor3f(0.0f, 1.0f, 1.0f);
			glVertex3f(0.5f, 0.5f, -0.5f);

			glColor3f(0.0f, 0.0f, 1.0f);
			glVertex3f(-0.5f, 0.5f, -0.5f);
		glEnd();

		glBegin(GL_QUADS);
			glColor3f(1.0f, 0.0f, 0.0f);
			glVertex3f(-0.5f, -0.5f, 0.5f);

			glColor3f(1.0f, 1.0f, 0.0f);
			glVertex3f(0.5f, -0.5f, 0.5f);

			glColor3f(0.0f, 1.0f, 0.0f);
			glVertex3f(0.5f, -0.5f, -0.5f);

			glColor3f(0.0f, 0.0f, 0.0f);
			glVertex3f(-0.5f, -0.5f,  -0.5f);
		glEnd();

	glPopMatrix();

	glFinish();
	SwapBuffers(wglGetCurrentDC());

	bBusy = FALSE;
}

BOOL CCubeView::OnEraseBkgnd(CDC* pDC)
{
	return TRUE;
}

void swap(void)
{
    SwapBuffers(wglGetCurrentDC());
}

void CCubeView::OnFileOpen() 
{
	// TODO: Add your command handler code here
	CString str;
	//CFileDialog dlg(TRUE, "CFG Files", "*.cfg");
	CFileDialog dlg2(TRUE, "BMP Files", "*.bmp");

	if (dlg2.DoModal() == IDOK) {
		status = NO_STATUS;
		
		str = dlg2.GetPathName();

		file_loaded = 1;
		LoadImageBMP((char *)LPCTSTR(str));
		
        int size_x, size_y;

		size_x = IMAGE_X;
		size_y = IMAGE_Y;
		TRACE("%d %d", size_x, size_y);

		GetParentFrame()->SetWindowPos(NULL, 0, 0, size_x+CHILDFRM_RIGHT_BAR_WIDTH, 
			size_y+CHILDFRM_TITLE_BAR_WIDTH, SWP_NOMOVE);
		GetParentFrame()->SetWindowText(LPCTSTR(str));
		AfxGetMainWnd()->SetWindowText("Paint");

		Invalidate(FALSE);

		//////////////////////////////////////////////
		GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);
		GLOBAL_CANNY_DONE = 0; 
		TRACE("cur_sigma = %f\n", cur_sigma);
		TRACE("max_grad2 = %f\n", max_grad2);
		/////////////////////////////////////////////////
		
	}	
	
	CClientDC dc(this);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
}

void CCubeView::LoadImageBMP(char *filename)
// Open BMP file
{
    int x, y;
	GLubyte r, g, b;
//	FILE *fp;

	//////////////////////////////////
	cur_file_name = filename;
	/////////////////////////////////

	LoadBMP((char *)LPCTSTR(filename), &memDC);
	IMAGE_X = m_Image.GetWidth();
	IMAGE_Y = m_Image.GetHeight();

	CClientDC dc(this);

	/// This is for double buffer
	bitmap2.DeleteObject();
	bitmap2.CreateCompatibleBitmap(&dc, IMAGE_X, IMAGE_Y);
	double_buffer.SelectObject(&bitmap2);

	// This is for memDC2
	bitmap3.DeleteObject();
	bitmap3.CreateCompatibleBitmap(&dc, IMAGE_X, IMAGE_Y);
	memDC2.SelectObject(&bitmap3);

	Dbuffer = new GLubyte [IMAGE_Y * IMAGE_X * 3];
	
	// Reverse the vertical for the OpenGL window
	for (y = IMAGE_Y - 1; y >= 0; y--) {
	//for (y = 0; y < IMAGE_Y; y++) {
		for (x = 0; x < IMAGE_X; x++) {
			COLORREF rgbcol = memDC.GetPixel(x, (IMAGE_Y-1)-y);
			r = (GLubyte)RGB_GETRED(rgbcol);
			g = (GLubyte)RGB_GETGREEN(rgbcol);
			b = (GLubyte)RGB_GETBLUE(rgbcol);
			Dbuffer[(y * IMAGE_X + x) * 3 + 0] = r;
			Dbuffer[(y * IMAGE_X + x) * 3 + 1] = g;
			Dbuffer[(y * IMAGE_X + x) * 3 + 2] = b;
		}
	}
}

void CCubeView::OnCurveFreehand()
{
	// TODO: Add your command handler code here
	status = FREEHAND_CURVE;
	ClearMemDC(&memDC);
	TRACE("status = FREEHAND_CURVE\n");
}

void CCubeView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	// TODO: Add your message handler code here and/or call default
	//key_input(nChar); // you can enable repeating key effect while pressing it down

	CView::OnKeyDown(nChar, nRepCnt, nFlags);
}



double cur_sigma;
double max_cur_sigma;
double hi_thres;
double max_hi_thres;
double lo_thres;
double max_lo_thres;
double factor1, factor2;
int	MASK_SIZE;
int MAX_MASK_SIZE, MIN_MASK_SIZE;

void CCubeView::OnEdgePaint()
{
	// TODO: Add your command handler code here
	//GetImage(IMAGE_X, IMAGE_Y, image, Dbuffer);
	if (file_loaded) {
		
		status = EDGE_PAINT;
	}
}

void CCubeView::OnEdgeBoost()
{
	// TODO: Add your command handler code here
		
	status = EDGE_TUNE;
	//OnGlobalCanny(dc, IMAGE_X, IMAGE_Y, gray, gray2);
}





void CCubeView::OnEdgeThinning()
{
	// TODO: Add your command handler code here
	GetImage(IMAGE_X, IMAGE_Y, image, Dbuffer);
	//max_grad = getGradient(IMAGE_X, IMAGE_Y, gradient, image);
	InitGradField(IMAGE_X, IMAGE_Y, gfield, image);
	//////////////////////////////////////////////////////////////////
	// This one is the best
	GVF2(IMAGE_X, IMAGE_Y, gfield, 100, 0.1); // for original picture
	///////////////////////////////////////
	///////////////////////////////////////////////////////
	CClientDC dc(this);
	////////////////////////////////////////////////////////////

	GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);
	//GetSobelGradient(IMAGE_X, IMAGE_Y, gray, gray2, G_mag);
	for (int i = 0; i < IMAGE_X; i++) {
		for (int j = 0; j < IMAGE_Y; j++) {
			G_mag[i][j] = (255-gray[i][j]) / 255.0;
		}
	}
	NonmaximaSuppression(IMAGE_X, IMAGE_Y, gfield, G_mag, gray2, 0.0);
	EdgeThinning(memDC, gray2); // morphological edge thinning
	
	DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray2);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
}



void CCubeView::OnEdgeGooch()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);
	GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);
	GetColorImage(IMAGE_X, IMAGE_Y, cmap, Dbuffer);
	//max_grad = getGradient(IMAGE_X, IMAGE_Y, gradient, image);
	//max_grad = getSobelGradient(IMAGE_X, IMAGE_Y, gradient, image);
	GaussSmooth5(gray, 0.7); // lena
	GetDOG2(gray, G_mag, gray, 1.0, 0.99);  // cameron diaz
	//GetDOG2(gray, G_mag, gray, 1.0, 0.995);  // cameron diaz
	//GetMedian(gray, 1.0, 1);
    Thresholding(IMAGE_X, IMAGE_Y, gray, 0.2); // DOG value 0.2 or below
	//////////////////////////////////////////////////////////
	BilateralColor(IMAGE_X, IMAGE_Y, cmap, 2.0, 0.1, 2);
	ColorQuantization(cmap, 10); 
	ConstructCartoonImage(cmap, gray, cmap);
	DrawColorImage(memDC, IMAGE_X, IMAGE_Y, cmap);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	TRACE("Gooch abstraction!\n");
}



void CCubeView::OnEdgeBilateral()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);
	GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);
	GetColorImage(IMAGE_X, IMAGE_Y, cmap, Dbuffer);
	BilateralColorSep(IMAGE_X, IMAGE_Y, cmap, 2.0, 10.0, 3);
	
	DrawColorImage(memDC, IMAGE_X, IMAGE_Y, cmap);
	//DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	TRACE("Bilateral color filter done!\n");
	CopyCmap2Membuffer(cmap, Dbuffer);

}




void CCubeView::OnToonWoodcut()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);
	GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);
	GaussSmoothSep(gray, 1.0); // lena
	//////////////////////////////////
	int image_x = gray.getRow();
	int image_y = gray.getCol();
	///////////////////////////////////
	double thres = 0.5;
	//double zone_size = 1.2;
	double zone_size = 1.5; // good for human facial portrait!
	imatrix bw(image_x, image_y);
	bw.copy(gray);
	Binarize(bw, thres); // zucker
	imatrix redzone(image_x, image_y);
	redzone.copy(bw);
	CreateRedZone(redzone, zone_size);
	///////////////////////////////////////////////
	DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, bw);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//////////////////////////////////////
	imatrix input(image_x, image_y);
	input.copy(gray);
	//////////////////////////////
	DrawGrayImage(memDC, image_x, image_y, gray);
	ETF e;
	e.init(image_x, image_y);
	e.set(input);
	//e.set(gray2);
	e.Smooth(4, 2); // works best for Yul Brynner!
	//e.Smooth(5, 3);
	//////////////////////////////////////////////////////
	//DrawETF_LIC(dc, e, 20);

	///*
	///////////////////////////////////////////////////////
	// Initial DoG or FDoG
	double tao = 0.99;
	//GetDOG2(gray, G_mag, gray, 1.0, 0.99);  
	GetFDoG(gray, e, 0.7, 0.7, tao); // Initial FDoG. Also works!
	Binarize(gray, 0.5); // zucker
	DrawCondInverseGrayImageRedZone(memDC, gray, bw, redzone);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	/////////////////////////////////////////////////////////
	//*/
	
	/////////////////////////////////////////////////////
	tao = 0.999;
	//////////////////
	for (int i = 0; i < 10; i++) {
		GaussSmoothSep(input, 1.0); 
		bw.copy(input);
		Binarize(bw, thres); 
		//////////////////////////////
		redzone.copy(bw);
		CreateRedZone(redzone, zone_size);
		///////////////////////////////////
		ConstructMergedImage(input, gray, gray); // merge with the original image
		GaussSmoothSep(gray, 1.1); // lena
		///////////////////////////////////////////////////////////////////////
		GetFDoG(gray, e, 0.7, 3.0, tao); // This might produce more realistic woodcut print!
		/////////////////////////////////////////////////////////////////
		Binarize(gray, 0.7); 
		DrawCondInverseGrayImageRedZone(memDC, gray, bw, redzone);
		dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
		TRACE("FDOG iteration %d complete.\n", i+1);
	}
	//*/

	TRACE("\nWoodcut printing completed!\n");
}

void CCubeView::OnToonShock()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);
	GetColorImage(IMAGE_X, IMAGE_Y, cmap, Dbuffer);

	int image_x = cmap.getRow();
	int image_y = cmap.getCol();
	//////////////////////////////
	GaussColSmoothSep(cmap, 1.0);
	CopyCol2GrayImage(image_x, image_y, cmap, gray); 
	DrawColorImage(memDC, image_x, image_y, cmap);
	dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);
	//////////////////////////////////////////////////

	ETF e;
	e.init(image_x, image_y);
	e.set(gray);
	//e.set(gray2);
	e.Smooth(4, 2); // works best for Yul Brynner!

	DrawETF_LIC(dc, e, 10);
	//DrawETFColor360(dc, e);
	//DrawETFColor180(dc, e);
	//DrawETFColorTensor(dc, e);
	
	double tau = 0.99;
	//tau = 0.995;
	double z1, z2;
	z1 = 0.05, z2 = 0.05;
	GetFlowColShockETF(cmap, e, gray, 1.0, 1.0, tau, z1, z2);
	DrawColorImage(memDC, IMAGE_X, IMAGE_Y, cmap);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	z1 = 0.05, z2 = 0.01;
	//tau = 0.99;
	for (int i = 0; i < 10; i++) {
		GaussColSmoothSep(cmap, 1.0);
		CopyCol2GrayImage(IMAGE_X, IMAGE_Y, cmap, gray); 
		GetFlowColShockETF(cmap, e, gray, 1.0, 3.0, tau, z1, z2);
		DrawColorImage(memDC, IMAGE_X, IMAGE_Y, cmap);
		dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
		TRACE("Toon-Shock iteration %d complete.\n", i+1);
	}
	////////////////////////////
	CopyCol2GrayImage(IMAGE_X, IMAGE_Y, cmap, gray); 
	tau = 0.99;
	GetColContrastEnhance(cmap, 0.5, 0.6); 
	DrawColorImage(memDC, IMAGE_X, IMAGE_Y, cmap);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	TRACE("\nEnd of F-Shock Filtering\n");
	//////////////////////////////////////////////////////////////////
}

void CCubeView::OnEdgeWeickert()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);

	GetColorImage(IMAGE_X, IMAGE_Y, cmap, Dbuffer);
	//////////////////////////////////////////////////////////
	GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);
	//////////////////////////////////
	int image_x = gray.getRow();
	int image_y = gray.getCol();

	///////////////////////////////////////////////////////
	imatrix sign(image_x, image_y);
	imatrix gray2(image_x, image_y);
	gray2.copy(gray);
	
	ETF e;
	//////////////
	e.init(image_x, image_y);
	//e.tensor(gray);
	//e.SmoothTensor(gray, 5.0);
	///*
	//////////////////////////////////////////////////////
	//DrawETF_LIC(dc, e, 10);
	//DrawETFColor360(dc, e);
	//DrawETFColor180(dc, e);
	//DrawETFColorTensor(dc, e);
	///////////////////////////////////////////////////////
	//*/
			
	for (int i = 0; i < 30; i++) {
		///////////////////////////////
		///////////////////////////////////////
		//e.init(image_x, image_y);
		e.tensor(gray);
		e.SmoothTensor(gray, 5.0);
		//e.set(gray);
		GaussSmoothSep(gray2, 4.0); // Weickert Mandrill: 4.0
		////////////////////////////////////////////////////////////////////////////
		GetFlowShockDoGETF2(gray2, e, sign, 0.9, 1.0); // Weickert fingerprint: 0.9
		GetWeickertUpwind(cmap, sign, 0.4, 1); // This works!
		DrawColorImage(memDC, image_x, image_y, cmap);
		dc.BitBlt(0, 0, image_x, image_y, &memDC, 0, 0, SRCCOPY);

		////////////////////////////////////////////////////////////////////
		CopyCol2GrayImage2(cmap, gray);  // used to get e
		CopyCol2GrayImage2(cmap, gray2);  // used to get sign (directional second derivative)
				
		TRACE("Weickert shock iteration %d\n", i+1);
	}

	////////////////////
	TRACE("Weickert shock filtering done!\n");
}

void CCubeView::OnEdgeEtf()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);
	//GetImage(IMAGE_X, IMAGE_Y, image, Dbuffer);
	GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);
	GaussSmoothSep(gray, 2.0);
	//////////////////////////////////////////////////
	int image_x = gray.getRow();
	int image_y = gray.getCol();
	DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray);
	//DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray2);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//////////////////////////////////////////////////

	ETF e;
	e.init(image_x, image_y);
	e.set(gray);
	//e.set(gray2);
	e.Smooth(4, 2);
	//////////////////////////////////////////////////////
	DrawETF_LIC(dc, e, 10);
	////////////////////////////////////////////////////////////
	TRACE("ETF construction done!\n");
	//*/
}


void CCubeView::OnEdgeTensor()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);

	GetColorImage(IMAGE_X, IMAGE_Y, cmap, Dbuffer);
	//////////////////////////////////////////////////////////
	GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);
	//////////////////////////////////
	int image_x = gray.getRow();
	int image_y = gray.getCol();

	///////////////////////////////////////////////////////
	imatrix sign(image_x, image_y);
	
	ETF e;
	//////////////
	e.init(image_x, image_y);
	e.tensor(gray);
	e.SmoothTensor(gray, 5.0);
	///*
	//////////////////////////////////////////////////////
	//DrawETF_LIC(dc, e, 10);
	//DrawETFColor360(dc, e);
	//DrawETFColor180(dc, e);
	//DrawETFColorTensor(dc, e);
	//DrawGradientFieldArrow(dc, e, 255, 0, 0);
	//DrawTangentFieldArrow(dc, gfield, 255, 0, 0);
	DrawETFArrow(dc, e, 255, 0, 0);
	///////////////////////////////////////////////////////
	//*/
			
	////////////////////
	TRACE("Structure Tensor done!\n");
}

void CCubeView::OnEdgeFbl()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);
	//GetImage(IMAGE_X, IMAGE_Y, image, Dbuffer);
	GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);
	GetColorImage(IMAGE_X, IMAGE_Y, cmap, Dbuffer);

	//GetPSNR(dc, cmap, "tileable64_noise2.bmp", "texture64_bl_2.0_20_5.bmp");

	int image_x = gray.getRow();
	int image_y = gray.getCol();
	imatrix input(image_x, image_y);
	input.copy(gray);

	////////////////////////////////////////
	DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	////////////////////////////////////////////////
	ETF e;
	e.init(image_x, image_y);
	e.set(input);
	//e.set(gray2);
	e.Smooth(10, 2);
	//////////////////////////////////////////////////////
	DrawETF_LIC(dc, e, 10);
	//DrawETFColor360(dc, e);
	//DrawETFColor180(dc, e);
	//DrawETFColorTensor(dc, e);
	
	GetFBL(cmap, e, 2.0, 2.0, 20.0, 10.0, 3); // small-scale texture
	//BilateralColorSep(image_x, image_y, cmap, 2.0, 10.0, 3);
	//BilateralColor(image_x, image_y, cmap, 1.0, 10, 5);
	//BilateralColorFull(cmap, 2.0, 25, 5);
	//BilateralColorFull(cmap, 2.0, 10, 3);
	
	//ColorQuantization(cmap, 4); 
	//ColorQuantizationLAB(cmap, 4); 
	//ColorQuantization(cmap, 10); 
	//ColorQuantization(cmap, 12); 
	DrawColorImage(memDC, IMAGE_X, IMAGE_Y, cmap);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//////////////////////////////////////////////////////////////////
	CopyCmap2Membuffer(cmap, Dbuffer);
}


void CCubeView::OnEdgeClditr()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);
	GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);
	int image_x = gray.getRow();
	int image_y = gray.getCol();
	imatrix input(image_x, image_y);
	input.copy(gray);

	//////////////////////////////
	DrawGrayImage(memDC, image_x, image_y, gray);
	//DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray2);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//InitGradField(IMAGE_X, IMAGE_Y, gfield, gray);
	//////////////////////////////////////////////////
	//GetSobelGradient2(IMAGE_X, IMAGE_Y, gray, gray2, G_mag); // used for ETF
	ETF e;
	e.init(image_x, image_y);
	e.set(input);
	e.Smooth(4, 2);
	//////////////////////////////////////////////////////
	///////////////////////////////////////////////////////
	// Initial DoG or FDoG
	double tao = 0.99;
	//GetDOG2(gray, G_mag, gray, 1.0, 0.99);  
	//GetDOG2(gray, G_mag, gray, 1.5, 0.99);  
	GetFDoG(gray, e, 0.9, 3.0, tao); // Initial FDoG. Also works!
	//GrayThresholding(gray, 0.7); // zucker
	Binarize(gray, 0.5); // 
	DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	/////////////////////////////////////////////////////////
	
	/////////////////////////////////////////////////////
	///*
	tao = 0.99;
	//double tao = 0.995;
	for (int i = 0; i < 4; i++) {
		ConstructMergedImage(input, gray, gray); // merge with the original image
		//ConstructMergedImageMult(input, gray, gray); // merge with the original image
		//GaussSmooth5(gray, 0.5); // leesy
		GaussSmoothSep(gray, 0.5); // lena
		//GetFlowDOG4(gray, gfield, gray, G_mag, 0.9, 2.0, tao); // Initial FDoG. Also works!
		//GetFDoG(gray, e, 0.7, 2.0, tao); // Initial FDoG. Also works!
		//GetFDoG(gray, e, 0.9, 2.0, tao); // Initial FDoG. Also works!
		GetFDoG(gray, e, 1.0, 2.0, tao); // Initial FDoG. Also works!
		//GetFDoG(gray, e, 3.0, 3.0, tao); // fingerprint
		//GetDoGLIC(gray, e, 1.0, 2.0, tao); // fingerprint
		//GetDoGLIC(gray, e, 1.0, 1.0, tao); // fingerprint
		Binarize(gray, 0.5); 
		//GrayThresholding(gray, 0.7); // zucker
		DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray);
		dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
		TRACE("FDOG iteration %d complete.\n", i+1);
	}
	//*/
	
	///////////////////////////////////////////
	DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//SaveBMPfromMemDC(memDC, IMAGE_X, IMAGE_Y, "_result.bmp");
	TRACE("Line drawing completed!\n");
	//printf("Line drawing completed!\n");
}

void CCubeView::OnEdgeFabstract()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);
	//GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);
	GetColorImage(IMAGE_X, IMAGE_Y, cmap, Dbuffer);
	BilateralColorSep(IMAGE_X, IMAGE_Y, cmap, 1.0, 10.0, 1);
	CopyCol2GrayImage(IMAGE_X, IMAGE_Y, cmap, gray); 

	//cmap.init(image_x, image_y);
	//merged.init(image_x, image_y);
	//gray2.copy(gray);
	GaussSmoothSep(gray, 1.0); // lena
	//GaussSmoothSep(gray, 2.0); // lena
	int image_x = gray.getRow();
	int image_y = gray.getCol();
	imatrix input(image_x, image_y);
	input.copy(gray);
	//////////////////////////////
	DrawGrayImage(memDC, image_x, image_y, gray);
	//DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray2);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//InitGradField(IMAGE_X, IMAGE_Y, gfield, gray);
	//////////////////////////////////////////////////
	//GetSobelGradient2(IMAGE_X, IMAGE_Y, gray, gray2, G_mag); // used for ETF
	ETF e;
	e.init(image_x, image_y);
	e.set(gray);
	e.Smooth(4, 2);
	///////////////////////////////////////////////////////
	// Initial DoG or FDoG
	double tao = 0.99;
	double lambda = 1.0;
	//GetDogSep(gray, 1.0, gray, tao);
	GetDoGLIC(gray, e, 1.0, 2.0, tao); // fingerprint
	Binarize(gray, 0.5); // zucker
	//GrayThresholding(gray, 0.7); // zucker
	//GrayThresholding(gray, 0.5); // zucker
	//Thresholding(IMAGE_X, IMAGE_Y, gray, 0.7); // zucker
	//Thresholding(IMAGE_X, IMAGE_Y, gray, 0.8); // zucker
	DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	/////////////////////////////////////////////////////////
	
	/////////////////////////////////////////////////////
	///*
	//tao = 0.98;
	tao = 0.99;
	//double tao = 0.995;
	for (int i = 0; i < 4; i++) {
		//ConstructMergedImage(input, gray, gray); // merge with the original image
		ConstructMergedImageMult(input, gray, gray); // merge with the original image
		//GaussSmooth5(gray, 0.5); // leesy
		GaussSmoothSep(gray, 0.5); // lena
		//GetFlowDOG4(gray, gfield, gray, G_mag, 0.9, 2.0, tao); 
		//GetFDoG(gray, e, 0.7, 2.0, tao); 
		//GetFDoG(gray, e, 0.9, 2.0, tao); 
		//GetFDoG(gray, e, 1.0, 0.5, tao); 
		GetFDoG(gray, e, 1.0, 3.0, tao); 
		//GetDoGLIC(gray, e, 1.0, 3.0, tao); 
		//GetDoGLIC(gray, e, 1.0, 1.0, tao); 
		Binarize(gray, 0.5); 
		//GrayThresholding(gray, 0.5); // zucker
		DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray);
		dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
		TRACE("FDOG iteration %d complete.\n", i+1);
	}
	//*/

	//GetFBL(cmap, e, 2.0, 2.0, 40.0, 40.0, 1); // 
	//GetFBL(cmap, e, 2.0, 0.5, 30.0, 20.0, 5); // 
	GetFBL(cmap, e, 2.0, 0.5, 40.0, 20.0, 5); // 
	//GetFBL(cmap, e, 2.0, 0.5, 20.0, 20.0, 5); // 
	//GetFBL(cmap, e, 2.0, 0.5, 20.0, 20.0, 5); // 
	//ColorQuantization(cmap, 4); 
	ColorQuantizationLAB(cmap, 4); 
	//ColorQuantizationLAB(cmap, 7); 
	ConstructMergedImageColor(cmap, gray, cmap); // merge with the original image
	
	///////////////////////////////////////////
	//DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray);
	DrawColorImage(memDC, IMAGE_X, IMAGE_Y, cmap);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//SaveBMPfromMemDC(memDC, IMAGE_X, IMAGE_Y, "_result.bmp");
	TRACE("FABTRACT completed!\n");
	//TRACE("Elapsed Time = %f\n", ElapsedTime());
}

void CCubeView::OnEdgeDog()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);
	GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);

	int image_x = gray.getRow();
	int image_y = gray.getCol();
	imatrix input(image_x, image_y);
	input.copy(gray);

	////////////////////////////////////////
	DrawGrayImage(memDC, image_x, image_y, input);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	
	///////////////////////////////////////////////////////
	// Initial DoG or FDoG
	double tao;
	
	//tao = 0.98; // noisy circle
	tao = 0.99;
	//tao = 1.0; // isolated points
	
	GetDogSep(gray, 1.0, gray, tao);
	Binarize(gray, 0.5); // 
	/////////////////////////////////////////////////////////
	
	///////////////////////////////////////////
	DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	SaveBMPfromMemDC(memDC, IMAGE_X, IMAGE_Y, "c:\\temp\\dog.bmp");
	TRACE("DoG Line drawing completed!\n");
}


void CCubeView::OnEdgeCld()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);
	GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);

	int image_x = gray.getRow();
	int image_y = gray.getCol();
	imatrix input(image_x, image_y);
	input.copy(gray);

	////////////////////////////////////////
	DrawGrayImage(memDC, image_x, image_y, input);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	ETF e;
	e.init(image_x, image_y);
	//e.set(gray);
	e.set(input);
	//e.set2(input);
	//e.set(gray2);
	e.Smooth(4, 2);

	///////////////////////////////////////////////////////
	// Initial DoG or FDoG
	double tao;
	
	//tao = 0.98; // noisy circle
	tao = 0.99;
	//tao = 1.0; // isolated points
	
	GetFDoG(gray, e, 1.0, 3.0, tao); 
	//GetFDoG(gray, e, 0.9, 1.0, tao); // Initial FDoG. Also works!
	//GrayThresholding(gray, 0.7); // zucker
	Binarize(gray, 0.5); // 
	/////////////////////////////////////////////////////////
	
	///////////////////////////////////////////
	DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, gray);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	SaveBMPfromMemDC(memDC, IMAGE_X, IMAGE_Y, "c:\\temp\\cld.bmp");
	TRACE("CLD Line drawing completed!\n");
}

void CCubeView::OnEdgeHybrid()
{
	// TODO: Add your command handler code here
	CClientDC dc(this);
	GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);

	int image_x = gray.getRow();
	int image_y = gray.getCol();
	imatrix input(image_x, image_y);
	input.copy(gray);

	////////////////////////////////////////
	DrawGrayImage(memDC, image_x, image_y, input);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	ETF e;
	e.init(image_x, image_y);
	//e.set(gray);

	e.set(input);
	//e.set2(input);
	//e.set(gray2);
	e.Smooth(4, 2);

	///////////////////////////////////////////////////////
	// Initial DoG or FDoG
	double tao;
	
	//tao = 0.98; // noisy circle
	tao = 0.99;
	char name[26];

	imatrix output(image_x, image_y);
	imatrix temp(image_x, image_y);

/*	output.copy(gray);
	GetDoGLIC(output, e, 0.9, 1.0, tao); // DoG + LIC
	Binarize(output, 0.5); // 
	DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, output);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	SaveBMPfromMemDC(memDC, IMAGE_X, IMAGE_Y, "C:\\temp\\hybrid.bmp");
*/
	lambda = 0.5;
	output.copy(gray);
	temp.copy(gray);

	double sigma = 1.0;
	double sigma2 = 3.0;
	imatrix edge(image_x, image_y);
	GetEdgeStructure(e, temp, edge, sigma, tao);
	//Binarize(edge, 0.5); // 
	DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, edge);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	SaveBMPfromMemDC(memDC, IMAGE_X, IMAGE_Y, "C:\\temp\\edge.bmp");
	
	GetDoGLIC2(output, e, edge, sigma, sigma2, tao); // DoG + LIC
	Binarize(output, 0.5); // 
	
/*	for(int i = 0; i < image_x; i++)
	{
		for(int j = 0; j < image_y; j++)
		{
			output[i][j] = min(output[i][j], edge[i][j]);
		}
	}
*/
	DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, output);
	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//sprintf(name, "C:\\temp\\hybrid%f.bmp", lambda);
	SaveBMPfromMemDC(memDC, IMAGE_X, IMAGE_Y, "C:\\temp\\adaptive.bmp");
/*
	for(lambda = 0.0; lambda <= 1.0; ) {
		output.copy(gray);
		GetDoGLIC2(output, e, 1.0, 3.0, tao, lambda); // DoG + LIC
		Binarize(output, 0.5); // 
		DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, output);
		dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
		sprintf(name, "C:\\temp\\hybrid%f.bmp", lambda);
		SaveBMPfromMemDC(memDC, IMAGE_X, IMAGE_Y, name);
//		if(lambda > 0.8) {
//			lambda += 0.005;
//		}else {
			lambda += 0.1;
//		}
	}
*/
	TRACE("Hybrid Line drawing completed!\n");
}

void CCubeView::OnKeyLeft()
{
	if(lambda > 0) {
	//	TRACE("Left!\n");
		lambda -= 0.1;
		if(lambda < 0) lambda = 0;
		CClientDC dc(this);

		GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);

		int image_x = gray.getRow();
		int image_y = gray.getCol();
		imatrix output(image_x, image_y);
		GetDoGLICAdjust(output, lambda); // DoG + LIC
		Binarize(output, 0.5); // 
		DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, output);
		dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//	sprintf(name, "C:\\temp\\hybrid%f.bmp", lambda);
	//	SaveBMPfromMemDC(memDC, IMAGE_X, IMAGE_Y, name);
	//	TRACE("Hybrid Line adjustment completed!\n");
	}
	TRACE("Lambda: %f\n", lambda);
}

void CCubeView::OnKeyRight()
{
	if(lambda < 1.0) {
		lambda += 0.1;
		if(lambda > 1.0) lambda = 1.0;
		CClientDC dc(this);

		GetGrayImage(IMAGE_X, IMAGE_Y, gray, Dbuffer);

		int image_x = gray.getRow();
		int image_y = gray.getCol();
		imatrix output(image_x, image_y);
		GetDoGLICAdjust(output, lambda); // DoG + LIC
		Binarize(output, 0.5); // 
		DrawGrayImage(memDC, IMAGE_X, IMAGE_Y, output);
		dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//	sprintf(name, "C:\\temp\\hybrid%f.bmp", lambda);
	//	SaveBMPfromMemDC(memDC, IMAGE_X, IMAGE_Y, name);
	}
	TRACE("Lambda: %f\n", lambda);
}

void CCubeView::OnSave()
{
	sprintf(name, "C:\\temp\\hybrid%f.bmp", lambda);
	SaveBMPfromMemDC(memDC, IMAGE_X, IMAGE_Y, name);
	TRACE("Saved\n");
}