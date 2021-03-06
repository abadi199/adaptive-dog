///////////////////////////////////////////
/// include files for avi write
#include <windows.h>
#include <windowsx.h>
#include <memory.h>
#include <mmsystem.h>
#include <vfw.h>

//#include "testview.h"
#include <fstream>
/////////////////////////////////////////////

#include "stdafx.h"
#include "Cube.h"

#include "CubeDoc.h"
#include "CubeView.h"

#include "Image.h"

//#include "gl\gl.h"
//#include "gl\glu.h"

#include "defines.h"
#include "globals.h"

//#include "mathclass.h"

//#include <math.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// definitions from render.c
#include <limits.h>

#define PI 3.1415926535
#define INFI 	10000

//#define FOCAL	200
//extern int FOCAL;

// definitions from CubeView.cpp
////////////////////////////////////////////////////////////////////

#define STEP	10

/////// globals from render.c
extern int cal3D;
extern point3D camera, Prevcamera, vanVec;
extern float camRot;
extern float camRot2;
extern float viewAngle;
extern int YZflag[MAX_TIP_FG_NUM];
extern int movie, fCounter;
extern char name[256];
extern char prefix[256];
extern short int *ibuf;

extern GLuint texName[MAX_TIP_FG_NUM];

//extern point2D TC[11][4];
extern int isVanish, fg_num;

static GLfloat winWidth, winHeight;
extern void swap(void);
extern int file_loaded;
extern GLubyte *Ibuffer;
extern GLubyte *Dbuffer;
extern GLubyte *Maskbuffer;
extern GLubyte *Backbuffer;

extern float x_size, y_size;
extern ObjPnt whichPoint(int obj, int x, int y);
extern point2D ImageCoord[7];

extern hCoord WC[7];
//extern hCoord FG[11][11];
extern hCoord FG[MAX_TIP_FG_NUM][MAX_FRAMES][MAX_TIP_FG_NUM];
extern hCoord static_FG[MAX_TIP_FG_NUM][MAX_FRAMES][MAX_TIP_FG_NUM];

extern point3D Wc, ray;
extern point2D Vc;
extern point3D v[MAX_TIP_FG_NUM];
extern float half_x, half_y;

extern void CalImageCoord();
extern void CalWorldCoord();
extern void renderNone(void);

//////// globals defined in CubeView.cpp
extern point2D vanish;
//extern point2D fg[11][11];
extern point2D movement;
extern int set_mode;
extern int isVanish, isPnt;
extern ObjPnt which;
extern int pcount;

extern int status;
extern int fg_num;

// image size
extern int IMAGE_X, IMAGE_Y;
// view size
extern int view_x, view_y;

// functions defined in other files
extern void renderInit(void);
//extern void renderInit(CDC& dc);
extern void renderReshape(int width, int height);
extern void renderScene();
extern char* fileLoad(char* filename);
extern void LoadRXYFile(CDC *, char *);
//extern void button_pressed(int button, int x, int y);
//extern void button_released(int button, int x, int y);
//extern void key_input(UINT nChar);
//extern void mouse_motion(UINT flag, int x, int y);


//BOOL SaveBMP( BITMAPINFO&, const char*, const char* );

#define WIDTHBYTES(bits)   (((bits) + 31) / 32 * 4)
#define DIB_HEADER_MARKER   ((WORD) ('M' << 8) | 'B')

BOOL SaveBMP(BITMAPINFO& bi, const char* img, const char* lpszFileName)
{
	FILE* file;
	BITMAPFILEHEADER bmfHdr;

	file = fopen( lpszFileName, "wb+" );

	bmfHdr.bfType = DIB_HEADER_MARKER;  // "BM"
	bmfHdr.bfSize=sizeof(BITMAPINFOHEADER)+bi.bmiHeader.biSizeImage+sizeof(BITMAPFILEHEADER);

	bmfHdr.bfReserved1 = 0;
	bmfHdr.bfReserved2 = 0;
	bmfHdr.bfOffBits=(DWORD)sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER);

	fwrite( (LPSTR)&bmfHdr, sizeof(BITMAPFILEHEADER),1,file);
	fwrite( &(bi.bmiHeader), sizeof(BITMAPINFOHEADER), 1, file );
	fwrite( img, bi.bmiHeader.biSizeImage, 1, file );

	fclose( file );
	return TRUE;
}

void bitmapInfo(BITMAPINFO& DIBInfo, int& dx, int& dy )
{
	//int ww = (viewer->wi/8)*8, hh = (viewer->he/8)*8; // width and height
	int ww = (IMAGE_X/8)*8, hh = (IMAGE_Y/8)*8; // width and height

	//dx=(ww+3)&(0xfffffffc); // why + 3 ?
	//dy=(hh+3)&(0xfffffffc); // why + 3 ?
	dx=(ww+0)&(0xfffffffc); // why + 3 ?
	dy=(hh+0)&(0xfffffffc); // why + 3 ?

	DIBInfo.bmiHeader.biSize		= 40;
	DIBInfo.bmiHeader.biWidth 		= dx;
	DIBInfo.bmiHeader.biHeight		= dy;
	DIBInfo.bmiHeader.biPlanes		= 1;
	DIBInfo.bmiHeader.biBitCount	= 24;
	DIBInfo.bmiHeader.biCompression	= BI_RGB;
	DIBInfo.bmiHeader.biSizeImage 	= dy*(DWORD)((dx*3+3)&~3);
	DIBInfo.bmiHeader.biXPelsPerMeter;
	DIBInfo.bmiHeader.biClrUsed		= 0;
	DIBInfo.bmiHeader.biClrImportant= 0;
}

void CCubeView::SaveFileBMP( const char* fn )
{
	char *img = NULL;
	HBITMAP hbmp;
	BITMAPINFO 	DIBInfo;
	int dx,dy;

	HDC mmdc;

	//CClientDC	dc(this);
	//memDC.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &dc, 0, 0, SRCCOPY);
	//dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	
	//draw();
	//GetWindowRect( viewer->hWnd, &rect );
	//viewer->g_pddsPrimary->GetDC( &gdc );

	//bitmapInfo( viewer, DIBInfo, dx, dy );
	bitmapInfo( DIBInfo, dx, dy );
	//hbmp= CreateDIBSection( gdc, (BITMAPINFO*)&DIBInfo, DIB_RGB_COLORS, (void**)&img, NULL, 0);
	hbmp = CreateDIBSection( memDC, (BITMAPINFO*)&DIBInfo, DIB_RGB_COLORS, (void**)&img, NULL, 0);

	//mmdc=CreateCompatibleDC(gdc);
	mmdc = CreateCompatibleDC(memDC);
	HGLOBAL ob = SelectObject(mmdc, hbmp); // the bmp object is now contained in mmdc
	//BitBlt(mmdc, 0, 0, dx, dy, gdc, rect.left, rect.top, SRCCOPY); // from gdc -> mmdc : the direction is important!
	//dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//BitBlt(mmdc, 0, 0, dx, dy, dc, rect.left, rect.top, SRCCOPY); // f
	//BitBlt(mmdc, 0, 0, dx, dy, memDC, rect.left, rect.top, SRCCOPY); // from memDC -> mmdc
	BitBlt(mmdc, 0, 0, dx, dy, memDC, 0, 0, SRCCOPY); // from memDC -> mmdc
	//viewer->g_pddsPrimary->ReleaseDC( gdc );

	if ( fn==NULL )	SaveBMP( DIBInfo, img, "TEST.BMP" );
	else			SaveBMP( DIBInfo, img, fn );
	//SaveBMP( DIBInfo, img, "TEST.BMP" );

	SelectObject(mmdc,ob);
	DeleteDC(mmdc);
	DeleteObject( hbmp );
}

//BOOL CImage::SaveBMP(LPCTSTR lpszFileName)
//BOOL CImage::LoadBMP(LPCTSTR lpszFileName)

void CCubeView::OnFileSaveBmp() 
{
	// TODO: Add your command handler code here
	//SaveFileBMP("test.bmp");
	//SaveFileBMP2(IMAGE_X, IMAGE_Y, "test.bmp");
	TRACE("I'm saving bmp\n");
	//SaveFileBMP2(IMAGE_X, IMAGE_Y, "a.bmp");
	SaveFileBMP2(IMAGE_X, IMAGE_Y, "test.bmp");
}

void SaveFileBMP2( int width, int height, const char* fn )
// Save memDC contents to file
{
	char *img = NULL;
	HBITMAP hbmp;
	BITMAPINFO 	DIBInfo;
	int dx,dy;

	HDC mmdc;

	//int ww = (width/8)*8, hh = (height/8)*8; // width and height
	int ww = width, hh = height; // width and height

	//dx=(ww+3)&(0xfffffffc); // why + 3 ?
	//dy=(hh+3)&(0xfffffffc); // why + 3 ?
	dx=(ww+0)&(0xfffffffc); 
	dy=(hh+0)&(0xfffffffc); 

	DIBInfo.bmiHeader.biSize		= 40;
	DIBInfo.bmiHeader.biWidth 		= dx;
	DIBInfo.bmiHeader.biHeight		= dy;
	DIBInfo.bmiHeader.biPlanes		= 1;
	DIBInfo.bmiHeader.biBitCount	= 24;
	DIBInfo.bmiHeader.biCompression	= BI_RGB;
	DIBInfo.bmiHeader.biSizeImage 	= dy*(DWORD)((dx*3+3)&~3);
	DIBInfo.bmiHeader.biXPelsPerMeter;
	DIBInfo.bmiHeader.biClrUsed		= 0;
	DIBInfo.bmiHeader.biClrImportant= 0;

	//bitmapInfo( DIBInfo, dx, dy );
	//hbmp= CreateDIBSection( gdc, (BITMAPINFO*)&DIBInfo, DIB_RGB_COLORS, (void**)&img, NULL, 0);
	hbmp = CreateDIBSection( memDC, (BITMAPINFO*)&DIBInfo, DIB_RGB_COLORS, (void**)&img, NULL, 0);

	//mmdc=CreateCompatibleDC(gdc);
	mmdc = CreateCompatibleDC(memDC);
	HGLOBAL ob = SelectObject(mmdc, hbmp); // the bmp object is now contained in mmdc
	//BitBlt(mmdc, 0, 0, dx, dy, gdc, rect.left, rect.top, SRCCOPY); // from gdc -> mmdc : the direction is important!
	//dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//BitBlt(mmdc, 0, 0, dx, dy, dc, rect.left, rect.top, SRCCOPY); // f
	//BitBlt(mmdc, 0, 0, dx, dy, memDC, rect.left, rect.top, SRCCOPY); // from memDC -> mmdc
	BitBlt(mmdc, 0, 0, dx, dy, memDC, 0, 0, SRCCOPY); // from memDC -> mmdc
	//viewer->g_pddsPrimary->ReleaseDC( gdc );

	if ( fn==NULL )	SaveBMP( DIBInfo, img, "TEST.BMP" );
	else			SaveBMP( DIBInfo, img, fn );
	//SaveBMP( DIBInfo, img, "TEST.BMP" );

	SelectObject(mmdc,ob);
	DeleteDC(mmdc);
	DeleteObject( hbmp );
}

void SaveBMPfromMemDC( CDC& dc, int width, int height, const char* fn )
{
	char *img = NULL;
	HBITMAP hbmp;
	BITMAPINFO 	DIBInfo;
	int dx,dy;

	HDC mmdc;

	//int ww = (width/8)*8, hh = (height/8)*8; // width and height
	int ww = width, hh = height; // width and height

	//dx=(ww+3)&(0xfffffffc); // why + 3 ?
	//dy=(hh+3)&(0xfffffffc); // why + 3 ?
	dx=(ww+0)&(0xfffffffc); 
	dy=(hh+0)&(0xfffffffc); 

	DIBInfo.bmiHeader.biSize		= 40;
	DIBInfo.bmiHeader.biWidth 		= dx;
	DIBInfo.bmiHeader.biHeight		= dy;
	DIBInfo.bmiHeader.biPlanes		= 1;
	DIBInfo.bmiHeader.biBitCount	= 24;
	DIBInfo.bmiHeader.biCompression	= BI_RGB;
	DIBInfo.bmiHeader.biSizeImage 	= dy*(DWORD)((dx*3+3)&~3);
	DIBInfo.bmiHeader.biXPelsPerMeter;
	DIBInfo.bmiHeader.biClrUsed		= 0;
	DIBInfo.bmiHeader.biClrImportant= 0;

	//bitmapInfo( DIBInfo, dx, dy );
	//hbmp= CreateDIBSection( gdc, (BITMAPINFO*)&DIBInfo, DIB_RGB_COLORS, (void**)&img, NULL, 0);
	hbmp = CreateDIBSection( dc, (BITMAPINFO*)&DIBInfo, DIB_RGB_COLORS, (void**)&img, NULL, 0);

	//mmdc=CreateCompatibleDC(gdc);
	mmdc = CreateCompatibleDC(dc);
	HGLOBAL ob = SelectObject(mmdc, hbmp); // the bmp object is now contained in mmdc
	//BitBlt(mmdc, 0, 0, dx, dy, gdc, rect.left, rect.top, SRCCOPY); // from gdc -> mmdc : the direction is important!
	//dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
	//BitBlt(mmdc, 0, 0, dx, dy, dc, rect.left, rect.top, SRCCOPY); // f
	//BitBlt(mmdc, 0, 0, dx, dy, memDC, rect.left, rect.top, SRCCOPY); // from memDC -> mmdc
	BitBlt(mmdc, 0, 0, dx, dy, dc, 0, 0, SRCCOPY); // from memDC -> mmdc
	//viewer->g_pddsPrimary->ReleaseDC( gdc );

	if ( fn==NULL )	SaveBMP( DIBInfo, img, "TEST.BMP" );
	else			SaveBMP( DIBInfo, img, fn );
	//SaveBMP( DIBInfo, img, "TEST.BMP" );

	SelectObject(mmdc,ob);
	DeleteDC(mmdc);
	DeleteObject( hbmp );
}

///*
int	CCubeView::LoadBMP(char* filename, CDC *memDC)
// Open BMP File
{
	int		m_nMul, m_nDiv;

	CClientDC dc(this);

	m_nMul = m_nDiv = 1;
	//	file_loaded = 1;
	if ( m_Image.Load(filename) == FALSE ) // the file does not exist
		return 0;

	//IMAGE_X = m_Image.GetWidth();
	bitmap.DeleteObject();
	bitmap.CreateCompatibleBitmap(&dc, m_Image.GetWidth(), m_Image.GetHeight());
	//memDC.CreateCompatibleDC(&dc);
	memDC->SelectObject(&bitmap);
	//	memDC.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &dc, 0, 0, SRCCOPY);
	//	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);

	CRect rcDIB;
	rcDIB.top = rcDIB.left = 0;
	rcDIB.right = m_Image.GetWidth();
	rcDIB.bottom = m_Image.GetHeight();
	CRect rcDest;
	rcDest = rcDIB;
	rcDest.left = rcDest.left * m_nMul / m_nDiv;
	rcDest.top = rcDest.top * m_nMul / m_nDiv;
	rcDest.right= rcDest.right * m_nMul / m_nDiv;
	rcDest.bottom= rcDest.bottom * m_nMul / m_nDiv;

		//m_Image.Draw(pDC->m_hDC,&rcDIB,&rcDest);
		// Draw the BMP image into a memDC!!!
		
		//m_Image.Draw(dc.m_hDC, &rcDIB, &rcDest);
		//memDC->BitBlt(0, 0, m_Image.GetWidth(), m_Image.GetHeight(), &dc, 0, 0, SRCCOPY);

	// Draw the BMP image on memDC
	m_Image.Draw(memDC->m_hDC, &rcDIB, &rcDest);
//		dc.BitBlt(0, 0, m_Image.GetWidth(), m_Image.GetHeight(), 
//			memDC, 0, 0, SRCCOPY);

	//OpenGL_InitMemoryDC();

	return 1; // file does exist
}
//*/

int	CCubeView::LoadBMP2(char* filename, CDC *memDC, int *width, int *height)
// Open BMP File and width and height
{
	int		m_nMul, m_nDiv;

	CClientDC dc(this);

	m_nMul = m_nDiv = 1;
	//	file_loaded = 1;
	if ( m_Image.Load(filename) == FALSE ) // the file does not exist
		return 0;

	//IMAGE_X = m_Image.GetWidth();
	bitmap.DeleteObject();
	bitmap.CreateCompatibleBitmap(&dc, m_Image.GetWidth(), m_Image.GetHeight());
	//memDC.CreateCompatibleDC(&dc);
	memDC->SelectObject(&bitmap);
	//	memDC.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &dc, 0, 0, SRCCOPY);
	//	dc.BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);

	CRect rcDIB;
	rcDIB.top = rcDIB.left = 0;
	rcDIB.right = m_Image.GetWidth();
	rcDIB.bottom = m_Image.GetHeight();
	CRect rcDest;
	rcDest = rcDIB;
	rcDest.left = rcDest.left * m_nMul / m_nDiv;
	rcDest.top = rcDest.top * m_nMul / m_nDiv;
	rcDest.right= rcDest.right * m_nMul / m_nDiv;
	rcDest.bottom= rcDest.bottom * m_nMul / m_nDiv;

		//m_Image.Draw(pDC->m_hDC,&rcDIB,&rcDest);
		// Draw the BMP image into a memDC!!!
			/*
		m_Image.Draw(dc.m_hDC, &rcDIB, &rcDest);
		memDC->BitBlt(0, 0, m_Image.GetWidth(), m_Image.GetHeight(), 
			&dc, 0, 0, SRCCOPY);
			*/

	// Draw the BMP image on memDC
	m_Image.Draw(memDC->m_hDC, &rcDIB, &rcDest);
//		dc.BitBlt(0, 0, m_Image.GetWidth(), m_Image.GetHeight(), 
//			memDC, 0, 0, SRCCOPY);
	*width = m_Image.GetWidth();
	*height = m_Image.GetHeight();

	return 1; // file does exist
}



/*
void CCubeView::OnFileInvertMask() 
{
	// TODO: Add your command handler code here
	CString str;
	CFileDialog dlg(TRUE, "TXT Files", "*.txt");

	if (dlg.DoModal() == IDOK) {
		str = dlg.GetPathName();
		//fileLoad((char *)LPCTSTR(str));
		//LoadFile((char *)LPCTSTR(str));

		InvertMASK((char *)LPCTSTR(str));
        
	}	
	
}
*/

