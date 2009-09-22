// cubeview.h : interface of the CCubeView class
//
// This is a part of the Microsoft Foundation Classes C++ library.
// Copyright (C) 1992-1998 Microsoft Corporation
// All rights reserved.
//
// This source code is only intended as a supplement to the
// Microsoft Foundation Classes Reference and related
// electronic documentation provided with the library.
// See these sources for detailed information regarding the
// Microsoft Foundation Classes product.
#include "gl\gl.h"
#include "gl\glu.h"

#include "Image.h"
#include "defines3.h"

class CCubeView : public CView
{
protected: // create from serialization only
	CCubeView();
	DECLARE_DYNCREATE(CCubeView)

// Attributes
public:
	CCubeDoc* GetDocument();

	CPalette    m_cPalette;
	CPalette    *m_pOldPalette;
	//CRect       m_oldRect;
	//float       m_fRadius;
	CClientDC   *m_pDC;

	CImage m_Image;
	int	test; // for testing access from dialog to CView

// Operations
public:
	void Init();
	void CreateRGBPalette(void);
	BOOL bSetupPixelFormat(void);
	unsigned char ComponentFromIndex(int i, UINT nbits, UINT shift);
	void DrawScene(void);

	//void OnGpaintOpenGL_MemoryDC();
	//void OpenGL_InitMemoryDC();
	//void OpenGL_Init(HDC& m_pDC);
	//void OpenGL_Init(CDC* memDC);
	//void OpenGL_Init(HDC& m_pDC);
	//void OpenGL_CreateRGBPalette(CDC* m_pDC);
	//unsigned char OpenGL_ComponentFromIndex(int i, UINT nbits, UINT shift);
	//bool OpenGL_SetDCPixelFormat(HDC hDC, DWORD dwFlags);
	//BOOL OpenGL_bSetupPixelFormat(CDC* m_pDC);
	int	LoadBMP(CDC& dc, char* filename, CDC *memDC);
	//void OpenGL_InitMemDC(int image_x, int image_y);
	

	void CCubeView::LoadImageBMP(char *filename);
	int	CCubeView::LoadBMP(char* filename, CDC *memDC);
	void CCubeView::SaveFileBMP( const char* fn );
	void CCubeView::OnFileSaveBmp();
	int	CCubeView::LoadBMP2(char* filename, CDC *memDC, int *width, int *height);
	void CCubeView::InitMemDCandBitmap();

	/*
	void CCubeView::OnLButtonDown(UINT nFlags, CPoint point); 
	void CCubeView::OnRButtonDown(UINT nFlags, CPoint point);
	void CCubeView::OnLButtonUp(UINT nFlags, CPoint point);
	void CCubeView::OnRButtonUp(UINT nFlags, CPoint point);
	void CCubeView::OnMouseMove(UINT nFlags, CPoint point);
	*/
	void CCubeView::key_input(UINT nChar);
	void CCubeView::button_pressed(UINT flag, int button, int x, int y);
	void CCubeView::LeftButtonPressed(int x, int y);
	void CCubeView::RightButtonPressed(UINT flag, int x, int y);
	void CCubeView::button_released(int button, int x, int y);
	void CCubeView::mouse_motion(UINT flag, int x, int y);



// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CCubeView)
	public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CCubeView();
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	BOOL    m_play;

// Generated message map functions
protected:
	//{{AFX_MSG(CCubeView)
	afx_msg void OnFilePlay();
	afx_msg void OnUpdateFilePlay(CCmdUI* pCmdUI);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnDestroy();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnTimer(UINT nIDEvent);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnFileOpen();
	afx_msg void OnPaintGeneticpainter();
	afx_msg void OnCurveFreehand();
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnKeyUp(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnEdgePaint();
	afx_msg void OnEdgeBoost();
	afx_msg void OnEdgeTvf();
	afx_msg void OnEdgeFlowedge();
	afx_msg void OnEdgeThinning();
	afx_msg void OnEdgeGabor();
	afx_msg void OnEdgeFlowBilateral();
	afx_msg void OnEdgeAbstract();
	afx_msg void OnEdgeGooch();
	afx_msg void OnEdgeBilateral();
	afx_msg void OnEdgeFlowbilateral2();
	afx_msg void OnEdgeAbstract2();
	afx_msg void OnEdgeFlowbilateral3();
	afx_msg void OnEdgeEtf();
	afx_msg void OnEdgeCld();
	afx_msg void OnToonShock();
	afx_msg void OnToonWoodcut();
	afx_msg void OnEdgeClditr();
	afx_msg void OnEdgeWeickert();
	afx_msg void OnEdgeTensor();
	afx_msg void OnEdgeFbl();
	afx_msg void OnEdgeFabstract();
	afx_msg void OnEdgeDog();
	afx_msg void OnEdgeHybrid();
	afx_msg void OnKeyLeft();
	afx_msg void OnKeyRight();
	afx_msg void OnSave();
};

#ifndef _DEBUG  // debug version in cubeview.cpp
inline CCubeDoc* CCubeView::GetDocument()
   { return (CCubeDoc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////
