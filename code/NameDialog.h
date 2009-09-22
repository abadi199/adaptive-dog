#pragma once
#include "afxwin.h"
#include "afxcmn.h"

// NameDialog dialog

class NameDialog : public CDialog
{
	DECLARE_DYNAMIC(NameDialog)

public:
	NameDialog(CWnd* pParent = NULL);   // standard constructor
	virtual ~NameDialog();

// Dialog Data
	enum { IDD = IDD_DIALOG1 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL OnInitDialog();
	//afx_msg void OnBnClickedButton1();
	bool m_Paint;
	afx_msg void OnBnClickedBoost();
	afx_msg void OnBnClickedPaint();
	//afx_msg void OnBnClickedConnect();
	//CEdit m_scale_val;
	CScrollBar m_Scale_Scroll;
	//CButton m_Paint_radio;
	//CButton m_Paint_radio2; // check box
	
	afx_msg void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	//CSliderCtrl m_slider;
	afx_msg void OnBnClickedGlobalCanny();
	int m_paint_group;
	int m_scale_group;
	afx_msg void OnBnClickedScale();
	afx_msg void OnBnClickedHthres();
	CScrollBar m_scrollbar_hthres;
	CScrollBar m_scrollbar_brush;
	afx_msg void OnBnClickedBrush();
	CScrollBar m_scrollbar_lthres;
	afx_msg void OnBnClickedLthres();
	afx_msg void OnBnClickedSample();
	afx_msg void OnBnClickedPick();
	afx_msg void OnBnClickedSamplePixels();
	afx_msg void OnBnClickedSegment();
	afx_msg void OnBnClickedNonmax();
};
