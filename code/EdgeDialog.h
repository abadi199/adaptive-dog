#pragma once


// EdgeDialog dialog

class EdgeDialog : public CDialog
{
	DECLARE_DYNAMIC(EdgeDialog)

public:
	EdgeDialog(CWnd* pParent = NULL);   // standard constructor
	virtual ~EdgeDialog();

// Dialog Data
	enum { IDD = IDD_DIALOG2 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
};
