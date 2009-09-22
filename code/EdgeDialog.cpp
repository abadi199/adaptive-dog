// EdgeDialog.cpp : implementation file
//

#include "stdafx.h"
#include "Cube.h"
#include "EdgeDialog.h"


// EdgeDialog dialog

IMPLEMENT_DYNAMIC(EdgeDialog, CDialog)
EdgeDialog::EdgeDialog(CWnd* pParent /*=NULL*/)
	: CDialog(EdgeDialog::IDD, pParent)
{
}

EdgeDialog::~EdgeDialog()
{
}

void EdgeDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(EdgeDialog, CDialog)
END_MESSAGE_MAP()


// EdgeDialog message handlers
