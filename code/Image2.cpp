// Image.cpp : implementation file
//

#include "stdafx.h"
#include "Image.h"
#include "ImagePixel.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE 
static char THIS_FILE[] = __FILE__;
#endif

/******************************************************
				������, �ʱ�ȭ �Լ�
******************************************************/

CImage::CImage()
{
	m_Size		 = CSize(1,1);
	m_hImage	 = NULL;
	m_pPal		 = NULL;
	m_hUndoImage = NULL;
}

CImage::CImage( CImage &Src )
{
	m_Size		 = CSize(1,1);
	m_hImage	 = NULL;
	m_pPal		 = NULL;
	m_hUndoImage = NULL;

	*this = Src;
}

CImage& CImage::operator=( CImage &Src )  // Right side is the argument.
{
	Free();
	m_hImage = (HDIB) ::CopyHandle( Src.GetHandle() );
	m_Size = Src.GetSize();
	m_pPal = new CPalette;
 	CreateDIBPalette();
	return (*this);
}

BOOL CImage::InitDIB(BOOL bCreatePalette)
{
	// �̹����� ����, ���� ũ�� ����
	LPSTR pDIB = (LPSTR)::GlobalLock((HGLOBAL) m_hImage);
	m_Size = CSize((int) ::DIBWidth(pDIB), (int) ::DIBHeight(pDIB));
	::GlobalUnlock((HGLOBAL) m_hImage);

	if(bCreatePalette)
	{
		if(m_pPal != NULL) delete m_pPal;
		// �ȷ�Ʈ ����
		m_pPal = new CPalette;
		if (CreateDIBPalette() == NULL)
		{
			// �ȷ�Ʈ�� ������ �ʴ� ���
			delete m_pPal;
			m_pPal = NULL;
			return FALSE;
		}
	}
	return TRUE;
}

void CImage::SetHandle(HANDLE hHandle)
{
	m_hImage = (HDIB)hHandle;
} 

BOOL CImage::Create(int width, int height, int depth)
{
    LPBITMAPINFOHEADER lpbi ;
	BYTE		*lpPal;
    DWORD       dwSizeImage;
    int         i;

    dwSizeImage = height*(DWORD)((width*depth/8+3)&~3);

	if(depth == 24)
		m_hImage= (HDIB)GlobalAlloc(GHND,sizeof(BITMAPINFOHEADER)+dwSizeImage);
    else
		m_hImage= (HDIB)GlobalAlloc(GHND,sizeof(BITMAPINFOHEADER)+dwSizeImage + 1024);

    if (m_hImage == NULL)
        return FALSE;

	lpbi = (LPBITMAPINFOHEADER)GlobalLock(m_hImage);
	lpbi->biSize            = sizeof(BITMAPINFOHEADER) ;
    lpbi->biWidth           = width;
    lpbi->biHeight          = height;
    lpbi->biPlanes          = 1;
    lpbi->biBitCount        = depth;
    lpbi->biCompression     = BI_RGB ;
    lpbi->biSizeImage       = dwSizeImage;
    lpbi->biXPelsPerMeter   = 0 ;
    lpbi->biYPelsPerMeter   = 0 ;
    lpbi->biClrUsed         = 0 ;
    lpbi->biClrImportant    = 0 ;

	lpPal = (BYTE *) lpbi;
	if (depth == 8)
	{
		lpbi->biClrUsed = 256;

		DWORD offDest = sizeof(BITMAPINFOHEADER);
		for(i = 0; i < 256; i++)
		{
			lpPal[offDest++] = (BYTE)i;
			lpPal[offDest++] = (BYTE)i;
			lpPal[offDest++] = (BYTE)i;
			lpPal[offDest++] = 0x00;
		}                  
	}

	InitDIB(FALSE);
	return TRUE;
}

BOOL CImage::CreateDIBPalette()
{
	LPLOGPALETTE lpPal;      // pointer to a logical palette
	HANDLE hLogPal;          // handle to a logical palette
	HPALETTE hPal = NULL;    // handle to a palette
	int i;                   // loop index
	WORD wNumColors;         // number of colors in color table
	LPSTR lpbi;              // pointer to packed-DIB
	LPBITMAPINFO lpbmi;      // pointer to BITMAPINFO structure (Win3.0)
	LPBITMAPCOREINFO lpbmc;  // pointer to BITMAPCOREINFO structure (old)
	BOOL bWinStyleDIB;       // flag which signifies whether this is a Win3.0 DIB
	BOOL bResult = FALSE;

	/* if handle to DIB is invalid, return FALSE */

	if (m_hImage == NULL)
	  return FALSE;

   lpbi = (LPSTR) ::GlobalLock((HGLOBAL) m_hImage);

   /* get pointer to BITMAPINFO (Win 3.0) */
   lpbmi = (LPBITMAPINFO)lpbi;

   /* get pointer to BITMAPCOREINFO (old 1.x) */
   lpbmc = (LPBITMAPCOREINFO)lpbi;

   /* get the number of colors in the DIB */
   wNumColors = ::DIBNumColors(lpbi);

   if (wNumColors != 0)
   {
		/* allocate memory block for logical palette */
		hLogPal = ::GlobalAlloc(GHND, sizeof(LOGPALETTE)
									+ sizeof(PALETTEENTRY)
									* wNumColors);

		/* if not enough memory, clean up and return NULL */
		if (hLogPal == 0)
		{
			::GlobalUnlock((HGLOBAL) m_hImage);
			return FALSE;
		}

		lpPal = (LPLOGPALETTE) ::GlobalLock((HGLOBAL) hLogPal);

		/* set version and number of palette entries */
		lpPal->palVersion = PALVERSION;
		lpPal->palNumEntries = (WORD)wNumColors;

		/* is this a Win 3.0 DIB? */
		bWinStyleDIB = IS_WIN30_DIB(lpbi);
		for (i = 0; i < (int)wNumColors; i++)
		{
			if (bWinStyleDIB)
			{
				lpPal->palPalEntry[i].peRed = lpbmi->bmiColors[i].rgbRed;
				lpPal->palPalEntry[i].peGreen = lpbmi->bmiColors[i].rgbGreen;
				lpPal->palPalEntry[i].peBlue = lpbmi->bmiColors[i].rgbBlue;
				lpPal->palPalEntry[i].peFlags = 0;
			}
			else
			{
				lpPal->palPalEntry[i].peRed = lpbmc->bmciColors[i].rgbtRed;
				lpPal->palPalEntry[i].peGreen = lpbmc->bmciColors[i].rgbtGreen;
				lpPal->palPalEntry[i].peBlue = lpbmc->bmciColors[i].rgbtBlue;
				lpPal->palPalEntry[i].peFlags = 0;
			}
		}

		/* create the palette and get handle to it */
		bResult = m_pPal->CreatePalette(lpPal);
		::GlobalUnlock((HGLOBAL) hLogPal);
		::GlobalFree((HGLOBAL) hLogPal);
	}

	::GlobalUnlock((HGLOBAL) m_hImage);

	return bResult;
}

/******************************************************
				�Ҹ���, ���� �Լ�
******************************************************/

void CImage::Free()
{
	if( m_hImage )
	{
		if( GlobalFree( m_hImage ) != NULL)
		{
			TRACE("Can't free handle in CImage::Free()");
		}
		m_hImage = NULL;
	}
	if( m_hUndoImage )
	{
		if( GlobalFree( m_hUndoImage ) != NULL)
		{
			TRACE("Can't free handle in CRawImage::Free()");
		}
		m_hUndoImage = NULL;
	}

	if(m_pPal != NULL)
	{
		delete m_pPal;
		m_pPal = NULL;
	}
}

/******************************************************
				�̹��� ������ ��� �Լ�
******************************************************/

int CImage::GetBitCount()
{
	if (m_hImage == NULL) return -1;
	LPBITMAPINFOHEADER lpbi;
	lpbi = (LPBITMAPINFOHEADER) ::GlobalLock((HGLOBAL) m_hImage );
	return lpbi->biBitCount;
}



/******************************************************
						�׸���
******************************************************/

BOOL CImage::Draw(HDC hDC, LPRECT lpDIBRect, LPRECT lpDCRect)
{
	LPSTR	lpDIBHdr;	// BITMAPINFOHEADER�� ����ų ������
	LPSTR	lpDIBBits;	// DIB ��Ʈ�� ����ų ������
	BOOL		bSuccess=FALSE;	 // Success/fail �÷���
	HPALETTE 	hPal=NULL;		 // DIB �ȷ�Ʈ
	HPALETTE 	hOldPal=NULL;	 // ���� �ȷ�Ʈ

	// �޸� ����
	lpDIBHdr  = (LPSTR) ::GlobalLock((HGLOBAL) m_hImage);
	// DIB ��Ʈ�� ����Ǿ� �ִ� ���� �ּҸ� ����
	lpDIBBits = ::FindDIBBits(lpDIBHdr);

	// �ȷ�Ʈ�� ��� DC�� ����
	if(m_pPal != NULL)
	{
		hPal = (HPALETTE) m_pPal->m_hObject;
		hOldPal = ::SelectPalette(hDC, hPal, TRUE);
	}

	::SetStretchBltMode(hDC, COLORONCOLOR);

	if ((RECTWIDTH(lpDCRect)  == RECTWIDTH(lpDIBRect)) &&
	   (RECTHEIGHT(lpDCRect) == RECTHEIGHT(lpDIBRect)))
		// ���� ũ��� �״�� ����ϴ� ���
		bSuccess = ::SetDIBitsToDevice(hDC, // hDC
			lpDCRect->left,		 			// DestX
			lpDCRect->top,		 			// DestY
			RECTWIDTH(lpDCRect),	 		// nDestWidth
			RECTHEIGHT(lpDCRect),			// nDestHeight
			lpDIBRect->left,		 		// SrcX
			(int)DIBHeight(lpDIBHdr) - lpDIBRect->top -	RECTHEIGHT(lpDIBRect),   		// SrcY
			0,                          	// nStartScan
			(WORD)DIBHeight(lpDIBHdr),  	// nNumScans
			lpDIBBits,                  	// lpBits
			(LPBITMAPINFO)lpDIBHdr,			// lpBitsInfo
			DIB_RGB_COLORS);				// wUsage
	else	// Ȯ�� �Ǵ� ����Ͽ� ����ϴ� ���
		bSuccess = ::StretchDIBits(hDC, 	// hDC
			lpDCRect->left,					// DestX
			lpDCRect->top,					// DestY
			RECTWIDTH(lpDCRect),			// nDestWidth
			RECTHEIGHT(lpDCRect),			// nDestHeight
			lpDIBRect->left,				// SrcX
			lpDIBRect->top,					// SrcY
			RECTWIDTH(lpDIBRect),			// wSrcWidth
			RECTHEIGHT(lpDIBRect),			// wSrcHeight
			lpDIBBits,						// lpBits
			(LPBITMAPINFO)lpDIBHdr,			// lpBitsInfo
			DIB_RGB_COLORS,					// wUsage
			SRCCOPY);						// dwROP

	// �޸� ������
   ::GlobalUnlock((HGLOBAL) m_hImage);
	// DC ����
	if (hOldPal != NULL) ::SelectPalette(hDC, hOldPal, TRUE);
	return bSuccess;
}

/******************************************************
				Undo ó���� ���� �Լ�
******************************************************/
BOOL CImage::Undo()
{
	HDIB hTemp;	
	
	// ����� �����ϸ� ����� �̹��� �ڵ��� ��ȯ
	if(m_hUndoImage)
	{
		hTemp = m_hImage;
		m_hImage = m_hUndoImage;
		m_hUndoImage = hTemp;
		return TRUE;
	}
	else
		return FALSE;
}

BOOL CImage::PrepareUndo()
{
	// �̹� ����� �����ϸ� �̸� ������
	if(m_hUndoImage) 
	{
		::GlobalFree((HGLOBAL)m_hUndoImage);
		m_hUndoImage = NULL;
	}

	// �̹����� ��°�� �����Ͽ� ����� ����
	if ((m_hUndoImage = (HDIB)::CopyHandle((HGLOBAL)m_hImage)) == NULL)
	{
		m_hUndoImage = NULL;
		return FALSE;
	}
	return TRUE;
}

/******************************************************
				���� �о����, �����ϱ�
******************************************************/

BOOL CImage::Save(LPCTSTR lpszFileName)
{
	CString filetype;
	filetype = lpszFileName;
	filetype.MakeUpper();

	if(filetype.Find(".BMP") > -1) return SaveBMP(lpszFileName);
	//else if(filetype.Find(".TIF") > -1) return SaveTIF(lpszFileName);
	//else if(filetype.Find(".GIF") > -1) return SaveGIF(lpszFileName);
	//else if(filetype.Find(".JPG") > -1) return SaveJPG(lpszFileName);
	else return FALSE;
}

BOOL CImage::Load(LPCTSTR lpszFileName)
{
	CString filetype;
	filetype = lpszFileName;
	filetype.MakeUpper();

	if(filetype.Find(".BMP") > -1) return LoadBMP(lpszFileName);
	//else if(filetype.Find(".TIF") > -1) return LoadTIF(lpszFileName);
	//else if(filetype.Find(".GIF") > -1) return LoadGIF(lpszFileName);
	//else if(filetype.Find(".JPG") > -1) return LoadJPG(lpszFileName);
	else return FALSE;
}

/******************************************************
  ���� �׷��� ������ �о���� ��ƾ�� ���� ���Ͽ� ����
 
	ImageFileBmp.cpp : BMP ���� (LoadBMP, SaveBMP)
	ImageFileGif.cpp : GIF ���� (LoadGIF, SaveGIF)
	ImageFileTif.cpp : TIFF ���� (LoadTIF, SaveTIF)
	ImageFileJpg.cpp : JPEG ���� (LoadJPG, SaveJPG)

******************************************************/


/******************************************************
				DIB�� ���õ� ���� �Լ�
******************************************************/

LPSTR WINAPI FindDIBBits(LPSTR lpbi)
{
	return (lpbi + *(LPDWORD)lpbi + ::PaletteSize(lpbi));
}


DWORD WINAPI DIBWidth(LPSTR lpDIB)
{
	LPBITMAPINFOHEADER lpbmi;  // pointer to a Win 3.0-style DIB
	LPBITMAPCOREHEADER lpbmc;  // pointer to an other-style DIB

	/* point to the header (whether Win 3.0 and old) */

	lpbmi = (LPBITMAPINFOHEADER)lpDIB;
	lpbmc = (LPBITMAPCOREHEADER)lpDIB;

	/* return the DIB width if it is a Win 3.0 DIB */
	if (IS_WIN30_DIB(lpDIB))
		return lpbmi->biWidth;
	else  /* it is an other-style DIB, so return its width */
		return (DWORD)lpbmc->bcWidth;
}


DWORD WINAPI DIBHeight(LPSTR lpDIB)
{
	LPBITMAPINFOHEADER lpbmi;  // pointer to a Win 3.0-style DIB
	LPBITMAPCOREHEADER lpbmc;  // pointer to an other-style DIB

	/* point to the header (whether old or Win 3.0 */

	lpbmi = (LPBITMAPINFOHEADER)lpDIB;
	lpbmc = (LPBITMAPCOREHEADER)lpDIB;

	/* return the DIB height if it is a Win 3.0 DIB */
	if (IS_WIN30_DIB(lpDIB))
		return lpbmi->biHeight;
	else  /* it is an other-style DIB, so return its height */
		return (DWORD)lpbmc->bcHeight;
}



WORD WINAPI PaletteSize(LPSTR lpbi)
{
   /* calculate the size required by the palette */
   if (IS_WIN30_DIB (lpbi))
	  return (WORD)(::DIBNumColors(lpbi) * sizeof(RGBQUAD));
   else
	  return (WORD)(::DIBNumColors(lpbi) * sizeof(RGBTRIPLE));
}



WORD WINAPI DIBNumColors(LPSTR lpbi)
{
	WORD wBitCount;  // DIB bit count

	/*  If this is a Windows-style DIB, the number of colors in the
	 *  color table can be less than the number of bits per pixel
	 *  allows for (i.e. lpbi->biClrUsed can be set to some value).
	 *  If this is the case, return the appropriate value.
	 */

	if (IS_WIN30_DIB(lpbi))
	{
		DWORD dwClrUsed;

		dwClrUsed = ((LPBITMAPINFOHEADER)lpbi)->biClrUsed;
		if (dwClrUsed != 0)
			return (WORD)dwClrUsed;
	}

	/*  Calculate the number of colors in the color table based on
	 *  the number of bits per pixel for the DIB.
	 */
	if (IS_WIN30_DIB(lpbi))
		wBitCount = ((LPBITMAPINFOHEADER)lpbi)->biBitCount;
	else
		wBitCount = ((LPBITMAPCOREHEADER)lpbi)->bcBitCount;

	/* return number of colors based on bits per pixel */
	switch (wBitCount)
	{
		case 1:
			return 2;

		case 4:
			return 16;

		case 8:
			return 256;

		default:
			return 0;
	}
}


/******************************************************
				Ŭ�����带 ���� ���� �Լ�
******************************************************/

HGLOBAL WINAPI CopyHandle (HGLOBAL h)
{
	if (h == NULL)
		return NULL;

	DWORD dwLen = ::GlobalSize((HGLOBAL) h);
	HGLOBAL hCopy = ::GlobalAlloc(GHND, dwLen);

	if (hCopy != NULL)
	{
		void* lpCopy = ::GlobalLock((HGLOBAL) hCopy);
		void* lp     = ::GlobalLock((HGLOBAL) h);
		memcpy(lpCopy, lp, dwLen);
		::GlobalUnlock(hCopy);
		::GlobalUnlock(h);
	}

	return hCopy;
}



