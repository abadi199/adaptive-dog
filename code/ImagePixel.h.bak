#ifndef _IMPTR_H_
#define _IMPTR_H_


/////////////////////////////////////////////////////////////////////////////
#include "Image.h"

typedef struct tagCOLOR
{
	BYTE B;
	BYTE G;
	BYTE R;
} COLOR;

typedef struct tagCOLOR* LPCOLOR;

class CColorPixel
{
public:
	CColorPixel();
	CColorPixel(COLOR);
	virtual ~CColorPixel();

	inline CColorPixel operator= (COLOR);
	inline CColorPixel operator >>(COLOR &);
	inline CColorPixel operator <<(COLOR);
	inline CColorPixel operator= (CColorPixel);
	inline CColorPixel operator* (CColorPixel);
	inline CColorPixel operator/ (CColorPixel);
	inline CColorPixel operator+ (CColorPixel);
	inline CColorPixel operator- (CColorPixel);
	inline CColorPixel operator*= (CColorPixel);
	inline CColorPixel operator/= (CColorPixel);
	inline CColorPixel operator+= (CColorPixel);
	inline CColorPixel operator-= (CColorPixel);

	inline CColorPixel operator= (int);
	inline CColorPixel operator* (int);
	inline CColorPixel operator/ (int);
	inline CColorPixel operator+ (int);
	inline CColorPixel operator- (int);
	inline CColorPixel operator*= (int);
	inline CColorPixel operator/= (int);
	inline CColorPixel operator+= (int);
	inline CColorPixel operator-= (int);
protected:
	int B, G, R;
};

CColorPixel CColorPixel::operator=(COLOR c)
{
	R = c.R;
	G = c.G;
	B = c.B;
	return *this;
}

CColorPixel CColorPixel::operator <<(COLOR c)
{
	R = c.R;
	G = c.G;
	B = c.B;
	return *this;
}

CColorPixel CColorPixel::operator >>(COLOR &c)
{
	c.R = (BYTE)__max(__min(255, R), 0);
	c.G = (BYTE)__max(__min(255, G), 0);
	c.B = (BYTE)__max(__min(255, B), 0);
	return *this;
}

CColorPixel CColorPixel::operator= (CColorPixel c)
{
	R = c.R;
	G = c.G;
	B = c.B;
	return *this;
}

CColorPixel CColorPixel::operator* (CColorPixel a)
{
	CColorPixel r;
	r.R = R*a.R;
	r.G = G*a.G;
	r.B = B*a.B;
	return r;
}

CColorPixel CColorPixel::operator/ (CColorPixel a)
{
	CColorPixel r;
	r.R = (int)((double)R/a.R+0.5);
	r.G = (int)((double)G/a.G+0.5);
	r.B = (int)((double)B/a.B+0.5);
	return r;
}

CColorPixel CColorPixel::operator+ (CColorPixel a)
{
	CColorPixel r;
	r.R = R+a.R;
	r.G = G+a.G;
	r.B = B+a.B;
	return r;
}

CColorPixel CColorPixel::operator- (CColorPixel a)
{
	CColorPixel r;
	r.R = R-a.R;
	r.G = G-a.G;
	r.B = B-a.B;
	return r;
}

CColorPixel CColorPixel::operator*= (CColorPixel a)
{
	R = R*a.R;
	G = G*a.G;
	B = B*a.B;
	return *this;
}

CColorPixel CColorPixel::operator/= (CColorPixel a)
{
	R = (int)((double)R/a.R+0.5);
	G = (int)((double)G/a.G+0.5);
	B = (int)((double)B/a.B+0.5);
	return *this;
}

CColorPixel CColorPixel::operator+= (CColorPixel a)
{
	R = R+a.R;
	G = G+a.G;
	B = B+a.B;
	return *this;
}

CColorPixel CColorPixel::operator-= (CColorPixel a)
{
	R = R-a.R;
	G = G-a.G;
	B = B-a.B;
	return *this;
}

CColorPixel CColorPixel::operator= (int a)
{
	R = a;
	G = a;
	B = a;
	return *this;
}

CColorPixel CColorPixel::operator* (int a)
{
	R = R*a;
	G = G*a;
	B = B*a;
	return *this;
}

CColorPixel CColorPixel::operator/ (int a)
{
	R = (int)((double)R/a+0.5);
	G = (int)((double)G/a+0.5);
	B = (int)((double)B/a+0.5);
	return *this;
}

CColorPixel CColorPixel::operator+ (int a)
{
	CColorPixel r;
	r.R = R+a;
	r.G = G+a;
	r.B = B+a;
	return r;
}

CColorPixel CColorPixel::operator- (int a)
{
	CColorPixel r;
	r.R = R-a;
	r.G = G-a;
	r.B = B-a;
	return r;
}

CColorPixel CColorPixel::operator*= (int a)
{
	R = R*a;
	G = G*a;
	B = B*a;
	return *this;
}

CColorPixel CColorPixel::operator/= (int a)
{
	R = (int)((double)R/a+0.5);
	G = (int)((double)G/a+0.5);
	B = (int)((double)B/a+0.5);
	return *this;
}

CColorPixel CColorPixel::operator+= (int a)
{
	R = R+a;
	G = G+a;
	B = B+a;
	return *this;
}

CColorPixel CColorPixel::operator-= (int a)
{
	R = R-a;
	G = G-a;
	B = B-a;
	return *this;
}

class CColorPixelPtr
{
public:
	LPCOLOR operator[](LONG index) {
		ASSERT(index>=0 && index<m_nHeight);
		return m_pPtr[index];
	}
	LPCOLOR operator[](int index) {
		ASSERT(index>=0 && index<m_nHeight);
		return m_pPtr[index];
	}
	CColorPixelPtr(CImage &Im);
	CColorPixelPtr(HDIB hHandle);
	virtual ~CColorPixelPtr();

protected:
	LPCOLOR *m_pPtr;
	HDIB m_hHandle;
	int m_nHeight;
};

class CPixel 
{
public:
	int I;
	CPixel();
	CPixel(int);
	virtual ~CPixel();

	inline CPixel operator>> (BYTE &);
	inline CPixel operator<< (BYTE);
	inline CPixel operator<< (int);
	inline operator BYTE() {return I;}
	inline operator int() {return (int)I;}
	inline CPixel operator= (int);
	inline CPixel operator* (int);
	inline CPixel operator/ (int);
	inline CPixel operator+ (int);
	inline CPixel operator- (int);
	inline CPixel operator*= (int);
	inline CPixel operator/= (int);
	inline CPixel operator+= (int);
	inline CPixel operator-= (int);
};

CPixel CPixel::operator >> (BYTE &i)
{
	i = (BYTE)__max(__min(255, I), 0);
	return *this;
}

CPixel CPixel::operator << (BYTE a)
{
	I = a;
	return *this;
}

CPixel CPixel::operator << (int a)
{
	I = a;
	return *this;
}

CPixel CPixel::operator= (int a)
{
	I = a;
	return *this;
}

CPixel CPixel::operator* (int a)
{
	CPixel r;
	r.I = I*a;
	return r;
}

CPixel CPixel::operator/ (int a)
{
	CPixel r;
	r.I = (int)((double)I/a+0.5);
	return r;
}

CPixel CPixel::operator+ (int a)
{
	CPixel r;
	r.I = I+a;
	return r;
}

CPixel CPixel::operator- (int a)
{
	CPixel r;
	r.I = I-a;
	return r;
}

CPixel CPixel::operator*= (int a)
{
	I *= a;
	return *this;
}

CPixel CPixel::operator/= (int a)
{
	I = (int)((double)I/a+0.5);
	return *this;
}

CPixel CPixel::operator+= (int a)
{
	I += a;
	return *this;
}

CPixel CPixel::operator-= (int a)
{
	I -= a;
	return *this;
}

class CPixelPtr
{
public:
	BYTE *operator[](LONG index) {
		ASSERT(index>=0 && index<m_nHeight);
		return m_pPtr[index];
	}
	BYTE *operator[](int index) {
		ASSERT(index>=0 && index<m_nHeight);
		return m_pPtr[index];
	}
	CPixelPtr(CImage &Im);
	CPixelPtr(HDIB hHandle);
	virtual ~CPixelPtr();

protected:
	HDIB m_hHandle;
	BYTE ** m_pPtr;
	int m_nHeight;
};

/////////////////////////////////////////////////////////////////////////////
#endif // _IMPTR_H_