#ifndef STROKE_H
#define STROKE_H

class Stroke {	
public:
	int	xc, yc;
	double	theta;
	int	w; // width
	int	st_x, st_y, end_x, end_y; // st points and ending points
	unsigned char		r, g, b; // default : private
	POINT P[4]; // for tile mosaic!
public:
	//float	xc, yc, theta, w, l;
	Stroke() { xc = yc = 0; theta = 0; w = 0; r = g = b = 0; }
	Stroke(int xc, int yc, int w, int st_x, int st_y, int end_x, int end_y, unsigned char r, unsigned char g, unsigned char b) 
	{ 
		this->xc = xc;
		this->yc = yc;
		this->theta = 0;
		this->w = w;	
		this->st_x = st_x;	
		this->st_y = st_y;
		this->end_x = end_x;
		this->end_y = end_y;
		this->r = r;
		this->g = g;
		this->b = b;
	}
	Stroke(int xc, int yc, int w, POINT* P, unsigned char r, unsigned char g, unsigned char b) 
	{ 
		this->xc = xc;
		this->yc = yc;
		this->w = w;	
		this->r = r;
		this->g = g;
		this->b = b;
		this->P[0].x = P[0].x;
		this->P[0].y = P[0].y;
		this->P[1].x = P[1].x;
		this->P[1].y = P[1].y;
		this->P[2].x = P[2].x;
		this->P[2].y = P[2].y;
		this->P[3].x = P[3].x;
		this->P[3].y = P[3].y;
	}
	void Draw(CDC& dc);
	bool DrawMosaic(CDC& dc, int overlap_test); 
	void DrawTileWhite(CDC& dc);
	bool DrawTileFast(CDC& dc, int overlap_test); 
	void SetValues(int xc2, int yc2, int w2, int st_x2, int st_y2, int end_x2, int end_y2, unsigned char r2, unsigned char g2, unsigned char b2);
	void SetValues(int xc2, int yc2, int w2, POINT* P, unsigned char r2, unsigned char g2, unsigned char b2);
	int GetXc() { return xc; }
	int GetYc() { return yc; }
	double GetTheta() { return theta; }
	int GetW() { return w; }
	//int	GetLL() { return ll; }
	//int GetRL() { return rl; }
	unsigned char GetR() { return r; }
	unsigned char GetG() { return g; }
	unsigned char GetB() { return b; }
};

#endif

