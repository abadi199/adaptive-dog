#ifndef _FIELD_H_
#define _FIELD_H_

//#define MAX_PTS	500
//#define MAX_CNT_PTS	500

//class vector;

#include "ivectormatrix.h"

struct Grad {
	double gx, gy, mag;
};


class Field {
private:
	int Nr, Nc;
	//double** p = new int*[rows];
    //vector p[3];
	//vector* p; // a sequence of row vectors
	Grad** p; 
	double max_grad;

    // addition & subtraction
    //friend matrix    operator+( matrix const&, matrix const& );
public:
    // constructors
	Field() 
    {
		Nr = 1, Nc = 1;
		//p = new vector[3];
		p = new Grad*[Nr];
		for(int i = 0; i < Nr; i++)
		   p[i] = new Grad[Nc];
		// Create 3x3 identity matrix
        p[0][0].gx=1.0; p[0][0].gy=1.0; p[0][0].mag=1.0;
        max_grad = 1.0;
    };
	Field(int i, int j) 
    {
		Nr = i, Nc = j;
		//p = new vector[i];
		p = new Grad*[Nr];
		for(i = 0; i < Nr; i++)
		   p[i] = new Grad[Nc];
		// Create ixj identity matrix
        //p[0][0]=1.0; p[0][1]=0.0; p[0][2]=0.0;
        //p[1][0]=0.0; p[1][1]=1.0; p[1][2]=0.0;
        //p[2][0]=0.0; p[2][1]=0.0; p[2][2]=1.0;
    };
	void ClearAll() {
		//for (int i = 0; i < Nr; i++) 
		//	delete[] p[i];
		delete[] p;
	}
	~Field() { ClearAll(); }
	////////////////////////////////////////////
	/// matrix[i] returns the l-value of p[i] (ith row vector)
	// Thus, matrix[i][j] returns the l-value of p[i][j]
	// Thus, matrix[i][j] = 1.0 is possible!
	//vector& operator[](int i) { return p[i]; };
	///////////////////////////////////////////////
	// Very important!
	Grad* operator[](int i) { return p[i]; };
	// This enables BOTH of the following expressions:
	// A[i][j] = 1.5;, tmp = A[i][j];
	//double*& operator[](int i) { return p[i]; }; This also works, but redundant!

	//doublegetValue( int i, int j ) const { return p[i].getValue( j );}
	Grad& get( int i, int j ) const { return p[i][j]; }
	int getRow() const { return Nr; }
	int getCol() const { return Nc; }
	void init(int i, int j) 
    {
		Nr = i, Nc = j;
		//p = new vector[i];
		//delete[] p;
		ClearAll();
		p = new Grad*[Nr];
		for(i = 0; i < Nr; i++)
		   p[i] = new Grad[Nc];
		max_grad = 1.0;
    };
	void copy(Field& s) 
    {
		//init(s.Nr, s.Nc);
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) {
				p[i][j].gx = s.p[i][j].gx;
				p[i][j].gy = s.p[i][j].gy;
				p[i][j].mag = s.p[i][j].mag;
			}
		max_grad = s.max_grad;
    };
	void zero()
	{
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) 
				p[i][j].gx = p[i][j].gy = p[i][j].mag = 0.0;
	}
	//void set(Image& image, Image& gradient); 
	void set(Image& image); 
	void set(imatrix& image); 
	Field& operator=( Field& s ) {
		//ClearAll();
		Nr = s.getRow(); Nc = s.getCol();
		init(Nr, Nc);
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) {
				p[i][j].gx = s[i][j].gx;
				p[i][j].gy = s[i][j].gy;
				p[i][j].mag = s[i][j].mag;
			}
		max_grad = s.max_grad;
		return *this;			
	}
	double GetMaxGrad() { return max_grad; }
	void normalize(); 
	
	//friend void ludcmp(matrix& a, int_vector& indx, double &d);
	//friend void lubksb(matrix& a, int_vector& indx, vector &b);
	//friend void inverse(matrix &a, matrix &inverse);
	//friend void merge_h(matrix &a, matrix &b, matrix &merged);
	//friend void split_v(matrix &merged, int a_row, matrix &a, matrix &b);
	//void print();
	//matrix inverse(); // This is the better version
	//matrix copy();
	//void copy(matrix& b);

	//friend vector operator*( matrix const&, vector const& );
	//friend matrix operator*( matrix const&, matrix const& );
	//friend void MatVecMult( matrix const&, vector const& , vector& );
	//friend void MatVecMultOverwrite( matrix const& a, vector& b);

};


#endif