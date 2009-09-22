#ifndef _CIMATRIX_H_
#define _CIMATRIX_H_

//class cimatrix;
//class ivector;

class cimatrix {
private:
	int Nr, Nc;
	//double** p = new int*[rows];
    //vector p[3];
	//vector* p; // a sequence of row vectors
	//int** p; 
	iRGB** p; 
	void delete_all() {
		for (int i = 0; i < Nr; i++) {
			delete[] p[i];
		}
		delete[] p;
	}
    // addition & subtraction
    //friend matrix    operator+( matrix const&, matrix const& );
public:
    // constructors
	cimatrix() 
    {
		Nr = 1, Nc = 1;
		//p = new vector[3];
		p = new iRGB*[Nr];
		for(int i = 0; i < Nr; i++)
		   p[i] = new iRGB[Nc];
		// Create 3x3 identity matrix
        p[0][0].r = 255; 
		p[0][0].g = 255; 
		p[0][0].b = 255; 
    };
	cimatrix(int i, int j) 
    {
		Nr = i, Nc = j;
		//p = new vector[i];
		
		p = new iRGB*[Nr];
		for(i = 0; i < Nr; i++)
		   p[i] = new iRGB[Nc];
		// Create ixj identity matrix
    };
	cimatrix(cimatrix& b) {// copy constructor
		// Without this copy constructor, the actual data pointed to by 'p' are not copied!
		// Especially when the class is passed as a function parameter
		Nr = b.Nr;
		Nc = b.Nc;
		//p = new vector[i];
		//delete[] p;
		p = new iRGB*[Nr];
		for (int i = 0; i < Nr; i++) {
			p[i] = new iRGB[Nc];
			for (int j = 0; j < Nc; j++) {
				p[i][j] = b[i][j];
			}
		}
	}
	void init(int i, int j) 
    {
		delete_all();
		Nr = i, Nc = j;
		//p = new vector[i];
		//delete[] p;
		p = new iRGB*[Nr];
		for(i = 0; i < Nr; i++)
		   p[i] = new iRGB[Nc];
    };

	~cimatrix()
	{
		delete_all();
	}
	////////////////////////////////////////////
	/// matrix[i] returns the l-value of p[i] (ith row vector)
	// Thus, matrix[i][j] returns the l-value of p[i][j]
	// Thus, matrix[i][j] = 1.0 is possible!
	//vector& operator[](int i) { return p[i]; };
	///////////////////////////////////////////////
	// Very important!
	iRGB* operator[](int i) { return p[i]; };
	// This enables BOTH of the following expressions:
	// A[i][j] = 1.5;, tmp = A[i][j];
	//double*& operator[](int i) { return p[i]; }; This also works, but redundant!

	//doublegetValue( int i, int j ) const { return p[i].getValue( j );}
	iRGB& get( int i, int j ) const { return p[i][j]; }
	int getRow() const { return Nr; }
	int getCol() const { return Nc; }
	
	void zero()
	{
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) { 
				p[i][j].r = 0;
				p[i][j].g = 0;
				p[i][j].b = 0;
			}
	}
	void copy(cimatrix& b)
	{
		init(b.Nr, b.Nc);
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) {
				p[i][j] = b.p[i][j];
			}
	}
	void print() {
		for (int i = 0; i < Nr; i++) {
			for (int j = 0; j < Nc; j++) {
				TRACE("  [%d][%d] = [%d,%d,%d]", i, j, p[i][j].r, p[i][j].g, p[i][j].b);
			}
			TRACE("\n");
		}
		TRACE("\n");
	}
};


////////////////////////////////////////
// Color (double) matrix

class cmatrix {
private:
	int Nr, Nc;
	//double** p = new int*[rows];
    //vector p[3];
	//vector* p; // a sequence of row vectors
	//int** p; 
	RGB** p; 
	void delete_all() {
		for (int i = 0; i < Nr; i++) {
			delete[] p[i];
		}
		delete[] p;
	}
    // addition & subtraction
    //friend matrix    operator+( matrix const&, matrix const& );
public:
    // constructors
	cmatrix() 
    {
		Nr = 1, Nc = 1;
		//p = new vector[3];
		p = new RGB*[Nr];
		for(int i = 0; i < Nr; i++)
		   p[i] = new RGB[Nc];
		// Create 3x3 identity matrix
        p[0][0].r = 255.; 
		p[0][0].g = 255.; 
		p[0][0].b = 255.; 
    };
	cmatrix(int i, int j) 
    {
		Nr = i, Nc = j;
		//p = new vector[i];
		
		p = new RGB*[Nr];
		for(i = 0; i < Nr; i++)
		   p[i] = new RGB[Nc];
		// Create ixj identity matrix
    };
	cmatrix(cmatrix& b) {// copy constructor
		// Without this copy constructor, the actual data pointed to by 'p' are not copied!
		// Especially when the class is passed as a function parameter
		Nr = b.Nr;
		Nc = b.Nc;
		//p = new vector[i];
		//delete[] p;
		p = new RGB*[Nr];
		for (int i = 0; i < Nr; i++) {
			p[i] = new RGB[Nc];
			for (int j = 0; j < Nc; j++) {
				p[i][j] = b[i][j];
			}
		}
	}
	void init(int i, int j) 
    {
		delete_all();
		Nr = i, Nc = j;
		//p = new vector[i];
		//delete[] p;
		p = new RGB*[Nr];
		for(i = 0; i < Nr; i++)
		   p[i] = new RGB[Nc];
    };

	~cmatrix()
	{
		delete_all();
	}
	////////////////////////////////////////////
	/// matrix[i] returns the l-value of p[i] (ith row vector)
	// Thus, matrix[i][j] returns the l-value of p[i][j]
	// Thus, matrix[i][j] = 1.0 is possible!
	//vector& operator[](int i) { return p[i]; };
	///////////////////////////////////////////////
	// Very important!
	RGB* operator[](int i) { return p[i]; };
	// This enables BOTH of the following expressions:
	// A[i][j] = 1.5;, tmp = A[i][j];
	//double*& operator[](int i) { return p[i]; }; This also works, but redundant!

	//doublegetValue( int i, int j ) const { return p[i].getValue( j );}
	RGB& get( int i, int j ) const { return p[i][j]; }
	int getRow() const { return Nr; }
	int getCol() const { return Nc; }
	
	void zero()
	{
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) { 
				p[i][j].r = 0.0;
				p[i][j].g = 0.0;
				p[i][j].b = 0.0;
			}
	}
	void copy(cmatrix& b)
	{
		init(b.Nr, b.Nc);
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) {
				p[i][j] = b.p[i][j];
			}
	}
	void print() {
		for (int i = 0; i < Nr; i++) {
			for (int j = 0; j < Nc; j++) {
				TRACE("  [%d][%d] = [%f,%f,%f]", i, j, p[i][j].r, p[i][j].g, p[i][j].b);
			}
			TRACE("\n");
		}
		TRACE("\n");
	}
};

#endif