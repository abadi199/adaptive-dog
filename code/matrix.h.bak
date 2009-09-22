#ifndef _MATRIX_H_
#define _MATRIX_H_

class vector;

class matrix {
private:
	int Nr, Nc;
	//double** p = new int*[rows];
    //vector p[3];
	//vector* p; // a sequence of row vectors
	double** p; 
	void delete_all() {
		for (int i = 0; i < Nr; i++) 
				delete[] p[i];
			delete[] p;
	}
    // addition & subtraction
    //friend matrix    operator+( matrix const&, matrix const& );
public:
    // constructors
	matrix() 
    {
		Nr = 3, Nc = 3;
		//p = new vector[3];
		p = new double*[Nr];
		for(int i = 0; i < Nr; i++)
		   p[i] = new double[Nc];
		// Create 3x3 identity matrix
        p[0][0]=1.0; p[0][1]=0.0; p[0][2]=0.0;
        p[1][0]=0.0; p[1][1]=1.0; p[1][2]=0.0;
        p[2][0]=0.0; p[2][1]=0.0; p[2][2]=1.0;
    };
	matrix(int i, int j) 
    {
		Nr = i, Nc = j;
		//p = new vector[i];
		p = new double*[Nr];
		for(i = 0; i < Nr; i++)
		   p[i] = new double[Nc];
		// Create ixj identity matrix
        //p[0][0]=1.0; p[0][1]=0.0; p[0][2]=0.0;
        //p[1][0]=0.0; p[1][1]=1.0; p[1][2]=0.0;
        //p[2][0]=0.0; p[2][1]=0.0; p[2][2]=1.0;
    };
	matrix(matrix& b) {// copy constructor
		Nr = b.Nr;
		Nc = b.Nc;
		//p = new vector[i];
		//delete[] p;
		p = new double*[Nr];
		for (int i = 0; i < Nr; i++) {
			p[i] = new double[Nc];
			for (int j = 0; j < Nc; j++) {
				p[i][j] = b[i][j];
			}
		}
	}
	~matrix() {
		delete_all();
	}
	////////////////////////////////////////////
	/// matrix[i] returns the l-value of p[i] (ith row vector)
	// Thus, matrix[i][j] returns the l-value of p[i][j]
	// Thus, matrix[i][j] = 1.0 is possible!
	//vector& operator[](int i) { return p[i]; };
	///////////////////////////////////////////////
	// Very important!
	double* operator[](int i) { return p[i]; };
	// This enables BOTH of the following expressions:
	// A[i][j] = 1.5;, tmp = A[i][j];
	//double*& operator[](int i) { return p[i]; }; This also works, but redundant!

	//doublegetValue( int i, int j ) const { return p[i].getValue( j );}
	double& get( int i, int j ) const { return p[i][j]; }
	int getRow() const { return Nr; }
	int getCol() const { return Nc; }
	void init(int i, int j) 
    {
		delete_all();
		Nr = i, Nc = j;
		//p = new vector[i];
		//delete[] p;
		p = new double*[Nr];
		for(i = 0; i < Nr; i++)
		   p[i] = new double[Nc];
    };
	void zero()
	{
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) 
				p[i][j] = 0.0;
	}
	void set_all(double x)
	{
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) 
				p[i][j] = x;
	}

	friend void ludcmp(matrix& a, int_vector& indx, double &d);
	friend void lubksb(matrix& a, int_vector& indx, vector &b);
	friend void inverse(matrix &a, matrix &inverse);
	friend void inverse2(matrix &a, matrix &inverse); // does not destroy original matrix a!!!
	friend void transpose(matrix &a, matrix &transpose);
	friend void merge_h(matrix &a, matrix &b, matrix &merged);
	friend void split_v(matrix &merged, int a_row, matrix &a, matrix &b);
	void print();
	matrix inverse(); // This is the better version
	matrix copy();
	void copy(matrix& b);
	
	friend vector operator*( matrix const&, vector const& );
	friend matrix operator*( matrix const&, matrix const& );
	friend void MatVecMult( matrix const&, vector const& , vector& );
	friend void MatVecMultOverwrite( matrix const& a, vector& b);
	friend void MatMatMult( matrix const& a, matrix const& b, matrix& c);
	friend void MatMatAdd( matrix const& a, matrix const& b, matrix& c);

};


#endif