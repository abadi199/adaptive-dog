#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <cmath>

class matrix;

class vector {
private:
	
public:
	int	N; 
	double* p; 
    vector() 
    {
		N = 1;
		p = new double[1];
        p[0]=1.0;
    };
	vector(int i) 
    {
		N = i;
		p = new double[N];
    };
	~vector()
	{
		delete[] p;
	}
	double& operator[](int i) { return p[i]; } 
	const double& operator[](int i) const { return p[i]; } 
	void zero() {
		for (int i = 0; i < N; i++)
			p[i] = 0.0;
	}
	void make_unit() { 
		double sum = 0.0;
		for (int i = 0; i < N; i++) {
			sum += p[i]*p[i];
		}
		sum = sqrt(sum);
		if (sum > 0.0) { 
			for (int i = 0; i < N; i++) {
				p[i] = p[i] / sum;
			}
		}
	}
	double norm() { 
		double sum = 0.0;
		for (int i = 0; i < N; i++) {
			sum += p[i]*p[i];
		}
		sum = sqrt(sum);
		return sum;
	}
	
	double get(int n) const { return p[n]; }
	int getMax() { return N; }
	void init(int i) {
		delete[] p;
		N = i;
		p = new double[N];
	}
	void print();
	void copy(vector& b);
	
	friend double operator*(const vector& a, const vector& b);
    friend vector    operator+( const vector&, const vector& );
    friend vector    operator-( vector const&, vector const& );
	friend vector    operator*( matrix const&, vector const& );
	friend void MatVecMult( matrix const&, vector const& , vector& );
	friend void MatVecMultOverwrite( matrix const& a, vector& b);
	friend void VecVecAdd( const vector& a, const vector& b, vector& c );
};


class int_vector {
private:
	int	N; // number of elements
	int* p; // a sequence of numbers

    
public:
    // constructors
    int_vector() 
    {
		N = 3;
		p = new int[3];
		// Create 3x1 unit vector
        p[0]=1; p[1]=0; p[2]=0;
    };
	int_vector(int i) 
    {
		N = i;
		p = new int[N];
		// Create 3x1 unit vector
    };
	~int_vector()
	{
		delete[] p;
	}
	////////////////////////////////////////////
	// vector[i] returns the l-value of p[i]
	// Thus, vector[i] = 1.0 is possible!
	// tmp = A[i]; is also possible!
	int& operator[](int i) { return p[i]; } 
	
	int get( int n ) const { return p[n]; }
	int getMax() { return N; }

	// negation
    //friend vector      operator-( vector const& );
    //friend unit_vector operator-( unit_vector const& );

    // addtion
    //friend vector&   operator+=( vector&, vector const& );
    
};



#endif