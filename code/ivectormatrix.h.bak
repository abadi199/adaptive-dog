#ifndef _IVECTORMATRIX_H_
#define _IVECTORMATRIX_H_

class imatrix;

class ivector {
private:
	
public:
	int	N; // number of elements
	int* p; // a sequence of numbers
    ivector() 
    {
		N = 3;
		p = new int[N];
        p[0]=1; p[1]=0; p[2]=0;
    };
	ivector(int i) 
    {
		N = i;
		p = new int[N];
    };
	ivector(ivector& b) {
		N = b.N;
		p = new int[N];
		for (int i = 0; i < N; i++) {
			p[i] = b[i];
		}
	}
	~ivector()
	{
		delete[] p;
	}
	int& operator[](int i) { return p[i]; } 
	const int& operator[](int i) const { return p[i]; } 
	
	int get(int n) const { return p[n]; }
	int getMax() { return N; }
	void init(int i) {
		N = i;
		delete[] p;
		p = new int[N];
	}
	void print() {
		for (int i = 0; i < N; i++) {
			TRACE("[%d] = %d  ", i, p[i]);
		}
		TRACE("\n");
	}
	void copy(ivector& b) {
		delete[] p;
		N = b.N;
		p = new int[N];
		for (int i = 0; i < N; i++) {
			p[i] = b[i];
		}
	}
};

class imatrix {
private:
	int Nr, Nc;
	int** p; 
	void delete_all() {
		for (int i = 0; i < Nr; i++) 
				delete[] p[i];
			delete[] p;
	}
public:
	imatrix() 
    {
		Nr = 1, Nc = 1;
		p = new int*[Nr];
		for(int i = 0; i < Nr; i++)
		   p[i] = new int[Nc];
        p[0][0]=1; 
    };
	imatrix(int i, int j) 
    {
		Nr = i, Nc = j;
		
		p = new int*[Nr];
		for(i = 0; i < Nr; i++)
		   p[i] = new int[Nc];
    };
	imatrix(imatrix& b) {
		Nr = b.Nr;
		Nc = b.Nc;
		p = new int*[Nr];
		for (int i = 0; i < Nr; i++) {
			p[i] = new int[Nc];
			for (int j = 0; j < Nc; j++) {
				p[i][j] = b[i][j];
			}
		}
	}
	void init(int i, int j) 
    {
		delete_all();
		Nr = i, Nc = j;
		p = new int*[Nr];
		for(i = 0; i < Nr; i++)
		   p[i] = new int[Nc];
    };

	~imatrix()
	{
		delete_all();
	}
	int* operator[](int i) { return p[i]; };

	int& get( int i, int j ) const { return p[i][j]; }
	int getRow() const { return Nr; }
	int getCol() const { return Nc; }
	
	void zero()
	{
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) 
				p[i][j] = 0;
	}
	void white()
	{
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) 
				p[i][j] = 255;
	}
	void copy(imatrix& b)
	{
		init(b.Nr, b.Nc);
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) 
				p[i][j] = b.p[i][j];
	}
};


#endif