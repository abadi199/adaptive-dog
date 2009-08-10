#ifndef _QMATRIX_H_
#define _QMATRIX_H_

// Matrix of queues
// Each queue contains a list of type T

template <class T> class qmatrix {
private:
	int Nr, Nc;
	deque<T>** p; // each element has a deque (not a pointer to deque)
	/////////////////////////////////////
	void delete_all() {
		for (int i = 0; i < Nr; i++) 
			delete[] p[i];
		delete[] p;
	}
	////////////////////////////////////
    // addition & subtraction
    //friend matrix    operator+( matrix const&, matrix const& );
public:
    // constructors
	qmatrix() 
    {
		Nr = 3, Nc = 3;
		//p = new vector[3];
		p = new deque<T>*[Nr];
		for(int i = 0; i < Nr; i++)
		   p[i] = new deque<T>[Nc];
		// Create 3x3 identity matrix
        //p[0][0]=1; p[0][1]=0; p[0][2]=0;
        //p[1][0]=0; p[1][1]=1; p[1][2]=0;
        //p[2][0]=0; p[2][1]=0; p[2][2]=1;
    };
	qmatrix(int i, int j) 
    {
		Nr = i, Nc = j;
		//p = new vector[i];
		
		p = new deque<T>*[Nr];
		for(i = 0; i < Nr; i++)
		   p[i] = new deque<T>[Nc];
		// Create ixj identity matrix
    };
	qmatrix(qmatrix& b) {// copy constructor
		// Without this copy constructor, the actual data pointed to by 'p' are not copied!
		// Especially when the class is passed as a function parameter
		Nr = b.Nr;
		Nc = b.Nc;
		//p = new vector[i];
		//delete[] p;
		p = new deque<T>*[Nr];
		for (int i = 0; i < Nr; i++) {
			p[i] = new deque<T>[Nc];
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
		p = new deque<T>*[Nr];
		for(i = 0; i < Nr; i++)
		   p[i] = new deque<T>[Nc];
    };
	~qmatrix()
	{
		delete_all();
	}
	/////////////////////////////////////////////////////////
	/// matrix[i] returns the l-value of p[i] (ith row vector)
	// Thus, matrix[i][j] returns the l-value of p[i][j]
	// Thus, matrix[i][j] = 1.0 is possible!
	//vector& operator[](int i) { return p[i]; };
	///////////////////////////////////////////////
	// Very important!
	deque<T>* operator[](int i) { return p[i]; }
	// This enables BOTH of the following expressions:
	// A[i][j] = 1.5; tmp = A[i][j];
	//double*& operator[](int i) { return p[i]; } // This also works, but redundant!

	//doublegetValue( int i, int j ) const { return p[i].getValue( j );}
	deque<T>& get( int i, int j ) const { return p[i][j]; }
	int getRow() const { return Nr; }
	int getCol() const { return Nc; }
	/////////////////////////////////
	void remove(int i, int j, T obj_ptr) { // remove a deque element containing obj_ptr at p[i][j]
		//deque<T>::iterator it;
		int k;
		if ( i > Nr-1 || i < 0 || j > Nc-1 || j < 0) {
			TRACE("qmatrix index error!\n");
			exit(1);
		}
		//TRACE("B_MAP REMOVE: [%d][%d]\n", i, j);
		//TRACE("b_map[%d][%d].size() = %d\n", i, j, p[i][j].size());

		for (k = 0; k < (signed)p[i][j].size(); k++) { // remove all instances! 
			if (p[i][j][k] == obj_ptr) {
				p[i][j].erase(p[i][j].begin() + k); 
				k--;
			}
		}
		/*
		// This is wrong! After erase(it), "it" is gone!
		it = p[i][j].begin();
		while (it != p[i][j].end()) {
			if (*it == obj_ptr) {
				if ((it+1) == p[i][j].end()) { // this was the last item
					p[i][j].erase(it); // remove all instances! 
					break;
				}
				else p[i][j].erase(it); // remove all instances! 
			}
			it++;
		}
		*/
	}
	void clear()
	{
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) 
				//p[i][j] = 0.0;
				p[i][j].clear();
	}
	void copy(imatrix& b)
	{
		init(b.Nr, b.Nc);
		for (int i = 0; i < Nr; i++) 
			for (int j = 0; j < Nc; j++) 
				p[i][j] = b.p[i][j];
	}
	void print() {
		for (int i = 0; i < Nr; i++) {
			for (int j = 0; j < Nc; j++) {
				//TRACE("  [%d][%d] = %d", i, j, p[i][j]);
			}
			TRACE("\n");
		}
		TRACE("\n");
	}
};


#endif