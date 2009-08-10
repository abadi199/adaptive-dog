#ifndef _IPEN_H_
#define _IPEN_H_

#include "MRBspline.h"

class IPen {
public:
	MRBspline curve[20];
	int	cc;

	IPen() {}
};

#endif

