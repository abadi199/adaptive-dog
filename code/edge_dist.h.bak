#pragma once

#define ABS(x) ( ((x)>0) ? (x) : (-(x)) )

class edge_dist
{
public:
	double r_mean, r_std;
	double g_mean, g_std;
	double b_mean, b_std;
	double gr_mean, gr_std;
	double a; // alpha

public:
	edge_dist() {
		r_mean = r_std = 0.0;
		g_mean = g_std = 0.0;
		b_mean = b_std = 0.0;
		gr_mean = gr_std = 0.0;
		a = 0.05;
	}
	void zero() {
		r_mean = r_std = 0.0;
		g_mean = g_std = 0.0;
		b_mean = b_std = 0.0;
		gr_mean = gr_std = 0.0;
	}
	void SetAlpha(double alpha) {
		a = alpha;
	}
	void print() {
		TRACE("r_mean = %.3f, r_std = %.5f\n", r_mean, r_std);
		TRACE("g_mean = %.3f, g_std = %.5f\n", g_mean, g_std);
		TRACE("b_mean = %.3f, b_std = %.5f\n", b_mean, b_std);
		TRACE("gr_mean = %.5f, gr_std = %.5f\n", gr_mean, gr_std);
	}
	void update(edge_dist& s) {
		if (r_mean == 0.0 && g_mean == 0.0 && b_mean == 0.0
			&& r_std == 0.0 && g_std == 0.0 && b_std == 0.0) { // first data
			r_mean = s.r_mean;
			g_mean = s.g_mean;
			b_mean = s.b_mean;
			gr_mean = s.gr_mean;
			r_std = s.r_std;
			g_std = s.g_std;
			b_std = s.b_std;
			gr_std = s.gr_std;
		}
		else {
			r_mean = (1-a) * r_mean + a * s.r_mean;
			g_mean = (1-a) * g_mean + a * s.g_mean;
			b_mean = (1-a) * b_mean + a * s.b_mean;
			gr_mean = (1-a) * gr_mean + a * s.gr_mean;
			r_std = (1-a) * r_std + a * ABS(s.r_mean - r_mean);
			g_std = (1-a) * g_std + a * ABS(s.g_mean - g_mean);
			b_std = (1-a) * b_std + a * ABS(s.b_mean - b_mean);
			gr_std = (1-a) * gr_std + a * ABS(s.gr_mean - gr_mean);
		}
		
		
		//TRACE("r_mean = %.1f, r_std = %.1f\n", r_mean, r_std);
		//TRACE("g_mean = %.1f, g_std = %.1f\n", g_mean, g_std);
		//TRACE("b_mean = %.1f, b_std = %.1f\n", b_mean, b_std);
		//TRACE("gr_mean = %.5f, gr_std = %.5f\n", gr_mean, gr_std);
	}
	void update(double r, double g, double b, double gr) {
		if (r_mean == 0.0 && g_mean == 0.0 && b_mean == 0.0
			&& r_std == 0.0 && g_std == 0.0 && b_std == 0.0) {
				TRACE("first data!\n");
				r_mean = r;
				g_mean = g;
				b_mean = b;
				gr_mean = gr;
		}				
		else {
			r_mean = (1-a) * r_mean + a * r;
			g_mean = (1-a) * g_mean + a * g;
			b_mean = (1-a) * b_mean + a * b;
			gr_mean = (1-a) * gr_mean + a * gr;
			r_std = (1-a) * r_std + a * ABS(r - r_mean);
			g_std = (1-a) * g_std + a * ABS(g - g_mean);
			b_std = (1-a) * b_std + a * ABS(b - b_mean);
			gr_std = (1-a) * gr_std + a * ABS(gr - gr_mean);
		}
		//TRACE("r_mean = %.1f, r_std = %.1f\n", r_mean, r_std);
		//TRACE("g_mean = %.1f, g_std = %.1f\n", g_mean, g_std);
		//TRACE("b_mean = %.1f, b_std = %.1f\n", b_mean, b_std);
		//TRACE("gr_mean = %.5f, gr_std = %.5f\n", gr_mean, gr_std);
	}
		//~edge_dist();
};
