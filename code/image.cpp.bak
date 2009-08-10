#include <cstdio>
#include "stdafx.h"
#include "Cube.h"

#include "CubeDoc.h"
#include "CubeView.h"

#include "gl\gl.h"
#include "gl\glu.h"

#include "globals.h"

int IMAGE_X = 640, IMAGE_Y = 480;

GLubyte FGtexture[MAX_TIP_FG_NUM][TEXTURE][TEXTURE][4];
GLubyte BGtexture[BACK_TEXTURE_WIDTH][TEXTURE][4];

//////////////////////////////////////////////
// Alpha texture map
GLubyte al_texture2[TEXTURE][TEXTURE];
GLubyte al_texture[TEXTURE][TEXTURE][4];
GLubyte al_texture3[4][TEXTURE][TEXTURE];

GLubyte *Dbuffer;
GLubyte *Maskbuffer;
GLubyte *Backbuffer;
GLubyte *brush_buffer;
GLubyte *Membuffer;
GLubyte *Doublebuffer;

char fgText[MAX_TIP_FG_NUM][255];
char fg_maskText[MAX_TIP_FG_NUM][255];
point2D fgSize[MAX_TIP_FG_NUM];
GLubyte *fgImage[MAX_TIP_FG_NUM];
extern int fg_num;

#define MIN(A, B) (((A) < (B))? (A) : (B))
#define MAX(A, B) (((A) > (B))? (A) : (B))

// r,g,b values are from 0 to 1
// h = [0,360], s = [0,1], v = [0,1]
//		if s == 0, then h = -1 (undefined)
void RGBtoHSV( float r, float g, float b, float *h, float *s, float *v )
{
	float min, max, delta;
	min = MIN( MIN(r, g), b );
	max = MAX( MAX(r, g), b );
	*v = max;				// v
	delta = max - min;
	if( max != 0 )
		*s = delta / max;		// s
	else {
		// r = g = b = 0		// s = 0, v is undefined
		*s = 0;
		*h = -1;
		return;
	}
	if( r == max )
		*h = ( g - b ) / delta;		// between yellow & magenta
	else if( g == max )
		*h = 2 + ( b - r ) / delta;	// between cyan & yellow
	else
		*h = 4 + ( r - g ) / delta;	// between magenta & cyan
	*h *= 60;				// degrees
	if( *h < 0 )
		*h += 360;
}

void LoadFGImage()
{
	int i,x,y, index;
	unsigned char r, g, b, a;
	float h, s, v;
	FILE *fp, *mkfp;

	for(i=1; i<=fg_num; i++){
		fgImage[i] = (GLubyte *)malloc(4*(int)fgSize[i].x*(int)fgSize[i].y*sizeof(GLubyte));
		fp = fopen(fgText[i], "rb");
		ASSERT(fp != NULL);	// ASSERT: if the file does not exist, program STOPS
		// Check if the corresponding mask image exists
		mkfp = fopen(fg_maskText[i], "rb");
		//ASSERT(mkfp != NULL);

		for (y = (int)fgSize[i].y - 1; y >= 0; y--) {
			for (x = 0; x < fgSize[i].x; x++) {
				index = (y*(int)fgSize[i].x+x)*4;
				r = getc(fp);
				g = getc(fp);
				b = getc(fp);
				if (mkfp != NULL) {
					a = getc(mkfp);
				}
				fgImage[i][index] = r;
				fgImage[i][index+1] = g;
				fgImage[i][index+2] = b;

				RGBtoHSV((float)r / 255, (float)g / 255, (float)b / 255,
						 &h, &s, &v);
				if (mkfp == NULL) // Maks file does not exist
					if (230 <= h && h <= 250 && s > 0.50 && v > 0.50) {
						// If the color is almost BLUE
						// Invisible: alpha = 0;
						fgImage[i][index+3] = 0;
					}
					else{
						// Visible: alpha = 1;
						fgImage[i][index+3] = 255;
					}
				else // MK_RAW file exists !!!
					fgImage[i][index+3] = a;

			}
		}

		fclose(fp);
		if (mkfp) fclose(mkfp);
	}
}

void makeTexture()
{
	int i, x, y;
	int index;
	unsigned char r, g, b;
	float h, s, v;
	int dx[] = { -1, 0, 1, -1, 0, 1, -1, 0, 1 };
	int dy[] = { -1, -1, -1, 0, 0, 0, 1, 1, 1 };


	gluScaleImage(GL_RGBA, IMAGE_X, IMAGE_Y, GL_UNSIGNED_BYTE, Backbuffer,
				  TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, BGtexture);

	for(i=1; i<=fg_num; i++){

		for (y = 0; y < fgSize[i].y; y++) {
			for (x = 0; x < fgSize[i].x; x++) {
				index = (y*(int)fgSize[i].x+x)*4;
				r = fgImage[i][index];
				g = fgImage[i][index+1];
				b = fgImage[i][index+2];

				RGBtoHSV((float)r / 255, (float)g / 255, (float)b / 255,
						 &h, &s, &v);
				if (180 <= h && h <= 280 && s > 0.50 && v > 0.50 &&
					fgImage[i][index + 3] == 0) {
					for (int dr = 0; dr < 8; dr++) {
						if (y + dy[dr] >= 0 && y + dy[dr] < fgSize[i].y &&
							x + dx[dr] >= 0 && x + dx[dr] < fgSize[i].x &&
							fgImage[i][((y+dy[dr])*(int)fgSize[i].x+x+dx[dr])*4 + 3] != 0) {
							fgImage[i][index] = fgImage[i][((y+dy[dr])*(int)fgSize[i].x+x+dx[dr])*4];
							fgImage[i][index+1] = fgImage[i][((y+dy[dr])*(int)fgSize[i].x+x+dx[dr])*4+1];
							fgImage[i][index+2] = fgImage[i][((y+dy[dr])*(int)fgSize[i].x+x+dx[dr])*4+2];
						}
					}
				}
			}
		}

		gluScaleImage(GL_RGBA, (int)fgSize[i].x, (int)fgSize[i].y, GL_UNSIGNED_BYTE, 
			fgImage[i], TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, FGtexture[i]);

	}
}

void OpenImage(char *filename, char *maskname, char *backname)
{
    int x, y;
	unsigned char r, g, b;
	FILE *fp;

	Dbuffer = new GLubyte [IMAGE_Y * IMAGE_X * 3];
	Maskbuffer = new GLubyte [IMAGE_Y * IMAGE_X * 3];
	Backbuffer = new GLubyte [IMAGE_Y * IMAGE_X * 4];
	
    fp = fopen(filename, "rb");
	ASSERT(fp != NULL);

	for (y = IMAGE_Y - 1; y >= 0; y--) {
		for (x = 0; x < IMAGE_X; x++) {
			r = getc(fp);
			g = getc(fp);
			b = getc(fp);
			Dbuffer[(y * IMAGE_X + x) * 3 + 0] = r;
			Dbuffer[(y * IMAGE_X + x) * 3 + 1] = g;
			Dbuffer[(y * IMAGE_X + x) * 3 + 2] = b;
		}
	}

	fclose(fp);

	fp = fopen(backname, "rb");
	ASSERT(fp != NULL);

	// Reverse the vertical order for the OpenGL window
	for (y = IMAGE_Y - 1; y >= 0; y--) {
		for (x = 0; x < IMAGE_X; x++) {
			r = getc(fp);
			g = getc(fp);
			b = getc(fp);
			Backbuffer[(y * IMAGE_X + x) * 4 + 0] = r;
			Backbuffer[(y * IMAGE_X + x) * 4 + 1] = g;
			Backbuffer[(y * IMAGE_X + x) * 4 + 2] = b;
			Backbuffer[(y * IMAGE_X + x) * 4 + 3] = 255;
		}
	}

	fclose(fp);

	LoadFGImage();
	makeTexture();
}


