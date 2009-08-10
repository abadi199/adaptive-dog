#include "stdafx.h"
#include "Cube.h"

#include "CubeDoc.h"
#include "CubeView.h"

#include "globals.h"

#include <cmath>
#include <climits>

#define PI 3.1415926535
#define INFI 	10000

#define TWOD	0
#define WIRE	1
#define RENDER	2

typedef struct{
	float upper;
	float xvalue;
	float rslope;
	int whichpol;
} Edge;

int FOCAL = 200;

int	image_mode;	
int cal3D;
point3D camera, Prevcamera, vanVec;
float camRot = 0.0;
float camRot2 = 0.0;
float viewAngle=90;
char name[256];
extern char prefix[256];

GLuint texName[MAX_TIP_FG_NUM];
GLuint tex_name;
GLuint tex_names[4];

extern point2D vanish;
extern point2D fg[MAX_TIP_FG_NUM][MAX_FRAMES][11];
extern point2D TC[MAX_TIP_FG_NUM][MAX_FRAMES][4];
extern int isVanish, fg_num;

static GLfloat winWidth, winHeight;
extern int file_loaded;
extern GLubyte *Dbuffer;
extern GLubyte *Maskbuffer;
extern GLubyte *Backbuffer;
extern GLubyte *Membuffer;
extern GLubyte *Doublebuffer;
extern float x_size, y_size;
extern int IMAGE_X, IMAGE_Y;
extern int view_mode;
point2D ImageCoord[7];

hCoord WC[7];
hCoord FG[MAX_TIP_FG_NUM][MAX_FRAMES][MAX_TIP_FG_NUM];

point3D Wc, ray;
point2D Vc;
float half_x, half_y;

int cur_t = 0;

/*
void renderInit()
{
	int i;

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL); 
	glDisable(GL_LIGHTING);
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0 );
	
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	
	glClearDepth(1.0);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	
	glGenTextures(1+fg_num, texName);
	glBindTexture(GL_TEXTURE_2D, texName[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, TEXTURE, TEXTURE, 0, GL_RGBA,
				 GL_UNSIGNED_BYTE, BGtexture);

	/////////////////////////////////////////////////////////////////////////
	/// This is the initial texture creation, not the actual binding.
	for(i=1; i<=fg_num; i++){
		glBindTexture(GL_TEXTURE_2D, texName[i]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
 		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, TEXTURE, TEXTURE, 0, GL_RGBA,
			   			         GL_UNSIGNED_BYTE, FGtexture[i]);
	}
	glLoadIdentity();
}
*/

void renderInit()
// Alpha texturemap
{
	//int i;
	
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL); 
	glDisable(GL_LIGHTING);
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0 );
	
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	
	glClearDepth(1.0);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	///////////////////////////////////////////////////////
	// Alpha texture map
	//glGenTextures(1, texName);
	glGenTextures(1, &tex_name);
	//glBindTexture(GL_TEXTURE_2D, texName[0]);
	glBindTexture(GL_TEXTURE_2D, tex_name);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, TEXTURE, TEXTURE, 0, GL_ALPHA, GL_UNSIGNED_BYTE, al_texture2);
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, TEXTURE, TEXTURE, 0, GL_RGBA, GL_UNSIGNED_BYTE, al_texture);
	
	glGenTextures(4, tex_names);
	for (int i = 0; i < 4; i++) {
		glBindTexture(GL_TEXTURE_2D, tex_names[i]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, TEXTURE, TEXTURE, 0, GL_ALPHA, GL_UNSIGNED_BYTE, al_texture3[i]);
		
	}
	glLoadIdentity();
}

#define RGB_GETRED(rgb)    ((rgb) & 0xff) 
#define RGB_GETGREEN(rgb)    (((rgb) >> 8) & 0xff) 
#define RGB_GETBLUE(rgb)    (((rgb) >> 16) & 0xff) 

void InitAlphaMap(CDC& dc)
{	
	GLubyte r, g, b;

	brush_buffer = new GLubyte [TEXTURE * TEXTURE];
	//brush_buffer = new GLubyte [TEXTURE * TEXTURE * 4];
	
	/////////////////////////////////////////////////////////////
	if ( !LoadBMP(dc, (char *)LPCTSTR("brush_texture64_3.bmp"), &memDC3) ) { // oil
	//if ( !LoadBMP(dc, (char *)LPCTSTR("brush_texture64_water3.bmp"), &memDC3) ) { // water (gaussian)
	//if ( !LoadBMP(dc, (char *)LPCTSTR("brush_texture64_water2.bmp"), &memDC3) ) { // water
	///////////////////////////////////////////////////////////////////
	//if ( !LoadBMP(dc, (char *)LPCTSTR("flower.bmp"), &memDC3) ) {
	//if ( !LoadBMP(dc, (char *)LPCTSTR("brushRGB3.bmp"), &memDC3) ) {
	//if ( !LoadBMP(dc, (char *)LPCTSTR("brush_texture64_2.bmp"), &memDC3) ) { // oil
	//if ( !LoadBMP(dc, (char *)LPCTSTR("_black_32.bmp"), &memDC3) ) {
	//if ( !LoadBMP(dc, (char *)LPCTSTR("pen_ink_32.bmp"), &memDC3) ) {
		TRACE("brushRGB.bmp does not exist!\n");
	}
	else 
		TRACE("brushRGB.bmp exists!\n");
	
	//dc.BitBlt(0, 0, TEXTURE, TEXTURE, &memDC3, 0, 0, SRCCOPY);
	
	///*
	for (int x = 0; x < TEXTURE; x++) {
		for (int y = 0; y < TEXTURE; y++) {
			r = (GLubyte)RGB_GETRED(memDC3.GetPixel(x, TEXTURE-1-y));
			g = (GLubyte)RGB_GETGREEN(memDC3.GetPixel(x, TEXTURE-1-y));
			b = (GLubyte)RGB_GETBLUE(memDC3.GetPixel(x, TEXTURE-1-y));
			brush_buffer[y * TEXTURE + x] = r;
		}
	}
	gluScaleImage(GL_ALPHA, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, brush_buffer, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, al_texture2);
	//gluScaleImage(GL_RGBA, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, brush_buffer, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, al_texture);
	
	if ( !LoadBMP(dc, (char *)LPCTSTR("_stroke_32_1.bmp"), &memDC3) ) 
		TRACE("brushRGB.bmp does not exist!\n");
	
	///*
	for (int x = 0; x < TEXTURE; x++) {
		for (int y = 0; y < TEXTURE; y++) {
			r = (GLubyte)RGB_GETRED(memDC3.GetPixel(x, TEXTURE-1-y));
			g = (GLubyte)RGB_GETGREEN(memDC3.GetPixel(x, TEXTURE-1-y));
			b = (GLubyte)RGB_GETBLUE(memDC3.GetPixel(x, TEXTURE-1-y));
			brush_buffer[y * TEXTURE + x] = r;
		}
	}
	gluScaleImage(GL_ALPHA, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, brush_buffer, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, 
		al_texture3[0]);
	//gluScaleImage(GL_RGBA, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, brush_buffer, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, al_texture);

	/*
	if ( !LoadBMP(dc, (char *)LPCTSTR("_stroke_32_2.bmp"), &memDC3) ) 
		TRACE("brushRGB.bmp does not exist!\n");
	
	///*
	for (int x = 0; x < TEXTURE; x++) {
		for (int y = 0; y < TEXTURE; y++) {
			r = (GLubyte)RGB_GETRED(memDC3.GetPixel(x, TEXTURE-1-y));
			g = (GLubyte)RGB_GETGREEN(memDC3.GetPixel(x, TEXTURE-1-y));
			b = (GLubyte)RGB_GETBLUE(memDC3.GetPixel(x, TEXTURE-1-y));
			brush_buffer[y * TEXTURE + x] = r;
		}
	}
	gluScaleImage(GL_ALPHA, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, brush_buffer, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, 
		al_texture3[1]);
	//gluScaleImage(GL_RGBA, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, brush_buffer, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, al_texture);

	if ( !LoadBMP(dc, (char *)LPCTSTR("_stroke_32_3.bmp"), &memDC3) ) 
		TRACE("brushRGB.bmp does not exist!\n");
	
	///*
	for (int x = 0; x < TEXTURE; x++) {
		for (int y = 0; y < TEXTURE; y++) {
			r = (GLubyte)RGB_GETRED(memDC3.GetPixel(x, TEXTURE-1-y));
			g = (GLubyte)RGB_GETGREEN(memDC3.GetPixel(x, TEXTURE-1-y));
			b = (GLubyte)RGB_GETBLUE(memDC3.GetPixel(x, TEXTURE-1-y));
			brush_buffer[y * TEXTURE + x] = r;
		}
	}
	gluScaleImage(GL_ALPHA, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, brush_buffer, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, 
		al_texture3[2]);
	
	if ( !LoadBMP(dc, (char *)LPCTSTR("_stroke_32_4.bmp"), &memDC3) ) 
		TRACE("brushRGB.bmp does not exist!\n");
	
	///*
	for (int x = 0; x < TEXTURE; x++) {
		for (int y = 0; y < TEXTURE; y++) {
			r = (GLubyte)RGB_GETRED(memDC3.GetPixel(x, TEXTURE-1-y));
			g = (GLubyte)RGB_GETGREEN(memDC3.GetPixel(x, TEXTURE-1-y));
			b = (GLubyte)RGB_GETBLUE(memDC3.GetPixel(x, TEXTURE-1-y));
			brush_buffer[y * TEXTURE + x] = r;
		}
	}
	gluScaleImage(GL_ALPHA, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, brush_buffer, TEXTURE, TEXTURE, GL_UNSIGNED_BYTE, 
		al_texture3[3]);
	*/
}



void renderReshape(int width, int height)
{
	float aspect;

	winWidth = (float)width;
	winHeight = (float)height;
	glViewport(0, 0, (int)winWidth, (int)winHeight);
	aspect = ((float) winWidth) / winHeight;

	glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	gluPerspective(100.0, aspect, 0.0, 30010.0);

	glMatrixMode(GL_MODELVIEW);
}

void renderNone(void)
{
  	static char message[] = "No Image loaded.";
  	static int width = 0;

	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glDisable(GL_LIGHTING);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, winWidth, 0, winHeight);
	glColor3f(1.0, 1.0, 1.0);
	glRasterPos2i(0, 0);
	glCallLists(sizeof(message), GL_UNSIGNED_BYTE, message);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

point3D normalize(point3D vec)
{
    point3D result;
    float size;

    size=(float)sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);

    result.x = vec.x/size;
    result.y = vec.y/size;
    result.z = vec.z/size;

    return result;
}

void CalImageCoord()
{
	half_x = (float)(IMAGE_X)/2;
	half_y = (float)(IMAGE_Y)/2;
	ImageCoord[1].x = -half_x;
	ImageCoord[1].y = half_y;
	ImageCoord[2].x = -half_x;
	ImageCoord[2].y = -half_y;
	ImageCoord[3].x = half_x;
	ImageCoord[3].y = half_y;
	ImageCoord[4].x = half_x;
	ImageCoord[4].y = -half_y;
	ImageCoord[5].x = -half_x;
	ImageCoord[5].y = vanish.y-half_y;
	ImageCoord[6].x = half_x;
	ImageCoord[6].y = vanish.y-half_y;
	
	cal3D = 0;
}

void CalWorldCoord()
{
	int i, j, k;
	float temp, rat;

	CalImageCoord();

	viewAngle = (float)((float)2*atan((float)(IMAGE_Y)/(2*FOCAL))*180/PI);
	
	WC[1].x = ImageCoord[1].x;
	WC[1].y = ImageCoord[1].y;
	WC[1].z = -(float)FOCAL;
	WC[1].w = 0.01f;
	
	WC[2].x = ImageCoord[2].x;
	WC[2].y = ImageCoord[2].y;
	WC[2].z = -(float)FOCAL;
	WC[2].w = 1;
	
	WC[3].x = ImageCoord[3].x;
	WC[3].y = ImageCoord[3].y;
	WC[3].z = -(float)FOCAL;
	WC[3].w = 0.01f;
	
	WC[4].x = ImageCoord[4].x;
	WC[4].y = ImageCoord[4].y;
	WC[4].z = -(float)FOCAL;
	WC[4].w = 1;
	
	WC[5].x = ImageCoord[5].x;
	WC[5].y = ImageCoord[5].y;
	WC[5].z = -(float)FOCAL;
	WC[5].w = 0.01f;
	
	WC[6].x = ImageCoord[6].x;
	WC[6].y = ImageCoord[6].y;
	WC[6].z = -(float)FOCAL;
	WC[6].w = 0.01f;

	vanVec.x = 0;
	vanVec.y = vanish.y-half_y;
	vanVec.z = -(float)FOCAL;
	vanVec = normalize(vanVec);

	if (fg_num) {
		temp = ( ImageCoord[1].x-ImageCoord[3].x ) * FOCAL * ( ImageCoord[6].y-ImageCoord[4].y );

		for ( i = 1; i <= fg_num; i++ ) {
			for ( k = 0; k < 1; k++ ) {

				    TRACE ("[%d] I'm HERE in YZflas == 0 !!!\n", i);
					FG[i][k][1].x = temp * ( fg[i][k][1].x-half_x );
					FG[i][k][1].y = temp * ( fg[i][k][1].y-half_y );
					FG[i][k][1].z = -FOCAL * temp; // -FOCAL means fg[i][k][1].z
					FG[i][k][1].w = temp * ( ImageCoord[6].y - (fg[i][k][1].y - half_y) )
							/ ( ImageCoord[6].y - ImageCoord[4].y );

					FG[i][k][2].x = temp * ( fg[i][k][2].x-half_x );
					FG[i][k][2].y = temp * ( fg[i][k][2].y-half_y );
					FG[i][k][2].z = -FOCAL * temp;
					FG[i][k][2].w = temp * ( ImageCoord[6].y - fg[i][k][2].y+half_y )
							/ ( ImageCoord[6].y - ImageCoord[4].y );

					for (j = 3; j <= fg[i][k][0].x; j++) { // fg[i][k][0].x contains the number of vertices
					// the other vertices
						if ( FG[i][k][1].z / FG[i][k][1].w < FG[i][k][2].z / FG[i][k][2].w ) {
							FG[i][k][j].z = ( (fg[i][k][j].x-fg[i][k][1].x) / (fg[i][k][2].x-fg[i][k][1].x) )
									* (FG[i][k][2].z-FG[i][k][1].z) + FG[i][k][1].z;
							rat = -FG[i][k][j].z / FOCAL;
							FG[i][k][j].x = ( fg[i][k][j].x-half_x ) * rat;
							FG[i][k][j].y = ( fg[i][k][j].y-half_y ) * rat;
							FG[i][k][j].w = ( (fg[i][k][j].x-fg[i][k][1].x) / (fg[i][k][2].x-fg[i][k][1].x) )
									* (FG[i][k][2].w-FG[i][k][1].w) + FG[i][k][1].w;
						}
						else {
						// if the 1st vertex lies farther than the 2nd vertex
							FG[i][k][j].z = ( (fg[i][k][2].x - fg[i][k][j].x) / (fg[i][k][2].x - fg[i][k][1].x) )
									* (FG[i][k][1].z- FG[i][k][2].z) + FG[i][k][2].z;
							rat = -FG[i][k][j].z/FOCAL;
							FG[i][k][j].x = ( fg[i][k][j].x-half_x ) * rat;
							FG[i][k][j].y = ( fg[i][k][j].y-half_y ) * rat;
							FG[i][k][j].w = ( ( fg[i][k][2].x - fg[i][k][j].x ) / ( fg[i][k][2].x - fg[i][k][1].x ) )
									* ( FG[i][k][1].w - FG[i][k][2].w ) + FG[i][k][2].w;
						}
					}
				
			} // end of for k
		} // end of for i
	} // end of if fg_num
	cal3D = 1; // The coordinate calculation has been already done once
}

void TIPrender(CDC* pDC)
{
	int i,j;
	GLint viewport[4];
	static int first_count = 0;
	
	if ( status==TWOD || status == NO_STATUS ) {
    	glGetIntegerv(GL_VIEWPORT, viewport);
   		glMatrixMode(GL_PROJECTION);
   		glLoadIdentity();

   		gluOrtho2D(0, viewport[2], 0, viewport[3]);

		/// this makes the image and the lines correctly appears
		glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

		/// Initializes the camera parameters
		camRot = 0.0;
		camRot2 = 0.0;
		camera.x = 0.0;
		camera.y = 0.0;
		camera.z = 0.0;

		POINT point = {0, 0};	// initial offset = (0, 0)
		
		if(status == TWOD || status == NO_STATUS) {  

			if (image_mode == OPENGL_MODE) {
				// For planar image !!!
				glRasterPos2i(0, 0);
				glDrawPixels(IMAGE_X, IMAGE_Y, GL_RGB, GL_UNSIGNED_BYTE, Dbuffer);
				//glClearColor(1.0, 1.0, 0.0, 1.0);
				//glClear(GL_COLOR_BUFFER_BIT);
				swap();
			}
			else if (image_mode == BITMAP_MODE) {
				// For panoramic image !!!
				pDC->BitBlt(0, 0, IMAGE_X, IMAGE_Y, &memDC, 0, 0, SRCCOPY);
				//
			}
						
				
		}
			
	}

	// draw here !! 
	if (status == TWOD) {
		glPointSize(3.0);
		glLineWidth(3.0);
		glColor4f(0.0, 1.0, 1.0, 1.0);
		if (isVanish) {
			// Drawing the vanishing line
			CalImageCoord();
			glBegin(GL_LINES);
			glVertex2f(ImageCoord[5].x+half_x, ImageCoord[5].y+half_y);
			glVertex2f(ImageCoord[6].x+half_x, ImageCoord[6].y+half_y);
			glEnd();
		}
		if (fg_num) {	// Drawing rectangles surrounding foreground objects
			for(i=1; i<=fg_num; i++){
				glBegin(GL_LINE_STRIP);
				glVertex2f(TC[i][cur_t][0].x, TC[i][cur_t][0].y);
				glVertex2f(TC[i][cur_t][1].x, TC[i][cur_t][1].y);
				glVertex2f(TC[i][cur_t][2].x, TC[i][cur_t][2].y);
				glVertex2f(TC[i][cur_t][3].x, TC[i][cur_t][3].y);
				glVertex2f(TC[i][cur_t][0].x, TC[i][cur_t][0].y);
				glEnd();
			}
			glColor4f(1.0, 1.0, 0.0, 1.0);
			for(i=1; i<=fg_num; i++) { 
				// Drawing polygons surrounding exact portion of foreground objects
				glBegin(GL_LINE_LOOP);
				for(j=1; j<=(int)fg[i][cur_t][0].x; j++) { // fg[i][0].x is the total vertex no
					glVertex2f(fg[i][cur_t][j].x, fg[i][cur_t][j].y);
				}
				glVertex2f(fg[i][cur_t][1].x, fg[i][cur_t][1].y);
				glEnd();
			}
		}
	}
	if (status == RENDER) {
		if (!cal3D) {	
			CalWorldCoord(); 
			renderInit(); 
		}
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    	glGetIntegerv(GL_VIEWPORT, viewport);
   		glMatrixMode(GL_PROJECTION);
   		glLoadIdentity();

		viewAngle = 68;
		gluPerspective(viewAngle, (float)((float)IMAGE_X/(float)(IMAGE_Y)), FOCAL, 30010.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(0, 0, 0, 0, 0, -1.0, 0.0, 1.0, 0.0);
		glTranslatef(-camera.x, -camera.y, -camera.z);
		glRotatef(-camRot, 0.0, 1.0, 0.0);
		glRotatef(-camRot2, 1.0, 0.0, 0.0);

		TRACE("camera.x = %f, camera.y = %f, camera.z = %f\n", camera.x, camera.y, camera.z);
		TRACE("camRot = %f\n", camRot);
		TRACE("camRot2 = %f\n", camRot2);
		TRACE("viewAngle = %f\n", viewAngle);

		glDisable(GL_LIGHTING);
		glEnable(GL_TEXTURE_2D);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		glBindTexture(GL_TEXTURE_2D, texName[0]);	// the background image texture
		glColor3f(1.0, 1.0, 1.0);
		
		///////////////////////////////////////////////////////////////
		/// Background model vertex coordinates for TIP
		glBegin(GL_QUAD_STRIP);
			glTexCoord2f(0.0, 1.0);
			glVertex4f(WC[1].x, WC[1].y, WC[1].z, WC[1].w);
			glTexCoord2f(1.0, 1.0);
			glVertex4f(WC[3].x, WC[3].y, WC[3].z, WC[3].w);

			glTexCoord2f(0.0, vanish.y / (float)(IMAGE_Y));
			glVertex4f(WC[5].x, WC[5].y, WC[5].z, WC[5].w);
			glTexCoord2f(1.0, vanish.y / (float)(IMAGE_Y));
			glVertex4f(WC[6].x, WC[6].y, WC[6].z, WC[6].w);
			glTexCoord2f(0.0, 0.0);
			glVertex4f(WC[2].x, WC[2].y, WC[2].z, WC[2].w);
			glTexCoord2f(1.0, 0.0);
			glVertex4f(WC[4].x, WC[4].y, WC[4].z, WC[4].w);
		glEnd();
			
		// These are for the foreground polygons
		for(i=1; i<=fg_num; i++){
			glBindTexture(GL_TEXTURE_2D, texName[i]); 
	
			glColor4f(1.0, 1.0, 1.0, 1.0);
			glBegin(GL_POLYGON);
			for(j=1; j<=(int)fg[i][cur_t][0].x; j++){
				glTexCoord4f( 
					(fg[i][cur_t][j].x-TC[i][cur_t][0].x) / (TC[i][cur_t][1].x-TC[i][cur_t][0].x) / FG[i][cur_t][j].w,
					(fg[i][cur_t][j].y-TC[i][cur_t][0].y) / (TC[i][cur_t][3].y-TC[i][cur_t][0].y) / FG[i][cur_t][j].w, 
							0, 1 / FG[i][cur_t][j].w );
				
				glVertex3f(FG[i][cur_t][j].x / FG[i][cur_t][j].w, 
					FG[i][cur_t][j].y / FG[i][cur_t][j].w, 
					FG[i][cur_t][j].z / FG[i][cur_t][j].w);
				
			}
			glEnd();
		}

		glDisable(GL_TEXTURE_2D);

	}
	if (status == PANORAMIC_NO_STATUS || status == PANORAMIC_SPECIFICATION || status == PANORAMIC_RENDER) {
		//TIPPanoramicDrawScene(pDC);
	}
	
	
}

void TIPDrawScene(CDC *pDC)
{
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

	if (!file_loaded) {
		renderNone(); 
		swap();
	}
	else
		TIPrender(pDC);

	swap(); 
}


void renderImage()
{
	glRasterPos2i(0, 0);
	glDrawPixels(IMAGE_X, IMAGE_Y, GL_RGB, GL_UNSIGNED_BYTE, Dbuffer);
}

