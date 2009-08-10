//#define CHILDFRM_TITLE_BAR_WIDTH 31 // windows 2000 style
//#define CHILDFRM_TITLE_BAR_WIDTH 38 // Windows XP style (MDI)
#define CHILDFRM_TITLE_BAR_WIDTH 88 // SDI
#define CHILDFRM_RIGHT_BAR_WIDTH 12

#define RGB_GETRED(rgb)    ((rgb) & 0xff) 
#define RGB_GETGREEN(rgb)    (((rgb) >> 8) & 0xff) 
#define RGB_GETBLUE(rgb)    (((rgb) >> 16) & 0xff)  

#define dist2(x1, y1, x2, y2) sqrt( (((double)x1)-((double)x2))*(((double)x1)-((double)x2)) + (((double)y1)-((double)y2))*(((double)y1)-((double)y2)) )
#define dist3(x, y, z) sqrt(((double)x)*((double)x) + ((double)y)*((double)y) + ((double)z)*((double)z))
#define dist5(x, y, r, g, b) sqrt(((double)x)*((double)x) + ((double)y)*((double)y) + ((double)r)*((double)r) + ((double)g)*((double)g) + ((double)b)*((double)b))


