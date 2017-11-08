/*=========================================================/
/  makePNG module for PX425 2017 assignment 2.             /
/  Writes a PNG of the simulation grid                     /
/                                                          /
/  N. Hine - October 2017                                  /
/  Originally by D. Quigley                                /
/=========================================================*/

#include <png.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "makePNG.h"

void writePNG(char *filename, double **grid, int height, int width){
  /*==================================================================/
  / Function to construct a PNG image from a 2D of data points in     /
  / range -1.0 to +1.0 representing the order parameter phi in a      /
  / Ginzgurg-Landau systems.  Needs to be linked against libpng.      /
  /-------------------------------------------------------------------/
  / Original B&W version by S. Brown - University of Warwick          /
  / Modified to use colour scales by G. Enstone - also Warwick        /
  /==================================================================*/
  int x,y;

  /* Grid limits */
  double min_grid = +1.0;
  double max_grid = -1.0;

  for (y=0; y<height;++y){
    for (x=0;x<width;++x){
      if ( grid[y][x] > max_grid ) max_grid = grid[y][x];
      if ( grid[y][x] < min_grid ) min_grid = grid[y][x];      
    }
  }
  if (fabs(max_grid)>0.0)
  {
    if (fabs(min_grid)>fabs(max_grid)) max_grid = -min_grid;
    else min_grid = -max_grid;
  }
  else min_grid = -max_grid;
  
  /* Open output file */
  FILE *fp = fopen(filename, "wb");
  png_byte **row_pointers;
  if (!fp) bork("Couldn't open %s for writing.\n",filename);

  /* Set colour scale */
  const int palette_size    = 40;    
  const int pixel_bit_depth = 8;
  png_color palette[palette_size];
  png_color colour;

  int colors[][3] =   { {71,46,230}, /* blue */
                       {56,46,230},
                       {46,50,230},
                       {46,65,230},
                       {46,80,230},
                       {46,95,230},
                       {46,109,230},
                       {46,124,230},
                       {46,139,230},
                       {46,154,230},
                       {46,169,230},
                       {46,183,230},
                       {46,198,230},
                       {46,213,230},
                       {46,228,230},
                       {46,230,217},
                       {46,230,202},
                       {46,230,187},
                       {46,230,172},
                       {46,230,158},
                       {46,230,143},
                       {46,230,128},
                       {46,230,113},
                       {46,230,99},
                       {46,230,84},
                       {46,230,69},
                       {46,230,54},
                       {52,230,46},
                       {67,230,46},
                       {82,230,46},
                       {97,230,46},
                       {111,230,46},
                       {126,230,46},
                       {141,230,46},
                       {156,230,46},
                       {170,230,46},
                       {185,230,46},
                       {200,230,46},
                       {215,230,46},
                       {230,230,46} } ;/* yellow */
  
  int icol;
  for (icol=0;icol<palette_size;icol++) {
    colour = (png_color){colors[icol][0],colors[icol][1],colors[icol][2]};
    palette[icol] = colour;
  }

  
  /* Set up the PNG file */
  png_structp png_ptr = png_create_write_struct
    (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) bork("Couldn't allocate PNG.\n");
  
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) bork("Couldn't allocate PNG info.\n");
  
  if (setjmp (png_jmpbuf (png_ptr))) bork("PNG error: long_jump\n");
  
  /* Set image attributes. */
  png_set_IHDR (png_ptr,
		info_ptr,
		width,
		height,
                pixel_bit_depth, 
		PNG_COLOR_TYPE_PALETTE,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
  
  /*png_set_packing(png_ptr); */
  png_set_PLTE(png_ptr, info_ptr, palette, palette_size);

  
  /* Initialise PNG rows */
  row_pointers = png_malloc (png_ptr, sizeof(png_byte *)*height);	
  for (y=0; y<height; ++y){
    png_byte *row = png_malloc(png_ptr, sizeof(png_byte)*width);
    row_pointers[y]=row;
    for (x=0;x<width;++x){
      row[x] = (png_byte) (((grid[y][x]-min_grid)/(max_grid-min_grid))*(double)palette_size); 
      if (grid[y][x] >= max_grid) row[x] = palette_size - 1;
      if (grid[y][x] <= min_grid) row[x] = 0;

    }
  }
  
  /* Output code */
  png_init_io(png_ptr, fp);
  png_set_rows(png_ptr, info_ptr, row_pointers);
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_PACKING, NULL);
  png_write_end(png_ptr,info_ptr);
  
  /* Tidy up */
  fclose(fp);
  for (y=0;y<height;y++){
    png_free(png_ptr, row_pointers[y]);
  }
  png_free(png_ptr, row_pointers);
  png_destroy_write_struct(&png_ptr,&info_ptr);
}


void bork(char *msg,...){
  va_list args;
  va_start(args,msg);
  vfprintf(stderr,msg,args);
  exit(EXIT_FAILURE);
}


