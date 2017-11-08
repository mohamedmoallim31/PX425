/*==========================================================//
//  ANSI-C code (optimised) for PX425 assignment 2 2017     //
//  Evolves a function H on a 2D grid via a discretisation  //
//  of the wave equation.                                   //
//                                                          //
//  d^2 u          d^2 u   d^2 u                            //
//  -----  =  c^2  ----- + -----                            //
//  d t^2          dx^2    dy^2                             //
//                                                          //
//  Original code created by D. Quigley                     //
//  Updated for 2017 by N. Hine                             //
//==========================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "makePNG.h"
#include "mt19937ar.h"

/* Function prototypes for memory management routines */
void allocate2d(double ***a,int num_rows,int num_cols);
void free2d(double ***a,int num_rows);

/* Function prototype for function which adds a new drop */
void addDrop(double ***u,double ***u_old,double dx,double dy,int Nx,int Ny,double height);

int main () {

  /* Function u on old new and current grid */
  double **u_new, **u_old, **u;

  /* Real variables (double precision) */
  double pudy,pudx;

  /* Number of grid points */
  int Nx = 480;
  int Ny = 640;

  /* Loop counters */
  int ix,iy,istep;

  /* Filename to which the grid is drawn */
  int  isnap=0;
  char filename[25];

  /*--------------*/
  /* Initial time */
  /*--------------*/
  clock_t t1 = clock();

  /*------------------------------------*/
  /* Initialise random number generator */
  /*------------------------------------*/
  unsigned long seed = 120492383972;
  init_genrand(seed);

  /*--------------------------*/
  /* Set grid 3x4 units       */
  /*--------------------------*/
  double dx = 3.0/(double)(Nx-1);
  double dy = 4.0/(double)(Ny-1);

  /*-----------------------------*/
  /* Set timestep and wave speed */
  /*-----------------------------*/
  double dt = 0.001;
  double c  = 1.0;

  /*--------------------------*/
  /* Number of steps to run   */
  /*--------------------------*/
  int nstep = 5000;

  /*--------------------------------------*/
  /* Allocate memory for a bunch of stuff */
  /*--------------------------------------*/
  allocate2d(&u,Nx,Ny);
  allocate2d(&u_old,Nx,Ny);
  allocate2d(&u_new,Nx,Ny);
      
  /* Add a drop */
  addDrop(&u,&u_old,dx,dy,Nx,Ny,0.02);

  /*------------------------------------*/
  /* Write an image of the initial grid */
  /*------------------------------------*/
  int stepincr = 10; 
  sprintf(filename,"snapshot%08d.png",isnap); 
  writePNG(filename,u,Nx,Ny); 
  isnap++;

  /* setup time */
  clock_t t2 = clock();
  printf("Setup time                    : %15.6f seconds\n",(double)(t2-t1)/(double)CLOCKS_PER_SEC);
  t1 = t2;

  /*===============================*/
  /* BEGIN SECTION TO BE OPTIMISED */
  /*===============================*/
 
  /*------------------------------------------*/
  /* Loop over the number of output timesteps */
  /*------------------------------------------*/
  double dx2 = dx*dx;
  double dy2 = dy*dy;
  double dt2 = dt*dt;
  double dInc = dt2/(dx2*dy2);
  for (istep=1;istep<nstep;istep++) {
    
    /* Randomly add new drops */
    if (genrand()<0.0005) {
      double height = genrand()*0.02;
      addDrop(&u,&u_old,dx,dy,Nx,Ny,height);
    }

    /*-----------------------*/
    /* Loop over grid points */
    /*-----------------------*/

    /*-----------------------------------------*/
    /* Loop over grid points away from boundary*/
    /*-----------------------------------------*/
    for(ix=1;ix<Nx-1;ix++) {
		for(iy=1;iy<Ny-1;iy++) {
			/* compute d2u/dy2 and d2u/dx2 away from boundaries*/
			pudy = (u[ix][iy+1] + u[ix][iy-1] - 2.0*u[ix][iy])*dx2;
			pudx = (u[ix+1][iy] + u[ix-1][iy] - 2.0*u[ix][iy])*dy2;
			/* compute new value at this grid point */
			u_new[ix][iy] = 2.0*u[ix][iy] - u_old[ix][iy] + dInc*(pudx+pudy); 
    	}
    }

    /*-----------------------------------------*/
    /* Loop over grid points ix==0 boundary    */
    /*-----------------------------------------*/
    pudx = (u[1][0] + u[Nx-1][0] - 2.0*u[0][0])*dy2;   /* lhs */
    pudy = (u[0][1] + u[0][Ny-1] - 2.0*u[0][0])*dx2;
    u_new[0][0] = 2.0*u[0][0] - u_old[0][0] + dInc*(pudx+pudy); 
    for(iy = 1; iy<Ny-1; iy++){
		pudx = (u[1][iy] + u[Nx-1][iy] - 2.0*u[0][iy])*dy2;   /* lhs */
		pudy = (u[0][iy+1] + u[0][iy-1] - 2.0*u[0][iy])*dx2;
		u_new[0][iy] = 2.0*u[0][iy] - u_old[0][iy] + dInc*(pudx+pudy); 
    }
    pudx = (u[1][Ny-1] + u[Nx-1][Ny-1] - 2.0*u[0][Ny-1])*dy2;   /* lhs */
    pudy = (u[0][0] + u[0][Ny-2] - 2.0*u[0][Ny-1])*dx2;   /* top    */
    u_new[0][Ny-1] = 2.0*u[0][Ny-1] - u_old[0][Ny-1] + dInc*(pudx+pudy); 

    /*-----------------------------------------*/
    /* Loop over grid points iy==0 boundary    */
    /*-----------------------------------------*/
    for(ix = 1; ix<Nx-1; ix++){
		pudy = (u[ix][1] + u[ix][Ny-1] - 2.0*u[ix][0])*dx2;   /* bottom */
		pudx = (u[ix+1][0] + u[ix-1][0] - 2.0*u[ix][0])*dy2;
		u_new[ix][0] = 2.0*u[ix][0] - u_old[ix][0] + dInc*(pudx+pudy); 
    }
    pudx = (u[0][0] + u[Nx-2][0] - 2.0*u[Nx-1][0])*dy2;   /* rhs */
    pudy = (u[Nx-1][1] + u[Nx-1][Ny-1] - 2.0*u[Nx-1][0])*dx2;   /* bottom */
    u_new[Nx-1][0] = 2.0*u[Nx-1][0] - u_old[Nx-1][0] + dInc*(pudx+pudy);

    /*--------------------------------------------*/
    /* Loop over grid points ix==Nx-1 boundary    */
    /*--------------------------------------------*/
    for(iy = 1; iy<Ny-1; iy++){
		pudx = (u[0][iy] + u[Nx-2][iy] - 2.0*u[Nx-1][iy])*dy2;   /* rhs */
		pudy = (u[Nx-1][iy+1] + u[Nx-1][iy-1] - 2.0*u[Nx-1][iy])*dx2;
		u_new[Nx-1][iy] = 2.0*u[Nx-1][iy] - u_old[Nx-1][iy] + dInc*(pudx+pudy); 
    }
    pudx = (u[0][Ny-1] + u[Nx-2][Ny-1] - 2.0*u[Nx-1][Ny-1])*dy2;   /* rhs */
    pudy = (u[Nx-1][0] + u[Nx-1][Ny-2] - 2.0*u[Nx-1][Ny-1])*dx2;   /* top    */
    u_new[Nx-1][Ny-1] = 2.0*u[Nx-1][Ny-1] - u_old[Nx-1][Ny-1] + dInc*(pudx+pudy);

    /*--------------------------------------------*/
    /* Loop over grid points iy==Ny-1 boundary    */
    /*--------------------------------------------*/
    for(ix = 1; ix<Nx-1; ix++){
	pudy = (u[ix][0] + u[ix][Ny-2] - 2.0*u[ix][Ny-1])*dx2;   /* top    */
        pudx = (u[ix+1][Ny-1] + u[ix-1][Ny-1] - 2.0*u[ix][Ny-1])*dy2;
	u_new[ix][Ny-1] = 2.0*u[ix][Ny-1] - u_old[ix][Ny-1] + dInc*(pudx+pudy); 
    }
        
    /* Shunt u into u_old */
    double **array = u_old;
    u_old = u;

    /* Shunt u_new into u */
    u = u_new;
    u_new = array;

    /*-----------------------------*/
    /* Snapshots of grid to file   */
    /* You are allowed to disable      */
    /* this if it causes a performance */
    /* hit.                            */
    /*-----------------------------*/
 /*   if (istep==isnap) {
        sprintf(filename,"snapshot%08d.png",isnap); 
        writePNG(filename,u,Nx,Ny); 
        if (stepincr==1){ isnap += stepincr; } 
        else if (istep>=1000) {isnap=isnap+1000;}
        else if (istep>=100) {isnap=isnap+100;}
        else { isnap *= stepincr;}
    }  
*/
  }

  /*=============================*/
  /* END SECTION TO BE OPTIMISED */
  /*=============================*/
  
  /* calculation time */
  t2 = clock();
  printf("Time taken for %8d steps : %15.6f seconds\n",nstep,(double)(t2-t1)/(double)CLOCKS_PER_SEC);
    
  /*----------------------------------*/
  /* Write an image of the final grid */
  /*----------------------------------*/
  sprintf(filename,"snapshot%08d.png",istep); 
  writePNG(filename,u,Nx,Ny); 

  /*--------------------------------------------*/
  /* Write final time-evolved solution to file. */
  /*--------------------------------------------*/
  FILE *fp = fopen("final_grid.dat","w");
  if (fp==NULL) printf("Error opening final_grid.dat for output\n");

  for(ix=0;ix<Nx-1;ix++) {
    for(iy=0;iy<Ny-1;iy++) {
      /* x and y at the current grid points */
      double x = dx*(double)ix;
      double y = dy*(double)iy;
      fprintf(fp,"%8.4f %8.4f %8.4e\n",x,y,u[ix][iy]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  /* Release memory */
  free2d(&u,Nx);
  free2d(&u_new,Nx);
  free2d(&u_old,Nx);
  
  return 0;
  
}
  

/*===============================================*/
/* Adds a new drop at a random point on the grid */ 
/*===============================================*/
void addDrop(double ***H,double ***H_old,double dx,double dy,int Nx,int Ny,double height) {
  
  const double rdecay = 8.0;
  const double kwave  = 40.0;

  double x,y,xc,yc,dxc,dyc;
  double Lx,Ly,r;
  int ix,iy;

  Lx = dx*(double)(Nx-1);
  Ly = dy*(double)(Ny-1);

  xc = genrand()*Lx;
  yc = genrand()*Ly;

  double **H_loc = *H;
  double **H_loc_old = *H_old;

  for (ix=0;ix<Nx;ix++) {

    /* x at the current grid point */
    x = dx*(double)ix;

    for (iy=0;iy<Ny;iy++) {
      
      /* y at the current grid point */
      y = dy*(double)iy;

      /* distance to drop centre */
      dxc = x-xc ; dxc -=  Lx*(int)(dxc/Lx+0.5);
      dyc = y-yc ; dyc -=  Ly*(int)(dyc/Ly+0.5);
        
      /* droplet function */
      r   = sqrt(dyc*dyc+dxc*dxc);

      H_loc[ix][iy]     += height*cos(kwave*r)*exp(-rdecay*r);
      H_loc_old[ix][iy] += height*cos(kwave*r)*exp(-rdecay*r);


    } /* ix */

  } /* iy */
 
} /* addDrop */



/*===========================================*/
/* Auxilliary routines for memory management */ 
/*===========================================*/
void allocate2d(double ***a,int Nx,int Ny) {

  double **b_loc; 

  b_loc = (double **)calloc(Nx,sizeof(double *));
  if (b_loc==NULL) printf("calloc error in allocate2d\n"); 

  int iy;
  for (iy=0;iy<Nx;iy++) {
    
    b_loc[iy] = (double *)calloc(Ny,sizeof(double));
    if (b_loc[iy]==NULL) printf("calloc error for row %d of %d in allocate2d\n",iy,Nx);

  }

  *a = b_loc;

}

void free2d(double ***a,int Nx) {

  int iy;

  double **b_loc = *a;

  /* Release memory */
  for (iy=0;iy<Nx;iy++) { 
    free(b_loc[iy]);
  }
  free(b_loc);
  *a = b_loc;

}









