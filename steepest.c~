#include <stdio.h> 
#include <stdlib.h>
#include <math.h>

/* function of f which is defined properly below  and it's derivatives */
double f(double x, double y);

double dxf(double x, double y);

double dyf(double x, double y);

int main(void) {

	/* define independent variables x and y*/
	double x0, y0;

	/* define t and gamma*/
	int t;
	double gamma = 0.01;


	/* read in input parameters */
	if(scanf("%lf %lf %d", &x0, &y0, &t) != 3){
		printf("error inputting value \n");
		return 1;
	};
	
	if(t > 4 || t < 1){
		printf("invalid value of t given");
		return 1;
	}


	if(x0>2.0 || x0<0.0 ||y0>2.0 || y0<0.0){
		printf("unsatisfactory initial conditions given");
		return 1;
	} 
	
	/*values required for convergence of algorithm*/
	double TOL = 1e-15;
	int MAX_ITER = 200;
	int counter = 1;
	
	/* set xN and yN to arbitrary values to start while loop*/
	double xN = x0 - (gamma*dxf(x0,y0));
	double yN = y0 - (gamma*dyf(x0,y0));
	

	while((fabs(f(xN,yN)-f(x0,y0))>TOL || fabs(x0 - xN)>TOL || fabs(y0-yN)>TOL) && counter<MAX_ITER){

		/* steepest descent in x */
		x0 = xN;
		xN = x0 - (gamma*dxf(x0,y0));

		/* steepest descent in y */
		y0 = yN;
		yN = y0 - (gamma*dyf(x0,y0));
		
		/*Apply counter*/
		counter = counter + 1;
	}
	
	if(counter == MAX_ITER){
		printf("Algorithm didn't converge");
	}else{
		printf("Minimum found at	%.15lf, %.15lf in %d iterations\n", xN, yN, counter);
	}
	
	return 0;
}
/*function of f*/
double f(double x, double y){
	double x2 = x*x;
	double x3 = x2*x;
	double y2 = y*y;
	double y3 = y2*y;
	return ((2*x3)-(3*x2) + 5)*((2*y3) -(3*y2) + 5);
}

/*derivative of f in x*/
double dxf(double x, double y){
	double y2 = y*y;
	double y3 = y2*y;
	return 6*((x*x) - x)*((2*y3)-(3*y2)+5);
}
/*derivative of f in y*/
double dyf(double x, double y){
	double x2 = x*x;
	double x3 = x2*x;
	return 6*((y*y) - y)*((2*x3)-(3*x2)+5);
}
