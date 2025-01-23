/*****************************************************************************
 * 2D MECHANIC-BASED MODEL FOR COLLECTIVE MIGRATION OF TANDEM PAIRS					 *
 * WITHOUT FLUID COMPUTATION    												                     *
 * A MORESE POTENTIAL AT THE INTERCELLULAR JUNCTION                          * 
 *                                                                           *
 * MDCK-LIKE DOUBLET SIMULATION                                              *                                                     
 *                                                                           *
 * NOTE:                                                                     * 
 * 		1. RELAXATION IS INCLUDED                                              *
 * 			 USE T_RELAX TO CONTROL THE END TIME FOR RELAXING TWO CELLS          *
 * 		2. CELL 1 = LEADER CELL, CELL 2 = TRAILER CELL                         *
 *                                                                           *
 * CREATED BY                                                                *
 * 			 Y.ZHANG 11/14/2024                                                  *
 ****************************************************************************/


#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits> 
#include <time.h> 
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <getopt.h>
#include <stdio.h>

using namespace std;

#include "structs_mem.h"

// TIME DISCRETIZATION
const double TMAX = 400.0;
const double T_RELAX = 10.0;
const double T_CELL2MOVE = 10.0;
const double tstep = 0.0001;
static int stop = floor(TMAX/tstep);
//static int stop = 10; /* RUNNING STATIC TESTS */

// CELL 1 PARAMETERS
const double k_cell_1 = 1000;
const double xi_1 = 0.05;
double swell_1 = 1.04;
const double pseudopod_critL_1 = 0.0; // critical length of pseudopod 

// CELL 2 PARAMETERS
const double k_cell_2 = 1000;
const double xi_2 = 0.05;
double swell_2 = 1.04; //0.4
const double pseudopod_critL_2 = 0.0; // critical length of pseudopod 

// CELL 1 MEMBRANE/CORTEX PARAMETERS
double k_mem_1 = 1.0; 					    // membrane spring stiffness (units of force)
double gamma_mem_1 = 1.0; 			    // resting tension in membrane (units of force) 
const double k_contact_1 = 100.0;		// contact constant (units of force/L^2)

// CELL 2 MEMBRANE/CORTEX PARAMETERS
double k_mem_2 = 1.0; 					    // membrane spring stiffness (units of force)
double gamma_mem_2 = 1.0; 			    // resting tension in membrane (units of force) 
const double k_contact_2 = 100.0; 	// contact constant (units of force/L^2)

// CELL 1 ACTIN POLYMERIZATION PARAMETERS
const int tip_poly_1 = 50;
const double setArcLength_1 = 10.0;

// CELL 2 ACTIN POLYMERIZATION PARAMETERS
const int tip_poly_2 = 50;
const double setArcLength_2 = 10.0;

// CELL 1 ADHESION BOND PARAMETERS
double adh_prop_1 = 1.0; 					
double crit_adh_1 = 10.5;            // this is a critical adhesion density!!
double kminus_bond_1 = 0.16;				
double a_bond_1 = 0.1;									
double e_bond_1 = 0.008;

// CELL 2 ADHESION BOND PARAMETERS
double adh_prop_2 = 1.0; 					
double crit_adh_2 = 10.5;           // this is a critical adhesion density!!
double kminus_bond_2 = 0.16;						
double a_bond_2 = 0.1;									
double e_bond_2 = 0.008;	

// INTERCELLULAR JUNCTION PARAMETERS
double junction_rep = 10.0;         // R
double k_junction = 5.0;            // A
double junction_rest_L = 1.0;       // r
double junction_dist = 3.0;         // a

// SURFACE PARAMETERS
const double channel_height = 5.0;
const double attach_dl0 = 0.4;

// Boilerplate
#define MAX_OUTFILE 162
#define OUTFILE1    "morse_cell1_all"
#define OUTFILE2		"morse_cell1_surface"
#define OUTFILE3    "morse_cell2_all"
#define OUTFILE4		"morse_cell2_surface"
char *outfile1;
char *outfile2;
char *outfile3;
char *outfile4;
bool dbg = false;

//---------------------------------------------------------------------------------------------------------//

// COMPUTE DIFFERENCE IN CLOCKTIME IN SECONDS
double diffclock(clock_t s, clock_t e) {
    double diffticks = s-e;
    double diffms = (diffticks)/CLOCKS_PER_SEC;
    return diffms;
}
// END COMPUTE DIFF IN CLOCKTIME IN SECONDS

// QSORT COMPARE FUNCTION (DECREASING ORDER SORT)
int cmpfunc (const void * x, const void * y)
{
	double xx = *(double*)x;
	double yy = *(double*)y;
  if (xx < yy) return 1;
  if (xx > yy) return -1;
	return 0;
}
// END QSORT COMPARE FUNCTION

// COMPUTE MAX
double computeMax(int Npts, double arr[]) {
	double max;
	max = arr[0];
	
	for(int i=1; i<Npts; i++) {
		if(arr[i]>max) {
			max = arr[i];
		}
	}
	
	return max;
}
// END COMPUTE MAX

// COMPUTE MAX INDEX
double computeMaxIdx(int Npts, double arr[]) {
	int idx = 0;
	double max;
	max = arr[0];
	
	for(int i=1; i<Npts; i++) {
		if(arr[i]>max) {
			max = arr[i];
			idx = i;
		}
	}
	
	return idx;
}
// END COMPUTE MAX INDEX

// COMPUTE MIN
double computeMin(int Npts, double arr[]) {
	double min;
	min = arr[0];
	
	for(int i=1; i<Npts; i++) {
		if(arr[i]<min) {
			min = arr[i];
		}
	}
	
	return min;
}
// END COMPUTE MIN

// COMPUTE AREA OF CELL
double computeArea(int Npts, vertex Nodes[]) {
	double area = 0.;	
	
	int i, j=Npts-1;
	for(i=0; i<Npts; i++) {
		area += -(Nodes[j].def.x + Nodes[i].def.x)*(Nodes[j].def.y-Nodes[i].def.y);
	  // depending on orientation of points this might return a negative area (assume clockwise arrangement)
		j = i; 
	}

  return area*0.5;
}
// END COMPUTE AREA

// COMPUTE REFERENCE VECTORS
void computeRefDists(int Npts, vertex Nodes[], double ref_dist_L[], double ref_dist_R[] ){
	int k, left, right;

	for(k=0; k<Npts; k++) {
	// setup neighbors and membrane length information
		if(k==0) { left = Npts-1; right = k+1; }
    else if(k==Npts-1) { left = k-1; right = 0; }
    else { left = k-1; right = k+1; }

    ref_dist_L[k] = sqrt( pow(Nodes[k].ref.x - Nodes[left].ref.x,2.0) + pow(Nodes[k].ref.y - Nodes[left].ref.y,2.0) );
    ref_dist_R[k] = sqrt( pow(Nodes[k].ref.x - Nodes[right].ref.x,2.0) + pow(Nodes[k].ref.y - Nodes[right].ref.y,2.0) );
  
  }
}
// END COMPUTE REFERENCE VECTORS

// COMPUTE VECTORS
void computeVectors(int Npts, vertex Nodes[], double def_dist_L[], double def_dist_R[], vector2 t_L[], vector2 t_R[], vector2 n[]){
	int k, left, right;

	for(k=0; k<Npts; k++) {
	// setup neighbors and membrane length information
		if(k==0) { left = Npts-1; right = k+1; }
    else if(k==Npts-1) { left = k-1; right = 0; }
    else { left = k-1; right = k+1; }

    def_dist_L[k] = sqrt( pow(Nodes[k].def.x - Nodes[left].def.x,2.0) + pow(Nodes[k].def.y - Nodes[left].def.y,2.0) );
    def_dist_R[k] = sqrt( pow(Nodes[k].def.x - Nodes[right].def.x,2.0) + pow(Nodes[k].def.y - Nodes[right].def.y,2.0) );
	  Nodes[k].characteristicL = 0.5*(def_dist_L[k] + def_dist_R[k]);

	  t_L[k].x = (Nodes[left].def.x - Nodes[k].def.x)/def_dist_L[k];
    t_L[k].y = (Nodes[left].def.y - Nodes[k].def.y)/def_dist_L[k];
    t_R[k].x = (Nodes[right].def.x - Nodes[k].def.x)/def_dist_R[k];
    t_R[k].y = (Nodes[right].def.y - Nodes[k].def.y)/def_dist_R[k];

    n[k].x = 0.5*(-t_R[k].y*def_dist_R[k] + t_L[k].y*def_dist_L[k])/Nodes[k].characteristicL;
    n[k].y = 0.5*(t_R[k].x*def_dist_R[k] - t_L[k].x*def_dist_L[k])/Nodes[k].characteristicL;

  }

}
// END COMPUTE VECTORS

// COMPUTE ELASTIC FORCES
void computeElasticForces(int Npts, vertex Nodes[], double ref_dist_L[], double ref_dist_R[], double def_dist_L[], double def_dist_R[], vector2 t_L[], vector2 t_R[], vector2 ef[], double k_mem, double gamma_mem) {
	int k, left, right;
	double ds;
	
	for(k=0; k<Npts; k++) {
		// setup neighbors and membrane length information
		if(k==0) { left = Npts-1; right = k+1; }
    else if(k==Npts-1) { left = k-1; right = 0; }
    else { left = k-1; right = k+1; }

    ds = Nodes[k].characteristicL;

    // compute elastic forces
    ef[k].x = ((gamma_mem + k_mem*((def_dist_L[k] - ref_dist_L[k])/ref_dist_L[k]))*t_L[k].x + (gamma_mem + k_mem*((def_dist_R[k]-ref_dist_R[k])/ref_dist_R[k]))*t_R[k].x)/ds;
    ef[k].y = ((gamma_mem + k_mem*((def_dist_L[k] - ref_dist_L[k])/ref_dist_L[k]))*t_L[k].y + (gamma_mem + k_mem*((def_dist_R[k]-ref_dist_R[k])/ref_dist_R[k]))*t_R[k].y)/ds;
	}
}
// END COMPUTE ELASTIC FORCES

// COMPUTE P0 TERM
double computeP0Term(int Npts, vertex Nodes[], double def_dist_L[], double def_dist_R[], vector2 n[], vector2 f[]) {
	int k, left, right;
	double p0_top, p0_bottom, p0;		
	
	p0_top = 0.0;
	p0_bottom = 0.0;	
	for(k=0; k<Npts; k++) {
    // setup neighbors and membrane length information
    if(k==0) { left = Npts-1; right = k+1; }
    else if(k==Npts-1) { left = k-1; right = 0; }
    else { left = k-1; right = k+1; }

    p0_top += (f[k].x/Nodes[k].characteristicL)*n[k].x*Nodes[k].characteristicL + (f[k].y/Nodes[k].characteristicL)*n[k].y*Nodes[k].characteristicL;
    p0_bottom += Nodes[k].characteristicL;
	}
 
  p0 = p0_top/p0_bottom;
	return p0;
}
// END COMPUTE P0 TERM

// COMPUTE PRESSURE CONTRIBUTION
void computePressureForces(int Npts, vertex Nodes[], vector2 n[], vector2 pf[], double p) {
	int k, left, right;

  // add pressure force
  for(k=0; k<Npts; k++) {
  	if(k==0) { left = Npts-1; right = k+1; }
    else if(k==Npts-1) { left = k-1; right = 0; }
    else { left = k-1; right = k+1; }

    // careful: these are force densities
   	pf[k].x = -p*n[k].x;
    pf[k].y = -p*n[k].y;
	}
}
// END COMPUTE PRESSURE CONTRIBUTION

// COMPUTE CONTACT FORCES
void computeContactForces(int Npts, vertex Nodes[], vector2 f_contact[], double k_contact) {
	int k;
	double v_dist;

	for(k=0; k<Npts; k++) {
  	if(Nodes[k].bottom == 1) {
    	v_dist = fabs(Nodes[k].def.y)-0.1;
			if(v_dist<0) {
				f_contact[k].y = k_contact*fabs(v_dist)*Nodes[k].def.y/fabs(Nodes[k].def.y);
			}
		}
	}
}
// END COMPUTE CONTACT FORCES

// COMPUTE ADHESION FORCES
void computeAdhesionForces(int Npts, vertex Nodes[], vector2 f_adh[], double Nbonds[], double kminusTerm[], double yTerm[], double formationRate[], double adh_prop, double kminus_bond, double a_bond, int cell_label) {
	int k;
	double dist;
	double Nbonds_prev[Npts], ft_prev[Npts];

	if(cell_label == 1){
	  for(k=0; k<Npts; k++) {
	  	Nbonds_prev[k] = 0.0;
	    ft_prev[k] = 0.0;
	    kminusTerm[k] = 0.0;
	    yTerm[k] = 0.0;

			if(Nodes[k].sticky == 1) {
				dist = sqrt( pow(Nodes[k].def.x - Nodes[k].adh_link.x,2) + pow(Nodes[k].def.y - Nodes[k].adh_link.y,2) );
				if(Nodes[k].front == 1){
	      	Nodes[k].k_adh = 2.0*adh_prop*Nbonds[k];
	      } else {
	      	Nodes[k].k_adh = adh_prop*Nbonds[k];
	      }        

	      f_adh[k].x = Nodes[k].k_adh*((dist-attach_dl0)/attach_dl0)*(Nodes[k].adh_link.x - Nodes[k].def.x)/dist;
	      f_adh[k].y = Nodes[k].k_adh*((dist-attach_dl0)/attach_dl0)*(Nodes[k].adh_link.y - Nodes[k].def.y)/dist;

			 	// Catch/slip (sharp)
	      Nbonds[k] += tstep*( formationRate[k]*(1.0-Nbonds[k]) - Nbonds[k]*kminus_bond*exp(-a_bond*adh_prop*fabs((dist-attach_dl0)/attach_dl0)) );

				kminusTerm[k] = kminus_bond*exp(-a_bond*adh_prop*fabs((dist-attach_dl0)/attach_dl0));
	      yTerm[k] = adh_prop*fabs((dist-attach_dl0)/attach_dl0);
			}
		}
	} else {
		for(k=0; k<Npts; k++) {
	  	Nbonds_prev[k] = 0.0;
	    ft_prev[k] = 0.0;
	    kminusTerm[k] = 0.0;
	    yTerm[k] = 0.0;

			if(Nodes[k].sticky == 1) {
				dist = sqrt( pow(Nodes[k].def.x - Nodes[k].adh_link.x,2) + pow(Nodes[k].def.y - Nodes[k].adh_link.y,2) );
	     	Nodes[k].k_adh = adh_prop*Nbonds[k];       

	      f_adh[k].x = Nodes[k].k_adh*((dist-attach_dl0)/attach_dl0)*(Nodes[k].adh_link.x - Nodes[k].def.x)/dist;
	      f_adh[k].y = Nodes[k].k_adh*((dist-attach_dl0)/attach_dl0)*(Nodes[k].adh_link.y - Nodes[k].def.y)/dist;

			 	// Catch/slip (sharp)
	      Nbonds[k] += tstep*( formationRate[k]*(1.0-Nbonds[k]) - Nbonds[k]*kminus_bond*exp(-a_bond*adh_prop*fabs((dist-attach_dl0)/attach_dl0)) );

				kminusTerm[k] = kminus_bond*exp(-a_bond*adh_prop*fabs((dist-attach_dl0)/attach_dl0));
	      yTerm[k] = adh_prop*fabs((dist-attach_dl0)/attach_dl0);
			}
		}

	}
}
// END COMPUTE ADHESION FORCES

// COMPUTE LEADER POLYMERIZATION FORCE
void computeLeaderPolymerizationForces(int Npts, vertex c_Nodes[], vector2 c_f_poly[], double c_f_correctPolymerization[], double c_Nbonds[], vector2 c_ef[], double c_prev_pfx[], double c_prev_pfy[], double setArcLength, double c_delta_tip, double c_xFront, double swell, double pseudopod_critL, double e_bond, int t, double tstep, double T_RELAX){
	int k, j;
	double c_load_force = 0.0;
	double c_count_ds = 0.0;
	//double c_delta_tip = 4.0;
	double c_frontSticky = 0.0;
	double c_sum_ds = 0.0;
	double c_sum_swell = 0.0;
	double c_current_arcLength = 0.0;
	double c_distFromSticky = 0.0;

	int c_count_tip = 0;
	int c_count_poly = 0;
	
	double yFront = 3.0;

	for(k=0; k<Npts; k++) {
		if(c_Nodes[k].def.x >= (c_xFront-c_delta_tip) && c_Nodes[k].def.y >= 0.5 && c_Nodes[k].def.y <= yFront) {
			c_load_force += sqrt(pow(c_ef[k].x+c_prev_pfx[k],2)+pow(c_ef[k].y+c_prev_pfy[k],2)); 
			c_count_tip++;
			c_count_ds += c_Nodes[k].characteristicL;
		}
	}
	c_load_force = c_load_force/c_count_tip;

	for(k=0; k<Npts; k++) {
    if(c_Nodes[k].def.x >= (c_xFront-c_delta_tip) && c_Nodes[k].def.y >= 0.5 && c_Nodes[k].def.y <= yFront) { 
        c_Nodes[k].front = 1;
        c_count_poly++;      
        c_f_poly[k].x = swell*exp(-6.59*c_load_force)+c_load_force;
    }
    c_sum_swell += c_f_poly[k].x*c_Nodes[k].characteristicL;

    // determine x-position of the front-most adhesion site
    if(c_Nodes[k].sticky == 1 && c_Nodes[k].def.x>c_frontSticky) {
      c_frontSticky = c_Nodes[k].def.x;
    }

    if(c_Nodes[k].sticky == 1) {
			c_sum_ds += c_Nodes[k].characteristicL;
		}
	}

	c_distFromSticky = c_xFront - c_frontSticky;
	// BALANCE POLYMERIZATION FORCES
	for(k=0; k<Npts; k++) { 
		// balance polymerization forces
		if(c_Nodes[k].sticky == 1) {
			c_f_correctPolymerization[k] = -c_sum_swell/c_sum_ds; // a force density
		}
	}	

	// STICK FRONT WITH DISTANCE CONDITION
	if(c_distFromSticky > pseudopod_critL) {
		// form new adhesions once protrusion is "formed" enough
		for(k=0; k<Npts; k++) {
			if(c_Nodes[k].def.x>c_frontSticky && c_Nodes[k].def.y<5.0 && c_Nodes[k].sticky == 0) { // Nodes[k].def.y<5.0 //before
				c_Nodes[k].sticky = 1;
				c_Nodes[k].bottom = 1;
				c_Nodes[k].adh_link.x = c_Nodes[k].def.x+1.0;
				c_Nodes[k].adh_link.y = -0.5;
				c_Nbonds[k] = tstep*e_bond;
			}
		}
	}

	// REFORM ADHESIONS SITES 
	if((t*tstep) > T_RELAX){
		for(k=0; k<Npts; k++) {
			if(c_Nbonds[k]==0 && c_Nodes[k].def.x<c_frontSticky && c_Nodes[k].def.y<0.5) {
				c_Nodes[k].sticky = 1;
				c_Nodes[k].bottom = 1;
				c_Nodes[k].adh_link.x = c_Nodes[k].def.x+1.0;
				c_Nodes[k].adh_link.y = -0.5;
				c_Nbonds[k] = tstep*e_bond;
			}
		}
	}

}
// END COMPUTE LEADER POLYMERIZATION FORCE


// COMPUTE TRAILER POLYMERIZATION FORCES
void computeTrailerPolymerizationForces(int Npts, vertex c_Nodes[], vector2 c_f_poly[], double c_f_correctPolymerization[], double c_Nbonds[], vector2 c_ef[], double c_prev_pfx[], double c_prev_pfy[], double c_prev_repx[], double c_prev_repy[], double setArcLength, double c_delta_tip, double c_xFront, double swell, double pseudopod_critL, double e_bond, double c1_xBack, double c2_junction_ymin, int t, double tstep, double T_RELAX){
	int k, j;
	double c_load_force = 0.0;
	double c_count_ds = 0.0;
	double c_frontSticky = 0.0;
	double c_sum_ds = 0.0;
	double c_sum_swell = 0.0;
	double c_current_arcLength = 0.0;
	double c_distFromSticky = 0.0;

	int c_count_tip = 0;
	int c_count_poly = 0;
	
	double yFront = 3.0;

	for(k=0; k<Npts; k++) {
		if(c_Nodes[k].def.x >= (c_xFront-c_delta_tip) && c_Nodes[k].def.y >= 0.5 && c_Nodes[k].def.y <= yFront) {
			c_load_force += sqrt(pow(c_ef[k].x+c_prev_pfx[k]+c_prev_repx[k],2)+pow(c_ef[k].y+c_prev_pfy[k]+c_prev_repy[k],2)); 
			c_count_tip++;
			c_count_ds += c_Nodes[k].characteristicL;
		}
	}
	c_load_force = c_load_force/c_count_tip;

	for(k=0; k<Npts; k++) {
    if(c_Nodes[k].def.x >= (c_xFront-c_delta_tip) && c_Nodes[k].def.y >= 0.5 && c_Nodes[k].def.y <= yFront) { 
        c_count_poly++;      
        c_f_poly[k].x = swell*exp(-6.59*c_load_force)+c_load_force;
    }
    c_sum_swell += c_f_poly[k].x*c_Nodes[k].characteristicL;

    // determine x-position of the front-most adhesion site
    if(c_Nodes[k].sticky == 1 && c_Nodes[k].def.x>c_frontSticky) {
      c_frontSticky = c_Nodes[k].def.x;
    }

    if(c_Nodes[k].sticky == 1) {
			c_sum_ds += c_Nodes[k].characteristicL;
		}
	}

	c_distFromSticky = c_xFront - c_frontSticky;
	// BALANCE POLYMERIZATION FORCES
	for(k=0; k<Npts; k++) { 
		// balance polymerization forces
		if(c_Nodes[k].sticky == 1) {
			c_f_correctPolymerization[k] = -c_sum_swell/c_sum_ds; // a force density
		}
	}	

	// STICK FRONT WITH DISTANCE CONDITION
	if(c_distFromSticky > pseudopod_critL) {
		// form new adhesions once protrusion is "formed" enough
		for(k=0; k<Npts; k++) {
			if(c_Nodes[k].def.x>c_frontSticky && c_Nodes[k].def.y<c2_junction_ymin) {
				if(c_Nodes[k].def.x>=c1_xBack){
					c_Nodes[k].sticky = 0;
					c_Nodes[k].bottom = 0;
					c_Nbonds[k] = 0.0;
				} else {
					c_Nodes[k].sticky = 1;
					c_Nodes[k].bottom = 1;
					c_Nodes[k].adh_link.x = c_Nodes[k].def.x+1.0;
					c_Nodes[k].adh_link.y = -0.5;
					c_Nbonds[k] = tstep*e_bond;
				}
			}
		}
	}

	// REFORM ADHESIONS SITES 
	if((t*tstep) > T_RELAX){
		for(k=0; k<Npts; k++) {
			if(c_Nbonds[k]==0 && c_Nodes[k].def.x<c_frontSticky && c_Nodes[k].def.y<0.5) {
				c_Nodes[k].sticky = 1;
				c_Nodes[k].bottom = 1;
				c_Nodes[k].adh_link.x = c_Nodes[k].def.x+1.0;
				c_Nodes[k].adh_link.y = -0.5;
				c_Nbonds[k] = tstep*e_bond;
			}
		}
	}
}

// END COMPUTE TRAILER POLYMERIZATION FORCES


// COMPUTE INTERCELLULAR JUNCTION FORCES
void computeJunctionForces(int Npts, vertex c1_Nodes[], vertex c2_Nodes[], vector2 c1_f_junction[], vector2 c2_f_junction[], int c2_count_poly, double k_junction, double junction_rest_L, double junction_rep, vector2 c1_n[], vector2 c2_n[], double c2_junction_ymin, double c2_sum_fpoly, double c1_sum_prev_junction, double c2_sum_prev_junction, vector2 c1_f_repulsion[], vector2 c2_f_repulsion[], int check, vector2 c2_f_poly[]){
	int k, j, c1_idx, c2_idx;
	double dist, dist1, dist2;
	double correctds;
	double c1_f_attract = 0.0;
	double c2_f_attract = 0.0; 
	double c1_f_repel = 0.0;
	double c2_f_repel = 0.0;

	int c1connection[c2_count_poly];
	int c2connection[c2_count_poly];
	int count = 0;
	int rep_count = 0;
	double c2_f_poly_max = 0.0;

	for(k=0; k<Npts; k++){
		if(c2_f_poly_max < c2_f_poly[k].x){
			c2_f_poly_max = c2_f_poly[k].x;
		}
		if(c2_Nodes[k].front == 1){
			dist2 = 1000.0;
			for(j=0; j<Npts; j++){
				if(c1_Nodes[j].back == 1 && c1_Nodes[j].junction == 0){
					if(c1_Nodes[j].def.y >= c2_junction_ymin){
						if(check == 0){
							dist1 = c1_Nodes[j].def.y - c2_Nodes[k].def.y;
						} else {
							dist1 = sqrt( pow(c1_Nodes[j].def.x - c2_Nodes[k].def.x,2) + pow(c1_Nodes[j].def.y - c2_Nodes[k].def.y,2) );
						}
						if(dist1 < dist2){
							dist2 = dist1;
							c1connection[count] = j; 
							c2connection[count] = k;
						}
					}
				}
			}
			c1_Nodes[c1connection[count]].junction = 1;
			count++;
		}
	}

	for(k=0; k<c2_count_poly; k++){
		c1_idx = c1connection[k];
		c2_idx = c2connection[k];

		dist = (c1_Nodes[c1_idx].def.x - c2_Nodes[c2_idx].def.x);

		if(dist > 0.0){
			c1_f_junction[c1_idx].x += junction_rep*exp(-fabs(dist/junction_rest_L)) - k_junction*exp(-fabs(dist/junction_dist));
			c2_f_junction[c2_idx].x += -(junction_rep*exp(-fabs(dist/junction_rest_L)) - k_junction*exp(-fabs(dist/junction_dist)));
		} else if(dist < 0.0) {
			c1_f_junction[c1_idx].x += -(junction_rep*exp(-fabs(dist/junction_rest_L)) - k_junction*exp(-fabs(dist/junction_dist)));
			c2_f_junction[c2_idx].x += junction_rep*exp(-fabs(dist/junction_rest_L)) - k_junction*exp(-fabs(dist/junction_dist));
		} else {
			c1_f_junction[c1_idx].x += 0.0;
			c2_f_junction[c2_idx].x += 0.0;
		}

	}
}
// END COMPUTE JUNCTION FORCES

//-------------------------------------------------------------------------------------------------//

// PROGRESSION IN TIME
void progress(int Npts, vertex c1_Nodes[], vertex c2_Nodes[]) {
	int j, k, l;

	// CELL 1 VECTORS
	vector2 c1_ef[Npts];	
  vector2 c1_f_poly[Npts];
	vector2 c1_pf[Npts];
	vector2 c1_f[Npts];
	vector2 c1_f_contact[Npts];
	vector2 c1_f_adh[Npts];
	vector2 c1_f_wall[Npts];
	vector2 c1_f_junction[Npts];
	vector2 c1_f_repulsion[Npts];
	vector2 c1_v[Npts];
	double c1_f_correctPolymerization[Npts];
	double c1_x_bond[Npts];
  double c1_ft_prev[Npts];
  double c1_Nbonds[Npts];

  vector2 c1_t_R[Npts];
  vector2 c1_t_L[Npts];
  vector2 c1_n[Npts];
  double c1_ref_dist_L[Npts];
  double c1_ref_dist_R[Npts];
  double c1_def_dist_L[Npts];
  double c1_def_dist_R[Npts]; 

  // CELL 2 VECTORS
	vector2 c2_ef[Npts];	
  vector2 c2_f_poly[Npts];
	vector2 c2_pf[Npts];
	vector2 c2_f[Npts];
	vector2 c2_f_contact[Npts];
	vector2 c2_f_adh[Npts];
	vector2 c2_f_wall[Npts];
	vector2 c2_f_junction[Npts];
	vector2 c2_f_repulsion[Npts];
	vector2 c2_v[Npts];
	double c2_f_correctPolymerization[Npts];
	double c2_x_bond[Npts];
  double c2_ft_prev[Npts];
  double c2_Nbonds[Npts];

  vector2 c2_t_R[Npts];
  vector2 c2_t_L[Npts];
  vector2 c2_n[Npts];
  double c2_ref_dist_L[Npts];
  double c2_ref_dist_R[Npts];
  double c2_def_dist_L[Npts];
  double c2_def_dist_R[Npts]; 

  // BOOLEAN CONTROL
  bool isDone = false;

  // file handling
  ofstream f1, f2, f3, f4;
	ifstream f01, f02, f03;
	f1.precision(16);
  f2.precision(16);
  f3.precision(16);
  f4.precision(16);

  // Open equalibrium data for cell 1
	f01.open("Configurations/C1_pre-Equil_1_0.1_0.14_0.01_10_spring_posrelax_d_0.5_10.txt");

	// Open equalibrium data for cell 2
	f02.open("Configurations/C2_pre-Equil_1_0.1_0.14_0.01_10_spring_posrelax_d_0.5_10.txt");

	// Open files for saving data
	f1.open(outfile1);
  f2.open(outfile2);
  f3.open(outfile3);
  f4.open(outfile4);

	// file headers with parameters
	f1 << "Membrane parameters N = " << Npts << " k_mem = " << k_mem_1 << " gamma_mem = " << gamma_mem_1 << " k_cell = " << k_cell_1 << " xi = " << xi_1 << endl;
  f1 << "Polymerization parameters swell = " << swell_1 << " pseudopod_critL = " << pseudopod_critL_1 << endl;
	f1 << "Adhesion bond parameters adh_prop = " << adh_prop_1 << " crit_adh = " << crit_adh_1 << " kminus = " << kminus_bond_1 << " a_bond = " << a_bond_1 << " kplus = " << e_bond_1 << endl;
	f1 << "Cell-cell junction parameters k_junction = " << k_junction << " junction_rest_L = " << junction_rest_L << " junction_dist = " << junction_dist << endl;
  f1 << "---------------" << endl;

	f2 << "Membrane parameters N = " << Npts << " k_mem = " << k_mem_1 << " gamma_mem = " << gamma_mem_1 << " k_cell = " << k_cell_1 << " xi = " << xi_1 << endl;
  f2 << "Polymerization parameters swell = " << swell_1 << " pseudopod_critL = " << pseudopod_critL_1 << endl;
	f2 << "Adhesion bond parameters adh_prop = " << adh_prop_1 << " crit_adh = " << crit_adh_1 << " kminus = " << kminus_bond_1 << " a_bond = " << a_bond_1 << " kplus = " << e_bond_1 << endl;
	f1 << "Cell-cell junction parameters k_junction = " << k_junction << " junction_rest_L = " << junction_rest_L << " junction_dist = " << junction_dist << endl;
  f2 << "---------------" << endl;

  f3 << "Membrane parameters N = " << Npts << " k_mem = " << k_mem_2 << " gamma_mem = " << gamma_mem_2 << " k_cell = " << k_cell_2 << " xi = " << xi_2 << endl;
  f3 << "Polymerization parameters swell = " << swell_2 << " pseudopod_critL = " << pseudopod_critL_2 << endl;
	f3 << "Adhesion bond parameters adh_prop = " << adh_prop_2 << " crit_adh = " << crit_adh_2 << " kminus = " << kminus_bond_2 << " a_bond = " << a_bond_2 << " kplus = " << e_bond_2 << endl;
  f3 << "Cell-cell junction parameters k_junction = " << k_junction << " junction_rest_L = " << junction_rest_L << " junction_dist = " << junction_dist << endl;
  f3 << "---------------" << endl;

	f4 << "Membrane parameters N = " << Npts << "k_mem = " << k_mem_2 << " gamma_mem = " << gamma_mem_2 << " k_cell = " << k_cell_2 << " xi = " << xi_2 << endl;
  f4 << "Polymerization parameters swell = " << swell_2 << " pseudopod_critL = " << pseudopod_critL_2 << endl;
	f4 << "Adhesion bond parameters adh_prop = " << adh_prop_2 << " crit_adh = " << crit_adh_2 << " kminus = " << kminus_bond_2 << " a_bond = " << a_bond_2 << " kplus = " << e_bond_2 << endl;
  f3 << "Cell-cell junction parameters k_junction = " << k_junction << " junction_rest_L = " << junction_rest_L << " junction_dist = " << junction_dist << endl;
  f4 << "---------------" << endl;

	// zero the force
  for(int k=0; k<Npts; k++) {
  	// CELL 1
		c1_Nodes[k].characteristicL = 0.0;
		c1_Nodes[k].force.x = 0.0;
  	c1_Nodes[k].force.y = 0.0;
		c1_f[k].x = 0.0;
		c1_f[k].y = 0.0;

		// CELL 2
		c2_Nodes[k].characteristicL = 0.0;
		c2_Nodes[k].force.x = 0.0;
  	c2_Nodes[k].force.y = 0.0;
		c2_f[k].x = 0.0;
		c2_f[k].y = 0.0;
  } 

  // compute initial areas
	j=Npts-1;
  double c1_A0 = 0.;
  double c2_A0 = 0.;
	for(int k=0; k<Npts; k++) {
    c1_A0 += -0.5*(c1_Nodes[j].ref.x + c1_Nodes[k].ref.x)*(c1_Nodes[j].ref.y-c1_Nodes[k].ref.y);
    c2_A0 += -0.5*(c2_Nodes[j].ref.x + c2_Nodes[k].ref.x)*(c2_Nodes[j].ref.y-c2_Nodes[k].ref.y);
    j = k;
  }

	// read in initial bond setup
	double c1;
	int counter = 0;
	while(f01 >> c1) {
		c1_Nbonds[counter] = c1;
		counter++;
	}  
	f01.close();
	double c2;
	counter = 0;
	while(f02 >> c2) {
		c2_Nbonds[counter] = c2;
		counter++;
	}  
	f02.close();

	for(int k=0; k<Npts; k++) {
		c1_x_bond[k] = 1.0;
		c1_ft_prev[k] = 0.0;

		c2_x_bond[k] = 1.0;
		c2_ft_prev[k] = 0.0;
	}

	// MARCHING IN TIME
	counter = 0;
	int check = 0;
	double c1_xFront, c1_xBack, c1_frontSticky, c1_distFromSticky, c1_lengthCell, c1_halfLength;
	double c1_p, c1_currentA, c1_p0;
	double c1_sum_swell, c1_sum_ds, c1_delta_tip, c1_load_force;
	double c1_prev_pfx[Npts], c1_prev_pfy[Npts]; // used for polymerization model 
	double c1_xPos[Npts], c1_negXPos[Npts], c1_kminusTerm[Npts], c1_yTerm[Npts];
	double c1_formationRate[Npts];
	double c1_prev_vx[Npts], c1_prev_vy[Npts];
	double c1_prev_repx[Npts], c1_prev_repy[Npts];

	double c2_xFront, c2_xBack, c2_frontSticky, c2_distFromSticky, c2_lengthCell, c2_halfLength;
	double c2_p, c2_currentA, c2_p0;
	double c2_sum_swell, c2_sum_ds, c2_delta_tip, c2_current_arcLength, c2_load_force;
	double c2_prev_pfx[Npts], c2_prev_pfy[Npts]; // used for polymerization model 
	double c2_xPos[Npts], c2_refxPos[Npts], c2_negXPos[Npts], c2_kminusTerm[Npts], c2_yTerm[Npts];
	double c2_formationRate[Npts];
	double c2_prev_vx[Npts], c2_prev_vy[Npts];
	double c1_sum_junction, c2_sum_junction; 
	double c1_sum_fpoly, c2_sum_fpoly; 
	double c1_sum_prev_junction, c2_sum_prev_junction;
	double c2_prev_repx[Npts], c2_prev_repy[Npts];

	for(int k=0; k<Npts; k++) {
		c1_prev_vx[k] = 0.0;
		c1_prev_vy[k] = 0.0;
		c1_prev_pfx[k] = 0.0;	
		c1_prev_pfy[k] = 0.0;
		c1_prev_repx[k] = 0.0;
		c1_prev_repy[k] = 0.0;		

		c2_prev_vx[k] = 0.0;
		c2_prev_vy[k] = 0.0;
		c2_prev_pfx[k] = 0.0;	
		c2_prev_pfy[k] = 0.0;
		c2_prev_repx[k] = 0.0;
		c2_prev_repy[k] = 0.0;	
	}
	
	c1_sum_junction = 0.0;
	c2_sum_junction = 0.0;

	c1_sum_fpoly = 0.0;
	c2_sum_fpoly = 0.0;

	c1_sum_prev_junction = 0.0;
	c2_sum_prev_junction = 0.0;

	// COMPUTE REFERENCE DISTANCES
	computeRefDists(Npts, c1_Nodes, c1_ref_dist_L, c1_ref_dist_R);
	computeRefDists(Npts, c2_Nodes, c2_ref_dist_L, c2_ref_dist_R);

	// START TIME STEPPING
	for(int t=0; t<stop; t++) {

		c1_sum_swell = 0.0;
 		c1_sum_ds = 0.0;
 		c2_sum_swell = 0.0;
 		c2_sum_ds = 0.0;

		// ZERO OUT STUFF
		for(int k=0; k<Npts; k++) {
			// CELL 1
			c1_ef[k].x = 0.0;
			c1_ef[k].y = 0.0;
			c1_f_poly[k].x = 0.0;
			c1_f_poly[k].y = 0.0;
      c1_f_contact[k].x = 0.0;
      c1_f_contact[k].y = 0.0;
      c1_f_adh[k].x = 0.0;
      c1_f_adh[k].y = 0.0;
      c1_pf[k].x = 0.0;
			c1_pf[k].y = 0.0;
			c1_f_wall[k].x = 0.0;
      c1_f_wall[k].y = 0.0;
      c1_f_junction[k].x = 0.0;
      c1_f_junction[k].y = 0.0;
      c1_f_repulsion[k].x = 0.0;
      c1_f_repulsion[k].y = 0.0;

      c1_t_R[k].x = 0.0;
      c1_t_R[k].y = 0.0;
  		c1_t_L[k].x = 0.0;
  		c1_t_L[k].y = 0.0;
  		c1_n[k].x = 0.0;
  		c1_n[k].y = 0.0;
  		c1_def_dist_L[k] = 0.0;
  		c1_def_dist_R[k] = 0.0; 
		
			c1_v[k].x = 0.0;
			c1_v[k].y = 0.0;
	
			c1_sum_junction = 0.0;
			c1_sum_fpoly = 0.0;
			c1_sum_prev_junction = 0.0;
			c1_load_force = 0.0;

			c1_f_correctPolymerization[k] = 0.0;
			c1_Nodes[k].front = 0;
			c1_Nodes[k].junction = 0;	
			c1_Nodes[k].back = 0;	
			c1_formationRate[k] = e_bond_1;
			c1_Nodes[k].flag = 0;

			c1_Nodes[k].bottom = 0; 
			if(c1_Nodes[k].def.y < 0.5) { c1_Nodes[k].bottom = 1; }

			c1_xPos[k] = c1_Nodes[k].def.x;
			c1_negXPos[k] = -c1_Nodes[k].def.x;

			// CELL 2
			c2_ef[k].x = 0.0;
			c2_ef[k].y = 0.0;
			c2_f_poly[k].x = 0.0;
			c2_f_poly[k].y = 0.0;
      c2_f_contact[k].x = 0.0;
      c2_f_contact[k].y = 0.0;
      c2_f_adh[k].x = 0.0;
      c2_f_adh[k].y = 0.0;
      c2_pf[k].x = 0.0;
			c2_pf[k].y = 0.0;
			c2_f_wall[k].x = 0.0;
      c2_f_wall[k].y = 0.0;
      c2_f_junction[k].x = 0.0;
      c2_f_junction[k].y = 0.0;
      c2_f_repulsion[k].x = 0.0;
      c2_f_repulsion[k].y = 0.0;

      c2_t_R[k].x = 0.0;
      c2_t_R[k].y = 0.0;
  		c2_t_L[k].x = 0.0;
  		c2_t_L[k].y = 0.0;
  		c2_n[k].x = 0.0;
  		c2_n[k].y = 0.0;
  		c2_def_dist_L[k] = 0.0;
  		c2_def_dist_R[k] = 0.0; 
		
			c2_v[k].x = 0.0;
			c2_v[k].y = 0.0;
	
			c2_sum_junction = 0.0;
			c2_sum_fpoly = 0.0;
			c2_sum_prev_junction = 0.0;
			c2_load_force = 0.0;

			c2_f_correctPolymerization[k] = 0.0;
			c2_Nodes[k].front = 0;	
			c2_Nodes[k].junction = 0;	
			c2_formationRate[k] = e_bond_2;
			c2_Nodes[k].flag = 0;

			c2_Nodes[k].bottom = 0; 
			if(c2_Nodes[k].def.y < 0.5) { c2_Nodes[k].bottom = 1; }

			c2_xPos[k] = c2_Nodes[k].def.x;
			c2_refxPos[k] = c2_Nodes[k].ref.x;
			c2_negXPos[k] = -c2_Nodes[k].def.x;
		}

		// KEEP TRACK OF OLD FRONT OF CELL 2
    if(t==0) {
      c2_frontSticky = computeMax(Npts,c2_refxPos);
    } else if (check==1) {
      c2_frontSticky = computeMax(Npts,c2_xPos);
      check = 0;
    }

		if(t==0) {
			for(int k=0; k<Npts; k++) {
				f1 << t*tstep << " " << c1_Nodes[k].def.x << " " << c1_Nodes[k].def.y << " " << c1_Nodes[k].ref.x << " " << c1_Nodes[k].ref.y << " " << c1_ef[k].x << " " << c1_ef[k].y << " " << c1_f_contact[k].x << " " << c1_f_contact[k].y << " " << c1_f_adh[k].x << " " << c1_f_adh[k].y << " " << c1_pf[k].x << " " << c1_pf[k].y << " " << c1_f[k].x << " " << c1_f[k].y << " " << c1_Nodes[k].characteristicL << " " << c1_f_correctPolymerization[k] << " " << c1_Nodes[k].flag << " " << c1_Nodes[k].sticky << " " << c1_Nbonds[k] << " " << c1_kminusTerm[k] << " " << c1_yTerm[k] << " " << c1_Nodes[k].front << " " << c1_f_poly[k].x << " " << c1_f_poly[k].y << " " << c1_f_junction[k].x << " " << c1_f_junction[k].y << " " << c1_f_repulsion[k].x << " " << c1_f_repulsion[k].y << endl;
    		f2 << t*tstep << " " << c1_Nodes[k].adh_link.x << " " << c1_Nodes[k].adh_link.y << " "<< c1_f_wall[k].x << " " << c1_f_wall[k].y << " " << c1_Nodes[k].sticky << endl;	

    		f3 << t*tstep << " " << c2_Nodes[k].def.x << " " << c2_Nodes[k].def.y << " " << c2_Nodes[k].ref.x << " " << c2_Nodes[k].ref.y << " " << c2_ef[k].x << " " << c2_ef[k].y << " " << c2_f_contact[k].x << " " << c2_f_contact[k].y << " " << c2_f_adh[k].x << " " << c2_f_adh[k].y << " " << c2_pf[k].x << " " << c2_pf[k].y << " " << c2_f[k].x << " " << c2_f[k].y << " " << c2_Nodes[k].characteristicL << " " << c2_f_correctPolymerization[k] << " " << c2_Nodes[k].flag << " " << c2_Nodes[k].sticky << " " << c2_Nbonds[k] << " " << c2_kminusTerm[k] << " " << c2_yTerm[k] << " " << c2_Nodes[k].front << " " << c2_f_poly[k].x << " " << c2_f_poly[k].y << " " << c2_f_junction[k].x << " " << c2_f_junction[k].y << " " << c2_f_repulsion[k].x << " " << c2_f_repulsion[k].y << endl;
    		f4 << t*tstep << " " << c2_Nodes[k].adh_link.x << " " << c2_Nodes[k].adh_link.y << " "<< c2_f_wall[k].x << " " << c2_f_wall[k].y << " " << c2_Nodes[k].sticky << endl;	
			}
		}

		// COMPUTE VECTORS
		computeVectors(Npts, c1_Nodes, c1_def_dist_L, c1_def_dist_R, c1_t_L, c1_t_R, c1_n);
		computeVectors(Npts, c2_Nodes, c2_def_dist_L, c2_def_dist_R, c2_t_L, c2_t_R, c2_n);

		// COMPUTE ELASTIC FORCES
		computeElasticForces(Npts, c1_Nodes, c1_ref_dist_L, c1_ref_dist_R, c1_def_dist_L, c1_def_dist_R, c1_t_L, c1_t_R, c1_ef, k_mem_1, gamma_mem_1);
		computeElasticForces(Npts, c2_Nodes, c2_ref_dist_L, c2_ref_dist_R, c2_def_dist_L, c2_def_dist_R, c2_t_L, c2_t_R, c2_ef, k_mem_2, gamma_mem_2);

		// DETERMINE LENGTH OF PSEUDOPOD FROM THE FRONT-MOST ADHESION SITE
		c1_xFront = computeMax(Npts, c1_xPos);
		c1_xBack = (-1.0)*computeMax(Npts, c1_negXPos);
		c2_xFront = computeMax(Npts, c2_xPos);
		c2_xBack = (-1.0)*computeMax(Npts, c2_negXPos);

		double c1_xStickyBack = 100000.0;
		for(int k=0; k<Npts; k++){
			if(c1_Nodes[k].sticky == 1 && c1_Nodes[k].def.x<c1_xStickyBack) {
      	c1_xStickyBack = c1_Nodes[k].def.x;
    	}
		}

		c1_lengthCell = c1_xFront - c1_xBack;
		c2_lengthCell = c2_frontSticky - c2_xBack;
		setArcLength_2 = 0.9*setArcLength_1;
		c1_halfLength = c1_xBack + 0.5*c1_lengthCell;
		c2_halfLength = c2_xBack + 0.25*c2_lengthCell;

		c1_delta_tip = 4.0;
		c2_delta_tip = 4.0;
		c2_current_arcLength = 0.0;
		double c2_junction_ymin = 0.0;
		int c2_count_junction = 0;

		// DETERMINE THE CONTACT REGION POINTS
		for(int k=0; k<Npts; k++){
			// CELL 1
			if(c1_Nodes[k].def.x <= c1_halfLength && c1_Nodes[k].def.y >= 0.5){
				c1_Nodes[k].back = 1;
			}

			// CELL 2
			if(c2_Nodes[k].def.x >= (c2_xFront-c2_delta_tip) && c2_Nodes[k].def.y >= 0.5 && c2_Nodes[k].def.y <= 3.0){
				c2_current_arcLength += c2_Nodes[k].characteristicL;
				if(c2_current_arcLength <= setArcLength_2) {
          	c2_Nodes[k].front = 1;
          	c2_count_junction++;

          if(c2_count_junction == 1){
          	c2_junction_ymin = c2_Nodes[k].def.y;
          } else {
          	if(c2_Nodes[k].def.y < c2_junction_ymin){
          		c2_junction_ymin = c2_Nodes[k].def.y;
          	}
          }
        }
			}

			if(c2_Nodes[k].def.x >= (c2_xFront-c2_delta_tip) && c2_Nodes[k].def.y > 3.0 && c2_Nodes[k].def.y < 6.0 ){
				c2_current_arcLength += c2_Nodes[k].characteristicL;
				if(c2_current_arcLength <= setArcLength_2) {
          c2_Nodes[k].front = 1;
          c2_count_junction++;
        }
      }
		}

		// COMPUTE POLYMERIZATION FORCES
		check = 1;
		if((t*tstep) > T_RELAX && (t*tstep) < T_CELL2MOVE){
			computeLeaderPolymerizationForces(Npts, c1_Nodes, c1_f_poly, c1_f_correctPolymerization, c1_Nbonds, c1_ef, c1_prev_pfx, c1_prev_pfy, setArcLength_1, c1_delta_tip, c1_xFront, swell_1, pseudopod_critL_1, e_bond_1, t, tstep, T_RELAX);
		} else if((t*tstep) >= T_CELL2MOVE){
			computeLeaderPolymerizationForces(Npts, c1_Nodes, c1_f_poly, c1_f_correctPolymerization, c1_Nbonds, c1_ef, c1_prev_pfx, c1_prev_pfy, setArcLength_1, c1_delta_tip, c1_xFront, swell_1, pseudopod_critL_1, e_bond_1, t, tstep, T_RELAX);
			computeTrailerPolymerizationForces(Npts, c2_Nodes, c2_f_poly, c2_f_correctPolymerization, c2_Nbonds, c2_ef, c2_prev_pfx, c2_prev_pfy, c2_prev_repx, c2_prev_repy, setArcLength_2, c2_delta_tip, c2_xFront, swell_2, pseudopod_critL_2, e_bond_2, c1_xStickyBack, c2_junction_ymin, t, tstep, T_RELAX);
		} else {
			check = 0;
		}

		// COMPUTE THE TOTAL POLYMERIZATION FORCES
		for(int k=0; k<Npts; k++){
			c1_sum_fpoly += sqrt( pow(c1_f_poly[k].x, 2) + pow(c1_f_poly[k].y, 2) );
			c2_sum_fpoly += sqrt( pow(c2_f_poly[k].x, 2) + pow(c2_f_poly[k].y, 2) );
		}

		// COMPUTE JUNCTION FORCES
		computeJunctionForces(Npts, c1_Nodes, c2_Nodes, c1_f_junction, c2_f_junction, c2_count_junction, k_junction, junction_rest_L, junction_rep, c1_n, c2_n, c2_junction_ymin, c2_sum_fpoly, c1_sum_prev_junction, c2_sum_prev_junction, c1_f_repulsion, c2_f_repulsion, check, c2_f_poly);

		// COMPUTE ADHESION FORCE DENSITIES
		computeAdhesionForces(Npts, c1_Nodes, c1_f_adh, c1_Nbonds, c1_kminusTerm, c1_yTerm, c1_formationRate, adh_prop_1, kminus_bond_1, a_bond_1, 1);
		computeAdhesionForces(Npts, c2_Nodes, c2_f_adh, c2_Nbonds, c2_kminusTerm, c2_yTerm, c2_formationRate, adh_prop_2, kminus_bond_2, a_bond_2, 2);

		// COMPUTE CONTACT FORCE DENSITIES
    computeContactForces(Npts, c1_Nodes, c1_f_contact, k_contact_1);
    computeContactForces(Npts, c2_Nodes, c2_f_contact, k_contact_2);

		// BALANCE POLYMERIZATION FORCES AND SUM UP ALL FORCES
		// RUPTURE ADHESIONS ABOVE DENSITY THRESHOLD
		for(int k=0; k<Npts; k++) {
			if(c1_Nbonds[k]<=0.000001) { c1_Nbonds[k] = 0.0; c1_Nodes[k].sticky = 0; c1_f_adh[k].x = 0.0; c1_f_adh[k].y = 0.0; }
			if(c2_Nbonds[k]<=0.000001) { c2_Nbonds[k] = 0.0; c2_Nodes[k].sticky = 0; c2_f_adh[k].x = 0.0; c2_f_adh[k].y = 0.0; }

			if(c1_Nbonds[k]>0.000001 && c1_Nodes[k].def.x<=c2_xFront+c2_delta_tip && sqrt(pow(c1_f_adh[k].x,2)+pow(c1_f_adh[k].y,2))/c1_Nbonds[k]>0.0*crit_adh_1){
				c1_Nbonds[k] = 0.0;
				c1_Nodes[k].sticky = 0;
				c1_f_adh[k].x = 0.0; 
				c1_f_adh[k].y = 0.0;
			} 

			if((tstep*t)>T_RELAX && (tstep*t)<T_CELL2MOVE){
				if(c1_Nbonds[k]>0.000001 && sqrt(pow(c1_f_adh[k].x,2)+pow(c1_f_adh[k].y,2))/c1_Nbonds[k]>crit_adh_1){
					c1_Nbonds[k] = 0.0;
					c1_Nodes[k].sticky = 0;
					c1_f_adh[k].x = 0.0; 
					c1_f_adh[k].y = 0.0;
				} 
			} else if((tstep*t)>=T_CELL2MOVE){
				if(c1_Nbonds[k]>0.000001 && sqrt(pow(c1_f_adh[k].x,2)+pow(c1_f_adh[k].y,2))/c1_Nbonds[k]>crit_adh_1){
					c1_Nbonds[k] = 0.0;
					c1_Nodes[k].sticky = 0;
					c1_f_adh[k].x = 0.0; 
					c1_f_adh[k].y = 0.0;
				} 

				if(c2_Nbonds[k]>0.000001 && sqrt(pow(c2_f_adh[k].x,2)+pow(c2_f_adh[k].y,2))/c2_Nbonds[k]>crit_adh_2){
					c2_Nbonds[k] = 0.0;
					c2_Nodes[k].sticky = 0;
					c2_f_adh[k].x = 0.0; 
					c2_f_adh[k].y = 0.0;
				} 
			}

			// compute total repulsion forces
			c1_sum_junction += sqrt(pow(c1_f_junction[k].x, 2)+pow(c1_f_junction[k].y, 2));
			c2_sum_junction += sqrt(pow(c2_f_junction[k].x, 2)+pow(c2_f_junction[k].y, 2));

			// compute total wall force densities
			c1_f_wall[k].x += -c1_f_contact[k].x-c1_f_adh[k].x;
      c1_f_wall[k].y += -c1_f_contact[k].y-c1_f_adh[k].y;

      c2_f_wall[k].x += -c2_f_contact[k].x-c2_f_adh[k].x;
      c2_f_wall[k].y += -c2_f_contact[k].y-c2_f_adh[k].y;

			// compute total forces (f labels forces not force densities) 
			c1_f[k].x = (c1_ef[k].x + c1_f_poly[k].x + c1_f_contact[k].x + c1_f_adh[k].x + c1_f_correctPolymerization[k] + c1_f_junction[k].x)*c1_Nodes[k].characteristicL;
      c1_f[k].y = (c1_ef[k].y + c1_f_contact[k].y + c1_f_adh[k].y + c1_f_junction[k].y)*c1_Nodes[k].characteristicL;	

      c2_f[k].x = (c2_ef[k].x + c2_f_poly[k].x + c2_f_contact[k].x + c2_f_adh[k].x + c2_f_correctPolymerization[k] + c2_f_junction[k].x)*c2_Nodes[k].characteristicL;
      c2_f[k].y = (c2_ef[k].y + c2_f_contact[k].y + c2_f_adh[k].y + c2_f_junction[k].y)*c2_Nodes[k].characteristicL;	
		}

		// COMPUTE PRESSURE CONTRIBUTION
		c1_p0 = computeP0Term(Npts, c1_Nodes, c1_def_dist_L, c1_def_dist_R, c1_n, c1_f); 
		c1_currentA = computeArea(Npts, c1_Nodes);
  	c1_p = c1_p0 + k_cell_1*log(c1_A0/c1_currentA);

  	c2_p0 = computeP0Term(Npts, c2_Nodes, c2_def_dist_L, c2_def_dist_R, c2_n, c2_f); 
		c2_currentA = computeArea(Npts, c2_Nodes);
  	c2_p = c2_p0 + k_cell_2*log(c2_A0/c2_currentA);

		computePressureForces(Npts, c1_Nodes, c1_n, c1_pf, c1_p);
		computePressureForces(Npts, c2_Nodes, c2_n, c2_pf, c2_p);

		for(int k=0; k<Npts; k++) {
      // careful: pf are force densities (need to convert them to forces)
      c1_f[k].x += c1_pf[k].x*c1_Nodes[k].characteristicL;
      c1_f[k].y += c1_pf[k].y*c1_Nodes[k].characteristicL;
    
			c1_prev_pfx[k] = c1_pf[k].x*c1_Nodes[k].characteristicL;
			c1_prev_pfy[k] = c1_pf[k].y*c1_Nodes[k].characteristicL;

			c2_f[k].x += c2_pf[k].x*c2_Nodes[k].characteristicL;
      c2_f[k].y += c2_pf[k].y*c2_Nodes[k].characteristicL;
    
			c2_prev_pfx[k] = c2_pf[k].x*c2_Nodes[k].characteristicL;
			c2_prev_pfy[k] = c2_pf[k].y*c2_Nodes[k].characteristicL;

			c1_prev_repx[k] = c1_f_junction[k].x*c1_Nodes[k].characteristicL;
			c1_prev_repy[k] = c1_f_junction[k].y*c1_Nodes[k].characteristicL;

			c2_prev_repx[k] = c2_f_junction[k].x*c2_Nodes[k].characteristicL;
			c2_prev_repy[k] = c2_f_junction[k].y*c2_Nodes[k].characteristicL;

			c1_sum_prev_junction += sqrt( pow(c1_prev_repx[k], 2) + pow(c1_prev_repy[k], 2) );
			if(c2_Nodes[k].def.y<=3.0){
				c2_sum_prev_junction += sqrt( pow(c2_prev_repx[k], 2) + pow(c2_prev_repy[k], 2) );
			}
		}

		// MOVE MEMBRANE POINTS
		for(int k=0; k<Npts; k++) {
			c1_Nodes[k].def.x = c1_Nodes[k].def.x + tstep*c1_f[k].x/(xi_1*c1_Nodes[k].characteristicL);
			c1_Nodes[k].def.y = c1_Nodes[k].def.y + tstep*c1_f[k].y/(xi_1*c1_Nodes[k].characteristicL);
   		c1_v[k].x = c1_f[k].x/(xi_1*c1_Nodes[k].characteristicL);
			c1_v[k].y = c1_f[k].y/(xi_1*c1_Nodes[k].characteristicL);

			c2_Nodes[k].def.x = c2_Nodes[k].def.x + tstep*c2_f[k].x/(xi_2*c2_Nodes[k].characteristicL);
			c2_Nodes[k].def.y = c2_Nodes[k].def.y + tstep*c2_f[k].y/(xi_2*c2_Nodes[k].characteristicL);
   		c2_v[k].x = c2_f[k].x/(xi_2*c2_Nodes[k].characteristicL);
			c2_v[k].y = c2_f[k].y/(xi_2*c2_Nodes[k].characteristicL);	

			c1_prev_vx[k] = c1_v[k].x;
			c1_prev_vy[k] = c1_v[k].y;
			
			c2_prev_vx[k] = c2_v[k].x;
			c2_prev_vy[k] = c2_v[k].y;
		}

		if(t%10000==0 && t!=0) {
			printf("%.3f, Cell 1: c1_sum_fpoly=%.5f, c1_sum_junction=%.5f, A=%.5f, Cell 2: c2_sum_fpoly=%.5f, c2_sum_junction=%.5f, A=%.5f\n",t*tstep,c1_sum_fpoly,c1_sum_junction,c1_currentA,c2_sum_fpoly,c2_sum_junction,c2_currentA);
			for(int k=0; k<Npts; k++) {
				f1 << t*tstep << " " << c1_Nodes[k].def.x << " " << c1_Nodes[k].def.y << " " << c1_Nodes[k].ref.x << " " << c1_Nodes[k].ref.y << " " << c1_ef[k].x << " " << c1_ef[k].y << " " << c1_f_contact[k].x << " " << c1_f_contact[k].y << " " << c1_f_adh[k].x << " " << c1_f_adh[k].y << " " << c1_pf[k].x << " " << c1_pf[k].y << " " << c1_f[k].x << " " << c1_f[k].y << " " << c1_Nodes[k].characteristicL << " " << c1_f_correctPolymerization[k] << " " << c1_Nodes[k].flag << " " << c1_Nodes[k].sticky << " " << c1_Nbonds[k] << " " << c1_kminusTerm[k] << " " << c1_yTerm[k] << " " << c1_Nodes[k].front << " " << c1_v[k].x << " " << c1_v[k].y << " " << c1_f_poly[k].x << " " << c1_f_poly[k].y << " " << c1_load_force << " " << c1_f_junction[k].x << " " << c1_f_junction[k].y << " " << c1_f_repulsion[k].x << " " << c1_f_repulsion[k].y << endl;  
				f2 << t*tstep << " " << c1_Nodes[k].adh_link.x << " " << c1_Nodes[k].adh_link.y << " "<< c1_f_wall[k].x << " " << c1_f_wall[k].y << " " << c1_Nodes[k].sticky << endl;	

				f3 << t*tstep << " " << c2_Nodes[k].def.x << " " << c2_Nodes[k].def.y << " " << c2_Nodes[k].ref.x << " " << c2_Nodes[k].ref.y << " " << c2_ef[k].x << " " << c2_ef[k].y << " " << c2_f_contact[k].x << " " << c2_f_contact[k].y << " " << c2_f_adh[k].x << " " << c2_f_adh[k].y << " " << c2_pf[k].x << " " << c2_pf[k].y << " " << c2_f[k].x << " " << c2_f[k].y << " " << c2_Nodes[k].characteristicL << " " << c2_f_correctPolymerization[k] << " " << c2_Nodes[k].flag << " " << c2_Nodes[k].sticky << " " << c2_Nbonds[k] << " " << c2_kminusTerm[k] << " " << c2_yTerm[k] << " " << c2_Nodes[k].front << " " << c2_v[k].x << " " << c2_v[k].y << " " << c2_f_poly[k].x << " " << c2_f_poly[k].y << " " << c2_load_force << " " << c2_f_junction[k].x << " " << c2_f_junction[k].y << " " << c2_f_repulsion[k].x << " " << c2_f_repulsion[k].y << endl;  
				f4 << t*tstep << " " << c2_Nodes[k].adh_link.x << " " << c2_Nodes[k].adh_link.y << " "<< c2_f_wall[k].x << " " << c2_f_wall[k].y << " " << c2_Nodes[k].sticky << endl;	
			}
		}
	}	

	f1.close();
  f2.close();
  f3.close();
  f4.close();
}

// DRAW START CONFIGURATION
void startCurve() {
  
  // file handling
  ifstream f1;
  ifstream f3;

  // CELL 1 DATA
	f1.open("Configurations/C1_MembraneN162_pre-relax_1_0.1_0.14_0.01_10_spring_posrelax_d_0.5_10.txt");
	cout << "cell 1 pre-relaxed profile opened." << endl;
	// determine length
  int Npoints = -1; 

  string b1;
  while( !f1.eof() ) { 
		getline(f1, b1); 
		Npoints++;
  }
  f1.close();

  cout << "Npoints = " << Npoints << endl;

  vertex c1_Nodes[Npoints];
  vertex c2_Nodes[Npoints];

  f1.open("Configurations/C1_MembraneN162_pre-relax_1_0.1_0.14_0.01_10_spring_posrelax_d_0.5_10.txt");
  // GET CELL 1 DATA
  int counter = 0;
	int M1 = 0;  // number of points on surface
  double c1, c2, c3, c4;
	while(f1 >> c1 >> c2 >> c3 >> c4) {
		c1_Nodes[counter].init.x = c1; 
    c1_Nodes[counter].init.y = c2; 
    c1_Nodes[counter].ref.x = c1; 
    c1_Nodes[counter].ref.y = c2; 
    c1_Nodes[counter].def.x = c3; 
    c1_Nodes[counter].def.y = c4; 
    c1_Nodes[counter].exterior = 0; 	

		// adhesion and surface information
		c1_Nodes[counter].adh_link.x = std::numeric_limits<double>::quiet_NaN();
		c1_Nodes[counter].adh_link.y = std::numeric_limits<double>::quiet_NaN();
		c1_Nodes[counter].bottom = 0;
		c1_Nodes[counter].top = 0;
		c1_Nodes[counter].sticky = 0;
		c1_Nodes[counter].flag = 0;
		c1_Nodes[counter].k_adh = 1.0;
		c1_Nodes[counter].junction = 0;

		// check bottom and introduce bottom wall and bottom channel 
		if(c1_Nodes[counter].init.y<0.5) {
			c1_Nodes[counter].bottom = 1;
			c1_Nodes[counter].adh_link.x = c1_Nodes[counter].ref.x;
     	c1_Nodes[counter].adh_link.y = -0.5;
			c1_Nodes[counter].sticky= 1; 
			M1++;
		}
		counter++;
  }
  f1.close();

  cout << "cell 1 counter = " << counter << endl;


  f3.open("Configurations/C2_MembraneN162_pre-relax_1_0.1_0.14_0.01_10_spring_posrelax_d_0.5_10.txt");
  // determine length
  Npoints = -1; 

  string b2;
  while( !f3.eof() ) { 
		getline(f3, b2); 
		Npoints++;
  }
  f3.close();
  cout << "Npoints = " << Npoints << endl;
	cout << "cell 2 pre-relaxed profile opened." << endl;
	f3.open("Configurations/C2_MembraneN162_pre-relax_1_0.1_0.14_0.01_10_spring_posrelax_d_0.5_10.txt");
  // GET CELL 2 DATA
  counter = 0;
	int M2 = 0;  // number of points on surface
  double c5, c6, c7, c8;
	while(f3 >> c5 >> c6 >> c7 >> c8) {
		c2_Nodes[counter].init.x = c5; 
    c2_Nodes[counter].init.y = c6; 
    c2_Nodes[counter].ref.x = c5; 
    c2_Nodes[counter].ref.y = c6; 
    c2_Nodes[counter].def.x = c7; 
    c2_Nodes[counter].def.y = c8; 
    c2_Nodes[counter].exterior = 0; 	

		// adhesion and surface information
		c2_Nodes[counter].adh_link.x = std::numeric_limits<double>::quiet_NaN();
		c2_Nodes[counter].adh_link.y = std::numeric_limits<double>::quiet_NaN();
		c2_Nodes[counter].bottom = 0;
		c2_Nodes[counter].top = 0;
		c2_Nodes[counter].sticky = 0;
		c2_Nodes[counter].flag = 0;
		c2_Nodes[counter].k_adh = 1.0;
		c2_Nodes[counter].junction = 0;

		// check bottom and introduce bottom wall and bottom channel 
		if(c2_Nodes[counter].init.y<0.5) {
			c2_Nodes[counter].bottom = 1;
			c2_Nodes[counter].adh_link.x = c2_Nodes[counter].ref.x;
     	c2_Nodes[counter].adh_link.y = -0.5;
			c2_Nodes[counter].sticky= 1; 
			M2++;
		}
		counter++;
  }
  f3.close();

  cout << "cell 2 counter = " << counter << endl;

	// progress points forward in time
  printf("%d points in each cell. Cell 1: %d points on the bottom / contact with surface. Cell 2: %d points on the bottom / contact with surface.\n",Npoints,M1,M2);
	progress(Npoints, c1_Nodes, c2_Nodes);
}

// Boilerplate: set up the filenames for the output files
void set_filenames(char *outsuffix){

  // allocate memory for output filenames
  outfile1 = new char[MAX_OUTFILE];
  outfile2 = new char[MAX_OUTFILE];
  outfile3 = new char[MAX_OUTFILE];
  outfile4 = new char[MAX_OUTFILE];

	// copy over the default filename for each file
  outfile1 = strcat(outfile1, OUTFILE1);
  outfile2 = strcat(outfile2, OUTFILE2);
  outfile3 = strcat(outfile3, OUTFILE3);
  outfile4 = strcat(outfile4, OUTFILE4);

	// append values of variables in filename
  stringstream c1_ss (stringstream::in | stringstream::out);
  c1_ss << "_" << gamma_mem_1 << "_" << k_mem_1 << "_" << adh_prop_1 << "_" << kminus_bond_1 << "_" << e_bond_1 << "_" << crit_adh_1 << "_" << swell_1 << "_" << k_junction << "_" << junction_rep << "_" << junction_rest_L << "_" << junction_dist << ".txt";
  std::string c1_str = c1_ss.str();
  outfile1 = strcat(outfile1, (char *) c1_str.c_str());

  // empty c1_ss
  c1_ss.str(std::string());
  // append all values using stringstream
  c1_ss << "_" << gamma_mem_1 << "_" << k_mem_1 << "_" << adh_prop_1 << "_" << kminus_bond_1 << "_" << e_bond_1 << "_" << crit_adh_1 << "_" << swell_1 << "_" << k_junction << "_" << junction_rep << "_" << junction_rest_L << "_" << junction_dist << ".txt";
  // convert to char* at the end
  c1_str = c1_ss.str();
  outfile2  = strcat(outfile2, (char *) c1_str.c_str());

  // append values of variables in filename
  stringstream c2_ss (stringstream::in | stringstream::out);
  c2_ss << "_" << gamma_mem_2 << "_" << k_mem_2 << "_" << adh_prop_2 << "_" << kminus_bond_2 << "_" << e_bond_2 << "_" << crit_adh_2 << "_" << swell_2 << "_" << k_junction << "_" << junction_rep << "_" << junction_rest_L << "_" << junction_dist << ".txt";
  std::string c2_str = c2_ss.str();
  outfile3 = strcat(outfile3, (char *) c2_str.c_str());

  // empty c2_ss
  c2_ss.str(std::string());
  // append all values using stringstream
  c2_ss << "_" << gamma_mem_2 << "_" << k_mem_2 << "_" << adh_prop_2 << "_" << kminus_bond_2 << "_" << e_bond_2 << "_" << crit_adh_2 << "_" << swell_2 << "_" << k_junction << "_" << junction_rep << "_" << junction_rest_L << "_" << junction_dist << ".txt";
  // convert to char* at the end
  c2_str = c2_ss.str();
  outfile4  = strcat(outfile4, (char *) c2_str.c_str());


  // printout the filename for debugging purposes
  if (dbg) {
    printf("Output1: %s\n", outfile1);
    printf("Output2: %s\n", outfile2);
    printf("Output3: %s\n", outfile3);
    printf("Output4: %s\n", outfile4);
	}
}
// END set_filenames

// Boilerplate: help
void help(){
  printf("\nUsage: dicty [-o <suffix>] [-d] -[h}\n");
  printf("options:\n");
  printf("\t-o <suffix>   Adds <suffix> string to all output filenames\n");
  printf("\t-d            Enables debugging messages\n");
  printf("\t-a <value>   Sets gamma_mem_1\n");
 	printf("\t-b <value> 	Sets k_mem_1\n");
	printf("\t-c <value>   Sets adh_prop_1\n");
	printf("\t-e <value>   Sets kminus_bond_1\n");
  printf("\t-f <value>   Sets e_bond_1\n");
  printf("\t-g <value>   Sets crit_adh_1\n");
	printf("\t-i <value>		Sets swell_1\n");
	printf("\t-j <value>   Sets gamma_mem_2\n");
 	printf("\t-k <value> 	Sets k_mem_2\n");
	printf("\t-l <value>   Sets adh_prop_2\n");
	printf("\t-m <value>   Sets kminus_bond_2\n");
  printf("\t-n <value>   Sets e_bond_2\n");
  printf("\t-p <value>   Sets crit_adh_2\n");
	printf("\t-q <value>		Sets swell_2\n");
  printf("\t-h            This message\n\n");
}
// END help

// MAIN
int main(int argc, char **argv) {
	int cmd;
  char *outsuffix = NULL;

  while ((cmd = getopt(argc, argv, "o:a:b:c:e:f:g:i:j:k:l:m:n:p:q:hd")) != -1) {
    switch(cmd) {
    case 'o':
      outsuffix = optarg;
      break;
    case 'h':
      help();
      return 0;
    case 'd':
      dbg = true;
      break;
    case 'a':
			gamma_mem_1 = atof(optarg);
      break;
		case 'b':
			k_mem_1 = atof(optarg);
			break;
		case 'c':
			adh_prop_1 = atof(optarg);
			break;
    case 'e':
      kminus_bond_1 = atof(optarg);
      break;
    case 'f':
      e_bond_1 = atof(optarg);
      break;
    case 'g':
      crit_adh_1 = atof(optarg);
      break;
    case 'i':
			swell_1 = atof(optarg);
			break;
		case 'j':
			gamma_mem_2 = atof(optarg);
      break;
		case 'k':
			k_mem_2 = atof(optarg);
			break;
		case 'l':
			adh_prop_2 = atof(optarg);
			break;
    case 'm':
      kminus_bond_2 = atof(optarg);
      break;
    case 'n':
      e_bond_2 = atof(optarg);
      break;
    case 'p':
      crit_adh_2 = atof(optarg);
      break;
    case 'q':
			swell_2 = atof(optarg);
			break;
		default:
      fprintf(stderr, "ERROR: unknown argument");
      return 1;
    }
  }
  
  // set the filenames for output
  set_filenames(outsuffix);

	printf("Two pill-shaped cells in channel simulation\n");
  clock_t begin = clock();
  startCurve();
  clock_t end = clock();
  printf("Total computation time (s): %.10f\n", double(diffclock(end,begin)));

  return 0;
}

	


