#ifndef STRUCTS_H
#define STRUCTS_H

// VECTOR STRUCT
typedef struct vector2{
  double x; // x-component
  double y; // y-component
} vector2;

// MEMBRANE/CORTEX VERTEX STRUCT
typedef struct vertex{
  vector2 init; 						// initial reference coords
  vector2 ref; 						// reference coords
  vector2 def; 						// deformation coords 
  vector2 force; 					// force
  int exterior; 					// 1 if this is a boundary point and 0 if this is an interior point
	double characteristicL; // characteristic length unit
	int bottom; 						// 1 if this is a 'bottom' point and 0 otherwise
	int top; 								// 1 if this a 'top' point and 0 otherwise
	vector2 adh_link;
	int sticky;							// 1 for an active adhesion bond and 0 otherwise
	int flag; 							// flag for whatever is needed
	int front; 							// flag for front of the cell
	double k_adh; 					// each vertex has its own adhesion constant (to do strain stiffening) 
	int junction;           // 1 if this is a point at cell-cell junction and 0 otherwise
	int back;               // flag for back of the cell
	int pt_rep;             // 1 if this is a point at the cell-cell junction experiencing repulsion forces and 0 otherwise
} vertex;
#endif
