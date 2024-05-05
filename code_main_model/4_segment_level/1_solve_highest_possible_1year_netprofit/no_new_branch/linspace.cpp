#include "parameters.h"

void linspace(double* grid, double min, double max, int length)
{
	
	double width=(max-min)/(length-1);
	
	for(int i=0; i<length; ++i){
		grid[i]=min+width*i;
	
	}
}

void nonlinspace(double* grid, double min, double middle, double max, int length, int length_left)
{
	
	double width_left=(middle-min)/(length_left-1);
	double width_right=(max-middle)/(length-length_left);
	
	for(int i=0; i<length_left; ++i){
		grid[i]=min+width_left*i;
	
	}
	for(int i=length_left; i<length; ++i){
		grid[i]=middle+width_right*(i-length_left+1);
	
	}
}
