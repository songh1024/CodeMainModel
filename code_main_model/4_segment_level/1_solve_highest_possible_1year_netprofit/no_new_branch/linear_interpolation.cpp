#include "parameters.h" 

/*double lin_interpo_2(int nx, int ny, double* gridx, double* gridy, double** wfun, double x, double y)
{
	double valuex1=0.0; //py
	double valuex2=0.0; //py+1
	double value=0.0;

	double widthx=gridx[1]-gridx[0]; // width of grid
	double widthy=gridy[1]-gridy[0]; // width of grid

	double minx=gridx[0];
	double maxx=gridx[nx-1];
	double miny=gridy[0];
	double maxy=gridy[ny-1];

	int position_x=(int)floor((x-minx)/widthx); // position_x is gonna be the left point
	int position_y=(int)floor((y-miny)/widthy); // position_y is gonna be the left point

	position_x=(position_x<=-1)? 0: position_x;
	position_x=(position_x>=nx-1)? (nx-2): position_x;
	position_y=(position_y<=-1)? 0: position_y;
	position_y=(position_y>=ny-1)? (ny-2): position_y;

	double weightx1=(x-gridx[position_x])/widthx;
	double weightx2=1-weightx1;

	valuex1=weightx1*wfun[position_x+1][position_y]+weightx2*wfun[position_x][position_y];
	valuex2=weightx1*wfun[position_x+1][position_y+1]+weightx2*wfun[position_x][position_y+1];

	double weighty1=(y-gridy[position_y])/widthy;
	double weighty2=1-weighty1;

	value=weighty1*valuex2+weighty2*valuex1;
    
	return value;
}*/


double lin_interpo_2_most_robust(int nx, int ny, double* gridx, double* gridy, double** vfun, double x, double y)
{
	double value_1=0.0;
	double value_2=0.0;
	double value=0.0;

	double minx=gridx[0];
	double miny=gridy[0];
	double maxx=gridx[nx-1];
	double maxy=gridy[ny-1];

	int position_x=0;
	int position_y=0;

	if (x>=maxx)
		position_x=nx-1;
	else if (x<minx)
		position_x=-1;

	if (y>=maxy)
		position_y=ny-1;
	else if (y<miny)
		position_y=-1;

	if (position_x==0)
	{
		while (gridx[position_x+1]<=x)
			++position_x;
	}

	if (position_y==0)
	{
		while (gridy[position_y+1]<=y)
			++position_y;
	}

	if ((position_x<=-1) && (position_y<=-1))
	{
		value_1=vfun[0][0]+(vfun[1][0]-vfun[0][0])*(x-gridx[0])/(gridx[1]-gridx[0]);
		value_2=vfun[0][1]+(vfun[1][1]-vfun[0][1])*(x-gridx[0])/(gridx[1]-gridx[0]);
		value=value_1+(value_2-value_1)*(y-gridy[0])/(gridy[1]-gridy[0]);
	}

	else if ((position_x<=-1) && (position_y<(ny-1)) && (position_y>-1))
	{
		value_1=vfun[0][position_y]+(vfun[1][position_y]-vfun[0][position_y])*(x-gridx[0])/(gridx[1]-gridx[0]);
		value_2=vfun[0][position_y+1]+(vfun[1][position_y+1]-vfun[0][position_y+1])*(x-gridx[0])/(gridx[1]-gridx[0]);
		value=value_1+(value_2-value_1)*(y-gridy[position_y])/(gridy[position_y+1]-gridy[position_y]);
	}

	else if ((position_x<=-1) && (position_y>=(ny-1)))
	{
		value_1=vfun[0][ny-2]+(vfun[1][ny-2]-vfun[0][ny-2])*(x-gridx[0])/(gridx[1]-gridx[0]);
		value_2=vfun[0][ny-1]+(vfun[1][ny-1]-vfun[0][ny-1])*(x-gridx[0])/(gridx[1]-gridx[0]);
		value=value_2+(value_2-value_1)*(y-gridy[ny-1])/(gridy[ny-1]-gridy[ny-2]); 
	}
	else if ((position_x<(nx-1)) && (position_x>-1) && (position_y<=-1))
	{
		value_1=vfun[position_x][0]+(vfun[position_x+1][0]-vfun[position_x][0])*(x-gridx[position_x])/(gridx[position_x+1]-gridx[position_x]);
		value_2=vfun[position_x][1]+(vfun[position_x+1][1]-vfun[position_x][1])*(x-gridx[position_x])/(gridx[position_x+1]-gridx[position_x]);
		value=value_1+(value_2-value_1)*(y-gridy[0])/(gridy[1]-gridy[0]);
	}

	else if ((position_x<(nx-1)) && (position_x>-1) && (position_y<(ny-1)) && (position_y>-1))
	{
		value_1=vfun[position_x][position_y]+(vfun[position_x+1][position_y]-vfun[position_x][position_y])*(x-gridx[position_x])/(gridx[position_x+1]-gridx[position_x]);
		value_2=vfun[position_x][position_y+1]+(vfun[position_x+1][position_y+1]-vfun[position_x][position_y+1])*(x-gridx[position_x])/(gridx[position_x+1]-gridx[position_x]);
		value=value_1+(value_2-value_1)*(y-gridy[position_y])/(gridy[position_y+1]-gridy[position_y]);
	}

	else if ((position_x<(nx-1)) && (position_x>-1) && (position_y>=(ny-1)))
	{
		value_1=vfun[position_x][ny-2]+(vfun[position_x+1][ny-2]-vfun[position_x][ny-2])*(x-gridx[position_x])/(gridx[position_x+1]-gridx[position_x]);
		value_2=vfun[position_x][ny-1]+(vfun[position_x+1][ny-1]-vfun[position_x][ny-1])*(x-gridx[position_x])/(gridx[position_x+1]-gridx[position_x]);
		value=value_2+(value_2-value_1)*(y-gridy[ny-1])/(gridy[ny-1]-gridy[ny-2]);
	}

	else if ((position_x>=(nx-1)) && (position_y<=-1))
	{
		value_1=vfun[nx-1][0]+(vfun[nx-1][0]-vfun[nx-2][0])*(x-gridx[nx-1])/(gridx[nx-1]-gridx[nx-2]);
		value_2=vfun[nx-1][1]+(vfun[nx-1][1]-vfun[nx-2][1])*(x-gridx[nx-1])/(gridx[nx-1]-gridx[nx-2]);
		value=value_1+(value_2-value_1)*(y-gridy[0])/(gridy[1]-gridy[0]);
	}

	else if ((position_x>=(nx-1)) && (position_y<(ny-1)) && (position_y>-1))
	{
		value_1=vfun[nx-1][position_y]+(vfun[nx-1][position_y]-vfun[nx-2][position_y])*(x-gridx[nx-1])/(gridx[nx-1]-gridx[nx-2]);
		value_2=vfun[nx-1][position_y+1]+(vfun[nx-1][position_y+1]-vfun[nx-2][position_y+1])*(x-gridx[nx-1])/(gridx[nx-1]-gridx[nx-2]);
		value=value_1+(value_2-value_1)*(y-gridy[position_y])/(gridy[position_y+1]-gridy[position_y]);
	}

	else if ((position_x>=(nx-1)) && (position_y>=(ny-1)))
	{
		value_1=vfun[nx-1][ny-2]+(vfun[nx-1][ny-2]-vfun[nx-2][ny-2])*(x-gridx[nx-1])/(gridx[nx-1]-gridx[nx-2]);
		value_2=vfun[nx-1][ny-1]+(vfun[nx-1][ny-1]-vfun[nx-2][ny-1])*(x-gridx[nx-1])/(gridx[nx-1]-gridx[nx-2]);
		value=value_2+(value_2-value_1)*(y-gridy[ny-1])/(gridy[ny-1]-gridy[ny-2]);
	}

	return value;
}


double lin_interpo_2(int nx, int ny, double* gridx, double* gridy, double** vfun, double x, double y)
{
	double value_1=0.0;
	double value_2=0.0;
	double value=0.0;

	double widthx=gridx[1]-gridx[0];
	double widthy=gridy[1]-gridy[0];

	double minx=gridx[0];
	double miny=gridy[0];
	double maxx=gridx[nx-1];
	double maxy=gridy[ny-1];

	int position_x=(int)floor((x-minx)/widthx);
	int position_y=(int)floor((y-miny)/widthy);

	if ((position_x<=-1) && (position_y<=-1))
	{
		value_1=vfun[0][0]+(vfun[1][0]-vfun[0][0])*(x-gridx[0])/widthx;
		value_2=vfun[0][1]+(vfun[1][1]-vfun[0][1])*(x-gridx[0])/widthx;
		value=value_1+(value_2-value_1)*(y-gridy[0])/widthy;
	}

	else if ((position_x<=-1) && (position_y<(ny-1)) && (position_y>-1))
	{
		value_1=vfun[0][position_y]+(vfun[1][position_y]-vfun[0][position_y])*(x-gridx[0])/widthx;
		value_2=vfun[0][position_y+1]+(vfun[1][position_y+1]-vfun[0][position_y+1])*(x-gridx[0])/widthx;
		value=value_1+(value_2-value_1)*(y-gridy[position_y])/widthy;
	}

	else if ((position_x<=-1) && (position_y>=(ny-1)))
	{
		value_1=vfun[0][ny-2]+(vfun[1][ny-2]-vfun[0][ny-2])*(x-gridx[0])/widthx;
		value_2=vfun[0][ny-1]+(vfun[1][ny-1]-vfun[0][ny-1])*(x-gridx[0])/widthx;
		value=value_2+((value_2-value_1)*(y-gridy[ny-1])/widthy); 
	}
	else if ((position_x<(nx-1)) && (position_x>-1) && (position_y<=-1))
	{
		value_1=vfun[position_x][0]+(vfun[position_x+1][0]-vfun[position_x][0])*(x-gridx[position_x])/widthx;
		value_2=vfun[position_x][1]+(vfun[position_x+1][1]-vfun[position_x][1])*(x-gridx[position_x])/widthx;
		value=value_1+(value_2-value_1)*(y-gridy[0])/widthy;
	}

	else if ((position_x<(nx-1)) && (position_x>-1) && (position_y<(ny-1)) && (position_y>-1))
	{
		value_1=vfun[position_x][position_y]+(vfun[position_x+1][position_y]-vfun[position_x][position_y])*(x-gridx[position_x])/widthx;
		value_2=vfun[position_x][position_y+1]+(vfun[position_x+1][position_y+1]-vfun[position_x][position_y+1])*(x-gridx[position_x])/widthx;
		value=value_1+(value_2-value_1)*(y-gridy[position_y])/widthy;
	}

	else if ((position_x<(nx-1)) && (position_x>-1) && (position_y>=(ny-1)))
	{
		value_1=vfun[position_x][ny-2]+(vfun[position_x+1][ny-2]-vfun[position_x][ny-2])*(x-gridx[position_x])/widthx;
		value_2=vfun[position_x][ny-1]+(vfun[position_x+1][ny-1]-vfun[position_x][ny-1])*(x-gridx[position_x])/widthx;
		value=value_2+((value_2-value_1)*(y-gridy[ny-1])/widthy);
	}

	else if ((position_x>=(nx-1)) && (position_y<=-1))
	{
		value_1=vfun[nx-1][0]+((vfun[nx-1][0]-vfun[nx-2][0])*(x-gridx[nx-1])/widthx);
		value_2=vfun[nx-1][1]+((vfun[nx-1][1]-vfun[nx-2][1])*(x-gridx[nx-1])/widthx);
		value=value_1+(value_2-value_1)*(y-gridy[0])/widthy;
	}

	else if ((position_x>=(nx-1)) && (position_y<(ny-1)) && (position_y>-1))
	{
		value_1=vfun[nx-1][position_y]+((vfun[nx-1][position_y]-vfun[nx-2][position_y])*(x-gridx[nx-1])/widthx);
		value_2=vfun[nx-1][position_y+1]+((vfun[nx-1][position_y+1]-vfun[nx-2][position_y+1])*(x-gridx[nx-1])/widthx);
		value=value_1+(value_2-value_1)*(y-gridy[position_y])/widthy;
	}

	else if ((position_x>=(nx-1)) && (position_y>=(ny-1)))
	{
		value_1=vfun[nx-1][ny-2]+((vfun[nx-1][ny-2]-vfun[nx-2][ny-2])*(x-gridx[nx-1])/widthx);
		value_2=vfun[nx-1][ny-1]+((vfun[nx-1][ny-1]-vfun[nx-2][ny-1])*(x-gridx[nx-1])/widthx);
		value=value_2+((value_2-value_1)*(y-gridy[ny-1])/widthy);
	}

	return value;

}

double lin_interpo_1(int nx, double* gridx, double* v, double x)   //linear interpolation using the value function
{
	double width=gridx[1]-gridx[0];
	if (x>=gridx[nx-1])
		return v[nx-1];
		//return v[nx-1]+(v[nx-1]-v[nx-2])/width*(x-gridx[nx-1]);
	else if (x<=gridx[0])
		return v[0];
		//return v[0]+(v[1]-v[0])/width*(x-gridx[0]);
	else
	{
		int x_index=0;
		x_index=(int)floor((x-gridx[0])/width);
		int low=x_index; 
		int high=x_index+1;
		
		double q=(gridx[high]-x)/(gridx[high]-gridx[low]);
		return v[low]*q+v[high]*(1-q);

	}
}