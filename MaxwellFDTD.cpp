// Solving Maxwell's equations using Finite Difference Time Domain FDTD
// Joshua D. John
// November 2021
//
// Compile it with:
//   g++ -o MaxwellFDTD MaxwellFDTD.cpp -lboost_iostreams -lboost_system -lboost_filesystem

#include <iostream>
#include <fstream>
//#include <string>
#include <iomanip>
#include <map>
#include <vector>
#include <cmath>

#include "gnuplot-iostream.h"

class Maxwell
{
	private:
	int NX{50}; // number of steps in x-direction
	int NY{50}; // number of steps in y-direction
	int NSize = NX*NY;
	int NSteps{100}; // number of time steps

	int NXMiddle = (int)NX/2.0; // middle point; location of pulse
	int NYMiddle = (int)NY/2.0;
	
	float* XGrid;
	float* YGrid;

	float DeltaX{0.01}; // cell size
	float DeltaT = (float)DeltaX/6e8;
	float Epsilon{8.8e-12};  
	float TZero = 20.0; // center of incident pulse
	float Width = 6.0; // width of the incident pulse
	float Time{0.0};
	
	float* E_z;
	float* H_x;
	float* H_y;
	
	public:
	Maxwell();
	~Maxwell();
	void IterateCPU();
	void IterateHetero();
	void Save(std::string fname);
	void Plot();
	
	
};

Maxwell::Maxwell()
{
	this->E_z = new float[this->NSize];
	this->H_x = new float[this->NSize];
	this->H_y = new float[this->NSize];
	
	this->XGrid = new float[this->NX];
	this->YGrid = new float[this->NY];
	
	// Initialize field to zero at all points
	for (size_t i{0}; i<this->NSize; i++)
	{
		this->E_z[i] = 0.f;
		this->H_x[i] = 0.f;
		this->H_y[i] = 0.f;	
	}
	
}

Maxwell::~Maxwell()
{
	
}

void Maxwell::IterateCPU()
{
	for (size_t t{0}; t < this->NSteps; t++)
	{
		this->Time++;
		std::cout << "Iteration step: " << t << std::endl;
		
		// calculate E_z field t
		for (size_t x{1}; x < this->NX; x++)
		{
			for (size_t y{1}; y < this->NY; y++)
			{
				int xy = (x*this->NX) + y;
				int xmy = ((x-1)*this->NX) + y;
				int xym = (x*this->NX) + (y-1);
				E_z[xy] = E_z[xy] + 0.5*(H_y[xy] - H_y[xmy] - H_x[xy] + H_x[xym]);
			}
		}
		
		// Gaussian pulse in the middle
		float inex = (float)(this->TZero - this->Time)/this->Width;
		float pulse = exp(-0.5*(inex*inex)); // Gaussian pulse
		int mid = this->NXMiddle*this->NX + this->NYMiddle;
		E_z[mid] = pulse;
		
		/*int low = 0.5*this->NXMiddle*this->NX + 0.5*this->NYMiddle;
		int up = 1.5*this->NXMiddle*this->NX + 1.5*this->NYMiddle;
		E_z[low] = pulse;
		E_z[up] = pulse;*/

		// Calculate H_x field
		for (size_t x{0}; x < this->NX-1; x++)
		{
			for (size_t y{0}; y < this->NY-1; y++)
			{
				int xy = (x*this->NX) + y;
				int xyp = (x*this->NX) + (y+1);
				H_x[xy] = H_x[xy] - 0.5*(E_z[xyp] - E_z[xy]);
			}
		}
		
		// Calculate H_y field
		for (size_t x{0}; x < this->NX-1; x++)
		{
			for (size_t y{0}; y < this->NY-1; y++)
			{
				int xy = (x*this->NX) + y;
				int xpy = ((x+1)*this->NX) + y;
				//int xyp = (x*this->NX) + (y+1);
				H_y[xy] = H_y[xy] + 0.5*(E_z[xpy] - E_z[xy]);
			}
		}
		
		std::string filename = "Iteration-" + std::to_string(t);
		this->Save(filename);
	}
}

void Maxwell::IterateHetero()
{
	
}

void Maxwell::Save(std::string fname)
{
	std::ofstream sfile (fname);
  	if (sfile.is_open())
  	{
  		for (size_t x {0}; x < this->NX; x++)
		{
			for (size_t y {0}; y < this->NY; y++)
			{
				int xy = (x*this->NX) + y;
				//this->surface.push_back(0.1f);
			
				float z = this->E_z[xy];
				sfile << x << "," << y << "," << z << "\n";
			}
		}
    		sfile.close();
	}
  	else std::cout << "Unable to open file";
}

void Maxwell::Plot()
{
	Gnuplot gp;
	std::vector<std::tuple<double, double, double>> surfaceGridPnts;
	
	for (size_t x {0}; x < this->NX; x++)
	{
		for (size_t y {0}; y < this->NY; y++)
		{
			int xy = (x*this->NX) + y;
			//this->surface.push_back(0.1f);
			
			float z = this->E_z[xy];
			surfaceGridPnts.push_back(std::make_tuple(x, y, z));
		}
	}
	gp << "set xrange [0:60]\nset yrange [0:60]\n";
	//gp << "set hidden3d\n";
	gp << "set pm3d interpolate 1,1\n";
	gp << "set dgrid3d 100,100 qnorm 2\n";
	gp << "set contour base\n";
	gp << "set cntrlabel  format '%8.3g' font ',7' start 5 interval 20\n";
	gp << "set cntrparam order 8\n";
	gp << "set cntrparam bspline\n";
	gp << "set style data lines\n";
	// gp << "set title 'contour of Sinc function'\n";
	// Data will be sent via a temporary file.  These are erased when you call
	// gp.clearTmpfiles() or when gp goes out of scope.  If you pass a filename
	// (e.g. "gp.file1d(pts, 'mydata.dat')"), then the named file will be created
	// and won't be deleted (this is useful when creating a script).
	// gp << "splot" << gp.file1d(xyz_pts_A) << "with lines title 'Gauss'," << std::endl;
	gp << "splot" << gp.file1d(surfaceGridPnts) << "with pm3d title 'Gauss'," << std::endl;
	

}

int main()
{
	Maxwell mx{};
	mx.IterateCPU();
	mx.Plot();
	return 0;
}
