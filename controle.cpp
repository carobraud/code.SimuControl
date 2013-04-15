#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>
#include <stdint.h>
#include <fstream>
#include <netcdfcpp.h>
#define cimg_debug         2

#include "../NetCDF.Tool/struct_parameter_NetCDF.h"
#include "../NetCDF.Tool/NetCDFinfo.h"
#include "../CImg.Tool/CImg_NetCDF.h"
#include "../CImg/CImg.h"

using namespace cimg_library;
using namespace std;


//! conversion integer to string
/**
 * \param [in] integer
 * \param [out] string
**/
string IntoStr(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

//! save in NetCDF format TODO: make class saveNetCDF(2C2D or 1C2D or 1C1D) + add time 
/**
 * \param [in] data        data 2C2D to save
 * \param [out] file_name   name of file to save
**/
template <typename Tdata>
int save_netcdf_2C2D(std::string fo, cimg_library::CImgList<Tdata>& imgList)
{
//dimension names
  vector<string> dim_names;
  string dim_time="dimt";
  dim_names.push_back("dimx");
  dim_names.push_back("dimy");
//variable names
  vector<string> var_names;
  var_names.push_back("u");
  var_names.push_back("v");
//unit names
  vector<string> unit_names;
  unit_names.push_back("m/s");
  unit_names.push_back("m/s");

// CREATE NC 
//  string fo="Test.nc";
  CImgListNetCDF<float> cimgList;
  //CImgList<Tdata> imgList(var_names.size(),nx,ny);
////file
  cout << "CImgListNetCDF::saveNetCDFFile(" << fo << ",...) return " 	<< cimgList.saveNetCDFFile((char*)fo.c_str()) << endl;
////dim
  cout << "CImgListNetCDF::addNetCDFDims(" << fo << ",...) return " 	<< cimgList.addNetCDFDims(imgList,dim_names,dim_time) << endl;
////var
  cout << "CImgListNetCDF::addNetCDFVar(" << fo << ",...) return " 	<< cimgList.addNetCDFVar(imgList,var_names,unit_names) << endl;
  //  for(int t=0;t<nt;++t)
  //  {
      cout << "CImgListNetCDF::addNetCDFData" << fo << ",...) return " 	<< cimgList.addNetCDFData(imgList) << endl;
      //  }
  cout << "*** SUCCESS writing example file " << fo << "!" << endl;
  return 0;
}

//! save in NetCDF format 
/**
 * \param [in] data        data 1C2D to save
 * \param [out] file_name   name of file to save
**/
template <typename Tdata>
int save_netcdf_1C2D(std::string fo, cimg_library::CImg<Tdata>& img)
{
//dimension names
  vector<string> dim_names;
  string dim_time="dimt";
  dim_names.push_back("dimx");
  dim_names.push_back("dimy");
//variable names
  string var_name="vortz";
//unit names
  string unit_name="m/s";


// CREATE NC 
//  string fo="Test.nc";
  CImgNetCDF<float> cimg;
  //CImgList<Tdata> imgList(var_names.size(),nx,ny);
////file
  cout << "CImgNetCDF::saveNetCDFFile(" << fo << ",...) return " 	<< cimg.saveNetCDFFile((char*)fo.c_str()) << endl;
////dim
  cout << "CImgNetCDF::addNetCDFDims(" << fo << ",...) return " 	<< cimg.addNetCDFDims(img,dim_names,dim_time) << endl;
////var
  cout << "CImgNetCDF::addNetCDFVar(" << fo << ",...) return " 	<< cimg.addNetCDFVar(img,var_name,unit_name) << endl;
  //  for(int t=0;t<nt;++t)
  //  {
      cout << "CImgNetCDF::addNetCDFData" << fo << ",...) return " 	<< cimg.addNetCDFData(img) << endl;
      //  }
  cout << "*** SUCCESS writing example file " << fo << "!" << endl;
  return 0;
}

//! GLOBAL VARIABLES IN time loop
/**/
int GcTime=0;
/**/
//! main function of the control program
/**
 * \param [in] nx,ny,u,v,x,y 
 * \param [out] nxc, xa, xb, xc, uc, vc 
 **/
extern "C" {
  void controle_ (int32_t &nx,int32_t &ny,float *x,float *y,float *u,float *v,int32_t &nxc,float *xc,float *uc,float *vc,float &xa,float &xb, int32_t *mask);}

void controle_ (int32_t &nx,int32_t &ny,float *x,float *y,float *u,float *v,int32_t &nxc,float *xc,float *uc,float *vc,float &xa,float &xb, int32_t *mask)
{
//! GLOBALE VARIABLES
/**/
  GcTime++;
  //  std::cout<< "Global variable: Compteur temps =" << GlobalCompteurTime << std::endl;
/**/
//! VARIABLES FOR CONTROL
/*
  TimeSave=;
  TimeInitialControl=;
*/
 
//! Extraction 2D vector field "data" for CImg
/*
* 
/**/
//  cout << "nx=" << nx <<"\t"<< "ny="<< ny <<endl; 
  //cimg_library::CImg<float> ux(nx,ny);
  //cimg_library::CImg<float> uy(nx,ny);
  cimg_library::CImgList<float> data(2,nx,ny);
  cimg_library::CImg<int> bump(nx,ny);
  
  int k=0;
  for (int j = 0; j<ny; j++) 
    {
      for (int i = 0; i<nx; i++) 
	{
	  data(0,i,j)=u[k];
	  data(1,i,j)=v[k];
	  bump(i,j)=mask[k];
	  k++;
	}
    }
  //data[0].display();
  //data[1].display();
  //data.print();


//! Computation 2D Vorticity field //input: 2D velocity vector field (u,v)
//!\todo: function to compute derivees and vorticity
// **** Variables 
  cimg_library::CImg<float> dvdx(nx-2,ny-2); //temp
  cimg_library::CImg<float> dudy(nx-2,ny-2); //temp
  cimg_library::CImg<float> vortz(nx-2,ny-2); //output
// **** 

  // Steps in X and Y directions
  //for (int i=0;i<ny-1;i++)
  //  {
    float dx=x[2]-x[1];
    float dy=y[2]-y[1];
    //  std::cout << "dy=" << dy << endl;
    // }

    int jcrop=0;
    for (int j=1;j<ny-1;j++)
      {
	int icrop=0; 
	for (int i=1;i<nx-1;i++)
	  {
	    //std::cout << "ICROP " << icrop << "\t JCROP "<< jcrop << endl;
	    //     std::cout << "COUCOU\tab i=" << i << endl;
	    //std::cout << "\t x=" << x << "\t y= "<< y << endl;
	    //std::cout << "dx=" << dx << "dy=" << dy << endl;
	    //std::cout << "V(i+1)=" << data[1](x+1,y) << endl;
	    // **** dvdx
	    dvdx(icrop,jcrop)=(data(1,i+1,j)-data(1,i-1,j))/(2*dx);
	    // **** dudy
	    dudy(icrop,jcrop)=(data(0,i,j+1)-data(0,i,j-1))/(2*dy);
	    // **** vortz
	    vortz(icrop,jcrop)=dvdx(icrop,jcrop)-dudy(icrop,jcrop);

	    //std::cout << "Vortz=" << vortz(icrop,jcrop) << endl;

	    icrop++;
	    //break;
	  }
	jcrop++;
      }
    
    //std::cout << ""  << endl;
    // vortz.display();
//*** End compute vorticity

//! Set control
//*** Find max bump 
//! \todo put it outside loop in time (no need to calculate it every time step)
/* 
   float ymax=0.;
    float xmax=0;
    int iymax=0;
    int ixmax=0;
      for(int j=0;j<ny;j++)
        {
	  for (int i=0;i<nx;i++)
	    {
	      if (bump(i,j) == 0) 
		{
		  if (y[j] >= ymax) 
		    {
		      ymax=y[j];		      
		      xmax=x[i];
		      iymax=j;
		      ixmax=i;
		      // std::cerr << "Summit x=" << xmax << "\t y=" << ymax << "\t Mask=" << bump(i,j) << endl;
		    }
		}
	    }
	}				       
      //std::cout << "Indice Summit i=" << ixmax << "\t j=" << iymax << endl;
      //std::cout << "Coordinate Summit x=" << xmax << "\t y=" << ymax << endl;
      
      //exit(0);
*/
//*** End Find max bump 

//*** Apply control relatively to submit (step fonction)
    //int ixa=ixmax;
    int ixa=50;
    int ixb=ixa+5;
    xa=x[ixa];
    xb=x[ixb];
      
    //    std::cout << "COUCOU" << endl;
    nxc=5;
    float ampl=0.1;
    //std::cout << "ixa=" << ixa << "\t xa=" << xa << endl;
      //std::cout << "ixb=" << ixb << "\t xb=" << xb << endl;
      //std::cout << "x(135)=" << x[135] << "\t x(140)=" << x[140] << endl;
      
      //!necessary to go back to fortran program
      float * temp;
      temp= new float [nxc];
      for(int i=ixa; i<ixb+1;i++)
	{
	  temp[i]=ampl;
	}
      //***************
      
      int TimeInitialControlUp=200;
      int TimeInitialControlDown=1500;
      if ((GcTime >= TimeInitialControlUp) && (GcTime <= TimeInitialControlDown))
	{
	  
	  for(int i=ixa; i<ixb+1;i++)
	    {
	      //! return to fortran program
	      *uc=temp[i];
	      *vc=temp[i];
	      *xc=x[i];
	      /*std::cout << "xc=" << *xc << endl;
		std::cout << "uc=" << *uc << endl;
		std::cout << "vc=" << *vc << endl;*/
	      xc++;
	      uc++;
	      vc++;
	    }
	}
      else
      	{
	  for(int i=ixa; i<ixb+1;i++)
	    {
	      *uc=0;
	      *vc=0;
	      *xc=x[i];
	      //std::cout << "x=" << x[i] << endl;
	      xc++;
	      uc++;
	      vc++;
	    }
	}

      delete[]temp;
 
	  /* Research of y at the bump (not necessary ...)	  		
		int j=0;
		
		for(int i=ixmax; i<ixmax+5;i++)
		{
		while(bump(i,j) == 0) 
		  {		    
		    //  std::cerr << "xmax=" << i << endl;
		    // std::cerr << "yborder1=" << j << endl;
		    // std::cout << "bump=" << bump(i,j) << endl;
		    j++;
		  } 

		//std::cerr << "yborder2=" << j << endl;
		data(0,ixmax+i,j)=ampl; 
		data(1,ixmax+i,j)=ampl; 
		}*/


	
      //      exit(0);

//! Find max vorticity

//! Extract horizontal line 

//! Locate max/min U or V at the vortex location AND at the bump location


//! Save 2D field in NetCDF file for each "TimeSave" step (function of "GlobalCompteurTime")
/*
* 
*
*/
  int TimeSave=100;
  if (GcTime%TimeSave == 0)
    {
      //std::cout <<" TIME TO SAVE!" << GlobalCompteurTime << std::endl;
      std::string number=IntoStr(GcTime);
      //cout << number <<endl;
      std::string file_uv="output/uv_"+number+".nc";
      std::string file_vortz="output/vortz_"+number+".nc";
      //cout << file_name <<endl;
      //std::string file_name="TEST.nc";
      save_netcdf_2C2D(file_uv,data); 
      //      vortz.print("VORTZ");
      save_netcdf_1C2D(file_vortz,vortz); 
    }
  //   exit(0);
 
}

