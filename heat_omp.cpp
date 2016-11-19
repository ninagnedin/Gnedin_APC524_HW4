/*
 * heat_omp.cpp
 *
 *  Created on: Nov 17, 2016
 *      Author: ninagnedin
 */

#include<iostream>
#include<fstream>
#include<stdio.h>
#include<string>
#include<stdlib.h>
#include<math.h>
#include<time.h>
using namespace std;


int main(int argc, char *argv[])
{
  clock_t timer;
  timer=clock();
  const double kappa=1;
  double avg=0;

  if(argc!=3){
      printf("Not enough inputs");
      exit(0);
  }

  const int size=atof(argv[1]);
  printf("Grid size: %c \n",size);
  const int thread_n=atof(argv[2]);

  double** grid = new double*[size];
  for (int i=0;i<size;i++){
     grid[i]=new double[size];
  }

  double** grid_n = new double*[size];
  for (int i=0;i<size;i++){
     grid_n[i]=new double[size];
  }

  double dx=M_PI/size;
  double dx2=pow(dx,2);
  const double dt=dx2/(8*kappa);
  const double time=0.5*pow(M_PI,2)/kappa;
  const double step_count=time/dt;
  //Boundary conditions
    for(int i=0;i<size;i++){
	    for(int j=0;j<size;j++){
		    grid[i][j]=0;
	    }
    }

  for (int i=0;i<size;i++){
    grid[i][0]=pow(cos(i*M_PI/double(size)),2);
    grid[i][size-1]=pow(sin(i*M_PI/double(size)),2);

  }


    for(int i=0;i<size;i++){
	    for(int j=0;j<size;j++){
		    grid_n[i][j]=grid[i][j];
	    }
    }


  for(int count=0;count<step_count;count++){
  #pragma omp parallel for num_threads(thread_n)
	  for(int i=1;i<size-1;i++){
		  for(int j=1;j<size-1;j++){
			  grid_n[i][j]=grid[i][j]+kappa*dt*(grid[i-1][j]+grid[i+1][j]+grid[i][j-1]+grid[i][j+1]-4*grid[i][j])/dx2;
		  }
	  }
  #pragma omp parallel for num_threads(thread_n)
    for(int i=1;i<size-1;i++){
	    grid_n[0][i]=grid[0][i]+kappa*dt*(grid[size-1][i]+grid[1][i]+grid[0][i-1]+grid[0][i+1]-4*grid[0][i])/dx2;
    }
#pragma omp parallel for num_threads(thread_n)
    for(int i=1;i<size-1;i++){
	    grid_n[size-1][i]=grid[size-1][i]+ kappa*dt*(grid[size-2][i]+grid[0][i]+grid[size-1][i-1]+grid[size-1][i+1]-4*grid[size-1][i])/dx2;
    }
#pragma omp parallel for num_threads(thread_n)
    for(int i=0;i<size;i++){
	    for(int j=0;j<size;j++){
		    grid[i][j]=grid_n[i][j];
	    }
    }

  }

#pragma omp parallel for num_threads(thread_n)
   for(int i=0;i<size;i++){
	    for(int j=0;j<size;j++){
		    avg+=grid[i][j];
	    }
    }

   timer=clock()-timer;
   printf("Total time taken is %f seconds \n ",float(timer)/CLOCKS_PER_SEC);
   printf("Average temperature over grid is: %f \n",avg/(size*size));

  char filename[50];
  sprintf(filename,"heatmap_omp_%d.txt",size);
  ofstream fout(filename);

    for(int i=0;i<size;i++){
	    for(int j=0;j<size;j++){
	      fout<< i<<" "<<j<<" "<< grid[i][j]<<endl;
	    }fout<<endl;
    }

    fout.close();
}




