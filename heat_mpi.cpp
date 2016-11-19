/*
 * heat_mpi.cpp
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
#include<mpi.h>
using namespace std;


int main(int argc, char *argv[])
{
  clock_t timer;
  timer=clock();
  const double kappa=1;
  double avg=0;

  if(argc!=2){
      printf("Not enough inputs");
      exit(0);
  }

  const int size=atof(argv[1]);
  printf("Grid size: %c \n",size);

  int n_pro;
  int rank;
  MPI_Status Stat[4];
  MPI_Request request[4];
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &n_pro);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int snp=size/n_pro;
  double** grid = new double*[snp+2];
  for (int i=0;i<size;i++){
      grid[i]=new double[size];
  }
  for(int i=0;i<snp + 2;i++){
      for(int j=0;j<size;j++){
	  grid[i][j]=0;
      }
  }

  for(int i=0;i<snp;i++){
      grid[i+1][size-1]=pow(sin((i + rank*snp)*M_PI/size),2);
  }

  for(int i=0;i<snp;i++){
      grid[i+1][0]=pow(cos((i + rank*snp)*M_PI/size),2);
  }

  double* ghost_fs=new double[size];
  double* ghost_bs=new double[size];
  double* ghost_fr=new double[size];
  double* ghost_br=new double[size];

  double** grid_n = new double*[snp+2];
  for (int i=0;i<snp+2;i++){
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
	  grid_n[i][j]=grid[i][j];
      }
  }

  for(int count=0;count<step_count;count++){
      for(int i=1;i<snp+1;i++){
	  for(int j=1;j<size-1;j++){
	      grid_n[i][j]=grid[i][j]+kappa*dt*(grid[i-1][j]+grid[i+1][j]+grid[i][j-1]+grid[i][j+1]-4*grid[i][j])/dx2;
	  }
      }
      for(int i=0;i<snp+2;i++){
	  for(int j=0;j<size;j++){
	      grid[i][j]=grid_n[i][j];
	  }
      }
      for(int i=0;i<size;i++){
	  ghost_bs[i]=grid[1][i];
	  ghost_fs[i]=grid[snp][i];
      }

      int rf=(rank+1)%n_pro;
      int rb=(rank-1)%n_pro;
      if(rb<0) rb=n_pro-1;

      int tag1=1;
      int tag2=2;

      MPI_Isend(&ghost_fs, size, MPI_DOUBLE, rf , tag1, MPI_COMM_WORLD, &request[0]);
      cout<<"Proc " <<rank<<" sending to  "<<rf<<" with tag "<<tag1<<endl;
      MPI_Isend(&ghost_bs, size, MPI_DOUBLE, rb , tag2, MPI_COMM_WORLD, &request[1]);
      MPI_Irecv(&ghost_fr, size, MPI_DOUBLE, rf , tag2, MPI_COMM_WORLD, &request[2]);
      MPI_Irecv(&ghost_br, size, MPI_DOUBLE, rb , tag1, MPI_COMM_WORLD, &request[3]);

      MPI_Waitall(4, request, Stat);
      printf("Received");

      for(int i=0;i<size;i++){
	  cout<<rank<<"in for loop "<<ghost_fr[i]<<" "<<ghost_br[i]<<endl;
	  grid[0][i]= ghost_fr[i];
	  grid[size/n_pro+1][i]=ghost_br[i];
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
      MPI_Finalize();
  }
}



