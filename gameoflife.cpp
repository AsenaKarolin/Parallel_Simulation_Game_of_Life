//Name: Asena Karolin Özdemir
//ID: 2018400366
//Compiling Status: Compiling
//Working Status: Working
//Periodical approach is used for the boundaries
//Checkered splits are implemented
#include "mpi.h"
#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#define N 360
using namespace std;
//This function sends the arrays to the neighbours on the 4 sides (north, south, west, east)
void sendTo4Sides(int norths[],int souths[],int wests[], int easts[], int len, int tag, int northnei, int southnei,int westnei, int eastnei){
	MPI_Send(norths, len, MPI_INT, northnei, tag, MPI_COMM_WORLD);
	MPI_Send(souths, len, MPI_INT, southnei, tag, MPI_COMM_WORLD);
	MPI_Send(wests, len, MPI_INT, westnei, tag, MPI_COMM_WORLD);
	MPI_Send(easts, len, MPI_INT, eastnei, tag, MPI_COMM_WORLD);
}
//This function receives the arrays from the neighbours on the 4 sides (south, north, east, west)
void receiveFrom4Sides(int southr[], int northr[], int westr[], int eastr[], int len, int tag, int southnei, int northnei, int eastnei, int westnei){
	MPI_Recv(southr, len, MPI_INT, southnei, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Recv(northr, len, MPI_INT, northnei, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Recv(eastr, len, MPI_INT, eastnei, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Recv(westr, len, MPI_INT, westnei, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}
//This function sends the corner points to the diagonal neighbours (northwest, northeast, southwest, southeast)
void sendToCorners(int *nws, int *nes, int *sws, int *ses, int tag, int nwnei, int nenei, int swnei, int senei){
	MPI_Send(nws, 1, MPI_INT, nwnei, tag, MPI_COMM_WORLD);
	MPI_Send(nes, 1, MPI_INT, nenei, tag, MPI_COMM_WORLD);
	MPI_Send(sws, 1, MPI_INT, swnei, tag, MPI_COMM_WORLD);
	MPI_Send(ses, 1, MPI_INT, senei, tag, MPI_COMM_WORLD);
}
//This function receives the corner points from the diagonal neighbours (southeast, southwest, northeast, northwest)
void receiveFromCorners(int *ser, int *swr, int *ner, int *nwr, int tag, int senei, int swnei, int nenei, int nwnei){
	MPI_Recv(ser, 1, MPI_INT, senei, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Recv(swr, 1, MPI_INT, swnei, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Recv(ner, 1, MPI_INT, nenei, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Recv(nwr, 1, MPI_INT, nwnei, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}
int main(int argc, char *argv[])
{
    //Initialization
    int i,j,k; //variables used to iterate in loops
    int bound=0,bound2=0; //boundaries for iterations in loops - used to split the input 
    int rank, size, next, prev, tag = 201;
    int iterations = atoi(argv[3]); //number of total iterations
    int it; //count of iterations
    int W,E,NO,S,NW,NE,SW,SE; //neighbours of each creaure (0 or 1) - NO stands for north, W stands for west, NW stands for northwest and so on
    int neigh; //count of neighbours of a creature (count of 1s)
    int grid[N][N]; //the total input
    // creates mpi
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //creates message to send the parts of the grid
    int squares = sqrt(size-1); //the number of splitted parts in each row and each column
    int len=N/squares; //the length of each side of a splitted part
    int message[len][len]; //message that is sent to each worker process from the manager process
    int temp[len][len]; //temporary array that is used to keep the information of the message during iterations of the game
    int norths[len],souths[len],easts[len],wests[len]; //arrays that are sent to neighbours - e.g. wests is the array that is sent to the west neighbour
    int northr[len],southr[len],eastr[len],westr[len]; //arrays that are received from neighbours - e.g. westr is the array that is received from the west neighbour
    int nws,nes,sws,ses,nwr,ner,swr,ser; //variables that are sent to and received from the neighbours in the diagonal - e.g. nws is the variable that is sent to the northwest neighbour and ner is the variable that is received from the northeast neighbour
    int northnei,southnei,eastnei,westnei,nwnei,swnei,nenei,senei; //the ranks of the neighbours for each process - e.g. northnei represents the rank of the north neighbour and senei represents the rank of the southeast neighbour
	if (rank==0){
		//reads from inputfile
		std::ifstream inputfile;
		inputfile.open(argv[1]);
		for (i=0; i<N; i++){
			for (j=0; j<N ; j++){
				inputfile >> grid[i][j];
			}
		}
		inputfile.close();
		//partitions the input, computes and sends a message to each process
		for (k=1; k<size; k++){
			//k represents the rank of the processes
			if(k%squares==0){
				bound2=(k/squares)-1;
			}
			else{
				bound2=k/squares;
			}
			for (i=len*bound2; i<len*(bound2+1); i++){
			//bound2 is the bound for the range of the rows that are send to the process with rank k
				if((k%squares)!=0){
					bound=(k%squares)-1;
				}
				else{
					bound=squares-1;
				}
				for (j=len*bound; j<len*(bound+1); j++){
				//bound is the bound for the range of the columns that are send to the process with rank k
					message[i%len][j%len] = grid[i][j];
				}
			}
			MPI_Send(message, len*len , MPI_INT, k, tag, MPI_COMM_WORLD);
		}
	}
	if (rank!=0) {
	//each process receives the message
		MPI_Recv(message, len*len, MPI_INT, 0, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		//the ranks of the neighbour processes are determined
		//sets northneighbour
		if (rank==squares){
			northnei=size-1;
		}
		else if (rank/squares!=0){
			northnei=rank-squares;
		}
		else{
			northnei=rank+size-1-squares;
		}
		//sets southneighbour
		if (rank+squares>size-1){
			southnei=rank-(size-1)+squares;
		}
		else{
			southnei=rank+squares;
		}
		//sets westneighbour
		if(rank%squares==1){
			westnei=rank+squares-1;
		}
		else{
			westnei=rank-1;
		}
		//sets eastneighbour
		if(rank%squares==0){
			eastnei=rank-squares+1;
		}
		else{
			eastnei=rank+1;
		}
		//sets northwestneighbour
		if(rank==1){
			nwnei=size-1;
		}
		else if(rank/squares==0){
			nwnei=rank+(size-1)-(squares+1);
		}
		else if(rank==squares){
			nwnei=rank+(size-1)-(squares+1);
		}
		else if(rank%squares==1){
			nwnei=rank-1;
		}
		else{
			nwnei=rank-squares-1;
		}
		//sets southwestneighbour
		if(rank==size-squares){
			swnei=squares;
		}
		else if(rank+squares>size-1){
			swnei=rank-size+squares;
		}
		else if(rank%squares==1){
			swnei=rank+(2*squares)-1;
		}
		else{
			swnei=rank+squares-1;
		}
		//sets northeastneighbour
		if(rank==squares){
			nenei=size-squares;
		}
		else if(rank/squares==0){
			nenei=rank+size-squares;
		}
		else if(rank%squares==0){
			nenei=rank-(2*squares)+1;
		}
		else{
			nenei=rank-squares+1;
		}
		//sets southeastneighbour
		if(rank==size-1){
			senei=1;
		}
		else if(rank+squares>size-1){
			senei=rank-(size-1)+squares+1;
		}
		else if(rank%squares==0){
			senei=rank+1;
		}
		else{
			senei=rank+squares+1;
		}
	}
	it=0;// iteration counter is set to 0
	while(it<iterations){ //the game is played as many times as specified by the variable 'iterations'
		if (rank!=0) {
			//each process creates the lines to be sent
			for (i=0; i<len; i++){
				norths[i]=message[0][i]; //sets the array that will be sent to the north neighbour
				souths[i]=message[len-1][i]; //sets the array that will be sent to the south neighbour
				wests[i]=message[i][0]; //sets the array that will be sent to the west neighbour
				easts[i]=message[i][len-1]; //sets the array that will be sent to the east neighbour
			}
			nws=message[0][0]; //sets the variable that will be sent to the northwest neighbour
			nes=message[0][len-1]; //sets the variable that will be sent to the northeast neighbour
			sws=message[len-1][0]; //sets the variable that will be sent to the southwest neighbour
			ses=message[len-1][len-1]; //sets the variable that will be sent to the southeast neighbour
		}
		//sends and receives the sides
		if(rank!=0 && (((rank%(2*squares)<=squares) && (rank%2==0) && (rank%(2*squares)!=0)) || ((rank%(2*squares)>squares) && (rank%2==1)))){
			//first the processes with even rank in odd numbered rows and the processes with odd rank in even numbered rows send their specified arrays to their north, south, west and east neighbours respectively
			//e.g. if the number of worker processes were 16, the processes with rank 2, 4, 5, 7, 10, 12, 13, 15 would send first
			sendTo4Sides(norths,souths,wests,easts,len,tag,northnei,southnei,westnei,eastnei);
			//the other processes receive these arrays and send their own arrays
			//then these processes receive the arrays of the other processes
			receiveFrom4Sides(southr,northr,westr,eastr,len,tag,southnei,northnei,eastnei,westnei);
			
		}
		if(rank!=0 && !(((rank%(2*squares)<=squares) && (rank%2==0) && (rank%(2*squares)!=0)) || ((rank%(2*squares)>squares) && (rank%2==1)))){
			//the other worker processes receive the arrays from their south, north, east, west neighbours respectively
			receiveFrom4Sides(southr,northr,westr,eastr,len,tag,southnei,northnei,eastnei,westnei);
			
			//then these processes send their specified arrays to their north, south, west and east neighbours respectively
			sendTo4Sides(norths,souths,wests,easts,len,tag,northnei,southnei,westnei,eastnei);
		}
		//sends and receives the corner points
		if(rank!=0 && (rank%(2*squares)<=squares) && (rank%(2*squares)!=0)){
			//first the processes in odd numbered rows send their corner points to their northwest, northeast, southwest, southeast neighbour respectively
			sendToCorners(&nws,&nes,&sws,&ses,tag,nwnei,nenei,swnei,senei);
			//the other processes receive these corner points and send their own corner points
			//then these processes receive the corner points of the other processes
			receiveFromCorners(&ser,&swr,&ner,&nwr,tag,senei,swnei,nenei,nwnei);
		}
		if(rank!=0 && !((rank%(2*squares)<=squares) && (rank%(2*squares)!=0))){
			//then the other processes receive the corner points from their southeast, southwest, northeast, northwest neighbours respectively
			receiveFromCorners(&ser,&swr,&ner,&nwr,tag,senei,swnei,nenei,nwnei);
			
			//then these processes send their corner points to their northwest, northeast, southwest, southeast neighbour respectively
			sendToCorners(&nws,&nes,&sws,&ses,tag,nwnei,nenei,swnei,senei);
		}
		if(rank!=0){
			//each process plays the game
			for (i=0; i<len; i++){
				for (j=0; j<len; j++){
					if(i!=0 && j!=0 && i!=len-1 && j!=len-1){
						//if the field is in the middle part of the game
						//uses the fields around it
						W=message[i][j-1];
						E=message[i][j+1];
						NO=message[i-1][j];
						S=message[i+1][j];
						NW=message[i-1][j-1];
						SW=message[i+1][j-1];
						NE=message[i-1][j+1];
						SE=message[i+1][j+1];
					}
					else if(i==0 && j!=0 && j!=len-1){
						//if the field is in the first row of the game but not on the corner
						//uses the fields around it and the array that was received from the north neighbour
						NO=northr[j];
						NW=northr[j-1];
						NE=northr[j+1];
						W=message[i][j-1];
						E=message[i][j+1];
						S=message[i+1][j];
						SW=message[i+1][j-1];
						SE=message[i+1][j+1];
					}
					else if(i==len-1 && j!=0 && j!=len-1){
						//if the field is in the last row of the game but not on the corner
						//uses the fields around it and the array that was received from the south neighbour
						W=message[i][j-1];
						E=message[i][j+1];
						NO=message[i-1][j];
						S=southr[j];
						NW=message[i-1][j-1];
						SW=southr[j-1];
						NE=message[i-1][j+1];
						SE=southr[j+1];
					}
					else if(j==0 && i!=0 && i!=len-1){
						//if the field is in the first column of the game but not on the corner
						//uses the fields around it and the array that was received from the west neighbour
						W=westr[i];
						E=message[i][j+1];
						NO=message[i-1][j];
						S=message[i+1][j];
						NW=westr[i-1];
						SW=westr[i+1];
						NE=message[i-1][j+1];
						SE=message[i+1][j+1];
					}
					else if(j==len-1 && i!=0 && i!=len-1){
						//if the field is in the last column of the game but not on the corner
						//uses the fields around it and the array that was received from the east neighbour
						W=message[i][j-1];
						E=eastr[i];
						NO=message[i-1][j];
						S=message[i+1][j];
						NW=message[i-1][j-1];
						SW=message[i+1][j-1];
						NE=eastr[i-1];
						SE=eastr[i+1];
					}
					else if(i==0 && j==0){
						//if the field is in the northwest corner of the game
						//uses the corner point received from the northwest neighbour and the arrays from the north and west neighbours, and the fields around it
						W=westr[i];
						E=message[i][j+1];
						NO=northr[j];
						S=message[i+1][j];
						NW=nwr;
						SW=westr[i+1];
						NE=northr[j+1];
						SE=message[i+1][j+1];
					}
					else if(i==0 && j==len-1){
						//if the field is in the northeast corner of the game
						//uses the corner point received from the northeast neighbour and the arrays from the north and east neighbours, and the fields around it
						W=message[i][j-1];
						E=eastr[i];
						NO=northr[j];
						S=message[i+1][j];
						NW=northr[j-1];
						SW=message[i+1][j-1];
						NE=ner;
						SE=eastr[i+1];
					}
					else if(j==0 && i==len-1){
						//if the field is in the southwest corner of the game
						//uses the corner point received from the southwest neighbour and the arrays from the south and west neighbours, and the fields around it
						W=westr[i];
						E=message[i][j+1];
						NO=message[i-1][j];
						S=southr[j];
						NW=westr[i-1];
						SW=swr;
						NE=message[i-1][j+1];
						SE=southr[j+1];
					}
					else if(i==len-1 && j==len-1){
						//if the field is in the southeast corner of the game
						//uses the corner point received from the southeast neighbour and the arrays from the south and east neighbours, and the fields around it
						W=message[i][j-1];
						E=eastr[i];
						NO=message[i-1][j];
						S=southr[j];
						NW=message[i-1][j-1];
						SW=southr[j-1];
						NE=eastr[i-1];
						SE=ser; 
					}
					neigh=W+E+NO+S+NE+NW+SE+SW; //calculates number of total neighbours that are alive
					//plays the game and stores the next stage of the game in the temp array
					if(message[i][j]==0 && neigh==3){
						temp[i][j]=1;
					}
					else if(message[i][j]==0 && neigh!=3){
						temp[i][j]=0;
					}
					else if(message[i][j]==1 && neigh<=3 && neigh>=2){
						temp[i][j]=1;
					}
					else {
						temp[i][j]=0;
					}
				}
			}
			//once the game is done the content of temp is copied into message again
			for (i=0; i<len; i++){
				for (j=0; j<len; j++){
					message[i][j]=temp[i][j];
				}
			}
		}
		it++; //increments iteration counter by 1
	}
	//game iterations are done
	if (rank!=0){
		//worker processes send their messages back to the manager process
		MPI_Send(message, len*len , MPI_INT, 0, tag, MPI_COMM_WORLD);
	}
	if (rank==0){
		//manager process receives the messages from the worker processes
		for (k=1; k<size; k++){
			//k represents the rank of the processes
			MPI_Recv(message, len*len, MPI_INT, k, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			if(k%squares==0){
				bound2=(k/squares)-1;
			}
			else{
				bound2=k/squares;
			}
			//bound2 is the bound for the range of the rows that are received from the process with rank k
			for (i=len*bound2; i<len*(bound2+1); i++){
				if((k%squares)!=0){
					bound=(k%squares)-1;
				}
				else{
					bound=squares-1;
				}
				//bound is the bound for the range of the columns that are received from the process with rank k
				for (j=len*bound; j<len*(bound+1); j++){
					//recreates the original grid using the messages
					grid[i][j]=message[i%len][j%len];
				}
			}
		}
		//writes to the output file
		std::ofstream outputfile;
		outputfile.open(argv[2]);
		for (i=0; i<N; i++){
			for (j=0; j<N ; j++){
				outputfile << grid[i][j] << " ";
			}
			outputfile << endl;
		}
		outputfile.close();
	}
	//finalizes MPI
	MPI_Finalize();
	return 0;
}
