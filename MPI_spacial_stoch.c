/**
 * Nathan Breedlove, Eliot Dixon, Kenny Tram
 * 
 * Compile: mpiCC -g -Wall -o spacial.out MPI_spacial_stoch.c
 * Run: mpiexec -n [# of processors] spacial.out
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <deque>

double TIME_MAX = 3000; //arbitrary time length the simulation runs for
int TEST_RANK = 1; //controls which processor the output prints information about
int WRITE_FILES = 0; //1 if producing visualization of output; 0 otherwise

struct States {
    int i;
    int j;
    int direction; //0, 1, 2, 3
    int receiver;
    int sender;
    double t; //time
};

void writeData(double t, int i, int j, int defuseDirection);
void writeInitial(int **x, int n, int my_rank, int comm_sz);
void writeFromHistory(std::deque<States> history, int n, int my_rank, int comm_sz);
int** seedCells(int max_Pop, int numb_Cells, int row_Size, int col_Size);
int getIdealPop(int **x, int n);
int getDefuseDirection(int i, int j, int n, int my_rank, int comm_sz);
int getDirectionNum(int i, int j, int n,int my_rank, int comm_sz);
void reactionsInCell(int i, int j, int** reactions);
void sendData(int row, int my_rank, int direction, int local_n, double t);
std::deque<States> rollback(int **x, std::deque<States> history, double goal_time, int my_rank, int comm_sz);
std::deque<States> updateFuture(std::deque<States> future, int culprit, double goal_time, int my_rank);

int main(){
	int comm_sz;
    int my_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	int n = 8;
    int local_n = n/comm_sz;
    int i, j, row, col, proceed;
	int rollFlag = 0;
	int flag = 0;
    double t = 0;
    double tau;
    double buf[3];
    double time_buf;
    double reactionNum; //random double between 0 and the sum of the propensity functions; decides which reaction takes place
    double rateCon = 1;	//rate constant
    double propSum = 0; //sum of cells' propensity functions
    double cumilSum = 0; //cumilitive sum of propensity functions
   
    MPI_Request request;
    MPI_Status status;

    srand(time(NULL) + my_rank);
    
    std::deque<States> history;
    std::deque<States> future;
    
    int **x = (int **)malloc(n* sizeof(int*)); //array holding all molecules
    double **propFunc = (double **)malloc(n* sizeof(double*)); // array holding propensity functions corresponding to each cell of array x
    int **reactions = (int **)malloc(n* sizeof(int*)); //array containing each cells number of reactions
   
	for (i = 0; i <  n; i++) {
		x[i] = (int *)malloc((local_n)* sizeof(int));
		reactions[i] = (int *)malloc((local_n) * sizeof(int));
		propFunc[i] = (double *)malloc((local_n) * sizeof(double));
	}

	for(i = 0; i < n; i++) { //Initializing x
		for(j = 0; j < local_n; j++) {
			x[i][j] = rand()%10; //my_rank*local_n+j;
		}
	}
	
	if(WRITE_FILES) {
	    writeInitial(x, n, my_rank, comm_sz);
	}
	
	while(t < TIME_MAX) {	
		proceed = 1;
		cumilSum = 0;
		propSum = 0;
		
		if(my_rank == TEST_RANK) {
			printf("Rank: %d at time %lf\n", my_rank, t);
			for(i = 0; i < n; i++) {
				for(j = 0; j < local_n; j++) {
					printf("%d\t", x[i][j]);
				}
				printf("\n");
			}
			printf("\n");
		}
		
		MPI_Iprobe(MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &rollFlag, MPI_STATUS_IGNORE);
		
		if(rollFlag != 0){
			MPI_Irecv(&time_buf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &request);
			
			if(my_rank == TEST_RANK) {
			    printf("Rank: %d ... Rank: %d rolled back, affecting my history and future up to time %lf\n", my_rank, status.MPI_SOURCE, time_buf);
			}
			
			history = rollback(x, history, time_buf, my_rank, comm_sz);
			t = history.back().t;
			
			future = updateFuture(future, status.MPI_SOURCE, time_buf, my_rank);
			
			if(my_rank == TEST_RANK) {
			    printf("Rank: %d at time %lf\n", my_rank, t);
		        for(i = 0; i < n; i++) {
			        for(j = 0; j < local_n; j++) {
				        printf("%d\t", x[i][j]);
			        }
			        printf("\n");
		        }
		        printf("\n");
			}
			
			rollFlag = 0;
		}
		
		MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &flag, &status);
	
		if(flag != 0) {
			MPI_Irecv(&buf, 3, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &request);  //changebuf
			
			//printf("Rank: %d Current Time: %.0lf ComFlag Recv at (%.0lf,%.0lf) from time %lf\n",my_rank,t, buf[0],buf[1], buf[2]);
			
			flag = 0;
			
			struct States recv_state;
			
			recv_state.i = buf[0];
			recv_state.j= buf[1];
			recv_state.direction = -2;
			recv_state.receiver = -1;
			recv_state.sender = status.MPI_SOURCE;
			recv_state.t = buf[2];
			
			if(buf[2] < t){
			    if(my_rank == TEST_RANK) {
				    printf("Rank: %d received something at [%d,%d] at PAST time %lf from Rank %d\n", my_rank, (int)buf[0], (int)buf[1], buf[2], status.MPI_SOURCE);
				}
				history = rollback(x, history, buf[2], my_rank, comm_sz);
				history.push_back(recv_state);
				
				t = history.back().t;

				x[(int)buf[0]][(int)buf[1]]++;
				
				if(my_rank == TEST_RANK) {
				    printf("Rank: %d now applying received event...\n", my_rank);
				    
				    printf("Rank: %d at time %lf\n", my_rank, t);
			        for(i = 0; i < n; i++) {
				        for(j = 0; j < local_n; j++) {
					        printf("%d\t", x[i][j]);
				        }
				        printf("\n");
			        }
			        printf("\n");
				}
			}
			else{
			    if(my_rank == TEST_RANK) {
				    printf("Rank: %d received something at [%d,%d] at FUTURE time %lf from Rank %d. Storing in FUTURE deque...\n", my_rank, (int)buf[0], (int)buf[1], buf[2], status.MPI_SOURCE);
				}
				
				if(future.empty()) {
				    future.push_back(recv_state);
				}
				else{
				    std::deque<States>::iterator it = future.begin();
				    
				    while(it < future.end()) {
				        struct States temp_state = *it;
                        if(temp_state.t < recv_state.t) {
                            ++it;
                        }
                        else {
                            future.insert(it, recv_state);
                            break;
                        }
                    }
                    
                    if(it == future.end()) {
                        future.insert(it, recv_state);
                    }
				}
        	    
        	    if(my_rank == 1) {
				    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	                printf("****************** Rank: %d, PRINTING FUTURE DEQUE\n\n", my_rank);
	                for(unsigned i = 0; i < future.size(); i++) {
                        struct States tempState = future.at(i);
                        printf("rank: %d, i: %d, j: %d, direction: %d, sender: %d, receiver: %d, time: %lf\n", my_rank, tempState.i, tempState.j, tempState.direction, tempState.sender, tempState.receiver, tempState.t);
                    }
                    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
                }
			}
		}
			
		for(i = 0; i < n; i++){
			for(j = 0; j < local_n; j++){ //creating propensity Functions
				propFunc[i][j] = (double)x[i][j] * getDirectionNum(i, j, n, my_rank, comm_sz) * rateCon;
				propSum += propFunc[i][j];
			}	
		}
		
		double timeRand = (double) rand()/RAND_MAX; //generate random value
		tau = log(1/timeRand)/(propSum*.001); //generate time step (time next reaction will take to complete)
		t += tau; //increment our proc's time 

		if(!future.empty()){
			struct States future_state = future.front();
			
			if(future_state.t < t){	//if future reaction occured within this step fire it
				if(my_rank == TEST_RANK) {
				    printf("!!! Rank: %d applying future event at [%d,%d] at time %lf from Rank %d !!!\n", my_rank, future_state.i, future_state.j, future_state.t, future_state.sender);
				}
				x[future_state.i][future_state.j] += 1;
				t = future_state.t;
				
				history.push_back(future_state);
				future.pop_front();
				
				proceed = 0;
			}
		}
		
		if(proceed) {
			struct States current_state;
			reactionNum = (double) rand()/RAND_MAX * propSum;
			
			for(i = 0; i < n; i++){
				for(j = 0; j < local_n; j++){
					cumilSum += propFunc[i][j];
					if(cumilSum >= reactionNum){	
						row = i;
						col = j;
						x[row][col] --;
						break;
					}
				}
				if(cumilSum >= reactionNum){
					break;
				}
			}
			
			current_state.i = row;
			current_state.j = col;
			current_state.receiver = -1;
			current_state.sender = -1;
			current_state.t = t;
			
			int defuse_direction = getDefuseDirection(row, col, n, my_rank, comm_sz);
			
			if((defuse_direction == 3 && col == local_n-1) ||
               (defuse_direction == 1 && col == 0)) {
                    if(defuse_direction == 3){
	                    current_state.receiver = my_rank+1;
                    }
                    else if(defuse_direction == 1){
	                    current_state.receiver = my_rank-1;
                    }
                    
	                sendData(row, my_rank, defuse_direction, local_n, t);
	                
                    defuse_direction = -1;
            }
			
			switch(defuse_direction){
				case(-1):
					break;
				case(0):
					x[row+1][col] += 1;
					break;
				case(1):
					x[row][col-1] += 1;
					break;
				case(2):
					x[row-1][col] += 1;
					break;
				case(3):
					x[row][col+1] += 1;
					break;
			}
			
			current_state.direction = defuse_direction;
			
			history.push_back(current_state);
			
			if(my_rank == TEST_RANK) {
                printf("Rank: %d ... [%d,%d] sent direction %d; completed at time %lf\n", my_rank, i, j, defuse_direction, t);
		    }
		}
		/*
		
		*/
	}
    
    if(WRITE_FILES) {
        writeFromHistory(history, n, my_rank, comm_sz);
    }
    
	MPI_Finalize();
	return 0;
}//MAIN



// =================================================================================
void writeData(double t, int i, int j, int defuseDirection){
	FILE *fp;
	fp = fopen("data.txt", "a");
	fprintf(fp,"%lf\t(%d,%d)\t%d\n",t,i,j,defuseDirection);
	fclose(fp);
}

// =================================================================================
void writeInitial(int **x, int n, int my_rank, int comm_sz) {
    FILE *file;
    char filename[20];
    sprintf(filename, "rank%02d_%06d.txt", my_rank, 0);
    file = fopen(filename, "w");
    
    int i, j;
    int local_n = n/comm_sz;
    
    if(file == NULL) {
        printf("ERROR: Output file not found. Terminating...\n");
	    exit(1);
    }
    else {
        fprintf(file, "%lf\n", 0.0); //prints time
        
        for(i = 0; i < n; i++) {
	        for(j = 0; j < local_n; j++) {
		        fprintf(file, "%d\t", x[i][j]);
	        }
	        fprintf(file, "\n");
        }
        fclose(file);
    }
}

// =================================================================================
void writeFromHistory(std::deque<States> history, int n, int my_rank, int comm_sz) {
    FILE *curr_file;
    FILE *prior_file;
    char filename[20];
    //char buffer[100];
    
    int i, j, row, col, direction;
    int local_n = n/comm_sz;
    int step = 1;
    
    int grid[n][local_n];
    
    for(unsigned k = 0; k < history.size(); k++) {
        sprintf(filename, "rank%02d_%06d.txt", my_rank, step-1);
        prior_file = fopen(filename, "r");
    
        sprintf(filename, "rank%02d_%06d.txt", my_rank, step);
        curr_file = fopen(filename, "w");
            
        struct States tempState = history.at(k);
        
	    if(prior_file == NULL) {
		    printf("ERROR: Prior output file not found. Terminating...\n");
		    exit(1);
	    }
	    else if(curr_file == NULL) {
		    printf("ERROR: Current output file not found. Terminating...\n");
		    exit(1);
	    }
	    else {
	        // *** Read and copy array from prior_file
	        fprintf(curr_file, "%lf\n", tempState.t); //prints time
	        
	        fscanf(prior_file, "%*[^\n]\n");
	        for(i = 0; i < n; i++) {
	            for(j = 0; j < local_n; j++) {
	                fscanf(prior_file, "%d", &grid[i][j]);
	            }
	        }
	        
		    // *** Print curr_file based on prior array and applying recorded move
		    row = tempState.i;
		    col = tempState.j;
		    direction = tempState.direction;
		    
		    switch(direction) {
		        case(-2): //RECEIVE
		            grid[row][col] += 1;
		            break;
			    case(-1): //SEND
			        grid[row][col] -= 1;
				    break;
			    case(0):
			        grid[row][col] -= 1;
				    grid[row+1][col] += 1;
				    break;
			    case(1):
			        grid[row][col] -= 1;
				    grid[row][col-1] += 1;
				    break;
			    case(2):
			        grid[row][col] -= 1;
				    grid[row-1][col] += 1;
				    break;
			    case(3):
			        grid[row][col] -= 1;
				    grid[row][col+1] += 1;
				    break;
		    }
		    
		    for(i = 0; i < n; i++) {
	            for(j = 0; j < local_n; j++) {
	                fprintf(curr_file, "%d\t" , grid[i][j]);
	            }
	            fprintf(curr_file, "\n");
	        }
		    
	        fclose(curr_file);
	        fclose(prior_file);
        }
	    step++;
    }
}

// =================================================================================
int getDefuseDirection(int i, int j, int n, int my_rank, int comm_sz){
    int validDirections[4];
    int index = 0;
    int direction;
    int local_n = n /comm_sz; 
    
    if(my_rank != 0 || j - 1 >= 0){
        validDirections[index] = 1;
        index++;
    }
    if(i - 1 >= 0){
        validDirections[index] = 2;
        index++;
    }
    if(my_rank != comm_sz-1 || j + 1 <= local_n-1){
        validDirections[index] = 3;
        index++;
    }
    if(i + 1 <= n-1){
        validDirections[index] = 0;
    }
    else{
        index--;
    }
    
    direction = validDirections[rand()%(index +1)];
   
    return direction;
} 

// =================================================================================
int getDirectionNum(int i, int j, int n,int my_rank, int comm_sz){
	int directions =4;
	int local_n = n/comm_sz; 
	if(j - 1 < 0 && my_rank == 0){
		directions--;
	}
	if(j + 1 > local_n-1 && my_rank == comm_sz -1){
		directions --;
	}
	if(i - 1 < 0){
		directions --;
	}
	if(i + 1 > n-1){
		directions --;
	}

	return directions;
}

// =================================================================================
int getIdealPop(int **x, int n) {
    int i, j;
    
    int sum = 0;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            sum += x[i][j];
        }
    }
    
    //Just for verifying math...
    printf("Total population is: %d\nTotal cells is: %d\n", sum, (i*j));
    
    return (sum / (i * j));
}


// =================================================================================
/*
function: creates a two dimensional array, simulating the grid our rections take place in
Parameters -
max_Pop: max population to generate within any given cell
numb_Cells: number of cells (randomly chosen within our grid) to generate population within (others have population of 0)
row_Size: number of rows
col_Size: number of columns
returns : grid_Array (pointer to the populated/seeded 2d grid)
*/
int** seedCells(int max_Pop, int numb_Cells, int row_Size, int col_Size){
	int i,j,k,numb;
	int** grid_Array = (int **)calloc((row_Size),sizeof(int*)); 			
	for (i = 0; i <  row_Size; i++){
		grid_Array[i] = (int *)calloc((col_Size),sizeof(int));
	}//create 2d array, intiialize all to zero

	k=0;
	while(k < numb_Cells){
		for (i = 0; i <  row_Size; i++){
			for(j = 0; j < col_Size; j++){
				double prob = (double)rand()/RAND_MAX; //gen prob 0.0-1.0
				printf("Something\n");
				if(prob <= (double)numb_Cells/(row_Size*col_Size) && grid_Array[i][j] == 0){
					numb = rand()%max_Pop;
					printf("NUM %d\n", numb);
					grid_Array[i][j] = numb;
					k++;
				}//if prob <= number_of_cells_we_want/grid_size AND this given cell is empty (0)

				if(k >= numb_Cells){
					break;
				}//if we've seeded numb_Cells, break

			}//nested for,navigate inner arrays

			if(k >= numb_Cells){
				break;
			}//if we've seeded numb_Cells, break

		}//navigate outer array
	}//while there are cells left to populate

	return grid_Array;
}//seedCells

// =================================================================================
void sendData(int row, int my_rank, int direction, int local_n, double t){
	int col = -1;
	int neighbor = -1;  
	MPI_Request request;
	
	if(direction == 3){
        neighbor = my_rank +1;
        col = 0;
    }
    else if(direction == 1){
        neighbor = my_rank -1;
        col = local_n -1;
    }
	
	double rowCol[3] = {(double)row, (double)col, t};
	
	/*
	if(my_rank == TEST_RANK) {MPI_STATUS_IGNORE
		printf("Rank: %d, sent something in direction %d to [%d,%d] of Rank %d!\n",my_rank, direction,row,col, neighbor);
	}
	*/
	
	MPI_Isend(rowCol, 3, MPI_DOUBLE, neighbor, 1, MPI_COMM_WORLD, &request);
}

// =================================================================================
std::deque<States> rollback(int **x, std::deque<States> history, double goal_time, int my_rank, int comm_sz) {  
	struct States temp_state;
	double current_time;
	int row, col, direction, receiver, sender;
	
	int left_neighbor = my_rank-1;
	double left_recv_time_buf = -1;
	double left_sent_time_buf = -1;
	
	int right_neighbor = my_rank+1;
	double right_recv_time_buf = -1;
	double right_sent_time_buf = -1;
	
	MPI_Request request;
	
	if(my_rank == TEST_RANK) {
	    printf("===================================================================================================\n");
		printf("$$$$$$$$$$$$$$$$$$ Rank: %d, ROLLING BACK, GOAL TIME: %lf\n\n", my_rank, goal_time);
	}
	
    while(!history.empty()) {
        temp_state = history.back();
        current_time = temp_state.t;
        
        if(current_time > goal_time) {
            direction = temp_state.direction;
            row = temp_state.i;
            col = temp_state.j;
		    sender = temp_state.sender;
		    receiver = temp_state.receiver;
		    
		    if(my_rank == TEST_RANK) {
			    printf("> Rank: %d, ROLLING BACK EVENT THAT OCCURED AT TIME: %lf\n", my_rank, current_time);
		        printf("\t> i: %d, j: %d, direction: %d, receiver: %d, sender: %d, time: %lf\n", temp_state.i, temp_state.j, temp_state.direction, temp_state.receiver, temp_state.sender, temp_state.t);
		    }
            
            switch(direction) {
                /*
                case(-3) //KILL SIGNAL RECEIVED
                    if(my_rank != 0) {
                        
                    
                    }
                
                */
			    case(-2): //RECEIVED
			        if(my_rank == TEST_RANK) {
			            printf("\t\t> It was a RECEIVE! Decrementing [%d,%d]...\n\n", row, col);
		            }
		            
				    x[row][col] -= 1;
				    
				    if(sender == right_neighbor) {
				        right_recv_time_buf = temp_state.t;
				    }
				    else {
				        left_recv_time_buf = temp_state.t;
				    }
				    
				    break;
		        case(-1): //SENT
		            if(my_rank == TEST_RANK) {
			            printf("\t\t> It was a SEND! Incrementing [%d,%d]...\n\n", row, col);
		            }
		            
		            x[row][col] += 1;
		            
		            if(receiver == right_neighbor) {
				        right_sent_time_buf = temp_state.t;
				    }
				    else {
				        left_sent_time_buf = temp_state.t;
				    }
				   
			        break;
		        case(0):
		            if(my_rank == TEST_RANK) {
			            printf("\t\t> It was a NORMAL REACTION! incrementing [%d,%d], decrementing [%d,%d]...\n\n", row, col, row+1, col);
		            }
		            
		            x[row][col] += 1;
			        x[row+1][col] -= 1;
			        break;
		        case(1):
		            if(my_rank == TEST_RANK) {
			            printf("\t\t> It was a NORMAL REACTION! incrementing [%d,%d], decrementing [%d,%d]...\n\n", row, col, row, col-1);
		            }
		            
		            x[row][col] += 1;
			        x[row][col-1] -= 1;
			        break;
		        case(2):
		            if(my_rank == TEST_RANK) {
			            printf("\t\t> It was a NORMAL REACTION! incrementing [%d,%d], decrementing [%d,%d]...\n\n", row, col, row-1, col);
		            }
		            
		            x[row][col] += 1;
			        x[row-1][col] -= 1;
			        break;
		        case(3):
		            if(my_rank == TEST_RANK) {
			            printf("\t\t> It was a NORMAL REACTION! incrementing [%d,%d], decrementing [%d,%d]...\n\n", row, col, row, col+1);
		            }
		            
		            x[row][col] += 1;
			        x[row][col+1] -=1;
			        break;
	        }//reverse events
	        
            history.pop_back();
        }//if current event's time > time to rollback to
        else {
            break;
        }
    }
    
	if(my_rank == TEST_RANK) {
	    if(history.empty()) {
	        current_time = 0;
	    }
	    
		printf("$$$$$$$$$$$$$$$$$$ Rank: %d, ROLLING BACK COMPLETE, TIME IS NOW %lf\n\n", my_rank, current_time);
		printf("===================================================================================================\n\n");
	}
	
    if(my_rank != 0) {
        if(left_recv_time_buf > -1 || left_sent_time_buf > -1) {
		    if(left_recv_time_buf == -1) { // use sent time
			    MPI_Isend(&left_sent_time_buf, 1, MPI_DOUBLE, left_neighbor, 3, MPI_COMM_WORLD, &request);
		    }
		    else if(left_sent_time_buf == -1) { //use recieve time
			    MPI_Isend(&left_recv_time_buf, 1, MPI_DOUBLE, left_neighbor, 3, MPI_COMM_WORLD, &request);
		    }
		    else {
			    if(left_recv_time_buf > left_sent_time_buf) { //use sent time
				    MPI_Isend(&left_sent_time_buf, 1, MPI_DOUBLE, left_neighbor, 3, MPI_COMM_WORLD, &request);
			    }
			    else if(left_sent_time_buf > left_recv_time_buf) { //use recv time
				    MPI_Isend(&left_recv_time_buf, 1, MPI_DOUBLE, left_neighbor, 3, MPI_COMM_WORLD, &request);
			    }
	        }
	    }
	}
	
	if(my_rank != comm_sz-1) {
	    if(right_recv_time_buf > -1 || right_sent_time_buf > -1) {
		    if(right_recv_time_buf == -1) { // use sent time
			    MPI_Isend(&right_sent_time_buf, 1, MPI_DOUBLE, right_neighbor, 3, MPI_COMM_WORLD, &request);
		    }
		    else if(right_sent_time_buf == -1) { //use recieve time
			    MPI_Isend(&right_recv_time_buf, 1, MPI_DOUBLE, right_neighbor, 3, MPI_COMM_WORLD, &request);
		    }
		    else {
			    if(right_recv_time_buf > right_sent_time_buf) { //use sent time
				    MPI_Isend(&right_sent_time_buf, 1, MPI_DOUBLE, right_neighbor, 3, MPI_COMM_WORLD, &request);
			    }
			    else if(right_sent_time_buf > right_recv_time_buf) { //use recv time
				    MPI_Isend(&right_recv_time_buf, 1, MPI_DOUBLE, right_neighbor, 3, MPI_COMM_WORLD, &request);
			    }
	        }
	    }
	}
    
    return history;
}//rollback

// =================================================================================
std::deque<States> updateFuture(std::deque<States> future, int culprit, double goal_time, int my_rank) {
    std::deque<States>::iterator it = future.end()-1;
    
    if(my_rank == TEST_RANK) {
		
		printf("===================================================================================================\n");
		printf("@@@@@@@@@@@@@@@@@@ Rank: %d, UPDATING FUTURE, GOAL TIME: %lf\n\n", my_rank, goal_time);
	}
	
    while(it >= future.begin()) {
        struct States temp_state = *it;
        
        if(temp_state.t >= goal_time) {
            if(temp_state.sender == culprit) {
                if(my_rank == TEST_RANK) {
			        printf("> Rank: %d, REMOVING FUTURE EVENT that was from Rank %d at time %lf\n", my_rank, temp_state.sender, temp_state.t);
                }
                future.erase(it);
            }
        }
        else {
            break;
        }
        
        --it;
    }
    
    if(my_rank == TEST_RANK) {
		printf("===================================================================================================\n\n");
	}
    
	return future;
}
