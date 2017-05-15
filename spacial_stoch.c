#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include <deque>

void writeFile(int **x, int n, int reaction_num);
int getDefuseDirection(int i, int j, int n);
void writeData(double t, int i, int j, int defuse_direction);
int getDirectionNum(int i, int j, int n);
int getIdealPop(int **x, int n);
void printGrid(int **x, int n);
int** seedCells(int max_Pop, int numb_Cells, int grid_Size);
int** seedFromInput(int n);
/*
struct States {
    int i;
    int j;
    int direction; //0, 1, 2, 3
    double t; //time
};
*/
main() {
    int animate = 1;
	int n = 50;
	double t = 0;
	double tau;
	double end = 50000;
	int reaction_count = 0;
	int i, j;
	double prop_sum = 0;      // sum of cells' propensity functions
	double reaction_num;      //random double between 0 and the sum of the propensity functions to decide which 
	int defuse_direction = -1; //random int between 0 and 3 decides which direction a molecule will defuse to
	double cumil_sum = 0;    ////cumilitive sum of propensity functions
	
	int **x = (int **)malloc((n)* sizeof(int*)); 				//array holding all molecules
	double **prop_func = (double **)malloc(n* sizeof(double*));  // array holding propensity functions corresponding to each cell of array x
	
	for (i = 0; i <  n; i++){
		x[i] = (int *)malloc((n)* sizeof(int));
		prop_func[i] = (double *)malloc(n* sizeof(double));
		
	}
    
    //std::deque<States> history;
    
	double rate_con =  2; // rate at which reactions occur
	srand(time(NULL));
	
	//printf("INITIAL GRID\n");
	//x = seedCells(5000, 1, n);
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			x[i][j] = 0;
		}
	}
	
	x[9][10] = 200;
	x[10][9] = 200;
	x[10][10] = 200;
	x[10][11] = 200;
	x[11][10] = 200;
	
	//x = seedFromInput(n);
	
	//printGrid(x, n);
	
	while(t < end){
		cumil_sum = 0;
		prop_sum = 0;
		
		//struct States curr_state;
		//curr_state.t = t;
		
		for(i = 0; i < n; i++) {
			for(j = 0; j < n; j++) {
				prop_func[i][j] = (double)x[i][j] * (double)getDirectionNum(i, j, n) * rate_con;
				prop_sum += prop_func[i][j];
			}
		}
		
		reaction_num = ((double)rand() / (double)RAND_MAX) * prop_sum;
		//printf("Reaction Number %lf\n", reaction_num);
		
		for(i = 0; i < n; i++){
		    for(j = 0; j < n; j++){
				cumil_sum += prop_func[i][j];
				
				if(reaction_num <= cumil_sum){	
						break;
				}
			}
			if(reaction_num <= cumil_sum){
					break;
			}
		}
		
        //curr_state.i = i;
        //curr_state.j = j;
        
		//0=south, 1=west, 2=north=, 3=east 
		defuse_direction = getDefuseDirection(i,j,n);
		x[i][j] -= 1;
		
		switch(defuse_direction){
		case(0):
			x[i+1][j] +=1;
			break;
		case(1):
		    x[i][j-1] += 1;
			break;
		case(2):
			x[i-1][j] += 1;
			break;
		case(3):
		    x[i][j+1] += 1;
			break;
		}
        
        //curr_state.direction = defuse_direction;
        //history.push_back(curr_state);
        
		double timeRand = (double) rand()/RAND_MAX;
		tau = log(1/timeRand)/(prop_sum*.001);
		t += tau;
		
		writeData(t, i, j, defuse_direction);
		
		if(animate) {
		    writeFile(x, n, reaction_count);
		}
		
		reaction_count++;
	}
	
	//printf("FINAL GRID\n");
	//printGrid(x, n);
    
    /*
    printf("PROP FUNCTIONS\n");
	for(j =0; j < n; j++){
			for(i = 0; i < n; i++){
				printf("%lf\t",prop_func[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	*/
	
	int ideal_pop = getIdealPop(x, n);
	printf("Number of Reactions: %d\n", reaction_count);
	printf("Ideal Population for Cell: ~%d\n", ideal_pop);
    
    /*
	for(unsigned i = 0; i < history.size(); i++) {
	    struct States temp_state = history.at(i);
        printf("i : %d\n", temp_state.i);
        printf("j : %d\n", temp_state.j);
        printf("direction : %d\n", temp_state.direction);
        printf("time : %lf\n", temp_state.t);
        printf("\n");
	}
    */
    
    return 0;
}

void writeFile(int **x, int n, int reaction) {
	int write_every = 500;
	FILE *file;
    char filename[20];
    int i, j;

    int max = getIdealPop(x, n) * 2;
    
    if(reaction % write_every == 0) {
        sprintf(filename, "output%09d.pgm", reaction);
        file = fopen(filename, "w");
        
	    if(file == NULL) {
		    printf("ERROR: Output file not found. Terminating...\n");
		    exit(1);
	    }
	    else {
		    fprintf(file, "%s\n%d %d\n%d\n", "P2", n, n, max);
		
		    for(i = 0; i < n; i++) {
			    for(j = 0; j < n; j++) {
				    int val = x[i][j];
				    if(val > max) {
				        val = max;
				    }
				    fprintf(file, "%d ", val);
			    }
			    fprintf(file, "\n");
		    }
		    fclose(file);
	    }
    }
}

void writeData(double t, int i, int j, int defuse_direction){
	FILE *fp;
	fp = fopen("data.txt", "a");
	fprintf(fp,"%lf\t(%d,%d)\t%d\n",t,i,j,defuse_direction);
	fclose(fp);
}

int getDefuseDirection(int i, int j, int n){
	int valid_directions[4];
	int index =0;
	int k;
	int direction;
	
	if(j - 1 >= 0){
		valid_directions[index] = 1;
		index++;
	}
	if(i - 1 >= 0){
		valid_directions[index] = 2;
		index++;
	}
	if(j + 1 <= n-1){
		valid_directions[index] = 3;
		index++;
	}
	if(i + 1 <= n-1){
		valid_directions[index] = 0;
	}
	else{
		index--;
	}
	
	direction = valid_directions[rand()%(index +1)];

	return direction;
}


int getDirectionNum(int i, int j, int n){
	int directions =0;
	if(j - 1 >= 0){
		directions++;
	}
	if(i - 1 >= 0){
		directions++;
	}
	if(j + 1 <= n-1){
		directions++;
	}
	if(i + 1 <= n-1){
		
		directions++;
	}
	return directions;
}


int getIdealPop(int **x, int n) {
    int i, j;
    
    int sum = 0;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            sum += x[i][j];
        }
    }
    
    //Just for verifying math...
    //printf("Total population is: %d\nTotal cells is: %d\n", sum, (n*n));
    
    return (sum / (n*n));
}

void printGrid(int **x, int n) {
    int i, j;
    
    for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			printf("%d\t",x[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

/*
function: creates a two dimensional array, simulating the grid our rections take place in
Parameters -
max_Pop: max population to generate within any given cell
numb_Cells: number of cells (randomly chosen within our grid) to generate population within (others have population of 0)
grid_Size: the size of our grid (grid_Size * grid_Size array)
returns : grid_Array (pointer to the populated/seeded 2d grid)
*/
int** seedCells(int max_Pop, int numb_Cells, int grid_Size){
	int i,j,k,numb;
	int** grid_Array = (int **)calloc((grid_Size),sizeof(int*)); 			
	for (i = 0; i <  grid_Size; i++){
		grid_Array[i] = (int *)calloc((grid_Size),sizeof(int));
	}//create 2d array, intiialize all to zero

	k=0;
	while(k < numb_Cells){
		for (i = 0; i <  grid_Size; i++){
			for(j = 0; j < grid_Size; j++){
				double prob = (double)rand()/RAND_MAX; //gen prob 0.0-1.0
				if(prob <= (double)numb_Cells/(grid_Size*grid_Size) && grid_Array[i][j] == 0){
					//numb = rand()%max_Pop;
					numb=max_Pop;
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
    /*
	for(i = 0; i < grid_Size; i++){
		for(j = 0; j < grid_Size; j++){
			printf("%d\t",grid_Array[i][j]);
		}//nested for
		printf("\n");
	}//for, print out seeded array
	printf("\n");
    */
	return grid_Array;
}//seedCells

int** seedFromInput(int n) {
	FILE *file;
	int i, j;
	
	int** array = (int **)calloc((n),sizeof(int*)); 			
	for (i = 0; i <  n; i++){
		array[i] = (int *)calloc((n),sizeof(int));
	}
	
	file = fopen("input.txt", "r");
	if(file == NULL) {
		printf("ERROR: Input file not found. Terminating...\n");
		exit(1);
	}
	else {
		for(i = 0; i < n; i++) {
			for(j = 0; j < n; j++) {
				fscanf(file, "%d", &array[i][j]);
			}
		}
		fclose(file);
	}
	
	return array;
}
