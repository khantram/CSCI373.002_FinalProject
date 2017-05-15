#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    int comm_sz = 4;
    int n = 8;
    int local_n = n/comm_sz;
    int i, j;
    
    double t = 0;
    double read_time = 0;
    double time_step = 10;
    double TIME_MAX = 3000;
    
    char final_file_name[20];
    char file_name[20];
    FILE *file, *final_file;
    int step_0 = 0, step_1 = 0, step_2 = 0, step_3 = 0, final_step = 0;;
    
    int grid[n][n];
    
    for(t = 0; t < TIME_MAX; t += time_step) {
        sprintf(final_file_name, "final_%06d.txt", final_step);
        final_file = fopen(final_file_name, "w");
        if(final_file == NULL) {
            printf("ERROR: %s not found. Terminating...\n", final_file_name);
	        exit(1);
        }
        else {
            
            //000000000000000000000000000000000000000000000000000000000000000000000000000
            while(1) {
                sprintf(file_name, "rank%02d_%06d.txt", 0, step_0);
                file = fopen(file_name, "r");
                
                if(file == NULL) {
                    printf("ERROR: %s not found. Terminating...\n", file_name);
	                exit(1);
                }
                else {
                    fscanf(file, "%lf", &read_time);
                }
                
                if(read_time < t) {
                    for(i = 0; i < n; i++) {
                        for(j = 0; j < 2; j++) {
                            fscanf(file, "%d", &grid[i][j]);
                        }
                    }
                }
                else {
                    fclose(file);
                    break;
                }
                fclose(file);
                step_0++;
            }
            
            
            
            /*
            //1111111111111111111111111111111111111111111111111111111111111111111111111111
            while(1) {
                sprintf(file_0_name, "rank%02d_%06d.txt", 1, step_1);
                file_1 = fopen(file_1_name, "r");
                
                if(file_0 == NULL) {
                    printf("ERROR: %s not found. Terminating...\n", file_1_name);
	                exit(1);
                }
                else {
                    fscanf(file_1, "%lf", &read_time);
                }
                
                
                if(read_time < t) {
                    for(i = 0; i < n; i++) {
                        for(j = 0; j < local_n; j++) {
                            fscanf(file_1, "%d", &grid[(i][1*local_n*j]);
                        }
                    }
                }
                else {
                    break;
                }
                
                step_1++;
            }
            fclose(file_1);
            
            //2222222222222222222222222222222222222222222222222222222222222222222222222222
            while(1) {
                sprintf(file_2_name, "rank%02d_%06d.txt", 2, step_2);
                file_2 = fopen(file_2_name, "r");
                
                if(file_2 == NULL) {
                    printf("ERROR: %s not found. Terminating...\n", file_2_name);
	                exit(1);
                }
                else {
                    fscanf(file_2, "%lf", &read_time);
                }
                
                if(read_time > t) {
                    step_2--;
                    break;
                }
                
                step_2++;
            }
            
            sprintf(file_2_name, "rank%02d_%06d.txt", 2, step_2);
            file_2 = fopen(file_2_name, "r");
            fscanf(file_2, "%*[^\n]\n");
            for(i = 0; i < n; i++) {
                for(j = (2*local_n); j < (3*local_n); j++) {
                    fscanf(file_2, "%d", &grid[i][j]);
                }
            }
            
            fclose(file_2);
            
            //3333333333333333333333333333333333333333333333333333333333333333333333333333333
            while(1) {
                sprintf(file_3_name, "rank%02d_%06d.txt", 3, step_3);
                file_3 = fopen(file_3_name, "r");
                
                if(file_3 == NULL) {
                    printf("ERROR: %s not found. Terminating...\n", file_3_name);
	                exit(1);
                }
                else {
                    fscanf(file_3, "%lf", &read_time);
                }
                
                if(read_time > t) {
                    step_3--;
                    break;
                }
                
                step_3++;
            }
            
            sprintf(file_3_name, "rank%02d_%06d.txt", 3, step_3);
            file_3 = fopen(file_3_name, "r");
            fscanf(file_3, "%*[^\n]\n");
            for(i = 0; i < n; i++) {
                for(j = (3*local_n); j < (4*local_n); j++) {
                    fscanf(file_3, "%d", &grid[i][j]);
                }
            }
            
            fclose(file_3);
            */
            //================================================================================
            printf("AAAAAAAAAAAAAAAAAAAHHHHH!!! t = %lf\n", t);
            
            for(i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    fprintf(final_file, "%d\t", grid[i][j]);
                }
                fprintf(final_file, "\n");
            }
        }
        final_step++;
        return 0;
    }
}
