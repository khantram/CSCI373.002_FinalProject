## CSCI 373.002 - Final Project (continued work)
### Parallelizing Stochastic Diffusion using MPI

Initially a group project with Nathan Breedlove and Eliot Dixon

@Kenny Tram hoping to continue the project


###### Notes
- spacial_stoch.c - Serial Implementation
- MPI_spacial_stoch.c - Parallel Implementation
- Notable animations of stochastic diffusion (generated in serial) included

###### To-do List
- Kill signals to assure that all processors "finish" when appropriate
- combine_output.c - A program to combine the output of each processor; synchronized by time
