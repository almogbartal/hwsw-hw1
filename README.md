# Prompt

generate me a kuramoto model firefly synchronization simulation in c that takes as an input parameter the number of fireflies to synchronize. every firefly should be specified by 2 floating point
  coordinates and they live on a [0,1]x[0,1] square. the fireflies should be placed randomly and uniformly inside the square, and the distance should be euclidean distance. different parameters, like
  coupling strength etc. can be defined arbitrarily or by looking at other projects. on the screen you should print every 1000 iterations mean phase, median phase and std.

  write it in a single .c file. dont optimize the program -- rather make it naive and very simple in terms of performance/optimization.

what can be a criterion to see that the simulation converges to synchronization? what can be printed on the cli to see this without visually rendering all the fireflies? implement it

also the firefly phases should also be initialized as uniformly random on [0:2pi]

make it stop whenever the simulation has converged


#final prompt

"please write a simple, single-threaded C program that simulates firefly synchronization using the Kuramoto model for 50,000 fireflies. each firefly lives on a [0,1]x[0,1] square and should be represented by a struct holding its ID, x and y coordinates, phase, and frequency. the coupling should be based on the Euclidean distance by comparing each firefly to all other fireflies in the simulation."
