# Prompt

generate me a kuramoto model firefly synchronization simulation in c that takes as an input parameter the number of fireflies to synchronize. every firefly should be specified by 2 floating point coordinates and they live on a [0,1]x[0,1] square. the fireflies should be placed randomly and uniformly inside the square, and the distance should be euclidean distance. different parameters, like coupling strength etc. can be defined arbitrarily or by looking at other projects. on the screen you should print every 1000 iterations mean phase, median phase and std. all fireflies should be in the same frequency.

write it in a single .c file. 

what can be a criterion to see that the simulation converges to synchronization? what can be printed on the cli to see this without visually rendering all the fireflies? implement it and print this parameter on the screen beside all the other parameters

also the firefly phases should also be initialized as uniformly random on [0:2pi]
make it stop whenever the simulation has converged. initialize the random number generator with value 42.

for every firefly make a struct with id, x,y and phase.

do it single threaded.

to change fireflies phase you should compare each firefly against all the rest
