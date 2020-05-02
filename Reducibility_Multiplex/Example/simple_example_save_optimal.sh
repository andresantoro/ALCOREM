
#!/bin/bash

#Launch the reducibility code using 100 iterations and 4 cores to compute it.
#Verbose activated & save the optimal aggregation in the current directory as edge lists!

./../reducibility_complexity nodes_unique.txt list_multiplex.txt 100 -p 4 -v -s

