
#!/bin/bash

#Launch the reducibility code using 100 iterations and 4 cores to compute it.
#Verbose activated & save all the aggregations in the current directory as edge lists!
#IMPORTANT: This option will create a lot of files in the current directory.

./../reducibility_complexity nodes_unique.txt list_multiplex.txt 100 -p 4 -v -w

