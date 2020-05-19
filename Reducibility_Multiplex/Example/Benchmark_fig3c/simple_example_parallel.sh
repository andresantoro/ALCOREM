#!/bin/bash

#Launch the reducibility code using 10 iterations, parallel computation (4 cores) and the verbose activated
./../../reducibility_complexity nodes_unique.txt list_multiplex.txt 10 -p 4 -v
#"**   OUTPUT: by default the algorithm returns the following info:                       **\n"
#"**   partition of the multiplex --> Number of layers, KC multiplex,                     **\n"
#"**   # links multiplex, KC aggregate, # links aggregate, complexity, quality function   **\n"
