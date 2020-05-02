# ALCOREM-0.1


This is **ALCOREM**(**AL**gorithmic **CO**mplexity and **RE**ducibility of **M**ultiplex networks), a C implementation of the algorithm for computing the algorithmic complexity and reducibility of a multiplex networks.

In the folders there are two different C codes, which respectively implement the algorithms to compute the algorithmic complexity of a multiplex network and the one for reducing the multiplex network. The algorithms are based on the procedure described in the paper:

A.Santoro, V. Nicosia, Algorithmic Complexity of Multiplex Networks, https://arxiv.org/abs/1903.08049

The code has been tested on OSX and Linux distributions using gcc/clang/icc compilers.


### Dependencies

The algorithm makes use of the following libraries: 
 - GMP (The GNU Multiple Precision Arithmetic Library - https://gmplib.org/ )
 - OpenMP
 - zlib (https://www.zlib.net/)
 
Therefore, before using the code, make sure you have installed it.  Otherwise, you can install easily as:

On Debian based systems:
```
sudo apt-get install libgmp3-dev
sudo apt-get install libomp-dev
sudo apt-get install zlib1g-dev
```

On OSX through brew:
```
brew install gmp
brew install zlib
brew install libomp
```

Sample bash scripts to launch the code with different options can be found in the "Example/" directory, where you will also find a sample real-world data set. This is the 13 layers London tube, originally provided in: 

M. De Domenico, A. Solé-Ribalta, S. Gómez, and A. Arenas, Navigability of interconnected networks under random failures. PNAS 111, 8351-8356 (2014)

Please consider citing this paper if you use this data set in a scientific work.

Another synthetic example is also reported in the directory "Example/". This include the synthetic multiplex network made of 10 ER random graphs reported in Fig 3c of the paper, i.e. five pairs of identical layers.


## License

This project is licensed under the GNU GPLv3 license - see the [LICENSE](LICENSE) file for details

------ 

## Python Package

My plan is to provide a python wrapper to better integrate the C code with python scripts in the near future.

