# ExpSustGRNs

## patrones.cc 
This library includes the source code on how the number of the distinct feedback loops patterns in the GRNs were calculated. That includes functions that count the number of self-activations, self-inhibitions, repricocal activations, reciprocal inhibitions, activations-inhibitions, and the 4 distinct kinds of 3-gene feedback loops. There are functions that count the total of positive or negative feedback loops in a GRN, as well as functions that calculate the number of genes on which cycles of positive interactions and mutual inhibitions converge (directly or indirectly). It also includes functions built to read a GRN from a source textfile and assign it to an object. This is crucial when analyzing the sample's random networks.

## tpanalysis.cc
It reads a network from a textfile ("originalnetwork.txt") which contains all the GRN interactions and calls the "patrones.cc" library functions to perform the patterns search. It returns the results in distinct textfiles ("act_inh_original.txt", "3ciclos_original.txt", "conjugados_original.txt", "total_original.txt").

## rnd\_tpanalysis.cc
It reads iteratively the sample's random networks (1,000) obtained from the original GRN and conducts the pattern analysis in each and every one of them. This is achieved with the same approach as in "tpanalysis.cc". In the end, it calculates the average and standard error of the 1,000-network sample values and returns these results in distinct textfiles ("act_inh_avgrandom.txt", "3ciclos_avgrandom.txt", "conjugados_avgrandom.txt", "total_avgrandom.txt"). In addition to this, the program returns other textfiles where the individual results for each network of the sample are contained ("act_inh_random.txt", "3ciclos_random.txt", "conjugados_random.txt", "total_random.txt").

## rng.R
It contains the source code of the functions that were developed to obtain the 1,000-network samples for each developmental GRN model in our dataset. We designed a modified version of the configuration model algorithm to generate networks that mantain the original GRN topological properties (size, number of both positive and negative interactions, degree distribution). The original GRN is read from the "originalnetwork.txt" file. Since all our GRN source files are named the same, it is crucial to have them in distinct directories. These are then loaded into the 'fuente' variable in the library. After running our program, it returns each of the randomized networks in individual files ("network1.txt", "network2.txt", ... "network1000.txt"). These networks are then examined with the "rnd_tpanalysis.cc" program.

