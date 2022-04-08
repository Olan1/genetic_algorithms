# genetic_algorithms

A genetic algorithm is an optimization algorithm inspired by the theory of
evolution and natural selection. It is what is known as a heuristic algorithm,
which is an algorithm designed to solve a problem faster and more efficiently
than traditional methods. This decrease in time however involves a trade-off
with accuracy and/or precision, leading to a near-optimal solution.

## GA_number_sort
This programme will implement a genetic algorithm that will be used to sort a
list of unique integers and will consist of 3 classes.

The first class, Gene, represents a single gene. It will take a single value as
an argument which will be stored as that genes value using the class
constructor. This class will overload the str magic method to return the genes
value, and the eq magic method to determine if the values of 2 different genes
equate. The gene class will also contain a method which returns the
genes value.

The second class, Chromosome, represents a group of genes. This class will take
the arguments length, which determines the number of genes it contains, and
randomised, which determines whether the genes are shuffled or not. The class
will construct a list of genes with values of 1 to length using its constructor
function. The str magic method will be overloaded to return the value of each
gene in the genes list. The contains magic method will be overloaded to
determine if a gene is already present in a chromosome. The class will also
contain a shuffle function, which shuffles the list of genes, a swap function
which picks two random genes in the list and swaps their locations, a slide
function which will select a random segment from the list of genes and move it
to another random location in the list, a rotate function which will select a
random segment from the list of genes and reverse its order, a mutate function
which will mutate a chromosome if a randomly generated number is below a
specified threshold using either the swap, slide, or rotate mutation functions,
and finally a getFitness function which determines the fitness of the
chromosome. The add magic method will be overloaded to add together two parent
chromosomes to create a new child chromosome.

The third class, Population, represents a group of chromosomes. This class will
take the arguments N, which is the number of chromosomes in the population, and
length, the number of genes in each chromosome. In the class constructor, a
list of chromosomes will be initialised and stored. The str magic method will
be overloaded to return the values of each gene, in each chromosome, in the
population. It will contain a getFittest function which will determine the
fittest chromosome in a list of chromosomes, and a tournamentSelect function,
which will randomly pick a subset of n chromosomes from the entire population
of chromosomes, and select the fittest chromosome from this subset. It will
contain an evolve function which will create a new generation of chromosomes
from the existing generation, and finally a solve function which evolve the
population over a specified number of generations, returning a list containing
the population instance for each generation.

The functionality of these classes will be demonstrated and tested.

## the_travelling_salesman

This programme will implement a genetic algorithm that will be used to solve
the Travelling Salesman Problem. The Travelling Salesman Problem is as follows:
Given a list of cities and the distances between each pair of cities, what is
the shortest route that can be taken so that each city is visited once?

This programme will consist of 3 classes.

The first class, City, represents a single gene. It will take a name, an x and
a y coordinate as arguments which will all be stored as instance variables
using the class constructor. This class will inherit from the Gene class in the
GA_number_sort file. A description of the Gene class can be found in this file.
The sub magic method will be overloaded to calculate the distance between 2
cities.

The second class, Route, represents a group of City objects, or genes. This
class will inherit from the Chromosome class in the GA_number_sort file, a
description of which can be found in this file. This class will take the
argument randomised, a boolean which determines whether the list of genes is
shuffled or not. The class constructor will create a list of City instances
with a specific name, x and y coordinates using a provided dictionary within
its constructor function. This class will override the Chromosome getFitness
function to determine the fitness of each route. This class will also have a
function getDistances which will calculate the distances between each city pair.

The third class, RoutePopulation, contains a group of Routes. This represents
the population of chromosomes in the genetic algorithm. This class will inherit
from the Population class in the GA_number_sort file, a description of which
can be found in this file. This class will take the argument N, which is the
number of Routes, or chromosomes, in the population. In the class constructor,
a list of Routes (chromosomes) will be initialised and stored. This class will
override the Population getOptimumFitness function to determine the fitness of
a sorted route.

The functionality of these classes will be demonstrated and tested below.
