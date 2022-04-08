"""
The Travelling Salesman

A genetic algorithm is an optimization algorithm inspired by the theory of
evolution and natural selection. It is what is known as a heuristic algorithm,
which is an algorithm designed to solve a problem faster and more efficiently
than traditional methods. This decrease in time however involves a trade-off
with accuracy and/or precision, leading to a near-optimal solution. This
programme will implement a genetic algorithm that will be used to solve the
Travelling Salesman Problem. The Travelling Salesman Problem is as follows:
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
"""

# Import libraries/modules
from GA_number_sort import Gene, Chromosome, Population
import numpy as np
from matplotlib import pyplot as plt
from itertools import permutations


class City(Gene):
    """
    A class used to represent a city. Each City instance will be a gene in the
    genetic algorithm. The city object will have a name, an x and a y
    coordinate representing its latitudal and longitudal coordinates. The City
    class inherits from the Gene class from the GA_number_sort.
    file.

    ...

    Args
    ----------
    name : TYPE : string
        DESCRIPTION: The name of the City
    x : TYPE : int/float
        DESCRIPTION: Optional -> Default is 0
                     The latitudal coordinate of the City
    y : TYPE : int/float
        DESCRIPTION: Optional -> Default is 0
                     The longitudal coordinate of the City

    Attributes
    ----------
    value : TYPE : string
        DESCRIPTION: The name of the City
    x : TYPE : int/float
        DESCRIPTION: Optional -> Default is 0
                     The latitudal coordinate of the City
    y : TYPE : int/float
        DESCRIPTION: Optional -> Default is 0
                     The longitudal coordinate of the City
    
    Methods
    -------
    __sub__(other)
        DESCRIPTION: Overload sub magic method to calculate the distance
                     between 2 cities.
    """
    
    def __init__(self, name, x = 0, y = 0):
        """
        Parameters
        ----------
        name : TYPE : string
            DESCRIPTION: The name of the city
        x : TYPE : int/float
            DESCRIPTION: Optional -> Default is 0
                         The latitudal coordinate of the City
        y : TYPE : int/float
            DESCRIPTION: Optional -> Default is 0
                         The longitudal coordinate of the City

        Returns
        -------
        None.

        """
        # Store the city name as the instance variable value
        self.value = name
        # Store x and y-coordinates as instance variables
        self.x = x
        self.y = y
    
    def __sub__(self, other):
        """
        Overload the sub magic method to calculate the distance between two
        cities

        ...
        
        Parameters
        ----------
        ther : TYPE : object
            DESCRIPTION: Instance of City

        Returns
        -------
        TYPE - distance
            DESCRIPTION - Distance between two cities
        """
        # Calculate distance between city pair using the pythagorean theorem
        # R^2 = (x2 - x1)^2 + (y2 - y1)^2 => R = sqrt((x2 - x1)^2 + (y2 - y1)^2)
        distance = np.sqrt((self.x - other.x)**2 + (self.y - other.y)**2)
        # Return calculated distance between 2 city
        return distance


class Route(Chromosome):
    """
    A class used to represent a route. The route is a chromosome in the genetic
    algorithm. Each route contains a list of City objects (genes). Each City
    object has a unique value (name), and an x and y coordinate corresponding
    to its latitudal and longitudal coordinates. The Route class inherits from
    the Chromosome class from the GA_number_sort file.

    ...

    Args
    ----------
    randomised : Type : boolean
        Description: Optional -> Default is True
                     Determines whether or not the list of genes is shuffled

    Attributes
    ----------
    chromosometype : Type : object
            DESCRIPTION: Object which will be the individual chromosomes in the 
                         population
    genetype : Type : object
            DESCRIPTION: Object which will be the individual genes in the 
                         chromosome
    genes : TYPE : List
        DESCRIPTION: A list of City objects

    Methods
    -------
    getFitness()
        DESCRIPTION: Calculate the fitness of the Route instance
    getDistances()
        DESCRIPTION: Calculate the distances between each city pair
    """
    
    def __init__(self, randomised = True):
        """
        Parameters
        ----------
        randomised : TYPE : boolean
            DESCRIPTION: Shuffle the list of genes if True.
                         Leave unshuffled if False.

        Returns
        -------
        None.

        """
        # Dictionary containing each city with its x and y-coordinate
        cities = {
                "Dublin":{'x': 0, 'y': 0},
                "Wicklow":{'x': 14.03508772, 'y':-37.71929825},
                "Wexford":{'x': -12.71929825, 'y': -113.5964912},
                "Cork":{'x':-145.1754386, 'y': -160.0877193},
                "Galway":{'x':182.4561404, 'y': -7.456140351},
                "Athlone":{'x': -110.9649123, 'y': 7.456140351},
                "Belfast":{'x' :17.98245614, 'y':138.1578947},
                "Letterkenny":{'x': -94.73684211, 'y': 175.877193},
                "Portlaois":{'x': -67.10526316, 'y':-32.89473684},
                "Longford":{'x': -99.56140351, 'y':39.9122807}
                }
        
        # Initialise genetype instance variable to equal the City class
        self.genetype = City
        # Initialise chromosometype instance variable to equal the Route class
        self.chromosometype = Route
        # Initialise a list to store each City instance 
        self.genes = [self.genetype(city, location['x'], location['y']) for city, location in cities.items()]
        # If randomised is set to True, shuffle the list of genes
        if randomised is True:
            self.shuffle()
    
    def getFitness(self):
        """
        Calculate the fitness of the Route instance

        ...

        Returns
        -------
        TYPE - int/float
            DESCRIPTION - Fitness of the Route instance
        """
        # Calculate the distances between each City pair
        distances = self.getDistances()
        # Calculate the absolute distance between all cities
        fitness = sum([abs(distance) for distance in distances])
        # Return fitness
        return fitness
    
    def getDistances(self):
        """
        Calculate the distances between each city pair

        ...

        Returns
        -------
        TYPE - list
            DESCRIPTION - Distances between each city pair
        """
        # Calculate distances between each city pair
        distances = [gene - self.genes[i-1] for i, gene in enumerate(self.genes) if i > 0]
        # Return calculated city pair distances
        return distances


class RoutePopulation(Population):
    """
    A class used to represent a population of routes. This represents a
    population of chromosomes in the genetic algorithm. This population is a
    list of Route objects. Each Route instance is composed of a list
    of City objects. Each City object has a unique name value, an x and y
    coordinate. The RoutePopulation class inherits from the Population class
    from the GA_number_sort file.

    ...

    Args
    ----------
    N : TYPE : int
        DESCRIPTION: Optional -> Default is 10
                     Number of Chromosome instances to create and store

    Attributes
    ----------
    chromosometype : Type : object
        DESCRIPTION: Object which will be the individual chromosomes in the 
                     population
    genetype : Type : object
        DESCRIPTION: Object which will be the individual genes in the 
                     chromosome
    chromosomes : TYPE : List
        DESCRIPTION: A list of N Route objects
    
    Methods
    -------
    getOptimumFitness():
        DESCRIPTION: Calculate the fitness of a sorted chromosome, i.e. the
                     shortest route
    """
    
    def __init__(self, N = 10):
        """
        Parameters
        ----------
        N : TYPE : Int
            DESCRIPTION: Optional -> Default is 10
                         Number of Route instances to create and store

        Returns
        -------
        None.

        """
        # Initialise genetype instance variable to equal the City class
        self.genetype = City
        # Initialise chromosometype instance variable to equal the Route class
        self.chromosometype = Route
        # Initialise list of Route objects
        self.chromosomes = [self.chromosometype() for i in range(N)]
    
    def getOptimumFitness(self):
        """
        Calculate the fitness of a sorted chromosome, i.e. the shortest route

        ...

        Returns
        -------
        TYPE - float
            DESCRIPTION - Fitness of a sorted route
        """
        # Calculate every possible permutation of routes and store as a list
        all_permutations = list(permutations(Route().genes))
        # Initialise a new RoutePopulation
        population = RoutePopulation()
        # Initialise the populations list of chromosomes to be an empty list
        population.chromosomes = []
        # Iterate through all permutations of routes
        for i, permutation in enumerate(all_permutations):
            # Create of population of equal size to the number of permutations
            population.chromosomes.append(Route())
            # Set the genes (cities) of each route equal to each permutation
            population.chromosomes[i].genes = list(permutation)
        # Determine the fittest route out of every permutation
        fittest = population.getFittest(population.chromosomes)
        # Calculate optimum fitness and return
        optimum_fitness = fittest.getFitness()
        return optimum_fitness


# Only run the following code if running the file directly
if __name__ == '__main__':
    # Test / Demonstrate functionality
    
    # City class tests
    print('############### City class tests ###############\n')
    # Initialise City instances and print to console
    city1 = City(name = 'Dublin', x = 1, y = 2)
    print(f'City 1: {city1}')
    
    # Test getValue() function and print to console
    city2 = City('Galway', 3, 4)
    city3 = City('Dublin', 1, 2)
    print(f'City 2: {city2.getValue()}')
    print(f'City 3: {city3.getValue()}\n')
    
    # Test eq operator and print to console
    print(f'City 1 == City 2: {city1 == city2}')
    print(f'City 2 == City 3: {city2 == city3}')
    print(f'City 1 == City 3: {city1 == city3}\n')
    
    # Test sub operator and print to console
    print(f'City 1 - 2 City 2: {city1 - city2}\n\n')
    
    
    # Route class tests
    print('############### Route class tests ###############\n')
    
    # Initilise unshuffled Route instance, test getDistances() and getFitness()
    route_unshuffled = Route(False)
    print(f'Unshuffled route: {route_unshuffled}')
    print(f'Unshuffled route distances: {route_unshuffled.getDistances()}')
    print(f'Unshuffled route fitness: {route_unshuffled.getFitness()}\n')
    
    # Initilise shuffled Route instance, test getFitness(), and print
    route_shuffled = Route()
    print(f'Shuffled route: {route_shuffled}')
    print(f'Shuffled route fitness: {route_shuffled.getFitness()}\n')
    
    # Test swap, slide, and rotate functions
    print(f'Unshuffled route: {Route(False)}')
    print(f'Swap mutation: {Route(False).swap()}')
    print(f'Slide mutation: {Route(False).slide()}')
    print(f'Rotate mutation: {Route(False).rotate()}\n')
    
    # Test mutate() function
    print(f'Unshuffled route: {Route(False)}')
    print(f'Mutate (swap) mutation: {Route(False).mutate(1, "swap")}')
    print(f'Mutate (slide) mutation: {Route(False).mutate(1, "slide")}')
    print(f'Mutate (rotate) mutation: {Route(False).mutate(1, "rotate")}\n')
    
    # Test crossover functionality
    print(f'route_unshuffled + route_unshuffled = {route_unshuffled + route_unshuffled}\n')
    print(f'route_unshuffled + route_shuffled = {route_unshuffled + route_shuffled}\n')\
    
    # Test contains functionality
    dublin = City('Dublin', 0, 0)
    london = City('London', 5, 5)
    print(f'Is Dublin in the route? {dublin in route_shuffled}')
    print(f'Is London in the route? {london in route_shuffled}\n\n')
    
    
    # RoutePopulation class tests
    print('############### RoutePopulation class tests ###############\n')
    
    # Initialise route population and print to console
    route_pop = RoutePopulation(10)
    print(f'Route population:\n{route_pop}')
    
    # Test getFittest() function
    print(f'Fittest route: {route_pop.getFittest(route_pop.chromosomes)}\n')
    
    # Test tournamentSelect() function
    print(f'Tournament select result: {route_pop.tournamentSelect()}\n')
    
    # Test evolve() function
    print(f'Evolve the route population 1 generation:\n{route_pop.evolve()}')
    
    # Test getOptimumFitness() function
    print(f'Fitness of a sorted route: {route_pop.getOptimumFitness()}\n')
    
    # Test / demonstrate solve() function
    population = RoutePopulation()
    # Specify number of generations to run simulation
    generations = 100
    # Specify mutation rate
    rate = 0.5
    # Run solve() for n generations and store results
    pop_evolution = population.solve(generations, rate)
    
    # Get the fitness of each route in each generation
    gen_fitness = [[route.getFitness() for route in gen.chromosomes] for gen in pop_evolution]
    # Get the fitness of the fittest route in each generation
    gen_fittest_fitness = [min(gen) for gen in gen_fitness]
    # Get the average fitness of each generation
    gen_avg_fitness = [int(sum(gen)/len(gen)) for gen in gen_fitness]
    
    # Get the final optimised solution and print to console
    optimised_solution = pop_evolution[-1].getFittest(pop_evolution[-1].chromosomes)
    print(f'Optimised final solution: {optimised_solution}')
    # Get the fitness of the final optimised solution and print to console
    final_fitness = gen_fittest_fitness[-1]
    print(f'Final solution fitness: {final_fitness}')
    # Calculate number of generations to find optimum solution and print to console
    num_generations = len(pop_evolution)
    print(f'Generations to find solution: {num_generations}')
    
    # Pyplot graph
    # Set x-axis label
    plt.xlabel('Generation')
    # Set y-axis label
    plt.ylabel('Fitness')
    # Set graph title
    plt.title('Generational Fitness')
    # Plot the fittest route fitness and average fitness versus generation
    plt.plot(range(num_generations), gen_fittest_fitness, label='Fittest route per generation')
    plt.plot(range(num_generations), gen_avg_fitness, label='Average fitness per generation')
    # Display legend
    plt.legend(loc='best')
