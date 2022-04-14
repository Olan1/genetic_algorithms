"""
GA Number Sort

A genetic algorithm is an optimization algorithm inspired by the theory of
evolution and natural selection. It is what is known as a heuristic algorithm,
which is an algorithm designed to solve a problem faster and more efficiently
than traditional methods. This decrease in time however involves a trade-off
with accuracy and/or precision, leading to a near-optimal solution. This
programme will implement a genetic algorithm that will be used to sort a list
of unique integers and will consist of 3 classes.

The first class, Gene, represents a single gene. It will take a single value as
an argument which will be stored as that genes value using the class
constructor. This class will overload the str magic method to return the genes
value, the eq magic method to determine if the values of 2 different genes
equate, and the sub magic method to determine the value of one gene subtracted
from the other. The gene class will also contain a method which returns the
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

The functionality of these classes will be demonstrated and tested below.
"""

# Import libraries/modules
import random as rm
import copy
from matplotlib import pyplot as plt


class Gene:
    """
    A class used to represent a gene. The Gene object has a numerical value

    ...

    Args
    ----------
    num : TYPE : int
        DESCRIPTION: Optional -> Default is 1
                     The numerical value of the gene

    Attributes
    ----------
    value : TYPE : int
        DESCRIPTION: Optional -> Default is 1
                     The numerical value of the gene

    Methods
    -------
    __eq__(other)
        DESCRIPTION: Overload the equality magic method the check if the values
                     of 2 gene objects are equal
    __str__()
        DESCRIPTION: Overload the string magic method to return the value of
                     the gene object as a string
    getValue()
        DESCRIPTION: Get the value of the Gene instance
    """

    def __init__(self, num = 1):
        """
        Parameters
        ----------
        num : TYPE : int
            DESCRIPTION: Optional -> Default is 1
                         The numerical value of the gene

        Returns
        -------
        None.

        """
        # Initialise instance variable value to equal provided num argument
        self.value = num

    def __eq__(self, other):
        """
        Overload the equality magic method to check if the values of 2 gene
        objects are equal

        ...

        Parameters
        ----------
        other : TYPE : object
            DESCRIPTION: Instance of Gene class

        Returns
        -------
        TYPE - boolean
            DESCRIPTION: True if values of each Gene instance are equal.
                         False if values of each Gene instance are unequal.
        """
        # Return True if both gene values are equal, else return False
        return self.value == other.value

    def __str__(self):
        """
        Overload the str magic method to return the value of the gene object as
        a string

        ...

        Returns
        -------
        TYPE - string
            DESCRIPTION - The value of the Gene instance as a string
        """
        # Return the value of the Gene instance as a string
        return str(self.value)

    def getValue(self):
        """
        Get the value of the Gene instance

        ...

        Returns
        -------
        TYPE - int
            DESCRIPTION - The value of the Gene instance
        """
        # Return the value of Gene instance
        return self.value
    
    def __sub__(self, other):
        """
        Overload the subtraction magic method to deduct the value of one gene
        object from another

        ...

        Parameters
        ----------
        other : TYPE : object
            DESCRIPTION: Instance of Gene class

        Returns
        -------
        TYPE - int
            DESCRIPTION - The value of one Gene instance subtracted from the
                          other
        """
        # Deduct value of one Gene instance from the other and return new value
        return self.value - other.value


class Chromosome:
    """
    A class used to represent a chromosome. A chromosome is a list of Gene
    objects. Each Gene object has a unique numerical value.

    ...

    Args
    ----------
    length : TYPE : int
        DESCRIPTION: Optional -> Default is 10
                     Number of Gene instances to create and store
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
        DESCRIPTION: A list of Gene objects with values ranging from 1 to
                     the specified length

    Methods
    -------
    __str__()
        DESCRIPTION: Overload the str magic method to return the values of each
                     Gene object in the genes list as a string
    __contains__(gene)
        DESCRIPTION: Overload the contains magic method to check if a Gene
                     instance is already present in the list of genes
    shuffle()
        DESCRIPTION: Randomly shuffles the list of genes
    swap()
        DESCRIPTION: Swaps the positions of 2 random elements in the genes list
    getFitness()
        DESCRIPTION: Calculate the fitness of the Chromosome instance
    slide()
        Description: Selects a random segment from the genes list and moves it
                     to a random location in the genes list
    rotate()
        Description: Reverse the order of a random segment of genes in the
                     genes list
    mutate(threshold, mutation)
        Description: Mutate the chromosome if a randomly generated number is
                     greater than the threshold
    __add__(parent2)
        Description: Overload the add magic method to produce a child chromosome
                     from the addition of two parent chromosomes
    """

    def __init__(self, length = 10, randomised = True):
        """
        Parameters
        ----------
        length : TYPE : Int
            DESCRIPTION: Optional -> Default is 10
                         Number of Gene instances to create and store in list
        randomised : Type : boolean
            Description: Optional -> Default is True
                         Determines if the list of genes is shuffled or not

        Returns
        -------
        None.

        """
        # Initialise chromosometype instance variable to equal Chromosome class
        self.chromosometype = Chromosome
        # Initialise genetype instance variable to equal Gene class
        self.genetype = Gene
        # Initialise list of Gene objects, providing Gene values 1 to length
        self.genes = [self.genetype(i+1) for i in range(length)]
        # If randomised is set to True, shuffle the list of genes
        if randomised is True:
            self.shuffle()

    def __str__(self):
        """
        Overload the str magic method to return the values of each Gene object
        in the genes list as a string

        ...

        Returns
        -------
        TYPE - str
            DESCRIPTION - List of instance values of Gene class

        """
        # Create list of each Gene objects value in the genes list
        # Return as a string
        return str([gene.value for gene in self.genes])

    def __contains__(self, gene):
        """
        Overload the contains magic method to check if a Gene instance is
        already present in the list of genes

        ...

        Parameters
        ----------
        gene : TYPE : object
            DESCRIPTION: Instance of Gene class

        Returns
        -------
        TYPE - boolean
            DESCRIPTION - True if gene appears in the list of values.
                          False if gene doesn't appear in the list of values.
        """
        # Return True if the provided gene is in the list of genes, else False
        return gene in self.genes

    def shuffle(self):
        """
        Randomly shuffles the list of genes

        ...

        Returns
        -------
        TYPE - Object
            DESCRIPTION - Updated instance of self
        """
        # Shuffle list of genes
        rm.shuffle(self.genes)
        # Return updated instance of self
        return self

    def swap(self):
        """
        Swaps the positions of 2 random elements in the genes list

        ...

        Returns
        -------
        TYPE - Object
            DESCRIPTION - Updated instance of self
        """
        # Pick 2 unique elements from list of genes
        choice1, choice2 = rm.sample(self.genes, 2)
        # Get the index of choice1 and choice 2 from the list of genes
        index1 = self.genes.index(choice1)
        index2 = self.genes.index(choice2)
        # Swap locations of choice1 and choice2
        self.genes[index1] = choice2
        self.genes[index2] = choice1
        # Return updated instance of self
        return self
    
    def slide(self):
        """
        Selects a random segment from the genes list and moves it to a random
        location in the genes list
        
        ...

        Returns
        -------
        TYPE - object
            DESCRIPTION - Updated instance of self

        """
        # Select 2 random indices from the genes list
        indices = rm.sample(range(len(self.genes)), 2)
        # Get the upper and lower indices and assign to variables
        min_i, max_i = sorted(indices)
        # Slice the genes list from min_i to max_i
        segment = self.genes[min_i : max_i]
        # Create a tuple of all possible indices the segment can be moved to
        start_range = tuple(range(0, len(self.genes) - len(segment)))
        # Select a random index to insert the segment
        insert_i = rm.choice(start_range)
        # Remove the segment genes from the genes list 
        self.genes = [gene for gene in self.genes if gene not in segment]
        # Insert the segment at the specified index
        self.genes[insert_i : insert_i] = segment
        # Return updated instance of self
        return self
    
    def rotate(self):
        """
        Reverse the order of a random segment of genes in the genes list
        
        ...

        Returns
        -------
        TYPE - object
            DESCRIPTION - Updated instance of self

        """
        # Select 2 random indices from the genes list
        indices = rm.sample(range(len(self.genes)), 2)
        # Get the upper and lower indices and assign to variables
        min_i, max_i = sorted(indices)
        # Reverse order of genes between min_i and max_i in the genes list
        self.genes[min_i : max_i] = reversed(self.genes[min_i : max_i])
        # Return updated instance of self
        return self
    
    def mutate(self, threshold = 0.5, mutation = 'swap'):
        """
        Mutate the chromosome if a randomly generated number is greater than
        the threshold
        
        ...

        Parameters
        ----------
        threshold : TYPE : float
            DESCRIPTION - Optional -> Default is 0.5
                          The rate at which mutations will occur
        mutation : TYPE : string
            DESCRIPTION - Optional -> Default is 'swap'
                          Corresponds to the type of mutation to perform
                          1. 'swap' - Swap the locations of 2 random genes in
                                      the genes list.
                          2. 'slide' - Move a random segment of genes to a
                                       random location in the genes list.
                          3. 'rotate' - Reverse the order of a random segment
                                        of genes in the genes list.

        Returns
        -------
        TYPE - object
            DESCRIPTION - Updated instance of self

        """
        # Dictionary containing the 3 possible mutation functions
        mutations = {'swap': self.swap,
                     'slide': self.slide,
                     'rotate': self.rotate}
        # Select the specified mutation function from the mutations dictionary
        mutate = mutations[mutation]
        # Generate a random number between 0 and 1
        r = rm.random()
        # If the random number is below the threshold, mutate the chromosome
        if r < threshold:
            mutate()
        # Return an updated instance of self
        return self

    def getFitness(self):
        """
        Calculate the fitness of the Chromosome instance

        ...

        Returns
        -------
        TYPE - int/float
            DESCRIPTION - Fitness of the Chromosome instance
        """
        # Initialise a fitness variable to 0
        fitness = 0
        # Iterate through each gene
        for i, gene in enumerate(self.genes):
            # Skip the first iteration
            if i == 0:
                continue
            # Calculate the absolute difference between the adjacent genes
            diff = abs(gene - self.genes[i-1])
            # Increment the fitness by the difference
            fitness += diff
        # Return the fitness
        return fitness
    
    def __add__(self, parent2):
        """
        Overload the add magic method to produce a child chromosome from the 
        addition of two parent chromosomes 
        
        ...

        Parameters
        ----------
        parent2 : TYPE : object
            DESCRIPTION - Instance of chromosome object

        Returns
        -------
        TYPE : object
            DESCRIPTION - Updated instance of self

        """
        # Pick 2 random indices from the parent 1 (self) genes list
        indices = rm.sample(range(len(self.genes)), 2)
        # Get the upper and lower indices and assign to variables
        min_i, max_i = sorted(indices)
        # Slice parent 1 genes list from min_i to max_i
        segment = self.genes[min_i : max_i]
        # Create child chromosome instance
        child = self.chromosometype(len(self.genes))
        # Replace child genes with parent2 genes that aren't in parent1 segment
        child.genes = [gene for gene in parent2.genes if gene not in segment]
        # Insert the segment at the specified index
        child.genes[min_i : min_i] = segment
        # Return updated instance of self
        return child


class Population:
    """
    A class used to represent a population of chromosomes. This population is a
    list of Chromosome objects. Each Chromosome instance is composed of a list
    of Gene objects. Each Gene object has a unique numerical value.

    ...

    Args
    ----------
    N : TYPE : int
        DESCRIPTION: Optional -> Default is 10
                     Number of Chromosome instances to create and store
    length : TYPE : int
        DESCRIPTION: Optional -> Default is 10
                     Number of Gene objects in each Chromosome instance

    Attributes
    ----------
    chromosometype : Type : object
        DESCRIPTION: Object which will be the individual chromosomes in the 
                     population
    genetype : Type : object
        DESCRIPTION: Object which will be the individual genes in the 
                     chromosome
    chromosomes : TYPE : List
        DESCRIPTION: A list of N Chromosome objects
    optimalfitness : TYPE : float
        DESCRIPTION: The fitness of a sorted chromosome

    Methods
    -------
    __str__():
        DESCRIPTION: Overloads the str magic method to return the values of
                     each Gene, in each Chromosome, in the population, as a
                     string
    getOptimumFitness():
        DESCRIPTION: Calculate the fitness of a sorted chromosome
    getFittest(chromosomes):
        DESCRIPTION: Determine the fittest Chromosome instance from a provided
                     list of Chromosomes
    tournamentSelect(n):
        DESCRIPTION: Create a subset of n random chromosomes from the list of
                     chromosomes and select the fittest chromosome in this
                     subset
    evolve(rate, n):
        DESCRIPTION: Create a new population generation of chromosomes from
                     the current generation
    solve(generations, rate, n):
        DESCRIPTION: Evolve the population over a specified number of
                     generations, with a specified mutation rate. Return a list
                     containing the population instance from each generation
                     once the most optimum solution has been found or the
                     number of generations has been reached
    """

    def __init__(self, N = 10, length = 10):
        """
        Parameters
        ----------
        N : TYPE : Int
            DESCRIPTION: Optional -> Default is 10
                         Number of Chromosome instances to create and store
        length : TYPE : Int
            DESCRIPTION: Optional -> Default is 10
                         Number of Gene objects in each Chromosome instance

        Returns
        -------
        None.

        """
        # Initialise chromosometype instance variable to equal Chromosome class
        self.chromosometype = Chromosome
        # Initialise genetype instance variable to equal Gene class
        self.genetype = Gene
        # Initialise list of Chromosome objects
        self.chromosomes = [self.chromosometype(length) for i in range(N)]
        # Store optimum fitness as an instance variable
        self.optimumfitness = self.getOptimumFitness()

    def __str__(self):
        """
        Overloads the str magic method to return the values of each Gene, in
        each Chromosome, in the population, as a string

        ...

        Returns
        -------
        pop_string : TYPE : str
            DESCRIPTION - The values of each Gene, in each Chromosome, in the
                          population, as a string

        """
        # Initialise an empty string variable to store the population values
        pop_string = ""
        # Iterate through each chromosome in the population
        for chromosome in self.chromosomes:
            # Create list of values for each Gene, in the current Chromosome
            chrom_vals = [gene.value for gene in chromosome.genes]
            # Append the chromosome values to the string and create a new line
            pop_string += str(chrom_vals) + '\n'
        # Return the population values as a string
        return pop_string
    
    def getOptimumFitness(self):
        """
        Calculate the fitness of a sorted chromosome

        ...

        Returns
        -------
        TYPE - int
            DESCRIPTION - Fitness of a sorted chromosome
        """
        return Chromosome(randomised = False).getFitness()
    
    def getFittest(self, chromosomes):
        """
        Determine the fittest Chromosome instance from a provided list of
        Chromosomes

        ...

        Parameters
        ----------
        chromosomes : TYPE : list
            DESCRIPTION: List of Chromosome instances

        Returns
        -------
        TYPE - object
            DESCRIPTION - Fittest Chromosome instance
        """
        # Calculate the fitness for each chromosome and add to the list
        fitness = [chromosome.getFitness() for chromosome in chromosomes]
        # Get the index of the fittest chromosome in the fitness list
        fittest_index = fitness.index(min(fitness))
        # Select and return the fittest chromosome from the chromosomes list
        fittest = chromosomes[fittest_index]
        return fittest

    def tournamentSelect(self, n = 5):
        """
        Create a subset of n random chromosomes from the list of chromosomes
        and select the fittest chromosome in this subset.

        ...

        Parameters
        ----------
        n : TYPE : int
            DESCRIPTION: Optional -> Default is 5
                         Number of chromosomes to randomly select from the
                         population

        Returns
        -------
        TYPE - object
            DESCRIPTION - Fittest Chromosome instance
        """
        # Randomly select n chromosomes from the list of chromosome objects
        subset = rm.sample(self.chromosomes, n)
        # Select and return the fittest chromosome in the subset
        selected = self.getFittest(subset)
        return selected
    
    def evolve(self, rate = 0.5, n = 5):
        """
        Create a new population generation of chromosomes from current
        generation

        ...

        Parameters
        ----------
        rate : TYPE : float
            DESCRIPTION: Optional -> Default is 0.5
                         Rate that mutations will occur in child chromosomes
        n : TYPE : int
            DESCRIPTION: Optional -> Default is 5
                         Number of chromosomes to randomly select from the
                         population in the tournamentSelect() function

        Returns
        -------
        TYPE - object
            DESCRIPTION - Updated instance of self
        """
        # Initialise empty list to store the new generation of chromosomes
        new_population = []
        # Iterate through the number of chromosomes in the current population
        for i in range(len(self.chromosomes)):
            # Use tournament select to select 2 parent chromosomes
            parent1 = self.tournamentSelect(n)
            parent2 = self.tournamentSelect(n)
            # Create a child chromosome from the 2 parent chromosomes
            child = parent1 + parent2
            # Select a random mutation
            mutation = rm.choice(['swap', 'slide', 'rotate'])
            # Mutate child randomly
            child.mutate(rate, mutation)
            # Add the child chromosome to the new generation population
            new_population.append(child)
        # Update chromosomes to equal the new generation population
        self.chromosomes = new_population
        # Return updated instance of self
        return self
    
    def solve(self, generations = 100, rate = 0.5, n = 5):
        """
        Evolve the population over a specified number of generations, with a
        specified mutation rate. Return a list containing the population
        instance from each generation once the most optimum solution has been
        found or the number of generations has been reached

        ...

        Parameters
        ----------
        generations : TYPE : int
            DESCRIPTION: Optional -> Default is 100
                         Number of generations to run genetic algorithm
        rate : TYPE : float
            DESCRIPTION: Optional -> Default is 0.5
                         Rate that mutations will occur in child chromosomes
        n : TYPE : int
            DESCRIPTION: Optional -> Default is 5
                         Number of chromosomes to randomly select from the
                         population in the tournamentSelect() function

        Returns
        -------
        pop_evolution : TYPE : list
            DESCRIPTION: List of each generations population
        """
        # Initialise an empty list to store each generations population
        pop_evolution = []
        # Iterate through each generation
        for i in range(generations):
            # Create a new generation population of chromosomes
            self.evolve(rate)
            # Create a deep copy of the current population instance
            generation = copy.deepcopy(self)
            # Append the new generation population to the list
            pop_evolution.append(generation)
            # Get the fittest chromosome in the current generation
            fittest = generation.getFittest(generation.chromosomes)
            # Get the fitness of the fittest chromosome
            fitness = fittest.getFitness()
            # If the fittest chromosome fitness is equal to the optimal fitness
            if fitness <= self.optimumfitness:
                # Return the list of generations
                return pop_evolution
        # Return the list of generations
        return pop_evolution


# Only run the following code if running the file directly
if __name__ == '__main__':
    
    # Demonstrate / Test Functionality
    
    # Gene class tests
    print('\n########## Gene class tests ##########\n')
    
    # Initialise gene and print to console
    gene1 = Gene(4)
    print('Gene 1 value: ' + str(gene1))
    
    # Initialise gene, test getValue() function, and print to console
    gene2 = Gene(14)
    print('Gene 2 value: ' + str(gene2.getValue()))
    
    # Initialise gene and print to console
    gene3 = Gene(4)
    print('Gene 3 value: ' + str(gene3.getValue()) + '\n')
    
    # Test equality magic method
    print('Gene 1 == Gene 2: ' + str(gene1 == gene2))
    print('Gene 1 == Gene 3: ' + str(gene1 == gene3))
    print('Gene 2 == Gene 3: ' + str(gene2 == gene3) + '\n')
    
    
    # Chromosome class tests
    print('\n\n########## Chromosome class tests ##########\n')
    
    # Initialise chromosome with randomised set to False, and print to console
    sorted_chromosome = Chromosome(10, False)
    print('Sorted chromosome: ' + str(sorted_chromosome) + '\n')
    
    # Initialise chromosome with randomised set to default(True), and print to console
    randomised_chromosome = Chromosome(10)
    print('Randomised chromosome: ' + str(randomised_chromosome) + '\n')
    
    # Initialise chromosome with randomised set to False
    # Test shuffle() function and print to console
    shuffled_chromosome = Chromosome(10, False).shuffle()
    print('Shuffled chromosome: ' + str(shuffled_chromosome) + '\n')
    
    # Test contains magic method
    print(f'Is Gene 1 in the chromosome {randomised_chromosome} ? ' + str(gene1 in randomised_chromosome))
    print(f'Is Gene 2 in the chromosome {randomised_chromosome} ? ' + str(gene2 in randomised_chromosome) + '\n')
    
    # Test getFitness() function
    print('Sorted chromosome fitness: ' + str(sorted_chromosome.getFitness()))
    print('Shuffled chromosome fitness: ' + str(shuffled_chromosome.getFitness()) + '\n')
    
    # Test the swap() function and print to console
    print('Original chromosome: ' + str(Chromosome(10, False)))
    print('Swapped mutation: ' + str(Chromosome(10, False).swap()) + '\n')
    
    # Test slide() mutation function and print to console
    print('Original chromosome: ' + str(Chromosome(10, False)))
    print('Slide mutation: ' + str(Chromosome(10, False).slide()) + '\n')
    
    # Test rotate() mutation function and print to console
    print('Original chromosome: ' + str(Chromosome(10, False)))
    print('Rotate mutation: ' + str(Chromosome(10, False).rotate()) + '\n')
    
    # Test mutate() mutation function and print to console
    # Set threshold to 1 so mutation occurs each time
    print('Original chromosome: ' + str(Chromosome(10, False)))
    print('Swap mutation: ' + str(Chromosome(10, False).mutate(threshold = 1, mutation = 'swap')))
    print('Slide mutation: ' + str(Chromosome(10, False).mutate(threshold = 1, mutation = 'slide')))
    print('Rotate mutation: ' + str(Chromosome(10, False).mutate(threshold = 1, mutation = 'rotate')) + '\n')
    
    # Test add magic operator and print to console
    chrom_1 = Chromosome(10, False)
    chrom_2 = Chromosome(10)
    print('Chromosome 1: ' + str(chrom_1))
    print('Chromosome 2: ' + str(chrom_2))
    print('Chromosome 1 + Chromosome 2 = ' + str(chrom_1 + chrom_2))
    
    
    # Population class tests
    print('\n\n########## Population class tests ##########\n')
    
    # Initialise instance of Population with default arguments and print to console
    population = Population(10)
    print('Population:\n' + str(population) + '\n')
    
    # Test getFittest() function on total population
    pop_fittest = population.getFittest(population.chromosomes)
    print('Fittest chromosome in population: ' + str(pop_fittest))
    # Get fitness of the fittest chromosome
    fittest_fitness = pop_fittest.getFitness()
    print('Fitness: ' + str(fittest_fitness) + '\n')
    
    # Test the tournamentSelect() function and print result to the console
    t_select = population.tournamentSelect(5)
    print('tournamentSelect result: ' + str(t_select))
    # Calculate the fitness of the tournamentSelect result
    print('Fitness: ' + str(t_select.getFitness()) + '\n')
    
    # Test evolve() function and print to console
    evolved_pop = population.evolve()
    print('Evolve population 1 generation: \n' + str(evolved_pop) +'\n')
    
    # Test / demonstrate solve() function
    # Specify number of generations to run simulation
    generations = 100
    # Specify mutation rate
    rate = 0.5
    # Run solve() for n generations and store results
    pop_evolution = population.solve(generations, rate)
    
    # Get the fitness of each chromosome in each generation
    gen_fitness = [[chrom.getFitness() for chrom in gen.chromosomes] for gen in pop_evolution]
    # Get the fitness of the fittest chromosome in each generation
    gen_fittest_fitness = [min(gen) for gen in gen_fitness]
    # Get the average fitness of each generation
    gen_avg_fitness = [int(sum(gen)/len(gen)) for gen in gen_fitness]
    
    # Get the final optimised solution and print to console
    optimised_solution = pop_evolution[-1].getFittest(pop_evolution[-1].chromosomes)
    print(f'Optimised final solution: {optimised_solution}')
    # Get the fitness of the final optimised solution and print to console
    final_fitness = gen_fittest_fitness[-1]
    print(f'Optimised solution fitness: {final_fitness}')
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
    # Plot the fittest chromosome fitness and average fitness versus generation
    plt.plot(range(num_generations), gen_fittest_fitness, label='Fittest chromosome per generation')
    plt.plot(range(num_generations), gen_avg_fitness, label='Average fitness per generation')
    # Display legend
    plt.legend(loc='best')
