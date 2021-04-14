def create_individual(bounds):
    '''
    Creates an individual with randomly chosen values from supplied parameter
    boundaries in array form (an n x 2 array where n is the number of
    parameters).
    Kind of old and not used in most things this module is used for.
    Should be modified and documented later.
    '''

    indiv = np.empty([len(bounds)])

    for n in range(0,len(bounds)):
        indiv[n] = random.uniform(bounds[n][0], bounds[n][1])

    return indiv

def create_population(size,bounds):

    '''
    Generates a population using the 'individual' function for each member of
    the population.
    Kind of old and not used in most things this module is used for.
    Should be modified and documented later.
    '''

    pop = np.empty([size, len(bounds)])

    for n in range(0,size):
        pop[n] = create_individual(bounds)

    return pop

def test_population(population,data,model):
    '''
    Tests each member of the population against the data. Returns a list of
    errors of each member of the population.
    Kind of old and not used in most things this module is used for.
    Should be modified and documented later.
    '''

    S = initialise_concs(model,data)

    errorlist = np.empty([len(population)])

    for n in range(0,len(population)):
        k = population[n]

        errorlist[n] = compare_model_data(k, data, S, model)

    return errorlist

def sort_population(population, errorlist):
    '''
    Sorts population based on an error list (low to high).
    Kind of old and not used in most things this module is used for.
    Should be modified and documented later.
    '''

    i = np.argsort(errorlist) # Get array indices from error list
    errorlist = errorlist[i] # Sort error list by its indices (lowest to highest)
    population = population[i,:] # Sort population by indices from error list (lowest to highest error)

    return population

def evolve(population, parambounds, retain = 0.25, random_select = 0.02, mutate_chance = 0.08):

    '''
    Evolves a sorted population of individuals. Random_select, mutate_chance
    and jump are compared to a random number generated between 0 and 1.
    Kind of old and not used in most things this module is used for.
    Should be modified and documented later.
    '''

    retain_length = int(len(population)*retain)

    parents = population[:retain_length]

    for indiv in population[retain_length:]:
        if random_select > random.random():
            parents = np.vstack((parents, indiv))

    replenish = len(population)-len(parents)

    children = np.empty([replenish, len(parambounds)])

    z = 0

    while z < replenish:

        malesele   = random.randint(0, len(parents)-1)
        femalesele = random.randint(0,len(parents)-1)
        child = np.empty([len(parambounds)])

        if malesele != femalesele:
            male   = parents[malesele]
            female = parents[femalesele]
            child[:int(len(parambounds)/2)]  = male[:int(len(parambounds)/2)]
            child[int(len(parambounds)/2):] = female[int(len(parambounds)/2):]
            children[z] = child
            z += 1

    for indiv in children: # Add mutations in children
        if mutate_chance > random.random():
            mutate_position = random.randint(0,len(parambounds)-1)
            indiv[mutate_position] = random.uniform(parambounds[mutate_position][0], parambounds[mutate_position][1])

    newpop = np.vstack((parents, children))

    return newpop

def evolution(data, model, bounds, size = 10, generations = 10):
    '''
    Parameter ranges from model.
    Kind of old and not used in most things this module is used for.
    Should be modified and documented later.
    '''

    pop = create_population(size,bounds)

    g = 0
    while g < generations:
        '''Screen population'''
        err = test_population(pop,data,model)

        '''Rank population'''
        sorted_pop = sort_population(pop, err)

        '''Evolve'''
        pop = evolve(pop, bounds, retain = 0.25, random_select = 0.02, mutate_chance = 0.08)
        g += 1

    '''Take the best parameter set after so many generations'''
    best_fit = sorted_pop[0]

    return  best_fit
