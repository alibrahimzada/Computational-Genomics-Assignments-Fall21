import argparse
import numpy as np
import time
import copy


def generate_data(args):
    '''
    creates the dna_strings.txt
    we dont use this part in our code, but we kept it for the case we need to produced a input file
    '''
    implanted_motif = np.random.choice(a=['a', 'c', 'g', 't'], size=args.l) # create original motif
    strings = []
    fw = open(args.input_file, 'w') # open the output file

    for iter in range(args.t):
        dna_str = np.random.choice(a=['a', 'c', 'g', 't'], size=args.n) # generate a random dna string of size n
        implanted_motif_copy = np.copy(implanted_motif)
        total_mutations = 0

        while True:

            if total_mutations == args.d: # if total mutations have reached d, then break
                break

            random_index = np.random.choice(a=args.l, size=1) # choose a random index
            random_bp = np.random.choice(a=['a', 'c', 'g', 't'], size=1) # choose a random base pair

            if implanted_motif[random_index] == random_bp: # if random base pair is same as the randomly chosen bp in dna, then continue
                continue

            implanted_motif_copy[random_index] = random_bp # create mutation
            total_mutations += 1

        random_dna_index = np.random.choice(a=np.arange(args.n - args.l + 1), size=1)[0] # choose a random index in the dna string
        dna_str[random_dna_index:random_dna_index + args.l] = implanted_motif_copy # implant the mutated motif in the dna string
        strings.append(dna_str)
        fw.write(''.join(implanted_motif_copy) + ' {}'.format(random_dna_index) + '\n')

    fw.write('original motif: {}\n'.format(''.join(implanted_motif)))

    for dna_str in strings:
        fw.write(''.join(dna_str) + '\n')

    fw.close()


def read_data(args):
    '''
    reads the input file and returns dna strings as an array
    '''
    lines = []
    with open(args.input_file) as fr: # read all lines of the input file
        lines = fr.readlines()

    print(lines[10].strip()) # print original motif
    dna_strings = []
    for dna_str in lines[11:]: # extract dna strings
        dna_strings.append(dna_str.strip())

    return dna_strings


def get_profile(best_motifs, randomized=True):
    '''
    creates a profile derived from best_motifs array.
    Also, we seperated gibbs and randomized search algorithms
    so, we could prevent any 0 chance probabilities for gibbs algorithm.
    '''
    profile = {} # a profile dictionary
    for motif in best_motifs:
        for i in np.arange(len(motif)):
            profile.setdefault(motif[i], [0] * len(motif)) # create the required space for a bp
            profile[motif[i]][i] += 1 # increment specific bp location

    for bp in profile:
        if randomized: # if randomized algorithm, then normalize
            profile[bp] = [count / len(best_motifs) for count in profile[bp]]
        else:  # if gibbs algorithm, then add 1 and then normalize
            profile[bp] = [(count + 1) / (len(best_motifs) + 4) for count in profile[bp]]

    return profile


def get_prob_score(profile, kmer):
    '''
    calculates the score of given kmer according to profile
    '''
    score = 1
    for i in np.arange(len(kmer)):
        score *= profile[kmer[i]][i] # calculate the prob of a kmer based on the profile

    return score


def calculate_prob_of_deleted_string(dna, profile, args):
    '''
    takes the k-mer and calculate its probability
    based on the profile and add to the prob_matrix dictionary
    '''
    prob_matrix = {}    # keys are the k-mers in the dna string, values equal the probability of the k-mer

    for i in np.arange(len(dna) - args.k + 1):
        kmer = dna[i:i + args.k] # get the current window kmer
        prob_score = get_prob_score(profile, kmer) # calculate its probability score
        prob_matrix[kmer] = prob_score

    return prob_matrix


def roll_dice(probabilities):
    '''
    this function takes probabilities dictionary as its argument and chooses the new k-mer randomly.
    but the k-mer with the highest probability is more likely to be chosen
    '''
    random_value = np.random.uniform(0, sum(probabilities.values()))

    total = 0
    for motif in probabilities:
        total += probabilities[motif]
        if total >= random_value:
            return motif


def get_likely_motifs(profile, dna, args):
    '''
    chooses motives that has highest profile scores from each dna string
    and return them as an array.
    '''
    prob_matrix = {}
    most_likely_motifs = []

    for i in np.arange(len(dna)):
        for j in np.arange(len(dna[i]) - args.k + 1):
            prob_matrix.setdefault(f'dna{i + 1}', ['', -1])
            kmer = dna[i][j:j + args.k]
            prob_score = get_prob_score(profile, kmer)
            if prob_score > prob_matrix[f'dna{i + 1}'][1]: # check if the prob. score is better than the current best kmer
                prob_matrix[f'dna{i + 1}'][0] = kmer
                prob_matrix[f'dna{i + 1}'][1] = prob_score

        most_likely_motifs.append(prob_matrix[f'dna{i + 1}'][0])

    return most_likely_motifs


def get_score(motifs):
    '''
    calculates the differences between motif set and
    returns that as score of motif set.
    '''
    score = 0
    profile = get_profile(motifs)

    for i in np.arange(len(profile['a'])): # dominant_bp = dominant base pair, column-wise
        dominant_bp = max(profile['a'][i], profile['g'][i], profile['c'][i], profile['t'][i])
        score += (1 - dominant_bp)

    return int(score * len(motifs))


def get_initial_motifs(dna, args):
    '''
    randomly chooses motives from given dna array.
    '''
    motifs = []
    random_index = np.random.choice(args.n - args.k + 1, 1)[0]
    for dna_str in dna:
        motifs.append(dna_str[random_index:random_index + args.k]) # store the randomly selected motif
        random_index = np.random.choice(args.n - args.k + 1, 1)[0] # choose a new random index for the next dna string
    return motifs


def gibbs_sampler(dna, args):
    '''
    gibbs algorithm
    '''
    motifs = get_initial_motifs(dna, args)  # get 10 motifs randomly
    best_motifs = copy.copy(motifs)         # assign the best motifs to the starting motifs.
    iter = 0

    # If the best score is not renewed 50 times in a row, the loop will end
    while iter < args.gibbs_iters:

        temp_motif = copy.copy(motifs)  # temporary list is used because 1 motif will be eliminated
        # the index of the to-be-eliminated motif is picked at random
        position_of_random_motif = np.random.randint(0, len(motifs))
        temp_motif.pop(position_of_random_motif)  # removes 1 randomly motif
        profile = get_profile(temp_motif, False)    # get profile of the motifs
        # probs of the k-mers in deleted dna string based on the profile are computed and stored in probabilities.
        probabilities = calculate_prob_of_deleted_string(dna[position_of_random_motif], profile, args)
        # randomly_chosen_motif equals to the motif is determined by throwing dice from the deleted dna string
        randomly_chosen_motif = roll_dice(probabilities)

        # new selected motif is added to the motifs
        motifs.pop(position_of_random_motif)
        motifs.insert(position_of_random_motif, randomly_chosen_motif)

        '''
        if the previous score is better than the best score, 
        the best score is synchronized to the previous score and the iter is reset,
        otherwise, iter is incremented by 1
        '''
        if get_score(motifs) < get_score(best_motifs):
            best_motifs = copy.copy(motifs)
            iter=0

        else:
            iter += 1

    return best_motifs, get_score(best_motifs)


def randomized_motif_search(dna, args):
    '''
    randomized motif search algorithm
    '''
    motifs = get_initial_motifs(dna, args) # get an initial set of motifs in the beginning
    best_motifs = copy.copy(motifs) # set the initial motifs as the best motifs

    while True:
        profile = get_profile(motifs) # get the motifs profile
        motifs = get_likely_motifs(profile, dna, args) # get the likely motifs based on profile

        if get_score(motifs) < get_score(best_motifs): # check if the score has been improved
            best_motifs = copy.copy(motifs)
        else:
            return best_motifs, get_score(best_motifs)


def main(args):
    if args.generate_data == 'True': # if generate data is True, then generate the random dataset
        generate_data(args)

    dna_strings = read_data(args) # read the data and store it in a data structure

    scores = []
    best_motifs = {}
    running_time = []

    '''
    runs algorithms according the given argument and compares their best scores of chosen algorithm and prints it
    '''
    for iter in np.arange(args.total_iter): # run the algorithms multiple times. i.e. 10
        start = time.time()
        motifs, score = 0, 0

        if args.algorithm == 'randomized': # check if algorithm is randomized motif search
            motifs, score = randomized_motif_search(dna_strings, args)

        elif args.algorithm == 'gibbs': # check if algorithm is gibbs sampling
            motifs, score = gibbs_sampler(dna_strings, args)

        scores.append(score) # store the best score of the run
        best_motifs.setdefault('best motif', [0, 100])
        if score < best_motifs['best motif'][1]: # check if the recent score beats the all time best score
            best_motifs['best motif'][0] = motifs # update motifs
            best_motifs['best motif'][1] = score # update score

        end = time.time()
        running_time.append(end - start) # store the elapsed time

    print('max score: {}\naverage score: {}\nbest score: {}'.format(max(scores), np.mean(scores), min(scores)))
    print('score: {}\nbest motifs: '.format(best_motifs['best motif'][1]))
    for motif in best_motifs['best motif'][0]:
        print(motif)
    print('average running time per iteration (s): {:.2f}\ntotal running time (s): {:.2f}'.format(np.mean(running_time),
                                                                                                  np.sum(running_time)))


def parse_args():
    parser = argparse.ArgumentParser("Solutions to Homework 1")
    parser.add_argument('-t', type=int, default=10, help='total number of dna strings')
    parser.add_argument('-l', type=int, default=10, help='implanted motif length')
    parser.add_argument('-d', type=int, default=4, help='total mutations in motif')
    parser.add_argument('-n', type=int, default=500, help='total length of a dna string')
    parser.add_argument('-k', type=int, default=10, help='length of consensus string')
    parser.add_argument('-gibbs_iters', type=int, default=50, help='total iterations of a gibbs run')    
    parser.add_argument('-total_iter', type=int, default=10, help='total runs for the greedy algorithm')
    parser.add_argument('-input_file', type=str, default='dna_strings.txt', help='input file')
    parser.add_argument('-generate_data', type=str, default='False', help='generate random dna strings')
    parser.add_argument('-algorithm', type=str, default='gibbs',
                        help='motif searching algorithm i.e. gibbs, randomized')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)
