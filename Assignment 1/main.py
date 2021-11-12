import argparse
import numpy as np
import time
import copy


def generate_data(args):
    implanted_motif = np.random.choice(a=['a', 'c', 'g', 't'], size=args.l)
    strings = []
    fw = open('dna_strings.txt', 'w')

    for iter in range(args.t):
        dna_str = np.random.choice(a=['a', 'c', 'g', 't'], size=args.n)

        implanted_motif_copy = np.copy(implanted_motif)
        total_mutations = 0

        while True:

            if total_mutations == args.d:
                break

            random_index = np.random.choice(a=args.l, size=1)
            random_bp = np.random.choice(a=['a', 'c', 'g', 't'], size=1)

            if implanted_motif[random_index] == random_bp:
                continue

            implanted_motif_copy[random_index] = random_bp
            total_mutations += 1

        random_dna_index = np.random.choice(a=np.arange(args.n - args.l + 1), size=1)[0]
        dna_str[random_dna_index:random_dna_index + args.l] = implanted_motif_copy
        strings.append(dna_str)
        fw.write(''.join(implanted_motif_copy) + ' {}'.format(random_dna_index) + '\n')

    fw.write('original motif: {}\n'.format(''.join(implanted_motif)))

    for dna_str in strings:
        fw.write(''.join(dna_str) + '\n')

    fw.close()


def read_data(args):
    lines = []
    with open(args.input_file) as fr:
        lines = fr.readlines()

    dna_strings = []
    for dna_str in lines[11:]:
        dna_strings.append(dna_str.strip())

    return dna_strings


def get_profile(best_motifs, randomized=True):
    profile = {}
    for motif in best_motifs:
        for i in np.arange(len(motif)):
            profile.setdefault(motif[i], [0] * len(motif))
            profile[motif[i]][i] += 1

    for bp in profile:
        if randomized:
            profile[bp] = [count / len(best_motifs) for count in profile[bp]]
        else:  # gibbs
            profile[bp] = [(count + 1) / (len(best_motifs) + 4) for count in profile[bp]]

            pass

    return profile


def get_prob_score(profile, kmer):
    score = 1
    for i in np.arange(len(kmer)):
        score *= profile[kmer[i]][i]

    return score


def calculate_prob_of_deleted_string(dna, profile, args):
    prob_matrix = {}    # keys are the k-mers in the dna string, values equal the probability of the k-mer

    '''
    takes the k-mer and calculate its probability 
    based on the profile and add to the prob_matrix dictionary
    '''
    for i in np.arange(len(dna) - args.k + 1):
        kmer = dna[i:i + args.k]
        prob_score = get_prob_score(profile, kmer)
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
    prob_matrix = {}
    most_likely_motifs = []

    for i in np.arange(len(dna)):
        for j in np.arange(len(dna[i]) - args.k + 1):
            prob_matrix.setdefault(f'dna{i + 1}', ['', -1])
            kmer = dna[i][j:j + args.k]
            prob_score = get_prob_score(profile, kmer)
            if prob_score > prob_matrix[f'dna{i + 1}'][1]:
                prob_matrix[f'dna{i + 1}'][0] = kmer
                prob_matrix[f'dna{i + 1}'][1] = prob_score

        most_likely_motifs.append(prob_matrix[f'dna{i + 1}'][0])

    return most_likely_motifs


def get_score(motifs):
    score = 0
    profile = get_profile(motifs)

    for i in np.arange(len(profile['a'])):
        dominant_bp = max(profile['a'][i], profile['g'][i], profile['c'][i], profile['t'][i])
        score += (1 - dominant_bp)

    return int(score * len(motifs))


def get_initial_motifs(dna, args):
    motifs = []
    random_index = np.random.choice(args.n - args.k + 1, 1)[0]
    for dna_str in dna:
        motifs.append(dna_str[random_index:random_index + args.k])
        random_index = np.random.choice(args.n - args.k + 1, 1)[0]
    return motifs


def gibbs_sampler(dna, args):
    motifs = get_initial_motifs(dna, args)  # get 10 motifs randomly
    best_motifs = copy.copy(motifs)         # assign the best motifs to the starting motifs.
    iter = 0

    # If the best score is not renewed 50 times in a row, the loop will end
    while iter < 50:

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
    motifs = get_initial_motifs(dna, args)
    best_motifs = copy.copy(motifs)

    while True:
        profile = get_profile(motifs)
        motifs = get_likely_motifs(profile, dna, args)

        if get_score(motifs) < get_score(best_motifs):
            best_motifs = copy.copy(motifs)
        else:
            return best_motifs, get_score(best_motifs)


def main(args):
    if args.generate_data == 'True':
        generate_data(args)

    dna_strings = read_data(args)

    scores = []
    best_motifs = {}
    running_time = []


    for iter in np.arange(args.total_iter):
        start = time.time()
        if args.algorithm == 'randomized':
            motifs, score = randomized_motif_search(dna_strings, args)
            scores.append(score)
            best_motifs.setdefault('best motif', [0, 100])
            if score < best_motifs['best motif'][1]:
                best_motifs['best motif'][0] = motifs
                best_motifs['best motif'][1] = score

        elif args.algorithm == 'gibbs':
            motifs, score = gibbs_sampler(dna_strings, args)
            scores.append(score)
            best_motifs.setdefault('best motif', [0, 100])
            if score < best_motifs['best motif'][1]:
                best_motifs['best motif'][0] = motifs
                best_motifs['best motif'][1] = score

        end = time.time()
        running_time.append(end - start)

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
    parser.add_argument('-total_iter', type=int, default=10, help='total runs for the greedy algorithm')
    parser.add_argument('-input_file', type=str, default='dna_strings.txt', help='input file')
    parser.add_argument('-generate_data', type=str, default='False', help='generate random dna strings')
    parser.add_argument('-algorithm', type=str, default='gibbs',
                        help='motif searching algorithm i.e. gibbs, randomized')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)
