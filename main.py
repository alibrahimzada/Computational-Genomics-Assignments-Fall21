import argparse
import numpy as np

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
        
        random_dna_index = np.random.choice(a=np.arange(args.n-args.l+1), size=1)[0]
        dna_str[random_dna_index:random_dna_index+args.l] = implanted_motif_copy
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

def get_profile(best_motifs):
    profile = {}
    for motif in best_motifs:
        for i in np.arange(len(motif)):
            profile.setdefault(motif[i], [0] * len(motif))
            profile[motif[i]][i] += 1

    for bp in profile:
        profile[bp] = [count / len(best_motifs) for count in profile[bp]]
    
    return profile

def get_prob_score(profile, kmer):
    score = 1
    for i in np.arange(len(kmer)):
        score *= profile[kmer[i]][i]

    return score

def get_likely_motifs(profile, dna, args):
    prob_matrix = {}
    most_likely_motifs = []

    for i in np.arange(len(dna)):
        for j in np.arange(len(dna[i])-args.k+1):
            prob_matrix.setdefault(f'dna{i+1}', ['', -1])
            kmer = dna[i][j:j+args.k]
            prob_score = get_prob_score(profile, kmer)
            if prob_score > prob_matrix[f'dna{i+1}'][1]:
                prob_matrix[f'dna{i+1}'][0] = kmer
                prob_matrix[f'dna{i+1}'][1] = prob_score
        
        most_likely_motifs.append(prob_matrix[f'dna{i+1}'][0])

    return most_likely_motifs

def get_score(motifs):
    score = 0
    profile = get_profile(motifs)

    for i in np.arange(len(profile['a'])):
        dominant_bp = max(profile['a'][i], profile['g'][i], profile['c'][i], profile['t'][i])
        score += (1 - dominant_bp)

    return int(score * len(motifs))

def randomized_motif_search(dna, args):
    random_index = np.random.choice(args.n-args.k+1, 1)[0]
    best_motifs = [dna_str[random_index:random_index+args.k] for dna_str in dna]

    while True:
        profile = get_profile(best_motifs)
        motifs = get_likely_motifs(profile, dna, args)

        if get_score(motifs) > get_score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs

def main(args):
    if args.generate_data == 'True':
        generate_data(args)

    dna_strings = read_data(args)

    for iter in np.arange(args.total_iter):
        if args.algorithm == 'randomized':
            best_motifs = randomized_motif_search(dna_strings, args)
            for m in best_motifs:
                print(m)

        break

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
    parser.add_argument('-algorithm', type=str, default='randomized', help='motif searching algorithm i.e. gibbs, randomized')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()    
    main(args)