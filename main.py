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

def main(args):
    if args.generate_data == 'True':
        generate_data(args)


def parse_args():
    parser = argparse.ArgumentParser("Solutions to Homework 1")
    parser.add_argument('-t', type=int, default=10, help='total number of dna strings')
    parser.add_argument('-l', type=int, default=10, help='implanted motif length')
    parser.add_argument('-d', type=int, default=4, help='total mutations in motif')
    parser.add_argument('-n', type=int, default=500, help='total length of a dna string')
    parser.add_argument('-generate_data', type=str, default='False', help='generate random dna strings')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()    
    main(args)