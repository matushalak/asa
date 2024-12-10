#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the SECOND Advanced Question of the Hidden Markov Models
    assignment in the Algorithms in Sequence Analysis course at the VU.

INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via
    corresponding Canvas assignment. Note this script is graded automatically,
    if and only if your "hmm.py" script succesfully implements Baum-Welch
    training!

AUTHOR:
    <Matus Halak 2724858>
"""

import os.path as op

from os import makedirs
from argparse import ArgumentParser, RawTextHelpFormatter
from hmm_utility import load_fasta, load_tsv, serialize
from hmm import viterbi



def parse_args():
    "Parses inputs from commandline and returns them as a Namespace object."

    parser = ArgumentParser(prog = 'python3 viterbi_training.py',
        formatter_class = RawTextHelpFormatter, description =
        '  Perform Viterbi training, given a set of sequences with A and E priors.\n\n'
        '  Example syntax:\n'
        '    python3 hmm.py seq.fasta A.tsv E.tsv -i 100 -o /viterbi_outputs'
        '    python3 hmm.py baumwelch in.fa priorA priorE -o ./outputs -i 1')

    # Positionals
    parser.add_argument('fasta', help='path to a FASTA formatted input file')
    parser.add_argument('transition', help='path to a TSV formatted transition matrix')
    parser.add_argument('emission', help='path to a TSV formatted emission matrix')

    # Optionals
    parser.add_argument('-o', dest='out_dir',
        help='path to a directory where output files are saved\n'
             '  (directory will be made if it does not exist)')
    parser.add_argument('-i', dest='max_iter', type=int, default=20,
        help='maximum number of iterations (default: 20 )')

    return parser.parse_args()

'''
Idea & Overview:
For every sequence, Viterbi algorithm allows us to obtain most likely path GIVEN A & E priors

Starting with SOME priors A & E, we run many iterations, where we try to RE-estimate A & E by
counting Transitions & Emissions in each most Vitebri path for each sequence

When we normalize these counts, we get a new estimate of A & E which we use as priors for next iteration
'''

def train_viterbi(X,A,E):
    #####################
    # START CODING HERE #
    #####################
    # Just like Baum Welch
    # Initialize a new (posterior) Transition and Emission matrix
    new_A = {}
    # rows -states
    for k in A:
        # columns - states
        new_A[k] = {l:0 for l in A[k]}
    
    new_E = {}
    # rows -states
    for k in E:
        # columns - emitted symbols
        new_E[k] = {s:0 for s in E[k]}

    # Get the state path of every sequence in X,
    # using the viterbi() function imported from hmm.py
    for seq in X:
        # state_path, Vitebri Prob, Vitebri Trellis for given training sequence
        state_path, _, _ = viterbi(seq, A, E)
        # Count the transitions and emissions for every state in the optimal state path
        for i, state in enumerate(state_path):
            if i == 0:
                last_state = 'B'
            else:
                last_state = state_path[i-1]
            # count transitions
            new_A[last_state][state] += 1
            # count emissions
            new_E[state][seq[i]] += 1
        # must include Transition to END state !!!
        last_state = state # use last state
        state = 'E'
        new_A[last_state][state] += 1

    # Normalize your row sums
    for fromState in new_A:
        row_sum = sum(new_A[fromState].values())
        if fromState != 'E':
            for toState in new_A[fromState]:            
                new_A[fromState][toState] /= row_sum

    for fromState in new_E:
        row_sum = sum(new_E[fromState].values())
        for emittedSymbol in new_E[fromState]:
            new_E[fromState][emittedSymbol] /= row_sum

    #####################
    #  END CODING HERE  #
    #####################
    # breakpoint()
    return new_A, new_E


def main(args = False):
    "Perform Viterbi training, given a set of sequences with A and E priors."
    
    # Process arguments and load specified files
    if not args: args = parse_args()

    set_X, labels = load_fasta(args.fasta) # List of sequences, list of labels
    A = load_tsv(args.transition) # Nested Q -> Q dictionary
    E = load_tsv(args.emission)   # Nested Q -> S dictionary
    
    i_max = args.max_iter
    
    #####################
    # START CODING HERE #
    #####################
    # Iterate until you've reached i_max or until your parameters have converged!
    # Note Viterbi converges discretely (unlike Baum-Welch), so you don't need to
    # track your Sum Log-Likelihood to decide this.
    for i in range(i_max):
        A_new, E_new = train_viterbi(set_X, A, E)
        
        if A_new == A and E_new == E:
            print(f'You have converged in iteration {i+1}!')
            break

        A, E = A_new, E_new
    #####################
    #  END CODING HERE  #
    #####################

    
    if args.out_dir:
        makedirs(args.out_dir, exist_ok=True) # Make sure the output directory exists.
        A_path = op.join(args.out_dir,'viterbi_posterior_A')
        with open(A_path,'w') as f: f.write(serialize(A))
        E_path = op.join(args.out_dir,'viterbi_posterior_E')
        with open(E_path,'w') as f: f.write(serialize(E))        



if __name__ == "__main__":
    main()