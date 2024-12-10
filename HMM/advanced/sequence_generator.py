#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the FIRST Advanced Question of the Hidden Markov Models
    assignment in the Algorithms in Sequence Analysis course at the VU.

INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via
    corresponding Canvas assignment. Note this script will be graded manually,
    if and only if your "hmm.py" script succesfully implements Baum-Welch
    training! Continuous Feedback will not be available for this script.

AUTHOR:
    <Matus Halak 2724858>
"""

from argparse import ArgumentParser, RawTextHelpFormatter
from hmm_utility import load_tsv
from numpy.random import choice



def parse_args():
    #####################
    # START CODING HERE #
    #####################
    # Implement a simple argument parser (WITH help documentation!) that parses
    # the information needed by main() from commandline. Take a look at the
    # argparse documentation, the parser in hmm_utility.py or align.py
    # (from the Dynamic Programming exercise) for hints on how to do this.

    parser = ArgumentParser(prog = 'python3 sequence_generator.py',
        formatter_class = RawTextHelpFormatter, description =
        '  Perform the specified algorithm, with given sequences and parameters.\n\n'
        '  Example syntax:\n'
        '    python3 sequence_generator.py A.tsv E.tsv'
        '    python3 hmm.py baumwelch in.fa priorA priorE -o ./outputs -i 1')

    # Positionals
    parser.add_argument('transition', help='path to a TSV formatted transition matrix')
    parser.add_argument('emission', help='path to a TSV formatted emission matrix')

    # Optionals
    parser.add_argument('-n', dest='n_seqs', type=int, default=1,
        help='How many sequences should be generated into fasta file, by default 1 sequence is generated!\n ')
    
    parser.add_argument('-o', dest='out_dir', default = 'sequences',
        help='path to a directory where output files are saved\n'
             '  (directory will be made if it does not exist)\n'
             '  (file names and contents depend on algorithm)')

    return parser.parse_args()
    
    #####################
    #  END CODING HERE  #
    #####################


def generate_sequence(A,E):
    #####################
    # START CODING HERE #
    #####################
    # Implement a function that generates a random sequence using the choice()
    # function, given a Transition and Emission matrix.
    
    # Look up its documentation online:
    # https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.random.choice.html
    # Start with start state
    state = 'B'
    sequence = ''        
    while True:
        to_states = list(A[state])
        transition_probs = list(A[state].values())   
        which_state = choice(a=to_states, p=transition_probs)
        # upon reaching end state break
        if which_state == 'E':
            break

        symbols = list(E[which_state])
        emission_probs = list(E[which_state].values())
        generated_symbol = choice(a=symbols, p=emission_probs)
        
        sequence += str(generated_symbol)
        state = which_state 
        # print(sequence)

    #####################
    #  END CODING HERE  #
    #####################
    
    return sequence



def main():
    args = parse_args()
    #####################
    # START CODING HERE #
    #####################
    # Uncomment and complete (i.e. replace '?' in) the lines below:
    N = args.n_seqs               # The number of sequences to generate
    out_file = args.out_dir        # The file path to which to save the sequences
    A = load_tsv(args.transition) # Transition matrix
    E = load_tsv(args.emission) # Emission matrix

    for seq_id in range(N):
        seq = generate_sequence(A, E)
        print(f'Sequence {seq_id}: ', seq)
        if out_file:
            fasta = out_file+'.fasta'
            with open(fasta, 'a') as f:
                f.write('>random_sequence_%i\n%s\n' % (seq_id,seq))
            
    if N == 0 and out_file:
        fasta = out_file+'.fasta'
        with open(fasta, 'a') as f:
            f.write('')

    #####################
    #  END CODING HERE  #
    #####################
    
if __name__ == "__main__":
    main()
