#!/usr/bin/env python3
# @matushalak
"""
DESCRIPTION:
    Template code for the BWT assignment in the Algorithms in Sequence Analysis course at the VU.
    
INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via corresponding Canvas assignment.

AUTHOR:
    <Matus Halak 2724858>
"""

import argparse

# Implement the following functions.

def read_fasta(filename):
    """Read a single sequence from the given FASTA file.
    The header of the sequence can be ignored, and if the file contains
    multiple sequences only the first one should be returned.
    If the file is not in correct FASTA format, an error message should
    be shown (e.g. by raising ValueError).
    """
    seq = []
    with open(filename,'r') as fa:
        for il, line in enumerate(fa):
            if il == 0 and line[0] != '>':
                raise ValueError('Error: Not Fasta format!')

            if line.startswith('>'):
                if not seq:
                    continue # skip first line with metadata
                else:
                    break # Stop if start of 2nd seq encountered
            else:
                seq.append(line.strip()) # add current line to sequence
    return ''.join(seq) # return as a string

def string_rotations(seq):
    """Return a list containing all rotations of the given sequence.
    The given sequence ends with a unique $ (appended in main).

    Example:
    >>> string_rotations('banana$')
    ['banana$', 'anana$b', 'nana$ba', 'ana$ban',
     'na$bana', 'a$banan', '$banana']
    """
    # return letters from index until end + letters until index
    return [seq[i:]+seq[:i] for i in range(len(seq))]

def bwt(rotations):
    """Return the Burrows-Wheeler Transform (BWT) of a sequence, given a list
    of its rotations.

    Example:
    >>> bwt(string_rotations('banana$'))
    'annb$aa'
    """
    # last letters of sorted btw list
    return ''.join([r[-1] for r in sorted(rotations)])

# TODO
def suffix_array_linear(string):
    # look up paper for linear suffix array algorithm & implement it
    pass

# TODO
def btw_linear(string):
    # use O(n) algorithm to construct suffix array
    pass

def rle(seq):
    """Return the Run-Length Encoding (RLE) of a string, as a string containing
    single characters alternating with digits.

    Example:
    >>> rle('annb$aa')
    'a1n2b1$1a2'
    """
    code = [seq[0]]
    count = 1
    for i, l in enumerate(seq):
        if seq[i-1] == l:
            count += 1
        else:
            count =1


    return ''.join(code)
    
def rle_invert(rle_seq):
    """Given a Run-Length Encoded string, return the original string in its
    uncompressed form.

    Example:
    >>> rle_invert('a1n2b1$1a2')
    'annb$aa'
    """
    raise NotImplementedError

def compute_rank_vector(bwt_seq):
    """Return the rank vector for the given BW-transformed string. The rank
    vector contains at each position i, the number of occurrences of bwt_seq[i]
    in positions before i; that is, the first A has rank 0, the third G has rank
    2, and so on.

    Example:
    >>> compute_rank_vector('annb$aa')
    [0, 0, 1, 0, 0, 1, 2]
    """
    raise NotImplementedError

# TODO
def better_rank_vector(bwt_seq):
    # compute the hash map later used for BWAlignment as well
    pass

def compute_f_map(bwt_seq):
    """Return, for the given BW-transformed string, a dictionary mapping each
    distinct character to its first occurrence in the first column of the
    Burrows-Wheeler matrix.

    Example:
    >>> compute_f_map('annb$aa')
    {'$': 0, 'a': 1, 'b': 4, 'n': 5}

    (The F-column for 'banana$' would be [$, a, a, a, b, n, n])
    """
    raise NotImplementedError

def bwt_invert(bwt_seq, rank, f_map):
    """Invert the Burrows-Wheeler Transform of a sequence, given the transformed
    sequence itself, the rank vector and the mapping of the F-column.

    Example:
    >>> seq = 'annb$aa'
    >>> bwt_invert(seq, compute_rank_vector(seq), compute_f_map(seq))
    'banana$'
    """
    raise NotImplementedError
        
# TODO
def bwa (subseq, seq):
    # align subseq to seq
    pass





# Code for testing:

def main():
    parser = argparse.ArgumentParser()
    file_or_stdin = parser.add_mutually_exclusive_group(required=True)
    file_or_stdin.add_argument('file', nargs='?', type=str, help='FASTA file to load')
    file_or_stdin.add_argument('--stdin', action='store_true', help='Read string as a single line from STDIN')
    file_or_stdin.add_argument(
            '--test',
            choices=['read_fasta', 'string_rotations', 'bwt',
                     'rle', 'rle_invert',
                     'compute_rank_vector', 'compute_f_map', 'bwt_invert',
                     'full'],
            help='Single function to test (codegrade)')
    args = parser.parse_args()

    if not args.test:
        rle_inv_seq = inv_seq = None
        try:
            # Calling the functions
            if args.stdin:
                print('Reading string from STDIN')
                seq = input() + '$'
            else:
                print('Reading ', args.file)
                seq = read_fasta(args.file) + '$'
            print('Sequence:')
            print(seq)

            rotations = string_rotations(seq)
            print('Computed rotations')

            bwt_seq = bwt(rotations)
            print('BWT of sequence:')
            print(bwt_seq)

            rle_seq = rle(bwt_seq)
            print('RLE of BWT:')
            print(rle_seq)

            rle_inv_seq = rle_invert(rle_seq)
            print('RLE inverted:')
            print(rle_inv_seq)

            rank = compute_rank_vector(rle_inv_seq)
            print('Computed rank vector')

            f_map = compute_f_map(rle_inv_seq)
            print('Computed F column map')

            inv_seq = bwt_invert(rle_inv_seq, rank, f_map)
            print('BWT inverted:')
            print(inv_seq)
        except NotImplementedError:
            print('(Exercise is unfinished.)')
        finally:
            # Check the resulting strings, even if not all tasks were finished
            if rle_inv_seq is not None:
                if rle_inv_seq == bwt_seq:
                    print('RLE inversion matches original BWT.')
                else:
                    print('RLE inversion does NOT match original BWT.')

            if inv_seq is not None:
                if inv_seq == seq:
                    print('BWT inversion matches original sequence.')
                else:
                    print('BWT inversion does NOT match original sequence.')
    else:
        # DO NOT CHANGE CODE BELOW -- NECESSARY FOR CODEGRADE

        inp = input()
        if args.test == 'read_fasta':
            print(read_fasta(inp))
        elif args.test == 'string_rotations':
            for rot in string_rotations(inp):
                print(rot)
        elif args.test == 'bwt':
            print(bwt(string_rotations(inp)))
        elif args.test == 'rle':
            print(rle(inp))
        elif args.test == 'rle_invert':
            print(rle_invert(inp))
        elif args.test == 'compute_rank_vector':
            for r in compute_rank_vector(inp):
                print(r)
        elif args.test == 'compute_f_map':
            for c, v in compute_f_map(inp).items():
                print(c, v)
        elif args.test == 'bwt_invert':
            ranks = compute_rank_vector(inp)
            f_map = compute_f_map(inp)
            print(bwt_invert(inp, ranks, f_map))
        elif args.test == 'full':
            inp = rle_invert(rle(bwt(string_rotations(inp))))
            print(bwt_invert(inp, compute_rank_vector(inp), compute_f_map(inp)))

if __name__ == '__main__':
    main()