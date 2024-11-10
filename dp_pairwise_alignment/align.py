#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the Dynamic Programming assignment in the Algorithms in Sequence Analysis course at the VU.
    
INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via corresponding Canvas assignment.

AUTHOR:
    Matus Halak - 2724858
"""



import argparse
import pickle



def parse_args():
    "Parses inputs from commandline and returns them as a Namespace object."

    parser = argparse.ArgumentParser(prog = 'python3 align.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Aligns the first two sequences in a specified FASTA\n'
        '  file with a chosen strategy and parameters.\n'
        '\n'
        'defaults:\n'
        '  strategy = global\n'
        '  substitution matrix = pam250\n'
        '  gap penalty = 2')
        
    parser.add_argument('fasta', help='path to a FASTA formatted input file')
    parser.add_argument('output', nargs='*', 
        help='path to an output file where the alignment is saved\n'
             '  (if a second output file is given,\n'
             '   save the score matrix in there)')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
        help='print the score matrix and alignment on screen', default=False)
    parser.add_argument('-s', '--strategy', dest='strategy',
        choices=['global','semiglobal','local'], default="global")
    parser.add_argument('-m', '--matrix', dest='substitution_matrix',
        choices=['pam250','blosum62','identity'], default='pam250')
    parser.add_argument('-g', '--gap_penalty', dest='gap_penalty', type=int,
        help='must be a positive integer', default=2)

    args = parser.parse_args()

    args.align_out = args.output[0] if args.output else False
    args.matrix_out = args.output[1] if len(args.output) >= 2 else False
                      # Fancy inline if-else statements. Use cautiously!
                      
    if args.gap_penalty <= 0:
        parser.error('gap penalty must be a positive integer')

    return args



def load_substitution_matrix(name):
    "Loads and returns the specified substitution matrix from a pickle (.pkl) file."
    # Substitution matrices have been prepared as nested dictionaries:
    # the score of substituting A for Z can be found with subst['A']['Z']
    # NOTE: Only works if working directory contains the correct folder and file!
    
    with open('substitution_matrices/%s.pkl' % name, 'rb') as f:
        subst = pickle.load(f)
    return subst
    
    

def load_sequences(filepath):
    "Reads a FASTA file and returns the first two sequences it contains."
    
    seq1 = []
    seq2 = []
    with open(filepath,'r') as f:
        for line in f:
            if line.startswith('>'):
                if not seq1:
                    current_seq = seq1
                elif not seq2:
                    current_seq = seq2
                else:
                    break # Stop if a 3rd sequence is encountered
            else:
                current_seq.append(line.strip())
    
    if not seq2:
        raise Exception('Error: Not enough sequences in specified FASTA file.')
    
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)
    return seq1, seq2



def align(seq1, seq2, strategy, substitution_matrix, gap_penalty):
    "Do pairwise alignment using the specified strategy and parameters."
    # This function consists of 3 parts:
    #
    #   1) Initialize a score matrix as a "list of lists" of the appropriate length.
    #      Fill in the correct values for the first row and column given the strategy.
    #        (local / semiglobal = 0  --  global = stacking gap penalties)
    #   2) Fill in the rest of the score matrix using Dynamic Programming, accounting
    #      for the selected alignment strategy, substitution matrix and gap penalty.
    #   3) Perform the correct traceback routine on your filled in score matrix.
    #
    # Both the resulting alignment (sequences with gaps and the corresponding score)
    # and the filled in score matrix are returned as outputs.
    #
    # NOTE: You are strongly encouraged to think about how you can reuse (parts of)
    #       your code between steps 2 and 3 for the different strategies!
    
    
    ### 1: Initialize
    M = len(seq1)+1
    N = len(seq2)+1
    score_matrix, pointers_matrix = [], []
    # empty matrix full of zeros
    for i in range(M):
        row, prow = [], []
        score_matrix.append(row)
        pointers_matrix.append(prow)
        for j in range(N):
            row.append(0)
            prow.append([])
    
    # initialize 1st row & col with gap penalties
    if strategy == 'global':
        #####################
        # START CODING HERE #
        #####################
        # first row
        score_matrix[0] = [-i*gap_penalty for i in range(N)]
        pointers_matrix[0] = [(0,i-1) if i > 0 else (0,0) for i in range(N)]
        # first column
        for j in range(1, M):
            score_matrix[j][0] = -j * gap_penalty
            pointers_matrix[j][0] = (j-1, 0)
    
    # SEMIGLOBAL
        # initialize 1st row & col with 0s == don't penalize gaps at the beginning of either sequence
        # 0s first row - X penalize gaps in seq1 / x / query; 
        # 0s first col - X penalize gaps in seq2 / y / reference
        # if you only want to allow gaps at beginning of one sequence but NOT another - initialize seq. w gaps with 0s but sequence without gaps same as global alignment
    # LOCAL - also initialize with zeros
    elif strategy in ('semiglobal', 'local'):
        minus = 1 if strategy == 'semiglobal' else 0
        # first row
        score_matrix[0] = [0 for _ in range(N)]
        pointers_matrix[0] = [(0,i-minus) if i > 0 else (0,0) for i in range(N)]
        # first column
        for j in range(1, M):
            score_matrix[j][0] = 0
            pointers_matrix[j][0] = (j-minus, 0)

        #####################
        #  END CODING HERE  #
        #####################

    ### 2: Fill in Score Matrix
    # same between global & semiglobal & local
    #####################
    # START CODING HERE #
    ##################### 
    def dp_function(ri,cj, X = seq1, Y = seq2,
                    MAT = score_matrix, SCORING = substitution_matrix, GAP = gap_penalty):
        # Seq1 = X
        # Seq2 = Y        
        neighbors = [MAT[ri][cj-1] - GAP, # left
                     MAT[ri-1][cj-1] + SCORING[X[ri-1]][Y[cj-1]], # diagonal
                     MAT[ri-1][cj] - GAP] # up
        
        best_nb = max(neighbors)

        neighbor_map = {
            0: (ri, cj-1), # low road
            1: (ri - 1, cj - 1), # diagonal
            2: (ri-1, cj)} # high road

        # stores all pointers for a given field
        pointers = {i:neighbor_map[i] for i, n in enumerate(neighbors) if n == best_nb}
        # high road take highest index
        pointers = pointers[max(pointers)]
        # if local check if sub / gap results in negative score
        if strategy ==  'local':
            # make sure that stop only at HARD zeros
            pointers = (ri,cj) if max(best_nb, 0) == 0 and 0 not in neighbors else pointers
            return max(best_nb, 0), pointers
        return best_nb, pointers
    
    for i in range(1,M):
        for j in range(1,N):
            score_matrix[i][j], pointers_matrix[i][j] = dp_function(i,j)    
    #####################
    #  END CODING HERE  #
    #####################   
    
    
    ### 3: Traceback
    # High Road: favor (mis)matches over gaps; otherwise randomly choose insertion / deletion
    #####################
    # START CODING HERE #
    #####################   
    # initialize with empty string
    aligned_seq1 = ''
    aligned_seq2 = ''
    
    if strategy == 'global':
        align_score = score_matrix[-1][-1] # alignment score from score matrix
        last_pointer = (M-1, N-1) # bottom right
        current_pointer = pointers_matrix[-1][-1]
    
    # different traceback
    elif strategy == 'semiglobal':
        # where you start from determines in which sequence you allow free gaps at the end
        # makes sense to allow free gaps at the end of the shorter sequence (but apparently not for this assignment)
        maxes = []
        max_val = score_matrix[-1][-1]
        for ic, val in enumerate(score_matrix[-1]):
            if val > max_val:
                max_val = val
                maxes = [(M-1, ic)]
            elif val == max_val:
                maxes.append((M-1, ic))
        for ir, val in enumerate([r[-1] for r in score_matrix]):
            if val > max_val:
                max_val = val
                maxes = [(ir, N-1)]
            elif val == max_val:
                maxes.append((ir, N-1))
        
        # high road - get absolute max from last row & last col, prioritizing lowest row index in last col
        start_index = sorted(maxes, key = lambda x: (x[1], -x[0]), reverse = True)[0]
        align_score = score_matrix[start_index[0]][start_index[1]]
        
        D = (start_index[0] - (M-1),
             start_index[1] - (N-1))
        end_gaps = 0 if D[0] == 0 else 1 # which sequence has gaps at the end
        no_gaps = abs(end_gaps - 1)
        n_gaps = abs(D[no_gaps])
        
        seqs = (seq1, seq2)
        aseqs = [aligned_seq1, aligned_seq2]
        aseqs[no_gaps] += seqs[no_gaps][-n_gaps:][::-1] if n_gaps > 0 else ''
        aseqs[end_gaps] += ''.join(['-' for _ in range(n_gaps)])

        # rest of traceback same as in global alignment
        # strings are immutable so need to do this for it to be used later
        aligned_seq1, aligned_seq2 = aseqs 
        last_pointer = start_index
        current_pointer = pointers_matrix[start_index[0]][start_index[1]]
    
    # traceback from max
    elif strategy ==  'local':
        # get all the max scores (stupid in native python, easy and faster in numpy)
        maxes = []
        max_val = score_matrix[0][0]
        for ir, rr in enumerate(score_matrix):
            for ic, val in enumerate(rr):
                if val > max_val:
                    max_val = val
                    maxes = [(ir,ic)]
                elif val == max_val:
                    maxes.append((ir,ic))

        # high road choice of optimal local alignment        
        start_row, start_col = sorted(maxes, key = lambda x: (x[1], -x[0]), reverse = True)[0]

        align_score = score_matrix[start_row][start_col]
        last_pointer = (start_row, start_col)
        current_pointer = pointers_matrix[start_row][start_col]

    # traverse pointers matrix
    while True:
        diff = (current_pointer[0] - last_pointer[0],
                current_pointer[1] - last_pointer[1])
        
        if last_pointer == current_pointer:
            break

        # both -1: diagonal
        if set(diff) == {-1}:
            aligned_seq1 += seq1[current_pointer[0]]
            aligned_seq2 += seq2[current_pointer[1]]
        # - row == high road, gap in Y
        elif diff[0] == -1:
            aligned_seq1 += seq1[current_pointer[0]]
            aligned_seq2 += '-'
        # - col == low road, gap in X
        elif diff[1] == -1:
            aligned_seq1 += '-'
            aligned_seq2 += seq2[current_pointer[1]]
        
        last_pointer = current_pointer
        current_pointer = pointers_matrix[current_pointer[0]][current_pointer[1]]

    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]
    #####################
    #  END CODING HERE  #
    #####################   
    alignment = (aligned_seq1, aligned_seq2, align_score)
    return (alignment, score_matrix)



def print_score_matrix(s1,s2,mat):
    "Pretty print function for a score matrix."
    
    # Prepend filler characters to seq1 and seq2
    s1 = '-' + s1
    s2 = ' -' + s2
    
    # Print them around the score matrix, in columns of 5 characters
    print(''.join(['%5s' % aa for aa in s2])) # Convert s2 to a list of length 5 strings, then join it back into a string
    for i,row in enumerate(mat):               # Iterate through the rows of your score matrix (and keep count with 'i').
        vals = ['%5i' % val for val in row]    # Convert this row's scores to a list of strings.
        vals.insert(0,'%5s' % s1[i])           # Add this row's character from s2 to the front of the list
        print(''.join(vals))                   # Join the list elements into a single string, and print the line.



def print_alignment(a):
    "Pretty print function for an alignment (and alignment score)."
    
    # Unpack the alignment tuple
    seq1 = a[0]
    seq2 = a[1]
    score = a[2]
    
    # Check which positions are identical
    match = ''
    for i in range(len(seq1)): # Remember: Aligned sequences have the same length!
        match += '|' if seq1[i] == seq2[i] else ' ' # Fancy inline if-else statement. Use cautiously!
            
    # Concatenate lines into a list, and join them together with newline characters.
    print('\n'.join([seq1,match,seq2,'','Score = %i' % score]))



def save_alignment(a,f):
    "Saves two aligned sequences and their alignment score to a file."
    with open(f,'w') as out:
        out.write(a[0] + '\n') # Aligned sequence 1
        out.write(a[1] + '\n') # Aligned sequence 2
        out.write('Score: %i' % a[2]) # Alignment score


    
def save_score_matrix(m,f):
    "Saves a score matrix to a file in tab-separated format."
    with open(f,'w') as out:
        for row in m:
            vals = [str(val) for val in row]
            out.write('\t'.join(vals)+'\n')
    


def main(args = False):
    # Process arguments and load required data
    if not args: args = parse_args()
    
    sub_mat = load_substitution_matrix(args.substitution_matrix)
    seq1, seq2 = load_sequences(args.fasta)

    # Perform specified alignment
    strat = args.strategy
    gp = args.gap_penalty
    alignment, score_matrix = align(seq1, seq2, strat, sub_mat, gp)

    # If running in "verbose" mode, print additional output
    if args.verbose:
        print_score_matrix(seq1,seq2,score_matrix)
        print('') # Insert a blank line in between
        print_alignment(alignment)
    
    # Save results
    if args.align_out: save_alignment(alignment, args.align_out)
    if args.matrix_out: save_score_matrix(score_matrix, args.matrix_out)



if __name__ == '__main__':
    main()