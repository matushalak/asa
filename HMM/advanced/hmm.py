#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the Hidden Markov Models assignment in the Algorithms in Sequence Analysis course at the VU.

INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via corresponding Canvas assignment.

AUTHOR:
    <Matus Halak 2724858>
"""

import os.path as op

from os import makedirs
from math import log10
from hmm_utility import parse_args, load_fasta, load_tsv, print_trellis, print_params, serialize


'''
Overview:
Goal: Hidden Markov Model that can recognise Domain (D) and Linker (L) regions in protein sequences.
Simplified AA alphabet H (Hydrophobic), P (Polar), C (Charged). 
Normal AA alphabet could easily be converted into the simplified alphabet, 
but of course a HMM model with the full alphabet would perform better.

4 States: 
Q = {B(egin state), D(omain), L(inker), E(nd state)}

Alphabet:
∑ = {H(ydrophobic), P(olar), C(harged)}

Transition probabilities between 4 states:
A1 (initial prior for those probabilities, needs to be trained to have a good model)

Emission probabilities of the D & L states:
E1 (initial prior for those emission probabilities, needs ot be trained to be a good HMM)
'''


def viterbi(X,A,E):
    """Given a single sequence, with Transition and Emission probabilities,
    return the most probable state path, the corresponding P(X), and trellis."""

    allStates = A.keys() # rows
    emittingStates = E.keys() # rows
    L = len(X) + 2 # 1 begin state + seq len + 1 end state

    # Initialize
    # states are rows, columns are positions in sequence
    V = {k:[0] * L for k in allStates} # The Viterbi trellis
    # first & last row = start & end states
    # first & last column = start & end
    V['B'][0] = 1.

    # Middle columns
    for i,s in enumerate(X):
        for l in emittingStates:
            terms = [V[k][i] * A[k][l] for k in allStates]
            # i+1 because 0th state is start state
            V[l][i+1] = max(terms) * E[l][s]

    # Last column
    for k in allStates:
        term = V[k][i+1] * A[k]['E'] 
        if term > V['E'][-1]:
            V['E'][-1] = term
            pi = k # Last state of the State Path

    # FOR VITERBI ONLY: Trace back the State Path
    l = pi
    i = L-2
    while i:
        i -= 1
        for k in emittingStates:
            if V[k][i] * A[k][l] * E[l][X[i]] == V[l][i+1]:
                pi = k + pi
                l = k
                break

    P = V['E'][-1] # The Viterbi probability: P(X,pi|A,E)
    return(pi,P,V) # Return the state path, Viterbi probability, and Viterbi trellis



def forward(X,A,E):
    """Given a single sequence, with Transition and Emission probabilities,
    return the Forward probability and corresponding trellis."""

    allStates = A.keys()
    emittingStates = E.keys()
    L = len(X) + 2

    # Initialize
    F = {k:[0] * L for k in allStates}
    F['B'][0] = 1

    #####################
    # START CODING HERE #
    #####################
    # HINT: The Viterbi and Forward algorithm are very similar! 
    # Adapt the viterbi() function to account for the differences.

    # Middle columns
    for i,s in enumerate(X):
        for l in emittingStates:
            terms = [F[k][i] * A[k][l] for k in allStates]
            F[l][i+1] = sum(terms) * E[l][s]

    # Last column
    # Forward probability of entire sequence GIVEN HMM model
    F['E'][-1] = sum(
        # sum of final column means sum of 2nd to last column, since last column is 'end state' column
        [F[k][-2]*A[k]['E'] for k in allStates]
                     )
    #####################
    #  END CODING HERE  #
    #####################
    # this is the P(X) that we use in Baum-Welch!!!
    P = F['E'][-1] # The Forward probability: P(X|A,E)
    return(P,F)



def backward(X,A,E):
    """Given a single sequence, with Transition and Emission probabilities,
    return the Backward probability and corresponding trellis."""

    allStates = A.keys()
    emittingStates = E.keys()
    L = len(X) + 2

    # Initialize
    B = {k:[0] * L for k in allStates} # The Backward trellis
    for k in allStates:
        # last column is initialized with last -> end transition probabilities for all states
        # here 2nd to last because last column reserved for end state
        B[k][-2] = A[k]['E']

    #####################
    # START CODING HERE #
    #####################
    # Remaining columns
    # normally iteration is from L-1 since last column filled at initialization
    # here, we iterate from L-3 since last column corresponds to L-2
    for i in range(L-3, -1, -1): # until -1 because we want 0 to be included
        for k in allStates:
            terms = [B[l][i+1] * A[k][l] * E[l][X[i]] for l in emittingStates]
            B[k][i] = sum(terms)
    #####################
    #  END CODING HERE  #
    #####################

    # Probability that we are in begin state & observed entire sequence
    # this is the P(X) that we use in Baum-Welch!!!
    P = B['B'][0] # The Backward probability -- should be identical to Forward!
    return(P,B)



def baumwelch(set_X,A,E):
    """Given a set of sequences X and priors A and E,
    return the Sum Log Likelihood of X given the priors,
    along with the calculated posteriors for A and E."""

    allStates = A.keys()
    emittingStates = E.keys()
    symbols = [s for s in E[list(E)[0]]]
    
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

    # Iterate through all sequences in X
    SLL = 0 # Sum Log-Likelihood
    # this is the jth training sequence
    # since A_kl and E_k(s) are just sums over j, can just += A{^j}_kl and += E{^j}_k(s) to A_kl and E_k(s)
    row_sumA, row_sumE = {a:0 for a in allStates}, {e:0 for e in emittingStates}
    for jseq, X in enumerate(set_X): # set of training sequences
        # this takes case of Sum-log likelihood
        P,F = forward(X,A,E)  # Save both the forward probability and the forward trellis
        pB,B = backward(X,A,E) # Forward P == Backward P, so only save the backward trellis
        try:
            SLL += log10(P)
        except ValueError:
            breakpoint()
        #####################
        # START CODING HERE #
        #####################

        # Inside the for loop: Expectation
        # Calculate the expected transitions and emissions for the sequence.
        # Add the contributions to your posterior matrices.
        # Remember to normalize to the sequence's probability P!
        
        for k in allStates:
            for l in emittingStates:                    
                # add current Ajkl to Akl
                # breakpoint()
                new_A[k][l] += sum(
                    [(F[k][i] * A[k][l] * E[l][symbol] * B[l][i+1]) / P  
                        for i, symbol in enumerate(X)]
                                    )
                if jseq == len(set_X) - 1:
                    row_sumA[k] += new_A[k][l]
            # END STATE!!! - essential for future iterations so that sequence can end!
            new_A[k]['E'] += (F[k][len(X)] * A[k]['E']) / P
            if jseq == len(set_X) - 1:
                row_sumA[k] += new_A[k]['E']

        for l in emittingStates:                    
                for s in symbols:
                    # add current Ejks to Eks
                    new_E[l][s] += sum(
                        [((F[l][i+1] * B[l][i+1]) / P if symbol == s else 0) 
                        for i, symbol in enumerate(X)]
                                        )
                    if jseq == len(set_X) - 1:
                        row_sumE[l] += new_E[l][s]
    # Outside the for loop: Maximization
    # Normalize row sums to 1 (except END STATE row in the Transition matrix!)
    # normalize to probabilities
    # breakpoint()
    for fromState in new_A:
        if fromState != 'E':
            for toState in new_A[fromState]:            
                new_A[fromState][toState] /= row_sumA[fromState]

    for fromState in new_E:
        for emittedSymbol in new_E[fromState]:
            new_E[fromState][emittedSymbol] /= row_sumE[fromState]
    #####################
    #  END CODING HERE  #
    #####################
    return(SLL,new_A,new_E)



def main(args = False):
    "Perform the specified algorithm, for a given set of sequences and parameters."
    
    # Process arguments and load specified files
    if not args: args = parse_args()

    cmd = args.command            # viterbi, forward, backward or baumwelch
    verbosity = args.verbosity
    set_X, labels = load_fasta(args.fasta)  # List of sequences, list of labels
    A = load_tsv(args.transition) # Nested Q -> Q dictionary
    E = load_tsv(args.emission)   # Nested Q -> S dictionary
    
    def save(filename, contents):
        if args.out_dir:
            makedirs(args.out_dir, exist_ok=True) # Make sure the output directory exists.
            path = op.join(args.out_dir,filename)
            with open(path,'w') as f: f.write(contents)
        # Note this function does nothing if no out_dir is specified!



    # VITERBI
    if cmd == 'viterbi':
        for j,X in enumerate(set_X): # For every sequence:
            # Calculate the most probable state path, with the corresponding probability and matrix
            Q, P, T = viterbi(X,A,E)

            # Save and/or print relevant output
            label = labels[j]
            save('%s.path' % label, Q)
            save('%s.matrix' % label, serialize(T,X))
            save('%s.p' % label, '%1.2e' % P)
            print('>%s\n Path = %s' % (label,Q))
            if verbosity: print(' Seq  = %s\n P    = %1.2e\n' % (X,P))
            if verbosity >= 2: print_trellis(T, X)
            


    # FORWARD or BACKWARD
    elif cmd in ['forward','backward']:
        if cmd == 'forward':
            algorithm = forward
        elif cmd == 'backward':
            algorithm = backward

        for j,X in enumerate(set_X): # For every sequence:
            # Calculate the Forward/Backward probability and corresponding matrix
            P, T = algorithm(X,A,E)

            # Save and/or print relevant output
            label = labels[j]
            save('%s.matrix' % label, serialize(T,X))
            save('%s.p' % label, '%1.2e' % P)
            if verbosity >= 2:
                print('\n>%s\n P = %1.2e\n' % (label,P))
                print_trellis(T, X)
            elif verbosity: print('>%-10s\tP = %1.2e' % (label,P))



    # BAUM-WELCH TRAINING
    elif cmd == 'baumwelch':
        # Initialize
        i = 1
        i_max = args.max_iter
        threshold = args.conv_thresh

        current_SLL, A, E = baumwelch(set_X,A,E)
        if verbosity: print('Iteration %i, prior SLL = %1.2e' % (i,current_SLL))
        if verbosity >= 2: print_params(A,E)
        
        last_SLL = current_SLL - threshold - 1 # Iterate at least once

        # Iterate until convergence or limit
        while i < i_max and current_SLL - last_SLL > threshold:
            i += 1
            last_SLL = current_SLL

            # Calculate the Sum Log-Likelihood of X given A and E,
            # and update the estimates (posteriors) for A and E.
            current_SLL, A, E = baumwelch(set_X,A,E)

            if verbosity: print('Iteration %i, prior SLL = %1.2e' % (i,current_SLL))
            if verbosity >= 2: print_params(A,E)

        converged = current_SLL - last_SLL <= threshold
        try:
            final_SLL = sum([log10(forward(X,A,E)[0]) for X in set_X])
        except ValueError:
            final_SLL = 0

        # Save and/or print relevant output
        save('SLL','%1.2e\t%i\t%s' % (final_SLL, i, converged))
        save('posterior_A',serialize(A))
        save('posterior_E',serialize(E))
        if verbosity: print('========================================\n')

        if converged:
            print('Converged after %i iterations.' % i)
        else:
            print('Failed to converge after %i iterations.' % i_max)

        if verbosity:
            print('Final SLL: %1.2e' % final_SLL)
            print('Final parameters:')
            print_params(A,E)



if __name__ == '__main__':
	main()