def encode(INPUT:str):
    # annoying but takes care of spaces
    add = '$' if ' ' not in INPUT else ' '
    if ' ' in INPUT:
        INPUT = INPUT.replace(' ', '$')
    INPUT += add
     
    shifts = []
    for i, l in enumerate(INPUT):
        shifts.append(
            INPUT[i:] + INPUT[:i]
            )
    # sorted on list of strings sorts them alphabetically!!!!
    shifts = sorted(shifts)
    # for s in shifts: print('\n', s)
    BWT =  (''.join(s[-1] for s in shifts),
            shifts.index(INPUT))
    print(BWT)
    return BWT

def decode(CODE:tuple[str, int]):
    last_col, row_original = CODE
    decoded = ''
    # The first column can be acquired by sorting the last column.
    first_col = sorted(last_col)
    #For every row of the table: Symbols in the first column follow on symbols in the last column, in the same way they do in the input string.
    lc_ranked = {l:[] for l in set(last_col)}
    fc_ranked = {l:[] for l in set(first_col)}

    for i, letter in enumerate(last_col):
        lc_ranked[letter].append(i)
    for i, letter in enumerate(first_col):
        fc_ranked[letter].append(i)
    
    F, L = first_col[0], last_col[0]
    for _ in range(len(last_col)):
        # print(F, L)
        if F == first_col[0]:
            decoded += L
            lc_index = lc_ranked[L].index(0)
        else:
            decoded += L
            lc_index = lc_ranked[L].index(F)

        # look for ELEMENT at same index in L, and find its FIRST INDEX in F
        F = fc_ranked[L][lc_index] # F is index now
        L = last_col[F]
    
    print(OUT := decoded.replace('$', ' ').strip()[::-1])
    return OUT

decode(encode('banana'))
decode(encode('bananabar'))
decode(encode('mississipi'))
decode(encode('Yellew Mellow'))
decode(encode('I saw the #sign and it 0opened up m% Eyes I'))