#cython: language_level=3
import re
from cpython cimport array
import array

IUPACnum = {}
IUPACcomp = {}
IUPACregex = {}

# IUPAC code for nucleotides and ambiguities
A = 1
C = 2
G = 4
T = 8
keys = [
    'A','C','G','T','U',
    'M','R','W','S','Y','K',
    'V','H','D','B','N','X'
]
complements = [
    'T','G','C','A','A',
    'K','Y','W','S','R','M',
    'B','D','H','V','N','X'
]
integers = [
    A,C,G,T,T,
    A|C,A|G,A|T,C|G,C|T,G|T,
    A|C|G,A|C|T,A|G|T,C|G|T,A|C|G|T,0
]
int_to_IUPAC = [
    'X','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'
]
# REGEX representation of IUPAC nucleotides
regex = [
    'A','C','G','T','T',
    '[AC]','[AG]','[AT]','[CG]','[CT]','[GT]',
    '[ACG]','[ACT]','[AGT]','[CGT]','.','-'
]

# Illumina quality scores
quality_scores = [
    '!','"','#','$','%','&',
    "'",'(',')','*','+',
    ',','-','.','/','0',
    '1','2','3','4','5',
    '6','7','8','9',':',
    ';','<','=','>','?',
    '@','A','B','C','D',
    'E','F','G','H','I','J'
]


for k,i,c,r in zip(keys, integers, complements, regex):
    IUPACregex[k] = r
    IUPACnum[k] = i
    IUPACnum[k.lower()] = i
    IUPACcomp[k] = c
    IUPACcomp[k.lower()] = c

def to_regex(sequence):
    """Converts an IUPAC-formatted string to a regex string"""
    return ''.join([IUPACregex[i] for i in sequence.upper()])

def _comp_array():
    comp_array =  [0]*256
    for a,b in zip(keys, complements):
        comp_array[ord(a)] = ord(b)
        comp_array[ord(a.lower())] = ord(b.lower())
    
    return comp_array

def _int_array():
    int_array =  [0]*256
    for k,i in zip(keys, integers):
        int_array[ord(k)] = i
        int_array[ord(k.lower())] = i
    
    return int_array

def _qual_array():
    qual_array =  [int(i*.33) for i in list(range(256))]
    for i,q in enumerate(quality_scores):
        qual_array[ord(q)] = i
        qual_array[ord(q.lower())] = i
    
    return qual_array

INT_ARRAY = array.array('i',_int_array())
QUAL_ARRAY = array.array('i',_qual_array())
COMP_ARRAY = array.array('i',_comp_array())

# All IUPAC triplets that unambiguously translate to a single amino acid
codons = [
    'ATT','ATC','ATA',                   # Isoleucine (I)
    'CTT','CTC','CTA','CTG','TTA','TTG', # Leucine (L)
    'GTT','GTC','GTA','GTG',             # Valine (V)
    'TTT','TTC',                         # Phenylalanine (F)
    'ATG',                               # Methionine (M)
    'TGT','TGC',                         # Cysteine (C)
    'GCT','GCC','GCA','GCG',             # Alanine (A)
    'GGT','GGC','GGA','GGG',             # Glycine (G)
    'CCT','CCC','CCA','CCG',             # Proline (P)
    'ACT','ACC','ACA','ACG',             # Threonine (T)
    'TCT','TCC','TCA','TCG','AGT','AGC', # Serine (S)
    'TAT','TAC',                         # Tyrosine (Y)
    'TGG',                               # Tryptophan (W)
    'CAA','CAG',                         # Glutamine (Q)
    'AAT','AAC',                         # Asparagine (N)
    'CAT','CAC',                         # Histidine (H)
    'GAA','GAG',                         # Glutamic acid (E)
    'GAT','GAC',                         # Aspartic acid (D)
    'AAA','AAG',                         # Lysine (K)
    'CGT','CGC','CGA','CGG','AGA','AGG', # Arginine (R)
    'TAA','TAG','TGA'                    # Stop codon (-)
]
# single-letter IUPAC amino acids for each entry in codons
aminos = [
    'I','I','I',
    'L','L','L','L','L','L',
    'V','V','V','V',
    'F','F',
    'M',
    'C','C',
    'A','A','A','A',
    'G','G','G','G',
    'P','P','P','P',
    'T','T','T','T',
    'S','S','S','S','S','S',
    'Y','Y',
    'W',
    'Q','Q',
    'N','N',
    'H','H',
    'E','E',
    'D','D',
    'K','K',
    'R','R','R','R','R','R',
    '-','-','-'
]

#def _comp_array():
#    nuc1 =  [0]*256
#    nuc2 =  [0]*256
#    nuc3 =  [0]*256
#    codon_array = [nuc1, nuc2, nuc3]
#    for codon,amino in zip(codons, aminos):
#        n1, n2, n3 = codon
#        codon_array[ord(n1)][ord(n2)][ord(n3)]
#        comp_array[ord(a)] = ord(b)
#        comp_array[ord(a.lower())] = ord(b)
#    
#    return 

# All codons with IUPAC ambiguities that code for a single amino acid
codons_ambiguous = [
    'ATY','ATW','ATM','ATD',             # Isoleucine (I)
    'CTM','CTR','CTW','CTS','CTY','CTK', # Leucine (L)
    'CTV','CTH','CTD','CTB','CTN',
    'YTA','YTG','TTR','YTR',
    'GTM','GTR','GTW','GTS','GTY','GTK', # Valine (V)
    'GTV','GTH','GTD','GTB','GTN',
    'TTY',                               # Phenylalanine (F)
    'TGY',                               # Cysteine (C)
    'GCM','GCR','GCW','GCS','GCY','GCK', # Alanine (A)
    'GCV','GCH','GCD','GCB','GCN',
    'GGM','GGR','GGW','GGS','GGY','GGK', # Glycine (G)
    'GGV','GGH','GGD','GGB','GGN',
    'CCM','CCR','CCW','CCS','CCY','CCK', # Proline (P)
    'CCV','CCH','CCD','CCB','CCN',
    'ACM','ACR','ACW','ACS','ACY','ACK', # Threonine (T)
    'ACV','ACH','ACD','ACB','ACN',
    'TCM','TCR','TCW','TCS','TCY','TCK', # Serine (S)
    'TCV','TCH','TCD','TCB','TCN','AGY',
    'TAY',                               # Tyrosine (Y)
    'CAR',                               # Glutamine (Q)
    'AAY',                               # Asparagine (N)
    'CAY',                               # Histidine (H)
    'GAR',                               # Glutamic acid (E)
    'GAY',                               # Aspartic acid (D)
    'AAR',                               # Lysine (K)
    'CGM','CGR','CGW','CGS','CGY','CGK', #Arginine (R)
    'CGV','CGH','CGD','CGB','CGN',
    'AGR','MGA','MGG','MGR',
    'TAR','TRA'                          # Stop codon (-)
]

# IUPAC amino acid that matches each element of codons_ambiguous
aminos_ambiguous = [
    'I','I','I','I',
    'L','L','L','L','L','L',
    'L','L','L','L','L',
    'L','L','L','L',
    'V','V','V','V','V','V',
    'V','V','V','V','V',
    'F',
    'C',
    'A','A','A','A','A','A',
    'A','A','A','A','A',
    'G','G','G','G','G','G',
    'G','G','G','G','G',
    'P','P','P','P','P','P',
    'P','P','P','P','P',
    'T','T','T','T','T','T',
    'T','T','T','T','T',
    'S','S','S','S','S','S',
    'S','S','S','S','S','S',
    'Y',
    'Q',
    'N',
    'H',
    'E',
    'D',
    'K',
    'R','R','R','R','R','R',
    'R','R','R','R','R',
    'R','R','R','R',
    '-','-'
]

quality_scores = [
    '!','"','#','$','%','&',
    "'",'(',')','*','+',
    ',','-','.','/','0',
    '1','2','3','4','5',
    '6','7','8','9',':',
    ';','<','=','>','?',
    '@','A','B','C','D',
    'E','F','G','H','I','J'
]

CODONhash = {}
for c,a in zip(codons, aminos):
    CODONhash[c] = a

AMBIGhash = {}
for c,a in zip(codons_ambiguous, aminos_ambiguous):
    AMBIGhash[c] = a

# Define a collection of Typed Memoryviews
cdef int[:] IUPACint = INT_ARRAY
cdef int[:] QUALint = QUAL_ARRAY
cdef int[:] COMPint = COMP_ARRAY

cpdef array.array nuc_to_int(str nuc_string, str qual_string='', int qualmask=12):
    """Converts a string of IUPAC-encoded nucleotides
    to an integer array."""
    global IUPACint
    global QUALint
    cdef:
        int s, qual
        Py_ssize_t i
        list nuc_as_int
    
    if qual_string == '':
        qual_string = 'J'*len(nuc_string)
    
    nuc_as_int = [IUPACint[s] for s in nuc_string]
    for i in range(len(nuc_string)):
        qual = QUALint[qual_string[i]]
        if qual <= qualmask:
            nuc_as_int[i] = 15
    
    return array.array('i',nuc_as_int)

cpdef bint is_homopolymer(str string, float threshold=0.8):
    """Returns whether a single character composes > threshold
    of a string."""
    cdef str n
    cdef int count_n, total_count, string_length, thresh_length
    
    string = string.upper()
    string_length = len(string)
    thresh_length = int(round(string_length * threshold))
    if string_length == 0:
        return True
    
    total_count = 0
    for n in ['A','T','G','C']:
        count_n = string.count(n)
        if count_n >= thresh_length:
            return True
        else:
            total_count += count_n
            if total_count > string_length - thresh_length:
                # Enough subthreshold nucleotides were found
                return False
    
    return False

cpdef double quality_score(str qual_string):
    """Takes an Illumina format quality string
    and returns the average score"""
    global QUALint
    if len(qual_string) == 0:
        return 0
    
    return float(sum([QUALint[s] for s in qual_string]))/len(qual_string)

cpdef str rc(str sequence):
    """Returns reverse complement of a nucleotide string."""
    global COMPint
    cdef:
        int s, str_len, i
        array.array rc_array
        str revcomp
    
    str_len = len(sequence)
    rc_array = array.array('i',[0]*str_len)
    i = str_len - 1
    for s in sequence:
        rc_array[i] = COMPint[s]
        i -= 1
    
    revcomp = ''.join(map(chr, rc_array))
    return revcomp

cpdef str complement(str sequence):
    """Returns complement of a nucleotide string."""
    global COMPint
    cdef:
        int s, str_len, i
        array.array comp_array
        str comp
    
    str_len = len(sequence)
    comp_array = array.array('i',[0]*str_len)
    i = 0
    for s in sequence:
        comp_array[i] = COMPint[s]
        i += 1
    
    comp = ''.join(map(chr, comp_array))
    return comp

cpdef int IUPACham(array.array a, array.array b, int stop_at=-1):
    """
    Returns the Hamming distance between two IUPAC numeric arrays
    """
    cdef int ham = -1
    if a == b: # Equality check to avoid unnecessary calculations
        return ham
    
    cdef int len_A, len_B, i
    cdef int[:] A = a
    cdef int[:] B = b
    len_A, len_B = len(A), len(B)
    if len_A != len_B:
        return ham
    
    if stop_at == -1: # No maximum distance was assigned
        stop_at = len_A
    
    ham = 0
    for i in range(len_A):
        x, y = A[i], B[i]
        if not x & y: # Bitwise-AND determines if two IUPAC characters match
            ham += 1
            if ham > stop_at: # Hamming distance has exceeded the maximum allowed
                return ham
    
    return ham

cpdef bint oligo_match(array.array a, array.array b, float mm_rate, int min_oligomer=8):
    """
    Returns a bool indicating if a sufficiently close match was found.
    For long oligos, allow a match <=mm_rate to any prefix of a >=half a's length
    """
    cdef int minmatch, max_ham, sub_ham
    cdef int ham = -1
    if a == b: # Equality check to avoid unnecessary calculations
        return True
    
    cdef int len_A, len_B, i
    cdef int[:] A = a
    cdef int[:] B = b
    len_A, len_B = len(A), len(B)
    if len_A != len_B:
        return False
    
    minmatch = min(len_A, min_oligomer)
    max_ham = int(len_A*mm_rate)
    ham = 0
    for i in range(len_A):
        x, y = A[i], B[i]
        if not x & y: # Bitwise-AND determines if two IUPAC characters match
            ham += 1
            if ham > max_ham:
                return False
        
        if i >= minmatch:
            sub_ham = int(i*mm_rate)
            if ham <= sub_ham:
                return True
    
    return ham <= max_ham

cpdef (int, int) best_sliding_fit(
        array.array query, array.array trim,
        int monomer=-1, int minmatch=5, double mm_rate=0.06,
        int maxlen=120):
    '''
    Slides trim along query and returns a double of
    (end_position, hamming_distance) for the longest possible match of trim in query.
    Allows for 3' extension of a monomer of unknown length
    '''
    global NULLvalues
    cdef:
        int q_len, t_len, best_pos, best_ham, min_oligomer, i, max_mismatch, oligo_counter, ocheck, n, a_i, spacer
        int l_a, l_b, ham, x, y, j
        int[:] a, b
    
    best_pos = 0
    best_ham = -1
    t_len = len(trim)
    q_len = len(query)
    if t_len == 0:
        return best_pos, best_ham
    
    cdef int[:] Q = query
    cdef int[:] T = trim
    min_oligomer = min(3, minmatch)
    i = minmatch
    while i <= min(q_len,maxlen):
        max_mismatch = int(mm_rate*min(i,t_len))
        spacer = i-t_len # Calculate whether the trim array is shorter than i
        if spacer > 0:
            a = Q[spacer:i] # Slice the query up to position i
        else:
            a = Q[:i]
        
        b = T[-i:]
        ham = -1
        l_a = len(a)
        l_b = len(b)
        if l_a == l_b:
            ham = 0
            for j in range(l_a):
                x = a[j]
                y = b[j]
                if not x & y: # Bitwise-AND determines if two IUPAC characters match
                    ham += 1
                    if ham > max_mismatch: # Hamming distance has exceeded the maximum allowed
                        break
        
        if ham == -1:
            print("WARNING: @{}: {} {} length mismatch".format(i, len(a),len(b)))
            return best_pos, best_ham
        
        if ham <= max_mismatch:
            if a[-1] & b[-1] > 0: # Require a terminal match
                if i == q_len:
                    return i, ham
                
                if monomer != -1: # Extend the match for as long as possible (<= max_mismatch)
                    oligo_counter = 1
                    ocheck = a[-oligo_counter]
                    while ocheck & monomer > 0 and oligo_counter < min_oligomer:
                        oligo_counter += 1
                        ocheck = a[-oligo_counter]
                    
                    best_pos = i
                    best_ham = ham
                    a_i = Q[i]
                    n = i
                    while ham <= max_mismatch: # Keep track of how many mismatches are encountered during monomer extension
                        if not a_i & monomer:
                            ham += 1
                            oligo_counter = 0
                        else: # If a_i matches the monomer, move best_pos up 1
                            oligo_counter += 1
                            if oligo_counter >= min_oligomer:
                                # best_pos = n+1
                                best_pos = n+1
                                best_ham = ham
                        
                        n += 1
                        max_mismatch = int(mm_rate*n)
                        if n == q_len:
                            return n, ham
                        
                        a_i = Q[n]
                    
                    if best_pos > i:
                        i = best_pos + 1
                else:
                    best_pos = i
                    best_ham = ham

        i += 1

    return best_pos, best_ham

cdef bint reads_are_reversed(int trimtype1, int trimtype2):
    '''
    If a pair of trims is concordant, they can either be in
    the sense (read1 forward, read2 reverse) or
    antisense (read2 forward, read1 reverse) orientation.
    Return boolean of whether the read pair is reversed.
    TRIMTYPES:
        -1=None, 0=S5, 1=E5, 2=S3, 3=E3
    '''
    if trimtype1 == -1:
        if trimtype2 == -1:
            return False
        elif trimtype2 == 0:
            return True
        elif trimtype2 == 1:
            return False
    elif trimtype1 == 0:
        return False
    elif trimtype1 == 1:
        return True
    elif trimtype2 == -1:
        if trimtype1 == -1:
            return False
        elif trimtype1 == 0:
            return False
        elif trimtype1 == 1:
            return True
    elif trimtype2 == 0:
        return True
    elif trimtype2 == 1:
        return False
    return False


cdef str trim_readstring(str readstring, int pos, int trimtype, bint qual=False, bint reverse=False):
    """Executes the trimming of a read given the type of trim (S5,S3,E5,E3)
    and the trimming position. Orients read in sense to the RNA strand.
    If pair, then the reverse complement is returned."""
    cdef str trim
    if trimtype == 0 or trimtype == 1: # Head trimming wins
        trim = readstring[pos:]
        if trimtype == 1: # Read is flipped
            if qual:
                trim = trim[::-1]
            else:
                trim = rc(trim)
    elif trimtype == 2 or trimtype == 3: # Tail trimming wins
        trim = readstring[:-pos]
        if trimtype == 2: # Read is flipped
            if qual:
                trim = trim[::-1]
            else:
                trim = rc(trim)
    else: # No possible trimming arrangements were provided
        if reverse:
            if qual:
                return readstring[::-1]
            else:
                return rc(readstring)
        else:
            return readstring
    
    return trim

cdef str get_umilabel(str readstring, int pos, int trimtype, str umi, (int, int) umi_range):
    """Extracts the UMI sequence from a specified position on the trimmed string"""
    cdef str umi_string
    umi_string = ""
    if trimtype == -1:
        return ""
    elif trimtype in [0,2] and umi != 'S':
        return ""
    elif trimtype in [1,3] and umi != 'E':
        return ""
    
    if trimtype == 0 or trimtype == 1: # Head trimming wins
        umi_string = readstring[:pos][umi_range[0]:umi_range[1]]
        if trimtype == 1: # Read is flipped
            umi_string = rc(umi_string)
    elif trimtype == 2 or trimtype == 3: # Tail trimming wins
        umi_string = readstring[-pos:][-umi_range[1]:-umi_range[0]]
        if trimtype == 2: # Read is flipped
            umi_string = rc(umi_string)
    else: # No possible trimming arrangements were provided
        return ""
    
    if len(umi_string) == 0:
        return ""
    else:
        return "_UMI="+umi_string

cdef str get_label(str readstring, int pos, int trimtype):
    """Returns a string labeling the end tag type and length."""
    cdef:
        str label = '',
        int str_len, trim_len
    if trimtype == -1:
        return label
    
    str_len = len(readstring)
    trim_len = str_len - pos
    if trimtype == 0 or trimtype == 2:
        label = 'S'
    elif trimtype == 1 or trimtype == 3:
        label = 'E'
    
    label += str(str_len-trim_len)
    return label

cdef bint is_flipped(int trimtype):
    return trimtype == 1 or trimtype == 2

cpdef (int, int, int) complementary_trim(
        array.array nucarray, int trimtype, 
        array.array S5array, int S5monomer, 
        array.array S3array, int S3monomer,
        array.array E5array, int E5monomer, 
        array.array E3array, int E3monomer, 
        int minstart, int minend, float mm_rate, int maxstart, int maxend):
    """Checks whether the reverse complement of a trimtype
    exists in a numeric nucleotide array.
    Returns a triple of (pos, ham, trimtype)"""
    cdef:
        int pos, ham
        int comptype
        array.array nucarray_rev
    
    pos = 0
    ham = -1
    comptype = -1
    if trimtype == 0:
        nucarray_rev = nucarray[::-1]
        pos,ham = best_sliding_fit(nucarray_rev, S3array, S3monomer, minstart, mm_rate, maxstart)
        comptype = 2
    elif trimtype == 2:
        pos,ham = best_sliding_fit(nucarray, S5array, S5monomer, minstart, mm_rate, maxstart)
        comptype = 0
    elif trimtype == 1:
        nucarray_rev = nucarray[::-1]
        pos,ham = best_sliding_fit(nucarray_rev, E3array, E3monomer, minend, mm_rate, maxend)
        comptype = 3
    elif trimtype == 3:
        pos,ham = best_sliding_fit(nucarray, E5array, E5monomer, minend, mm_rate, maxend)
        comptype = 1

    if ham == -1:
        comptype = -1
    
    return pos, ham, comptype

cpdef str collapse_reads(str string1, str string2, str qual1, str qual2, double mm_rate=0.06, int qualmask=12):
    """Given two strings of the same or differing lengths,
    give the consensus model, resolving ambiguities if possible.
    Ex: collapse_reads('AGGCNNTC','ARGCTATCGGCAATA')
    AGGCTATCGGCAATA
    """
    cdef:
        int l1, l2, bigL, smallL, short_len, nextnuc, a, b, ldiff
        array.array longer, shorter, collapsed
        str longstring
        int i, l, s, mismatches
        double mm_allowed
    
    mismatches = 0
    l1, l2 = len(string1), len(string2)
    if l1 >= l2:
        bigL, smallL = l1, l2
        longer = nuc_to_int(string1, qual1, qualmask)
        shorter = nuc_to_int(string2, qual2, qualmask)
        longstring = string1
    else:
        bigL, smallL = l2, l1
        longer = nuc_to_int(string2, qual2, qualmask)
        shorter = nuc_to_int(string1, qual1, qualmask)
        longstring = string2

    if len(shorter) == 0:
        # shorter necessarily collapses into longer
        return longstring
    
    collapsed = array.array('i',[0]*bigL)
    mm_allowed = mm_rate*smallL
    # (1) Left-justify
    for i in range(bigL): # Get a consensus of all shared nucleotides
        if i >= smallL:
            a = 15
        else:
            a = shorter[i]
        
        b = longer[i]
        nextnuc = a & b
        if nextnuc == 0:
            mismatches += 1
            if mismatches > mm_allowed:
                # The strings can't be collapsed
                if bigL == smallL:
                    return ''

                break

        collapsed[i] = nextnuc

    if nextnuc == 0: # Left-justified alignment failed
        # (2) Right-justify
        for i in range(bigL): # Get a consensus of all shared nucleotides
            ldiff = bigL-smallL
            if i >= ldiff:
                a = shorter[i-ldiff]
            else:
                a = 15
            
            b = longer[i]
            nextnuc = a & b
            if nextnuc == 0:
                mismatches += 1
                if mismatches > mm_allowed:
                    # The strings can't be collapsed
                    if bigL == smallL:
                        return ''
                    
                    return ''
                

            collapsed[i] = nextnuc
    
    return ''.join([int_to_IUPAC[i] for i in collapsed])

def terminal_trim(
        str mate1, str qual1,
        str mate2, str qual2, 
        array.array S5array, int S5monomer, 
        array.array S3array, int S3monomer, 
        array.array E5array, int E5monomer, 
        array.array E3array, int E3monomer, 
        str strand, int minstart, int minend, int minlen, double minqual, int qualmask, float mm_rate,
        str umi, (int,int) umi_range, int maxstart, int maxend):
    '''
    Trims the most well-supported adapter sequence(s) at the end(s)
    of RNA seq reads, based on the expected structure of the double-stranded cDNA.
      S5--------------E3
      S3--------------E5

    The following adapter pairs are possible:
      None, S5|E5
      S5,   None|S3|E5
      S3,   S5
      E5,   None|S5|E3
      E3,   E5

    Returns trimmed (mate1, mate2, label).
    '''
    cdef:
        int pos1, pos2, ham1, ham2, trimtype1, trimtype2
        int E5pos1, E5ham1, E3pos1, E3ham1, S5pos1, S5ham1, S3pos1, S3ham1
        array.array mate1array, mate2array, mate1array_rev
        str trim1, qtrm1, label1, trim2, qtrm2, label2, umilabel
        bint flipped1, flipped2, reverse
    
    # Convert input strings to numeric IUPAC arrays
    mate1array = nuc_to_int(mate1, qual1, qualmask)
    mate2array = nuc_to_int(mate2, qual2, qualmask)
    umilabel = ""
    reverse = strand == 'reverse'
    # (1) Calculate 2 possible forward trims on mate1
    pos1 = 0
    ham1 = len(mate1array)
    trimtype1 = -1
    if strand != 'reverse': # Check for 5P adapter (sense orientation)
        S5pos1,S5ham1 = best_sliding_fit(mate1array, S5array, S5monomer, minstart, mm_rate, maxstart)
        if S5ham1 != -1: # A match was found
            pos1 = S5pos1
            ham1 = S5ham1
            trimtype1 = 0
    
    if strand != 'forward': # Check for 3P adapter (antisense orientation)
        # (2) Check if E5 is a better match
        E5pos1,E5ham1 = best_sliding_fit(mate1array, E5array, E5monomer, minend, mm_rate, maxend)
        if E5ham1 != -1:
            if (E5pos1 > pos1) or (E5pos1 == pos1 and E5ham1 < ham1): # Improved match
                pos1 = E5pos1
                ham1 = E5ham1
                trimtype1 = 1
    
    if len(mate2) == 0: # SINGLE-END; check for tail matches
        trimtype2 = -1
        mate1array_rev = mate1array[::-1]
        if strand != 'forward': # Check for RC of 5P adapter (antisense orientation)
            S3pos1,S3ham1 = best_sliding_fit(mate1array_rev, S3array, S3monomer, minstart, mm_rate, maxstart)
            if S3ham1 != -1:
                if (S3pos1 > pos1) or (S3pos1 == pos1 and S3ham1 < ham1) or trimtype1 == 1: # Improved match
                    if trimtype1 == 1: # Can be shared with first trim
                        pos2 = S3pos1
                        ham2 = S3ham1
                        trimtype2 = 2
                    else: # Replaces first trim
                        pos1 = S3pos1
                        ham1 = S3ham1
                        trimtype1 = 2
                        pos2 = 0
                        ham2 = -1
                        trimtype2 = -1
        if strand != 'reverse': # Check for RC of 3P adapter (sense orientation)
            E3pos1,E3ham1 = best_sliding_fit(mate1array_rev, E3array, E3monomer, minend, mm_rate, maxend)
            if E3ham1 != -1:
                if (E3pos1 > pos1) or (E3pos1 == pos1 and E3ham1 < ham1) or trimtype1 == 0: # Improved match
                    if trimtype1 == 0: # Can be shared with first trim
                        pos2 = E3pos1
                        ham2 = E3ham1
                        trimtype2 = 3
                    else: # Replaces first trim
                        pos1 = E3pos1
                        ham1 = E3ham1
                        trimtype1 = 3
                        pos2 = 0
                        ham2 = -1
                        trimtype2 = -1
        
        trim1 = trim_readstring(mate1, pos1, trimtype1, qual=False, reverse=reverse)
        qtrm1 = trim_readstring(qual1, pos1, trimtype1, qual=True, reverse=reverse)
        umilabel = get_umilabel(mate1, pos1, trimtype1, umi, umi_range)
        label1 = get_label(mate1, pos1, trimtype1)
        flipped1 = is_flipped(trimtype1)
        if trimtype2 != -1: # A matching trim was found
            if trimtype1 == 0: # S5 + E3 trim
                trim1 = trim1[:-pos2]
                qtrm1 = qtrm1[:-pos2]
                label1 = label1+'E'+str(pos2)
            elif trimtype1 == 1: # E5 + S3 trim
                trim1 = trim1[pos2:]
                qtrm1 = qtrm1[pos2:]
                label1 = 'S'+str(pos2)+label1
        
        if len(trim1) < minlen: # trim1 length fail
            return('', '', None, None, label1+umilabel) # Return empty read with label
        else:
            qualscore1 = quality_score(qtrm1)
            if qualscore1 < minqual: # trim1 quality fail
                return('', '', None, None, label1+umilabel) # Return empty read with NO label
        
        return (trim1, qtrm1, None, None, label1+umilabel) # Checks passed, return trimmed read
    else: # PAIRED-END: mate2 exists, find the best mate2 trim
        # (3) Calculate 2 possible forward trims on mate2
        pos2 = 0
        ham2 = len(mate2array)
        trimtype2 = -1
        if strand != 'forward': # Check for mate2 5P adapter (antisense orientation of mate1)
            S5pos2,S5ham2 = best_sliding_fit(mate2array, S5array, S5monomer, minstart, mm_rate, maxstart)
            if S5ham2 != -1: # A match was found
                pos2 = S5pos2
                ham2 = S5ham2
                trimtype2 = 0

        if strand != 'reverse': # Check for mate2 3P adapter (sense orientation of mate1)
            E5pos2,E5ham2 = best_sliding_fit(mate2array, E5array, E5monomer, minend, mm_rate, maxend)
            if E5ham2 != -1:
                if (E5pos2 > pos2) or (E5pos2 == pos2 and E5ham2 < ham2): # Improved match
                    pos2 = E5pos2
                    ham2 = E5ham2
                    trimtype2 = 1

        # Generate the trimmmed sequences that represent the best match for mate1 and mate2
        trim1 = trim_readstring(mate1, pos1, trimtype1, qual=False, reverse=reverse)
        qtrm1 = trim_readstring(qual1, pos1, trimtype1, qual=True, reverse=reverse)
        label1 = get_label(mate1, pos1, trimtype1)
        umilabel = get_umilabel(mate1, pos1, trimtype1, umi, umi_range)
        flipped1 = is_flipped(trimtype1)
        
        trim2 = trim_readstring(mate2, pos2, trimtype2, qual=False, reverse=reverse)
        qtrm2 = trim_readstring(qual2, pos2, trimtype2, qual=True, reverse=reverse)
        label2 = get_label(mate2, pos2, trimtype2)
        if len(umilabel) == 0:
            umilabel = get_umilabel(mate2, pos2, trimtype2, umi, umi_range)
        
        flipped2 = is_flipped(trimtype2)
        
        # Determine which of the two trims above is better
        if trimtype1 == -1 and trimtype2 == -1: # No trimming was performed at all
            return (mate1, qual1, mate2, qual2, ''+umilabel)
        
        if pos1 > pos2 or (pos1 == pos2 and ham1 <= ham2): # Mate1's trim is the best supported (wins ties)
            # Check for complement of mate1's trim
            posC, hamC, trimtypeC = complementary_trim(mate2array, trimtype1, S5array, S5monomer, S3array, S3monomer, E5array, E5monomer, E3array, E3monomer, minstart, minend, mm_rate, maxstart, maxend)
            if trimtypeC != -1: # The complement putatively exists
                trimC = trim_readstring(mate2, posC, trimtypeC, qual=False, reverse=reverse)
                qtrmC = trim_readstring(qual2, posC, trimtypeC, qual=True, reverse=reverse)
                if len(umilabel) == 0:
                    umilabel = get_umilabel(mate2, posC, trimtypeC, umi, umi_range)
                
                labelC = get_label(mate2, posC, trimtypeC)
                flippedC = is_flipped(trimtypeC)
                
                collapsed = collapse_reads(trim1, trimC, qtrm1, qtrmC, mm_rate, qualmask)
                if collapsed != '': # Complementary trimming is compatible, output solo read
                    if len(qtrm1) == len(collapsed): # mate1 was the longer sequence
                        return (collapsed, qtrm1, None, None, labelC+umilabel)
                    else: # mate2 was the longer sequence
                        return (collapsed, qtrmC, None, None, label1+umilabel)
                else: # Pick a solo output with the higher quality score
                    qualscore1 = quality_score(qtrm1)
                    qualscoreC = quality_score(qtrmC)
                    if len(trim1) < minlen or qualscore1 < minqual: # trim1 fail
                        if len(trimC) < minlen or qualscoreC < minqual: # trimC fail
                            return('', '', None, None, '') # Return empty read
                        else: # trimC pass
                            return(trimC, qtrmC, None, None, labelC+umilabel)
                    elif len(trimC) < minlen or qualscoreC < minqual: # trim2 pass, trimC fail
                        return(trim1, qtrm1, None, None, label1)
                    else: # both pass
                        if qualscore1 >= qualscoreC:
                            return(trim1, qtrm1, None, None, label1+umilabel)
                        else:
                            return(trimC, qtrmC, None, None, labelC+umilabel)
            else: # Could not find complementary trim
                if len(trim1) < minlen or quality_score(qtrm1) < minqual: # Trimmed length is too short to be usable
                    return(mate2, qual2, None, None, ''+umilabel) # Return empty read
            
            # Check if the original trim was compatible)
            if trimtype1 == trimtype2: # Incompatible two-headed or two-tailed trim
                trim2, qtrm2, label2, trimtype2 = mate2, qual2, '', -1
        else: # Mate2's trim is the best supported.
            # Check for complement of mate2's trim
            posC, hamC, trimtypeC = complementary_trim(mate1array,trimtype2, S5array, S5monomer, S3array, S3monomer, E5array, E5monomer, E3array, E3monomer, minstart, minend, mm_rate, maxstart, maxend)
            if trimtypeC != -1: # The complement putatively exists
                trimC = trim_readstring(mate1, posC, trimtypeC, qual=False, reverse=reverse)
                qtrmC = trim_readstring(qual1, posC, trimtypeC, qual=True, reverse=reverse)
                if len(umilabel) == 0:
                    umilabel = get_umilabel(mate1, posC, trimtypeC, umi, umi_range)
                
                labelC = get_label(mate1, posC, trimtypeC)
                flippedC = is_flipped(trimtypeC)
                
                collapsed = collapse_reads(trim2, trimC, qtrm2, qtrmC, mm_rate, qualmask)
                if collapsed != '': # Complementary trimming is compatible, output solo read
                    if len(qtrm2) == len(collapsed): # mate1 was the longer sequence
                        return (collapsed, qtrm2, None, None, labelC+umilabel)
                    else: # mate1 was the longer sequence
                        return (collapsed, qtrmC, None, None, label2+umilabel)
                else: # Pick the higher quality solo-mapper
                    qualscore2 = quality_score(qtrm2)
                    qualscoreC = quality_score(qtrmC)
                    if len(trim2) < minlen or qualscore2 < minqual: # trim2 fail
                        if len(trimC) < minlen or qualscoreC < minqual: # trimC fail
                            return('', '', None, None, '') # Return empty read
                        else: # trimC pass
                            return(trimC, qtrmC, None, None, labelC+umilabel)
                    elif len(trimC) < minlen or qualscoreC < minqual: # trim2 pass, trimC fail
                        return(trim2, qtrm2, None, None, label2+umilabel)
                    else: # both pass
                        if qualscore2 >= qualscoreC:
                            return(trim2, qtrm2, None, None, label2+umilabel)
                        else:
                            return(trimC, qtrmC, None, None, labelC+umilabel)
            else: # Could not find complementary trim
                if len(trim2) < minlen or quality_score(qtrm2) < minqual: # Trimmed length is too short to be usable
                    return(mate1, qual1, None, None, ''+umilabel) # Return empty read
            
            # Check if the original trim was compatible)
            if trimtype1 == trimtype2: # Incompatible two-headed or two-tailed trim
                trim1, qtrm1, label1, trimtype1 = mate1, qual1, '', -1

    if 'E' in label1:
        label = label2+label1
    else:
        label = label1+label2
    
    reverse = reads_are_reversed(trimtype1,trimtype2)
    if reverse:
        # mates 1 and 2 must be swapped
        trim1, qtrm1, label1, flipped1, trim2, qtrm2, label2, flipped2 = trim2, qtrm2, label2, flipped2, trim1, qtrm1, label1, flipped1

    if 'S' in label and 'E' in label: # Special case: mate1 S, mate2 E
        # Because trim_readstring() flipped both mates, one must be flipped back
        trim2, qtrm2 = rc(trim2), qtrm2[::-1]

    if label:
        if label[0] == 'E':
            trim2, qtrm2 = rc(trim2), qtrm2[::-1]
            if int(label[1:]) == len(mate1) and len(trim1) == len(mate1):
                # One fully trimmed and one fully untrimmed read exist. This is uninformative as an end
                label = ''

    return (trim1, qtrm1, trim2, qtrm2, label+umilabel)

# Importing a FASTA file as a genome object
cpdef import_genome(str genome_FASTA, str split_on=' ', bint keep_case=True, bint indexed=False):
    """Reads FASTA file to a dict."""
    cdef:
        str rawline, line, chromname, chromstring, firstchar, index
        long linelen, start_position, trimmed_length, untrimmed_length
        list current_lines = []
        dict genome = {}
        bint length_mismatch
    
    chromname = 'none'
    chromstring = ''
    trimmed_length = untrimmed_length = start_position = 0
    length_mismatch = False
    if indexed:
        index = open(genome_FASTA+'.fai','r').read()
    else:
        index = ''
    
    genome_file = open(genome_FASTA)
    rawline = genome_file.readline()
    while rawline:
        line = rawline.rstrip()
        linelen = len(line)
        if linelen == 0:
            rawline = genome_file.readline()
            continue
        
        firstchar = line[0]
        if firstchar == '>':
            if chromname != 'none':
                chromstring = ''.join(current_lines)
                genome[chromname] = chromstring
                if not indexed:
                    index += '{}\t{}\t{}\t{}\t{}\n'.format(chromname, len(chromstring), start_position, trimmed_length, untrimmed_length)
            
            chromname = line[1:].split(split_on)[0]
            start_position = genome_file.tell()
            trimmed_length = 0
            length_mismatch = False
            chromstring = ''
            current_lines = []
            rawline = genome_file.readline()
            continue
        elif not indexed:
            if length_mismatch:
                print("Indexing error: [{}] line length mismatch. Indexing ignored.".format(genome_file.tell()))
                indexed = True
            
            if trimmed_length == 0:
                trimmed_length = len(line)
                untrimmed_length = len(rawline)
            elif trimmed_length != len(line):
                length_mismatch = True
        
        if not keep_case:
            line = line.upper()
        
        current_lines.append(line)
        rawline = genome_file.readline()

    chromstring = ''.join(current_lines)
    genome[chromname] = chromstring
    if not indexed:
        index += '{}\t{}\t{}\t{}\t{}\n'.format(chromname, len(chromstring), start_position, trimmed_length, untrimmed_length)
    
    genome_file.close()
    return genome, index


def generate_softbridges(dict genome_dict, int minlen, int maxlen):
    """From a genome dict, yields one BED12 line for each start/stop
    of a softmasked region of the FASTA file, demarcated by lowercase letters.
    """
    cdef:
        bint soft_toggle
        str chromname, chromstring, c, out_string
        int current_pos, start_pos, end_pos, sb_length
        bint lowercase
        Py_ssize_t i
    
    soft_toggle = False
    current_pos = start_pos = end_pos = sb_length = 0
    for chromname in sorted(list(genome_dict.keys())):
        soft_toggle = False
        current_pos = 0
        start_pos = 0
        end_pos = 0
        chromstring = genome_dict[chromname]
        for i in range(len(chromstring)):
            c = chromstring[i]
            lowercase = c.islower()
            if lowercase:
                if not soft_toggle: # Start a softbridge
                    soft_toggle = True
                    start_pos = current_pos
            elif soft_toggle: # End a softbridge
                soft_toggle = False
                end_pos = current_pos
                sb_length = end_pos - start_pos
                if sb_length >= minlen and sb_length <= maxlen:
                    out_string = '{}\t{}\t{}\t.\t0.01\t.\t0\t0\t204,204,180\t1\t{}\t0\t0.01\tsoftbridge\t..\n'.format(
                        chromname, start_pos, end_pos, sb_length
                    )
                    yield out_string
            
            current_pos += 1


def number_chromosomes(genome):
    """Returns a sorted index of chromosome starting positions in a genome."""    
    running_count = 0
    chromosome_number = {}
    for c in sorted(genome.keys()):
        chromosome_number[c] = running_count
        running_count += 1
    
    return chromosome_number

def translate(codon):
    """Looks up a nucleotide triplet in the codon hashtables."""
    if len(codon) == 3:
        c = CODONhash.get(codon, '?')
    else:
        return ''
    
    if c == '?':
        return AMBIGhash.get(codon,'X')
    else:
        return c

def longest_orf(sequence,allow_truncation=True):
    """Locates the longest open reading frame.
    
    Outputs a triple of:
        amino acid sequence (string)
        start:stop positions in nucleotide sequence (0-indexed int)
        Is the codon full or truncated? (bool)
    """
    sequence = sequence.upper()
    frame  = {}
    frame[0] = [sequence[i:(i+3)] for i in range(0,len(sequence),3)]
    frame[1] = [sequence[i:(i+3)] for i in range(1,len(sequence),3)]
    frame[2] = [sequence[i:(i+3)] for i in range(2,len(sequence),3)]
    orf = ''
    span = []
    stopless = False
    start_met = re.compile('^.*?(M.*)$')
    for f in [0,1,2]:
        translation = ''.join([translate(i) for i in frame[f]])
        potential_orfs = translation.split('-')
        stopless_list = [False]*(len(potential_orfs)-1)+[True]
        if not allow_truncation:
            potential_orfs = potential_orfs[:-1]
            stopless = stopless_list[:-1]
        for p,s in zip(potential_orfs,stopless_list):
            M_match = start_met.match(p)
            if M_match:
                o = M_match.groups()[0]
                if len(o) > len(orf):
                    orf = o
                    if s:
                        stopless = True
                        span = [i*3+f for i in re.search(
                            o,translation).span()]
                    else:
                        stopless = False
                        span = [i*3+f for i in re.search(
                            o+'-',translation).span()]
    return (orf,span,stopless)
