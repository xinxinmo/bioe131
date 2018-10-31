#!/usr/bin/python

"""

Needleman-Wunsch Aligner
Bioengineering 131/231, Fall 2018

Command-line script to read in the contents of a multi-FASTA containing two
sequences and a score matrix, perform a global alignment, and print the
resulting alignment matrix and optimal alignment to STDOUT.

"""

import os
import sys

class NWAligner:
    def __init__(self, score_matrix_fname):
        self.score_matrix, self.gap_penalty = self.load_score_matrix(score_matrix_fname)

    @staticmethod
    def load_score_matrix(fname):
        """
        Input: (String) A path to a scoring matrix file.
        Output: (Dictionary) A nested dictionary where keys are strings
                and elements are scores as integers.
    
        Example:
    
        >>> matrix, gap_penalty = NWAligner.load_score_matrix('/home/bioe131/BLOSUM62')
        >>> matrix['A']['A']
        4
        >>> matrix['W']['W']
        11
        >>> gap_penalty
        -4

        """

        score_matrix = {}
        gap_penalty = None
        name = []

        with open(fname) as fp:
            for line_num, line in enumerate(fp):
                # ignore comments in matrix file
                if line.startswith("#"):
                    continue

                ### TODO ###
                # Parse matrix file line-by-line and load into nested dictionaries.
                # Last line of matrix contains the gap penalty which must be pulled
                # out and returned.
                lineList = line.split()
                if len(lineList) > 1:
                    # if the line is not the last line
                    if line_num == 0:
                        # if at the first line
                        name = lineList
                    else:
                        # if at the line with numbers
                        score_matrix.setdefault(name[line_num-1], {})
                        for i in range(len(lineList)):
                            score_matrix[name[line_num-1]][name[i]] = int(lineList[i])
                else:
                    # if the line is the last line
                    gap_penalty = int(line)
        return score_matrix, gap_penalty

    @staticmethod
    def load_FASTA(fname):
        """
        Input: (String) A path to a FASTA file containing exactly two sequences.
        Output: (List) A list containing two strings: one for each sequence.

        Example:

        >>> seqs = NWAligner.load_FASTA('example.fa')
        >>> seqs[0]
        'YAADSKATPGNPAFHQDEIFLARIAFIYQMWDGGQLKLIDYAPHHVMCEE'
        >>> seqs[1]
        'WVGQPNMKVQHWSNMKACCVKFITWTFIAPEKHACKWTETAYQADCDIIW'
        >>> len(seqs)
        2

        """

        seqs = []

        ### TODO ###
        # Load FASTA file and return list of sequences.
        # Throw an error if there are more than two sequences in the file.
        # open the file
        infile = open(fname)
        counter = 0
        for line in infile.readlines():
            if line.startswith(">"):
                counter += 1
                continue
            else:
                if counter > 2:
                    raise ValueError('There are more than two sequences in the file!')
                    break
                seqs.append(line[:-1])    
        infile.close()
        return seqs

    def align(self, seq_x, seq_y, print_matrix = False):
        """
        Input: (Strings) Two sequences to be aligned (seq_x and seq_y).
               (Boolean) If print_matrix is True, print the dynamic programming
                         matrix before traceback.
        Output: (Tuple of strings) Two sequences, aligned.

        Example:

        >>> aligner = NWAligner('BLOSUM62')
        >>> seqs = aligner.load_FASTA('example.fa')
        >>> aligner.align(seqs[0], seqs[1])
        ('YAAD-SKATPGNPAF---HQDEIF--L-AR--IA-FIYQM-WDGGQLK-LIDYAPH-HVM-C---E-------E---',
         'W---VGQ--P-N--MKVQH----WSNMKA-CCV-KFI---TW------TFI--APEKH--ACKWTETAYQADCDIIW')

        """

        ###
        ### INITIALIZATION
        ###

        # create two empty matrices with sizes based on the input sequences.
        # one contains the dynamic programming matrix, the other contains
        # pointers we'll use during traceback
        matrix = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]
        pointers = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]

        ### TODO ###
        # Fill the top row of the matrix with scores for gaps.
        for i in range(1, len(seq_x) + 1):
            if matrix[0]:
                matrix[0][i] = -10 * i
        # Fill the first column of the matrix with scores for gaps.
        for i in range(1, len(seq_y) + 1):
            if matrix[i]:
                matrix[i][0] = -10 * i
        
        ###
        ### RECURSION
        ###

        # fill the dynamic programming and pointer matrices
        '''The dynamic programming matrix has been filled with the max score, each evaluated using the gap penalty or the score 
        matrix. The pointer matrix is filled in with 1s in certain positions based off of which score was determined to be the 
        max. For example if score_diag is the max, the diagonal position[x-1][y-1] in the pointer matrix will be assigned the 
        value 1. This will be our method of traceback - to follow the 1s back to the origin.'''
        ### TODO ###
                # Take the maximum score of three possibilities:
                #   1) The element in the matrix diagonal from this one
                #      plus the score of an exact match
                #   2) The element to the left plus a gap penalty
                #   3) The element above plus a gap penalty
                # ... and set the current element (matrix[x][y]) equal to that
                #
                # Keep track of which of these choices you made by setting
                #   the same element (i.e., pointers[x][y]) to some value that
                #   has meaning to you.
        for x in range(1, len(seq_x) + 1):
            for y in range(1, len(seq_y) + 1):
                match_score = self.score_matrix[seq_x[x - 1]][seq_y[y - 1]]
                score_diag = matrix[x - 1][y - 1] + match_score
                score_left = matrix[x][y - 1] + self.gap_penalty
                score_above = matrix[x - 1][y] + self.gap_penalty
                matrix[x][y] = max(score_diag, score_left, score_above)
                
                if matrix[x][y] == score_diag:
                    pointers[x - 1][y - 1] = 1
                elif matrix[x][y] == score_left:
                    pointers[x][y - 1] = 1
                elif matrix[x][y] == score_above:
                    pointers[x - 1][y] = 1
        pointers[x][y] = 1
        # print the dynamic programming matrix
        #print_matrix = True
        if print_matrix:
            for x in range(len(seq_x) + 1):
                print " ".join(map(lambda i: str(int(i)), matrix[x]))
        ###
        ### TRACEBACK
        ###

        # starting from the bottom right corner, follow the pointers back
        x, y = len(seq_x), len(seq_y)

        # fill these lists with the aligned sequences
        align_x = []
        align_y = []
        
        '''Goes through the pointer matrix which has been filled with 1s according to if that particular position has been 
        chosen for the score of the current position. By matching the values at different positions to that of the  
        current position, we append elements of the sequence into the alignment accordingly.'''
        while x > 0 or y > 0:
            move = pointers[x][y]
            if pointers[x-1][y-1] == move:
                align_x.append(seq_x[x-1])
                align_y.append(seq_y[y-1])
                x -= 1
                y -= 1
            elif pointers[x-1][y] == move:
                align_x.append(seq_x[x-1])
                align_y.append('-')
                x -= 1
            elif pointers[x][y-1] == move:
                align_x.append('-')
                align_y.append(seq_y[y-1])
                y -= 1
        
        # flip the alignments, as they're reversed
        return ("".join(align_x[::-1]), "".join(align_y[::-1]))

###                                      ###
### NO NEED TO EDIT CODE BELOW THIS LINE ###
###                                      ###

if __name__ == '__main__':
    def usage():
        print 'usage: %s matrixfilename stringfilename'
        sys.exit(1)

    if len(sys.argv) != 3:
        usage()

    for fname in sys.argv[1:]:
        if not os.path.isfile(fname):
            print 'Can not open %s' % (fname,)
            usage()

    aligner = NWAligner(sys.argv[1])
    seqs = aligner.load_FASTA(sys.argv[2])
    result = aligner.align(seqs[0], seqs[1])

    print('>seq1\n%s\n>seq2\n%s' % (result[0], result[1]))
