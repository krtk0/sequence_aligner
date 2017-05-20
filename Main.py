from Substitution_Matrix import pam120, pam250, blosum62
from Aligners import Global_Aligner, Local_Aligner

# from terminaltables import AsciiTable


def get_sequences(file):
    """
    retreiving secuences from .fasta files to a tuple
    (sequences_catalog_dictionary, sequences_list)
    """
    prot_file = open(file, 'r')
    _protein_sequences = []
    sequence = ''
    number = 0
    sequences_dictionary = {}
    for line in prot_file:
        if line[0] != '>':
            sequence += line
        elif len(sequence) != 0:
            _protein_sequences.append(sequence.replace('\n', ''))
            sequence = ''
            sequences_dictionary[number] = line.replace('\n', '')
            number += 1
    _protein_sequences.append(sequence.replace('\n', ''))
    return (sequences_dictionary, _protein_sequences)

''' --- Get the sequences of amino acids --- '''
protein_sequences = get_sequences('protein-sequences.fasta')
ww_sequences = get_sequences('WW-sequence.fasta')


''' --- Global Aligner test ---'''
global_test = Global_Aligner(ww_sequences[1][0], ww_sequences[1][1], pam250, -8, -3)
# print(AsciiTable(global_test.get_global_matrix()).table)
# print(global_test.get_global_matrix())
print('Global alignment test: ', global_test.get_global_alignment())


''' --- Local Aligner test --- '''
local_test = Local_Aligner(protein_sequences[1][0], protein_sequences[1][3], pam120, -10, -2)
# print(AsciiTable(local_test.get_local_matrix()).table)
# print(local_test.get_local_matrix())
print('Local alignment test: ', local_test.get_local_alignment())