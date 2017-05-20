import requests


class Substitution_Matrix(object):
    """
    Class Substitution_Matrix is created to parse and store
    Substitution Matrixes (PAM, Blosum) from mbio.ncsu.edu
    and check substitution values of amino acids
    """

    def __init__(self, matrix_name):
        self.matrix_name = matrix_name
        self.url = 'http://www.mbio.ncsu.edu/bioedit/tables/' + self.matrix_name
        self.matrix_full = requests.get(self.url).text

    def _get_substitution_matrix_list(self):
        """
        transforming the matrix of amino acids
        to a nested-list format
        """
        acids_list = self.matrix_full.split('\n')
        while acids_list[0][0] != 'A':
            acids_list.pop(0)
        acids_list.pop()
        for i in range(len(acids_list)):
            acids_list[i] = acids_list[i].split()
        return acids_list

    def _get_substitution_matrix_dict(self):
        """
        dictionary (Acid : serial number in the matrix) serves
        for convenient access to substitution values
        """
        acids_dictionary = {}
        acids_list = self._get_substitution_matrix_list()
        val = 1
        for acid in acids_list:
            acids_dictionary[acid[0]] = val
            val += 1
        return acids_dictionary

    def get_substitution_value(self, acid1, acid2):
        """
        allows to check alignment value for two amino acids
        based on chosen substitution matrix
        """
        acids_list = self._get_substitution_matrix_list()
        acids_dictionary = self._get_substitution_matrix_dict()
        pos1 = acids_dictionary[acid1] - 1
        pos2 = acids_dictionary[acid2]
        result = int(acids_list[pos1][pos2])
        # print('\n' + acid1 + ' <=> ' + acid2 + ': ' + str(result))
        return result

    def __str__(self):
        return self.matrix_full


pam120 = Substitution_Matrix('PAM120')
pam250 = Substitution_Matrix('PAM250')
blosum62 = Substitution_Matrix('BLOSUM62')
