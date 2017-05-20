class Global_Aligner(object):
    """
    Class Global_Aligner is created to build a global alignment
    scoring matrix and make a global alignment of two sequences
    of amino acids with linear or affine gap penalty.
    (If linear one needed do not fill in the last argument)
    """

    def __init__(self, seq1, seq2, substitution_matrix, penalty, ext_penalty=None):
        self.seq1 = seq1
        self.seq2 = seq2
        self.penalty = penalty
        self.ext_penalty = ext_penalty
        self.substitution_matrix = substitution_matrix

    def get_global_matrix(self):
        """building a global alignment scoring matrix"""
        align_matrix = []

        # filling in 2 first fows
        head = []
        head.append('')
        head.append('')
        for acid in range(len(self.seq1)):
            head.append(self.seq1[acid])
        align_matrix.append(head)

        second = [''] * (len(self.seq1) + 2)
        second[1] = 0

        '''filling in rest matrix'''
        if self.ext_penalty is None:

            '''Linear gap penalty case'''
            for i in range(2, len(second)):
                second[i] = second[i - 1] + self.penalty
                # second[i] = 0   # in case of semiglobal alignment; if run comment previous line
            align_matrix.append(second)

            for acid in range(len(self.seq2)):
                new_line = [''] * (len(self.seq1) + 2)
                new_line[0] = self.seq2[acid]
                new_line[1] = align_matrix[-1][1] + self.penalty
                # new_line[1] = 0   # in case of semiglobal alignment; if run comment previous line
                align_matrix.append(new_line)

            for j in range(2, len(self.seq2) + 2):
                for i in range(2, len(self.seq1) + 2):
                    match = align_matrix[j - 1][i - 1] + self.substitution_matrix.get_substitution_value(
                        align_matrix[0][i], align_matrix[j][0])
                    delete = align_matrix[j - 1][i] + self.penalty
                    insert = align_matrix[j][i - 1] + self.penalty
                    align_matrix[j][i] = max(match, delete, insert)

        else:
            '''Affine gap penalty case'''
            second[2] = self.penalty
            for i in range(3, len(second)):
                second[i] = second[i - 1] + self.ext_penalty
                # second[i] = 0   # in case of semiglobal alignment; if run comment previous line
            align_matrix.append(second)

            third = [''] * (len(self.seq1) + 2)
            third[0] = self.seq2[0]
            third[1] = self.penalty
            align_matrix.append(third)

            for acid in range(1, len(self.seq2)):
                new_line = [''] * (len(self.seq1) + 2)
                new_line[0] = self.seq2[acid]
                new_line[1] = align_matrix[-1][1] + self.ext_penalty
                # new_line[1] = 0   # in case of semiglobal alignment; if run comment previous line
                align_matrix.append(new_line)

            gap_check = 0
            for j in range(2, len(self.seq2) + 2):
                for i in range(2, len(self.seq1) + 2):
                    match = align_matrix[j - 1][i - 1] + self.substitution_matrix.get_substitution_value(
                        align_matrix[0][i], align_matrix[j][0])
                    if gap_check == 0:
                        delete = align_matrix[j - 1][i] + self.penalty
                        insert = align_matrix[j][i - 1] + self.penalty
                    else:
                        delete = align_matrix[j - 1][i] + self.ext_penalty
                        insert = align_matrix[j][i - 1] + self.ext_penalty
                    align_matrix[j][i] = max(match, delete, insert)

                    if max(match, delete, insert) is match:
                        gap_check = 0
                    else:
                        gap_check = 1

        return align_matrix

    def get_global_alignment(self):
        """
        global alignment of two sequences. Returns a tuple
        (align-seq1, align-seq2, align-score)
        """
        alignment_matrix = self.get_global_matrix()

        # alignment_matrix_table = AsciiTable(alignment_matrix)
        # print(alignment_matrix_table.table)

        alignment_seq1 = ''
        alignment_seq2 = ''

        i = len(self.seq1) + 1
        j = len(self.seq2) + 1

        alignment_score = alignment_matrix[j][i]

        if self.ext_penalty is None:
            '''Linear gap penalty case'''
            while i > 1 and j > 1:
                val = alignment_matrix[j][i]  # current value
                val_diag = alignment_matrix[j - 1][i - 1]  # value from the left top diagonal
                val_up = alignment_matrix[j - 1][i]  # value from the top
                val_left = alignment_matrix[j][i - 1]  # value from the left

                if val == val_diag + self.substitution_matrix.get_substitution_value(alignment_matrix[0][i],
                                                                                     alignment_matrix[j][0]):
                    alignment_seq1 = alignment_matrix[0][i] + alignment_seq1
                    alignment_seq2 = alignment_matrix[j][0] + alignment_seq2
                    i -= 1
                    j -= 1
                elif val == val_left + self.penalty:
                    alignment_seq1 = alignment_matrix[0][i] + alignment_seq1
                    alignment_seq2 = "-" + alignment_seq2
                    i -= 1
                elif val == val_up + self.penalty:
                    alignment_seq1 = "-" + alignment_seq1
                    alignment_seq2 = alignment_matrix[j][0] + alignment_seq2
                    j -= 1

            while i > 1:
                alignment_seq1 = alignment_matrix[0][i] + alignment_seq1
                alignment_seq2 = "-" + alignment_seq2
                i -= 1

            while j > 1:
                alignment_seq1 = "-" + alignment_seq1
                alignment_seq2 = alignment_matrix[j][0] + alignment_seq2
                j -= 1

        else:
            '''Affine gap penalty case'''
            while i > 1 and j > 1:
                val = alignment_matrix[j][i]
                val_diag = alignment_matrix[j - 1][i - 1]
                val_up = alignment_matrix[j - 1][i]
                val_left = alignment_matrix[j][i - 1]

                if val == val_diag + self.substitution_matrix.get_substitution_value(alignment_matrix[0][i],
                                                                                     alignment_matrix[j][0]):
                    alignment_seq1 = alignment_matrix[0][i] + alignment_seq1
                    alignment_seq2 = alignment_matrix[j][0] + alignment_seq2
                    i -= 1
                    j -= 1
                elif val == val_left + self.penalty or val_left + self.ext_penalty:
                    alignment_seq1 = alignment_matrix[0][i] + alignment_seq1
                    alignment_seq2 = "-" + alignment_seq2
                    i -= 1
                elif val == val_up + self.penalty or val_up + self.ext_penalty:
                    alignment_seq1 = "-" + alignment_seq1
                    alignment_seq2 = alignment_matrix[j][0] + alignment_seq2
                    j -= 1

            while i > 1:
                alignment_seq1 = alignment_matrix[0][i] + alignment_seq1
                alignment_seq2 = "-" + alignment_seq2
                i -= 1

            while j > 1:
                alignment_seq1 = "-" + alignment_seq1
                alignment_seq2 = alignment_matrix[j][0] + alignment_seq2
                j -= 1

        # for tests — uncomment line below
        # print(alignment_seq1 + '\n' + alignment_seq2 + '\n' + str(alignment_score) + '\n')

        return (alignment_seq1, alignment_seq2, alignment_score)


class Local_Aligner(object):
    """
    Class Local_Aligner is created to build a local alignment scoring matrix
    and make a local alignment of two sequences of amino acids
    with linear or affine gap penalty.
    (If linear one needed do not fill in the last argument)
    """

    def __init__(self, seq1, seq2, substitution_matrix, penalty, ext_penalty=None):
        self.seq1 = seq1
        self.seq2 = seq2
        self.penalty = penalty
        self.ext_penalty = ext_penalty
        self.substitution_matrix = substitution_matrix

    def get_local_matrix(self):
        """building a global alignment scoring matrix"""
        align_matrix = []
        for i in range(len(self.seq2) + 2):
            new_line = [''] * (len(self.seq1) + 2)
            align_matrix.append(new_line)

        '''filling in 2 first rows'''
        for i in range(len(self.seq1)):
            align_matrix[0][i + 2] = self.seq1[i]
            align_matrix[1][i + 2] = 0

        '''filling in 2 first columns'''
        for i in range(len(self.seq2)):
            align_matrix[i + 2][0] = self.seq2[i]
            align_matrix[i + 2][1] = 0

        '''filling in rest matrix'''
        align_matrix[1][1] = 0

        if self.ext_penalty is None:
            '''Linear gap penalty case'''
            for j in range(2, len(self.seq2) + 2):
                for i in range(2, len(self.seq1) + 2):
                    match = align_matrix[j - 1][i - 1] + self.substitution_matrix.get_substitution_value(
                        align_matrix[0][i], align_matrix[j][0])
                    delete = align_matrix[j - 1][i] + self.penalty
                    insert = align_matrix[j][i - 1] + self.penalty
                    align_matrix[j][i] = max(match, delete, insert, 0)

        else:
            '''Affine gap penalty case'''
            for j in range(2, len(self.seq2) + 2):
                for i in range(2, len(self.seq1) + 2):
                    match = align_matrix[j - 1][i - 1] + self.substitution_matrix.get_substitution_value(
                        align_matrix[0][i], align_matrix[j][0])
                    delete = align_matrix[j - 1][i] + self.penalty
                    insert = align_matrix[j][i - 1] + self.penalty
                    align_matrix[j][i] = max(match, delete, insert)

        return align_matrix

    def _get_max_element(self):
        """
        returns a tuple (the biggest value; its position)
        of local alignment scoring matrix
        """
        maximum = 0
        maximum_pos = ()
        alignment_matrix = self.get_local_matrix()
        # result_arr = []

        for j in range(2, len(self.seq2) + 2):
            for i in range(2, len(self.seq1) + 2):
                if alignment_matrix[j][i] > maximum:
                    maximum = alignment_matrix[j][i]
                    maximum_pos = (j, i)

        return (maximum, maximum_pos)

    def get_local_alignment(self):
        """
        local alignment of two sequences. Returns a tuple
        (align-seq1, align-seq2, align-score)
        """
        max_element = self._get_max_element()[0]
        alignment_matrix = self.get_local_matrix()
        alignment_seq1 = ''
        alignment_seq2 = ''

        i = self._get_max_element()[1][1]
        j = self._get_max_element()[1][0]

        if self.ext_penalty is None:
            '''Linear gap penalty case'''
            while i > 1 and j > 1 and alignment_matrix[j][i] != 0:
                val = alignment_matrix[j][i]
                val_diag = alignment_matrix[j - 1][i - 1]
                val_up = alignment_matrix[j - 1][i]
                val_left = alignment_matrix[j][i - 1]

                if val == val_diag + self.substitution_matrix.get_substitution_value(alignment_matrix[0][i],
                                                                                     alignment_matrix[j][0]):
                    alignment_seq1 = alignment_matrix[0][i] + alignment_seq1
                    alignment_seq2 = alignment_matrix[j][0] + alignment_seq2
                    i -= 1
                    j -= 1
                elif val == val_left + self.penalty:
                    alignment_seq1 = alignment_matrix[0][i] + alignment_seq1
                    alignment_seq2 = "-" + alignment_seq2
                    i -= 1
                elif val == val_up + self.penalty:
                    alignment_seq1 = "-" + alignment_seq1
                    alignment_seq2 = alignment_matrix[j][0] + alignment_seq2
                    j -= 1

            while i > 1 and alignment_matrix[j][i] != 0:
                alignment_seq1 = alignment_matrix[0][i] + alignment_seq1
                alignment_seq2 = "-" + alignment_seq2
                i -= 1

            while j > 1 and alignment_matrix[j][i] != 0:
                alignment_seq1 = "-" + alignment_seq1
                alignment_seq2 = alignment_matrix[j][0] + alignment_seq2
                j -= 1

        else:
            '''Affine gap penalty case'''
            while i > 1 and j > 1 and alignment_matrix[j][i] != 0:
                val = alignment_matrix[j][i]
                val_diag = alignment_matrix[j - 1][i - 1]
                val_up = alignment_matrix[j - 1][i]
                val_left = alignment_matrix[j][i - 1]

                if val == val_diag + self.substitution_matrix.get_substitution_value(alignment_matrix[0][i],
                                                                                     alignment_matrix[j][0]):
                    alignment_seq1 = alignment_matrix[0][i] + alignment_seq1
                    alignment_seq2 = alignment_matrix[j][0] + alignment_seq2
                    i -= 1
                    j -= 1
                elif val == val_left + self.penalty or val_left + self.ext_penalty:
                    alignment_seq1 = alignment_matrix[0][i] + alignment_seq1
                    alignment_seq2 = "-" + alignment_seq2
                    i -= 1
                elif val == val_up + self.penalty or val_up + self.ext_penalty:
                    alignment_seq1 = "-" + alignment_seq1
                    alignment_seq2 = alignment_matrix[j][0] + alignment_seq2
                    j -= 1

            while i > 1 and alignment_matrix[j][i] != 0:
                alignment_seq1 = alignment_matrix[0][i] + alignment_seq1
                alignment_seq2 = "-" + alignment_seq2
                i -= 1

            while j > 1 and alignment_matrix[j][i] != 0:
                alignment_seq1 = "-" + alignment_seq1
                alignment_seq2 = alignment_matrix[j][0] + alignment_seq2
                j -= 1

        # print(alignment_seq1 + '\n' + alignment_seq2 + '\n' + str(max_element) + '\n')

        return (alignment_seq1, alignment_seq2, max_element)
