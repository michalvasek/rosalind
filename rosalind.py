import sys


#  Title: Counting DNA Nucleotides
#  URL:   http://rosalind.info/problems/dna/
def dna(dataset_file):

    acgt = {'A' : 0,
            'C' : 0,
            'G' : 0,
            'T' : 0 }

    with open(dataset_file, 'r') as file:
        for line in file:
            for char in line:
                if char in acgt:
                    acgt[char] += 1

    result = '%d %d %d %d' % (acgt['A'], acgt['C'], acgt['G'], acgt['T']) 

    return result


if __name__ == '__main__':

    try:
        problem_name = sys.argv[1]
        dataset_file = 'datasets/rosalind_' + problem_name + '.txt'
        output_file = 'output.txt'

        result = locals()[problem_name](dataset_file)
        with open(output_file, 'w') as file:
            file.write(result)

        print('Solved problem name: ' + problem_name)
        print('Input dataset file: ' + dataset_file)
        print('Generated output file: ' + output_file)

    except Exception as e:
        print('An exception has occured:', e)
        sys.exit(1)