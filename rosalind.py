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


#  Title: Transcribing DNA into RNA
#  URL:   http://rosalind.info/problems/rna/
def rna(dataset_file):

    result = ''

    with open(dataset_file, 'r') as file:
        for line in file:
            for char in line:
                char = 'U' if char == 'T' else char
                result += char

    return result


#  Title: Complementing a Strand of DNA
#  URL:   http://rosalind.info/problems/revc/
def revc(dataset_file):

    result = ''
    complement = {'A' : 'T',
                  'T' : 'A',
                  'C' : 'G',
                  'G' : 'C' }

    with open(dataset_file, 'r') as file:
        for line in file:
            for char in line[-1::-1]:
                if char in complement:
                    char = complement[char]
                    result += char
    result += '\n'

    return result


#  Title: Rabbits and Recurrence Relations
#  URL:   http://rosalind.info/problems/fib/
def fib(dataset_file):

    adults, children = 1, 0

    with open(dataset_file, 'r') as file:
        for line in file:
            numbers = [int(n) for n in line.split()]
            n = numbers[0]
            k = numbers[1]

    for month in range(1, n):
        adults, children = adults + children, adults * k

    result = str(adults)

    return result


#  Title: Computing GC Content
#  URL:   http://rosalind.info/problems/gc/
def gc(dataset_file):

    dna = ''
    dataset = []
    result = ''

    with open(dataset_file, 'r') as file:
        for line in file:
            if line[0] == '>':
                if dna:
                    dataset.append({'label' : label,
                                    'dna'   : dna,
                                    'gc'    : None})
                    dna = ''
                label = line[1:-1]
            else:
                dna += line[:-1]
        dataset.append({'label' : label,
                        'dna'   : dna,
                        'gc'    : None})

    for item in dataset:
        gc_count = 0
        for char in item['dna']:
            if char in ['G', 'C']:
                gc_count += 1
        item['gc'] = gc_count / len(item['dna'])

    max = 0
    for item in dataset:
        if item['gc'] > max:
            result = item['label'] + '\n' + str(round(item['gc'] * 100, 3))
            max = item['gc']

    return result


#  Title: Counting Point Mutations
#  URL:   http://rosalind.info/problems/hamm/
def hamm(dataset_file):

    with open(dataset_file, 'r') as file:
        s = file.readline()[:-1]
        t = file.readline()[:-1]

    distance = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            distance += 1

    result = str(distance)

    return result


#  Title: Translating RNA into Protein
#  URL:   http://rosalind.info/problems/prot/
def prot(dataset_file):

    result = ''
    stop_codons = ['UAA', 'UAG', 'UGA']
    codons = {'AAA' : 'K', 'AAC' : 'N', 'AAG' : 'K', 'AAU' : 'N',
              'ACA' : 'T', 'ACC' : 'T', 'ACG' : 'T', 'ACU' : 'T',
              'AGA' : 'R', 'AGC' : 'S', 'AGG' : 'R', 'AGU' : 'S',
              'AUA' : 'I', 'AUC' : 'I', 'AUG' : 'M', 'AUU' : 'I',
              'CAA' : 'Q', 'CAC' : 'H', 'CAG' : 'Q', 'CAU' : 'H',
              'CCA' : 'P', 'CCC' : 'P', 'CCG' : 'P', 'CCU' : 'P',
              'CGA' : 'R', 'CGC' : 'R', 'CGG' : 'R', 'CGU' : 'R',
              'CUA' : 'L', 'CUC' : 'L', 'CUG' : 'L', 'CUU' : 'L',
              'GAA' : 'E', 'GAC' : 'D', 'GAG' : 'E', 'GAU' : 'D',
              'GCA' : 'A', 'GCC' : 'A', 'GCG' : 'A', 'GCU' : 'A',
              'GGA' : 'G', 'GGC' : 'G', 'GGG' : 'G', 'GGU' : 'G',
              'GUA' : 'V', 'GUC' : 'V', 'GUG' : 'V', 'GUU' : 'V',
              'UAC' : 'Y', 'UAU' : 'Y', 'UCA' : 'S', 'UCC' : 'S',
              'UCG' : 'S', 'UCU' : 'S', 'UGC' : 'C', 'UGG' : 'W',
              'UGU' : 'C', 'UUA' : 'L', 'UUC' : 'F', 'UUG' : 'L',
              'UUU' : 'F' }

    with open(dataset_file, 'r') as file:
        s = file.readline()[:-1]

    for i in range(s.find('AUG'), len(s), 3):
        codon = s[i:i+3]
        if codon in stop_codons:
            break
        result += codons[codon]
    
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