import sys

from math import factorial
from itertools import permutations


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


def reverse_complement(dna):

    result = ''
    complement = {'A' : 'T',
                  'T' : 'A',
                  'C' : 'G',
                  'G' : 'C' }

    for char in dna[::-1]:
        result += complement[char]

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


#  Title: Finding a Motif in DNA
#  URL:   http://rosalind.info/problems/subs/
def subs(dataset_file):

    result = []

    with open(dataset_file, 'r') as file:
        s = file.readline()[:-1]
        t = file.readline()[:-1]

    length = len(t)
    for i in range(len(s) - length + 1):
        if s[i:i + length] == t:
            result.append(i + 1)

    result = ' '.join(str(n) for n in result)

    return result


#  Title: Mendel's First Law
#  URL:   http://rosalind.info/problems/iprb/
def iprb(dataset_file):

    result = ''

    with open(dataset_file, 'r') as file:
        for line in file:
            numbers = [int(n) for n in line.split()]

    dom = numbers[0]
    hetero = numbers[1]
    rec = numbers[2]

    individuals = (dom + hetero + rec)
    n_all = individuals * individuals - individuals
    n_100 = individuals * individuals - (hetero + rec) * (hetero + rec) - dom
    n_75  = hetero * hetero - hetero
    n_50  = 2 * hetero * rec

    result = (1.00 * n_100 + 0.75 * n_75 + 0.50 * n_50) / n_all
    result = str(round(result, 5))

    return result


#  Title: Locating Restriction Sites
#  URL:   http://rosalind.info/problems/revp/
def revp(dataset_file):

    result = ''
    dna = ''
    min, max = 4, 12

    with open(dataset_file, 'r') as file:
        for line in file:
            if line[0] == '>':
                continue
            else:
                dna += line[:-1]
    
    len_dna = len(dna)
    for i in range(len_dna):
        for n in range(min, max + 1):
            if i + n > len_dna:
                continue
            segment = dna[i:i + n]
            revc = reverse_complement(segment)
            if segment == revc:
                result += '%d %d\n' % (i + 1, n)

    return result


#  Title: Mortal Fibonacci Rabbits
#  URL:   http://rosalind.info/problems/fibd/
def fibd(dataset_file):

    generation = list()
    children, adults = 0, 0
    generation.append({'children' : children, 'adults' : adults})
    children, adults = 1, 0
    generation.append({'children' : children, 'adults' : adults})

    with open(dataset_file, 'r') as file:
        for line in file:
            numbers = [int(n) for n in line.split()]
            m = numbers[0]
            n = numbers[1]
    
    for month in range(2, m + 1):
        children = generation[month - 1]['adults']
        adults = generation[month - 1]['children']
        adults += generation[month - 1]['adults']
        if month > n:
            adults -= generation[month - n]['children']
        generation.append({'children' : children, 'adults' : adults})

    children = generation[month]['adults']
    adults = generation[month]['children']
    result = str(children + adults)

    return result


#  Title: Enumerating Gene Orders
#  URL:   http://rosalind.info/problems/perm/
def perm(dataset_file):

    result = ''

    with open(dataset_file, 'r') as file:
        n = int(file.readline()[:-1])

    result += str(factorial(n)) + '\n'

    integers = [str(i) for i in range(1, n + 1)]
    perms = permutations(integers)
    for perm in perms:
        result += (' '.join(perm)) + '\n'

    return result


def load_dna(dataset_file):

    dna = ''
    dataset = []

    with open(dataset_file, 'r') as file:
        for line in file:
            if line[0] == '>':
                if dna:
                    dataset.append({'label' : label,
                                    'dna'   : dna})
                    dna = ''
                label = line[1:-1]
            else:
                dna += line[:-1]
        dataset.append({'label' : label,
                        'dna'   : dna})

    return dataset


#  Title: Overlap Graphs
#  URL:   http://rosalind.info/problems/grph/
def grph(dataset_file):

    result = ''
    nodes = load_dna(dataset_file)

    for node in nodes:
        node.update({'prefix' : node['dna'][:3],
                     'suffix' : node['dna'][-3:]})

    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if i == j:
                continue
            if nodes[i]['suffix'] == nodes[j]['prefix']:
                result += '%s %s\n' % (nodes[i]['label'], nodes[j]['label'])

    return result


#  Title: Calculating Protein Mass
#  URL:   http://rosalind.info/problems/prtm/
def prtm(dataset_file):

    result = 0
    weights = {'A' : 71.03711,  'C' : 103.00919, 'D' : 115.02694,
               'E' : 129.04259, 'F' : 147.06841, 'G' : 57.02146,
               'H' : 137.05891, 'I' : 113.08406, 'K' : 128.09496,
               'L' : 113.08406, 'M' : 131.04049, 'N' : 114.04293,
               'P' : 97.05276,  'Q' : 128.05858, 'R' : 156.10111,
               'S' : 87.03203,  'T' : 101.04768, 'V' : 99.06841,
               'W' : 186.07931, 'Y' : 163.06333 }

    with open(dataset_file, 'r') as file:
        protein = file.readline()[:-1]

    for aa in protein:
        result += weights[aa]

    result = str(round(result, 3))
    return result


def load_numbers(dataset_file):

    with open(dataset_file, 'r') as file:
        for line in file:
            numbers = [int(n) for n in line.split()]

    return numbers


#  Title: Calculating Expected Offspring
#  URL:   http://rosalind.info/problems/iev/
def iev(dataset_file):

    result = 0
    couples = [{'genotypes' : 'AA-AA', 'dominant' : 1.00},
               {'genotypes' : 'AA-Aa', 'dominant' : 1.00},
               {'genotypes' : 'AA-aa', 'dominant' : 1.00},
               {'genotypes' : 'Aa-Aa', 'dominant' : 0.75},
               {'genotypes' : 'Aa-aa', 'dominant' : 0.50},
               {'genotypes' : 'aa-aa', 'dominant' : 0.00}]

    counts = load_numbers(dataset_file)
    for i, couple in enumerate(couples):
        couple.update({'count' : counts[i]})

    for couple in couples:
        result += couple['dominant'] * couple['count'] * 2
    result = str(result)

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