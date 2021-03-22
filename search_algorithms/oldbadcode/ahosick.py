import time #time.time()
import itertools

import ahocorasick 
from memory_profiler import memory_usage #memory_usage()
from Bio import SeqIO
from primer3 import calcTm 

from search_algorithms.searchprimers import Set_Search


def GetPrimers(aut, len_primer, GC, Tm):
    dict_nucl = 'AGCT'

    for el in itertools.product(dict_nucl, repeat=len_primer):
        primer = ''.join(el)

        primer_Tm = calcTm(primer)        
        gc_count = (el.count('G') + el.count('C'))
        gc_count = gc_count/len_primer*100
        
        if (primer_Tm >= Tm[0] and primer_Tm <= Tm[1]) and (gc_count >= GC[0] and gc_count <= GC[1]):
            aut.add_word(primer, [primer, len_primer, str(gc_count), str(primer_Tm)])
        
    return aut

def GetGenom(genom):
    genom2 = ''

    for i in genom:
        if i == 'A':
            genom2 += 'T'
        elif i == 'G':
            genom2 += 'C'
        elif i == 'C':
            genom2 += 'G'
        elif i == 'T':
            genom2 += 'A'
    
    return genom2

def Primer_Search(genom_file, diapazon, len_primer, GC, Tm):
    results = []

    genom_file = SeqIO.read(genom_file, 'fasta') #считывание генома
    genom = str(genom_file.seq)

    if diapazon[1] != 'all':
        genom = genom[diapazon[0]-1: int(diapazon[1])]
    else:
        genom = genom[diapazon[0]-1:]

    genom_2 = GetGenom(genom)

    A = ahocorasick.Automaton() #создание автомата (дерева)

    for len_p in range(len_primer[0], len_primer[1]+1):
        A = GetPrimers(A, len_p, GC, Tm)

    if A:
        A.make_automaton() #конвертируем дерево в автомат, чтобы включить поиск по Ахо-Корасику
        results_g1 = list(A.iter(genom, ignore_white_space=True)) #список индексов и ключей
        results_g2 = list(A.iter(genom_2, ignore_white_space=True))

        primers = Set_Search(results_g1, results_g2)
        
        for nabor in primers:
            results.append(nabor)

        if len(results) >= 10:
            return results
    
    return results
#logs = open(r'text_files/logs.txt', 'w') #создание файла логов

#время алгоритма пересений ничтожно мало, поэтому им можно пренебречь
#logs.write('Длина генома: ' + str(len(genom)) + '\n') #запись логов
#logs.write('Время работы: ' + str(round(time.time()-start_time, 6)) + ' seconds' + '\n')
#logs.write('Время работы алгоритма Ахо-Корасик: ' + str(round(aho_time, 6)) + ' seconds' + '\n') 
#logs.write('Память программы: ' + str(memory_usage())[1:-5] + ' Megabytes')  

#logs.close() #закрытие файла
