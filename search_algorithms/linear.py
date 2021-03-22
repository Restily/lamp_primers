import time #time.time()
import itertools

from Bio import SeqIO
from primer3 import calcTm 

from search_algorithms.searchprimers import Set_Search


def GetPrimers(n_ind, genom, len_primer, GC, Tm, gnum):
    """ Нахождение праймеров в геноме"""
    results = []
    len_g = len(genom)

    for i in range(n_ind, len_g-len_primer): #нахождение праймеров
        if gnum == 1:
            primer = genom[i:i+len_primer]
        else:
            primer = genom[i:i+len_primer][::-1]

        primer_Tm = calcTm(primer)      
        gc_count = (primer.count('G') + primer.count('C'))
        gc_count = gc_count/len_primer*100
        
        if (primer_Tm >= Tm[0] and primer_Tm <= Tm[1]) and (gc_count >= GC[0] and gc_count <= GC[1]):
            results.append([i+len_primer-1, [primer, len_primer, str(gc_count), str(primer_Tm)]])
    
    return results

def GetGenom(genom):
    """ Нахождение 5'-3' цепочки """
    comp = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
    
    return ''.join([comp[nucl] if nucl in comp else nucl for nucl in genom])

def Primer_Search(genom_file, diapazon, len_primer, GC, Tm):
    """ На вход подается 5'-3' цепочка, производится поиск праймеров и их сборка по наборам"""
    results_n1 = []
    results_n2 = []

    genom_file = SeqIO.read(genom_file, 'fasta') #считывание генома
    genom = str(genom_file.seq)

    n_ind = diapazon[0] - 1
    if diapazon[1] != 'all': #преобразуем геном в необходимую длину
        genom = genom[:int(diapazon[1])]

    genom_2 = GetGenom(genom) #цепочка 5'-3'

    for len_p in range(len_primer[0], len_primer[1]+1): #нахождение праймеров
        results_n1 += GetPrimers(n_ind, genom, len_p, GC, Tm, gnum=1)
        results_n2 += GetPrimers(n_ind, genom_2, len_p, GC, Tm, gnum=2)

    #сортировка результатов поиска
    results_n1.sort(key=lambda index: index[0])
    results_n2.sort(key=lambda index: index[0])

    #нахождение наборов
    primers = Set_Search(results_n1, results_n2) 

    return primers

#logs = open(r'text_files/logs.txt', 'w') #создание файла логов

#время алгоритма пересений ничтожно мало, поэтому им можно пренебречь
#logs.write('Длина генома: ' + str(len(genom)) + '\n') #запись логов
#logs.write('Время работы: ' + str(round(time.time()-start_time, 6)) + ' seconds' + '\n')
#logs.write('Время работы алгоритма Ахо-Корасик: ' + str(round(aho_time, 6)) + ' seconds' + '\n') 
#logs.write('Память программы: ' + str(memory_usage())[1:-5] + ' Megabytes')  

#logs.close() #закрытие файла

#Primer_Search('text_files/seq.fasta', [1, 1000], [22, 22], [50, 60], [55, 65])
