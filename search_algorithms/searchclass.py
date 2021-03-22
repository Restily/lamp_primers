import time
import numpy as np

class Search():
    F3_F2 = [0, 60]
    F2_F1 = [20, 40]
    F1_B1c = [0, 30] 
    st = [False]*5
    f = [1, 2, 0, 1, 2]

    def __init__(self, results_53, results_35):  
        self.results_53 = results_53
        self.results_35 = results_35
        self.lens = [len(results_35), len(results_53)]
        self.search_primers = []

    def get_primers(self):
        return self.search_primers

    def nucls(self, result1, result2):
        """ Проверка, что весь набор лежит в диапазоне 280 нуклеотидов """
        if result2[0] - result2[1][1] - result1[0] <= 280:
            return True
        return False

    def check(self, primer1, result1, result2, dist):
        """ Проверка праймеров на условия (расстояние) """
        nuclsnum = self.nucls(primer1, result2)
        if not nuclsnum:
            return False

        if result2[0] - result2[1][1] - result1[0] > dist[0]:
            if result2[0] - result2[1][1] - result1[0] <= dist[1]:
                if result1[1][0] == result2[1][0] or abs(result1[1][1] - result2[1][1]) > 3: 
                    return 2
                else:
                    return True
            else:
                return False
        
        return 2

    def start_search_set(self):
        """ Запуск создания наборов """
        for f_ind in range(self.lens[0]):
            self.st = [False]*5
            self.set_primers(sc=1, indexes=[f_ind, 0, 0, 0, 0, 0])

    def set_primers(self, sc, indexes):
        """ Основной алгоритм """
        for ind in range(self.f[sc-1], self.lens[sc % 2]):
            if sc % 2 == 0:
                step = self.check(self.results_53[indexes[0]], self.results_53[indexes[sc-1]], self.results_35[ind], self.F2_F1)
            else:
                if sc == 1:
                    step = self.check(self.results_53[indexes[0]], self.results_53[indexes[sc-1]], self.results_53[ind], self.F3_F2)
                elif sc == 3:
                    step = self.check(self.results_53[indexes[0]], self.results_35[indexes[sc-1]], self.results_53[ind], self.F1_B1c)
                else:
                    step = self.check(self.results_53[indexes[0]], self.results_35[indexes[sc-1]], self.results_35[ind], self.F3_F2)

            if step == 1:
                if not self.st[sc-1]:
                    self.f[sc-1] = ind
                    self.st[sc-1] = True
                
                indexes[sc] = ind

                if sc == 5:
                    self.search_primers.append(indexes)
                else:
                    self.set_primers(sc+1, indexes)
            elif not step:
                return None
        
        return None