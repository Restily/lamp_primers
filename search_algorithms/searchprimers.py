def step(primer1, result1, result2, ras):
    """ Проверка праймеров на условия (температуара, длина, расстояние) """
    nuclsnum = nucls(primer1, result2)
    if not nuclsnum:
        return False

    if result2[0] - result2[1][1] - result1[0] > ras[0]:
        if result2[0] - result2[1][1] - result1[0] <= ras[1]:
            if result1[1][0] == result2[1][0] or abs(result1[1][1] - result2[1][1]) > 3: 
                return 'not yet'
            else:
                return True
        else:
            return False
    
    return 'not yet'

def nucls(result1, result2): #проверка, что весь набор лежит в диапазо
    if result2[0] - result2[1][1] - result1[0] <= 280:
        return True
    return False

def Set_Search(results_g53, results_g35): #results_g53 - 5'-3', results_g35 - 3'-5'
    srch_primers = []
    f0 = 0
    f1 = 1
    f2 = 2
    d1 = 0
    d2 = 1
    d3 = 2
    F3_F2 = [0, 60]
    F2_F1 = [20, 40]
    F1_B1c = [0, 30]
    len_g35 = len(results_g35)
    len_g53 = len(results_g53)    

    for s0 in range(f0, len_g53): 
        st1 = False
        st2 = False
        st3 = False
        st4 = False
        st5 = False

        for s1 in range(f1, len_g53): #step 1
            step1 = step(results_g53[s0], results_g53[s0], results_g53[s1], F3_F2)

            if step1 == True and st1 == False:
                f1 = s1
                st1 = True
            
            if step1 == False:
                break

            if step1 != 'not yet':
                for s2 in range(d1, len_g35): #step 2
                    step2 = step(results_g53[s0], results_g53[s1], results_g35[s2], F2_F1)
                    
                    if step2 == True and st2 == False:
                        d1 = s2
                        st2 = True
                    
                    if step2 == False:
                        break
                    
                    if step2 != 'not yet':
                        for s3 in range(f2, len_g53): #step3
                            step3 = step(results_g53[s0], results_g35[s2], results_g53[s3], F1_B1c)
                            
                            if step3 == True and st3 == False:
                                f2 = s3
                                st3 = True
                            
                            if step3 == False:
                                break
                            
                            if step3 != 'not yet':
                                for s4 in range(d2, len_g35): #step 4
                                    step4 = step(results_g53[s0], results_g53[s3], results_g35[s4], F2_F1)
                                    
                                    if step4 == True and st4 == False:
                                        d2 = s4
                                        st4 = True
                                    
                                    if step4 == False:
                                        break
                                    
                                    if step4 != 'not yet':
                                        for s5 in range(d3, len_g35): #step 5
                                            step5 = step(results_g53[s0], results_g35[s4], results_g35[s5], F3_F2)
                                            
                                            if step5 == True and st5 == False:
                                                d3 = s5
                                                st5 = True
                                            
                                            if step5 == False:
                                                break
                                            
                                            if step5 != 'not yet':
                                                srch_primers.append([results_g53[s0], results_g53[s1], results_g35[s2], results_g53[s3], results_g35[s4], results_g35[s5]]) #'F3 ', 'F2 ', 'F1c', 'B1c', 'B2 ', 'B3 '

    return srch_primers


