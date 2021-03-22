import time

import PySimpleGUI as gui

from search_algorithms.linear import Primer_Search


layout = [
    [gui.Text('Файл генома: ', background_color='gray'), gui.InputText(), #интефейс приложения
                         gui.FileBrowse(button_text='Добавить файл')],
    [gui.Text('Диапазон генома: ', background_color='gray'), gui.InputText(size=(7,2), default_text='1'), gui.Text(' : ', background_color='gray'), gui.InputText(size=(7,2), default_text='all'),
    gui.Text('\t       ' + 'Длина праймеров:          ', background_color='gray'), gui.InputText(size=(7,2), default_text='20'), gui.Text(' : ', background_color='gray'), gui.InputText(size=(7,2), default_text='22')], 
    [gui.Text('GC-состав: '  + '          ', background_color='gray'), gui.InputText(size=(7,2), default_text='40'), gui.Text(' : ', background_color='gray'), gui.InputText(size=(7,2), default_text='60'), gui.Text('%', background_color='gray'),
    gui.Text('\t' + 'Температура плавления: ', background_color='gray'), gui.InputText(size=(7,2), default_text='55'), gui.Text(' : ', background_color='gray'), gui.InputText(size=(7,2), default_text='65'), gui.Text('°C', background_color='gray')
    ],
    [gui.Output(size=(150, 30), key='-OUT-')],
    [gui.Submit(button_text='Начать поиск'), gui.Exit('Выход из программы')]
]

window = gui.Window('Поиск праймеров', layout, background_color='gray', auto_size_buttons=True) #создание окна

while True:
    event, values = window.read() #чтение окна
    
    if event in (None, 'Cancel'): #выход из цикла при условии ошибки или выхода пользователя
        break

    if event == 'Начать поиск':
        start_time = time.time()
        
        window.FindElement('-OUT-').Update('')  
        print('Поиск запущен. Пожалуйста, подождите его завершения, либо нажмите Escape, чтобы остановить выполенение программы.')
        
        genom_file = values[0] #считывание генома
        diapazon = [int(values[1]), values[2]]
        len_primer = [int(values[3]), int(values[4])]
        GC = [float(values[5]), float(values[6])]
        Tm = [float(values[7]), float(values[8])]

        results = Primer_Search(genom_file, diapazon, len_primer, GC, Tm) #вызов основной функции
        
        window.FindElement('-OUT-').Update('')
        
        print('Праймер | Индексы вхождений | GC-состав | Температура плаления \n') #вывод наборов
        if results:
            names = ['F3 ', 'F2 ', 'F1c', 'B1c', 'B2 ', 'B3 ']
            for ind, nabor in enumerate(results):
                if ind == 11:
                    break
                print('Набор #' + str(ind+1))
                for ind_p, primer in enumerate(nabor):
                    print(names[ind_p] + '  ' + primer[1][0] + '  ' + str(primer[0]-primer[1][1]+2) + ' -> ' + str(primer[0]+1) + '  ' + primer[1][2][:5] + '%  ' + primer[1][3][:5])
                print('\nFIP  ' + nabor[2][1][0] + '  ' + nabor[1][1][0])
                print('BIP  ' + nabor[3][1][0] + '  ' + nabor[4][1][0])
                print('\n')
        else:
            print('Не найдено наборов праймеров')
        
        print('Время выполнения программы: ' + str(time.time()-start_time))

window.close() #завершение программы