
'''
В этом модуле описан класс ReadWrite, методы которого позволяют считывать или 
записывать файлы в разных форматах.

Методы класса:

\\textcolor{magenta}{r\_pot\_3\_claster(self, file\_name)}

Этот метод считывает информацию о потенциале (по сути сам потенциал), реализованную нами 
в виде кластерного разложения конфигурационной энергии с точным учетом трехчастичных взаимодействий.
Такая информация записана в файлах типа x\_par.in, x\_par.out.

На вход метод требует имя файла с потенциалом file\_name.

Данные сохраняет в ряде атрибутов экземпляра:

self.if\_angle\_pot  -  значения 'T' or 'F'

self.n\_par,self.n\_bas,self.n\_sort 

self.Rmin,self.Rmin\_a 

self.r\_cut1 

self.mult\_e 

self.Rmax 

self.Rmax\_a 

self.n\_sp\_fi, self.n\_sp\_a, self.n\_sp\_g, self.n\_sp\_emb 

self.R\_sp\_fi

self.R\_sp\_a

self.R\_sp\_g

self.R\_sp\_emb

self.emb - список, содержащий списки self.n\_sp\_emb значений  функций emb

self.bas - список, содержащий списки self.n\_sp\_a значений  функций bas + значение производной в тчк. self.Rmin\_a 
для каждой bas функции записанное в конце каждого списка.

self.g - словарь, ключи которого есть индексы str(pq) (см. формулу разложения конфигурационной энергии),
а значение для каждого ключа есть список self.n\_sp\_g значений  функции g(pq) + 1 непонятное число записанное в конце списка.

self.fi - список, содержащий self.n\_sp\_fi значений  функции  fi + значение производной в тчк. self.Rmin 
записанное в конце списка.


\\textcolor{magenta}{w\_columns(self, *args, file\_name)}

Данный метод позволяет записывать произвольное количество переданных ему списков или
кортежей данных в столбцы, т.е. число столбцов == числу переданных списков или кортежей данных.

На вход метод требует:

file\_name - имя файла, в который будут записаны столбцы данных (передается только как ключевой аргумент)

*args - произвольное количество позиционных аргументов (списков, кортежей с данными)
 (передается перед ключевым аргументом file\_name).
 
\\textcolor{magenta}{r\_columns(self, file\_name)}

Данный метод позволяет считывать произвольное количество столбцов данных типа float 
из файла file\_name, при этом число столбцов данных вычисляется по числу столбцов в первой строке файла 
file\_name.

На вход метод требует только имя файла с данными.

Считанные данные записываются в атрибут класса self.columns.

self.columns - это список, содержащий списки, каждый из которых есть данные из соответствующего столбца
из файла file\_name.


\\textcolor{magenta}{r\_xyz\_format(self,name\_file)}

Этот метод считывает данные из файла формата ".xyz".

На вход нужно дать методу имя файла name\_file имеющего ".xyz" формат.

Считанные данные хранятся в следующих переменных:

self.n\_at

self.info

self.r\_at

self.type\_at


\\textcolor{magenta}{r\_abinit\_out(self, file\_name)}

Этот метод считывает данные из выходных (.out) файлов Abinit.

Методу на вход подается только имя выходного файла Abinit file\_name.
Вся информация о считанных данных помещается в словарь self.ab\_out, ключи которого
есть считанные данные. Доступны ключи:

'natom', 'acell', 'etotal', 'fcart', 'rprim', 'xangst', 'xred', 'time\_cpu', 'strten'.

Каждому ключу соответствует список значений параметров. 

Все энергии записываются в eV, силы в eV/Angstrom, acell в ангстремах.

Отметим, что по ключу 'strten' доступен список списков по 6 чисел - компонент тензора напряжений (для данного расчета), компоненты заданы в GPa, последовательность компонент в списке следующая: $\sigma_{11}$, $\sigma_{22}$, $\sigma_{33}$, 
$\sigma_{32}$, $\sigma_{31}$, $\sigma_{21}$.


\\textcolor{magenta}{r\_rv\_at(self, file\_at)}

Этот метод считывает данные записанные в формате rv\_at (такие данные может записывать программа wr\_rv\_at.f), 
т.е. этот формат содержит информацию о решетке включающую помимо типов атомов, их масс и др. еще и
координаты атомов и скорости их движения.

На вход метод требует только имя файла file\_at записанного в формате rv\_at.

На выходе метод записывает в экземпляр класса словарь:

self.r\_rv\_at\_dict, ключи которого n\_at, n\_mark\_at, size, a\_lattice3, r\_at, v\_at, mass\_at, i\_sort\_at,  num\_at\_r, mark\_green, mark\_at. 


\\textcolor{magenta}{w\_rv\_at(self, file\_at, rv\_at\_dict)}

Этот метод записывает данные в формате rv\_at (такие данные также может записывать программа wr\_rv\_at.f)

На вход метод требует имя файла file\_at, в который нужно записать данные о сверхрешетке, и словарь rv\_at\_dict,
в котором обязательно должны быть ключи n\_at, n\_mark\_at, size, a\_lattice3, r\_at, v\_at, mass\_at, i\_sort\_at, num\_at\_r, mark\_green, mark\_at, т.е. все данные которые нужно записать в указанном формате.


\\textcolor{magenta}{r\_rf\_at(self, file\_at)}

Этот метод считывает данные записанные в формате rf\_at (такие данные может записывать программа wr\_rf\_at.f), 
т.е. этот формат содержит информацию о решетке включающую помимо типов атомов, их масс и др. еще и
координаты атомов и силы, действующие на атомы.

На вход метод требует только имя файла file\_at записанного в формате rf\_at.

На выходе метод записывает в экземпляр класса словарь:

self.r\_rf\_at\_dict, ключи которого 'n\_at', 'size', 'a\_lattice3', 'r\_at', 'f\_at', 'mass\_at', 'i\_sort\_at', 'de\_pot', 
здесь de\_pot - энергия считываемой конфигурации относительно некоторого начала отчета, например идеальной решетки.


\\textcolor{magenta}{w\_rf\_at(self, file\_at, rf\_at\_dict)}

Этот метод записывает данные в формате rf\_at (такие данные также может записывать программа wr\_rf\_at.f)

На вход метод требует имя файла file\_at, в который нужно записать данные о сверхрешетке, и словарь rf\_at\_dict,
в котором обязательно должны быть ключи 'n\_at', 'size', 'a\_lattice3', 'r\_at', 'f\_at', 'mass\_at', 'i\_sort\_at', 'de\_pot',
 т.е. все данные которые нужно записать в указанном формате, здесь de\_pot - энергия записываемой конфигурации относительно некоторого начала отчета, например идеальной решетки.
 
 
\\textcolor{magenta}{r\_seat(self, file\_seat)}

Этот метод считывает данные из файла, записанные в seat формате (аналог программы wr\_seat.f). 

На вход методу нужно подать имя файла file\_seat с данными в seat формате.

Метод записывает в экземпляр класса словарь self.seat\_dict, ключи словаря:

r\_main  - список списков ([0-2][0-2]) компонентов векторов трансляций суперячейки;

n\_seat - число атомов в ячейке;

seat - список списков ([0 - (n\_seat-1)][0-2]) абсолютных координат атомов в ангстремах;

sort\_seat - список сортов атомов;

mass\_seat - список масс атомов;

a\_lattice3 - список постоянных решетки (редко используется);




\\textcolor{magenta}{w\_seat(self, file\_seat, seat\_dict)}

Этот метод записывает данные в файл в seat формате (аналог программы wr\_seat.f). 

На вход методу нужно подать имя файла file\_seat и словарь seat\_dict, в котором должны быть следующие ключи: 

r\_main, n\_seat, seat, sort\_seat, mass\_seat, a\_lattice3. 

Подробное описание ключевых параметров смотри в методе \\textcolor{blue}{r\_seat} этого класса.



\\textcolor{magenta}{r\_mol\_static\_format(self, name\_file)}

Этот метод считывает данные в формате молекулярной статики (его пишет программа molecul\_static\_m.py)

На вход методу нужно подать имя файла name\_file, который содержит данные в указанном формате.

На выходе, метод создает словарь self.mol\_static\_dict, ключи которого - это имена файлов с конфигурациями, а значения - это словари, со следующими ключами: 'n\_at', 'a\_lattice', 'e\_at', 'e\_llena'.



\\textcolor{magenta}{r\_alfa\_file(self, name\_file)}

Этот метод позволяет считывать .alfa файл, т.е. файл который содержит данные о величинах относительных деформаций, 
которые были использованы для некоторой матрицы деформации при деформировании некоторой кристаллической решетки.

Такие файлы обычно записываются в используемых методах, где применяется матрица деформации, запись осуществляется методом \\textcolor{magenta}{w\_alfa\_file(self, *args, alfa, acell, info, name\_file)} из данного класса. 

.alfa файл нужен при обработке результатов расчетов (т.к. необходимо знать величину деформации для конкретной полученной конфигурации).


На вход методу нужно подать имя .alfa файла name\_file.


На выходе, создается словарь self.alfa\_dict, доступны следующие ключи:

'info' - строка, с информацией о проведенной деформации (например, здесь может содержаться вид матрицы деформации);

'acell' - список трех чисел, постоянных решетки в ангстремах (при нулевой деформации);

'alfa\_list' - список значений alfa относительных деформаций;

'etc\_list' - список списков, дополнительных данных содержащихся в .alfa файле, здесь каждый вложеный список - это данные из столбцов, следующих по порядку за первым столбцом alfa деформаций, соответственно длина списка 'etc\_list' равна количеству столбцов в .alfa файле - 1.




\\textcolor{magenta}{w\_alfa\_file(self, *args, alfa, acell, info, name\_file)}

Метод записывает .alfa файл.

На вход методу нужно подать:

alfa - список значений относительных деформаций;

acell - список трех чисел, постоянных решетки в ангстремах для нулевой деформации;

info - строка, содержащая информацию о проведенной деформации, например вид матрицы деформации, а также информацию о списках из args. Замечание: в info не должно быть целой строки только с числами (иначе файл при считывании методом \\textcolor{magenta}{r\_alfa\_file(self, name\_file)}  считается неправильно);

name\_file - имя alfa файла, в который запишутся данные (желательно в имя включать расширение .alfa).

args - произвольное количество списков или кортежей с дополнительными данными (без относительных деформаций), но длина каждого списка должна быть равна длине списка alfa;


На выходе метод записывает .alfa файл в name\_file, в котором первый столбец с данными - это значения относительных деформаций alfa.






Данный сценарий только импортируемый.

'''


class ReadWrite: 
    def r_pot_3_claster(self, file_name):
        f = open(file_name)
        # считывание общей информации: число точек, аргументы R, Rmin, Rmax и др. 
        self.if_angle_pot = f.readline().split()[0]
        self.n_par,self.n_bas,self.n_sort = list(map(int, f.readline().split()[0:3]))
        self.Rmin,self.Rmin_a = list(map(float, f.readline().split()[0:2]))
        for i in 1,2,3,4:
            v = float(f.readline().split()[0])
            if i==1: self.r_cut1 = v
            elif i==2: self.mult_e = v
            elif i==3: self.Rmax = v
            elif i==4: self.Rmax_a = v
        self.n_sp_fi, self.n_sp_a, self.n_sp_g, self.n_sp_emb = list(map(int, f.readline().split()[0:4]))
        self.R_sp_fi = list(map(float, f.readline().split()))        
        self.R_sp_a = list(map(float, f.readline().split()))
        self.R_sp_g = list(map(float, f.readline().split()))
        self.R_sp_emb = list(map(float, f.readline().split()))
        assert (len(self.R_sp_fi)==self.n_sp_fi and 
                len(self.R_sp_a)==self.n_sp_a and
                len(self.R_sp_g)==self.n_sp_g and
                len(self.R_sp_emb)==self.n_sp_emb), ' Длина одного из списков R не равна числу точек n'
        # считывание значений функций в точках для сплайнов
        # emb functions
        self.emb = [0 for i in range(self.n_bas)]
        for i in range(self.n_bas): 
            self.emb[i]=list(map(float,f.readline().split()))
            assert len(self.emb[i])==self.n_sp_emb, ' Длина списка emb['+str(i)+'] не равна числу точек n_sp_emb'
        # bas functions
        self.bas = [0 for i in range(self.n_bas)]
        for i in range(self.n_bas): 
            self.bas[i]=list(map(float,f.readline().split()))
            assert (len(self.bas[i])-1)==self.n_sp_a, ' Длина списка bas['+str(i)+'] не равна числу точек n_sp_a'
        # g functions
        if self.if_angle_pot=='T':
            self.g = dict()
            for i1 in range(self.n_bas):
                for i2 in range(i1+1):
                    self.g[str(i1+1)+str(i2+1)] = list(map(float,f.readline().split()))
                    assert ((len(self.g[str(i1+1)+str(i2+1)])-1)==self.n_sp_g, 
                    ' Длина списка g[',str(i1+1)+str(i2+1),'] не равна числу точек n_sp_g')
        # fi function
        self.fi=list(map(float,f.readline().split()))
        assert (len(self.fi)-1)==self.n_sp_fi, ' Длина списка fi не равна числу точек n_sp_fi'
        
        f.close()
        
    def w_columns(self, *args, file_name):
        f = open(file_name, 'w')
        imax = 0
        for i in args:
            if len(i)>imax: imax=len(i)
        for i in range(imax):
            for j in range(len(args)):
                try: s = '{0:15.8f}     '.format(args[j][i])
                except ValueError: s = '{0:}     '.format(args[j][i])
                except IndexError: 
                    s=20*' '
                    print('Количество значений в переданных списках разное!!!!')
                f.write(s)
            f.write('\n')
        f.close()
        
        
    def r_columns(self, file_name):
        f = open(file_name)
        ii=0
        for i in f:
            cur = map(float,i.split())
            if ii==0: columns = [[] for j in range(len(i.split()))]
            kk=0
            for j in cur: 
                columns[kk].append(j)
                kk+=1
            ii+=1
        self.columns = columns
        f.close()
        
    
    def r_xyz_format(self,name_file):
        f = open(name_file)
        self.n_at = int(f.readline())
        self.info = f.readline()
        self.r_at = []; self.type_at = []
        for i in range(self.n_at):
            s = f.readline().split()
            assert len(s)==4, 'len(s)!=4'
            self.r_at.append( list(map(float,s[1:])) )
            self.type_at.append(s[0])
        f.close()
        assert self.n_at==len(self.r_at), 'self.n_at!=len(self.r_at)'
        


    def r_abinit_in(self, file_name):
        # создаем словарь из списка параметров abinit
        ab_in = dict()
        ab_in['natom']=[]
        ab_in['ntypat']=[]
        ab_in['typat']=[]
        ab_in['znucl']=[]
     
        f1 = open(file_name)
        for line in f1:
                          
            if 'natom' in line:
                ab_in['natom'].append( int(line.split()[1]) )

            elif 'ntypat' in line:
                ab_in['ntypat'].append(int(line.split()[1]))
                                    
            elif 'typat' in line:
                  line = line.split()
                  typat = [int(i) for i in line[1:]]
                  ab_in['typat'].append(typat)
                 
            elif 'znucl' in line:
                  line = line.split()
                  znucl = [int(i) for i in line[1:]]
                  ab_in['znucl'].append(znucl)
                 
          
        f1.close()
        self.ab_in = ab_in

    def r_abinit_out(self, file_name):
        # создаем словарь из списка параметров abinit
        ab_out = dict()
        ab_out['natom']=[]
        ab_out['acell']=[]
        ab_out['etotal']=[]
        ab_out['fcart']=[]
        ab_out['rprim']=[]
        ab_out['xangst']=[]
        ab_out['xred']=[]
        ab_out['time_cpu'] = []
        ab_out['strten'] = []
        ab_out['typat'] = ''
        
    
        f1 = open(file_name)
        for line in f1:
            if 'END DATASET(S)' in line:
                new_line = True
                while True:
                
                    if new_line==True: line = f1.readline()
                    else: ...
                    
                    if 'natom' in line:
                        ab_out['natom'].append( int(line.split()[1]) )
                        new_line = True
        
                    elif 'typat' in line:
                        if 'ntypat' in line: ab_out['typat'] += ''
                        elif not 'ntypat' in line: 
                            typat_list = line.split()[1:]
                            for i in typat_list: ab_out['typat'] += i+' '
                        new_line = True

                    elif 'acell' in line:
                        line = line.split()
                        acell = [0.5291772108*float(i) for i in line[1:-1]]
                        ab_out['acell'].append(acell)
                        new_line = True
                    
                    elif 'etotal' in line: 
                        ab_out['etotal'].append(27.2113845*float(line.split()[1]))
                        new_line = True
                    
                    elif 'fcart' in line:
                        fcart = []
                        fcur = [(27.2113845/0.5291772108)*float(i) for i in line.split()[1:]]
                        fcart.append(fcur)
                        while True:
                            line = f1.readline()
                            try: 
                                fcur = [(27.2113845/0.5291772108)*float(i) for i in line.split()]
                                if len(fcur)==0: float('a')
                                fcart.append(fcur)
                            except ValueError: 
                                new_line = False
                                ab_out['fcart'].append(fcart)
                                break
                    
                    elif 'xangst' in line:
                        xangst = []
                        xacur = [float(i) for i in line.split()[1:]]
                        xangst.append(xacur)
                        while True:
                            line = f1.readline()
                            try: 
                                xacur = [float(i) for i in line.split()]
                                if len(xacur)==0: float('a')
                                xangst.append(xacur)
                            except ValueError: 
                                new_line = False
                                ab_out['xangst'].append(xangst)
                                break                    

                    elif 'xred' in line:
                        xred = []
                        xrcur = [float(i) for i in line.split()[1:]]
                        xred.append(xrcur)
                        while True:
                            line = f1.readline()
                            try: 
                                xrcur = [float(i) for i in line.split()]
                                if len(xrcur)==0: float('a')
                                xred.append(xrcur)
                            except ValueError: 
                                new_line = False
                                ab_out['xred'].append(xred)
                                break
                                
                    elif 'rprim' in line:
                        rprim = []
                        cur = [float(i) for i in line.split()[1:]]
                        rprim.append(cur)
                        for j in 1,2:
                            line = f1.readline()
                            cur = [float(i) for i in line.split()]
                            rprim.append(cur)
                        ab_out['rprim'].append(rprim)
                        new_line = True
                        
                    elif 'strten' in line:
                        strten = []
                        cur = [float(i) for i in line.split()[1:]]
                        strten.extend(cur)
                        line = f1.readline()
                        cur = [float(i) for i in line.split()]
                        strten.extend(cur)
                        # перевод компонент в ГПа (сейчас они в Ha/Bohr**3)
                        strten = [27.2113845*160.217653/0.5291772108**3*i for i in strten]   # GPa
                        ab_out['strten'].append(strten)    # последовательность компонент: 11, 22, 33, 32, 31, 21
                        new_line = True
                        
                       
                        
                    elif 'Proc.   0 individual time (sec): cpu=' in line:
                        cpu1 = float(line.split()[7])
                        ab_out['time_cpu'].append('1 cpu time = '+str(round(cpu1/3600, 4))+' hours ('+str(round(cpu1/(3600*24), 4))+' days)')
                        new_line = True
                        
                    elif '+Overall time at end (sec) : cpu=' in line:
                        cpu1 = float(line.split()[7])
                        ab_out['time_cpu'].append('Full cpu time = '+str(round(cpu1/3600, 4))+' hours ('+str(round(cpu1/(3600*24), 4))+' days)')
                        break                       
                        
                    else: new_line=True
                    
        if len(ab_out['rprim'])==0: ab_out['rprim']=[[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]]
        if len(ab_out['xangst'])==0 and ab_out['natom']==[1] : ab_out['xangst']=[[[0.,0.,0.]]]
        if len(ab_out['xred'])==0 and ab_out['natom']==[1] : ab_out['xred']=[[[0.,0.,0.]]]
                    
        f1.close()
        self.ab_out = ab_out                            
                    


    def r_rv_at(self, file_at):
        f = open(file_at)
        line = f.readline().split()
        n_at = int(line[0]); n_mark_at = int(line[1])
        line = f.readline().split()
        size = [float(i) for i in line]
        line = f.readline().split()
        a_lattice3 = [float(i) for i in line]
        fline = f.read().split()
        f.close()
        iw=0
        r_at=[]
        for i in range(n_at):
            x=float(fline[iw]); iw+=1
            y=float(fline[iw]); iw+=1
            z=float(fline[iw]); iw+=1
            r_at.append([x,y,z])

        v_at=[]
        for i in range(n_at):
            vx=float(fline[iw]); iw+=1
            vy=float(fline[iw]); iw+=1
            vz=float(fline[iw]); iw+=1
            v_at.append([vx,vy,vz])

        mass_at=[]
        for i in range(n_at):
            mass=float(fline[iw]); iw+=1
            mass_at.append(mass)

        i_sort_at=[]
        for i in range(n_at):
            i_sort=int(fline[iw]); iw+=1
            i_sort_at.append(i_sort)

        num_at_r=[]
        for i in range(n_at):
            num_at=int(fline[iw]); iw+=1
            num_at_r.append(num_at)

        mark_green=[]
        for i in range(n_at):
            mark=fline[iw]; iw+=1
            mark_green.append(mark)

        mark_at=[]
        for i in range(n_at):
            mark_at.append([])
            for j in range(n_mark_at): mark_at[i].append(j)
            
        for j in range(n_mark_at):
            for i in range(n_at):
                mark_at[i][j] = fline[iw]; iw+=1
            
        # Создаем словарь данных
        self.r_rv_at_dict = dict( zip(['n_at', 'n_mark_at', 'size', 'a_lattice3', 'r_at', 'v_at', 'mass_at', 'i_sort_at',
                                    'num_at_r' , 'mark_green', 'mark_at' ], 
                                  [n_at, n_mark_at, size, a_lattice3, r_at, v_at, mass_at, i_sort_at,
                                   num_at_r, mark_green, mark_at]) )                        
                


          
            
    def w_rv_at(self, file_at, rv_at_dict):
        f = open(file_at,'w')
        f.write(str(rv_at_dict['n_at'])+'    '+str(rv_at_dict['n_mark_at'])+' \n')
        f.write(str(rv_at_dict['size'][0])+'   '+str(rv_at_dict['size'][1])+'   '+str(rv_at_dict['size'][2])+' \n')
        f.write(str(rv_at_dict['a_lattice3'][0])+'   '+str(rv_at_dict['a_lattice3'][1])+'   '+str(rv_at_dict['a_lattice3'][2])+' \n')
        for i in range(rv_at_dict['n_at']): f.write(str(rv_at_dict['r_at'][i][0])+ ' '*4+str(rv_at_dict['r_at'][i][1])+ ' '*4+str(rv_at_dict['r_at'][i][2])+ ' '*4)
        f.write('\n')
        for i in range(rv_at_dict['n_at']): f.write(str(rv_at_dict['v_at'][i][0])+ ' '*4+str(rv_at_dict['v_at'][i][1])+ ' '*4+str(rv_at_dict['v_at'][i][2])+ ' '*4)
        f.write('\n')        
        for i in range(rv_at_dict['n_at']): f.write(str(rv_at_dict['mass_at'][i])+'   ')
        f.write('\n') 
        for i in range(rv_at_dict['n_at']): f.write(str(rv_at_dict['i_sort_at'][i])+'  ')
        f.write('\n') 
        for i in range(rv_at_dict['n_at']): f.write(str(rv_at_dict['num_at_r'][i])+'  ')
        f.write('\n') 
        for i in range(rv_at_dict['n_at']): f.write(str(rv_at_dict['mark_green'][i])+' ')
        f.write('\n')
        for i in range(rv_at_dict['n_mark_at']):
            for j in range(rv_at_dict['n_at']): f.write(str(rv_at_dict['mark_at'][j][i])+' ')
            f.write('\n') 
        f.close()            
            
            



    def r_rf_at(self, file_at):
        f = open(file_at)
        line = f.readline().split()
        n_at = int(line[0])
        line = f.readline().split()
        size = [float(i) for i in line[0:3]]; de_pot = float(line[3])
        line = f.readline().split()
        a_lattice3 = [float(i) for i in line]        
        fline = f.read().split()
        f.close()
        
        iw=0
        r_at=[]
        for i in range(n_at):
            x=float(fline[iw]); iw+=1
            y=float(fline[iw]); iw+=1
            z=float(fline[iw]); iw+=1
            r_at.append([x,y,z])        
        f_at=[]
        for i in range(n_at):
            fx=float(fline[iw]); iw+=1
            fy=float(fline[iw]); iw+=1
            fz=float(fline[iw]); iw+=1
            f_at.append([fx,fy,fz])
        mass_at=[]
        for i in range(n_at):
            mass=float(fline[iw]); iw+=1
            mass_at.append(mass)

        i_sort_at=[]
        for i in range(n_at):
            i_sort=int(fline[iw]); iw+=1
            i_sort_at.append(i_sort)            

        # Создаем словарь данных
        self.r_rf_at_dict = dict( zip(['n_at', 'size', 'a_lattice3', 'r_at', 'f_at', 'mass_at', 'i_sort_at', 'de_pot'], 
                                      [n_at, size, a_lattice3, r_at, f_at, mass_at, i_sort_at, de_pot]) ) 

       


    def w_rf_at(self, file_at, rf_at_dict):
        f = open(file_at,'w')
        f.write(str(rf_at_dict['n_at'])+' \n')
        f.write(str(rf_at_dict['size'][0])+'   '+str(rf_at_dict['size'][1])+'   '+str(rf_at_dict['size'][2])+ '      '+str(rf_at_dict['de_pot'])+' \n')
        f.write(str(rf_at_dict['a_lattice3'][0])+'   '+str(rf_at_dict['a_lattice3'][1])+'   '+str(rf_at_dict['a_lattice3'][2])+' \n')
        for i in range(rf_at_dict['n_at']): f.write(str(rf_at_dict['r_at'][i][0])+ ' '*4+str(rf_at_dict['r_at'][i][1])+ ' '*4+str(rf_at_dict['r_at'][i][2])+ ' '*4)
        f.write('\n')
        for i in range(rf_at_dict['n_at']): f.write(str(rf_at_dict['f_at'][i][0])+ ' '*4+str(rf_at_dict['f_at'][i][1])+ ' '*4+str(rf_at_dict['f_at'][i][2])+ ' '*4)
        f.write('\n')
        for i in range(rf_at_dict['n_at']): f.write(str(rf_at_dict['mass_at'][i])+'   ')
        f.write('\n') 
        for i in range(rf_at_dict['n_at']): f.write(str(rf_at_dict['i_sort_at'][i])+'  ')
        f.write('\n') 
        f.close()     



    
    def r_seat(self, file_seat):
        f1 = open(file_seat)
        r_main = []
        for i in range(3):
            r_main.append( [float(j) for j in f1.readline().split()] )
        n_seat = int(f1.readline())
        sort_seat = [int(j) for j in f1.readline().split()]
        coord = [float(j) for j in f1.readline().split()]
        ii = 0; seat=[]
        for i in range(n_seat):
            at = []
            for j in range(3):
                at.append( coord[ii] )
                ii+=1
            seat.append( at )
        mass_seat = [float(j) for j in f1.readline().split()]
        a_lattice3 = [float(j) for j in f1.readline().split()]
        f1.close()
        # Объединяем считанные данные в словарь
        self.seat_dict = dict(r_main=r_main, n_seat=n_seat, sort_seat=sort_seat,
                                 seat=seat, mass_seat=mass_seat, a_lattice3=a_lattice3)
        
        
        
    def w_seat(self, file_seat, seat_dict):
        f1 = open(file_seat, 'w')
        for i in range(3):
            f1.writelines( [str(j)+'   ' for j in seat_dict['r_main'][i]] )
            f1.write('\n')
        f1.write(str(seat_dict['n_seat'])+'\n')
        f1.writelines([str(j)+'  ' for j in seat_dict['sort_seat']])
        f1.write('\n')
        for i in range(seat_dict['n_seat']):
            f1.writelines([str(j)+'  ' for j in seat_dict['seat'][i]])
        f1.write('\n')
        f1.writelines([str(j)+'  ' for j in seat_dict['mass_seat']])
        f1.write('\n')        
        f1.writelines([str(j)+'   ' for j in seat_dict['a_lattice3']])
        f1.write('\n')
        
    
    
    
    def r_mol_static_format(self, name_file):
        mol_static_dict = dict()
        f1 = open(name_file)
        for i in 0,1: f1.readline()
        for line in f1:
            line1 = line.split()
            mol_static_dict[line1[0]] = {}
            mol_static_dict[line1[0]]['n_at'] = int(line1[1])
            mol_static_dict[line1[0]]['a_lattice'] = [float(i) for i in line1[2:5]]
            mol_static_dict[line1[0]]['e_at'] = float(line1[5])
            mol_static_dict[line1[0]]['e_llena'] = mol_static_dict[line1[0]]['e_at']*mol_static_dict[line1[0]]['n_at']
        f1.close()
        self.mol_static_dict = mol_static_dict
            
      
      
            
    def r_alfa_file(self, name_file):
    
        info = ''  # некоторая информация об alfa файле (как правило содержит матрицу деформации)
        etc_list = []; alfa_list = [] # списки, содержащие значения alfa относительных деформаций и проч. данные
        for i in open(name_file):
            line = i.split()
            if len(line)==0: continue
            if 'acell' in i: acell = [float(j) for j in line[1:4]]
            else:
                try:
                    l_cur = [float(j) for j in line]
                    if len(l_cur)>1: etc_list.append( l_cur[1:] )
                    alfa_list.append( l_cur[0] )
                except ValueError: info+=i
        # переделаем self.etc_list
        try: l1 = [[] for i in range(len(etc_list[0]))]
        except IndexError: l1=[]
        for j in etc_list:
            for i in range(len(l1)): l1[i].append(j[i])
        etc_list = l1
        # сохраняем все считанные данные в словарь self.alfa_dict
        self.alfa_dict = dict(zip(['info','etc_list','alfa_list','acell'],[info,etc_list,alfa_list,acell]))
            
            
    
    def w_alfa_file(self, *args, alfa, acell, info, name_file):
        self.w_columns(alfa, *args, file_name=name_file)
        f1 = open(name_file)
        lines = f1.readlines()
        f1.close()
        f1 = open(name_file, 'w')
        f1.write('acell    '+str(acell[0])+'   '+str(acell[1])+'   '+str(acell[2])+'\n\n')
        f1.write(info+'\n\n')
        
        f1.writelines(lines)
        f1.close()
        
        
            
            
            
            
            
            
         
            
            
        
        
        
        
        
       
