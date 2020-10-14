##### Please load python/3.5.2 for running in terminal
##### A version without using numpy library

##### Step 1 #####
# Ask user for input of a multiple sequence alignment file

import sys
import math #for calculation functions
input_file = sys.argv[1] #read file

seq_align_name = []
seq_align = []
with open(input_file) as fp:
    temp = ''
    for line in fp:
        if line[0] != '\n': 
            seq_align_name.append(line[:line.index("\t")])
            line = line.replace('\n','')
            seq_align.append(line[line.index("\t")+1:])

##### Step 2 #####
# Compute BLOSUM matrix

### 2.1 Observed frequency
a_letters = 0
c_letters = 0
d_letters = 0
e_letters = 0
f_letters = 0
g_letters = 0
h_letters = 0
i_letters = 0
k_letters = 0
l_letters = 0
m_letters = 0
n_letters = 0
p_letters = 0
q_letters = 0
r_letters = 0
s_letters = 0
t_letters = 0
v_letters = 0
w_letters = 0
y_letters = 0
gaps = 0
for i in seq_align: 
    a_letters = a_letters + i.count('A')
    c_letters = c_letters + i.count('C')
    d_letters = d_letters + i.count('D')
    e_letters = e_letters + i.count('E')
    f_letters = f_letters + i.count('F')
    g_letters = g_letters + i.count('G')
    h_letters = h_letters + i.count('H')
    i_letters = i_letters + i.count('I')
    k_letters = k_letters + i.count('K')
    l_letters = l_letters + i.count('L')
    m_letters = m_letters + i.count('M')
    n_letters = n_letters + i.count('N')
    p_letters = p_letters + i.count('P')
    q_letters = q_letters + i.count('Q')
    r_letters = r_letters + i.count('R')
    s_letters = s_letters + i.count('S')
    t_letters = t_letters + i.count('T')
    v_letters = v_letters + i.count('V')
    w_letters = w_letters + i.count('W')
    y_letters = y_letters + i.count('Y')
total_letters = a_letters+c_letters+d_letters+e_letters+f_letters+g_letters+h_letters+i_letters+k_letters+l_letters+m_letters+n_letters+p_letters+q_letters+r_letters+s_letters+t_letters+v_letters+w_letters+y_letters

obsFreq_dict = {'A':a_letters / total_letters,
                'C':c_letters / total_letters,
                'D':d_letters / total_letters,
                'E':e_letters / total_letters,
                'F':f_letters / total_letters,
                'G':g_letters / total_letters,
                'H':h_letters / total_letters,
                'I':i_letters / total_letters,
                'K':k_letters / total_letters,
                'L':l_letters / total_letters,
                'M':m_letters / total_letters,
                'N':n_letters / total_letters,
                'P':p_letters / total_letters,
                'Q':q_letters / total_letters,
                'R':r_letters / total_letters,
                'S':s_letters / total_letters,
                'T':t_letters / total_letters,
                'V':v_letters / total_letters,
                'W':w_letters / total_letters,
                'Y':y_letters / total_letters}

### 2.2 Expected mutations
# If same letter: obsFreq*obsFreq
# If diff letter: 2*x_obsFreq*y_obsFreq
### 2.3 Observed mutations
# a) Pairwise letter alignments --> denominator
# b) count numerator as the sum of the pair
# c) Observed mutations
### 2.4 M(ij) values

# 2.3(a) 
ct = 0
pairLetAlign = 0
while ct < len(seq_align[0]): #for each column in all alignments
    temp = "" #collect the column as a string 
    for each in seq_align: 
        temp = temp + each[ct]
    length = len(temp) - temp.count('-')
    if length >= 2: 
        comb = math.factorial(length) / ( math.factorial(length-2)*math.factorial(2) )  
        pairLetAlign = pairLetAlign + comb
    ct = ct+1

# Other small steps
aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
blosum_res = {}
for i in aa: 
    for j in aa:
        newi = int(aa.index(i))
        newj = int(aa.index(j))
        if newi == newj: #if same letter/AA
            # 2.3(b) 
            count = 0
            numerator = 0
            while count < len(seq_align[0]): #for each column in all alignments
                temp = "" 
                for each in seq_align: 
                    temp = temp + each[count] #collect the column as a string 
                num_i = temp.count(i)
                if num_i >= 2: 
                    numerator = numerator + math.factorial(num_i) / ( math.factorial(num_i-2)*math.factorial(2) )  
                count = count+1
            obsMut = numerator / pairLetAlign
            # 2.2
            expMut = obsFreq_dict[i]*obsFreq_dict[j]
            #2.4
            if (expMut != 0 and obsMut/expMut == 0): 
                mij = "NaN"
            else: 
                mij = round( 2 * math.log(obsMut/expMut,2),3 )
            #Collecting results
            blosum_res.update({i+j: mij})  
        elif newi < newj: #if different letter/AA
            # 2.3(b) 
            count = 0
            numerator = 0
            while count < len(seq_align[0]): #for each column in all alignments
                temp = "" 
                for each in seq_align: 
                    temp = temp + each[count] #collect the column as a string 
                num_i = temp.count(i)
                num_j = temp.count(j)
                numerator = numerator + num_i*num_j
                count = count+1
            obsMut = numerator / pairLetAlign
            #2.2
            expMut = 2*obsFreq_dict[i]*obsFreq_dict[j]
            #2.4
            if obsMut/expMut == 0: 
                mij = "NaN"
            else: 
                mij = round( 2 * math.log(obsMut/expMut,2),3 )
            #Collecting results
            blosum_res.update({i+j: mij}) 

##### Step 3 #####
# print BLOSUM matrix on the terminal

matrix = ["\t"+'A'+"\t"+'C'+"\t"+'D'+"\t"+'E'+"\t"+'F'+"\t"+'G'+"\t"+'H'+"\t"+'I'+"\t"+'K'+"\t"+'L'+"\t"+'M'+"\t"+'N'+"\t"+'P'+"\t"+'Q'+"\t"+'R'+"\t"+'S'+"\t"+'T'+"\t"+'V'+"\t"+'W'+"\t"+'Y']
for i in aa: 
    rows = i
    for j in aa:
        newi = int(aa.index(i))
        newj = int(aa.index(j))
        if newi <= newj: 
            inp = i+j
        else: 
            inp = j+i
        rows = rows + "\t" + str(blosum_res[inp])
    matrix.append(rows)

for i in matrix: 
    print(i)


