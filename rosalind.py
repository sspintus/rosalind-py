def sum_odds(start, upto):
    sum=0
    if start%2!=1:  start=start+1
    for i in range(start, upto+1, 2): sum = sum+i
    return sum

def even_lines(file_name):
    text=''
    file=open(file_name, 'r')
    counter=1
    for line in file:
        if counter%2==0:
            text=text+line
        counter+=1
    return text

def str2dict(string):
    words=string.split(' ')
    dwords={}
    for w in words:
        if(w not in dwords):
            dwords[w] = 0
        dwords[w]+=1
    for key, value in dwords.items():
        print(key,value)

def count_nt(dna_seq):
    atgc={}
    for n in dna_seq:
        if(n not in atgc):
            atgc[n]=0
        atgc[n]+=1
    print(atgc['A'],atgc['C'],atgc['G'],atgc['T'])

def dna2rna(dna_seq):
    rna_seq=''
    for n in dna_seq:
        if(n=='T'):
            n='U'
        rna_seq+=n
    return rna_seq

def complement(dna_seq):
    rev_cmpl=''
    for n in dna_seq[::-1]:
        match n:
            case 'A':
                pair='T'
            case 'T':
                pair='A'
            case 'G':
                pair='C'
            case 'C':
                pair='G'
        rev_cmpl+=pair
    return rev_cmpl

def gc(dna_seq):
    if len(dna_seq)==0:
        return 0
    gc=0
    for nt in dna_seq:
        if(nt=='G' or nt=='C'):
            gc+=1;
    gc*=100
    gc/=len(dna_seq)
    return gc

def max_gc(fasta):
    is_end_record=False
    is_eof=False
    seq=''
    max_gc=0
    max_id=''
    f=open(fasta)
    while(True):
        line=next(f,'$')
        line=line.rstrip()
        if(len(line)==0):
            line='$'
        if(line[0]=='>' or line[0]=='$'):
            seq_gc=gc(seq)
            seq=''
            if(seq_gc > max_gc):
                max_gc=seq_gc
                max_id=seq_id
            if(line[0]=='>'):
                seq_id=line[1:]
            else:
                break
        else:
            seq+=line
    f.close()
    return(max_id+'\n'+str(max_gc))

def hamming(reads_filename):
    f=open(reads_filename)
    reads=[f.readline().rstrip(),f.readline().rstrip()]
    f.close
    distance=0
    for i in range(0,len(reads[0])):
        if(reads[0][i]!=reads[1][i]):
            distance+=1
    return distance

def fibonacci_array(months, progeny):
    pairs=[1,1]
    month=1
    while(month < months-1):
        month+=1
        pairs.append(pairs[month-2]*progeny+pairs[month-1])
        print(month+1, pairs[-1])
    return pairs[-1]

def fibonacci_12(months, progeny):
    F1=1
    F2=1
    month=2
    while(month<months):
        month=month+1
        F=F1*progeny+F2
        F1=F2
        F2=F
    return F

def fa2gc(fasta):
    f=open(fasta)
    line=f.readline().rstrip()
    seq_id=line[1:]
    line=f.readline().rstrip()
    seq=line
    max_id=''
    max_gc=0
    while(line):
        while(line and line[0] != '>'):
            seq+=line
            line=f.readline().rstrip()
        seq_gc=100*(seq.count('G')+seq.count('C'))/len(seq)
        seq=''
        if(seq_gc > max_gc):
            max_gc=seq_gc
            max_id=seq_id
        if(line):
            seq_id=line[1:]
            line=f.readline().rstrip()
    f.close
    return max_id,max_gc

def codetable2dic(code_table):
    code_dic={}
    f=open(code_table)
    for code_line in f:
        code_list=code_line.rstrip().split()
        for i in range(0,len(code_list),2):
            code_dic[code_list[i]]=code_list[i+1]
    f.close
    return code_dic

def rna2aa(rna_seq, code_table):
    aa_seq=''
    code_dic=codetable2dic(code_table)
    for i in range(0,len(rna_seq),3):
        codon=rna_seq[i]+rna_seq[i+1]+rna_seq[i+2]
        aa=code_dic[codon]
        if(aa=='Stop'):
            break
        aa_seq+=aa
    return aa_seq

def substring(data_file):
    f=open(data_file)
    target=f.readline().rstrip()
    query= f.readline().rstrip()
    f.close()
    loc=[]
    for shift in range(0,len(target)-len(query)):
        is_match=True
        for i in range(0,len(query)):
            if(query[i] != target[i+shift]):
                is_match=False
                break
        if(is_match):
            loc.insert(len(loc),str(shift+1)) #0-based to 1-based
    return " ".join(loc)
