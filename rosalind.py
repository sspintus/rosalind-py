import math

def expect_offspring(gt_freqs):
#    AA-AA
#    AA-Aa
#    AA-aa
#    Aa-Aa
#    Aa-aa
#    aa-aa
    probs_A=[1.0, 1.0, 1.0, 0.75, 0.5, 0]
    return 2*sum([p*f for p,f in zip(probs_A, list(map(float,gt_freqs.split())))])

def overlap_graph_dic_noset(fasta, overlap=3):
    f=open(fasta)
    entries={}
    seq_id=''
    for line in f:
        line=line.rstrip()
        if(line[0]=='>'):
            seq_id=line.lstrip('>')
            entries[seq_id]=''
        else:
            entries[seq_id]+=line
    f.close()
    for seq_id1 in entries.keys():
        for seq_id2 in entries.keys():
            if(seq_id1 == seq_id2):
                continue
            if(entries[seq_id1][-3:] == entries[seq_id2][:3]):
                print(seq_id1+' '+seq_id2)

def overlap_graph_dic(fasta, overlap=3):
    f=open(fasta)
    entries={}
    seq_id=''
    for line in f:
        line=line.rstrip()
        if(line[0]=='>'):
            seq_id=line.lstrip('>')
            entries[seq_id]=''
        else:
            entries[seq_id]+=line
    f.close()
    edges=set()
    for seq_id1 in entries.keys():
        for seq_id2 in entries.keys():
            if(seq_id1 == seq_id2):
                continue
            if(entries[seq_id1][-3:] == entries[seq_id2][:3]):
                edges.add(seq_id1+' '+seq_id2)
            if(entries[seq_id2][-3:] == entries[seq_id1][:3]):
                edges.add(seq_id2+' '+seq_id1)
    for e in edges:
        print(e)

def fibonacci_array(months, progeny):
    pairs=[1,1]
    month=1
    while(month < months-1):
        month+=1
        pairs.append(pairs[month-2]*progeny+pairs[month-1])
        print(month+1, pairs[-1])
    return pairs[-1]

def fibonacci_12(months, progeny):
    if(months==1 or months==2):
        return 1
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

def mendel(dom,het,rec):
    return 1-(((rec+het/2)**2-rec-(het/4))/((dom+het+rec)*(dom+het+rec-1)))

def consensus(data_file):
    consensus=''
    profile={'A':[],'C':[],'G':[],'T':[]}
    f=open(data_file)
    seq_id=f.readline().rstrip()
    line=f.readline().rstrip()
    seq=''
    is_first_entry=True
    seq_index=0
    while(line):
        while(line and line[0] != '>'):
            seq+=line
            line=f.readline().rstrip()
        if(is_first_entry):
            is_first_entry=False
            for nt in 'ACGT':
                profile[nt]=[0]*len(seq)
        for i in range(0,len(seq)):
            profile[seq[i]][i]+=1
        seq=''
        line=f.readline().rstrip()
    f.close()
    consensus=''

    for i in range(0,len(profile['A'])):
        max_nt=''
        max_score=-1
        for nt in 'ACGT':
            if(profile[nt][i] > max_score):
                max_score = profile[nt][i]
                max_nt=nt
        consensus+=max_nt
    profile_txt=''

    for nt in 'ACGT':
        profile_txt=profile_txt+nt+': '+' '.join(list(map(str,profile[nt])))+'\n'
    return consensus+'\n'+profile_txt

def fibonacci_mortal(month, lifespan):
    newborn=[1,0]
    mature=[0,1]
    fib=[1,1]
    for m in range(2, month):
        newborn.append(mature[m-1])
        mature.append(mature[m-1]+newborn[m-1])
        if(m > lifespan - 1):
            mature[-1]-=newborn[m-lifespan]
        fib.append(mature[-1]+newborn[-1])
    return fib[-1]

def protein_mass(aa_seq,rmasses='masses.txt'):
    f=open(rmasses)
    rm_dic={}
    for line in f:
        entry=line.split()
        rm_dic[line.split()[0]]=float(line.split()[1])
    f.close()
    pmass=0.0
    for aa in aa_seq:
        pmass+=rm_dic[aa]
    return pmass

def random_strings(nt_seq, gc_seq):
    probs=[]
    pN={}
    for p in list(map(float,gc_seq.split())):
        pN['G']=pN['C']=math.log10(p)-math.log10(2)
        pN['A']=pN['T']=math.log10(1-p)-math.log10(2)
        lp=sum(list(map(lambda x: pN[x],nt_seq)))
        probs.append(str("%.3f" % lp))
    return " ".join(probs)
        
