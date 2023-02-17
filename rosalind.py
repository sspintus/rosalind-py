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

