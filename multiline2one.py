def multi2one(x):
    fr=open('x', 'r')
    fw=open('out.fasta', 'w')
    seq={}
    for line in fr:
        if line.startswith('>'):    #判断字符串是否以‘>开始’
            name=line.split()[0]    #以空格为分隔符，并取序列为0的项。
            seq[name]=''
        else:
            seq[name]+=line.replace('\n', '')
    fr.close()
    for i in seq.keys():
        fw.write(i)
        fw.write('\n')
        fw.write(seq[i])
        fw.write('\n')
    
    fr.close()
