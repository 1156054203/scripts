#!/usr/bin/env python3.5
from Bio import SeqIO
from Bio.Blast import NCBIWWW,NCBIXML
import argparse
import os,re,sys,logging,threading,gzip,time
from queue import Queue

def parse_args(): 
    parser = argparse.ArgumentParser(usage='\n     python %(prog)s [options]')
    parser.add_argument('fastq',metavar='fastq',help='input fastq file')
    parser.add_argument('-o',dest='outdir',metavar='outdir',default='.',help="output directory [default:%(default)s]")
    parser.add_argument('-f',dest='outfile',metavar='outfile',default='species_stat.txt',help="output file [default:%(default)s]")
    parser.add_argument('-n',dest='seqnum',metavar='seqnum',type=int,default=1000,help='abstract how many reads in each fastq [default:%(default)s]')
    parser.add_argument('-t',dest='thread',metavar='thread',type=int,default=8,help="thread number [default:%(default)s]")
    parser.add_argument('-e',dest='evalue',metavar='evalue',default="1e-5",help="evalue [default:%(default)s]")
    parser.add_argument('-m',dest='maxlimit',metavar='maxlimit',type=int,default=1,help="Number of hitlist_size to logger.debug [default:%(default)s]")
    return parser.parse_args()

def createlogger(name,filename):
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(filename)
    formatter = logging.Formatter('[%(asctime)s] %(name)s-%(levelname)s-%(message)s',datefmt='%Y-%m-%d %X')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

def main(fastq,seqnum,thread,evalue,maxlimit,outdir,outfile):
    threadlist=[]
    recordlist=[]
    filehand = gzip.open(fastq, "rt") if os.path.splitext(fastq)[1] == ".gz" or os.path.splitext(fastq)[1] == ".gzip" else open(fastq, 'r')
    recorditer = SeqIO.parse(filehand, "fastq")
    recordnum=0
    while recordnum < seqnum:
        try:
            logger.info('Reading %s record %d ...' %(fastq, recordnum))
            record = next(recorditer)
            recordlist.append(record)
            recordnum += 1
        except StopIteration:
            break
    filehand.close()

    # Number of threads
    for i in range(thread):
        t = threading.Thread(target=runblast,args=(i,evalue,maxlimit))
        threadlist.append(t)
        t.start() 
    # put data in queue
    for i in range(len(recordlist)):
        q.put(recordlist[i])
    q.join()

    for i in range(thread):
        q.put(None)

    for t in threadlist:
        t.join()
   
    logger.info("all threads done.")
    
    statfh = open(os.path.join(outdir,outfile), 'w')
    totalnum = sum(speciesDict.values())
    for species in speciesDict.keys():
        statfh.write('%s\t%d\t%.2f%%\n' %(species, speciesDict[species], speciesDict[species]/totalnum*100))
    statfh.close()
    
def runblast(i,evalue,maxlimit):
     while True:
         record = q.get()
         if record is None: break
         logger.info("Thread %d blast for record %s start." % (i,record.id))
         trynum = 0
         while trynum < 10:
             trynum += 1
             try:
                 blasttmp = NCBIWWW.qblast("blastx", "nr", record.format("fasta"), expect=evalue, hitlist_size=maxlimit)
                 if blasttmp: break
             except ValueError as err:
                 blastError = err
                 time.sleep(10)
         if trynum >=10: 
             logger.error("blast for record %s has been tried 10 times, still get error!" %(record.id))
             logger.error(blastError)
             continue
         ## parse blast result
         blastRecord = NCBIXML.read(blasttmp)
         blasttmp.close()
         if len(blastRecord.alignments):
             searchobj = re.search(" \[([\w\.\(\)\\-\[\]]+ +[^\]|\[]+)\]", blastRecord.alignments[0].title)
             if searchobj:
                 species = searchobj.group(1)
         else:
             species = "No Hit"

         logger.info("Thread %d blast for record %s done." %(i,record.id))
         threadlock.acquire()
         speciesDict[species] = speciesDict[species] + 1 if species in speciesDict else 1
         threadlock.release()
         q.task_done()

if __name__=='__main__':
    args=parse_args()
    if args.outdir:
        outdir=args.outdir
        outdir = os.path.abspath(outdir)
        if not os.path.exists(outdir): os.makedirs(outdir)
    if args.outfile:
        outfile=args.outfile
    if args.seqnum:
        seqnum=args.seqnum
    if args.thread:
        thread=args.thread
    if args.evalue:
        evalue=args.evalue
    if args.maxlimit:
        maxlimit=args.maxlimit
    
    logger=createlogger('blast',os.path.join(outdir,'run.log'))
    threadlock = threading.Lock()
    q = Queue(thread)
    speciesDict={}
    main(args.fastq,seqnum,thread,evalue,maxlimit,outdir,outfile)
