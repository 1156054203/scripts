#!/usr/bin/env python3
#codingL:utf-8
# Author: Haining.Wang
# Common Use Function
import os
import re
import sys
import time
import xlrd
import googletrans

# run cmd, had tested
def run_cmd(input_cmd=None):
    '''
    blat.run_cmd
    :param input_cmd:Cmd Info;
    :return: None;
    '''
    if input_cmd == None:
	    print('Input Cmd == None!')
	    sys.exit()
    print('Run Cmd:\n\t%s;' % input_cmd)
    print( 'Start:%s;' % now_time() );sys.stdout.flush()
    os.system(input_cmd)
    print('Cmd End');
    print( 'End:%s;' % now_time() );sys.stdout.flush()

# get time, had tested
def now_time():
    '''
    common.now_time()
    :return: 20171212_132430
    '''
    n_time = time.strftime('%Y%m%d_%H%M%S',time.localtime())
    return n_time

# argv api
def argv(*argv_list):
    '''
    argv1+argv2(+argv3)(+argv4)
    :param argv_list: input+file/fold+file_name/path+a(append)
    :return: argv_dict: { argv1:{argv2:argv3}, argv1:{argv2:[argv3-1,argv3-2]} }
    '''
    argv_dict = {}
    #print(argv_list);sys.exit()
    for argv_info in argv_list:
        #print(argv_info);sys.exit()
        info_list = re.split('\+',argv_info.strip())
        info_list[0] = info_list[0].lower()
        # check
        if len(info_list) < 2 or len(info_list) > 4:
            print('Error! common.argv: Input argv must name1+name2+name3(+name4)!')
            sys.exit()
        # input+file/fold+...
        if re.search('^input$',info_list[0]):
            info_list[1] = info_list[1].lower()
            if re.search('^file$',info_list[1]):
                if check_info(info_list[2],'file'):
                    pass
            elif re.search('^fold$',info_list[1]):
                pass
        # add suffix: a:append;r: replace;
        if re.search('^(r|a)$',info_list[-1],re.IGNORECASE):
            info_list[-1] = info_list[-1].lower()
        else:
            info_list.append('r')
        if len(info_list) == 3:
            info_list[1] = re.sub('\*',' ',info_list[1])
            if info_list[-1] == 'r':
                argv_dict[info_list[0]] = info_list[1]
            else:
                if info_list[0] not in argv_dict:
                    argv_dict[info_list[0]] = [info_list[1]]
                else:
                    argv_dict[info_list[0]].append(info_list[1])
        else:
            info_list[2] = re.sub('\*',' ',info_list[2]) # -m*2 -> -m 2;
            if info_list[0] not in argv_dict:
                if info_list[3] == 'a':
                    argv_dict.update({ info_list[0]:{info_list[1]:[info_list[2]]} })
                else:
                    argv_dict.update({info_list[0]: {info_list[1]: info_list[2]}})
            elif info_list[0] in argv_dict and info_list[1] not in argv_dict[info_list[0]]:
                if info_list[3] == 'a':
                    argv_dict[info_list[0]].update({ info_list[1]:[info_list[2]] })
                else:
                    argv_dict[info_list[0]].update({ info_list[1]: info_list[2] })
            else:
                if info_list[3] == 'a':
                    argv_dict[info_list[0]][info_list[1]].append(info_list[2])
                else:
                    argv_dict[info_list[0]][info_list[1]] = info_list[2]
    return argv_dict

# argv_dict merge config_dict
def config_set(argv_dict={}):
    '''
    config_set: merge argv_dict and config_set
    :param argv_dict:
    :return: argv_dict
    '''
    if check_info(argv_dict,'dict') == False:
        print('common.config_set: argv_dict Not Dict')
        sys.exit()
    config_dict = {
        # db
        'seqdir' : '/', # gfclient/gfpcr
        'host' : 'localhost', # gfclient/gfpcr/gfserver
        'bit' : '/online/home/wanghn/data/Genome/H_sapiens/Hg19/hg19.2bit',
        'ref' : '/online/home/wanghn/data/Genome/H_sapiens/Hg19/hg19.fa',
        'base': {
            'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
        },
        'psl': {
            'header': [
                'match', 'mismatch', 'rep.match', 'N\'s',
                'Q gap count', 'Q_gap_bases', 'T_gap_count', 'T_gap_bases', 'strand',
                'Q_name', 'Q_size', 'Q_start', 'Q_end', 'T_name', 'T_size', 'T_start', 'T_end',
                'block_count', 'blockSizes', 'qStarts', 'tStarts',
            ],
        },
        'pool': '2',
        'samtools_flag' : {
            0 : 'read paired', 1 : 'read mapped in proper pair',
            2 : 'read unmapped', 3 : 'mate unmapped',
            4 : 'read reverse strand', 5 : 'mate reverse strand',
            6 : 'first in pair', 7 : 'second in pair',
            8 : 'not primary alignment', 9 : 'read fails platform/vendor quality checks',
            10 : 'read is PCR or optical duplicate', 11 : 'supplementary alignment',
        },
        'aa_name' : {
            'Ala' : 'A', 'A' : 'Ala', 'Arg' : 'R', 'R' : 'Arg',
            'Asn' : 'N', 'N' : 'Asn', 'Asp' : 'D', 'D' : 'Asp',
            'Cys' : 'C', 'C' : 'Cys', 'Gln' : 'Q', 'Q' : 'Gln',
            'Glu' : 'E', 'E' : 'Glu', 'Gly' : 'G', 'G' : 'Gly',
            'His' : 'H', 'H' : 'His', 'Ile' : 'I', 'I' : 'Ile',
            'Leu' : 'L', 'L' : 'Leu', 'Lys' : 'K', 'K' : 'Lys',
            'Met' : 'M', 'M' : 'Met', 'Phe' : 'F', 'F' : 'Phe',
            'Pro' : 'P', 'P' : 'Pro', 'Ser' : 'S', 'S' : 'Ser',
            'Thr' : 'T', 'T' : 'Thr', 'Trp' : 'W', 'W' : 'Trp',
            'Tyr' : 'Y', 'Y' : 'Tyr', 'Val' : 'V', 'V' : 'Val',
        },

        # soft
        'gfserver' : {
            'gfserver' : 'gfServer',
            '-tileSize' : '-tileSize=6',
        },
        'gfpcr' : {
            'gfpcr' : 'gfPcr',
            '-minPerfect' : '-minPerfect=10',
            '-minGood' : '-minGood=10',
            '-out' : '-out=psl',
        },
        'gfclient' : {
            'gfclient' : 'gfClient',
            '-minScore' : '-minScore=0',
            '-minIdentity' : '-minIdentity=80',
        },
        'vt' : {
            'vt':'vt',
        },
        'samtools' : {
            'samtools' : 'samtools',
        },
        'vep' : {
            'vep' : 'vep', '--fork' : '--fork 2', '--format' : '--format vcf', '--vcf' : '--vcf',
            '--sift' : '--sift b', '--polyphen' : '--polyphen b',
            '--vcf_info_field' : '--vcf_info_field ANN', '--terms' : '--terms SO',
            '--hgvs' : '--hgvs', '--hgvsg' : '--hgvsg', '--canonical' : '--canonical',
            '--af' : '--af', '--af_1kg' : '--af_1kg', '--af_esp' : '--af_esp',
            '--af_gnomad' : '--af_gnomad', '--offline' : '--offline',
        },
    }
    # merge argv_dict and config_dict
    for test_key in config_dict.keys():
        if test_key not in argv_dict:
            argv_dict[test_key] = config_dict[test_key]
        else:
            if isinstance(config_dict[test_key],(dict)):
                for par_name in config_dict[test_key].keys():
                    if par_name not in argv_dict[test_key]:
                        argv_dict[test_key][par_name] = config_dict[test_key][par_name]
    return argv_dict

# check file/fold/exist
def check_info(input_info=None,check_type=None):
    if input_info == None:
        print('common.check_info: Input_info == None!')
        return ()
    if check_type == None:
        check_type = 'exist'
    if check_type == 'int' or check_type == 'float':
        if isinstance(input_info,(int,float)) == False:
            print('%s Not Int/Float!' % str(input_info))
            return False
    else:
        #print(check_type);sys.exit()
        if re.search('^exist$',check_type,re.IGNORECASE):
            if os.path.exists(os.path.abspath(input_info)) == None:
                print('%s Not Exist!' % input_info)
                return False
        elif re.search('^file$',check_type,re.IGNORECASE):
            if os.path.isfile(os.path.abspath(input_info)) == None:
                print('%s Not File!' % input_info)
                return False
        elif re.search('^fold$',check_type,re.IGNORECASE):
            if os.path.isfile(os.path.abspath(input_info)) == None:
                print('%s Not Fold!' % input_info)
                return False
        elif re.search('^str$',check_type,re.IGNORECASE):
            if isinstance(input_info,(str)) == False:
                print('%s Not Str!' % input_info)
                return False
        elif re.search('^dict$',check_type,re.IGNORECASE):
            if isinstance(input_info,(dict)) == False:
                print('Not Dict!');print(input_info);
                return False
        elif re.search('^list$',check_type,re.IGNORECASE):
            if isinstance(input_info,(list)) == False:
                print('Not List!'); print(input_info)
                return False
        else:
            print('Not Find Check Type:%s!' % check_type)
            sys.exit()

    return True

# merge regions, had test
def merge_region(*region_infos):
    '''
    merge_region: merge regions
    :param region_infos: [ [chr1,pos1,pos2],[chr2,pos1,pos2],...... ]
    :return: region_list: [ [],[], ]
    '''
    if region_infos == None:
        print('Input_infos must list[ [chr1,pos1,pos2],[chr2,pos1,pos2] ]!')
        print(region_infos);
        sys.exit()
    # sort chrom,pos1,pos2
    region_infos = sorted( region_infos,key=lambda x:(x[0],int(x[1]),int(x[2])) )
    # start merge
    merge_regions = []
    for info_list in region_infos:
        if merge_regions:
            find_s = 0
            for merge_num,merge_list in enumerate(merge_regions):
                if merge_list[0] == info_list[0]: # chrom eq
                    if ( int(merge_list[1]) > int(info_list[2]) 
                    or int(merge_list[2]) < int(info_list[1]) ) == False:
                        merge_regions[merge_num] = [ merge_list[0],merge_list[1],info_list[2] ]
                        find_s = 1
                        break
            if find_s == 0:
                merge_regions.append(info_list)
        else:
            merge_regions = [info_list]
    
    return merge_regions

# merge xls
def merge_xls(*input_list):
    '''
    merge_xls: merge xls_list into merge_xls.xls
    :param input_list: name1.xls(x) name2.xls(x);
    :return: output_file: merge_xls.xls;
    '''
    head_dict = {}
    output_fold = os.getcwd()
    output_file = os.path.join(output_fold,'merge_xls.xls')
    output_h = open(output_file,'w')
    for input_xls in input_list:
        if os.path.exists(input_xls):
            with xlrd.open_workbook(input_xls) as xls_h:
                for table_data in xls_h.sheets():
                    nrow_num = table_data.nrows; #ncol_num = table_data.ncols
                    for row_num in range(nrow_num):
                        infos = table_data.row_values(row_num)
                        if head_dict:
                            if infos[0].lower() in head_dict:
                                head_info = ':'.join( sorted( list(map( lambda x:x.lower(),infos)) ) )
                                if head_info != ':'.join( sorted(head_dict.keys()) ): # head not same
                                    print('Head Not Same!');print(input_xls);print(head_dict)
                                    sys.exit()
                            else: # output
                                output_h.write( '\t'.join( list(map(lambda x:str(x),infos)) )+'\n' )
                        elif not head_dict: # head infos
                            head_dict = { head_value.lower():head_num for head_num,head_value in enumerate(infos) }
                            output_h.write( '\t'.join(infos)+'\n' )
        else:
            print('Input File-%s Not xls!' % input_xls);
            sys.exit()
    output_h.close()

    return output_file

# google tran eng -> chs
def google_tran(input_text):
    # google
    google_tran = googletrans.client.Translator(service_urls=['translate.google.cn'])
    # tran
    tran_text = ''
    if len(input_text) > 0:
        if len(input_text) > 5000:
            total_line = ''
            for s_line in re.split('(\.|\,)\s+',input_text):
                if len(total_line)+len(s_line) > 4500:
                    tran_text += google_tran.translate(total_line, dest='zh-CN').text
                    total_line = s_line
                else:
                    total_line += ' '+s_line
            tran_text += google_tran.translate(total_line, dest='zh-CN').text
        else:
            tran_text = google_tran.translate(input_text, dest='zh-CN').text
    return tran_text

if __name__ == '__main__':
    if re.match('now_time|now|time',sys.argv[1],re.IGNORECASE):
        print(now_time())
    elif re.match('merge_xls',sys.argv[1],re.IGNORECASE):
        xls_list = []
        if os.path.isdir(sys.argv[2]):
            input_fold = os.path.abspath(sys.argv[2].strip())
            for input_file in os.popen('dir %s' % input_fold):
                input_name = re.split('\s+',input_file.strip())[-1]
                input_file = os.path.join(input_fold,input_name)
                xls_list.append(input_file)
            print( merge_xls(*xls_list) )
        else:
            print( merge_xls(*sys.argv[2:]) )
    elif re.match('merge_region',sys.argv[1],re.IGNORECASE):
        test = [ ['chr1','25','97'],['chr2','808','808'],['chr1','25','184'],['chr2','777','808'] ]
        print(merge_region(test))
    else:
        print('Use:')
        print('\tpython *py time/now_time;')
        print('\tpython *py merge_xls name1.xls(x) name2.xls(x) ......;')
    print('End')

