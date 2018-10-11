#!/usr/bin/env python3
#coding:utf-8
import os
import re
import copy
import sys
import gzip

if os.path.exists('/path'):
    sys.path.append('/path')
else:
    sys.path.append('.')
import common
import vep

# 染色体规范化
def chr_format(chrom,genomic_version='hg19'):
    '''
    vcf.chr_format-染色体规范化: chr23/MT -> chrX/M
    :param argv_infos: chrom,genomic_version
    :return: chr1/X/Y/MT
    '''
    #1 check
    if not chrom :
        print('chrom Not Set!');
        sys.exit()
    #2 format
    chrom = re.split('chr',chrom,re.IGNORECASE)[-1]
    if re.search('^(hg19|grch)$',genomic_version,re.IGNORECASE):
        if re.search('[^a-zA-Z]',chrom) or re.search('^(MT|X|Y)$',chrom,re.IGNORECASE):
            chrom = 'chr'+chrom
            if re.search('^(chrMT)$',chrom,re.IGNORECASE):
                chrom = 'chrM'
            elif re.search('^(chr23|chr24)$',chrom,re.IGNORECASE):
                chrom = 'chrX' if chrom == 'chr23' else 'chrY'
    else:
        print('Not Set Genome-%s!' % genomic_version);

    return chrom

# vcf多位点分割
def decompose(*argv_infos):
    '''
    vcf.vcf_decompose-多位点分割: OLD_MULTIALLELIC=1:3759889:TA/TAA/TAAA/T;
    :param argv_infos: input+file+*.vcf
    :return: *.vcf
    '''
    #1 set argv and check input/output file
    argv_dict = common.argv(*argv_infos)
    if 'input' not in argv_dict or 'file' not in argv_dict['input']:
        print('Input File Not Set!'); print(argv_infos)
        sys.exit()
    argv_dict['input']['file'] = os.path.abspath(argv_dict['input']['file'])
    input_h = open(argv_dict['input']['file'],'r')
    if re.search('\.gz$',argv_dict['input']['file'],re.IGNORECASE):
        input_h = gzip.open(argv_dict['input']['file'],'rt')
    if 'output' not in argv_dict or 'file' not in argv_dict['output']:
        input_name = os.path.basename(argv_dict['input']['file'])
        input_fold = os.path.dirname(argv_dict['input']['file'])
        output_name = re.sub( '\.vcf'+re.split('\.vcf',input_name)[-1],'_decompose.vcf.gz',input_name )
        argv_dict.update({ 'output':{ 'file':os.path.join(input_fold,output_name) } })
    output_h = open(argv_dict['output']['file'],'w')
    if re.search('\.gz$',argv_dict['output']['file']):
        output_h = gzip.open(argv_dict['output']['file'],'wt')
    #2 read
    for line in input_h:
        if re.match('#',line) == None:
            infos = re.split('\t',line.strip())
            allele_list = re.split(',',infos[4])
            if len(allele_list) == 1:
                output_h.write('\t'.join(infos)+'\n')
            else:
                de_alleles = copy.deepcopy(allele_list)
                de_alleles.insert(0,infos[3])
                for allele_info in allele_list:
                    new_infos = copy.deepcopy(infos)
                    decompose_str = 'OLD_MULTIALLELIC=%s:%s:%s' % \
                                    (infos[0],infos[1],'/'.join(de_alleles))
                    new_infos[7] += ';'+decompose_str+';'
                    new_infos[4] = allele_info
                    output_h.write('\t'.join(new_infos)+'\n')
        else:
            output_h.write(line)
    input_h.close()
    output_h.close()

    return argv_dict['output']['file']

# vcf位置矫正
def normalize(*argv_infos):
    '''
    vcf.vcf_normal-vcf位置矫正: pos normal
    :param argv_infos: input+file+*.vcf vt+-o+-o*...vcf vt+-r+-r*...fa
    :return: output_file: *_normal.vcf;
    '''
    #1 set argv
    argv_dict = common.argv(*argv_infos)
    argv_dict = common.config_set(argv_dict)
    #2 input/output
    if 'input' not in argv_dict or 'file' not in argv_dict['input']:
        print('Input File Not Set!'); print(argv_infos)
        sys.exit()
    argv_dict['input']['file'] = os.path.abspath(argv_dict['input']['file'])
    if 'output' not in argv_dict or 'file' not in argv_dict['output']:
        input_name = os.path.basename(argv_dict['input']['file'])
        input_fold = os.path.dirname(argv_dict['input']['file'])
        output_name = re.sub( '\.vcf'+re.split('\.vcf',input_name)[-1],'_normalize.vcf',input_name )
        argv_dict.update({ 'output':{'file':os.path.join(input_fold,output_name)} })
    #3 ref
    if 'vt' not in argv_dict or '-r' not in argv_dict['vt']:
        if 'vt' not in argv_dict:
            argv_dict.update({ 'vt':{'-r':'-r %s' % argv_dict['ref']} })
        else:
            argv_dict['vt'].update({'-r':'-r %s' % argv_dict['ref']})
    if '-o' not in argv_dict['vt']:
        argv_dict['vt'].update({ '-o':'-o %s' % argv_dict['output']['file'] })
    #4 vt normal
    par_str = ' '.join([
        argv_dict['vt'][tag_key] for tag_key in argv_dict['vt'].keys() if tag_key != 'vt'
    ])
    normal_cmd = ' '.join([
        argv_dict['vt']['vt'], 'normalize', par_str,
        argv_dict['input']['file'],
        '>',re.split('\s+',argv_dict['vt']['-o'])[-1]+'.log',
    ])
    common.run_cmd(normal_cmd)

    return argv_dict['output']['file']

# xls header format
def _head_format(header_line=None):
    '''
    vcf.head_format: xls header format
    :param header_line: chr\tpos\t.....
    :return: format_list,format_dict: ['chr','pos',....]
    '''
    if header_line == None:
        print('Xls header Line Not Find!'); print(header_line)
        sys.exit()
    header_list = re.split( '\t',header_line.strup('\n') )
    format_dict = []
    for info_num,info_value in enumerate(header_list):
        info_value = info_value.strip().lower()
        # tumor/normal_id
        if re.search('^(sample_id|tumor_id)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'sample_id'; format_dict['sample_id'] = info_num
        elif re.search('^(normal_id)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'normal_id'; format_dict['normal_id'] = info_num
        elif re.search('^(gene|symbol)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'symbol'; format_dict['symbol'] = info_num
        elif re.search('^(gene_id)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'gene_id'; format_dict['gene_id'] = info_num
        elif re.search('^(mutation_type|consequence)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'consequence'; format_dict['consequence'] = info_num
        elif re.search('^(Variant_Type)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'variant_type'; format_dict['variant_type'] = info_num
        elif re.search('^(chr|chrom)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'chr'; format_dict['chr'] = info_num
        elif re.search('^(pos|start)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'pos'; format_dict['pos'] = info_num
        elif re.search('^(ref|reference)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'ref'; format_dict['ref'] = info_num
        elif re.search('^(alt|alternate)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'alt'; format_dict['alt'] = info_num
        elif re.search('^(hgvs)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'hgvs'; format_dict['hgvs'] = info_num;
        elif re.search('^(hgvsc|hgvs\.c)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'hgvsc'; format_dict['hgvsc'] = info_num;
        elif re.search('^(hgvsg|hgvs\.g)$',info_value,re.IGNORECASE):
            format_dict[info_num] = 'hgvsg'; format_dict['hgvsg'] = info_num;
        else:
            format_dict[info_num] = info_value

    return format_dict

# 将注释后的vcf转化为xls
def vcf_to_xls(*argv_infos):
    # 1 set
    argv_dict = common.argv(*argv_infos)
    # 2 check input
    if 'input' not in argv_dict or 'file' not in argv_dict['input']:
        print('Input File Not Set!'); print(argv_infos)
        sys.exit()
    argv_dict['input']['file'] = os.path.abspath(argv_dict['input']['file'])
    input_h = open(argv_dict['input']['file'],'r',encoding='utf-8')
    if re.search('\.gz$',argv_dict['input']['file'],re.IGNORECASE):
        input_h = gzip.open(argv_dict['input']['file'],'rt',encoding='utf-8')
    if 'output' not in argv_dict or 'file' not in argv_dict['output']:
        input_name = os.path.basename(argv_dict['input']['file'])
        input_path = os.path.dirname(argv_dict['input']['file'])
        output_name = re.sub( '\.vcf'+re.split('\.vcf',input_name)[-1],'_format.xls',input_name )
        argv_dict.update({ 'output':{'file':os.path.join(input_path,output_name)} })
    output_h = open(argv_dict['output']['file'],'w',encoding='utf-8')
    # 3 read
    ann_dict = {}
    format_dict = {}
    output_ann_tag = [
        'Consequence','IMPACT','SYMBOL','Gene',
        'Feature_type','Feature','BIOTYPE','EXON','HGVS',
    ]
    for line in input_h:
        if re.match('#',line): # header
            if re.match('#CHROM',line): # output
                infos = re.split('\t', line.strip())
                output_h.write('\t'.join(['CHROM','POS','REF','ALT']))
                output_h.write('\t'+'\t'.join(output_ann_tag))
                if len(infos) >= 10:
                    for sample_id in infos[9:]:
                        output_h.write('\t%s_Fre\t%s_DP\t%s_AD\t%s_ADF\t%s_ADR' % [sample_id]*5)
                output_h.write('\n')
            elif re.match('##INFO=<ID=ANN,Number=',line): # vep annot header
                line = re.sub('(\s+)|(\">)','',line.strip())
                ann_list = re.split( '\|',re.split(':',line)[-1] )
                ann_dict = { ann_tag.lower():ann_num for ann_num,ann_tag in enumerate(ann_list) }
        else:
            infos = re.split('\t',line.strip())
            if re.search(',',infos[4]):
                print('Please Use Decompose!');
                sys.exit()
            alt = infos[4]; alt_num = 0; total_num = 2; multi_list = []
            output_h.write('\t'.join([ infos[0],infos[1],infos[3],infos[4], ]))
            for info_str in re.split(';',infos[7]):
                if re.match('ANN=',info_str):
                    parse_ann_list = vep.ann_parse(info_str,ann_dict)
                    for ann_tag in output_ann_tag:
                        if len(parse_ann_list) == 0:
                            output_h.write('\tNA')
                        if ann_tag.lower() in ann_dict:
                            for info_num,info_list in enumerate(parse_ann_list):
                                if info_num == 0:
                                    output_h.write('\t%s' % info_list[ann_dict[ann_tag.lower()]])
                                else:
                                    output_h.write(',%s' % info_list[ann_dict[ann_tag.lower()]])
                        else:
                            if ann_tag.lower() == 'hgvs':
                                for info_num,info_list in enumerate(parse_ann_list):
                                    hgvsp = info_list[ann_dict['hgvsp']] \
                                        if len(info_list[ann_dict['hgvsp']]) != 0 else 'NA'
                                    hgvsc = info_list[ann_dict['hgvsc']] \
                                        if len(info_list[ann_dict['hgvsc']]) != 0 else 'NA'
                                    hgvsg = info_list[ann_dict['hgvsg']] \
                                        if len(info_list[ann_dict['hgvsg']]) != 0 else 'NA'
                                    hgvs = '|'.join([ hgvsp,hgvsc,hgvsg ])
                                    if info_num == 0:
                                        output_h.write('\t%s' % hgvs)
                                    else:
                                        output_h.write(',%s' % hgvs)
                            else:
                                print('Ann_Tag-%s Not Ann_Dict!' % ann_tag); print(ann_dict)
                                sys.exit()
                elif re.match('OLD_VARIANT',info_str):
                    alt = re.split('\/',info_str)[-1]
                elif re.match('OLD_MULTIALLELIC',info_str):
                    multi_list = re.split('\/',re.split(':',info_str)[-1])
                    total_num = len(multi_list)
            for multi_num,multi_alt in enumerate(multi_list):
                if multi_alt == alt:
                    alt_num = multi_num-1
            if not format_dict and len(infos) >= 9:
                format_dict = {
                    format_tag : format_num for format_num,format_tag in enumerate(re.split(':',infos[8]))
                }
            if len(infos) >= 10:
                for sample_info in infos[9:]:
                    info_list = re.split(':',sample_info)
                    dp = 0; ad = 0; adf = 0; adr = 0; fre = 0
                    if 'DP' in format_dict:
                        dp = int(info_list[format_dict['DP']])
                    for depth_tag in ['AD','ADF','ADR']:
                        if depth_tag in format_dict:
                            depth_list = re.split(',',info_list[format_dict[depth_tag]])
                            depth = int(depth_list[alt_num])
                            if total_num == len(depth_list):
                                depth = int(depth_list[alt_num+1])
                            if depth_tag == 'AD':
                                ad = depth
                            elif depth_tag == 'ADF':
                                adf = depth
                            else:
                                adr = depth
                    if 'DP4' in format_dict:
                        dp4_list = re.split(',',info_list[format_dict['DP4']])
                        adf = int(dp4_list[2*(alt_num+1)])
                        adr = int(dp4_list[2*(alt_num+1)+1])
                        ad = adf+adr
                    fre = ad/dp if dp != 0 else 0
                    output_h.write('\t%.4f\t%d\t%d\t%d\t%d' % (fre,dp,ad,adf,adr))
            output_h.write('\n')
    output_h.close()
    input_h.close()
    return argv_dict['output']

# xls to maf
def to_maf(*argv_infos):
    '''
    vcf.xls_to_maf
    :param argv_infos: input+file+...
    :return: maf_file
    '''
    # 1 set and check
    argv_dict = common.argv(*argv_infos)
    if 'input' not in argv_dict or 'file' not in argv_dict['input']:
        print('Input File Not Set!'); print(argv_infos)
        sys.exit()
    argv_dict['input']['file'] = os.path.abspath(argv_dict['input']['file'])
    if 'output' not in argv_dict or 'file' not in argv_dict['output']:
        argv_dict['output']['file'] = re.sub( re.split('\.',argv_dict['input']['file'])[-1],'maf',argv_dict['input']['file'] )

    # 2 maf_dict
    maf_list = [
        'Hugo_Symbol', 'Entrez_Gene_Id', 'NCBI_Build',
        'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification',
        'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
        'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode',
        'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2',
        'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2',
        'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2',
        'Verification_Status', 'Validation_Status', 'Mutation_Status',
        'Sequence_Source', 'Validation_Method',
    ]
    maf_header_dict = {
        'Hugo_Symbol'.lower():'symbol', 'Entrez_Gene_Id'.lower():'gene_id', 'NCBI_Build'.lower():'genome',
        'Chromosome':'chr', 'Start_Position':'pos', 'End_Position':'end_pos', 'Strand':'strand',
        'Variant_Classification':'consequence', 'Variant_Type':'variant_type',
        'Reference_Allele':'ref', 'Tumor_Seq_Allele1':'alt', 'Tumor_Seq_Allele2':'alt',
        'Tumor_Sample_Barcode':'sample_id', 'Matched_Norm_Sample_Barcode':'normal_id',
        'Match_Norm_Seq_Allele1':'ref', 'Match_Norm_Seq_Allele2':'ref',
        'Tumor_Validation_Allele1':'alt', 'Tumor_Validation_Allele2':'alt',
        'Match_Norm_Validation_Allele1':'ref', 'Match_Norm_Validation_Allele2':'ref',
    }

    # 3 read and output
    output_h = open(argv_dict['output']['file'],'w')
    output_h.write( '\t'.join(maf_list)+'\n' )
    head_dict = {}
    sample_id = re.split('_|\.',os.path.basename(argv_dict['input']['file']))[0]
    with open(argv_dict['input']['file'],'r') as input_h:
        for line in input_h:
            infos = re.split('\t',line.strip('\n'))
            if  len(infos) >0:
                if not head_dict:
                    head_dict = _head_format(line.strip('\n'))
                else:
                    for tag_name in maf_list:
                        tag_name = tag_name.lower()
                        if tag_name in maf_header_dict:
                            match_name = maf_header_dict[tag_name].lower()
                            if match_name not in head_dict:
                                if match_name == 'genome':
                                    output_h.write('\thg19')
                                elif match_name == 'strand':
                                    output_h.write('\t+')
                                elif match_name == 'end_pos':
                                    if 'ref' in head_dict and 'alt' in head_dict and 'pos' in head_dict:
                                        ref = infos[head_dict['ref']]; alt = infos[head_dict['alt']]
                                        pos = int(infos[head_dict['pos']]); end_pos = pos
                                        if re.search('[^ATCG]',ref,re.IGNORECASE) == None and \
                                        re.search('[^ATGC]',alt,re.IGNORECASE) == None:
                                            if len(ref) > len(alt):
                                                end_pos = pos+abs(len(ref)-len(alt))
                                            output_h.write('\t%d' % end_pos)
                                        else:
                                            output_h.write('\tNA')
                                    else:
                                        print('Error! Ref or Alt or Pos not IN Head dict!'); print(head_dict)
                                        sys.exit()
                                elif match_name == 'variant_type':
                                    ref = infos[head_dict['ref']]; alt = infos[head_dict['alt']]
                                    if re.search('[^ATGC]',ref,re.IGNORECASE) or re.search('[^ATCG]]',alt,re.IGNORECASE):
                                        output_h.write('\tNA')
                                    else:
                                        v_type = 'SNP'
                                        if len(ref) == len(alt):
                                            if len(ref) == 2:
                                                v_type = 'DNP'
                                            elif len(ref) == 3:
                                                v_type = 'TNP'
                                            elif len(re) == 4:
                                                v_type = 'ONP'
                                        else:
                                            v_type = 'INS' if len(ref) < len(alt) else 'DEL'
                                        output_h.write('\t%s' % v_type)
                                elif match_name == 'sample_id':
                                    output_h.write('\t%s' % sample_id)
                                else:
                                    print('Error!%s Not Find in Header!' % match_name);print(head_dict)
                                    sys.exit()
                            else:
                                if match_name == 'symbol':
                                    output_h.write(infos[head_dict[match_name]])
                                else:
                                    output_h.write('\t%s' % infos[head_dict[match_name]])
                        else:
                            output_h.write('\tNone')
                    output_h.write('\n')
    input_h.close()
    output_h.close()

    return argv_dict['output']['file']

if __name__ == '__main__':
    if re.search('to_maf',sys.argv[1],re.IGNORECASE):
        print( to_maf(*sys.argv[2:]) )
    elif re.search('^decompose$',sys.argv[1],re.IGNORECASE):
        print( decompose(*sys.argv[2:]) )
    elif re.search('^normalize$',sys.argv[1],re.IGNORECASE):
        print( normalize(*sys.argv[2:]) )
    elif re.search('^vcf_to_xls$',sys.argv[1],re.IGNORECASE):
        print(vcf_to_xls(*sys.argv[2:]))
    else:
        print('Ex:')
        print('\tpython3 *py to_maf input+file+...;')
        print('\tpython3 *py vcf_to_xls input+file+...;')
    sys.exit()
