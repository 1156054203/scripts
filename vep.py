#!/usr/bin/env python3
import os
import re
import gzip
import sys
if os.path.exists('/path'):
    sys.path.append('/path')
else:
    sys.path.append('.')
import common

# vep注释相关信息
def _vep_info():
    vep_info_dict = {
        'impact'.lower() : {
            'HIGH'.lower() : 4,
            'MODERATE'.lower() : 3,
            'LOW'.lower() : 2,
            'MODIFIER'.lower() : 1,
        }
    }
    return vep_info_dict

# vep注释命令
def annot(*argv_infos):
    '''
    vep.annot -> *_vep.vcf
    :param argv_infos: input+file+*
    :return:
    '''
    print('Vep Run:'); print(common.now_time()); sys.stdout.flush()
    #1 set argv
    argv_dict = common.argv(*argv_infos)
    argv_dict = common.config_set(argv_dict)
    #2 input/output
    if '--input_file' not in argv_dict['vep']:
        if 'input' not in argv_dict or 'file' not in argv_dict['input']:
            print('Input File Not Set!'); print(argv_infos)
            sys.exit()
        elif 'input' in argv_dict and 'file' in argv_dict['input']:
            argv_dict['input']['file'] = os.path.abspath(argv_dict['input']['file'])
            argv_dict['vep'].update({ '--input_file':'--input_file %s' % argv_dict['input']['file'] })
    if '--output_file' not in argv_dict['vep']:
        if 'output' not in argv_dict or 'file' not in argv_dict['output']:
            input_name = os.path.basename(re.split('\s+',argv_dict['vep']['--input_file'])[-1])
            input_fold = os.path.dirname(re.split('\s+',argv_dict['vep']['--input_file'])[-1])
            output_name = re.sub( '\.'+re.split('\.vcf.+',input_name)[-1],'_vep.vcf',input_name )
            argv_dict['vep'].update({ '--output_file':'--output_file %s' % os.path.join(input_fold,output_name) })
        elif 'output' in argv_dict and 'file' in argv_dict['output']:
            argv_dict['output']['file'] = os.path.abspath(argv_dict['output']['file'])
            argv_dict['vep'].update({ '--output_file':'--output_file %s' % argv_dict['output']['file'] })
    #3 run_cmd
    vep_par_str = ' '.join([
        argv_dict['vep'][tag_name] for tag_name in argv_dict['vep'].keys() if tag_name.lower() != 'vep'
    ])
    vep_cmd = ' '.join([
        argv_dict['vep']['vep'],vep_par_str,
        '>',re.split('\s+',argv_dict['vep']['--output_file'])[-1]+'.log','2>&1',
    ])
    common.run_cmd(vep_cmd)
    print('Vep Run End!'); print(common.now_time()); sys.stdout.flush()

    return re.split('\s+',argv_dict['vep']['--output_file'])[-1]

# 解析vep注释并获取最重要注释
def ann_parse(ann_line=None,ann_dict=None):
    if ann_line == None or ann_dict == None:
        print('Ann Info Not Set or Ann Head Not Set!')
        sys.exit()
    argv_dict = common.config_set({})
    parse_ann_list = []
    vep_info_dict = _vep_info()
    ann_list = re.split(',',re.split('=',ann_line.strip())[-1])
    for ann_str in ann_list:
        info_list = re.split('\|',ann_str)
        gene_name = info_list[ann_dict['Gene'.lower()]]
        gene_symbol = info_list[ann_dict['SYMBOL'.lower()]]
        gene_impact = info_list[ann_dict['IMPACT'.lower()]]
        if len(gene_impact) == 0:
            continue
        # hgvsp ENSG.. -> Gene_Name, and AA3-AA1
        if 'hgvsp'.lower() in ann_dict and len(info_list[ann_dict['hgvsp'.lower()]]) > 0:
            info_list[ann_dict['hgvsp'.lower()]] = re.sub(
                re.split(':',info_list[ann_dict['hgvsp'.lower()]])[0], gene_symbol,
                info_list[ann_dict['hgvsp'.lower()]]
            )
            info_list[ann_dict['hgvsp'.lower()]] = re.sub(
                '\%3D','=',info_list[ann_dict['hgvsp'.lower()]]
            )
            for aa_name1 in argv_dict['aa_name'].keys():
                if len(aa_name1) != 1:
                    info_list[ann_dict['hgvsp'.lower()]] = re.sub(
                        aa_name1,argv_dict['aa_name'][aa_name1],
                        info_list[ann_dict['hgvsp'.lower()]]
                    )
        if not parse_ann_list:
            parse_ann_list.append(info_list)
        else:
            if len(gene_impact) == 0:
                print(info_list);sys.exit()
            impact_num = vep_info_dict['impact'][gene_impact.lower()]
            list_num = vep_info_dict['impact'][parse_ann_list[0][ann_dict['impact'.lower()]].lower()]
            if impact_num > list_num:
                parse_ann_list = []; parse_ann_list.append(info_list)
            elif impact_num == list_num:
                parse_ann_list.append(info_list)
    return parse_ann_list

if __name__ == '__main__':
    if re.search('^annot$',sys.argv[1],re.IGNORECASE):
        print( annot(*sys.argv[2:]) )
    else:
        print('Ex:')
        print('\tpython3 *py annot input+file+*.vcf ......;')


