#!/usr/bin/env python3
#coding:utf-8
import os
import re
import sys
import gzip

def exac(exac_file = None):
    if exac_file == None:
        print('Exac File Not Set!');
        sys.exit()
    exac_file = os.path.abspath(exac_file.strip())
    exac_h = ''
    if re.search('\.gz$',exac_file):
        exac_h = gzip.open(exac_file,'rt',encoding='utf-8')
    else:
        exac_h = open(exac_file,'r',encoding='utf-8')
    input_name = os.path.basename(exac_file)
    output_name = re.sub( '.vcf'+re.split('\.vcf',input_name)[-1],'.xls',input_name )
    output_file = os.path.join(os.path.dirname(exac_file),output_name)
    output_h = open(output_file,'w',encoding='utf-8')
    output_dict = {
        'AF':1, 'DP':1,
        'AC':1, 'AC_Adj':1, #'AC_Hemi':1, 'AC_Het':1, 'AC_Hom':1, 'AC_MALE':1, 'AC_FEMALE':1,
        'AN':1, 'AN_Adj':1, #'AN_Hemi':1, 'AN_Het':1, 'AN_Hom':1, 'AN_MALE':1, 'AN_FEMALE':1,
        'AC_AFR':1, 'AC_AMR':1, 'AC_EAS':1, 'AC_FIN':1, 'AC_NFE':1, 'AC_OTH':1, 'AC_SAS':1,
        'AN_AFR':1, 'AN_AMR':1, 'AN_EAS':1, 'AN_FIN':1, 'AN_NFE':1, 'AN_OTH':1, 'AN_SAS':1,
        #'Hemi_AFR':1, 'Hemi_AMR':1, 'Hemi_EAS':1, 'Hemi_FIN':1, 'Hemi_NFE':1, 'Hemi_OTH':1, 'Hemi_SAS':1,
        #'Het_AFR': 1, 'Het_AMR': 1, 'Het_EAS': 1, 'Het_FIN': 1, 'Het_NFE': 1, 'Het_OTH': 1, 'Het_SAS': 1,
        #'Hom_AFR': 1, 'Hom_AMR': 1, 'Hom_EAS': 1, 'Hom_FIN': 1, 'Hom_NFE': 1, 'Hom_OTH': 1, 'Hom_SAS': 1,
        'ESP_AF_GLOBAL':1, 'ESP_AC':1,
        'KG_AF_GLOBAL':1, 'KG_AC':1,
    }
    output_h.write('\t'.join( ['chr','pos','ref','alt'] ))
    for tag_name in sorted(output_dict.keys()):
        output_h.write('\t%s' % tag_name)
    output_h.write('\n')
    for line in exac_h:
        if re.match('#',line) == None:
            infos = re.split('\t',line.strip())
            temp_dict = {}
            total_len = 1
            target_num = 0
            de_str = list(filter( lambda x:re.match('OLD_MULTIALLELIC',x),re.split(';',infos[7]) ))
            var_str = list(filter( lambda x:re.match('OLD_VARIANT',x),re.split(';',infos[7]) ))
            if len(de_str) >= 1:
                de_dict = { t_base: t_num for t_base,t_num in enumerate(re.split('/',re.split(':',de_str[0])[-1])) }
                total_len = len(de_dict)
                target_base = infos[4] if len(var_str) == 0 else re.split('/',re.split(':',var_str[0])[-1])[-1]
                if target_base in de_dict:
                    target_num = de_dict[target_base]
            for info_str in re.split(';',infos[7]):
                info_list = re.split('=',info_str)
                if info_list[0] in output_dict:
                    detail_list = re.split(',',info_list[1])
                    #print('test');print(target_num);print(total_len);print(info_str)
                    if len(detail_list) == total_len:
                        temp_dict[info_list[0]] = detail_list[target_num]
                    elif len(detail_list)+1 == total_len:
                        temp_dict[info_list[0]] = detail_list[target_num-1]
                    elif re.search('(^AN)|(\_AN)|DP|(^AF$)',info_list[0]):
                        temp_dict[info_list[0]] = info_list[-1]
                    else:
                        print('Not Find Pos');print(info_str);
                        print(de_str); print(var_str); print(target_num)
                        sys.exit()
            output_h.write('\t'.join( [infos[0],infos[1],infos[3],infos[4]] ))
            for tag_name in sorted(output_dict.keys()):
                temp_dict[tag_name] = 'Unknown' if tag_name not in temp_dict else temp_dict[tag_name]
                output_h.write('\t%s' % temp_dict[tag_name])
            output_h.write('\n')
    return output_file

if __name__ == '__main__':
    print( exac(sys.argv[1]) )
    print('End')







