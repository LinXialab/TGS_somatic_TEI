import os
import pandas as pd


samtools = "/biosoft/samtools-1.9/samtools"
bedtools = "/biosoft/bedtools-2.30.0"
sniffles2 = "/biosoft/sniffles2/bin/sniffles"
minimap2 = "/biosoft/minimap2-2.17/minimap2"


hg38_fa = "/nanopore/ref/hg38_mainChr.fa"
tr_bed = "/nanopore/ref/human_GRCh38_no_alt_analysis_set.trf.bed"
repeat = tr_bed
repeatmasker_dir = '/nanopore/ref/repeat'


work_dir = "/nanopore/NewSample"

ms_dir = os.path.join(work_dir, '02.minimap2_sniffles')
coverage_dir = os.path.join(work_dir,'03.coverage')
somatic_dir = os.path.join(work_dir, "04.Somatic_minimap2_2.4")
Iris_dir = os.path.join(work_dir, 'Iris_test')
Somatic_TEI_dir = os.path.join(work_dir, 'somatic_INS_TEI_dir')

all_file_list = ['Sample1']



## 20230925
def reannotation_region_and_ins_region(file_csv_x,TYPE_repeat):
    final_decision = pd.DataFrame()
    final_decision['address'] = [x for x in file_csv_x['clean_decision'].split('|')]
    df1 = final_decision['address'].str.split('_', expand=True)
    df1 = df1.rename(columns={0: 'start', 1: 'end', 2: 'TYPE'})
    df1['TYPE'] = df1['TYPE'].str.split('/', expand=True)[0]
    df1_Simple_repeat = df1[df1['TYPE'] == TYPE_repeat]
    df1_Simple_repeat['start'] = df1_Simple_repeat['start'].astype(int)
    df1_Simple_repeat['end'] = df1_Simple_repeat['end'].astype(int)
    df1_Simple_repeat['anno_len'] = df1_Simple_repeat['end'] - df1_Simple_repeat['start']
    df1_Simple_repeat = df1_Simple_repeat.sort_values(by='start')
    df1_Simple_repeat['finally_decision']  =df1_Simple_repeat.apply(lambda x: 1 if ((x['start'] <= 52) or (x['end'] >= (48 + int(file_csv_x['total_len']) - 100))) and ((x['start'] <= (48 + int(file_csv_x['total_len']) - 100)) and (x['end'] >= 52)) else 0, axis=1)
    df1_Simple_repeat_1 = df1_Simple_repeat[df1_Simple_repeat['finally_decision'] == 1]
    return df1_Simple_repeat_1



def reannotation_region_TE_judgement_TE(file_csv_x):
    final_decision = pd.DataFrame()
    final_decision['address'] = [x for x in file_csv_x['clean_decision'].split('|')]
    df1 = final_decision['address'].str.split('_', expand=True)
    df1 = df1.rename(columns={0: 'start', 1: 'end', 2: 'TYPE'})
    df1['TYPE'] = df1['TYPE'].str.split('/', expand=True)[0]
    df1_Simple_repeat = df1[((df1['TYPE'] != 'Low') &(df1['TYPE'] != 'Simple') & (df1['TYPE'] != 'Tandem')  & (df1['TYPE'] != 'Satellite')) ]
    # df1_Simple_repeat = df1
    df1_Simple_repeat['start'] = df1_Simple_repeat['start'].astype(int)
    df1_Simple_repeat['end'] = df1_Simple_repeat['end'].astype(int)
    df1_Simple_repeat['anno_len'] = df1_Simple_repeat['end'] - df1_Simple_repeat['start']
    df1_Simple_repeat = df1_Simple_repeat.sort_values(by='start')
    print(df1_Simple_repeat)
    reslute = 'NA'
    if df1_Simple_repeat.shape[0] == 0:
        reslute = 'de_novo|te'
    else:
        df1_Simple_repeat['finally_decision'] = df1_Simple_repeat.apply(lambda x: 1 if ((x['start'] <= 50  and  50 <= x['end'] <= (int(file_csv_x['total_len']) - 50))
                                                                                        or (50 <= x['start']<=(int(file_csv_x['total_len'])-50) and (int(file_csv_x['total_len']) - 50) <= x['end'] <= int(file_csv_x['total_len']) )
                                                                                        or ((50 <= x['start']<=(int(file_csv_x['total_len'])-50)) and 50 <= x['end'] <= (int(file_csv_x['total_len']) - 50)))  else 0, axis=1)  ## 注释区域跨过或者再INS之中
        df1_Simple_repeat_1 = df1_Simple_repeat[df1_Simple_repeat['finally_decision'] == 1]
        if df1_Simple_repeat_1.shape[0] > 0:
            if df1_Simple_repeat_1.shape[0] == 1:
                reslute = 'TE|solo|' + df1_Simple_repeat_1['TYPE'].values[0]
            else:
                df1['TYPE'] = df1.apply(lambda x: 'Low_complex' if x['TYPE'] == 'Low' else 'Simple_repeat' if x['TYPE'] == 'Simple' else 'Tandem_repeat' if x['TYPE'] == 'Tandem' else x['TYPE'], axis=1)
                TE_complex = '|'.join(df1['TYPE'].values.tolist())  # complex  加上TD的类型
                reslute = 'TE|complex|' + TE_complex
        else:
            reslute = 'de_novo|te'
    return reslute



def reannotation_region_TE_judgement_TD(file_csv_x):
    final_decision = pd.DataFrame()
    final_decision['address'] = [x for x in file_csv_x['clean_decision'].split('|')]
    df1 = final_decision['address'].str.split('_', expand=True)
    df1 = df1.rename(columns={0: 'start', 1: 'end', 2: 'TYPE'})
    df1['TYPE'] = df1['TYPE'].str.split('/', expand=True)[0]
    df1_Simple_repeat = df1[((df1['TYPE'] != 'Low') &(df1['TYPE'] != 'Simple') & (df1['TYPE'] != 'Tandem')  & (df1['TYPE'] != 'Satellite')) ]
    df1_Simple_repeat['start'] = df1_Simple_repeat['start'].astype(int)
    df1_Simple_repeat['end'] = df1_Simple_repeat['end'].astype(int)
    df1_Simple_repeat['anno_len'] = df1_Simple_repeat['end'] - df1_Simple_repeat['start']
    df1_Simple_repeat = df1_Simple_repeat.sort_values(by='start')
    print(df1_Simple_repeat)
    reslute = 'NA'
    if df1_Simple_repeat.shape[0] == 0:
        reslute = 'de_novo|td'
    else:
        df1_Simple_repeat['finally_decision'] = df1_Simple_repeat.apply(lambda x: 1 if ((x['start'] <= 50  and  50 <= x['end'] <= (int(file_csv_x['total_len']) - 50))
                                                                                        or (50 <= x['start']<=(int(file_csv_x['total_len'])-50) and (int(file_csv_x['total_len']) - 50) <= x['end'] <= int(file_csv_x['total_len']) )
                                                                                        or ((50 <= x['start']<=(int(file_csv_x['total_len'])-50)) and 50 <= x['end'] <= (int(file_csv_x['total_len']) - 50)))  else 0, axis=1)  ## 注释区域跨过或者再INS之中
        df1_Simple_repeat_1 = df1_Simple_repeat[df1_Simple_repeat['finally_decision'] == 1]
        if df1_Simple_repeat_1.shape[0] > 0:
            if df1_Simple_repeat_1.shape[0] == 1:
                reslute = 'TE|solo|' + df1_Simple_repeat_1['TYPE'].values[0]
            else:
                df1['TYPE'] = df1.apply(lambda x: 'Low_complex' if x['TYPE'] == 'Low' else 'Simple_repeat' if x['TYPE'] == 'Simple' else 'Tandem_repeat' if x['TYPE'] == 'Tandem' else x['TYPE'], axis=1)
                TE_complex = '|'.join(df1['TYPE'].values.tolist())   ## complex TD
                reslute = 'TE|complex|' + TE_complex
        else:
            reslute = 'de_novo|td'
    return reslute



## 20231017 TEI type
def reannotation_region_TE(file_csv_x):
    final_decision = pd.DataFrame()
    final_decision['address'] = [x for x in file_csv_x['clean_decision'].split('|')]
    df1 = final_decision['address'].str.split('_', expand=True)
    df1 = df1.rename(columns={0: 'start', 1: 'end', 2: 'TYPE'})
    df1['TYPE'] = df1['TYPE'].str.split('/', expand=True)[0]
    df1_Simple_repeat = df1
    if df1_Simple_repeat.shape[0] == 1:
        reslute_TE= 'TE|solo|'+df1_Simple_repeat['TYPE'].values[0]
    else:
        TE_complex = '|'.join(df1_Simple_repeat['TYPE'].values.tolist())
        reslute_TE = 'TE|complex|' + TE_complex
    return reslute_TE




def get_coordinate_and_type_cor(file_csv_x, reslute):
    coordinate = 'Unmask'
    type_cor = 'Unmask'
    if 'TD' in reslute:
        coordinate = file_csv_x['trf_cor']
        type_cor = file_csv_x['trf_seq']
        if coordinate == 0 or coordinate == '0' or type_cor == 0 or type_cor == '0':
            coordinate = file_csv_x['sdust_cor']
            type_cor = file_csv_x['trf_ratio']
    elif 'TE' in reslute:
        coordinate = file_csv_x['coordinate']
        type_cor = file_csv_x['type_cor']
        if coordinate == 0 or coordinate == '0':
            coordinate = file_csv_x['sdust_cor']
            type_cor = file_csv_x['trf_ratio']
    elif 'de_novo' in reslute:
        coordinate = 'Unmask'
        type_cor = 'Unmask'
    return coordinate, type_cor



def get_TE_TD_de_novo(file_csv_x):
    # final_decision = pd.DataFrame()
    reslute = 'None'
    TE_goon = 'None'
    print(file_csv_x)
    if file_csv_x["total_len"] != 'total_len':
        if 'Simple_repeat' in file_csv_x['clean_decision']:
            TYPE_repeat = 'Simple'
            df1_Simple_repeat_1 = reannotation_region_and_ins_region(file_csv_x,TYPE_repeat)
            if df1_Simple_repeat_1.shape[0] >0 :
                reslute =  'TD|Simple_repeat_expasion'
            else:
                TE_goon = 'TE_judge'
                reslute = 'de_novo|Simple_repeat'
        elif 'Tandem_repeat' in file_csv_x['clean_decision']:
            TYPE_repeat = 'Tandem'
            df1_Simple_repeat_1 = reannotation_region_and_ins_region(file_csv_x,TYPE_repeat)
            if df1_Simple_repeat_1.shape[0] >0 :
                reslute =  'TD|Tandem_repeat_expasion'
            else:
                TE_goon = 'TE_judge'
                reslute = 'de_novo|Tandem_repeat'
        elif  'Satellite'  in file_csv_x['clean_decision']:
            TYPE_repeat = 'Satellite'
            df1_Simple_repeat_1 = reannotation_region_and_ins_region(file_csv_x,TYPE_repeat)
            if df1_Simple_repeat_1.shape[0] >0 :
                reslute =  'TD|Satellite_expasion'
            else:
                TE_goon = 'TE_judge'
                reslute = 'de_novo|Satellite'
        elif  'Low_complex'  in file_csv_x['clean_decision']:
            TYPE_repeat = 'Low'
            df1_Simple_repeat_1 = reannotation_region_and_ins_region(file_csv_x,TYPE_repeat)
            if df1_Simple_repeat_1.shape[0] >0 :
                reslute =  'TD|Low_complex_expasion'
            else:
                TE_goon = 'TE_judge'
                reslute = 'de_novo|Low_complex_expasion'
        elif 'Unmask' in file_csv_x['clean_decision']:
            reslute = 'de_novo|Unmask'
            # print('de_novo')
        elif('Low_complex'  not  in file_csv_x['clean_decision'])  and  ('Tandem_repeat'  not  in file_csv_x['clean_decision'])  and ('Simple_repeat' not in file_csv_x['clean_decision']) and ('Unmask' not in file_csv_x['clean_decision'])and ('Satellite' not in file_csv_x['clean_decision']):
            reslute = reannotation_region_TE_judgement_TE(file_csv_x)
        if TE_goon == 'TE_judge':
            reslute = reannotation_region_TE_judgement_TD(file_csv_x)
    coordinate, type_cor = get_coordinate_and_type_cor(file_csv_x, reslute)
    print(reslute, coordinate, type_cor)
    return (reslute, coordinate, type_cor)



for sampel in all_file_list:
    file_path_01 = '{Iris_dir}/ALL_reannotation/{sampel}'.format(sampel=sampel,Iris_dir=Iris_dir)
    file = os.path.join(file_path_01, 'all.INS.sdust.trf.replaced.cor.type.final.tsv')
    file_csv_01 = pd.read_csv(file, sep='\t')
    file_csv_01.loc[:, 'TE_TD_de_novo_type'] = 0
    file_csv_01.loc[:, 'TE_TD_de_novo_type_coordinate'] = 0
    file_csv_01.loc[:, 'TE_TD_de_novo_type_type_cor'] = 0
    file_csv_01 = file_csv_01.loc[:,~file_csv_01.columns.duplicated()]
    file_csv_01['TE_TD_de_novo_type'] ,file_csv_01['TE_TD_de_novo_type_coordinate'],file_csv_01['TE_TD_de_novo_type_type_cor'] = zip(*file_csv_01.apply(lambda x: get_TE_TD_de_novo(x), axis=1))
    file_csv_01.to_csv(os.path.join(file_path_01, 'all.INS.sdust.trf.replaced.cor.type.TE_TD_de_novo_type.tsv'),
                       sep='\t')






for sample in all_file_list:
    ins_csv = '{Iris_dir}/ALL_reannotation/{sampel}/all.INS.sdust.trf.replaced.cor.type.TE_TD_de_novo_type.tsv'.format(sampel=sample,Iris_dir=Iris_dir)
    all_ins_1 = pd.read_csv(ins_csv, sep='\t')
    all_ins_1 = all_ins_1[all_ins_1['dup_region'] != 'dup_region']
    all_ins  = pd.DataFrame()
    all_ins[['sample', 'old_type', 'svid']] = all_ins_1.apply(lambda x: split_row(x), axis=1)
    all_ins[['chr', 'start_2', 'end_2']] = all_ins_1['dup_region'].str.split('_', expand=True)
    all_ins['start'] = all_ins.apply(lambda x:min(int(x['start_2']),int(x['end_2'])),axis=1)
    all_ins['end'] = all_ins.apply(lambda x:max(int(x['start_2']),int(x['end_2'])),axis=1)
    all_ins['TE_TD_de_novo_type'] = all_ins_1['TE_TD_de_novo_type']
    all_ins['TE_TD_de_novo_type_2'] = all_ins_1['TE_TD_de_novo_type'].apply(lambda x:x.split('|')[0])
    all_ins['coordinate'] = all_ins_1['TE_TD_de_novo_type_coordinate']
    # all_ins['trf_seq'] = all_ins_1['trf_seq']
    all_ins['type_cor'] = all_ins_1['TE_TD_de_novo_type_type_cor']
    subdir = os.path.join(Somatic_TEI_dir, sample)
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    polish_reslut_path = '/nanopore/NewSample/02.polish_4_times_andsniffles_test/{sample}/{sample}_vcf_region_df_only_positive_polish.csv'.format(sample=sample)
    if os.path.exists(polish_reslut_path):
        sample_sv = pd.read_csv(polish_reslut_path, sep='\t')
        sample_sv = sample_sv.dropna(subset=['chr'])
        sample_sv = sample_sv.loc[:,~sample_sv.columns.duplicated()]
        sample_sv['somSV'] = sample_sv['somSV'].apply(lambda x: 'True_somatic_SV' if x == 'True_somatic_SV' or x == 'True_somatic_SV.support_reads_less_half' or x == 'True_somatic_SV.support_reads_over_half' else x)
        sample_sv['svid'] = sample_sv['sv_id']
        sample_sv = sample_sv[['svid', 'dir_name', 'sv_length', 'somSV']]
        sample_sv = sample_sv[sample_sv['somSV'] == 'True_somatic_SV']
        print(sample)
        sample_ins = all_ins[all_ins['sample'] == sample]
        sample_ins[sample_ins['TE_TD_de_novo_type_2'] == 'de_nove'] = sample_ins[
            sample_ins['TE_TD_de_novo_type_2'] == 'de_nove'].replace('de_nove', 'de_novo')
        for sv in ['TD', 'TE', 'de_novo','all']:
            if sv == 'all':
                subins = sample_ins
            else:
                subins = sample_ins[sample_ins['TE_TD_de_novo_type_2']== sv]
            if subins.shape[0] != 0:
                subsv = sample_sv[sample_sv['svid'].isin(subins['svid'])]
                subdf = pd.merge(subins, subsv, on=['svid'])
                subdf1 = subdf[['chr', 'start', 'end', 'TE_TD_de_novo_type','coordinate','type_cor','svid', 'dir_name', 'sv_length', 'somSV']]
                subdf1 = subdf1.sort_values(by=['chr', 'start', 'end'])  # sort by 'chr_x','start_x','end_x'
                subdf1.to_csv(os.path.join(subdir, '_'.join([sample, 'tumor_somatic', sv, 'ins.csv'])), sep='\t', index=False)
                subdf_bk_start = subdf1.copy(deep=True)
                subdf_bk_start['end'] = subdf_bk_start['start']
                subdf_bk_start['start'] = subdf_bk_start['start'].map(int) - 1
                subdf_bk_start['start'] = subdf_bk_start['start'].astype(str)
                subdf_bk_end = subdf1.copy(deep=True)
                subdf_bk_end['start'] = subdf_bk_end['end']
                subdf_bk_end['end'] = subdf_bk_end['end'].map(int) + 1
                subdf_bk_end['end'] = subdf_bk_end['end'].astype(str)
                subdf_bk = pd.concat([subdf_bk_start, subdf_bk_end], axis=0)
                subdf_bk = subdf_bk.sort_values(by=['chr', 'start', 'end'])
                subdf_bk.to_csv(os.path.join(subdir, '_'.join([sample, 'tumor_somatic', sv, 'breakpoint']) + '.bed'),index=0, header=0, sep='\t')
                subdf1.to_csv(os.path.join(subdir, '_'.join([sample, 'tumor_somatic', sv]) + '.bed'), sep='\t', index=0, header=0)


















