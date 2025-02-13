import os
import pysam
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool

racon = "/biosoft/racon/build/bin/racon"
samtools = "/biosoft/samtools-1.9/samtools"   ##1.9
minimap2 = "/biosoft/minimap2-2.17/minimap2"    ##2.17
sniffles = "/biosoft/sniffles2/bin/sniffles"
seqkit = '/biosoft/seqkit'
bedtools = "/biosoft//bedtools-2.30.0"



work_dir = "/nanopore/NewSample"

ms_dir = os.path.join(work_dir, '02.minimap2_sniffles')
coverage_dir = os.path.join(work_dir,'03.coverage')
somatic_dir = os.path.join(work_dir, "04.Somatic_minimap2_2.4")




def get_blood_bam(ParamList):
    '''
    return SV whether should reverse polish
    input:
    ParamList:
        work_dir
        tumor_bam_path
        blood_bam_path
        vcf_region_df: vcf_region_df
        index:
    output:
        record {sampel_name_x}_blood.bam
    '''
    # region, blood_bam_path, wok_dir, dir_name
    wok_dir, tumor_bam_path, blood_bam_path, vcf_region_df, index = ParamList
    line_x = vcf_region_df.loc[index]
    dir_name = line_x['dir_name']
    region = line_x['region']
    dir_path = os.path.join(wok_dir, dir_name)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    bam_out_region = os.path.join(wok_dir, dir_name, '{sampel_name_x}_region_blood.bam'.format(sampel_name_x=dir_name))
    bam_out = os.path.join(wok_dir, dir_name, '{sampel_name_x}_blood.bam'.format(sampel_name_x=dir_name))
    os.system(
        "{samtools} view -@ 10 -h -b {bam_ont} {region} > {bam_out}".format(samtools=samtools, bam_ont=blood_bam_path,
                                                                            region=region, bam_out=bam_out_region))
    cmd = "({samtools} view -H {inbam}; {samtools} view {inbam} | grep 'blood')| {samtools} view -b - -o {out}".format(
        samtools=samtools, inbam=bam_out_region,out=bam_out)
    os.system(cmd)
    os.system("rm {bam_out_region}".format(bam_out_region=bam_out_region))
    os.system("{samtools} index {bam_out}".format(samtools=samtools, bam_out=bam_out))


def get_tumor_support_reads_bam(ParamList):
    '''
    return SV whether should reverse polish
    input:
    ParamList:
    work_dir
    tumor_bam_path
        blood_bam_path
        vcf_region_df: vcf_region_df
        index:
        # wok_dir: bamFile dir
        # SVinfo: SVinfo df
        # index: SV PID
    output:
       {sampel_name_x}_tumor.bam  {sampel_name_x}_tumor_support_reads.fq
    '''
    # region, tumor_bam_path, wok_dir, line_x
    wok_dir, tumor_bam_path, blood_bam_path, vcf_region_df, index = ParamList
    line_x = vcf_region_df.loc[index]
    dir_name = line_x['dir_name']
    region = line_x['region']
    print(dir_name)
    dir_path = os.path.join(wok_dir, dir_name)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    bam_out_region = os.path.join(wok_dir, dir_name, '{sampel_name_x}_region_tumor.bam'.format(sampel_name_x=dir_name))
    bam_out = os.path.join(wok_dir, dir_name, '{sampel_name_x}_tumor.bam'.format(sampel_name_x=dir_name))
    region_tumor_bam_parh_x_fq = os.path.join(dir_path, '{sampel_name_x}_tumor.fq'.format(sampel_name_x=dir_name))
    os.system(
        "{samtools} view -@ 10 -h -b {bam_ont} {region} > {bam_out}".format(samtools=samtools, bam_ont=tumor_bam_path,
                                                                            region=region, bam_out=bam_out_region))
    cmd = "({samtools} view -H {inbam}; {samtools} view {inbam}| grep 'tumor')| {samtools} view -b - -o {out}".format(
        samtools=samtools, inbam=bam_out_region,out=bam_out)
    os.system(cmd)
    os.system("rm {bam_out_region}".format(bam_out_region=bam_out_region))
    os.system("{samtools} index {bam_out}".format(samtools=samtools, bam_out=bam_out))
    os.system("{samtools} fastq {bam_out} > {fq_out}".format(samtools=samtools, bam_out=bam_out,
                                                             fq_out=region_tumor_bam_parh_x_fq))  # Get tumor'fastq




def GetSoftClip(bam):
    '''
    return soft clip location for each read
    input:
        opened bam file handel "pysam.AlignmentFile()"
    output:
        record format: [quaryName, referenceTag, lsoftclip_starts, lsoftclip_ends, rsoftclip_start, rsoftclip_end]
    '''
    record = []
    for reads in bam.fetch():
        if reads.is_secondary or reads.is_supplementary or reads.is_unmapped:
            # print('JJJJ')
            continue
        else:
            cigar = reads.get_aligned_pairs()
            # cigar for soft clip on the left and right of aligned pairs
            ReadLen = cigar[-1][0]
            A_Match = np.array([cigar[i] for i in range(len(cigar)) if cigar[i][1] != None])
            l_softsite, l_slen = A_Match[0][1], A_Match[0][0]
            r_softsite, r_slen = A_Match[-1][1], ReadLen - A_Match[-1][0]
            record.append((reads.qname, reads.qlen, reads.reference_name, l_softsite, l_slen, r_softsite, r_slen))
    return (record)





def polish_four_times(work_dir_x, blood_sequence_fq_x, sampel_name_x):
    '''
   return soft clip location for each read
   input:
       {sampel_name_x}_blood.fq
   output:
       {sampel_name_x}.polished.round6.fasta
    #    '''
    for i in [1, 2, 3, 4,5,6]:
        print('{sampel_name_x} in {next_n} times polish'.format(sampel_name_x=sampel_name_x, next_n=i))
        work_dir = work_dir_x
        previous_sequence = work_dir + '/{sampel_name_x}.polished.round{last_n}.fasta'.format(
            sampel_name_x=sampel_name_x,
            last_n=i - 1)
        Up_to_data_sam = work_dir + '/{sampel_name_x}.polished.R{next_n}.sam'.format(sampel_name_x=sampel_name_x,
                                                                                     next_n=i)
        Up_to_data_bam = work_dir + '/{sampel_name_x}.polished.R{next_n}.bam'.format(sampel_name_x=sampel_name_x,
                                                                                     next_n=i)
        Up_to_data_sequence = work_dir + '/{sampel_name_x}.polished.round{next_n}.fasta'.format(
            sampel_name_x=sampel_name_x,
            next_n=i)
        if i == 1:
            previous_sequence = '/nanopore/ref/hg38_mainChr.fa'
        os.system(
            '{minimap2} --MD -ax map-ont -t 10 {previous_sequence} {blood_sequence_fq}  > {Up_to_data_sam}'.format(minimap2=minimap2,
                previous_sequence=previous_sequence, blood_sequence_fq=blood_sequence_fq_x,
                Up_to_data_sam=Up_to_data_sam))  # map blood fastq with last time fasta
        os.system(
            '{racon} --threads 10  {blood_sequence} {Up_to_data_sam} {previous_sequence}  > {Up_to_data_sequence}'.format(racon=racon,
                blood_sequence=blood_sequence_fq_x, Up_to_data_sam=Up_to_data_sam, previous_sequence=previous_sequence,
                Up_to_data_sequence=Up_to_data_sequence))  # racon sam with blood fasta (fist time is ref)
        # del polish filter
        if i in [1, 2, 3, 4,5,6]:
            os.system('rm {Up_to_data_sam}'.format(Up_to_data_sam=Up_to_data_sam))
        if i in [2, 3, 4,5,6]:
            os.system('rm {previous_sequence}'.format(previous_sequence=previous_sequence))  # del before fasta
            print('del {previous_sequence}'.format(previous_sequence=previous_sequence))
        if i == 6:
            os.system('samtools faidx {Up_to_data_sequence}'.format(Up_to_data_sequence=Up_to_data_sequence))






def polish_4_times_blood(wok_dir, dir_name):
    '''
   return soft clip location for each read
   input:
       wok_dir, dir_name
   output:
       {sampel_name_x}_blood.fq {sampel_name_x}_blood.bam {sampel_name_x}_blood.bam.bai
   '''
    # wok_dir, SVinfo, index = ParamList
    dir_path = os.path.join(wok_dir, dir_name)
    region_blood_bam_parh_x = os.path.join(dir_path, '{sampel_name_x}_blood.bam'.format(sampel_name_x=dir_name))
    region_blood_bam_parh_x_fq = os.path.join(dir_path, '{sampel_name_x}_blood.fq'.format(sampel_name_x=dir_name))
    os.system(
        ''' %s fastq %s > %s ''' % (samtools, region_blood_bam_parh_x, region_blood_bam_parh_x_fq))  # bam to fastq
    work_dir_polish = os.path.join(dir_path, 'polish')
    if not os.path.exists(work_dir_polish):
        os.makedirs(work_dir_polish)
    # 4 times polish
    dir_name_cancer_type = dir_name + '_blood'
    # print(dir_name_cancer_type)
    polish_four_times(work_dir_polish, region_blood_bam_parh_x_fq, dir_name_cancer_type)







def minimap2_and_sniffle(wok_dir, dir_name):
    '''
    return  无
    input:
        wok_dir, dir_name
    output:
        record format: [quaryName, referenceTag, lsoftclip_starts, lsoftclip_ends, rsoftclip_start, rsoftclip_end]
    '''
    print(dir_name + 'Start run minimap2 and sniffles')
    dir_path = os.path.join(wok_dir, dir_name)
    tumor_fq_path = os.path.join(dir_path, '{sampel_name_x}_tumor.fq'.format(sampel_name_x=dir_name))
    blood_fq_path = os.path.join(dir_path, '{sampel_name_x}_blood.fq'.format(sampel_name_x=dir_name))
    merge_fq_path = os.path.join(dir_path, '{sampel_name_x}_merge.fq'.format(sampel_name_x=dir_name))
    os.system('cat {tumor_fq_path} {blood_fq_path} > {merge_fq_path} ' .format(tumor_fq_path =tumor_fq_path ,blood_fq_path =blood_fq_path,merge_fq_path=merge_fq_path))
    blood_fa_path = dir_path + '/polish' + '/{sampel_name_x}_blood.polished.round{next_n}.fasta'.format(
        sampel_name_x=dir_name,
        next_n=6)
    dir_minimap_path = os.path.join(dir_path, 'polish')
    if not os.path.exists(dir_minimap_path):
        os.makedirs(dir_minimap_path)
    out_bam = os.path.join(dir_minimap_path, '{sampel_name_x}_ref_blood_polish.bam'.format(sampel_name_x=dir_name))
    out_vcf = os.path.join(dir_minimap_path, '{sampel_name_x}_ref_blood_polish.vcf'.format(sampel_name_x=dir_name))
    cmd = '{minimap2} --MD -ax map-ont -t 10 {ref} {fastq}  | /NAS/wg_fzt/software/samtools-1.9/samtools sort -O BAM -o {out_bam}'.format(
        minimap2=minimap2,
        ref=blood_fa_path, fastq=merge_fq_path, out_bam=out_bam)  # 比较得到bam文件
    os.system(cmd)
    cmd_2 = '{samtools} index {out_bam}'.format(samtools=samtools, out_bam=out_bam)
    os.system(cmd_2)
    cmd_3 = "{sniffles} -i {bam} -v {out_file} --minsupport 3 --cluster-binsize 500 --minsvlen 50  --phase -t 10 --mapq 40 " \
            "--min-alignment-length 1000 --output-rnames --allow-overwrite --long-ins-length 100000 " \
            "--reference {ref}".format(sniffles = sniffles,bam=out_bam, out_file=out_vcf,ref=blood_fa_path)
    os.system(cmd_3)







def if_reverse_polish_by_softclipe(ParamList):
    '''
    return SV whether should reverse polish
    input:
    ParamList:
        work_dir
        tumor_bam_path
        blood_bam_path
        vcf_region_df: vcf_region_df
        index:
        # wok_dir: bamFile dir
        # SVinfo: SVinfo df
        # index: SV PID
    output:
       {sampel_name_x}_tumor.bam  {sampel_name_x}_tumor_support_reads.fq
    '''
    wok_dir, tumor_bam_path, blood_bam_path, vcf_region_df, index = ParamList
    line_x = vcf_region_df.loc[index]
    dir_name = line_x['dir_name']
    print(dir_name)
    dir_path = os.path.join(wok_dir, dir_name)
    dir_minimap_path = os.path.join(dir_path, 'polish')
    out_vcf_blood = os.path.join(dir_minimap_path, '{sampel_name_x}_ref_blood_polish.vcf'.format(sampel_name_x=dir_name))
    dir_minimap_path_2  = os.path.join(dir_path, 'polish_tumor')
    out_vcf_tumor = os.path.join(dir_minimap_path_2, '{sampel_name_x}_ref_tumor_polish.vcf'.format(sampel_name_x=dir_name))
    # print(dir_name)
    dir_path = os.path.join(wok_dir, dir_name)
    region_blood_bam_parh_x = os.path.join(dir_path, '{sampel_name_x}_blood.bam'.format(sampel_name_x=dir_name))
    bam = pysam.AlignmentFile(region_blood_bam_parh_x,'r',check_sq=False) # read blood bam
    softclipDf = pd.DataFrame.from_records(GetSoftClip(bam)) # get readsoftclipe information
    if len(softclipDf) == 0:
        print('blood all reads MAPQ = 0，can not in forward polish，') ##
        print('need reverse polish')
        return ((dir_name, 'reverse_polish'))
    softclipDf.columns = ['readID', 'readLen', 'chrom', 'l_softClip', 'l_softlen', 'r_softClip', 'r_softlen']
    softclipDf = softclipDf.sort_values(['chrom', 'l_softClip'], ascending=True) # sort by chrom and softclipe's location
    sv_start = int(line_x['start']) - 100
    sv_end = int(line_x['end']) + 100
    # filter softclipe'location in SV breakpoint +- 100bp and softclipe length >500bp reads
    softclipDf_condidata = softclipDf.loc[((softclipDf['l_softClip'] >= sv_start) & (softclipDf['l_softClip'] <= sv_end) & (
            softclipDf['l_softlen'] >= 500)) | ((softclipDf['r_softClip'] >= sv_start) & (
            softclipDf['r_softClip'] <= sv_end) & (softclipDf['r_softlen'] >= 500)), :]
    if softclipDf_condidata.shape[0] >= 0.5 * softclipDf.shape[0]: # softclipe'reads>0.5*reads
        print('softclipe reads over 0.5 all reads')
        print('need reverse polish')
        return ((dir_name, 'reverse_polish'))
    else:
        if not os.path.exists (out_vcf_blood):
            polish_4_times_blood(wok_dir, dir_name) # polish 6
            minimap2_and_sniffle(wok_dir,dir_name) # minimap2 and sniffle
            print('softclipe reads less 0.5 all reads')
            print('not need reverse polish')
        return ((dir_name, 'only_positive_polish'))



def ParseSV(vcf_path,line_x):
    with open(vcf_path,'r') as fin:
        records = [x.strip().split("\t") for x in fin.readlines() if not re.search('##', x)]
    vcfDf = pd.DataFrame.from_records(records[1:])
    if len(vcfDf) == 0:
        print('vcf is empty，not somatic SV ')
        sub = {}
        return sub
    else:
        vcfDf.columns = records[0]  ## 第0行为headr
        #### BreakPoint1
        vcfDf['CHR1'] = vcfDf['#CHROM']
        vcfDf['BRKP1'] = vcfDf['POS']
        vcfDf['PRECISE'] = vcfDf['INFO'].apply(lambda x:x.split(";")[0])
        ### BreakPoint2
        vcfDf['CHR2'] = vcfDf.apply(lambda x:x['INFO'].split("CHR2=")[-1].split(";")[0] if re.search("CHR2=",x['INFO']) else x['#CHROM'],axis=1)
        vcfDf['BRKP2'] = vcfDf.apply(lambda x:x['INFO'].split("END=")[-1].split(";")[0] if re.search("END=",x['INFO']) else re.match('[0-9]+',x['ALT'].split(":")[-1]).group(0),axis=1)
        vcfDf['SVType'] = vcfDf['INFO'].apply(lambda x:x.split("SVTYPE=")[-1].split(";")[0])
        #### +/- strands
        vcfDf['STRANDS'] = vcfDf['INFO'].apply(lambda x:x.split("STRAND=")[-1].split(";")[0])
        vcfDf['STRAND1'] = vcfDf['STRANDS'].apply(lambda x:x[0])
        vcfDf['STRAND2'] = vcfDf['STRANDS'].apply(lambda x:x[1] if len(x)==2 else x[0])
        ### SV length
        vcfDf['SVLen'] = vcfDf['INFO'].apply(lambda x:np.abs(int(x.split("SVLEN=")[-1].split(";")[0])) if re.search("SVLEN=",x) else 0)
        vcfDf['STD_LEN'] = vcfDf['INFO'].apply(lambda x:x.split("STDEV_LEN=")[-1].split(";")[0] if re.search("SVLEN=",x) else 0)
        vcfDf['STD_POS'] = vcfDf['INFO'].apply(lambda x:x.split("STDEV_POS=")[-1].split(";")[0])
        vcfDf['SNIF_COVERAGE'] = vcfDf['INFO'].apply(lambda x:x.split("COVERAGE=")[-1].split(";")[0])
        vcfDf['SNIF_AF'] = vcfDf['INFO'].apply(lambda x:x.split("AF=")[-1].split(";")[0])
        vcfDf['SNIF_GT'] = vcfDf['SAMPLE'].apply(lambda x:x.split(":")[0])
        vcfDf['SNIF_GQ'] = vcfDf['SAMPLE'].apply(lambda x:x.split(":")[1])
        vcfDf['SNIF_DR'] = vcfDf['SAMPLE'].apply(lambda x:x.split(":")[2])
        vcfDf['SNIF_DV'] = vcfDf['SAMPLE'].apply(lambda x:x.split(":")[3])
        #### Support Reads
        vcfDf['RNAMES'] = vcfDf['INFO'].apply(lambda x:x.split("RNAMES=")[-1].split(";")[0])
        ### Extract Tumor / Normal Support Reads
        vcfDf['TR'] = vcfDf.apply(lambda x: [r for r in x['RNAMES'].split(",") if re.search('tumor',r)], axis=1)
        vcfDf['NR'] = vcfDf.apply(lambda x: [r for r in x['RNAMES'].split(",") if re.search('blood',r)], axis=1)
        vcfDf['Support'] = vcfDf.apply(lambda x: x['TR']+x['NR'], axis=1)
        vcfDf['TS'] = vcfDf['TR'].apply(lambda x: len(x))
        vcfDf['NS'] = vcfDf['NR'].apply(lambda x: len(x))
        vcfDf['TR'] = vcfDf.apply(lambda x: ",".join([r for r in x['RNAMES'].split(",") if re.search('tumor',r)]), axis=1)
        vcfDf['NR'] = vcfDf.apply(lambda x: ",".join([r for r in x['RNAMES'].split(",") if re.search('blood',r)]), axis=1)
        ### remove mito SV
        vcfDf = vcfDf.loc[(vcfDf['CHR1']!='chrM') & (vcfDf['CHR2']!='chrM')]
        vcfDf['Patient'] = line_x['sample_name']
        vcfDf['PID'] = vcfDf.apply(lambda x:"-".join([x['Patient'],x['ID'],]),axis=1)
        ##ReadsOI
        vcfDf_ROI = vcfDf
        #### out
        sub = pd.DataFrame(vcfDf_ROI,columns=['Patient','PID','ID','SVType','SVLen',
                                              'CHR1','BRKP1','CHR2', 'BRKP2',
                                              'STRANDS','STRAND1','STRAND2',
                                              'ALT','FILTER','PRECISE',
                                              'STD_LEN','STD_POS','SNIF_COVERAGE','SNIF_AF','SNIF_GT','SNIF_GQ','SNIF_DR','SNIF_DV',
                                              'TR','NR','Support','TS','NS'])
        return(sub)




def if_somatic_sv_by_sniffles_vcf(ParamList):
    '''
    return SV whether should reverse polish
    input:
    ParamList:
        wok_dir: bamFile dir
        SVinfo: SVinfo df
        index: SV PID
    output:
        record {'NONE_UNkonw_vcf_is_empty'|'NONE_UNkonw_vcf_is_empty'|'NONE_UNkonw_wrong_last_step'|'False_somatic_SV'|'True_somatic_SV'}
    '''
    wok_dir, tumor_bam_path, blood_bam_path, vcf_region_df, index = ParamList
    line_x = vcf_region_df.loc[index]
    dir_name = line_x['dir_name']
    dir_path = os.path.join(wok_dir, dir_name)
    vcf_path = os.path.join(dir_path, 'polish', '{sampel_name_x}_ref_blood_polish.vcf'.format(sampel_name_x=dir_name))
    if not os.path.exists(vcf_path):
        return ((dir_name, 'NONE_UNkonw_wrong_last_step'))
    else:
        sub = ParseSV(vcf_path, line_x)
        ### Determine if the VCF file is empty
        if len(sub) == 0 :
            return ((dir_name, 'NONE_UNkonw_vcf_is_empty' ))
        else:
            number_same_chr_1 = sub[sub['CHR1'] == line_x['chr']].shape[0]
            number_same_chr_2 = sub[sub['CHR2'] == line_x['chr']].shape[0]
            ### Determine whether the CHR1 of the new vcf and the previous SV are on the same chromosome
            if number_same_chr_1 == 0 and number_same_chr_2 == 0:
                return ((dir_name, 'NONE_UNkonw_wrong_different_chr' ))
            else:
                ###  Extracting Somatic SV and its Information
                SVTable_Som = sub.loc[(sub['NS']==0)&(sub['TS']>=3)] ## support reads over 3
                SVTable_Som = SVTable_Som.loc[(SVTable_Som['CHR1'] != 'chrY') & (SVTable_Som['CHR2'] != 'chrY')]
                SVTable_Som['SVLen'] = SVTable_Som['SVLen'].astype(int)
                SVTable_INDEL  = SVTable_Som.loc[(SVTable_Som['SVType'] == 'INS') & (SVTable_Som['SVLen'] <= 50)]
                SVTable_INDEL = SVTable_INDEL.append(SVTable_Som.loc[(SVTable_Som['SVType'] == 'DEL') & (SVTable_Som['SVLen'] <= 50)])
                SVTable_Som = SVTable_Som.append(SVTable_INDEL).drop_duplicates(subset=['PID'],keep=False)
                SVTable_Som.index = SVTable_Som['PID']
                ### Determine if there is a somatic SV after Polish 2023/4/17更改
                if len(SVTable_Som) == 0:
                    return ((dir_name, 'False_somatic_SV'))
                else:
                    SVTable_Som.to_csv(os.path.join(dir_path, 'polish', '{sampel_name_x}_ref_blood_polish_somatic_SV.csv'.format(sampel_name_x=dir_name)),index=False)
                    SVTable_polish = line_x['support_reads_name'].split(',')
                    SVTable_polish = set(SVTable_polish)
                    SVTable_Som['Support'] = SVTable_Som['Support'].apply(lambda x: set(x))
                    SVTable_Som['Support'] = SVTable_Som['Support'].apply(lambda x: len(x.intersection(SVTable_polish)))
                    SVTable_Som['support'] = SVTable_Som['Support'].astype(int)
                    SVTable_Som['Support'] = SVTable_Som['Support'].apply(lambda x: x/max(len(SVTable_polish),x))
                    SVTable_Som['Support'] = SVTable_Som['Support'].apply(lambda x: True if x > 0.5 else False)
                    SVTable_Som['Support'] = SVTable_Som['Support'].astype(str)
                    SVTable_Som['Support'] = SVTable_Som['Support'].apply(lambda x: 'True_somatic_SV' if x == 'True' else 'False_somatic_SV')
                    if 'True_somatic_SV' in SVTable_Som['Support'].values:
                        return ((dir_name, 'True_somatic_SV.support_reads_over_half'))
                    else:
                        return ((dir_name, 'True_somatic_SV.support_reads_less_half'))



def make_new_df_from_vcf(path_vcf,sample_name):
    vcf_df = pd.read_csv(path_vcf,sep='\t',header=None,comment='#')
    vcf_region_df = pd.DataFrame()
    vcf_region_df['chr'] = vcf_df[0]
    vcf_region_df['start'] = vcf_df[1].apply(lambda x: str(int(x) - 50))
    vcf_region_df['end'] = vcf_df.apply(lambda x:str(abs(int(x[7].split(';')[3].split('=')[1])) + 50) if x[7].split('SVTYPE=') [-1].split(";")[0] != 'BND' else int(x[1]) + 50,axis=1)
    vcf_region_df['sv_id'] = vcf_df[2]
    vcf_region_df['sample_name'] = sample_name
    vcf_region_df['region'] = vcf_region_df['chr'].astype(str) + ':' + vcf_region_df['start'].astype(str) + '-' + vcf_region_df['end'].astype(str)
    vcf_region_df['dir_name'] = vcf_region_df['sample_name'] + '_' + vcf_region_df['sv_id'] + '_' + vcf_region_df['chr'] + '_' + vcf_region_df['start'].astype(str)
    vcf_region_df['support_reads_name'] = vcf_df[7].apply(lambda x: x.split('RNAMES=')[-1].split(';')[0])
    vcf_region_df['sv_type'] = vcf_df[7].str.split(';').str[1].str.split('=').str[1]
    vcf_region_df['sv_length'] = vcf_df[7].apply(lambda x:np.abs(int(x.split("SVLEN=")[-1].split(";")[0])) if re.search("SVLEN=",x) else 0)
    vcf_region_df.index = np.array(vcf_region_df['dir_name']) #set ndex as dir_name
    return vcf_region_df


def get_merge_bam_path(sample_name):
    tumor_bam = os.path.join(somatic_dir, 'Cancer_cells', sample_name,
                                      '%s_Cancer_cells_minimap2_sorted_tag.bam' % sample_name)
    sampelid_blood = sample_name
    blood_bam = os.path.join(somatic_dir, 'blood', sampelid_blood,
                                      '%s_blood_minimap2_sorted_tag.bam' % sampelid_blood)
    return tumor_bam, blood_bam


all_file_list = ['Sample1']

for sample_name in all_file_list:
    path_vcf = '/nanopore/NewSample/Iris_test/{OSCC5}/{OSCC5}_iris_out.vcf'.format(OSCC5=sample_name)
    work_dir = '/nanopore/NewSample/02.polish_4_times_andsniffles_test/{OSCC5}/'.format(OSCC5=sample_name)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    tumor_bam_path, blood_bam_path = get_merge_bam_path(sample_name)
    vcf_region_df = make_new_df_from_vcf(path_vcf,sample_name)
    ParamList = [(work_dir,tumor_bam_path,blood_bam_path,vcf_region_df, index) for index in vcf_region_df.index]
    print('Start tumor support reads bam')
    P = Pool(60)
    P.map(get_tumor_support_reads_bam, ParamList) #get tumorsupport sv bam and fastq
    P.close()
    del P
    print('Over  tumor support reads bam')
    print('Start blood support reads bam')
    P = Pool(60)
    P.map(get_blood_bam, ParamList) #get blood sv bam and fastq
    P.close()
    del P
    print('Over blood support reads bam')
    print('Start polish')
    P = Pool(60)
    RecordList = P.map(if_reverse_polish_by_softclipe, ParamList)
    P.close()
    del P
    print('polish Over')
    print('Out reslute')
    JudgementDf = pd.DataFrame.from_records(RecordList)
    JudgementDf.columns = ['dir_name', 'reverse_polish']
    JudgementDf.index = np.array(JudgementDf['dir_name'])
    ComDf = pd.concat([vcf_region_df, JudgementDf], axis=1) # merge dataframe,dir_name as index
    ComDf.to_csv(work_dir + '/{OSCC5}_the_resluct_of_polish.csv'.format(OSCC5 = sample_name),index=False,header=True,sep='\t')
    print(work_dir + '/{OSCC5}_the_resluct_of_polish.csv'.format(OSCC5 = sample_name))
    print('Complete polish and Out reslute')
    plt.title('{OSCC5} blood_polish result'.format(OSCC5 = sample_name))
    ComDf['reverse_polish'].value_counts().plot(kind='pie',autopct='%.2f%%')
    plt.savefig(work_dir + '/{OSCC5}_blood_polish_result.png'.format(OSCC5 = sample_name))
    plt.close()
    vcf_region_df_only_positive_polish  = ComDf[ComDf['reverse_polish'] == 'only_positive_polish']
    print('Start judge somatic SV Ture or False')
    P = Pool(60)
    RecordList_vcf = P.map(if_somatic_sv_by_sniffles_vcf, ParamList)
    P.close()
    del P
    print('Complete judge somatic SV Ture or False')
    JudgementDf_vcf = pd.DataFrame.from_records(RecordList_vcf)
    JudgementDf_vcf.columns = ['dir_name', 'somSV']
    JudgementDf_vcf.index = np.array(JudgementDf_vcf['dir_name'])
    ComDf_vcf = pd.concat([vcf_region_df_only_positive_polish, JudgementDf_vcf], axis=1)
    ComDf_vcf = ComDf_vcf.dropna(subset=['chr'])
    ComDf_vcf = ComDf_vcf.loc[:,~ComDf_vcf.columns.duplicated()]
    ComDf_vcf.to_csv(work_dir + '/{OSCC5}_vcf_region_df_only_positive_polish.csv'.format(OSCC5 = sample_name),index=False,header=True,sep='\t')
    print(work_dir + '/{OSCC5}_vcf_region_df_only_positive_polish.csv'.format(OSCC5 = sample_name))
    plt.title('{OSCC5} blood_only_positive_polish result'.format(OSCC5 = sample_name))
    plt.pie(ComDf_vcf['somSV'].value_counts(),labels=ComDf_vcf['somSV'].value_counts().index,autopct='%1.1f%%')
    plt.savefig(work_dir + '/{OSCC5}_blood_only_positive_polish_result.png'.format(OSCC5 = sample_name),bbox_inches='tight')
    plt.close()
    os.system('rm -rf {work_dir}/{OSCC5}_Sniffles2.*'.format(work_dir=work_dir,OSCC5 = sample_name)) # del filter






