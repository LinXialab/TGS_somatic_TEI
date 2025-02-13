import os
import os
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from collections import Counter
import scipy.stats as stats


samtools = "/biosoft/samtools-1.9/samtools"
bedtools = "/biosoft/bedtools-2.30.0"
sniffles2 = "/biosoft/sniffles2/bin/sniffles"
minimap2 = "/biosoft/minimap2-2.17/minimap2"


hg38_fa = "/nanopore/ref/hg38_mainChr.fa"
tr_bed = "/nanopore/ref/human_GRCh38_no_alt_analysis_set.trf.bed"
repeat = tr_bed
repeatmasker_dir = '/nanopore/ref/repeat'


work_dir_0 = "/nanopore/NewSample"
ms_dir = os.path.join(work_dir_0, '02.minimap2_sniffles')
coverage_dir = os.path.join(work_dir_0,'03.coverage')
somatic_dir = os.path.join(work_dir_0, "04.Somatic_minimap2_2.4")

workdir =  os.path.join(work_dir_0, '2_std_test')




def extract_region_bam(region, rnames, inbam, out):
    print('extract_region_bam')
    inbam_bai_path ='.'.join([inbam,'bai'])
    if not os.path.exists(inbam_bai_path):
        os.system('{samtools} index -@ 20 {inbam}'.format(samtools=samtools, inbam=inbam))
    if rnames:
        os.system("({samtools} view -H {inbam}; {samtools} view {inbam} {region}|grep -f {qnames})|"
                  "{samtools} view -b -@ 10 - -o {out}".format(
            samtools=samtools, inbam=inbam, region=region, qnames=rnames, out=out))
        if os.path.getsize(out) == 0:
            region1 = region.split(':')[0]+':'+str(int(region.split(':')[1].split('-')[0])-2500)+'-'+str(int(region.split(':')[1].split('-')[0])+2500)
            os.system("({samtools} view -H {inbam}; {samtools} view {inbam} {region}|grep -f {qnames})|"
                  "{samtools} view -b -@ 10 - -o {out}".format(
            samtools=samtools, inbam=inbam, region=region1, qnames=rnames, out=out))
    else:
        os.system("{samtools} view -@ 10 -h -b {inbam} {region} > {out}".format(
         samtools=samtools, inbam=inbam, region=region, out=out))
    os.system("{samtools} index {out}".format(samtools=samtools,out=out))
    print("{samtools} index {out}".format(samtools=samtools,out=out))


def big_del(svid, svlen, reads_list, region, tumorbam, index,intersect_info):
    print('big_del_wrong')
    print(svid)
    tumorsvBam= os.path.join(subdir, '_'.join([sampleid,svid,'tumor-sv.bam']))
    region1 = region[0] + ':' + str(int(region[1]) - 2500) + '-' + str(int(region[2]) + 2500)
    vr_txt = os.path.join(subdir, '%s_%s_variant_ont_reads.txt' % (sampleid,svid))
    with open(vr_txt, 'w') as out:
        for vr in reads_list:
            print(vr)
            if 'blood' in vr:
                print(vr)
            else:
                out.write(vr.split('tumor')[1] + '\n')
    if not os.path.exists(tumorsvBam):
        extract_region_bam(region1, vr_txt, tumorbam, tumorsvBam)
    if not os.path.exists(tumorsvBam+'.bai'):
        os.system("{samtools} index {tumorsvBam}".format(samtools=samtools, tumorsvBam=tumorsvBam))
    tumorsv_pysam = pysam.AlignmentFile(tumorsvBam, 'rb') #,ignore_truncation=True
    reads = tumorsv_pysam.fetch()
    if os.path.exists(os.path.join(subdir, '_'.join([sampleid,svid,'tumor_supportReads_del.bed']))):
        os.system('rm '+ os.path.join(subdir, '_'.join([sampleid,svid,'tumor_supportReads_del.bed'])))
    sv_dic = os.path.join(subdir, '_'.join([sampleid,svid,'tumor_supportReads_del.bed']))
    reads_alignment = {}
    for AlignedSegment in reads:
        query_n = AlignedSegment.query_name
        ref_n = AlignedSegment.reference_name
        if ref_n == region[0]:
            if query_n not in reads_alignment:
                reads_alignment[query_n] =[ref_n,AlignedSegment.get_blocks()[-1][1]]
            # else:
            #     reads_alignment[query_n] = [ref_n, AlignedSegment.get_blocks()[-1][1]]
                reads_alignment[query_n].append(AlignedSegment.get_blocks()[0][0])
    reads_alignment1 = pd.DataFrame.from_dict(reads_alignment,orient='index')
    print('1')
    print(reads_alignment1)
    svlen_list = np.abs(reads_alignment1[2]-reads_alignment1[1])
    print('2')
    print(svlen_list)
    reads_alignment1['svlen']=svlen_list
    print('3')
    print('reads_alignment1')
    reads_alignment1.to_csv(sv_dic,sep='\t',index=None,header=None)
    with open(os.path.join(subdir,sampleid+'_'+str(svid)+'_test.bed'),'w') as out:
            out.write('\t'.join([region[0],str(int(region[1])-100), str(int(region[2]+100))])+'\n')
    intersect_info = os.popen(bedtools+' intersect -a %s -b %s -f 0.4 -wa' % (os.path.join(subdir,
                                                                                           '_'.join([sampleid,svid,'tumor_supportReads_del.bed'])),
                                                                                  os.path.join(subdir,sampleid+'_'+str(svid)+'_test.bed'))).read().strip().split('\n')
    svlen_list = []
    if intersect_info[0]:
        for info in intersect_info:
            svlen_list.append(info.split('\t')[-1])
        svlen_list = list(map(lambda x: int(x),[y for y in svlen_list]))
        svlen_std = np.std(svlen_list)
        Q1=np.percentile(svlen_list, 25)
        Q3=np.percentile(svlen_list, 75)
        IQR=Q3-Q1
        outlier_step = IQR*1.5
        if index:
            if svlen > np.mean(svlen_list)+2*svlen_std or svlen< np.mean(svlen_list)-2*svlen_std:
                std_test = 'abnormal_bk'
            else:
                std_test = 'normal_bk'
            if svlen<(Q1-outlier_step) or svlen>(Q3+outlier_step):
                iqr_test = 'abnormal_bk'
            else:
                iqr_test='normal_bk'
            return([std_test,iqr_test])
        else:
            return(svlen_list)
    else:
        if index:
            return(['abnormal_bk','abnormal_bk'])
        else:
            return('Null')


# In[262]:


def check_whether_abnormalLen_deletion(tumorsv_bam,svid, sampleid, region,svlen,subdir, reads_list, index):
    inbam_bai_path = '.'.join([tumorsv_bam, 'bai'])
    if not os.path.exists(inbam_bai_path):
        os.system("{samtools} index {out}".format(samtools=samtools, out=tumorsv_bam))
    print("check_whether_abnormalLen_deletion")
    print(tumorsv_bam)
    tumorsv_pysam = pysam.AlignmentFile(tumorsv_bam, 'rb') #, ignore_truncation=True
    reads = tumorsv_pysam.fetch()
    if os.path.exists(os.path.join(subdir, '_'.join([sampleid,svid,'tumor_supportReads_del.bed']))):
        os.system('rm '+ os.path.join(subdir, '_'.join([sampleid,svid,'tumor_supportReads_del.bed'])))
    sv_dic = open(os.path.join(subdir, '_'.join([sampleid,svid,'tumor_supportReads_del.bed'])),'w')
    for AlignedSegment in reads:
        query_n = AlignedSegment.query_name
        ref_n = AlignedSegment.reference_name
        if ref_n == region[0]:
            cigartuples = AlignedSegment.cigartuples
            aligned_pairsv = AlignedSegment.get_aligned_pairs(matches_only=False, with_seq=False)
            cigar_loc = 0
            query_loc = 0
            for cigar in cigartuples:
                if cigar[0] != 2:
                    query_loc += cigar[1]
                if cigar[0] != 5:  # supplementary is hard clip
                    cigar_loc += cigar[1]
                if cigar[0] == 2 and cigar[1] >= 30:  # DEL
                    query_sv_end = query_loc
                    query_sv_start = query_loc - 1
                    align_sv_end = cigar_loc
                    align_sv_start = cigar_loc - 1
                    sv_len = cigar[1]
                    sv_ref_start = aligned_pairsv[align_sv_start][1] + 1
                    sv_ref_end = aligned_pairsv[align_sv_end-1][1] + 1
                    if AlignedSegment.is_reverse:
                        qlen = AlignedSegment.query_length
                        query_sv_start_r = qlen - query_sv_end
                        query_sv_end_r = qlen - query_sv_start
                        sv_dic.write('\t'.join(map(lambda x: str(x), [y for y in [ref_n, sv_ref_start,sv_ref_end,
                                                                                  sv_len, query_n + '_reverse', query_sv_start_r, query_sv_end_r,
                                                                                  align_sv_start, align_sv_end]]))+'\n')
                    else:
                        sv_dic.write('\t'.join(map(lambda x: str(x), [y for y in[ref_n, sv_ref_start, sv_ref_end,
                                                                                 sv_len, query_n, query_sv_start, query_sv_end,
                                                                                 align_sv_start, align_sv_end]]))+'\n')
    tumorsv_pysam.close()
    sv_dic.close()
    with open(os.path.join(subdir,sampleid+'_'+str(svid)+'_test.bed'),'w') as out:
        out.write('\t'.join(map(lambda x:str(x),[y for y in region]))+'\n')
    intersect_info = os.popen(bedtools+' intersect -a %s -b %s -f 0.4 -wa' % (os.path.join(subdir, '_'.join([sampleid,svid,'tumor_supportReads_del.bed'])),
                                                                              os.path.join(subdir,sampleid+'_'+str(svid)+'_test.bed'))).read().strip().split('\n')
    svlen_list = []
    print(intersect_info[0])
    if intersect_info[0]:
        for info in intersect_info:
            svlen_list.append(info.split('\t')[3])
        svlen_list = list(map(lambda x: int(x),[y for y in svlen_list]))
        svlen_std = np.std(svlen_list)
        Q1=np.percentile(svlen_list, 25)
        Q3=np.percentile(svlen_list, 75)
        IQR=Q3-Q1
        outlier_step = IQR*1.5
        if index:
            if svlen > np.mean(svlen_list)+2*svlen_std or svlen< np.mean(svlen_list)-2*svlen_std:
                std_test = 'abnormal_bk'
            else:
                std_test = 'normal_bk'
            if svlen<(Q1-outlier_step) or svlen>(Q3+outlier_step):
                iqr_test = 'abnormal_bk'
            else:
                iqr_test='normal_bk'
            return([std_test,iqr_test])
        else:
            return(svlen_list)
    else:
        return(big_del(svid, svlen, reads_list, region, tumorbam, index,intersect_info))





def check_whether_abnormalLen_insertion(tumorsv_bam,svid, sampleid, region,svlen,subdir,index):
    inbam_bai_path = '.'.join([tumorsv_bam, 'bai'])
    if not os.path.exists(inbam_bai_path):
        os.system("{samtools} index {out}".format(samtools=samtools, out=tumorsv_bam))
    print("check_whether_abnormalLen_insertion")
    print(tumorsv_bam)
    tumorsv_pysam = pysam.AlignmentFile(tumorsv_bam, 'rb') #, ignore_truncation=True
    reads = tumorsv_pysam.fetch()
    if os.path.exists(os.path.join(subdir, '_'.join([sampleid,svid,'tumor_supportReads_ins.bed']))):
        os.system('rm '+ os.path.join(subdir, '_'.join([sampleid,svid,'tumor_supportReads_ins.bed'])))
    sv_dic = open(os.path.join(subdir, '_'.join([sampleid,svid,'tumor_supportReads_ins.bed'])),'w')
    for AlignedSegment in reads:
        query_n = AlignedSegment.query_name
        ref_n = AlignedSegment.reference_name
        if ref_n == region[0]:
            cigartuples = AlignedSegment.cigartuples
            aligned_pairsv = AlignedSegment.get_aligned_pairs(matches_only=False, with_seq=False)
            cigar_loc = 0
            query_loc = 0
            for cigar in cigartuples:
                if cigar[0] != 2:
                    query_loc += cigar[1]
                if cigar[0] != 5:  # supplementary is hard clip
                    cigar_loc += cigar[1]
                if cigar[0] == 1 and cigar[1] >= 30:  # INS
                    query_sv_end = query_loc
                    query_sv_start = query_loc - cigar[1]
                    align_sv_end = cigar_loc
                    align_sv_start = cigar_loc - cigar[1]
                    sv_len = cigar[1]
                    sv_ref_start = aligned_pairsv[align_sv_start - 1][1] + 1
                    sv_ref_end = aligned_pairsv[align_sv_end][1] + 1
                    if AlignedSegment.is_reverse:
                        qlen = AlignedSegment.query_length
                        query_sv_start_r = qlen - query_sv_end
                        query_sv_end_r = qlen - query_sv_start
                        sv_dic.write('\t'.join(map(lambda x: str(x), [y for y in [ref_n, sv_ref_start,sv_ref_end,
                                                                                  sv_len, query_n + '_reverse', query_sv_start_r, query_sv_end_r,
                                                                                  align_sv_start, align_sv_end]]))+'\n')
                    else:
                        sv_dic.write('\t'.join(map(lambda x: str(x), [y for y in[ref_n, sv_ref_start, sv_ref_end,
                                                                                 sv_len, query_n, query_sv_start, query_sv_end,
                                                                                 align_sv_start, align_sv_end]]))+'\n')
    tumorsv_pysam.close()
    sv_dic.close()
    with open(os.path.join(subdir,sampleid+'_'+str(svid)+'_test.bed'),'w') as out:
        out.write('\t'.join([region[0],str(int(region[1])-100), str(int(region[2]+100))])+'\n')
    intersect_info = os.popen(bedtools+' intersect -a %s -b %s -f 0.4 -wa' % (os.path.join(subdir, '_'.join([sampleid,svid,'tumor_supportReads_ins.bed'])),
                                                                              os.path.join(subdir,sampleid+'_'+str(svid)+'_test.bed'))).read().strip().split('\n')
    svlen_list = []
    if intersect_info[0]:
        for info in intersect_info:
            svlen_list.append(info.split('\t')[3])
        svlen_list = list(map(lambda x: int(x),[y for y in svlen_list]))
        svlen_std = np.std(svlen_list)
        Q1=np.percentile(svlen_list, 25)
        Q3=np.percentile(svlen_list, 75)
        IQR=Q3-Q1
        outlier_step = IQR*1.5
        if index:
            if svlen > np.mean(svlen_list)+2*svlen_std or svlen< np.mean(svlen_list)-2*svlen_std:
                std_test = 'abnormal_bk'
            else:
                std_test = 'normal_bk'
            if svlen<(Q1-outlier_step) or svlen>(Q3+outlier_step):
                iqr_test = 'abnormal_bk'
            else:
                iqr_test='normal_bk'
            return([std_test,iqr_test])
        else:
            return(svlen_list)
    else:
        if index:
            return(['abnormal_bk','abnormal_bk'])
        else:
            return('Null')






all_file_list = ['Sample1']

for sampleid  in all_file_list:
        subdir = os.path.join(workdir, sampleid)
        if not os.path.exists(subdir):
            os.makedirs(subdir)
        os.chdir(subdir)
        print(subdir)
        B_ID = sampleid.split("C")[0]
        mergebed = "{ms_dir}/{sample}/sniffles2/{sample}_merge_minimap2_sniffles_v2.bed".format(ms_dir=ms_dir,
            sample=sampleid)
        tumorbed = "{ms_dir}/{sample}/sniffles2/{sample}_merge_minimap2_sniffles_v2_tumor_somatic.bed".format(ms_dir=ms_dir,
            sample=sampleid)
        tumorbam = "{somatic_dir}/Cancer_cells/{sample}/{sample}_Cancer_cells_minimap2_sorted_tag.bam".format(somatic_dir=somatic_dir,
            sample=sampleid)
        bloodbam = "{somatic_dir}/blood/{sample}/{sample}_blood_minimap2_sorted_tag.bam".format(somatic_dir=somatic_dir,
            sample=B_ID)
        vcf = "{ms_dir}/{sample}/sniffles2/{sample}_merge_minimap2_sniffles_v2.vcf".format(ms_dir=ms_dir,
            sample=sampleid)
        usecols1 = [0, 1, 2, 6, 7, 8, 10, 12, 13, 14, 15, 16]
        df = pd.read_csv(mergebed, sep='\t', header=None)
        sample_state = []
        for i in list(df.index):
            if 'blood' in df.loc[i, 10]:
                if 'tumor' in df.loc[i, 10]:
                    sample_state.append('TB')
                else:
                    sample_state.append('B')
            else:
                sample_state.append('T')
        df['sample_state'] = sample_state
        somaticsv_df = df[df['sample_state'] == 'T']
        sv_list = list(somaticsv_df[7])
        newsv_list = {}
        with open(vcf) as vfile:
            for line in vfile.readlines():
                if '#' not in line:
                    comment = line.split('\t')
                    # print(comment[2])
                    if comment[2] in sv_list:
                        # print(comment[2])
                        for com in comment:
                            if 'STDEV_POS' in com:
                                com1 = comment
                                std = line.split('STDEV_POS=')[1].split(';')[0].split('\t')[0]
                                if float(std) >= 150:
                                    newsv_list[comment[2]] = 'abnormal_std'
                                else:
                                    newsv_list[comment[2]] = 'normal_std'
        newsv_list1 = pd.DataFrame.from_dict(newsv_list, orient='index')
        newsv_list2 = newsv_list1.sort_index()
        somaticsv_df1 = somaticsv_df.sort_values(by=[7])
        std_state = []
        for i in range(0, somaticsv_df1.shape[0]):
            std_state.append(newsv_list2.loc[somaticsv_df1.iloc[i, 7], 0])
        somaticsv_df1['stdpos'] = std_state
        somaticsv_df1 = somaticsv_df1[somaticsv_df1[8] < 1000000]
        print(somaticsv_df1.shape[0])
        repeat_overlapped_bed = os.path.join(subdir, '_'.join([sampleid, 'merge', 'sv', 'repeatmasker.bed']))
        os.system('cat %s/*_repeat.bed |' % repeatmasker_dir + bedtools + ' intersect -a %s -b - -wa -wb -f 0.5 > %s' % (
        mergebed, repeat_overlapped_bed))
        repeat_df = pd.read_csv(repeat_overlapped_bed, sep='\t', header=None,
                                usecols=usecols1,
                                names=['chr', 'start', 'end', 'svtype', 'svid', 'svlen', 'svreads', 'r_chr', 'r_start',
                                       'r_end', 'pattern', 'repeat_type'])
        repeat_df = repeat_df[repeat_df['svlen'] < 1000000]
        sample_state = []
        for i in list(repeat_df.index):
            if 'blood' in repeat_df.loc[i, 'svreads']:
                if 'tumor' in repeat_df.loc[i, 'svreads']:
                    sample_state.append('TB')
                else:
                    sample_state.append('B')
            else:
                sample_state.append('T')
        repeat_df['index'] = repeat_df['r_chr'] + '-' + repeat_df['r_start'].map(str) + '-' + repeat_df['r_end'].map(
            str) + '-' + repeat_df['repeat_type']
        repeat_df['sample_state'] = sample_state
        repeat_df.head()
        repeat_df1 = repeat_df[repeat_df['svtype'].isin(['INS', 'DEL'])]
        repeat_t_b = {}
        for repeat in list(set(repeat_df1['index'])):
            subdf = repeat_df1[repeat_df1['index'] == repeat]
            svtypes = dict(Counter(subdf['svtype']))
            for sv in svtypes:
                subdat = subdf[subdf['svtype'] == sv]
                if subdat.shape[0] > 1:
                    reads = subdat.iloc[0].to_list()[6].split(',')
                    sample_states = list(set(subdat['sample_state']))
                    print(sample_states)
                    if 'T' in sample_states:
                        if len(set(sample_states) - {'T'}) > 0:
                            # print(subdat)
                            tumor_specifc = subdat[subdat['sample_state'] == 'T']
                            region = [subdat.loc[list(subdat.index)[0], 'chr'], min(subdat['start']), max(subdat['end'])]
                            for i in list(tumor_specifc.index):
                                if tumor_specifc.loc[i, 'svid'] not in repeat_t_b:
                                    sub_tumor_specific_svreads = tumor_specifc.loc[i, 'svreads'].split(',')
                                    sub_svid = tumor_specifc.loc[i, 'svid']
                                    sub_svlen = tumor_specifc.loc[i, 'svlen']
                                    vr_txt = os.path.join(subdir, '%s_%s_variant_ont_reads.txt' % (sampleid, sub_svid))
                                    with open(vr_txt, 'w') as out:
                                        for vr in sub_tumor_specific_svreads:
                                            out.write(vr.split('tumor')[1] + '\n')
                                    region1 = region[0] + ':' + str(int(region[1]) - 2500) + '-' + str(
                                        int(region[2]) + 2500)
                                    tumorBam = os.path.join(subdir, '_'.join([sampleid, sub_svid, 'tumor_ALL.bam']))
                                    tumorsvBam = os.path.join(subdir, '_'.join([sampleid, sub_svid, 'tumor-sv.bam']))
                                    bloodBam = os.path.join(subdir, '_'.join([sampleid, sub_svid, 'blood.bam']))
                                    if not os.path.exists(tumorsvBam) or os.path.getsize(tumorsvBam) < 1000:
                                        extract_region_bam(region1, vr_txt, tumorbam, tumorsvBam)
                                    if not os.path.exists(tumorBam) or os.path.getsize(tumorBam) < 1000:
                                        extract_region_bam(region1, '', tumorbam, tumorBam)
                                    if not os.path.exists(bloodBam) or os.path.getsize(bloodBam) < 1000:
                                        print('bloodBam_extract_region_bam')
                                        extract_region_bam(region1, '', bloodbam, bloodBam)
                                        extract_region_bam(region1, '', bloodbam, bloodBam)
                                    if not os.path.exists(bloodBam) or os.path.getsize(bloodBam) < 1000:
                                        print('bloodBam_extract_region_bam')
                                        extract_region_bam(region1, '', bloodbam, bloodBam)
                                        extract_region_bam(region1, '', bloodbam, bloodBam)
                                    if sv == 'INS':
                                        t_svlen_list = check_whether_abnormalLen_insertion(tumorsvBam, sub_svid, sampleid,
                                                                                           region, sub_svlen, subdir, False)
                                        b_svlen_list = check_whether_abnormalLen_insertion(bloodBam, sub_svid, sampleid,
                                                                                           region, sub_svlen, subdir, False)
                                    else:
                                        t_svlen_list = check_whether_abnormalLen_deletion(tumorsvBam, sub_svid, sampleid,
                                                                                          region, sub_svlen, subdir, reads,
                                                                                          False)
                                        b_svlen_list = check_whether_abnormalLen_deletion(bloodBam, sub_svid, sampleid,
                                                                                          region, sub_svlen, subdir, reads,
                                                                                          False)
                                    if 'Null' == b_svlen_list:
                                        p = 0.0000000000005
                                    elif 'Null' == t_svlen_list:
                                        print('error ' + sub_svid)
                                    else:
                                        t1 = len(set(t_svlen_list))
                                        t2 = len(set(b_svlen_list))
                                        if (t1 == t2 == 1) and (set(t_svlen_list) == set(b_svlen_list)):
                                            p = 0.0000000000005
                                        else:
                                            stat, p = stats.mannwhitneyu(t_svlen_list, b_svlen_list,
                                                                         alternative='two-sided')
                                        # stat,p = stats.mannwhitneyu(t_svlen_list, b_svlen_list, alternative='two-sided')
                                    if p < 0.05:
                                        repeat_t_b[tumor_specifc.loc[i, 'svid']] = 'repeat_dif'
                                    else:
                                        repeat_t_b[tumor_specifc.loc[i, 'svid']] = 'repeat_same'
                        else:
                            repeat_t_b[subdat.iloc[0]['svid']] = 'repeat'
                            # repeat_t_b[tumor_specifc.loc[i, 'svid']] = 'repeat'
                else:
                    if 'T' in subdat['sample_state']:
                        repeat_t_b[subdat.iloc[0]['svid']] = 'repeat'
        new_repeat_df = repeat_df[repeat_df['sample_state'] == 'T']
        new_repeat_df1 = new_repeat_df.sort_values(by=['svid'])
        svlen_state = []
        for i in list(new_repeat_df1.index):
            svid = new_repeat_df1.loc[i, 'svid']
            if svid in repeat_t_b:
                svlen_state.append(repeat_t_b[svid])
            else:
                svlen_state.append('repeat')
        new_repeat_df1['svlen_state'] = svlen_state
        Counter(new_repeat_df1['svlen_state'])
        somaticsv_df1 = somaticsv_df1[[0, 1, 2, 6, 7, 8, 10, 'stdpos']]
        somaticsv_df1.columns = ['chr', 'start', 'end', 'svtype', 'svid', 'svlen', 'svreads', 'stdpos']
        somaticsv_df1 = somaticsv_df1.sort_values(by=['svid'])
        svlen_state = []
        for i in list(somaticsv_df1['svid']):
            if i in list(new_repeat_df1['svid']):
                svlen_state.append(new_repeat_df1[new_repeat_df1['svid'] == i]['svlen_state'].to_list()[0])
            else:
                svlen_state.append('None')
        somaticsv_df1['repeat_svlen_state'] = svlen_state
        bk_state = {}
        iqr_state = {}
        n = 0
        for i in list(somaticsv_df1.index):
            print(n)
            svtype = somaticsv_df1.loc[i, 'svtype']
            svid = somaticsv_df1.loc[i, 'svid']
            if svtype in ['INS', 'DEL']:
                svlen = somaticsv_df1.loc[i, 'svlen']
                tumorsvBam = os.path.join(subdir, '_'.join([sampleid, svid, 'tumor-sv.bam']))
                region = somaticsv_df1.loc[i, ['chr', 'start', 'end']]
                region1 = region[0] + ':' + str(int(region[1]) - 2500) + '-' + str(int(region[2]) + 2500)
                reads = somaticsv_df1.loc[i, 'svreads'].split(',')
                vr_txt = os.path.join(subdir, '%s_%s_variant_ont_reads.txt' % (sampleid, svid))
                with open(vr_txt, 'w') as out:
                    for vr in reads:
                        out.write(vr.split('tumor')[1] + '\n')
                if not os.path.exists(tumorsvBam):
                    extract_region_bam(region1, vr_txt, tumorbam, tumorsvBam)
                if not os.path.exists(tumorsvBam + '.bai'):
                    os.system('samtools index %s' % tumorsvBam)
                if svtype == 'INS':
                    # print(svid)
                    bkState = check_whether_abnormalLen_insertion(tumorsvBam, svid, sampleid, region.to_list(), svlen,
                                                                  subdir, True)
                    bk_state[svid] = bkState[0]
                    iqr_state[svid] = bkState[1]
                else:
                    # print(svid)
                    bkState = check_whether_abnormalLen_deletion(tumorsvBam, svid, sampleid, region.to_list(), svlen,
                                                                 subdir, reads, True)
                    bk_state[svid] = bkState[0]
                    iqr_state[svid] = bkState[1]
            else:
                bk_state[svid] = 'None'
                iqr_state[svid] = 'None'
                # print(svid)
            n += 1
        bk_state1 = []
        for svid in list(somaticsv_df1['svid']):
            bk_state1.append(bk_state[svid])
        iqr_state1 = []
        for svid in list(somaticsv_df1['svid']):
            iqr_state1.append(iqr_state[svid])
        somaticsv_df1['bk_state'] = bk_state1
        somaticsv_df1['iqr_state'] = iqr_state1
        somaticsv_df1.to_csv(os.path.join(subdir, '_'.join([sampleid, 'somatic_sv_state-20221229.csv'])), sep=',', index=0)
        data = pd.read_csv(os.path.join(subdir, '_'.join([sampleid, 'somatic_sv_state-20221229.csv']))) #'/NAS/wg_xialin/pancancer/check_somaticsv/GBM4/test/GBM4_somatic_sv_state-20221229.csv'
        data1 = data[data['stdpos'] == 'normal_std']
        data2 = data1[data1['repeat_svlen_state'] != 'repeat_same']
        nonrepeat = data2[data2['repeat_svlen_state'] == 'None']
        repeat = data2[data2['repeat_svlen_state'] != 'None']
        nonrepeat_filter = nonrepeat[(nonrepeat['bk_state'] != 'abnormal_bk') | (nonrepeat['iqr_state'] != 'abnormal_bk')]
        final_abnormal = nonrepeat_filter.append(repeat)
        final_abnormal.to_csv(os.path.join(subdir, '_'.join([sampleid, 'somatic_sv_state-20221230-final.csv'])), sep=',',
                              index=0)



