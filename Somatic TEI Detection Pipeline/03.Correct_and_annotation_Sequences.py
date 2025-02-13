import pandas as pd
import os
import pysam
from multiprocessing import Pool


racon = "/biosoft/racon/build/bin/racon"
samtools = "/biosoft/samtools-1.9/samtools"   ##1.9
minimap2 = "/biosoft/minimap2-2.17/minimap2"    ##2.17
sniffles = "/biosoft/sniffles2/bin/sniffles"
seqkit = '/biosoft/seqkit'
bedtools = "/biosoft//bedtools-2.30.0"


hg38_fa = "/nanopore/ref/hg38_mainChr.fa"
tr_bed = "/nanopore/ref/human_GRCh38_no_alt_analysis_set.trf.bed"


work_dir = "/nanopore/NewSample"

ms_dir = os.path.join(work_dir, '02.minimap2_sniffles')
coverage_dir = os.path.join(work_dir,'03.coverage')
somatic_dir = os.path.join(work_dir, "04.Somatic_minimap2_2.4")


workdir_0 = os.path.join(work_dir, 'Iris_test')


############################################## Correct Sequences ##############################################

def variant_ont_qnames(vr_txt_x, ont_bam_x, region_x, ont_bam_region_x, vr_bam_x, vr_fq_x,vr_fq_x_2):
    # candidate region bam
    os.system("{samtools} view -@ 10 -h -b {bam_ont} {region} > {bam_out}".format(
        samtools=samtools, bam_ont=ont_bam_x, region=region_x, bam_out=ont_bam_region_x))
    # variant bam
    os.system("({samtools} view -H {inbam}; {samtools} view {inbam}|grep -f {qnames})|"
              "{samtools} view -b - -o {out}".format(
                samtools=samtools, inbam=ont_bam_region_x, qnames=vr_txt_x, out=vr_bam_x))
    os.system('{samtools} index {bam}'.format(samtools=samtools, bam=vr_bam_x))
    # # bam to fastq, forward
    os.system(''' %s view -F 2048 %s | awk 'BEGIN {FS="\\t"} {print "@"$1"\\n"$10"\\n+\\n"$11}' - > %s ''' % (
        samtools, vr_bam_x, vr_fq_x))
    os.system(''' %s view -F 2048 %s | awk 'BEGIN {FS="\\t"} {print "@"$1"\\t"$5"\\t%s\\t"$10"\\n"}' - > %s ''' % (
        samtools, vr_bam_x,region_x,vr_fq_x_2))


def find_fq_MQ_max_reads(vr_fq_x):
    df = pd.read_csv(vr_fq_x, sep='\t', names=["reads_name", "MQ", "breakpoint", "sequence"])
    col = "MQ"
    max_x = df.loc[df[col].astype(float).idxmax()]
    df_2 = pd.DataFrame(columns=["reads_name", "MQ", "breakpoint", "sequence"])
    for idx, row in df.iterrows():
        if row.MQ >= max_x.MQ:
            df_2 = df_2.append(row)
    reads_name_list = []
    for idx, row in df_2.iterrows():
        reads_name_list.append(row.reads_name.split('@')[1])
    return reads_name_list

def DUP_pysam(region_file,reads_name,end,star,out_fa,out_only_DUP_fa,sv_region):
    n = 0
    query_star = 0
    reads_query_star_end = {}
    reads_name_list2 = []
    reads_name_list3 = []
    for AlignedSegmennt in region_file:
        n = n + 1
        aligned_pairsv = AlignedSegmennt.get_aligned_pairs(matches_only=False, with_seq=False)
        for pair1 in aligned_pairsv:
            query_star1_2 = 0
            if pair1[1] and pair1[0]:
                if pair1[1] == int(star):
                    query_star1_2 = pair1[0]
                    # print(query_star1_2)
                    if (AlignedSegmennt.query_length - int(query_star1_2)) > (int(end) - int(star)):
                        reads_name_list3.append(reads_name) ##across DUP region
                    break
        if query_star1_2 == 0:
            query_star = 0
        else:
            if query_star < int(query_star1_2):
                query_star = query_star1_2
    fa_df = pd.read_csv(out_fa, sep="\t", names=['reads_name', 'sequence'])
    for index, row in fa_df.iterrows():
        if len(row.sequence) > 3:  ## over 3bp
            line1 = '>%s' % reads_name + '\n'
            line2 = row.sequence[int(query_star):int(query_star) + int(sv_region)]
            newline = line1 + line2
            with open(out_only_DUP_fa, 'w') as DUP_fa:
                DUP_fa.write(newline)  # DUP in fasta
            reads_name_list2.append(reads_name)
            reads_query_star_end[reads_name] = [int(query_star), int(query_star) + int(sv_region)]
        else:
            reads_name_list2 = []
            # print(fa_df)
            # print(reads_name)
    return reads_query_star_end, reads_name_list3, reads_name_list2


def INS_pysam(region_file,reads_name,out_fa,out_only_DUP_fa,svlen):
    print('star to  INS_pysam')
    query_star = 0
    query_end = 0
    reads_query_star_end = {}
    reads_name_list2 = []
    reads_name_list3 = []
    for AlignedSegmennt in region_file:
        cigartuples = AlignedSegmennt.cigartuples
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
                # if AlignedSegmennt.is_reverse:
                query_star = query_sv_start
                query_end = query_sv_end
    fa_df = pd.read_csv(out_fa, sep="\t", names=['reads_name', 'sequence'])
    for index, row in fa_df.iterrows():
        if len(row.sequence) > 3:
            line1 = '>%s' % reads_name + '\n'
            line2 = row.sequence[int(query_star):int(query_end)]
            newline = line1 + line2
            with open(out_only_DUP_fa, 'w') as DUP_fa:
                DUP_fa.write(newline)
            reads_name_list2.append(reads_name)
            reads_name_list3.append(reads_name)
            reads_query_star_end[reads_name] = [int(query_star), int(query_end)]
        else:
            reads_name_list2 = []
    print(reads_query_star_end, reads_name_list3, reads_name_list2)
    return reads_query_star_end, reads_name_list3, reads_name_list2


def greads_query_star_end(reads_name_list,vr_dir,svloc,star,end,vr_bam,svlen):
    print(svloc)
    print('star to greads_query_star_end')
    for reads_name in reads_name_list:
        out_bam = os.path.join(vr_dir, '%s_%s_variant_ont_reads.bam' % (svloc, reads_name))
        out_fa = os.path.join(vr_dir, '%s_%s_variant_ont_reads.fastq' % (svloc, reads_name))
        out_only_DUP_fa = os.path.join(vr_dir, '%s_%s_DUP_only_reads.fasta' % (svloc, reads_name))
        os.system("({samtools} view -H {inbam}; {samtools} view {inbam}|grep  {qnames})|"
                  "{samtools} view -b - -o {out}".format(
            samtools=samtools, inbam=vr_bam, qnames=reads_name, out=out_bam)) # in vr_bam get
        os.system(''' %s view -F 2048  %s | awk 'BEGIN {FS="\\t"} {print "@"$1"\\t"$10"\\n"}' - > %s ''' % (
            samtools, out_bam,  out_fa))
        os.system("{samtools} index {out}".format(samtools=samtools, out=out_bam))
        bam_file = pysam.AlignmentFile(out_bam, 'rb')
        region_file = bam_file.fetch()
        sv_region = int(end) - int(star)
        if svloc.split('.')[-2] == 'DUP':
            reads_query_star_end_x, reads_name_list3_x, reads_name_list2_x = DUP_pysam(region_file, reads_name, end, star, out_fa, out_only_DUP_fa,sv_region)
            return reads_query_star_end_x, reads_name_list3_x, reads_name_list2_x
        elif svloc.split('.')[-2] == 'INS':
            print('INS')
            reads_query_star_end_x, reads_name_list3_x, reads_name_list2_x = INS_pysam(region_file, reads_name,out_fa, out_only_DUP_fa,svlen)
            return reads_query_star_end_x, reads_name_list3_x, reads_name_list2_x



def get_line_dup(sample_name, chr, star, sv_type, end, reads, svloc,svlen):
    bam_merge = get_merge_bam_path(sample_name)
    workdir = os.path.join(workdir_0,sample_name)
    assembly_dir = os.path.join(workdir, 'assembly_DUP', '%s' % sample_name)
    if not os.path.exists(assembly_dir):
        os.makedirs(assembly_dir)
    sv_dir = os.path.join(assembly_dir, svloc)
    if not os.path.exists(sv_dir):
        os.mkdir(sv_dir)
    vr_dir = os.path.join(sv_dir, 'variant_reads')
    if not os.path.exists(vr_dir):
        os.mkdir(vr_dir)
    vr_txt = os.path.join(vr_dir, '%s_variant_ont_reads.txt' % svloc)
    with open(vr_txt, 'w') as out:
        for vr in reads:
            if 'tumor' in vr:
                out.write(vr + '\n')
    # region bam
    region = "%s:%s-%s" % (chr, int(star) - 1000, int(end) + 1000)
    ont_bam_region = os.path.join(vr_dir, '%s_flk1000_region.bam' % svloc)
    vr_bam = os.path.join(vr_dir, '%s_variant_ont_reads.bam' % svloc)
    vr_fq = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)
    vr_fq_2 = os.path.join(vr_dir,
                           '%s_variant_ont_reads_for_Iris.fastq' % svloc)
    variant_ont_qnames(vr_txt, bam_merge, region, ont_bam_region, vr_bam, vr_fq,
                       vr_fq_2)
    if os.path.getsize(vr_fq_2) > 0:
        reads_name_list = find_fq_MQ_max_reads(vr_fq_2)
    else:
        ont_bam_region = os.path.join(vr_dir, '%s_flk10000_region.bam' % svloc)
        region = "%s:%s-%s" % (chr, int(star) - 10000, int(end) + 10000)
        variant_ont_qnames(vr_txt, bam_merge, region, ont_bam_region, vr_bam, vr_fq,
                           vr_fq_2)
        if os.path.getsize(vr_fq_2) > 0:
            reads_name_list = find_fq_MQ_max_reads(vr_fq_2)
        else:
            return 0
    greads_query_star_end_list, across_DUP_reads_name_list, across_DUP_reads_name_list2 = greads_query_star_end(reads_name_list, vr_dir, svloc, star, end, vr_bam,svlen)
    across_DUP_reads_name_list3 = set(across_DUP_reads_name_list).intersection(across_DUP_reads_name_list2)
    if not across_DUP_reads_name_list3: # across_DUP_reads_name_list3 is empty
        # print(svloc)
        line_dup = '*'
        return line_dup
    else:
        out_only_DUP_fa = os.path.join(vr_dir, '%s_%s_DUP_only_reads.fasta' % (
            svloc, list(across_DUP_reads_name_list3)[0]))
        for_racon_fa = os.path.join(vr_dir, '%s_for_racon_Assembly.fasta' % (
            svloc))  ## first across reads
        # print(for_racon_fa)
        with open(out_only_DUP_fa, 'r') as only_DUP_fa:
            for line in only_DUP_fa:
                if line[0] != '>':
                    line_dup = line.strip()
                    return line_dup



def get_line_dup_all_reads(sample_name, chr, star, sv_type, end, reads, svloc,svlen):
    bam_merge = get_merge_bam_path(sample_name)
    workdir = os.path.join(workdir_0,sample_name)
    assembly_dir = os.path.join(workdir, 'assembly_DUP', '%s' % sample_name)
    if not os.path.exists(assembly_dir):
        os.makedirs(assembly_dir)
    # print(line)
    sv_dir = os.path.join(assembly_dir, svloc)
    if not os.path.exists(sv_dir):
        os.mkdir(sv_dir)
    vr_dir = os.path.join(sv_dir, 'variant_reads')
    if not os.path.exists(vr_dir):
        os.mkdir(vr_dir)
    vr_txt = os.path.join(vr_dir, '%s_variant_ont_reads.txt' % svloc)
    with open(vr_txt, 'w') as out:
        for vr in reads:
            out.write(vr + '\n')
    # region bam
    region = "%s:%s-%s" % (chr, int(star) - 1000, int(end) + 1000)
    ont_bam_region = os.path.join(vr_dir, '%s_flk1000_region.bam' % svloc)
    vr_bam = os.path.join(vr_dir, '%s_variant_ont_reads.bam' % svloc)
    vr_fq = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)
    vr_fq_2 = os.path.join(vr_dir,
                           '%s_variant_ont_reads_for_Iris.fastq' % svloc)
    variant_ont_qnames(vr_txt, bam_merge, region, ont_bam_region, vr_bam, vr_fq,
                       vr_fq_2)
    if os.path.getsize(vr_fq_2) > 0:
        reads_name_list = find_fq_MQ_max_reads(vr_fq_2)
    else:
        ont_bam_region = os.path.join(vr_dir, '%s_flk10000_region.bam' % svloc)
        region = "%s:%s-%s" % (chr, int(star) - 10000, int(end) + 10000)
        variant_ont_qnames(vr_txt, bam_merge, region, ont_bam_region, vr_bam, vr_fq,
                           vr_fq_2)
        if os.path.getsize(vr_fq_2) > 0:
            reads_name_list = find_fq_MQ_max_reads(vr_fq_2)
        else:
            return 0
    greads_query_star_end_list, across_DUP_reads_name_list, across_DUP_reads_name_list2 = greads_query_star_end(reads_name_list, vr_dir, svloc, star, end, vr_bam,svlen)
    across_DUP_reads_name_list3 = set(across_DUP_reads_name_list).intersection(across_DUP_reads_name_list2)
    if not across_DUP_reads_name_list3: # across_DUP_reads_name_list3 is empty
        # print(svloc)
        line_dup = '*'
        return line_dup
    else:
        out_only_DUP_fa = os.path.join(vr_dir, '%s_%s_DUP_only_reads.fasta' % (
            svloc, list(across_DUP_reads_name_list3)[0]))
        for_racon_fa = os.path.join(vr_dir, '%s_for_racon_Assembly.fasta' % (
            svloc))
        # print(for_racon_fa)
        with open(out_only_DUP_fa, 'r') as only_DUP_fa:
            for line in only_DUP_fa:
                if line[0] != '>':
                    line_dup = line.strip()
                    return line_dup


def get_line_ins(sample_name, chr, star, sv_type, end, reads, svloc,svlen):
    bam_merge = get_merge_bam_path(sample_name)
    workdir = os.path.join(workdir_0, sample_name)
    assembly_dir = os.path.join(workdir, 'assembly_INS', '%s' % sample_name)
    if not os.path.exists(assembly_dir):
        os.makedirs(assembly_dir)
    # print(line)
    sv_dir = os.path.join(assembly_dir, svloc)
    if not os.path.exists(sv_dir):
        os.mkdir(sv_dir)
    vr_dir = os.path.join(sv_dir, 'variant_reads')
    if not os.path.exists(vr_dir):
        os.mkdir(vr_dir)
    vr_txt = os.path.join(vr_dir, '%s_variant_ont_reads.txt' % svloc)
    with open(vr_txt, 'w') as out:
        for vr in reads:
            if 'tumor' in vr:
                out.write(vr + '\n')
    # region bam
    region = "%s:%s-%s" % (chr, int(star) - 10000, int(end) + 10000)
    ont_bam_region = os.path.join(vr_dir, '%s_flk10000_region.bam' % svloc)
    vr_bam = os.path.join(vr_dir, '%s_variant_ont_reads.bam' % svloc)
    vr_fq = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)  ### 组装
    vr_fq_2 = os.path.join(vr_dir,
                           '%s_variant_ont_reads_for_Iris.fastq' % svloc)  ### 组装
    variant_ont_qnames(vr_txt, bam_merge, region, ont_bam_region, vr_bam, vr_fq,
                       vr_fq_2)
    if os.path.getsize(vr_fq_2) > 0:
        reads_name_list = find_fq_MQ_max_reads(vr_fq_2)
    else:
        line_dup = '*'
        return line_dup
    greads_query_star_end_list, across_DUP_reads_name_list, across_DUP_reads_name_list2 = {},{},{}
    greads_query_star_end_list, across_DUP_reads_name_list, across_DUP_reads_name_list2 = greads_query_star_end(reads_name_list, vr_dir, svloc, star, end, vr_bam,svlen)
    if not greads_query_star_end_list:
        print('greads_query_star_end_list is empty')
    elif not across_DUP_reads_name_list:
        print('across_DUP_reads_name_list is empty')
    elif not across_DUP_reads_name_list2:
        print('across_DUP_reads_name_list2 is empty')
    across_DUP_reads_name_list3 = set(across_DUP_reads_name_list).intersection(across_DUP_reads_name_list2)
    if not across_DUP_reads_name_list3:
        # print(svloc)
        line_dup = '*'
        return line_dup
    else:
        out_only_DUP_fa = os.path.join(vr_dir, '%s_%s_DUP_only_reads.fasta' % (svloc, list(across_DUP_reads_name_list3)[0]))
        with open(out_only_DUP_fa, 'r') as only_DUP_fa:
            for line in only_DUP_fa:
                if line[0] != '>':
                    line_dup = line.strip()
                    return line_dup



def get_line_ins_all_reads(sample_name, chr, star, sv_type, end, reads, svloc,svlen):
    bam_merge = get_merge_bam_path(sample_name)
    workdir = os.path.join(workdir_0, sample_name)
    assembly_dir = os.path.join(workdir, 'assembly_INS', '%s' % sample_name)
    if not os.path.exists(assembly_dir):
        os.makedirs(assembly_dir)
    # print(line)
    sv_dir = os.path.join(assembly_dir, svloc)
    if not os.path.exists(sv_dir):
        os.mkdir(sv_dir)
    vr_dir = os.path.join(sv_dir, 'variant_reads')
    if not os.path.exists(vr_dir):
        os.mkdir(vr_dir)
    vr_txt = os.path.join(vr_dir, '%s_variant_ont_reads.txt' % svloc)
    with open(vr_txt, 'w') as out:
        for vr in reads:
            out.write(vr + '\n')
    # region bam
    region = "%s:%s-%s" % (chr, int(star) - 10000, int(end) + 10000)
    ont_bam_region = os.path.join(vr_dir, '%s_flk10000_region.bam' % svloc)
    vr_bam = os.path.join(vr_dir, '%s_variant_ont_reads.bam' % svloc)
    vr_fq = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)
    vr_fq_2 = os.path.join(vr_dir,
                           '%s_variant_ont_reads_for_Iris.fastq' % svloc)
    variant_ont_qnames(vr_txt, bam_merge, region, ont_bam_region, vr_bam, vr_fq,
                       vr_fq_2)
    if os.path.getsize(vr_fq_2) > 0:
        reads_name_list = find_fq_MQ_max_reads(vr_fq_2)
    else:
        line_dup = '*'
        return line_dup
    greads_query_star_end_list, across_DUP_reads_name_list, across_DUP_reads_name_list2 = {},{},{}
    greads_query_star_end_list, across_DUP_reads_name_list, across_DUP_reads_name_list2 = greads_query_star_end(reads_name_list, vr_dir, svloc, star, end, vr_bam,svlen)
    if not greads_query_star_end_list:
        print('greads_query_star_end_list is empty')
    elif not across_DUP_reads_name_list:
        print('across_DUP_reads_name_list is empty')
    elif not across_DUP_reads_name_list2:
        print('across_DUP_reads_name_list2 is empty')
    across_DUP_reads_name_list3 = set(across_DUP_reads_name_list).intersection(across_DUP_reads_name_list2)
    if not across_DUP_reads_name_list3:
        # print(svloc)
        line_dup = '*'
        return line_dup
    else:
        out_only_DUP_fa = os.path.join(vr_dir, '%s_%s_DUP_only_reads.fasta' % (svloc, list(across_DUP_reads_name_list3)[0]))
        with open(out_only_DUP_fa, 'r') as only_DUP_fa:
            for line in only_DUP_fa:
                if line[0] != '>':
                    line_dup = line.strip()
                    return line_dup



def change_DUP_in_Iris_vcf(Iris_in_vcf_x,Iris_in_vcf_DUP_x,sample_name):
    with open(Iris_in_vcf_DUP_x, 'w') as Iris_in_vcf_DUP_line:
        with open(Iris_in_vcf_x, 'r') as Iris_in_vcf_line:
            for line in Iris_in_vcf_line:
                if line[0] != '#':
                    annotation = line.strip().split('\t')
                    chr = annotation[0]
                    star = annotation[1]
                    sv_id = annotation[2]
                    print(sv_id)
                    information = annotation[7]
                    infor = information.split(';')
                    sv_type = infor[1].split('=')[1]
                    svlen = infor[2].split('=')[1]
                    end = infor[3].split('=')[1]
                    reads = (infor[5].split('=')[1]).split(',')
                    svloc = '_'.join([sample_name, chr, star, end, sv_type, sv_id])
                    if sv_type == 'DUP' :
                        if int(svlen) <= 10000: ## 参考Jasmine的参数 只把1万以下的dup的序列提出来
                            line_dup = get_line_dup(sample_name, chr, star, sv_type, end, reads, svloc, svlen)
                            if line_dup == 0:
                                print(line)
                            else:
                                id = sv_id.split('.')
                                infor = annotation[7].split(';')
                                infor_split = infor[1].split('=')
                                infor_split[1] = 'INS'
                                infor[1] = '='.join([str(i) for i in infor_split])
                                annotation[7] = ';'.join([str(i) for i in infor])
                                # ((annotation[7].split(';'))[1].split('='))[1] == 'INS'
                                annotation[2] = '.'.join([str(i) for i in id])
                                annotation[4] = line_dup
                                line2 = '\t'.join([str(i) for i in annotation]) + '\n'
                                # Iris_in_vcf_DUP_line.write(line2)
                                if line_dup:
                                    if len(str(line_dup)) > 3:
                                        Iris_in_vcf_DUP_line.write(line2)
                    elif annotation[4] == '<INS>':
                        line_dup = get_line_ins(sample_name, chr, star, sv_type, end, reads, svloc,svlen)
                        annotation[4] = line_dup
                        line = '\t'.join([str(i) for i in annotation]) + '\n'
                        # Iris_in_vcf_DUP_line.write(line)
                        if line_dup:
                            if len(str(line_dup)) > 3:
                                Iris_in_vcf_DUP_line.write(line)
                    else:
                        Iris_in_vcf_DUP_line.write(line)
                else:
                    Iris_in_vcf_DUP_line.write(line)


def change_DUP_in_Iris_vcf_all_reads(Iris_in_vcf_x,Iris_in_vcf_DUP_x,sample_name):
    with open(Iris_in_vcf_DUP_x, 'w') as Iris_in_vcf_DUP_line:
        with open(Iris_in_vcf_x, 'r') as Iris_in_vcf_line:
            for line in Iris_in_vcf_line:
                if line[0] != '#':
                    annotation = line.strip().split('\t')
                    chr = annotation[0]
                    star = annotation[1]
                    sv_id = annotation[2]
                    print(sv_id)
                    information = annotation[7]
                    infor = information.split(';')
                    sv_type = infor[1].split('=')[1]
                    svlen = infor[2].split('=')[1]
                    end = infor[3].split('=')[1]
                    reads = (infor[5].split('=')[1]).split(',')
                    svloc = '_'.join([sample_name, chr, star, end, sv_type, sv_id])
                    if sv_type == 'DUP' :
                        if int(svlen) <= 10000: ## in Jasmine only get dup seq less 1W
                            line_dup = get_line_dup_all_reads(sample_name, chr, star, sv_type, end, reads, svloc, svlen)
                            if line_dup == 0:
                                print(line)
                            else:
                                id = sv_id.split('.')
                                infor = annotation[7].split(';')
                                infor_split = infor[1].split('=')
                                infor_split[1] = 'INS'
                                infor[1] = '='.join([str(i) for i in infor_split])
                                annotation[7] = ';'.join([str(i) for i in infor])
                                # ((annotation[7].split(';'))[1].split('='))[1] == 'INS'
                                annotation[2] = '.'.join([str(i) for i in id])
                                annotation[4] = line_dup
                                line2 = '\t'.join([str(i) for i in annotation]) + '\n'
                                # Iris_in_vcf_DUP_line.write(line2)
                                if line_dup:
                                    if len(str(line_dup)) > 3:
                                        Iris_in_vcf_DUP_line.write(line2)
                    elif annotation[4] == '<INS>': ##
                        line_dup = get_line_ins_all_reads(sample_name, chr, star, sv_type, end, reads, svloc,svlen)
                        annotation[4] = line_dup
                        line = '\t'.join([str(i) for i in annotation]) + '\n'
                        if line_dup:
                            if len(str(line_dup)) > 3:
                                Iris_in_vcf_DUP_line.write(line)
                    else:
                        Iris_in_vcf_DUP_line.write(line)
                else:
                    Iris_in_vcf_DUP_line.write(line)


def IRSI_function(sample_name):
    print(sample_name)
    vcf_file = get_vcf_file_path(sample_name)
    CSV_file = get_CSV_file_path(sample_name)
    reads_in = get_merge_bam_path(sample_name)
    vcf_out = os.path.join('{sample2}_iris_out.vcf'.format(sample2=sample_name))
    log_out = os.path.join(workdir_0,
                           '{sample}/{sample2}_iris_out.log'.format(sample=sample_name,sample2=sample_name))
    out_dir = os.path.join(workdir_0,'{sample}/'.format(sample=sample_name))
    Iris_result_file = os.path.join(workdir_0 , '{sample_name}/resultsstore.txt'.format(
        sample_name=sample_name))
    workdir = os.path.join(workdir_0,"{sample}".format(sample=sample_name))
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    all_sv_file = pd.read_csv(CSV_file, sep=',')
    all_sv_id_list = []
    for svid in list(set(all_sv_file['svid'])):
        svtype = svid.split('.')
        if svtype[1] in ['INS', 'DUP']:
            all_sv_id_list.append(svid)
    Iris_in_vcf = os.path.join(workdir, '{CCA1}_INS_difference_sv_id_for_Iris.vcf'.format(CCA1=sample_name))
    grep_sv_id = '#'
    for sv_id_1 in all_sv_id_list:
        grep_sv_id = grep_sv_id + '|' + sv_id_1
    os.system(
        "grep -E -w '{grep_sv_id}' {vcf_file} > {Iris_in_vcf} ".format(grep_sv_id=grep_sv_id, vcf_file=vcf_file,Iris_in_vcf=Iris_in_vcf))
    Iris_in_vcf_DUP = os.path.join(workdir, '{CCA1}_INS_and_DUP_for_Iris.vcf'.format(CCA1=sample_name))
    change_DUP_in_Iris_vcf(Iris_in_vcf, Iris_in_vcf_DUP, sample_name)
    txt = os.popen(
        "cd {workdir} && iris genome_in=/nanopore/ref/hg38_mainChr.fa vcf_in={Iris_in_vcf} reads_in={reads_in} vcf_out={vcf_out} out_dir = {out_dir} genome_buffer=1000 --keep_long_variants threads = 20".format(
            workdir = workdir,sample=sample_name, Iris_in_vcf=Iris_in_vcf_DUP, reads_in=reads_in, out_dir=out_dir,
            vcf_out=vcf_out, log_out=log_out)).read()
    with open(log_out, 'w') as log_out_line:
        log_out_line.write(txt)


def get_merge_bam_path(sample_name):
    bam_merge = os.path.join(somatic_dir, 'Cancer_cells', sample_name,
                             '%s_Cancer_cells_minimap2_sorted_tag.bam' % sample_name)
    return bam_merge


def get_CSV_file_path(sample_name):
    CSV_dir =  '{work_dir}/2_std_test'.format(work_dir=work_dir)
    CSV_path = os.path.join(CSV_dir,sample_name,"{sample}_somatic_sv_state-20221230-final.csv".format(sample=sample_name))
    return CSV_path



def get_vcf_file_path(sample_name):
    vcf_dir =  '{ms_dir}/{sample_name}/sniffles2'.format(sample_name=sample_name,ms_dir=ms_dir)
    vcf_file = os.path.join(vcf_dir,"{sample}_merge_minimap2_sniffles_v2.vcf".format(sample=sample_name))
    return vcf_file





Pool_core = 10

all_file_list = ['Sample1']

pools =Pool(Pool_core)

for  sample_name  in all_file_list:
    pools.apply_async(IRSI_function,(sample_name,))

pools.close()
pools.join()
del pools






############################################## Annocation Sequences ##############################################


def fa_hg_38(region):
    fa = pysam.FastaFile(hg38_fa)
    seq = fa.fetch(region=region)
    return seq


def Repeat_masker(sample_name):
    Iris_out_vcf_file = os.path.join(workdir_0,'{CCA_test2}/{CCA_test2}_iris_out.vcf'.format(CCA_test2=sample_name))
    with open (Iris_out_vcf_file,'r') as Iris_out_vcf:
        for line in Iris_out_vcf:
            if line[0] != '#':
                annotation = line.strip().split('\t')
                chr = annotation[0]
                star = annotation[1]
                sv_id = annotation[2]
                svtype  = sv_id.split('.')[1]
                sequence = annotation[4]
                sv_len = len(sequence[1::])
                end = int(star) + sv_len
                region_star = chr + ':' + str(int(star)-50) + '-' + str(star)
                region_end = chr + ':' + str(end) + '-' + str(int(end)+50)
                seq_star = fa_hg_38(region_star)
                seq_end = fa_hg_38(region_end)
                seq = seq_star + sequence[1::] + seq_end
                sv_id2 = '_'.join([sample_name, chr, str(star), str(int(star) + 1), 'INS',
                                   sv_id])
                sv_id3 = '_'.join([sample_name, chr, 'INS', sv_id, 'srm'])
                fa_for_rm = os.path.join(workdir_0,
                                         'RepeatMasker/{CCA1}/'.format(CCA1=sample_name), sv_id2,
                                         '{sv_id}_srm2hg38_INS_sequence.fasta'.format(sv_id=sv_id2))
                fa_too_long = os.path.join(workdir_0,'RepeatMasker/{CCA1}/'.format(CCA1=sample_name), sv_id2,
                                           '{sv_id}_srm2hg38_INS_long_sequence.fasta'.format(sv_id=sv_id2))
                out_dir = os.path.join(workdir_0,
                                       'RepeatMasker/{CCA1}/'.format(CCA1=sample_name), sv_id2,
                                       'RepeatMasker')
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                newline1 = '>%s' % sv_id3 + '\n'
                newline2 = seq + '\n'
                newline3 = newline1 + newline2
                if len(sequence) < 20000:
                    with open(fa_for_rm, 'w') as fa2:
                        fa2.write(newline3)
                else:
                    with open(fa_too_long, 'w') as fa3:
                        fa3.write(newline3)
                    newline2 = sequence[0:10000] + '\n'
                    newline3 = newline1 + newline2
                    with open(fa_for_rm, 'w') as fa2:
                        fa2.write(newline3)
                os.system('/biosoft/RepeatMasker/RepeatMasker -species human -engine RMBlast -q -parallel 4 %s -dir %s' % (fa_for_rm, out_dir))


pools =Pool(Pool_core)
for  sample_name  in all_file_list:
    pools.apply_async(Repeat_masker,(sample_name,))

pools.close()
pools.join()
del pools




for sample_name in all_file_list:
    reannotation_path = '{workdir_0}/ALL_reannotation/{PAAD}'.format(workdir_0=workdir_0,PAAD = sample_name)
    if not os.path.exists(reannotation_path):
        os.mkdir(reannotation_path)
    reannotation_path_bulid_txt = os.path.join(reannotation_path,'ALL_reannotation_path_bulid.txt')
    reannotation_path_bulid = open(reannotation_path_bulid_txt,'w')
    Iris_out_vcf_file = os.path.join(workdir_0,'{CCA_test2}/{CCA_test2}_iris_out.vcf'.format(CCA_test2=sample_name))
    with open (Iris_out_vcf_file,'r') as Iris_out_vcf:
        for line in Iris_out_vcf:
            if line[0] != '#':
                annotation = line.strip().split('\t')
                chr = annotation[0]
                star = annotation[1]
                sv_id = annotation[2]
                sequence = annotation[4]
                if sequence != '<DUP>' :
                    sv_id2 = '_'.join([sample_name, chr, str(star), str(int(star) + 1), 'INS', sv_id])
                    sv_id3 = '_'.join([sample_name, chr, 'INS', sv_id,'srm'])
                    out_dir = os.path.join(workdir_0,'RepeatMasker/{CCA1}/{sv_id2}'.format(CCA1=sample_name,sv_id2=sv_id2))
                    hg38_location = '_'.join([chr, str(star), str(int(star) + 1)])
                    line = '\t'.join([sv_id2, out_dir, hg38_location]) + '\n'
                    # print(line)
                    reannotation_path_bulid.write(line)
    reannotation_path_bulid.close()
    os.system('cp -r /nanopore/script_dir/script {reannotation_path}'.format(reannotation_path=reannotation_path)) #script is in github snakemake script,need downlode
    with open (os.path.join(reannotation_path,'script','Snakefile'),'w') as reannotation_path_bulid:
        with open ('/nanopore/script_dir/script/Snakefile','r') as reannotation_path_bulid2:
            n = 0
            for line in reannotation_path_bulid2:
                n += 1
                if n== 4:
                    line = '''workdir: "{reannotation_path}" '''.format(reannotation_path = reannotation_path) + '\n'
                reannotation_path_bulid.write(line)
    os.system('cd {reannotation_path}/script && snakemake --cores 10 '.format(reannotation_path=reannotation_path))






