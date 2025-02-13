
from multiprocessing import Pool
import os


samtools = "/biosoft/samtools-1.9/samtools"
bedtools = "/biosoft/bedtools-2.30.0"
sniffles2 = "/biosoft/sniffles2/bin/sniffles"
minimap2 = "/biosoft/minimap2-2.17/minimap2"


hg38_fa = "/nanopore/ref/hg38_mainChr.fa"
tr_bed = "/nanopore/ref/human_GRCh38_no_alt_analysis_set.trf.bed"

work_dir = "/nanopore/NewSample"
program_dir = os.path.join(work_dir, '01.program') # sh script is in sh_script in github
ms_dir = os.path.join(work_dir, '02.minimap2_sniffles')
coverage_dir = os.path.join(work_dir,'03.coverage')
somatic_dir = os.path.join(work_dir, "04.Somatic_minimap2_2.4")

dir_list = [program_dir, ms_dir, somatic_dir,coverage_dir]
for dir in dir_list:
    if not os.path.exists(dir):
        os.mkdir(dir)

Pool_core = 10

######################### Step1 map and add tag #################

def bam_tag(script_tag,fq,new_ID,somatic_sample_dir,hg38_fa,stdout,stderr):
    os.system('/bin/bash {script} {a} {b} {c} {d} 1 > {out} 2 > {err}'.format(
            script=script_tag, a=fq, b=new_ID, c=somatic_sample_dir, d=hg38_fa,out=stdout, err=stderr))


all_file_list = ['Sample1']
TB_list = ['blood','tumor']


pools = Pool(Pool_core)
for sampleid in all_file_list:
    for sampletype in TB_list:
        if sampletype == 'tumor':
            fq = os.path.join('/nanopore/NewSample/01.fastq_dir',sampletype,'%s_%s.cat.fastq.gz' % (sampleid,sampletype) )
            tag = 'tumor'
        if sampletype == 'blood':
            fq = os.path.join('/nanopore/NewSample/01.fastq_dir',sampletype,'%s_%s.cat.fastq.gz' % (sampleid,sampletype) )
            tag = 'blood'
        script_tag = os.path.join(program_dir, "ONT_bam_add_tag_minimap2_%s.sh" % tag)
        somatic_sample_dir = os.path.join(somatic_dir, tag, sampleid)
        if not os.path.exists(somatic_sample_dir):
            os.makedirs(somatic_sample_dir)
        new_ID = '%s_%s' % (sampleid, tag)
        stdout = os.path.join(somatic_sample_dir, '%s_add_tag.o' % new_ID)
        stderr = os.path.join(somatic_sample_dir, '%s_add_tag.e' % new_ID)
        pools.apply_async(bam_tag, args=(script_tag,fq,new_ID,somatic_sample_dir,hg38_fa,stdout,stderr))


pools.close()
pools.join()
del pools



######################### Step2 Merge bam and call SV #################

for sampleid in all_file_list:
    sniffles2_sample_dir = os.path.join(ms_dir, sampleid, 'sniffles2')
    if os.path.exists(sniffles2_sample_dir):
        os.system("rm -rf %s" % sniffles2_sample_dir)
    os.makedirs(sniffles2_sample_dir)
    sv_vcf2 = os.path.join(sniffles2_sample_dir, '%s_merge_minimap2_sniffles_v2.vcf' % sampleid)
    bam_sort_tag_tumor = os.path.join(somatic_dir, 'Cancer_cells', sampleid,
                                      '%s_Cancer_cells_minimap2_sorted_tag.bam' % sampleid)
    sampelid_blood = sampleid
    bam_sort_tag_blood = os.path.join(somatic_dir, 'blood', sampelid_blood, '%s_blood_minimap2_sorted_tag.bam' % sampelid_blood)
    bam_merge_path = os.path.join(somatic_dir, 'Merge_bam', sampleid)
    bam_merge = os.path.join(somatic_dir, 'Merge_bam', sampleid,
                             '%s_minimap2_sorted_blood_tumor_merge.bam' % sampleid)
    if os.path.exists(bam_merge):
        os.system("rm -rf %s" % bam_merge)
    # os.mkdir(sniffles2_sample_dir)
    os.mkdir(bam_merge_path)
    # script
    script_ms = os.path.join(ms_dir, sampleid, '%s_bam_merge_sniffles.2nd.sh' % sampleid)
    with open(script_ms, 'w') as out:
        out.write("#! /bin/bash" + '\n')
        out.write('''echo "$(date) 1. Start to merge bam: %s" ''' % sampleid + '\n')
        out.write("{samtools} merge -@ 25 -h {bam_b} {out_bam} {bam_b} {bam_t}".format(
            samtools=samtools, bam_b=bam_sort_tag_blood, bam_t=bam_sort_tag_tumor, out_bam=bam_merge) + '\n')
        out.write('''/biosoft/samtools-1.9/samtools index -@ 15 %s \n''' % bam_merge)
        out.write('''echo "$(date) 1. Finish to merge bam: %s" ''' % sampleid + '\n')
        out.write('''echo "$(date) 2. Start to sniffles2: %s" ''' % sampleid + '\n')
        out.write('''export PATH="/data/fs01/wangzf/software/anaconda3/bin:$PATH" \n''')
        out.write("source activate nanoplot \n")
        out.write(
            "{sniffles} -i {bam_sort} -v {vcf} --tandem-repeats {tr} -t 4 --minsupport 1 --mapq 10 --min-alignment-length 1000 "
            "--output-rnames --allow-overwrite --long-ins-length 100000 --reference {ref}".format(
                sniffles=sniffles2, bam_sort=bam_merge, vcf=sv_vcf2, tr=tr_bed,
                ref=hg38_fa) + '\n')  # --minsvlen default 35bp
        out.write('''echo "$(date) 2. Finish to sniffles2: %s" ''' % sampleid + '\n')
        out.write('''echo "$(date) 3. Start to Dele merbam: %s" ''' % sampleid + '\n')
        out.write(" rm {out_bam} " .format(out_bam=bam_merge)+'\n')
        out.write('''echo "$(date) 3. Finish to Dele merbam: %s" ''' % sampleid + '\n')



def bam_tag_2(script_tag,stdout,stderr):
    os.system('/bin/bash {script}  1 > {out} 2 > {err}'.format(
            script=script_tag, out=stdout, err=stderr))


for sampleid in all_file_list:
    script_ms = os.path.join(ms_dir,sampleid, '%s_bam_merge_sniffles.2nd.sh' % sampleid)
    stdout = script_ms.replace(".sh", ".o")
    stderr = script_ms.replace(".sh", ".e")
    bam_tag_2(script_ms, stdout, stderr)


######################### Step3 Filter SV #################



########################################################################################
# Somatic post-process:
########################################################################################
def sniffles2_vcftobed(vcf, bed, support):
    """vcf to bed
    header: chrom start end (chr start end) svtype id length RE RNAMES IMPRECISE/PRECISE STD STRAND
    1-based to 0-based
    """
    bedout = open(bed, 'w')
    with open(vcf, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            else:
                content = line.split()
                chr1 = content[0]
                start = content[1]
                svid = content[2]
                info = content[7].split(';')
                fil = info[0]
                svtype = info[1].split('=')[1]
                if chr1 not in ['chrM', 'chrY']:
                    if svtype == 'BND':
                        # if ']' in content[4]:
                        #     content4 = content[4].split(']')
                        # else:
                        #     content4 = content[4].split('[')
                        # end = content4[1].split(':')[1]
                        end = content[4].replace('N', '').split(':')[1][:-1]
                        chr2 = info[-2].split('=')[1]
                        if chr2 not in ['chrM', 'chrY']:
                            svlen = 1
                            re = info[2].split('=')[1]
                            if int(re) >= support:
                                rnames = info[3].split('=')[1]
                                # std = info[-1].split('=')[1]
                                # suptype = info[-4].split('=')[1]
                                anotation = ['TRA', svid, str(abs(int(svlen))), re, rnames, fil]
                                newline1 = '\t'.join(
                                    [chr1, str(int(start) - 1), start, chr2, str(int(end) - 1), end] + anotation) + '\n'
                                newline2 = '\t'.join(
                                    [chr2, str(int(end) - 1), end, chr1, str(int(start) - 1), start] + anotation) + '\n'
                                newline = newline1 + newline2
                                bedout.write(newline)
                    else:
                        svlen = info[2].split('=')[1]
                        re = info[4].split('=')[1]
                        if int(re) >= support:
                            if abs(int(svlen)) >= 50:
                                rnames = info[5].split('=')[1]
                                # std = info[-1].split('=')[1]
                                # suptype = info[-4].split('=')[1]
                                end = info[3].split('=')[1]
                                anotation = ["NA","NA","NA",svtype, svid, str(abs(int(svlen))), re, rnames, fil]
                                newline = '\t'.join([chr1, str(int(start) - 1), end] + anotation) + '\n'
                                bedout.write(newline)
                    bedout.flush()
    bedout.close()


for sampleid in all_file_list:
    sniffles2_sample_dir = os.path.join(ms_dir, sampleid, 'sniffles2')
    sv_vcf = os.path.join(sniffles2_sample_dir, '%s_merge_minimap2_sniffles_v2.vcf' % sampleid)
    sv_bed = sv_vcf.replace(".vcf", ".bed")
    sniffles2_vcftobed(sv_vcf, sv_bed, 3)
    sv_bed_tumor_somatic = sv_bed.replace(".bed", "_tumor_somatic.bed")
    sv_bed_blood_somatic = sv_bed.replace(".bed", "_blood_somatic.bed")
    with open(sv_bed_tumor_somatic, 'w') as out:
        with open(sv_bed, 'r') as f:
            for line in f:
                c = line.strip().split('\t')
                blood_c = 0
                tumor_c = 0
                # if c[6] != "TRA":
                #     variant_reads = c[7]
                # else:
                variant_reads = c[10]
                variant_reads_list = variant_reads.split(',')
                all_c = len(variant_reads_list)
                for reads in variant_reads_list:
                    if reads[:5] == "tumor":
                        tumor_c += 1
                    elif reads[:5] == "blood":
                        blood_c += 1
                if blood_c == 0:
                    out_list = c + list(map(str, [all_c, tumor_c, blood_c]))
                    out.write('\t'.join(out_list) + '\n')
                    out.flush()
    with open(sv_bed_blood_somatic, 'w') as out:
        with open(sv_bed, 'r') as f:
            for line in f:
                c = line.strip().split('\t')
                blood_c = 0
                tumor_c = 0
                # if c[6] != "TRA":
                #     variant_reads = c[7]
                # else:
                variant_reads = c[10]
                variant_reads_list = variant_reads.split(',')
                all_c = len(variant_reads_list)
                for reads in variant_reads_list:
                    if reads[:5] == "tumor":
                        tumor_c += 1
                    elif reads[:5] == "blood":
                        blood_c += 1
                if tumor_c == 0:
                    out_list = c + list(map(str, [all_c, tumor_c, blood_c]))
                    out.write('\t'.join(out_list) + '\n')
                    out.flush()


def sniffles_svtype(bed):
    """select RE>=3
    delete chrY,chrM
    """
    svtype_dic = {}
    with open(bed, 'r') as f:
        for line in f:
            content = line.strip().split('\t')
            chrom = content[0]
            if chrom != 'chrM' and chrom != 'chrY':
                svtype_x = content[6]
                if svtype_x not in svtype_dic:
                    svtype_file = bed.split('.bed')[0] + '_%s.bed' % svtype_x
                    svtype_dic[svtype_x] = open(svtype_file, 'w')
                newline = '\t'.join(content) + '\n'
                svtype_dic[svtype_x].write(newline)
    for svtype_x in svtype_dic:
        svtype_dic[svtype_x].close()


for sampleid in all_file_list:
    sniffles2_sample_dir = os.path.join(ms_dir, sampleid, 'sniffles2')
    sv_bed_somatic = os.path.join(sniffles2_sample_dir, '%s_merge_minimap2_sniffles_v2_tumor_somatic.bed' % sampleid)
    sniffles_svtype(sv_bed_somatic)
















