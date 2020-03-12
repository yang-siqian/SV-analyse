#!/usr/bin/env python
import os
import sys
import random
import time
from datetime import datetime
import threading

usage = '''Program:  smrtsv (SV detection from pacbio sequencing)
Version:  0.0.1 
Usage:    python3 smrtmax.py <snakefile.config>
Contact:  YSQ '''


root_dir = ""
out_dir = ""
reference = ""
process_name = ""
runtimekey = ""
n_thread = 20

#align
fasta_file = ""
#detect/assemble
alignments = ""
#assemble
# candidates = ""
#call
local_assembly_alignments = ""

TimeFormat = '%m/%d/%Y %H:%M:%S'


def myprint(string):
    print('[' + datetime.now().strftime(TimeFormat) + ']' + string)
    return


def parse_config_file(config_file):
    global root_dir
    global out_dir
    global reference
    global process_name
    global n_thread
    global fasta_file
    global alignments
    global candidates
    global local_assembly_alignments

    root_dir = os.path.split(os.path.abspath(__file__))[0]

    with open(config_file, 'r') as config_file_fp:
        while 1:
            line = config_file_fp.readline()
            if not line:
                break
            line = line.strip()
            if '#' in line:
                line = line.split('#')[0]
            if len(line) == 0:
                continue
            if '=' not in line:
                continue
            line = line.split('=')
            if len(line) < 2:
                continue
            key = line[0].strip().lower()
            value = '='.join(line[1:])
            value = value.strip()
            if (not key) or (not value):
                continue

            if key == 'out_dir':
                out_dir = os.path.abspath(value)
            elif key == 'reference':
                reference = os.path.abspath(value)
            elif key == 'process_name':
                process_name = value
            elif key == 'n_thread' and int(value) > 0:
                n_thread = int(value)

            #align
            elif key == 'fasta_file':
                fasta_file = os.path.abspath(value)

            # detect
            elif key == 'alignments':
                alignments = os.path.abspath(value)

            # assemble
            elif key == 'candidates':
                candidates = os.path.abspath(value)

            # call
            elif key == 'local_assembly_alignments':
                local_assembly_alignments = os.path.abspath(value)



    if reference == "":
        myprint('reference is required')
        sys.exit()

    elif process_name == "":
        myprint('process name is required!')
        sys.exit()

    elif out_dir == "":
        myprint('out_dir is required!')
        sys.exit()

    elif process_name == "align" and fasta_file == "":
        myprint('run_align: fasta_file is required!')
        sys.exit()

    elif process_name == "detect" and alignments == "":
        myprint('run_detect: alignments is required!')
        sys.exit()

    elif process_name == "assemble" and candidates == "" and alignments == "":
        myprint('run_assemble: candidates and alignments is required!')
        sys.exit()

    elif process_name == "call" and local_assembly_alignments == "":
        myprint('run_call: local_assembly_alignments is required!')
        sys.exit()



def main():
    global runtimekey
    if len(sys.argv) < 2:
        myprint(usage)
        sys.exit()

    myprint('program started')
    myprint('reading snakefile.config file')
    config_file = sys.argv[1]
    parse_config_file(config_file)

    random.seed()
    runtimekey = random.randint(10000000000, 99999999999)

    myprint('running tasks: %s' % process_name)

    if process_name == "index":
        i = Index(out_dir, reference, process_name)
        i.run_index_cmd()
    elif process_name == "align":
        a = Align(out_dir, reference, process_name, n_thread, fasta_file)
        a.run_align_cmd()
    elif process_name == "detect":
        d = Detect(root_dir, out_dir, reference, alignments)
        d.run_detect_cmd()
    elif process_name == "assemble":
        # 1.创建regions_by_contig、mhap_assembly、local_assemblies、临时文件夹tmp
        regions_by_contig = os.path.join(out_dir, 'regions_by_contig')
        check_dir_exist_rm(regions_by_contig)

        mhap_assembly = os.path.join(out_dir, 'mhap_assembly')
        check_dir_exist_rm(mhap_assembly)

        local_assemblies = os.path.join(out_dir, 'local_assemblies')
        check_dir_exist_rm(local_assemblies)

        task_list = []

        # 2.按行读取candidates.bed, contig=chrI.bed, chr_bed为contig的绝对路径，将相同的contig写入同一个chr_bed文件中
        with open(candidates, 'r') as fc:
            for line in fc.readlines():
                # chr = chrI, start起始位置，end结束位置
                chr = line.strip().split()[0]
                start = line.strip().split()[1]
                end = line.strip().split()[2]

                stand_region = chr + ':' + start + '-' + end

                contig = line.strip().split()[0] + '.bed'
                chr_bed = os.path.join(regions_by_contig, contig)
                with open(chr_bed, 'a') as fb:
                    fb.write(line)

                task_list.append(stand_region)
        n = len(task_list)
        print(n)
        for a in range(0, n, 5):
            t1 = Assemble(a, task_list[a], root_dir, out_dir, reference, alignments)
            print(task_list[a])
            t1.start()
            if a + 1 < n:
                t2 = Assemble(a + 1, task_list[a+1], root_dir, out_dir, reference, alignments)
                t2.start()
            if a + 2 < n:
                t3 = Assemble(a + 2, task_list[a+2], root_dir, out_dir, reference, alignments)
                t3.start()
            if a + 3 < n:
                t4 = Assemble(a + 3, task_list[a+3], root_dir, out_dir, reference, alignments)
                t4.start()
            if a + 4 < n:
                t5 = Assemble(a + 4, task_list[a + 4], root_dir, out_dir, reference, alignments)
                t5.start()


            t1.join()
            t2.join()
            t3.join()
            t4.join()
            t5.join()


        list_filename_list = os.listdir(local_assemblies)

        local_assemblies_list = []
        for base_filename in list_filename_list:
            bam_name = base_filename.replace('list.txt', 'bam')
            out_bam = os.path.join(local_assemblies, bam_name)

            # 3.生成local_assemblies目录下的local_assembly_alignments.chrI.bam
            list_filename = os.path.join(local_assemblies, base_filename)
            input_fai = reference + '.fai'

            # a14 collect_assembly_alignments
            cmd = "while read file; do sed 's/\/0_[0-9]\+//' $file; done < %s | samtools view -Sbu -t %s - | bamleftalign -f %s | samtools sort -O bam -T %s -o %s" % (list_filename, input_fai, reference, list_filename, out_bam)
            os.system(cmd)
            step = 14
            finish_cmd = 'echo %d_finish' % step
            only_check_exist(out_bam, finish_cmd)

            # 将生成的local_assembly_alignments.chrI.bam添加到local_assemblies_list列表中
            local_assemblies_list.append(out_bam)

        # 4.定义local_assembly_alignments.bam
        assembly_alignments = os.path.join(out_dir, 'local_assembly_alignments.bam')

        local_assemblies_str = ' '.join(local_assemblies_list)

        # 5.如果local_assemblies_list列表长度大于1，执行merge_cmd,否则执行view_cmd
        if len(local_assemblies_list) > 1:
            merge_cmd = "samtools merge %s " % assembly_alignments + local_assemblies_str
            os.system(merge_cmd)
        else:
            view_cmd = "samtools view -b -o %s " % assembly_alignments + local_assemblies_str
            os.system(view_cmd)

        # 6.给local_assembly_alignments.bam建索引
        index_cmd = "samtools index %s" % assembly_alignments
        os.system(index_cmd)




    elif process_name == "call":
        c = Call(root_dir, out_dir, reference, alignments, n_thread, local_assembly_alignments)
        c.run_call_cmd()



    finish_file = '%s/%s.finish' % (out_dir, process_name)
    if os.path.exists(finish_file):
        myprint('program finished')



def only_check_exist(file, cmd):
    if os.path.exists(file):
        os.system(cmd)
    else:
        myprint("Error! %s is not exist!" % file)
        sys.exit()


def check_file_exist(file, cmd):
    if os.path.exists(file) and os.stat(file).st_size > 0:
        os.system(cmd)
    else:
        myprint("Error! %s is not exist or size is zero!" % file)
        sys.exit()


def check_dir_exist(dir):
    if os.path.exists(dir):
        myprint("%s has been exist!" % dir)
    else:
        os.mkdir(dir)


def check_dir_exist_rm(dir):
    if os.path.exists(dir):
        os.system("rm -rf %s/*" % dir)
    else:
        os.mkdir(dir)



class Index(object):
    def __init__(self, out_dir, reference, process_name):
        self.reference = reference
        self.out_dir = out_dir
        self.process_name = process_name


    def run_index_cmd(self):

        ref_path = os.path.dirname(self.reference)
        ref_name = os.path.basename(self.reference)
        cd_cmd = "cd %s \n" % ref_path

        fai_cmd = cd_cmd + "samtools faidx %s;" % self.reference
        check_file_exist(self.reference, fai_cmd)

        out_fai = "%s.fai" % self.reference
        out_sa = ref_name + '.sa'
        sa_cmd = cd_cmd + "sawriter %s %s ;" % (out_sa, self.reference)
        check_file_exist(out_fai, sa_cmd)


        out_ctab = ref_name + '.ctab'
        sa_file = "%s.sa" % self.reference
        ctab_cmd = cd_cmd + "printTupleCountTable %s > %s ;" % (self.reference, out_ctab)
        check_file_exist(sa_file, ctab_cmd)

        ctab_file = "%s.ctab" % self.reference
        finish_cmd = 'echo "program finished" > %s/%s.finish' % (self.out_dir, self.process_name)
        check_file_exist(ctab_file, finish_cmd)



class Align(object):
    def __init__(self, out_dir, reference, process_name, n_thread, fasta_file):
        self.out_dir = out_dir
        self.reference = reference
        self.process_name = process_name
        self.n_thread = n_thread
        self.fasta_file = fasta_file


    def run_align_cmd(self):

        # 1.定义ref_sa、ref_ctab、alignments文件夹和out_bam文件绝对路径
        ref_sa = self.reference + '.sa'
        ref_ctab = self.reference + '.ctab'

        alignments_path = os.path.join(self.out_dir, 'alignments')
        check_dir_exist_rm(alignments_path)

        out_bam = os.path.join(alignments_path, '0.bam')

        # 2.判断ref_sa、ref_ctab、fasta_file、reference是否存在
        if os.path.exists(ref_sa) and os.path.exists(ref_ctab) and os.path.exists(self.fasta_file) and os.path.exists(self.reference) and os.stat(ref_sa).st_size > 0 and os.stat(ref_ctab).st_size > 0 and os.stat(self.reference).st_size > 0:
            cmd1 = 'cd %s ;blasr %s %s -unaligned /dev/null -out /dev/stdout -sam -sa %s -ctab %s -nproc %d -clipping subread -bestn 2 -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 50 | samtools view -F 0x4 -hS - | samtools sort -@ 1 -O bam -T 0 -o %s -;' % (alignments_path, self.fasta_file, self.reference, ref_sa, ref_ctab, self.n_thread, out_bam) + '\n'
            os.system(cmd1)
        else:
            myprint("required file is not all exist!")
            sys.exit()

        # 3.out_bam是否生成
        cmd2 = 'cd %s ; samtools index %s;' % (alignments_path, out_bam) + '\n'
        check_file_exist(out_bam, cmd2)


        # 4.out_fai是否生成
        out_fai = "%s.fai" % out_bam
        finish_cmd = 'echo "program finished" > %s/%s.finish' % (self.out_dir, self.process_name)
        check_file_exist(out_fai, finish_cmd)



class Detect(object):
    def __init__(self, root_dir, out_dir, reference, alignments):
        self.root_dir = root_dir
        self.out_dir = out_dir
        self.reference = reference
        self.alignments = alignments


        self.min_clipping = 500
        self.bin_size = 500
        self.merge_distance = 500
        self.slop = 10000
        self.mapping_quality = 0
        self.min_length = 50
        self.min_support = 5
        self.max_support = 100
        self.min_coverage = 5
        self.max_coverage = 100
        self.assembly_window_size = 60000
        self.assembly_window_slide = 20000
        self.min_hardstop_support = 11
        self.max_candidate_length = 60000


    def run_detect_cmd(self):
        self.find_gaps_in_aligned_reads()
        self.find_hardstops_in_aligned_reads()
        self.calculate_coverage_per_batch()
        self.identify_gaps_in_reference_assembly()
        self.classify_gaps_in_aligned_reads("deletion")
        self.classify_gaps_in_aligned_reads("insertion")
        self.merge_coverage_per_batch()
        self.merge_gap_support_from_aligned_reads("deletion")
        self.merge_gap_support_from_aligned_reads("insertion")
        self.annotate_coverage_of_merged_gap_support("deletion")
        self.annotate_coverage_of_merged_gap_support("insertion")
        self.collect_hardstops()
        self.create_hardstop_breakpoints()
        self.filter_candidates("deletion")
        self.filter_candidates("insertion")
        self.combine_filtered_candidates()
        self.create_genomic_bins()
        self.count_hardstops_per_genomic_bin()
        self.filter_and_merge_adjacent_hardstop_bins()
        self.merge_filtered_candidates()
        self.create_windows_for_tiled_assemblies()
        self.merge_filtered_candidates_with_tiled_windows()
        self.annotate_coverage_for_candidates()
        self.get_regions()


    # d1
    def find_gaps_in_aligned_reads(self):
        """
        Parse CIGAR string of aligned reads for insertions and deletions.

        """
        gaps_in_aligned_reads = os.path.join(self.out_dir, "gaps_in_aligned_reads")
        out_bed = os.path.join(gaps_in_aligned_reads, '0.bed')
        check_dir_exist(gaps_in_aligned_reads)

        if os.path.exists(self.alignments) and os.path.exists(self.reference) and os.stat(self.alignments).st_size > 0 and os.stat(self.reference).st_size > 0:
            cmd = 'samtools view -F 0x4 -q %d %s | python %s/scripts/PrintGaps.py %s /dev/stdin --tsd 0 --condense 20 > %s' % (
            self.mapping_quality, self.alignments, self.root_dir, self.reference, out_bed)
            os.system(cmd)
        else:
            myprint("alignments or reference is not exist!")
            sys.exit()

        step = 1
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d2
    def find_hardstops_in_aligned_reads(self):
        """
        Parse CIGAR string of aligned reads for clipped alignments.

        """

        hardstops_in_aligned_reads = os.path.join(self.out_dir, 'hardstops_in_aligned_reads')
        out_bed = os.path.join(hardstops_in_aligned_reads, '0.bed')
        check_dir_exist(hardstops_in_aligned_reads)

        cmd = "%s/scripts/mcst/hardstop %s %d %d %s " % (
            self.root_dir, self.alignments, self.mapping_quality, self.min_clipping, out_bed) + '\n'
        cmd += "sort -k 1,1 -k 2,2n -o %s %s;" % (out_bed, out_bed) + '\n'
        os.system(cmd)

        step = 2
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d3
    def calculate_coverage_per_batch(self):
        """
        Calculate coverage from each batch.

        """
        coverage = os.path.join(self.out_dir, 'coverage')
        out_bed = os.path.join(coverage, '0.bed')
        check_dir_exist(coverage)

        cmd = "%s/scripts/mcst/coverage %s -in %s;" % (self.root_dir, out_bed, self.alignments) + '\n'
        os.system(cmd)

        step = 3
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d4
    def identify_gaps_in_reference_assembly(self):
        """
        Find gap bases in the reference assembly to exclude from hardstop collection.

        """
        out_bed = os.path.join(self.out_dir, 'gaps_in_reference_assembly.bed')

        cmd = 'python %s/scripts/find_fasta_gaps.py %s > %s;' % (self.root_dir, self.reference, out_bed) + '\n'
        os.system(cmd)

        step = 4
        finish_cmd = 'echo %d_finish' % step
        only_check_exist(out_bed, finish_cmd)


    # d5
    def classify_gaps_in_aligned_reads(self, type):
        """
         Classify insertions and deletions into their own output files.

        """
        gaps_in_aligned_reads = os.path.join(self.out_dir, 'gaps_in_aligned_reads')
        input_bed = os.path.join(gaps_in_aligned_reads, '0.bed')
        aligned_reads_dir = os.path.join(self.out_dir, "aligned_reads_%s/%s" % (type, type))
        out_bed = os.path.join(aligned_reads_dir, '0.bed')
        check_dir_exist(aligned_reads_dir)

        cmd = "awk '$4 == \"%s\"' %s | sort -k 1,1 -k 2,2n > %s;" % (type, input_bed, out_bed) + '\n'
        check_file_exist(input_bed, cmd)

        step = 5
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # d6
    def merge_coverage_per_batch(self):
        """
        Collect coverages from all alignments.

        """
        coverage = os.path.join(self.out_dir, 'coverage')
        input_bed = os.path.join(coverage, '0.bed')
        out_bed = os.path.join(self.out_dir, 'coverage.bed')

        cmd = "paste %s | awk '" % input_bed + 'OFS="\\t" {{ sum = 0; for (i = 4; i <= NF; i++) {{ if (i % 4 == 0) {{ sum += $i }} }} print $1,$2,$3,sum }}\'' + '| sort -k 1,1 -k 2,2n > %s;' % out_bed + '\n'
        check_file_exist(input_bed, cmd)

        step = 6
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d7
    def merge_gap_support_from_aligned_reads(self, type):
        """
        Merge gap support for each type of event.

        """
        aligned_reads_dir = os.path.join(self.out_dir, 'aligned_reads_%s/%s' % (type, type))
        input_bed = os.path.join(aligned_reads_dir, '0.bed')
        out_bed = os.path.join(self.out_dir, 'merged_support_for_%s.bed' % type)

        cmd = "set -o pipefail; sort -k 1,1 -k 2,2n -m %s | python %s/scripts/PrintGapSupport.py /dev/stdin /dev/stdout | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4n -k 5,5n -k 6,6 -k 7,7 -k 8,8 -k 9,9 > %s;" % (
            input_bed, self.root_dir, out_bed) + '\n'
        check_file_exist(input_bed, cmd)

        step = 7
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # d8
    def annotate_coverage_of_merged_gap_support(self, type):
        """
        Annotate merged gap support with alignment coverage.

        """

        input_support = os.path.join(self.out_dir, 'merged_support_for_%s.bed' % type)
        input_coverage = os.path.join(self.out_dir, 'coverage.bed')
        out_bed = os.path.join(self.out_dir, 'coverage_and_merged_support_for_%s.bed' % type)

        cmd = "set -e; set -o pipefail; bedtools intersect -a %s -b %s -sorted -wao | awk 'OFS=\"\\t\" {{ if ($13 == \".\") {{ $13 = 0 }} print }}' | cut -f 1-6,8- | groupBy -i stdin -g 1,2,3,4,5,6,7,8 -c 12 -o mean > %s;" % (input_support, input_coverage, out_bed) + '\n'
        os.system(cmd)

        step = 8
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # d9
    def collect_hardstops(self):
        """
        Collect gaps in one command

        """
        input_bed = os.path.join(self.out_dir, 'hardstops_in_aligned_reads/0.bed')
        out_bed = os.path.join(self.out_dir, 'hardstops.bed')

        cmd = "sort -k 1,1 -k 2,2n -m %s > %s" % (input_bed, out_bed) + '\n'
        os.system(cmd)

        step = 9
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d10
    def create_hardstop_breakpoints(self):
        """
        Create hardstop breakpoints from hardstop locations (either left, right, or
        both). If a read is clipped on the left, use its start position as the
        breakpoint. If it is clipped on the right, use its end position. If a read is
        clipped from both sides, print two separate breakpoints using these same rules
        for left and right breakpoints.

        """
        input_bed = os.path.join(self.out_dir, 'hardstops.bed')
        out_bed = os.path.join(self.out_dir, 'hardstop_breakpoints.bed')

        cmd = "awk 'OFS=\"\\t\" {{ if ($6 == \"left\") {{ print $1,$2,$2 + 1,$4,$7,\"left\" }} else if ($6 == \"right\") {{ print $1,$3 - 1,$3,$4,$8,\"right\" }} else if ($6 == \"both\") {{ print $1,$2,$2 + 1,$4,$7,\"left\"; print $1,$3 - 1,$3,$4,$8,\"right\" }} }}' %s | sort -k 1,1 -k 2,2n > %s;" % (
            input_bed, out_bed) + '\n'
        os.system(cmd)

        step = 10
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d11
    def filter_candidates(self, type):
        """
        Filter candidates by support and coverage.

        """
        input_bed = os.path.join(self.out_dir, 'coverage_and_merged_support_for_%s.bed' % type)
        out_bed = os.path.join(self.out_dir, 'filtered_candidates_for_%s.bed' % type)

        cmd = "awk '$4 >= %d && $5 >= %d && $5 <= %d && $9 >= %d && $9 <= %d' %s > %s.tmp" % (
            self.min_length, self.min_support, self.max_support, self.min_coverage,
            self.max_coverage, input_bed, out_bed) + '\n'
        cmd += "bedtools merge -i %s.tmp -d 1 -c 6,4,5 -o distinct,mean,mean > %s" % (out_bed, out_bed) + '\n'
        cmd += 'rm -f %s.tmp' % out_bed + '\n'
        os.system(cmd)

        step = 11
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # d12
    def combine_filtered_candidates(self):
        """
        Summarize filtered candidates by event attributes.

        """
        input_del = os.path.join(self.out_dir, 'filtered_candidates_for_deletion.bed')
        input_ins = os.path.join(self.out_dir, 'filtered_candidates_for_insertion.bed')
        out_bed = os.path.join(self.out_dir, 'filtered_candidates.tab')

        cmd = 'sort -k 1,1 -k 2,2n -m %s %s | cut -f 1-4,6 > %s;' % (input_del, input_ins, out_bed) + '\n'
        os.system(cmd)

        step = 12
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d13
    def create_genomic_bins(self):
        input_ref = '%s.fai' % self.reference
        out_bed = os.path.join(self.out_dir, 'hardstop_bins.bed')

        cmd = 'bedtools makewindows -g %s -w %d | sort -k 1,1 -k 2,2n > %s;' % (input_ref, self.bin_size, out_bed) + '\n'
        os.system(cmd)

        step = 13
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d14.
    def count_hardstops_per_genomic_bin(self):
        """
        Count hardstops per genomic bin.

        """
        input_hardstops = os.path.join(self.out_dir, 'hardstop_breakpoints.bed')
        input_bins = os.path.join(self.out_dir, 'hardstop_bins.bed')
        out_bed = os.path.join(self.out_dir, 'hardstops_per_bin.bed')

        cmd = "bedtools intersect -a %s -b %s -sorted -c > %s;" % (input_bins, input_hardstops, out_bed) + '\n'
        os.system(cmd)

        step = 14
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d15
    def filter_and_merge_adjacent_hardstop_bins(self):
        """
        Filter hardstops to exclude reference gaps and small SV candidates. Then merge
        adjacent bins that pass filters.

        """
        input_hardstops = os.path.join(self.out_dir, "hardstops_per_bin.bed")
        input_gaps = os.path.join(self.out_dir, "gaps_in_reference_assembly.bed")
        input_small_svs = os.path.join(self.out_dir, "filtered_candidates.tab")
        out_bed = os.path.join(self.out_dir, 'merged_hardstops_per_bin.bed')

        cmd = "awk '$4 > %d' %s | bedtools window -w 1000 -a stdin -b %s -v | bedtools window -w 1000 -a stdin -b %s -v | bedtools merge -i stdin -d 1 > %s;" % (
            self.min_hardstop_support, input_hardstops, input_gaps, input_small_svs, out_bed) + '\n'
        os.system(cmd)

        step = 15
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d16
    def merge_filtered_candidates(self):
        """
        Merge filtered candidates.

        """
        input_ref = '%s.fai' % self.reference
        input_tab = os.path.join(self.out_dir, 'filtered_candidates.tab')
        input_bed = os.path.join(self.out_dir, 'merged_hardstops_per_bin.bed')
        out_bed = os.path.join(self.out_dir, 'assembly_candidates.bed')

        cmd = "set -o pipefail; cut -f 1-3 %s %s | sort -k 1,1 -k 2,2n | bedtools merge -i stdin -d %d | bedtools slop -i stdin -g %s -b %d > %s;" % (
            input_tab, input_bed, self.merge_distance, input_ref, self.slop, out_bed) + '\n'
        os.system(cmd)

        step = 16
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d17
    def create_windows_for_tiled_assemblies(self):
        """
        Windows for tiled assemblies.

        """
        input_ref = '%s.fai' % self.reference
        out_bed = os.path.join(self.out_dir, 'windows_for_tiled_assembly.bed')

        cmd = "bedtools makewindows -g %s -w %d -s %d | sort -k 1,1 -k 2,2n > %s;" % (
            input_ref, self.assembly_window_size, self.assembly_window_slide, out_bed) + '\n'
        os.system(cmd)

        step = 17
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d18
    def merge_filtered_candidates_with_tiled_windows(self):
        """
        Merge filtered candidates with tiled windows.

        """
        input_bed_1 = os.path.join(self.out_dir, 'assembly_candidates.bed')
        input_bed_2 = os.path.join(self.out_dir, 'windows_for_tiled_assembly.bed')
        out_bed = os.path.join(self.out_dir, 'assembly_candidates_and_windows.bed')

        cmd = "sort -k 1,1 -k 2,2n -m %s %s > %s;" % (input_bed_1, input_bed_2, out_bed) + '\n'
        os.system(cmd)

        step = 18
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d19
    def annotate_coverage_for_candidates(self):
        """
        Annotate assembly candidates with coverage.

        """
        input_candidates = os.path.join(self.out_dir, 'assembly_candidates_and_windows.bed')
        input_coverage = os.path.join(self.out_dir, '%s/coverage.bed')
        out_bed = os.path.join(self.out_dir, '%s/assembly_candidates_with_coverage.bed')

        cmd = "bedtools intersect -a %s -b %s -sorted -wao | awk 'OFS=\"\\t\" {{ if ($7 == \".\") {{ $7 = 0 }} print }}' | groupBy -i stdin -g 1,2,3 -c 7 -o mean > %s;" % (
            input_candidates, input_coverage, out_bed) + '\n'
        os.system(cmd)

        step = 19
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # d20
    def get_regions(self):
        """
        Annotate assembly candidates with coverage.

        """
        input_bed = os.path.join(self.out_dir, "assembly_candidates_with_coverage.bed")
        out_bed = os.path.join(self.out_dir, 'candidates.bed')

        cmd = "awk '$4 >= %s && $4 <= %s && $3 - $2 <= %s' %s > %s" % (
            self.min_coverage, self.max_coverage, self.max_candidate_length, input_bed, out_bed) + '\n'
        os.system(cmd)

        step = 20
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)



class Assemble(threading.Thread):
    def __init__(self, task_id, stand_region, root_dir, out_dir, reference, alignments):
        super(Assemble, self).__init__()
        self.stand_region = stand_region
        self.task_id = task_id
        self.root_dir = root_dir
        self.out_dir = out_dir
        self.reference = reference
        self.n_thread = 4
        self.alignments = alignments

        self.mapping_quality = 30
        self.read_length = 1000
        self.partitions = 50
        self.max_runtime = "10m"
        self.celera_spec = "/share/home/chuanlex/xieshangqian/software/smrtsv/pacbio_variant_caller/celera/pacbio.local_human.spec"
        self.bin_dir = "/share/home/chuanlex/xieshangqian/software/smrtsv/pacbio_variant_caller/bin"
        self.region_size = 60000
        self.alignment_parameters = "-affineAlign -affineOpen 8 -affineExtend 0 -bestn 1 -maxMatch 30 -sdpTupleSize 13"




    def run(self):

        self.run_assemble_cmd(self.task_id, self.stand_region)



    def run_assemble_cmd(self, task_id, stand_region):
        print('%s------>%d' % (stand_region, task_id))
        myprint('%s is assembling!' % stand_region)
        mhap_assembly = os.path.join(self.out_dir, 'mhap_assembly')
        local_assemblies = os.path.join(self.out_dir, 'local_assemblies')

        region = stand_region.replace(':', '-')
        chr = region.split('-')[0]
        chr_dir = os.path.join(mhap_assembly, chr)
        if os.path.exists(chr_dir):
            myprint('%s has been exists!' % chr_dir)
        else:
            os.mkdir(chr_dir)

        # 创建region_dir(mhap_assembly/chrI/chrI-100000-160000)
        region_dir = os.path.join(chr_dir, region)
        check_dir_exist(region_dir)

        # 定义local_assemblies目录下的local_assembly_alignments.chrI.bam
        base_filename = 'local_assembly_alignments.{}.list.txt'.format(chr)
        list_filename = os.path.join(local_assemblies, base_filename)
        tmp = os.path.join(self.out_dir, 'tmp_%s' % task_id)
        if os.path.exists(tmp):
            os.system('rm -rf %s/*' % tmp)
        else:
            os.mkdir(tmp)

        # 1,2,3,4,5,6,7,8,,9,10,11,12,13
        self.create_reference_region(tmp, region)
        self.extract_reference_sequence(tmp)
        self.get_reads(tmp, self.stand_region)
        self.convert_reads_to_bam(tmp)
        self.convert_reads_to_fasta(tmp)
        self.convert_reads_to_fastq(tmp)
        self.assemble_reads(tmp, region_dir, region)  # 耗时
        self.index_assembly(tmp)
        self.map_reads_to_assembly(tmp)
        self.convert_assembly_alignments_to_hdf5(tmp)
        self.quiver_assembly(tmp, region)
        self.trim_consensus(tmp)
        self.align_consensus_to_reference_region(tmp, region_dir)
        # 如果mhap_assembly/chrI/chrI-100000-160000/consensus_reference_alignment.sam文件已生成，
        # 将绝对路径写进local_assemblies/local_assembly_alignments.chrI.list.txt文件中
        consensus_sam = os.path.join(region_dir, 'consensus_reference_alignment.sam')
        if os.path.exists(consensus_sam):
            with open(list_filename, 'a') as fl:
                fl.write(consensus_sam + '\n')
        os.system('rm -rf %s' % tmp)


    # a1
    def create_reference_region(self, tmp, region):
        out_bed = os.path.join(tmp, "reference_region.bed")

        cmd = 'echo %s | sed "s/-/\t/g" > %s;' % (region, out_bed) + '\n'
        os.system(cmd)

        step = 1
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # a2
    def extract_reference_sequence(self, tmp):
        region_bed = os.path.join(tmp, "reference_region.bed")
        out_fasta = os.path.join(tmp, "reference_region.fasta")

        cmd = 'bedtools getfasta -fi %s -bed %s -fo %s;' % (self.reference, region_bed, out_fasta) + '\n'
        os.system(cmd)

        step = 2
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_fasta, finish_cmd)


    # a3
    def get_reads(self, tmp, stand_region):
        out_sam = os.path.join(tmp, "reads.sam")

        cmd = 'samtools view -H %s > %s;' % (self.alignments, out_sam) + '\n'
        cmd += 'samtools view -q %d %s %s >> %s;' % (
            self.mapping_quality, self.alignments, stand_region, out_sam) + '\n'
        os.system(cmd)

        step = 3
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_sam, finish_cmd)


    # a4
    def convert_reads_to_bam(self, tmp):
        input_sam = os.path.join(tmp, "reads.sam")
        out_bam = os.path.join(tmp, "reads.bam")

        cmd = 'samtools view -buS %s > %s;' % (input_sam, out_bam) + '\n'
        os.system(cmd)

        step = 4
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bam, finish_cmd)


    # a5
    def convert_reads_to_fasta(self, tmp):
        input_bam = os.path.join(tmp, "reads.bam")
        out_fasta = os.path.join(tmp, "reads.fasta")

        cmd = "samtools view %s | awk '{{ print \">\"$1; print $10 }}' | python2 %s/scripts/FormatFasta.py --fakename > %s;" % (input_bam, self.root_dir, out_fasta) + '\n'
        os.system(cmd)

        step = 5
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_fasta, finish_cmd)


    # a6
    def convert_reads_to_fastq(self, tmp):
        input_fasta = os.path.join(tmp, "reads.fasta")
        out_fastq = os.path.join(tmp, "reads.fastq")

        cmd = 'python2 %s/scripts/FastaToFakeFastq.py %s %s;' % (self.root_dir, input_fasta, out_fastq) + '\n'
        os.system(cmd)

        step = 6
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_fastq, finish_cmd)


    # a7
    def assemble_reads(self, tmp, region_dir, region):
        input_fastq = os.path.join(tmp, "reads.fastq")
        out_fasta = os.path.join(tmp, "assembly.fasta")
        assembly_log = os.path.join(self.out_dir, "assembly.log")
        region_log = os.path.join(region_dir, "assembly.log")
        assembly_output = os.path.join(tmp, "local/9-terminator/asm.ctg.fasta")
        unitig_output = os.path.join(tmp, "local/9-terminator/asm.utg.fasta")

        step = 7
        assembly_exists = False
        try:
            cmd1 = "cd %s;timeout %s %s/PBcR -threads %d -length %d -partitions %d -l local -s %s -fastq %s genomeSize=%d assembleMinCoverage=5 &> %s" % (tmp,self.max_runtime, self.bin_dir, self.n_thread, self.read_length, self.partitions, self.celera_spec, input_fastq, self.region_size, region_log) + '\n'
            cmd1 += "echo %d_try" % step + '\n'
            os.system(cmd1)

        except:
            cmd2 = 'echo -e "%s\tassembly_crashed" >> %s' % (region, assembly_log) + '\n'
            cmd2 += "echo %d_except" % step + '\n'
            os.system(cmd2)

        if os.path.exists(assembly_output) and os.stat(assembly_output).st_size > 0:
            cmd3 = "cat %s > %s" % (assembly_output, out_fasta) + '\n'
            cmd3 += 'echo -e "%s\tassembly_exists" >> %s' % (region, assembly_log) + '\n'
            assembly_exists = True
            cmd3 += "echo %d_if" % step + '\n'
            os.system(cmd3)

        elif os.path.exists(unitig_output) and os.stat(unitig_output).st_size > 0:
            cmd4 = "cat %s > %s" % (unitig_output, out_fasta) + '\n'
            cmd4 += 'echo -e "%s\tunitig_assembly_exists" >> %s' % (region, assembly_log) + '\n'
            assembly_exists = True
            cmd4 += "echo %d_elif" % step + '\n'
            os.system(cmd4)

        else:
            cmd5 = 'echo -e "%s\tno_assembly_exists" >> %s' % (region, assembly_log) + '\n'
            cmd5 += "echo %d_else" % step + '\n'
            os.system(cmd5)

        # Create an empty assembly for failed regions.
        if not assembly_exists:
            cmd6 = 'echo -e ">%s\nN" > %s' % (region, out_fasta) + '\n'
            cmd6 += "echo %d_if_not" % step + '\n'
            os.system(cmd6)

        finish_cmd = 'echo %d_finish' % step
        only_check_exist(out_fasta, finish_cmd)


    # a8
    def index_assembly(self, tmp):
        input_fasta = os.path.join(tmp, "assembly.fasta")
        out_fai = "%s.fai" % input_fasta

        cmd = "samtools faidx %s;" % input_fasta + '\n'
        os.system(cmd)

        step = 8
        finish_cmd = 'echo %d_finish' % step
        only_check_exist(out_fai, finish_cmd)


    # a9
    def map_reads_to_assembly(self, tmp):
        input_fastas = "%s/reads.fasta %s/assembly.fasta" % (tmp, tmp)
        out_sam = os.path.join(tmp, "alignment.sam")

        cmd = "blasr %s -sam -bestn 1 -out /dev/stdout -nproc %d | samtools view -h -F 0x4 -S - > %s;" % (
            input_fastas, self.n_thread, out_sam) + '\n'
        os.system(cmd)

        step = 9
        finish_cmd = 'echo %d_finish' % step
        only_check_exist(out_sam, finish_cmd)


    # a10
    def convert_assembly_alignments_to_hdf5(self, tmp):
        input_sam = os.path.join(tmp, "alignment.sam")
        out_h5 = os.path.join(tmp, "alignment.cmp.h5")

        cmd = "cp %s %s;" % (input_sam, out_h5) + '\n'
        os.system(cmd)

        step = 10
        finish_cmd = 'echo %d_finish' % step
        only_check_exist(out_h5, finish_cmd)


    # a11
    def quiver_assembly(self, tmp, region):
        assembly_log = os.path.join(self.out_dir, "assembly.log")
        input_fasta = os.path.join(tmp, "assembly.fasta")
        out_fasta = os.path.join(tmp, "consensus.fasta")

        cmd = "echo -e '%s\tquiver_failed' >> %s" % (region, assembly_log) + '\n'
        cmd += "cat %s > %s" % (input_fasta, out_fasta) + '\n'
        os.system(cmd)

        step = 11
        finish_cmd = 'echo %d_finish' % step
        only_check_exist(out_fasta, finish_cmd)


    # a12
    def trim_consensus(self, tmp):
        input_fasta = os.path.join(tmp, "consensus.fasta")
        out_fasta = os.path.join(tmp, "consensus.trimmed.fasta")

        cmd = " python2 %s/scripts/trim_lowercase.py %s %s;" % (self.root_dir, input_fasta, out_fasta) + '\n'
        os.system(cmd)

        step = 12
        finish_cmd = 'echo %d_finish' % step
        only_check_exist(out_fasta, finish_cmd)


    # a13
    def align_consensus_to_reference_region(self, tmp, region_dir):
        consensus_fasta = os.path.join(tmp, "consensus.trimmed.fasta")
        region_fasta = os.path.join(tmp, "reference_region.fasta")
        out_sam = os.path.join(region_dir, "consensus_reference_alignment.sam")

        cmd = 'blasr %s %s -clipping subread -out /dev/stdout -sam %s | samtools view -q %d - | awk \'OFS="\\t" {{ sub(/:/, "-", $3); num_of_pieces=split($3, pieces, "-"); $3 = pieces[1]; $4 = pieces[2] + $4; print }}\' | sed \'s/RG:Z:\w\+\\t//\' > %s;' % (consensus_fasta, region_fasta, self.alignment_parameters, self.n_thread, out_sam) + '\n'
        os.system(cmd)

        step = 13
        finish_cmd = 'echo %d_finish' % step
        only_check_exist(out_sam, finish_cmd)







class Call(object):
    def __init__(self, root_dir, out_dir, reference, alignments, n_thread, local_assembly_alignments):
        self.root_dir = root_dir
        self.out_dir = out_dir
        self.reference = reference
        self.alignments = alignments
        self.n_thread = n_thread
        self.local_assembly_alignments = local_assembly_alignments

        self.tsd_length = 20
        self.indel_pack_distance = 20
        self.reference_window = 5000
        self.min_contig_length = 40000
        self.window = 20
        self.overlap = 0.5
        self.bin_dir = "/share/home/chuanlex/xieshangqian/software/smrtsv/pacbio_variant_caller/bin"
        self.sample = "UCSF_Yeast9464"
        self.species = "yeast"



    def run_call_cmd(self):
        sv_calls = os.path.join(self.out_dir, "sv_calls")
        check_dir_exist(sv_calls)
        indel_calls = os.path.join(self.out_dir, "indel_calls")
        check_dir_exist(indel_calls)

        d = Detect(self.root_dir, self.out_dir, self.reference, self.alignments)
        d.calculate_coverage_per_batch()
        d.merge_coverage_per_batch()
        
        # 1,2,3,4,5,6,7,8
        self.find_calls_by_gaps_in_alignments()
        self.find_indel_gaps_in_alignments()
        self.find_inversions()
        self.calculate_coverage_from_assembled_contigs()
        self.identify_calls_by_type("insertion")
        self.identify_calls_by_type("deletion")
        self.merge_inversions()
        self.convert_inversion_bed_to_vcf()
        self.tile_contigs_from_alignments()

        # 9,10,11,12,13,14,15,16,17(ins)
        self.create_sv_fasta("insertion")
        self.trf_mask_sv_fasta("insertion")
        self.repeatmask_sv_fasta("insertion")
        self.annotate_sv_calls_with_repeatmasker("insertion")
        self.annotate_sv_calls_with_trf("insertion")
        self.summarize_calls_by_repeat_type("insertion")
        self.filter_indel_gaps_by_tiling_path()
        self.split_indels_by_type("insertion")
        self.collect_summarized_sv_calls_within_type("insertion")

        # 9,10,11,12,13,14,16,17(del)
        self.create_sv_fasta("deletion")
        self.trf_mask_sv_fasta("deletion")
        self.repeatmask_sv_fasta("deletion")
        self.annotate_sv_calls_with_repeatmasker("deletion")
        self.annotate_sv_calls_with_trf("deletion")
        self.summarize_calls_by_repeat_type("deletion")
        self.split_indels_by_type("deletion")
        self.collect_summarized_sv_calls_within_type("deletion")

        # 18,19,20,21,22,23,24,25,26(ins)
        self.collect_all_summarized_sv_calls()
        self.convert_sv_bed_to_vcf()
        self.remove_indel_events_in_homopolymers("insertion")
        self.filter_indel_events_by_size("insertion")
        self.annotate_coverage_of_pacbio_reads_for_indels("insertion")
        self.annotate_coverage_of_assembled_contigs_for_indels("insertion")
        self.cut_indel_events_by_columns("insertion")
        self.annotate_support_from_assembled_contigs_for_indels("insertion")
        self.combine_annotations_for_indel_type("insertion")

        # 20,21,22,23,24,25,26(del)
        self.remove_indel_events_in_homopolymers("deletion")
        self.filter_indel_events_by_size("deletion")
        self.annotate_coverage_of_pacbio_reads_for_indels("deletion")
        self.annotate_coverage_of_assembled_contigs_for_indels("deletion")
        self.cut_indel_events_by_columns("deletion")
        self.annotate_support_from_assembled_contigs_for_indels("deletion")
        self.combine_annotations_for_indel_type("deletion")

        # 27,28,29
        self.call_indels()
        self.convert_indel_bed_to_vcf()
        self.call_variants()



    # c1
    def find_calls_by_gaps_in_alignments(self):
        out_bed = os.path.join(self.out_dir, "sv_calls/gaps.bed")

        cmd = "samtools view -h %s | %s/scripts/PrintGaps.py %s /dev/stdin --qpos --condense %d --tsd %d | sort -k 1,1 -k 2,2n > %s;" % (
            self.local_assembly_alignments, self.root_dir, self.reference, self.indel_pack_distance, self.tsd_length, out_bed) + '\n'
        os.system(cmd)

        step = 1
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # c2
    def find_indel_gaps_in_alignments(self):
        out_bed = os.path.join(self.out_dir, "indel_calls/gaps.bed")
        self.indel_pack_distance = 0

        cmd = "samtools view %s | %s/scripts/PrintGaps.py %s /dev/stdin --minLength 30 --maxLength 50 --context 6 --removeAdjacentIndels --onTarget --minContigLength %d --condense %d --outFile %s;" % (
            self.local_assembly_alignments, self.root_dir, self.reference, self.min_contig_length,
            self.indel_pack_distance, out_bed) + '\n'
        os.system(cmd)

        step = 2
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # c3
    def find_inversions(self):
        out_bed = os.path.join(self.out_dir, "sv_calls/inversions.bed")

        cmd = "samtools view %s | %s/scripts/mcst/screenInversions /dev/stdin %s %s -w %d -r --noClip -j %d;" % (self.local_assembly_alignments, self.root_dir, self.reference, out_bed, self.reference_window,
            self.n_thread) + '\n'
        os.system(cmd)

        step = 3
        finish_cmd = 'echo %d_finish' % step
        only_check_exist(out_bed, finish_cmd)


    # c4
    def calculate_coverage_from_assembled_contigs(self):
        out_bed = os.path.join(self.out_dir, "assembled_contigs.depth.bed")

        cmd = "bedtools bamtobed -i %s | %s/scripts/BedIntervalsToDepth.py /dev/stdin %s --out /dev/stdout | sort -k 1,1 -k 2,2n > %s;" % (self.local_assembly_alignments, self.root_dir, self.reference, out_bed) + '\n'
        os.system(cmd)

        step = 4
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # c5
    def identify_calls_by_type(self, type):
        input_bed = os.path.join(self.out_dir, "sv_calls/gaps.bed")
        out_bed = os.path.join(self.out_dir, "sv_calls/calls.%s.bed" % type)

        if type == "insertion":
            call_comparison_action = "window"
        else:
            call_comparison_action = "intersect"

        cmd = 'awk \'$4 == "%s" && index($6, "N") == 0\' %s | awk \'OFS="\\t" {{ if ("%s" == "insertion") {{ $3=$2 + 1 }} print }}\' | python %s/scripts/cluster_calls.py --window %d --reciprocal_overlap %f /dev/stdin %s | awk \'OFS="\\t" {{ if ("%s" == "insertion") {{ $3=$2 + $5 }} print }}\' | sort -k 1,1 -k 2,2n | while read line; do set -- $line; coverage=`samtools view -c %s $1:$2-$3`; echo -e "$line\\t$coverage"; done > %s;' % (type, input_bed, type, self.root_dir, self.window, self.overlap, call_comparison_action, type,self.local_assembly_alignments, out_bed) + '\n'
        os.system(cmd)

        step = 5
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # c6
    def merge_inversions(self):
        inversions = os.path.join(self.out_dir, "sv_calls/inversions.bed")
        out_bed = os.path.join(self.out_dir, "sv_calls/merged_inversions.bed")

        cmd = 'sort -k 1,1 -k 2,2n %s | bedtools merge -i stdin -d 0 -c 4 -o count | while read line; do set -- $line; coverage=`samtools view -c %s $1:$2-$3`; echo -e "$line\\t$coverage"; done | awk \'OFS="\\t" {{ print $1,$2,$3,"inversion",$4,$5 }}\' > %s' % (inversions, self.local_assembly_alignments, out_bed) + '\n'
        os.system(cmd)

        step = 6
        finish_cmd = 'echo %d_finish' % step
        only_check_exist(out_bed, finish_cmd)


    # c7
    def convert_inversion_bed_to_vcf(self):
        input_bed = os.path.join(self.out_dir, "sv_calls/merged_inversions.bed")
        out_vcf = os.path.join(self.out_dir, "inversions.vcf")

        cmd = "%s/scripts/variants_bed_to_vcf.py %s %s %s %s inversion;" % (
            self.root_dir, input_bed, self.reference, out_vcf, self.sample) + '\n'
        os.system(cmd)

        step = 7
        finish_cmd = 'echo %d_finish' % step
        only_check_exist(out_vcf, finish_cmd)



    # c8
    def tile_contigs_from_alignments(self):
        out_bed = os.path.join(self.out_dir, "tiling_contigs.tab")

        cmd = "samtools view -h %s | %s/scripts/TilingPath.py /dev/stdin > %s;" % (
            self.local_assembly_alignments, self.root_dir, out_bed) + '\n'
        os.system(cmd)

        step = 8
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # c9
    def create_sv_fasta(self, type):
        input_bed = os.path.join(self.out_dir, "sv_calls/calls.%s.bed" % type)
        out_bed = os.path.join(self.out_dir, "sv_calls/%s/%s.fasta" % (type, type))
        out_dir = os.path.join(self.out_dir, "sv_calls/%s" % type)
        check_dir_exist(out_dir)

        cmd = "%s/scripts/GapBedToFasta.py %s %s;" % (self.root_dir, input_bed, out_bed) + '\n'
        os.system(cmd)
        
        step = 9
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # c10
    def trf_mask_sv_fasta(self, type):
        input_bed = os.path.join(self.out_dir, "sv_calls/%s/%s.fasta" % (type, type))
        out_trf = os.path.join(self.out_dir, "sv_calls/%s/rm/%s.fasta.trf" % (type, type))
        out_dir = os.path.join(self.out_dir, "sv_calls/%s/rm" % type)
        check_dir_exist(out_dir)

        cmd = "cd %s; %s/trf %s 2 7 7 80 10 20 500 -m -ngs -h > %s;" % (self.out_dir, self.bin_dir, input_bed, out_trf) + '\n'
        os.system(cmd)

        step = 10
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_trf, finish_cmd)


    # c11
    def repeatmask_sv_fasta(self, type):
        input_bed = os.path.join(self.out_dir, "sv_calls/%s/%s.fasta" % (type, type))
        out = os.path.join(self.out_dir, "sv_calls/%s/rm/%s.fasta.out" % (type, type))

        cmd = """RepeatMasker -species "%s" -dir `dirname %s` -xsmall -no_is -s -pa %s %s;""" % (
            self.species, out, self.n_thread, input_bed) + '\n'
        os.system(cmd)

        step = 11
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out, finish_cmd)


    # c12
    def annotate_sv_calls_with_repeatmasker(self, type):
        calls = os.path.join(self.out_dir, "sv_calls/calls.%s.bed" % type)
        repeats = os.path.join(self.out_dir, "sv_calls/%s/rm/%s.fasta.out" % (type, type))
        masked_fasta = os.path.join(self.out_dir, "sv_calls/%s/rm/%s.fasta.masked" % (type, type))
        out_bed = os.path.join(self.out_dir, "sv_calls/all_annotated.%s.bed" % type)

        cmd = "%s/scripts/AnnotateGapBed.py %s %s %s %s;" % (self.root_dir, calls, out_bed, repeats, masked_fasta) + '\n'
        os.system(cmd)

        step = 12
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # c13
    def annotate_sv_calls_with_trf(self, type):
        input_bed = os.path.join(self.out_dir, "sv_calls/all_annotated.%s.bed" % type)
        input_trf = os.path.join(self.out_dir, "sv_calls/%s/rm/%s.fasta.trf" % (type, type))
        out_bed = os.path.join(self.out_dir, "sv_calls/all_annotated_with_trf.%s.bed" % type)

        cmd = "%s/scripts/AnnotateWithTRF.py %s %s %s;" % (self.root_dir, input_bed, input_trf, out_bed) + '\n'
        os.system(cmd)

        step = 13
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # c14
    def summarize_calls_by_repeat_type(self, type):
        input_bed = os.path.join(self.out_dir, "sv_calls/all_annotated_with_trf.%s.bed" % type)
        out_dir = os.path.join(self.out_dir, "sv_calls/summarized_%s" % type)
        check_dir_exist(out_dir)

        cmd = "cd %s; awk '$20 > 0.8' %s > %s/TRF.bed; " % (self.out_dir, input_bed, out_dir) + '\n'

        cmd += "awk '$20 <= 0.8' %s > sv_calls/not_trf.bed; " % input_bed + '\n'

        cmd += "%s/scripts/FixMasked.py 17 < sv_calls/not_trf.bed | awk '$18 < 0.7' > %s/NotMasked.bed; " % (
            self.root_dir, out_dir) + '\n'

        cmd += "%s/scripts/FixMasked.py 17 < sv_calls/not_trf.bed | awk '$18 >= 0.7' > sv_calls/repeat.bed; " % self.root_dir + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py sv_calls/repeat.bed --prefix AluY --minPrefix 1 --maxPrefix 1 --maxNotPrefix 0 --maxSTR 0 --remainder %s/1.bed > %s/AluY.simple.bed; " % (
            self.root_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/1.bed --prefix AluS --minPrefix 1 --maxPrefix 1 --maxNotPrefix 0 --maxSTR 0 --remainder %s/2.bed > %s/AluS.simple.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/2.bed --minSTR 1 --maxNotPrefix 0 --remainder %s/4.bed > %s/STR.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/4.bed --prefix L1HS --maxNotPrefix 0 --remainder %s/5.bed > %s/L1HS.simple.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/5.bed --prefix Alu --minPrefix 1 --maxNotPrefix 0 --maxSTR 0 --remainder %s/6.bed > %s/Alu.Mosaic.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/6.bed --prefix Alu --minSTR 1 --minPrefix 1 --maxNotPrefix 0 --remainder %s/7.bed > %s/Alu.STR.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/7.bed --prefix ALR --minPrefix 1 --maxNotPrefix 0 --remainder %s/8.bed > %s/ALR.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/8.bed --prefix SVA --minPrefix 1 --maxNotPrefix 0 --remainder %s/9.bed > %s/SVA.simple.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/9.bed --prefix HERV --minPrefix 1 --maxNotPrefix 0 --remainder %s/10.bed > %s/HERV.simple.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/10.bed --prefix L1P --minPrefix 1 --maxNotPrefix 0 --remainder %s/11.bed > %s/L1P.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/11.bed --prefix BSR/Beta --minPrefix 1 --maxNotPrefix 0 --remainder %s/12.bed > %s/Beta.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/12.bed --prefix HSAT --minPrefix 1 --maxNotPrefix 0 --remainder %s/13.bed > %s/HSAT.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/13.bed --prefix MER --minPrefix 1 --maxNotPrefix 0 --remainder %s/14.bed > %s/MER.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/14.bed --prefix L1 --minPrefix 1 --maxNotPrefix 0 --remainder %s/15.bed > %s/L1.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/15.bed --prefix LTR --minPrefix 1 --maxNotPrefix 0 --remainder %s/16.bed > %s/LTR.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "%s/scripts/PrintUniqueEvents.py %s/16.bed --max 1 --remainder %s/17.bed > %s/Singletons.bed; " % (
            self.root_dir, out_dir, out_dir, out_dir) + '\n'

        cmd += "mv -f %s/17.bed %s/Complex.bed; " % (out_dir, out_dir) + '\n'
        cmd += 'echo 14_%s_finish' % type + '\n'
        os.system(cmd)


    # c15
    def filter_indel_gaps_by_tiling_path(self):
        input_bed = os.path.join(self.out_dir, "indel_calls/gaps.bed")
        input_tab = os.path.join(self.out_dir, "tiling_contigs.tab")
        out_bed = os.path.join(self.out_dir, "indel_calls/gaps.tiled.bed")
        log = os.path.join(self.out_dir, "indel_calls/gaps.tiled.log")
        cmd = "%s/scripts/FilterGapsByTilingPath.py %s %s > %s 2> %s;" % (
            self.root_dir, input_bed, input_tab, out_bed, log) + '\n'
        os.system(cmd)

        step = 15
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # c16
    def split_indels_by_type(self, type):
        input_bed = os.path.join(self.out_dir, "indel_calls/gaps.tiled.bed")
        out_bed = os.path.join(self.out_dir, "indel_calls/%s/gaps.bed" % type)
        out_dir = os.path.join(self.out_dir, "indel_calls/%s" % type)
        check_dir_exist(out_dir)

        cmd = "grep %s %s | sort -k 1,1 -k 2,2n > %s;" % (type, input_bed, out_bed) + '\n'
        os.system(cmd)

        step = 16
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # c17
    def collect_summarized_sv_calls_within_type(self, type):
        input_dir = os.path.join(self.out_dir, "sv_calls/summarized_%s" % type)
        out_bed = os.path.join(self.out_dir, "sv_calls/repeat_classified_%s.bed" % type)
        cmd = "for file in %s/*.bed; do repeat_type=`basename ${{file/.bed/}} | sed 's/\./_/g'`; awk -v repeat_type=$repeat_type 'OFS=\"\\t\" {{ print $0,repeat_type }}' $file; done | sort -k 1,1 -k 2,2n | uniq > %s;" % (input_dir, out_bed) + '\n'
        os.system(cmd)

        step = 17
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # c18
    def collect_all_summarized_sv_calls(self):
        input_bed = "%s/sv_calls/repeat_classified_insertion.bed %s/sv_calls/repeat_classified_deletion.bed" % (self.out_dir, self.out_dir)
        out_bed = os.path.join(self.out_dir, "sv_calls/sv_calls_with_repeats.bed")

        cmd = "sort -k 1,1 -k 2,2n -m %s | uniq > %s;" % (input_bed, out_bed) + '\n'
        os.system(cmd)

        step = 18
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_bed, finish_cmd)


    # c19
    def convert_sv_bed_to_vcf(self):
        input_bed = os.path.join(self.out_dir, "sv_calls/sv_calls_with_repeats.bed")
        out_vcf = os.path.join(self.out_dir, "sv_calls.vcf")

        cmd = "%s/scripts/variants_bed_to_vcf.py %s %s %s %s sv;" % (
            self.root_dir, input_bed, self.reference, out_vcf, self.sample) + '\n'
        os.system(cmd)

        step = 19
        finish_cmd = 'echo %d_finish' % step
        check_file_exist(out_vcf, finish_cmd)


    # c20
    def remove_indel_events_in_homopolymers(self, type):
        out_dir = os.path.join(self.out_dir, "indel_calls/%s" % type)
        input_bed = os.path.join(out_dir, "gaps.bed")
        out_bed = os.path.join(out_dir, "gaps_without_homopolymers.bed")

        base_cmd = "awk '$11 == \"F\"' %s | cut -f 1,2,3,5,6,8 | sort -k 1,1 -k 2,2n -k 4,4n | %s/scripts/PrintSNVSupport.py /dev/stdin /dev/stdout" % (
            input_bed, self.root_dir)
        if type == "insertion":
            cmd = "%s | %s/scripts/BedMod.py --leftjustify 1 /dev/stdin %s" % (base_cmd, self.root_dir, out_bed) + '\n'
            os.system(cmd)
        else:
            cmd = "%s > %s" % (base_cmd, out_bed) + '\n'
            os.system(cmd)

        step = 20
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # c21
    def filter_indel_events_by_size(self, type):
        input_bed = os.path.join(self.out_dir, "indel_calls/%s/gaps_without_homopolymers.bed" % type)
        out_bed = os.path.join(self.out_dir, "indel_calls/%s/gaps_2bp_or_more_without_homopolymers.bed" % type)
        cmd = "awk '$4 >= 2' %s | bedtools groupby -c 6 -o max -full > %s;" % (input_bed, out_bed) + '\n'
        os.system(cmd)

        step = 21
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # c22
    def annotate_coverage_of_pacbio_reads_for_indels(self, type):
        input_bed = os.path.join(self.out_dir, "indel_calls/%s/gaps_2bp_or_more_without_homopolymers.bed" % type)
        coverage = os.path.join(self.out_dir, "coverage.bed")
        out_txt = os.path.join(self.out_dir, "indel_calls/%s/read_coverage.txt" % type)

        cmd = "cut -f 1-3 %s " % input_bed
        cmd += "| bedtools intersect -a stdin -b %s -loj -sorted | sed 's/\t\./\t0/g' | " % coverage
        cmd += "bedtools groupby -c 7 -o mean -full | cut -f 8 | "
        cmd += "awk '{{ printf(\"%2.2f\\n\", $1) }}' "
        cmd += " > %s;" % out_txt + '\n'
        os.system(cmd)

        step = 22
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_txt, finish_cmd)


    # c23
    def annotate_coverage_of_assembled_contigs_for_indels(self, type):
        input_bed = os.path.join(self.out_dir, "indel_calls/%s/gaps_2bp_or_more_without_homopolymers.bed" % type)
        input_depth = os.path.join(self.out_dir, "assembled_contigs.depth.bed")
        out_txt = os.path.join(self.out_dir, "indel_calls/%s/assembled_contigs_coverage.txt" % type)

        cmd = "cut -f 1-3 %s | bedtools intersect -a stdin -b %s -loj -sorted | sed 's/\\t\./\\t0/g' | bedtools groupby -c 7 -o max -full | cut -f 8 > %s;" % (
            input_bed, input_depth, out_txt) + '\n'
        os.system(cmd)

        step = 23
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_txt, finish_cmd)

    # c24
    def cut_indel_events_by_columns(self, type):
        input_bed = os.path.join(self.out_dir, "indel_calls/%s/gaps_2bp_or_more_without_homopolymers.bed" % type)
        out_bed = os.path.join(self.out_dir, "indel_calls/%s/filtered_gaps.bed" % type)

        cmd = "cut -f 1-5 %s > %s;" % (input_bed, out_bed) + '\n'
        os.system(cmd)

        step = 24
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # c25
    def annotate_support_from_assembled_contigs_for_indels(self, type):
        input_bed = os.path.join(self.out_dir, "indel_calls/%s/gaps_2bp_or_more_without_homopolymers.bed" % type)
        out_txt = os.path.join(self.out_dir, "indel_calls/%s/support_by_assembled_contigs.txt" % type)

        cmd = "cut -f 6 %s > %s;" % (input_bed, out_txt) + '\n'
        os.system(cmd)

        step = 25
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_txt, finish_cmd)


    # c26
    def combine_annotations_for_indel_type(self, type):
        input_bed = os.path.join(self.out_dir, "indel_calls/%s/filtered_gaps.bed" % type)
        input_txt1 = os.path.join(self.out_dir, "indel_calls/%s/support_by_assembled_contigs.txt" % type)
        input_txt2 = os.path.join(self.out_dir, "indel_calls/%s/assembled_contigs_coverage.txt" % type)
        input_txt3 = os.path.join(self.out_dir, "indel_calls/%s/read_coverage.txt" % type)
        out_tab = os.path.join(self.out_dir, "indel_calls/%s.tab" % type)

        cmd = 'paste %s %s %s %s | awk \'OFS="\\t" {{ print $0,"%s" }}\' > %s;' % (
            input_bed, input_txt1, input_txt2, input_txt3, type, out_tab)
        os.system(cmd)

        step = 26
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_tab, finish_cmd)


    # c27
    def call_indels(self):
        input_tab = "%s/indel_calls/insertion.tab %s/indel_calls/deletion.tab" % (self.out_dir, self.out_dir)
        out_bed = os.path.join(self.out_dir, "indel_calls.bed")

        cmd = "sort -k 1,1 -k 2,2n %s > %s;" % (input_tab, out_bed) + '\n'
        os.system(cmd)

        step = 27
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_bed, finish_cmd)


    # c28
    def convert_indel_bed_to_vcf(self):
        input_bed = os.path.join(self.out_dir, "indel_calls.bed")
        out_vcf = os.path.join(self.out_dir, "indel_calls.vcf")

        cmd = "%s/scripts/variants_bed_to_vcf.py %s %s %s %s indel;" % (
            self.root_dir, input_bed, self.reference, out_vcf, self.sample) + '\n'
        os.system(cmd)

        step = 28
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_vcf, finish_cmd)


    # c29
    def call_variants(self):
        input_vcf1 = os.path.join(self.out_dir, "sv_calls.vcf")
        input_vcf2 = os.path.join(self.out_dir, "indel_calls.vcf")
        input_vcf3 = os.path.join(self.out_dir, "inversions.vcf")
        out_vcf = os.path.join(self.out_dir, "variants.vcf")

        cmd = 'grep "^##" %s | sed \'/INFO/d\' > %s; ' % (input_vcf1, out_vcf) + '\n'
        cmd += 'grep -h "^##" %s %s %s | grep INFO | sort | uniq >> %s; ' % (
            input_vcf1, input_vcf2, input_vcf3, out_vcf) + '\n'
        cmd += 'grep -h "^#CHROM" %s >> %s; ' % (input_vcf1, out_vcf) + '\n'
        cmd += "sed '/^#/d' %s %s %s | sort -k 1,1 -k 2,2n >> %s;" % (input_vcf1, input_vcf2, input_vcf3, out_vcf) + '\n'
        os.system(cmd)

        step = 29
        finish_cmd = 'echo %d_%s_finish' % (step, type)
        check_file_exist(out_vcf, finish_cmd)



if __name__ == '__main__':
    main()
