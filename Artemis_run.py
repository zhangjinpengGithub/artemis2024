import subprocess
import time
import os,sys
import argparse
from multiprocessing import Pool
import pandas as pd
###edit zhangjp 20240407
def run_command(cmd):
        return_code = subprocess.call(cmd, shell=True)
def run_artemis_kmers(sample,bam, outdir,cpu):
    print("\n")
    os.makedirs(outdir, exist_ok=True)
    cmd=f"samtools fastq {bam} -@ {cpu} > {outdir}/{sample}.fq"
    run_command(cmd)
    cmd=f"/dssg/home/zhangjp/bin/jellyfish count  --mer-len 24 --size 3G --threads {cpu} --canonical --quality-start=32 --if /dssg/home/zhangjp/soft/ARTEMIS/kmers_all.fasta --output {outdir}/{sample}_kmers_all.jellyfish {outdir}/{sample}.fq"
    run_command(cmd)
    
def a2_count_reads_denom(sample,bam,outdir):
    file_path=f'{outdir}/{sample}_kmers_all.jellyfish'
    if os.path.exists(file_path):
        print("kmers_all.jellyfish is found")
    else:
        output_file=outdir+'R1.jellyfish.nofound'
        with open(output_file, 'a+') as f1:
            f1.write("kmers_all_R1.jellyfish is not found" + '\n')
        error_callback("kmers_all_R1.jellyfish")
    cmd=f"samtools view -c -q 30 -F 3844 -@ 8 {bam} > {outdir}/{sample}.aligned_reads"
    run_command(cmd)
def a3_fast_aggregate(sample,Feature_file,outdir,cpu):
    runcmd=f"python /dssg/home/zhangjp/soft/ARTEMIS/newscript/TempAnalysis/artemis.py --jellyfish_R1 {outdir}/{sample}_kmers_all.jellyfish  --outdir {outdir}  --kmers_type_path /dssg/home/zhangjp/soft/ARTEMIS/final_kmers  --cpu {cpu} --jellyfish /dssg/home/zhangjp/bin/jellyfish --Feature_file {Feature_file}  --SampleID {sample} --aligned_reads {outdir}/{sample}.aligned_reads"
    print(runcmd)
    run_command(runcmd)
    if os.path.exists(Feature_file):
        run_command(f"rm -rf {outdir}")
def error_callback(e):
    print(f"An exception occurred: {e}")
    raise SystemExit  # 可以选择抛出一个异常来中断程序

def run_script():
    parser = argparse.ArgumentParser(description="This is a script for generating Artemis features")
    
    parser.add_argument("--sample", type=str, help="sample name.")
    parser.add_argument("--outdir", type=str, help="Path to the output directory.")
    parser.add_argument("--Feature_file", type=str, help="Path to the Feature file name.")
    parser.add_argument("--bam", type=str, help="Path to the  bam file.")
    parser.add_argument("--cpu", type=int,help="Number of CPUs used.")
    args = parser.parse_args()
    
    start_time = time.time()
    run_artemis_kmers(args.sample, args.bam, args.outdir,args.cpu)
    a2_count_reads_denom(args.sample,args.bam,args.outdir)
    a3_fast_aggregate(args.sample,args.Feature_file,args.outdir,args.cpu)    
    end_time = time.time()
    run_time = end_time - start_time
    print(f"Run time: {run_time} seconds")

if __name__ == "__main__":
    run_script()
