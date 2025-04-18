import time
import os,sys,re
import argparse
import pandas as pd
from joblib import Parallel, delayed
import subprocess
#def run_command(cmd):
#        return_code = subprocess.call(cmd, shell=True)
def fast_aggregate(sample,ref_type,outdir,cpu,jellyfish,R1_jellyfish):
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outdir+"/features", exist_ok=True)
    os.makedirs(outdir+"/Raw_features", exist_ok=True)
    ref_files = sorted([f for f in os.listdir(ref_type) if f.endswith(".fasta")])
    n_jobs = int(cpu) -1
    results = Parallel(n_jobs=n_jobs, verbose=0)(delayed(process_ref_file)(ref_file, outdir, sample, ref_type,jellyfish,R1_jellyfish) for ref_file in ref_files)
def process_ref_file(ref_file, outdir, sample,ref_type,jellyfish,R1_jellyfish): 
    #jellyfish=args.jellyfish
    ref_name = ref_file.replace(".fasta", "")
    refpath=ref_type+"/"+ref_file
    r1=f"{jellyfish} query {R1_jellyfish} -s '{refpath}'  | awk "+" '{ sum += $2 } END { print sum }'"
    r1out=subprocess.getoutput(r1)
    val=int(r1out)
    output_file=outdir+'/Raw_features/'+f'{ref_name}.txt'
    with open(output_file, 'w+') as f1:
        f1.write(f'{ref_name} {val} {sample}'+"\n")
def run_merge(sample,ref_type,outdir):
    ref_files = sorted([f for f in os.listdir(ref_type) if f.endswith(".fasta")])
    flag=0
    for i in ref_files:
        ref_name = i.replace(".fasta", "")
        output_file=outdir+'/Raw_features/'+f'{ref_name}.txt'
        df=pd.read_csv(output_file,header=None,sep=" ")
        if(flag==0):
            dfall=df
            flag=1
        else:
            dfall=pd.concat([dfall,df])
    dfall.to_csv(outdir+"/features/"+f'{sample}.txt',sep=" ",header=0,index=False)

def run_featureMatrix(sample,outdir,out_file,aligned_reads):
    dat=pd.read_csv(outdir+'/features/'+f'{sample}.txt',header=None,sep=" ")
    dat.columns=["feature","count","SampleID"]
    counts=int(pd.read_csv(aligned_reads,header=None).iloc[-1,0])/1000000
    dat['n_norm']=dat['count'].map(lambda x:int(x)/counts)
    dat['feature']=dat['feature'].map(lambda x:"ARTEMIS."+x)
    dat['feature'] = [re.sub("[^a-zA-Z0-9._]", "_", col) for col in dat['feature']]
    dat=dat.sort_values("feature")[['feature','n_norm']].T.reset_index(drop=True)
    features_SampleID=['SampleID']+[sample]
    dat.insert(0, '0', features_SampleID)
    dat.to_csv(out_file,index=False,header=None)
    
def run_script():
    parser = argparse.ArgumentParser(description="This is a script for generating Artemis features")

    # 添加参数
    parser.add_argument("--jellyfish_R1", type=str, help="Path to the _kmers_all_R1.jellyfish file.")
    parser.add_argument("--outdir", type=str, help="Path to the output directory.")
    parser.add_argument("--kmers_type_path", type=str, help="Path to  the  1280 types of repetitive   file.")
    parser.add_argument("--cpu", type=int,help="Number of CPUs used.")
    parser.add_argument("--jellyfish", type=str,help="Path to the  jellyfish software ")
    parser.add_argument("--Feature_file", type=str,help="Path to the Output Artemis feature file ")
    parser.add_argument("--SampleID", type=str,help="Please enter the sample name ")
    parser.add_argument("--aligned_reads", type=str,help="Path to the aligned_reads file ")
    # 解析命令行参数
    args = parser.parse_args()
    # 调用主函数
    start_time = time.time()
    fast_aggregate(args.SampleID,args.kmers_type_path,args.outdir,args.cpu,args.jellyfish,args.jellyfish_R1)
    run_merge(args.SampleID,args.kmers_type_path,args.outdir)
    run_featureMatrix(args.SampleID,args.outdir,args.Feature_file,args.aligned_reads)
    end_time = time.time()
    run_time = end_time - start_time
    print(f"Run time: {run_time} seconds")

if __name__ == "__main__":
    run_script()

