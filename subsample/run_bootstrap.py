import os
import sys
import subprocess
import multiprocessing
import time
import argparse
import shutil

# ================= CONFIGURATION =================
RAW_BAM_DIR = "/data1/guoliang/ZhuoLJ/1.depth_threshold_depth_bootstrap"
BASE_POOL_DIR = "bootstrap_results"
OUT_DIR = "bootstrap_results_py" 

REF_DIR = "/data1/guoliang/ZhuoLJ/ref"
BED_DIR = "/data1/guoliang/ZhuoLJ/getgene/share_getgene_24"
MAP_FILE = "sample_ref_map.txt"

TARGET_MAX_COV = 30
TARGET_COVS = [5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20] 
ITERATIONS = 100

# 满足你的要求：每个任务 8 核，50 个并发
THREADS_PER_TASK = 5
MAX_WORKERS = 30 
EXTRACT_WORKERS = 10

SHM_TMP_DIR = "/dev/shm/guoliang_bootstrap_tmp"
SHM_MIN_FREE_GB = 300
# =================================================

# 全局字典：主进程加载一次，子进程共享内存
GLOBAL_4D_SITES = {}

def get_disk_free_gb(path):
    try:
        st = os.statvfs(path)
        return (st.f_bavail * st.f_frsize) / (1024**3)
    except Exception:
        return 0

def get_single_base_median_cov(bam_file):
    if not os.path.exists(bam_file):
        return 0
    cmd = (
        f"samtools depth -q 20 -Q 30 {bam_file} 2>/dev/null | "
        f"awk '{{print $3}}' | "
        f"sort -n | "
        f"awk '{{a[NR]=$1}} END {{if(NR>0) print a[int(NR/2)]; else print 0}}'"
    )
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.stdout.strip().isdigit():
        return int(res.stdout.strip())
    return 0

def extract_worker(args):
    sample, ref = args
    raw_bam = f"{RAW_BAM_DIR}/{sample}.chr.bam"
    base_30x_bam = f"{BASE_POOL_DIR}/{sample}_30X.bam"
    
    print(f"[{sample}] Scanning median coverage...", end="", flush=True)
    median_cov = get_single_base_median_cov(raw_bam)
    
    if median_cov == 0:
        print(f" FAILED to get depth. Skipping.", flush=True)
        return False, sample
        
    frac = TARGET_MAX_COV / median_cov
    
    if frac > 1.0:
        print(f" Depth: {median_cov}X < {TARGET_MAX_COV}X. SKIPPING sample.", flush=True)
        if os.path.exists(base_30x_bam): os.remove(base_30x_bam)
        if os.path.exists(f"{base_30x_bam}.bai"): os.remove(f"{base_30x_bam}.bai")
        return False, sample
        
    print(f" Depth: {median_cov}X, Fraction: {frac:.4f}", flush=True)
    print(f"[{sample}] Extracting base pool...", flush=True)
    
    cmd = (
        f"samtools view -u -q 30 -f 0x2 -F 0x900 -s {frac:.4f} {raw_bam} | "
        f"samtools sort -@ {THREADS_PER_TASK} -m 500M -o {base_30x_bam} - && "
        f"samtools index -@ {THREADS_PER_TASK} {base_30x_bam}"
    )
    subprocess.run(cmd, shell=True)
    success = os.path.exists(base_30x_bam)
    if success:
        print(f"[{sample}] Base extraction complete.", flush=True)
    else:
        print(f"[{sample}] Base extraction failed.", flush=True)
    return success, sample

def extract_individual_30x_pools(my_samples):
    need_extract = []
    for sample, ref in my_samples:
        bam = f"{BASE_POOL_DIR}/{sample}_30X.bam"
        bai = f"{bam}.bai"
        if not (os.path.exists(bam) and os.path.exists(bai)):
            need_extract.append((sample, ref))
            
    if need_extract:
        print(f"\n[STEP 1] Extracting pools for {len(need_extract)} missing samples...", flush=True)
        with multiprocessing.Pool(processes=EXTRACT_WORKERS) as pool:
            for success, sample in pool.imap_unordered(extract_worker, need_extract):
                pass
    else:
        print("\n[STEP 1] All assigned base pools already exist. Skipping extraction.")

def process_task(task_args):
    sample, ref_prefix, cov, i = task_args
    base_30x_bam = f"{BASE_POOL_DIR}/{sample}_30X.bam"
    ref_fa = f"{REF_DIR}/{ref_prefix}.fna"
    interval_file = f"{BED_DIR}/{ref_prefix}.intervals.txt"
    tmp_res = f"{OUT_DIR}/tmp_res_{sample}_{cov}.txt"

    if get_disk_free_gb(SHM_TMP_DIR) < SHM_MIN_FREE_GB:
        with open(tmp_res, "a") as f:
            f.write(f"{cov}X\t{i}\tNA\n")
        return

    # 从全局内存字典拿 4D 位点
    target_positions = GLOBAL_4D_SITES.get(ref_prefix, set())

    tmp_bam = f"{SHM_TMP_DIR}/tmp_{sample}_{cov}_{i}_{os.getpid()}.bam"
    frac = cov / TARGET_MAX_COV
    seed = (os.getpid() * 32768 + i) % 2147483647
    samtools_seed = f"{seed}.{int(frac * 10000)}"

    subprocess.run(f"samtools view -@ {THREADS_PER_TASK} -b -s {samtools_seed} {base_30x_bam} > {tmp_bam} 2>/dev/null", shell=True)
    
    het_val = "NA"
    if os.path.exists(tmp_bam) and os.path.getsize(tmp_bam) > 0:
        subprocess.run(f"samtools index -@ {THREADS_PER_TASK} {tmp_bam} 2>/dev/null", shell=True)
        
        # 🌟 核心修复：在 query 中加上 %CHROM 和 %POS 让 Python 能拿到坐标
        cmd = (
            f"bcftools mpileup -B -I -A -Ou -q 30 -Q 20 -f {ref_fa} "
            f"-R {interval_file} -a FORMAT/DP --threads {THREADS_PER_TASK} {tmp_bam} 2>/dev/null | "
            f"bcftools call -cO u 2>/dev/null | "
            f"bcftools query -f '%CHROM\t%POS\t[%GT,%DP\\n]' 2>/dev/null"
        )
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        output = res.stdout.strip()

        sites = 0
        het = 0
        for line in output.split('\n'):
            if not line: continue
            cols = line.split('\t')
            if len(cols) < 3: continue
            
            chrom = cols[0]
            pos = cols[1]
            gt_dp = cols[2].split(',')
            if len(gt_dp) < 2: continue
            
            gt, dp = gt_dp[0], int(gt_dp[1]) if gt_dp[1].isdigit() else 0
            
            # 1. 必须命中 4D 位点字典 (完美对齐 Perl)
            if f"{chrom}:{pos}" not in target_positions:
                continue

            # 2. 判定深度和基因型
            if dp >= 5 and gt != ".":
                sites += 1
                if gt == "0/1" or gt == "1/0":
                    het += 1
                    
        if sites > 0:
            het_val = f"{het / sites:.6f}"
            
    if os.path.exists(tmp_bam): os.remove(tmp_bam)
    if os.path.exists(f"{tmp_bam}.bai"): os.remove(f"{tmp_bam}.bai")

    with open(tmp_res, "a") as f:
        f.write(f"{cov}X\t{i}\t{het_val}\n")

def print_progress(completed, total, start_time, node_id, stage=""):
    percent = completed / total
    bar_length = 30
    filled_length = int(bar_length * percent)
    bar = '#' * filled_length + '-' * (bar_length - filled_length)
    elapsed = time.time() - start_time
    eta = (elapsed / completed) * (total - completed) if completed > 0 else 0
    eta_str = time.strftime("%H:%M:%S", time.gmtime(eta))
    elapsed_str = time.strftime("%H:%M:%S", time.gmtime(elapsed))
    
    sys.stdout.write(f"\rNode {node_id} {stage} |{bar}| {completed}/{total} ({percent*100:.1f}%) | Elapsed: {elapsed_str} | ETA: {eta_str}")
    sys.stdout.flush()
    if completed == total:
        sys.stdout.write("\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--node-id", type=int, required=True, help="Current node ID")
    parser.add_argument("--total-nodes", type=int, default=1, help="Total number of nodes")
    args = parser.parse_args()

    os.makedirs(BASE_POOL_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(SHM_TMP_DIR, exist_ok=True)

    print("==================================================")
    print(f"Starting Bootstrap Pipeline on Node {args.node_id}")
    print("==================================================")
    
    if not os.path.exists(MAP_FILE):
        print(f"[FATAL] {MAP_FILE} not found!")
        return

    # 物理行号交替分配
    all_samples = []
    with open(MAP_FILE, "r") as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3: continue
            all_samples.append((parts[0], parts[2]))

    my_samples = []
    for idx, (sample, ref) in enumerate(all_samples):
        if idx % args.total_nodes == args.node_id:
            my_samples.append((sample, ref))

    print(f"Node {args.node_id} is assigned {len(my_samples)} samples.")
    print("Assigned list:")
    for s, r in my_samples:
        print(f"  - {s} ({r})")
    print("--------------------------------------------------")

    extract_individual_30x_pools(my_samples)

    samples_info = []
    for sample, ref in my_samples:
        bam = f"{BASE_POOL_DIR}/{sample}_30X.bam"
        if os.path.exists(bam) and os.path.exists(f"{bam}.bai"):
            samples_info.append((sample, ref))
        else:
            print(f"[WARNING] No base pool found for {sample}. It will not be processed.")

    if not samples_info:
        print("\n[ERROR] No samples available for bootstrap. Exiting.")
        return

    # 主进程一次性把 4D 位点读进内存
    print("\n[STEP 2] Pre-loading 4D site dictionaries into memory...")
    for sample, ref in samples_info:
        if ref in GLOBAL_4D_SITES:
            continue
        four_d_file = f"{BED_DIR}/{ref}.4Dsites.txt"
        if os.path.exists(four_d_file):
            sites_set = set()
            with open(four_d_file, 'r') as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split("\t") 
                    if len(parts) >= 2:
                        sites_set.add(f"{parts[0]}:{parts[1]}")
            GLOBAL_4D_SITES[ref] = sites_set
            print(f"  -> Loaded {len(sites_set)} sites for {ref}", flush=True)
        else:
            print(f"  -> ⚠️ Warning: {four_d_file} not found!", flush=True)

    print("\n[STEP 3] Scanning existing results for resuming...")
    tasks = []
    for sample, ref in samples_info:
        res_file = f"{OUT_DIR}/{sample}_bootstrap_het.txt"
        finished_tasks = set()
        if os.path.exists(res_file):
            with open(res_file, "r") as rf:
                next(rf)
                for line in rf:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3 and parts[2].strip() != "":
                        cov_str, i_str = parts[0], parts[1]
                        if cov_str.endswith('X') and i_str.isdigit():
                            finished_tasks.add((int(cov_str[:-1]), int(i_str)))
        else:
            with open(res_file, "w") as rf:
                rf.write("Coverage\tIteration\tHeterozygosity\n")
                
        for cov in TARGET_COVS:
            tmp_res = f"{OUT_DIR}/tmp_res_{sample}_{cov}.txt"
            if os.path.exists(tmp_res):
                with open(tmp_res, "r") as tf:
                    for line in tf:
                        parts = line.strip().split('\t')
                        if len(parts) >= 3 and parts[2].strip() != "":
                            cov_str, i_str = parts[0], parts[1]
                            if cov_str.endswith('X') and i_str.isdigit():
                                finished_tasks.add((int(cov_str[:-1]), int(i_str)))

        for cov in TARGET_COVS:
            for i in range(1, ITERATIONS + 1):
                if (cov, i) not in finished_tasks:
                    tasks.append((sample, ref, cov, i))
    
    print(f"Remaining micro-tasks to run on Node {args.node_id}: {len(tasks)}")
    
    if len(tasks) == 0:
        print("\nAll tasks for this node are already finished! Exiting.")
        return

    print(f"\n[STEP 4] Running {len(tasks)} micro-tasks with {MAX_WORKERS} workers...")
    start_time = time.time()
    
    with multiprocessing.Pool(processes=MAX_WORKERS) as pool:
        completed = 0
        for _ in pool.imap_unordered(process_task, tasks):
            completed += 1
            if completed % 10 == 0 or completed == len(tasks):
                print_progress(completed, len(tasks), start_time, args.node_id, stage="Calculating")

    print(f"\nAll micro-tasks completed on Node {args.node_id} in {time.time() - start_time:.2f}s.")
    
    print(f"\n[STEP 5] Merging and sorting results...")
    for sample, ref in samples_info:
        res_file = f"{OUT_DIR}/{sample}_bootstrap_het.txt"
        for cov in TARGET_COVS:
            tmp_res = f"{OUT_DIR}/tmp_res_{sample}_{cov}.txt"
            if os.path.exists(tmp_res):
                lines = open(tmp_res).readlines()
                valid_lines = [l for l in lines if len(l.strip().split('\t')) >= 3 and l.strip().split('\t')[2].strip() != ""]
                valid_lines.sort(key=lambda x: int(x.split('\t')[1]))
                with open(res_file, "a") as f:
                    f.writelines(valid_lines)
                os.remove(tmp_res)
                
        if os.path.exists(res_file):
            lines = open(res_file).readlines()
            if len(lines) > 1:
                header, data = lines[0], lines[1:]
                unique_data = {}
                for line in data:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        key = (parts[0], parts[1])
                        unique_data[key] = line
                final_data = list(unique_data.values())
                final_data.sort(key=lambda x: (int(x.split('\t')[0].replace('X', '')), int(x.split('\t')[1])))
                with open(res_file, "w") as f:
                    f.write(header)
                    f.writelines(final_data)
        print(f"Merged {sample}")

    print(f"\nNode {args.node_id} pipeline finished successfully!")

if __name__ == "__main__":
    main()
