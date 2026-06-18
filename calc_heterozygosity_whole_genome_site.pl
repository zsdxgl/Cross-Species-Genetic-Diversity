#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max sum);
use File::Temp qw(tempfile);
use File::Basename;

# ============================================================
# 参数解析
# ============================================================
my ($vcf_file, $out_file, $total_sites, $threads, $help);
GetOptions(
    'vcf=s'          => \$vcf_file,
    'out=s'          => \$out_file,
    'total-sites=i'  => \$total_sites,
    'threads=i'      => \$threads,
    'help'           => \$help,
) or die "Usage: $0 --vcf <vcf.gz> --out <output.txt> --total-sites <N> --threads <N>\n";

if ($help || !$vcf_file || !$total_sites) {
    die "Usage: $0 --vcf <vcf.gz> --out <output.txt> --total-sites <N> --threads <N>\n";
}

$out_file ||= "heterozygosity.txt";
$threads ||= 1;

# ============================================================
# 0. 检查并生成 VCF 索引
# ============================================================
print STDERR "[0/6] Checking VCF index...\n";

my $index_file = $vcf_file . ".tbi";
if (!-e $index_file) {
    print STDERR "      Index file not found, creating...\n";
    my $cmd_index = "bcftools index -t $vcf_file 2>/dev/null";
    system($cmd_index) == 0 or die "Failed to create index for $vcf_file\n";
    print STDERR "      Index created successfully\n";
} else {
    print STDERR "      Index file exists\n";
}

# ============================================================
# 第一步：获取样本名
# ============================================================
print STDERR "[1/6] Getting sample names from VCF...\n";
my @samples = `bcftools query -l $vcf_file 2>/dev/null`;
chomp @samples;
my $n_samples = scalar @samples;
die "No samples found in $vcf_file" if $n_samples == 0;
print STDERR "      Found $n_samples samples\n";

# ============================================================
# 第二步：读取所有位点坐标并写入临时文件
# ============================================================
print STDERR "[2/6] Reading all site positions from VCF...\n";

my ($fh_all, $all_sites_file) = tempfile(UNLINK => 1, SUFFIX => ".all");
binmode $fh_all;

my $cmd1 = "bcftools query -f '%CHROM\\t%POS\\n' $vcf_file 2>/dev/null";
open(my $fh1, "$cmd1 |") or die "Cannot run bcftools query: $!";

my $line_count = 0;
while (<$fh1>) {
    chomp;
    $line_count++;
    my ($chr, $pos) = split /\t/;
    print $fh_all "$chr\t$pos\n";
    if ($line_count % 10000000 == 0) {
        print STDERR "      Read $line_count sites from VCF\n";
    }
}
close $fh1;
close $fh_all;
print STDERR "      Total sites in VCF: $line_count\n";

# ============================================================
# 第三步：使用 shuf 命令随机采样
# ============================================================
print STDERR "[3/6] Randomly sampling $total_sites sites from VCF using shuf...\n";

my $sample_file = "sampled_sites_$$.txt";
my $cmd_shuf = "shuf -n $total_sites $all_sites_file > $sample_file";
system($cmd_shuf) == 0 or die "shuf failed: $!";
print STDERR "      Sampled $total_sites sites\n";

# 读取采样位点
open(my $fh_sample, '<', $sample_file) or die "Cannot open $sample_file: $!";
my @sampled_sites = <$fh_sample>;
chomp @sampled_sites;
close $fh_sample;
unlink $sample_file;  # 删除临时采样文件

# ============================================================
# 第四步：第一遍扫描，提取采样位点的 DP 值计算中位数
# ============================================================
print STDERR "[4/6] First pass: extracting DP values for sampled sites...\n";

# 为每个样本创建临时文件存储 DP 值
my @temp_files_dp;
my @temp_filenames_dp;
for (my $i = 0; $i < $n_samples; $i++) {
    my ($fh, $filename) = tempfile(UNLINK => 1, SUFFIX => ".dp");
    $temp_files_dp[$i] = $fh;
    $temp_filenames_dp[$i] = $filename;
}

# 创建采样位点列表文件（用于 bcftools query -R）
my ($fh_sample_list, $sample_list_file) = tempfile(UNLINK => 1, SUFFIX => ".samples");
binmode $fh_sample_list;
foreach my $site (@sampled_sites) {
    print $fh_sample_list "$site\n";
}
close $fh_sample_list;

# 使用 bcftools query 配合 -R 提取采样位点的 DP 值
my $cmd2 = "bcftools query -R $sample_list_file -f '[%DP\\t]\\n' $vcf_file 2>/dev/null";
open(my $fh2, "$cmd2 |") or die "Cannot run bcftools query: $!";

my $line_count2 = 0;
while (<$fh2>) {
    chomp;
    $line_count2++;
    my @dps = split /\t/;
    pop @dps if $dps[-1] eq "";
    next if @dps < $n_samples;
    
    for (my $i = 0; $i < $n_samples; $i++) {
        if (defined $dps[$i] && $dps[$i] =~ /^\d+$/) {
            print {$temp_files_dp[$i]} $dps[$i], "\n";
        }
    }
    
    if ($line_count2 % 100000 == 0) {
        print STDERR "      Processed $line_count2 sites for DP collection\n";
    }
}
close $fh2;
print STDERR "      Total sites processed for DP: $line_count2\n";

# 关闭临时文件
for (my $i = 0; $i < $n_samples; $i++) {
    close $temp_files_dp[$i];
}

# ============================================================
# 第五步：计算每个样本的 DP 中位数和阈值
# ============================================================
print STDERR "[5/6] Calculating DP medians and thresholds...\n";

my @dp_low;
my @dp_high;
my @dp_median;

for (my $i = 0; $i < $n_samples; $i++) {
    open(my $fh, '<', $temp_filenames_dp[$i]) or die "Cannot open $temp_filenames_dp[$i]: $!";
    my @dps;
    while (<$fh>) {
        chomp;
        next unless /^\d+$/;
        push @dps, $_;
    }
    close $fh;
    
    my $median = 0;
    if (@dps) {
        my @sorted = sort { $a <=> $b } @dps;
        if (@sorted % 2 == 0) {
            $median = ($sorted[@sorted/2 - 1] + $sorted[@sorted/2]) / 2;
        } else {
            $median = $sorted[(@sorted-1)/2];
        }
        $dp_median[$i] = $median;
        $dp_low[$i] = 5;  # 最低阈值固定为 5
        $dp_high[$i] = $median * 2;
        print STDERR "      Sample $samples[$i]: median=$median low=5 high=$dp_high[$i]\n";
    } else {
        $dp_low[$i] = 5;
        $dp_high[$i] = 1e9;
        print STDERR "      Sample $samples[$i]: no DP data, using default threshold\n";
    }
}

# ============================================================
# 第六步：第二遍扫描，使用阈值计算杂合度
# ============================================================
print STDERR "[6/6] Second pass: calculating heterozygosity with DP filter...\n";

my @sample_totals = (0) x $n_samples;
my @sample_hets = (0) x $n_samples;

my $cmd3 = "bcftools query -f '[%GT\\t%DP\\t]\\n' $vcf_file 2>/dev/null";
open(my $fh3, "$cmd3 |") or die "Cannot run bcftools query: $!";

my $line_count3 = 0;
while (<$fh3>) {
    chomp;
    $line_count3++;
    my @vals = split /\t/;
    pop @vals if $vals[-1] eq "";
    next if @vals < $n_samples * 2;
    
    for (my $i = 0; $i < $n_samples; $i++) {
        my $gt = $vals[$i*2];
        my $dp = $vals[$i*2+1];
        next unless $gt =~ /^[0-9][\/\|][0-9]$/;
        next unless $dp =~ /^\d+$/ && $dp >= 5;
        next unless $dp >= $dp_low[$i] && $dp <= $dp_high[$i];
        
        $sample_totals[$i]++;
        if ($gt =~ /^0[\/\|]1$/ || $gt =~ /^1[\/\|]0$/) {
            $sample_hets[$i]++;
        }
    }
    
    if ($line_count3 % 10000000 == 0) {
        print STDERR "      Processed $line_count3 sites for heterozygosity\n";
    }
}
close $fh3;
print STDERR "      Total sites processed for heterozygosity: $line_count3\n";

# ============================================================
# 第七步：输出结果
# ============================================================
open(my $out, '>', $out_file) or die "Cannot write to $out_file: $!";
print $out "Sample\tTotal_Sites\tHeterozygous_Sites\tHeterozygosity\tDP_Median\n";
for (my $i = 0; $i < $n_samples; $i++) {
    my $total = $sample_totals[$i] || 0;
    my $het = $sample_hets[$i] || 0;
    my $ratio = $total > 0 ? sprintf("%.6f", $het/$total) : "NA";
    my $median = $dp_median[$i] || "NA";
    print $out join("\t", $samples[$i], $total, $het, $ratio, $median), "\n";
}
close $out;
print STDERR "Done. Results written to $out_file\n";
