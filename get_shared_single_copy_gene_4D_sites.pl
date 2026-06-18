#!/usr/bin/perl -w
use strict;
use List::MoreUtils ':all';
use List::Util qw(min max);

# ============================================================================
# 脚本功能：
#   1. 读取每个物种的染色体长度文件 (.chr.length,例如anguilla_Japonica.chr.length)
#   2. 读取每个物种的 BUSCO 完整表格 (.full_table.tsv,例如anguilla_Japonica.full_table.tsv)
#   3. 提取所有 "Complete" 的单拷贝基因，过滤边缘和重叠基因
#   4. 找出所有物种共有的共享单拷贝基因
#   5. 输出每个物种共享基因的区间文件 (.intervals.txt，例如anguilla_Japonica.intervals.txt)，按染色体和起始位置排序
#   6. 输出每个物种共享基因内部的4D位点文件 (.4Dsites.txt,例如anguilla_Japonica.4Dsites.txt)，按染色体和位置排序
# ============================================================================

# ----------------------------------------------------------------------------
# 统一物种列表（24个物种）
# ----------------------------------------------------------------------------
my @species = (
"anguilla_Japonica","bahaba_taipingensis","channa_argus","cheilinus_undulatus","ctenopharyngodon_idella","cynoglossus_semilaevis","dicentrarchus_labrax","hypophthalmichthys_nobilis","larimichthys_crocea","lateolabrax_maculatus","megalobrama_amblycephala","micropterus_salmoides","misgurnus_anguillicaudatus","oreochromis_niloticus","paralichthys_olivaceus","scophthalmus_maximus","siniperca_chuatsi","siniperca_roulei","tachysurus_fulvidraco","takifugu_rubripes","trachinotus_ovatus","salmo_salar","cyprinus_carpio","oncorhynchus_mykiss"
);

# 初始化数据结构
my (%chr, %pos, @genes4spall);
my (%overlap_genes);

# ----------------------------------------------------------------------------
# 排序函数：用于区间文件（按第一列字符顺序，第二列数字顺序）
# ----------------------------------------------------------------------------
sub sort_by_chr_pos {
    my ($a_str, $b_str) = @_;
    return 0 if not defined $a_str or not defined $b_str;
    my ($chr_a, $pos_a) = split /\t/, $a_str;
    my ($chr_b, $pos_b) = split /\t/, $b_str;
    my $chr_cmp = ($chr_a cmp $chr_b);
    return $chr_cmp if $chr_cmp != 0;
    return ($pos_a <=> $pos_b);
}

# ----------------------------------------------------------------------------
# 排序函数：用于4D位点文件（按第一列字符顺序，第二列数字顺序）
# ----------------------------------------------------------------------------
sub sort_4dsite_by_chr_pos {
    my ($a_str, $b_str) = @_;
    return 0 if not defined $a_str or not defined $b_str;
    my ($chr_a, $pos_a) = split /\t/, $a_str;
    my ($chr_b, $pos_b) = split /\t/, $b_str;
    my $chr_cmp = ($chr_a cmp $chr_b);
    return $chr_cmp if $chr_cmp != 0;
    return ($pos_a <=> $pos_b);
}

# ----------------------------------------------------------------------------
# 第一步：读取每个物种的染色体长度文件
# ----------------------------------------------------------------------------
foreach my $sp (@species) {
    my $file = $sp . ".chr.length";
    open(my $fh, '<', $file) or die "无法打开 $file: $!";
    while (my $line = <$fh>) {
        chomp $line;
        my @col = split /\t/, $line;
        $chr{$sp}{$col[0]} = $col[1];
    }
    close $fh;
    print STDERR "完成读取 $file\n";
}

# ----------------------------------------------------------------------------
# 第二步：处理每个物种的 BUSCO 完整表格
# ----------------------------------------------------------------------------
foreach my $sp (@species) {
    my $file = $sp . ".full_table.tsv";
    open(my $fh, '<', $file) or die "无法打开 $file: $!";
    
    my @all_intervals;
    my %overlap_in_sp = ();  # 每个物种独立的重叠标记
    
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^#/;
        if ($line =~ /Complete/) {
            my @col = split /\t/, $line;
            
            # 检查染色体是否存在
            unless (exists $chr{$sp}{$col[2]}) {
                warn "警告: $sp 的染色体 $col[2] 未定义，跳过基因 $col[0]\n";
                next;
            }
            
            if ($col[3] >= 1000 && $col[3] <= ($chr{$sp}{$col[2]} - 1000) &&
                $col[4] >= 1000 && $col[4] <= ($chr{$sp}{$col[2]} - 1000)) {
                
                my ($start, $end) = ($col[3] < $col[4]) ? ($col[3], $col[4]) : ($col[4], $col[3]);
                
                my $interval = {
                    chr   => $col[2],
                    start => $start,
                    end   => $end,
                    gene  => $col[0]
                };
                
                push @all_intervals, $interval;
                push @genes4spall, $col[0];
            }
        }
    }
    close $fh;
    
    # ------------------------------------------------------------------------
    # 第三步：检测并标记重叠的基因（同一染色体上位置重叠）
    # ------------------------------------------------------------------------
    # 按染色体分组
    my %intervals_by_chr;
    foreach my $interval (@all_intervals) {
        push @{$intervals_by_chr{$interval->{chr}}}, $interval;
    }
    
    # 对每条染色体检测重叠
    foreach my $chr_name (keys %intervals_by_chr) {
        my @sorted = sort { $a->{start} <=> $b->{start} } @{$intervals_by_chr{$chr_name}};
        
        for (my $k = 0; $k < @sorted; $k++) {
            my $cur = $sorted[$k];
            for (my $m = $k + 1; $m < @sorted; $m++) {
                my $next = $sorted[$m];
                last if $next->{start} > $cur->{end};
                if ($next->{start} < $cur->{end}) {
                    $overlap_in_sp{$cur->{gene}} = 1;
                    $overlap_in_sp{$next->{gene}} = 1;
                }
            }
        }
    }
    
    # ------------------------------------------------------------------------
    # 第四步：保存非重叠基因的位置信息
    # ------------------------------------------------------------------------
    foreach my $interval (@all_intervals) {
        next if exists $overlap_in_sp{$interval->{gene}};
        $pos{$sp}{$interval->{gene}} = join("\t", $interval->{chr}, $interval->{start}, $interval->{end});
    }
    
    my $total_genes = scalar @all_intervals;
    my $overlap_count = scalar keys %overlap_in_sp;
    my $kept_count = $total_genes - $overlap_count;
    print STDERR "物种 $sp: 总基因数 = $total_genes, 重叠基因数 = $overlap_count, 保留基因数 = $kept_count\n";
}

# ----------------------------------------------------------------------------
# 第五步：识别所有物种共有的共享单拷贝基因
# ----------------------------------------------------------------------------
my @unique_genes = uniq @genes4spall;
my @shared_genes;

foreach my $gene (@unique_genes) {
    my $found_in_all = 1;
    foreach my $sp (@species) {
        unless (exists $pos{$sp}{$gene}) {
            $found_in_all = 0;
            last;
        }
    }
    if ($found_in_all) {
        push @shared_genes, $gene;
    }
}

print STDERR "\n共享基因数量: ", scalar @shared_genes, "\n";

if (scalar @shared_genes == 0) {
    print STDERR "错误：未找到任何共享基因，请检查输入文件。\n";
    exit 1;
}

# ----------------------------------------------------------------------------
# 第六步：输出每个物种的共享基因区间和4D位点（按染色体和位置排序）
# ----------------------------------------------------------------------------
foreach my $sp (@species) {
    # ========== 6.1 处理基因区间文件 ==========
    my @interval_lines;
    foreach my $gene (@shared_genes) {
        if (exists $pos{$sp}{$gene}) {
            push @interval_lines, $pos{$sp}{$gene};
        }
    }
    
    # 排序
    my @sorted_intervals = sort sort_by_chr_pos @interval_lines;
    
    # 输出区间文件
    my $interval_out = $sp . ".intervals.txt";
    open(my $out_fh, '>', $interval_out) or die "无法创建 $interval_out: $!";
    foreach my $line (@sorted_intervals) {
        print $out_fh $line, "\n";
    }
    close $out_fh;
    print STDERR "完成写入 $interval_out（已排序，共 ", scalar @sorted_intervals, " 行）\n";
    
    # ========== 6.2 处理4D位点文件 ==========
    my @site_lines;
    foreach my $gene (@shared_genes) {
        my $site_file = $sp . ".reseqtools_out/" . $gene . ".4Dsite.gz";
        
        if (-e $site_file) {
            open(my $in_fh, "gunzip -cd $site_file |") or warn "无法打开 $site_file: $!";
            while (my $line = <$in_fh>) {
                chomp $line;
                next unless $line =~ /\S/;
                my @col = split /\t/, $line;
                push @site_lines, $col[0] . "\t" . $col[1];
            }
            close $in_fh;
        } else {
            warn "警告: 文件不存在 $site_file\n";
        }
    }
    
    # 排序
    my @sorted_sites = sort sort_4dsite_by_chr_pos @site_lines;
    
    # 输出4D位点文件
    my $site_out = $sp . ".4Dsites.txt";
    open(my $site_fh, '>', $site_out) or die "无法创建 $site_out: $!";
    foreach my $line (@sorted_sites) {
        print $site_fh $line, "\n";
    }
    close $site_fh;
    print STDERR "完成写入 $site_out（已排序，共 ", scalar @sorted_sites, " 个位点）\n";
}

print STDERR "\n所有处理完成！\n";
