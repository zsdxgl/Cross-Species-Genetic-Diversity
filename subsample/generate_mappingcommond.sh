#!/bin/bash

DATA_PATH="/data1/guoliang/ZhuoLJ/1.genomescope"
REF_PATH="/data1/guoliang/ZhuoLJ/ref"
WORK_PATH="/data1/guoliang/ZhuoLJ/1.depth_threshold_depth_bootstrap"

mkdir -p "$WORK_PATH"
cd "$WORK_PATH" || exit 1

# 写入样本与参考基因组的映射表
cat << EOF > sample_ref_map.txt
SRR12443991	Micropterus_nigricans	micropterus_salmoides
SRR11300853	Siniperca_chuatsi	siniperca_chuatsi
SRR27213130	Paralichthys_olivaceus	paralichthys_olivaceus
SRR28605936	Megalobrama_amblycephala	megalobrama_amblycephala
SRR27396995	Ctenopharyngodon_idella	ctenopharyngodon_idella
SRR8197016	Larimichthys_crocea	larimichthys_crocea
SRR11842918	Oreochromis_niloticus	oreochromis_niloticus
SRR19974267	Anguilla_japonica	anguilla_Japonica
SRR19215814	Channa_argus	channa_argus
SRR36321447	Takifugu_rubripes	takifugu_rubripes
SRR28761201	Hypophthalmichthys_nobilis	hypophthalmichthys_nobilis
DRR932747	Misgurnus_anguillicaudatus	misgurnus_anguillicaudatus
CRR1985732	Trachinotus_anak	trachinotus_ovatus
SRR9039887	Scophthalmus_maximus	scophthalmus_maximus
SRR31802955	Cynoglossus_semilaevis	cynoglossus_semilaevis
SRR36613568	Oncorhynchus_mykiss	oncorhynchus_mykiss
SRR6080033	Lateolabrax_maculatus	lateolabrax_maculatus
SRR7817079	Tachysurus_fulvidraco	tachysurus_fulvidraco
SRR6019467	Salmo_salar	salmo_salar
SRR8380261	Cyprinus_carpio	cyprinus_carpio
ERR15903826	Dicentrarchus_labrax	dicentrarchus_labrax
SRR31143828	Siniperca_roulei	siniperca_roulei
SRR11496422	Cheilinus_undulatus	cheilinus_undulatus
EOF

# 清空旧的命令文件
> work.pro1.sh
> work.pro2.sh
> work.pro3.sh

SKIP_COUNT=0
GEN_COUNT=0

while IFS=$'\t' read -r sample latin ref_prefix; do

    # 最终产物是 .chr.bam，已存在则跳过
    if [ -s "${WORK_PATH}/${sample}.chr.bam" ]; then
        echo "⏭️  跳过已完成的样本: $sample (已存在 .chr.bam)"
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi

    REF_FILE="${ref_prefix}.fna"
    CHR_FILE="${REF_PATH}/${ref_prefix}.chr.txt"

    if [ ! -f "${REF_PATH}/${REF_FILE}" ]; then
        echo "❌ 严重错误: 找不到参考基因组 ${REF_PATH}/${REF_FILE}，请检查!"
        exit 1
    fi
    if [ ! -f "${CHR_FILE}" ]; then
        echo "❌ 严重错误: 找不到染色体列表 ${CHR_FILE}，请检查!"
        exit 1
    fi

    # =====================================================================
    # 工作流 1：全参考比对 + 坐标排序
    # =====================================================================
    echo "set -o pipefail; bwa mem -t 32 -M -Y -R \"@RG\\tID:${sample}\\tSM:${sample}\" ${REF_PATH}/${REF_FILE} ${DATA_PATH}/${sample}_1.fastq.gz ${DATA_PATH}/${sample}_2.fastq.gz 2> ${WORK_PATH}/${sample}.bwa.log | samtools sort -@ 32 -m 1G -T ${WORK_PATH}/${sample}.sorttmp -o ${WORK_PATH}/${sample}.sorted.bam - 2> ${WORK_PATH}/${sample}.sort.log" >> work.pro1.sh

    # =====================================================================
    # 工作流 2：全参考去重，保留 D.bam
    # =====================================================================
    echo "set -o pipefail;ulimit -n 65000 ; sambamba markdup -t 32 --tmpdir=${WORK_PATH} --sort-buffer-size=2000000 -r ${WORK_PATH}/${sample}.sorted.bam ${WORK_PATH}/${sample}.D.bam 2>${WORK_PATH}/${sample}.dup.log" >> work.pro2.sh
    #echo "rm -f ${WORK_PATH}/${sample}.sorted.bam" >> work.pro2.sh

    # =====================================================================
    # 工作流 3：去重后按染色体过滤
    #   -q 30    : MAPQ >= 30
    #   -f 0x2   : 只保留 proper pair
    #   -F 0x900 : 排除 secondary (0x100) 和 supplementary (0x800)
    # 过滤后重新排序并建索引
    # =====================================================================
    echo "set -o pipefail; CHR_LIST=\$(cut -f1 ${CHR_FILE} | paste -sd' ' -); samtools view -u -q 30 -f 0x2 -F 0x900 ${WORK_PATH}/${sample}.D.bam \${CHR_LIST} | samtools sort -@ 16 -m 2G -T ${WORK_PATH}/${sample}.chrsorttmp -o ${WORK_PATH}/${sample}.chr.bam - && samtools index ${WORK_PATH}/${sample}.chr.bam" >> work.pro3.sh

    GEN_COUNT=$((GEN_COUNT + 1))
done < sample_ref_map.txt

echo "--------------------------------------------------"
echo "✅ 校验通过！"
echo "   - 跳过了 $SKIP_COUNT 个已完成的样本。"
echo "   - 为 $GEN_COUNT 个新样本生成了命令。"
echo "   - work.pro1.sh / work.pro2.sh / work.pro3.sh 已更新。"
