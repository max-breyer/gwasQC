SAMPLES = ["HapMap_3_r3_1"]

# change the input of this rule to decide what snakemake intends to generate
rule target:
	input: 
		expand("{sample}_P8.bim", sample=SAMPLES),
        expand("{sample}_dropSummary.txt", sample=SAMPLES),
        #expand("{sample}_dropSummary.txt", sample=SAMPLES)

# exercise 1
rule removeChromosomeZero:
	input:
		expand("{sample}.bim", sample=SAMPLES)
	output:
		"{sample}_QCdrop.SNPs"
	shell:
		"awk '{{if($1==0) print}}' {wildcards.sample}.bim > {wildcards.sample}_QCdrop.SNPs"

rule minorAlleleFreqCheck:
	input:
		expand("{sample}.bed", sample=SAMPLES),
		expand("{sample}.bim", sample=SAMPLES),
		expand("{sample}.fam", sample=SAMPLES),
		expand("{sample}_QCdrop.SNPs", sample=SAMPLES)
	output:
		"{sample}_P2.bed",
		"{sample}_P2.bim",
		"{sample}_P2.fam",
		"{sample}_P2.log"
	shell:
		"plink --bfile {wildcards.sample} --exclude {wildcards.sample}_QCdrop.SNPs --maf 0.01 --make-bed --out {wildcards.sample}_P2"

# exercise 2
rule dropQCsubjects:
	input:
		expand("{sample}_P2.bed", sample=SAMPLES),
		expand("{sample}_P2.bim", sample=SAMPLES),
		expand("{sample}_P2.fam", sample=SAMPLES)
	output:
		"{sample}_P3.bed",			 
		"{sample}_P3.bim",
		"{sample}_P3.fam",
		"{sample}_P3.log",
		"{sample}_QC_subjects.bed",	 
		"{sample}_QC_subjects.bim",
		"{sample}_QC_subjects.fam",
        "{sample}_QC_subjects.log"
		"{sample}_P2_QC.subjects"
	shell:
		"""
		grep NA1 {wildcards.sample}_P2.fam > {wildcards.sample}_P2_QC.subjects
		plink --bfile {wildcards.sample}_P2 --remove {wildcards.sample}_P2_QC.subjects --make-bed --out {wildcards.sample}_P3
		plink --bfile {wildcards.sample}_P2 --keep {wildcards.sample}_P2_QC.subjects --make-bed --out {wildcards.sample}_QC_subjects
		"""

# exerise 3
rule collectNonAutosomalSNP:
	input:
		expand("{sample}_P2.bim", sample=SAMPLES)
	output:
		"{sample}_nonAUT.SNPs"
	shell:
		"awk '{{if($1>22) print}}' {wildcards.sample}_P2.bim > {wildcards.sample}_nonAUT.SNPs"

rule excludeNonAutosomalSNP:
	input:
		expand("{sample}_P3.bed", sample=SAMPLES),
		expand("{sample}_P3.bim", sample=SAMPLES),
		expand("{sample}_P3.fam", sample=SAMPLES),
		expand("{sample}_nonAUT.SNPs", sample=SAMPLES)
	output:
		"{sample}_P3.1.log",
		"{sample}_P4.bed",
		"{sample}_P4.bim",
		"{sample}_P4.fam",
		"{sample}_P4.log"
	run:
		shell("plink --bfile {wildcards.sample}_P2 --mind 0.05 --exclude {wildcards.sample}_nonAUT.SNPs --make-bed --out {wildcards.sample}_P3.1")
		shell("plink --bfile {wildcards.sample}_P3 --keep {wildcards.sample}_P3.1.fam --make-bed --out {wildcards.sample}_P4")
#	shell:
#		"""
#		plink --bfile {wildcards.sample}_P2 --mind 0.05 --exclude {wildcards.sample}_nonAUT.SNPs --make-bed --out {wildcards.sample}_P3.1
#		plink --bfile {wildcards.sample}_P3 --keep {wildcards.sample}_P3.1.fam --make-bed --out {wildcards.sample}_P4
#		rm -rf {wildcards.sample}.P3.1*
#		"""

rule excludeMinorAlleleFrequency:
	input:
		expand("{sample}_P4.bed", sample=SAMPLES),
		expand("{sample}_P4.bim", sample=SAMPLES),
		expand("{sample}_P4.fam", sample=SAMPLES)
	output:
		"{sample}_P5.log",
		"{sample}_P5.bed",
		"{sample}_P5.bim",
		"{sample}_P5.fam"
	shell:
		"plink --bfile {wildcards.sample}_P4 --geno 0.05 --make-bed --out {wildcards.sample}_P5"

# exercise 4
rule checkSex:
	input:
		expand("{sample}_P5.bed", sample=SAMPLES),
		expand("{sample}_P5.bim", sample=SAMPLES),
		expand("{sample}_P5.fam", sample=SAMPLES)
	output:
		"{sample}_P5.sexcheck",
		"{sample}_P5.sexfail"
	shell:
		"""
		plink --bfile {wildcards.sample}_P5 --maf 0.2 --check-sex --out {wildcards.sample}_P5 
		awk '{{if($5=="PROBLEM") print $1, $2}}' {wildcards.sample}_P5.sexcheck > {wildcards.sample}_P5.sexfail
		"""

rule excludeSex:
	input:
		expand("{sample}_P5.bed", sample=SAMPLES),
		expand("{sample}_P5.bim", sample=SAMPLES),
		expand("{sample}_P5.fam", sample=SAMPLES),
		expand("{sample}_P5.sexcheck", sample=SAMPLES),
		expand("{sample}_P5.sexfail", sample=SAMPLES)
	output:
		"{sample}_P6.bed",
		"{sample}_P6.bim",
		"{sample}_P6.fam",
		"{sample}_P6.log"
	shell:
		"plink --bfile {wildcards.sample}_P5 --remove {wildcards.sample}_P5.sexfail --make-bed --out {wildcards.sample}_P6"
# exercise 5
rule checkRelatedness:
	input:
		expand("{sample}_P6.bed", sample=SAMPLES),
		expand("{sample}_P6.bim", sample=SAMPLES),
		expand("{sample}_P6.fam", sample=SAMPLES)
	output:
		"{sample}_P6.prune.in",
		"{sample}_P6.prune.out",
	shell:
		"plink --bfile {wildcards.sample}_P6 --maf 0.1 --exclude {wildcards.sample}_nonAUT.SNPs --indep-pairwise 80 8 0.15 --out {wildcards.sample}_P6"

rule excludeRelatedness:
	input:
		expand("{sample}_P6.bed", sample=SAMPLES),
		expand("{sample}_P6.bim", sample=SAMPLES),
		expand("{sample}_P6.fam", sample=SAMPLES),
		expand("{sample}_P6.prune.in", sample=SAMPLES)
	output:
		"{sample}_P6.IBD.bed",
		"{sample}_P6.IBD.bim",
		"{sample}_P6.IBD.fam",
		"{sample}_P6.IBD.log"
	shell:
		"plink --bfile {wildcards.sample}_P6 --extract {wildcards.sample}_P6.prune.in --make-bed --out {wildcards.sample}_P6.IBD"

rule checkHeterozygosity:
	input:
		expand("{sample}_P6.IBD.bed", sample=SAMPLES),
		expand("{sample}_P6.IBD.bim", sample=SAMPLES),
		expand("{sample}_P6.IBD.fam", sample=SAMPLES)
	output:
		"{sample}_het_check.het",
		"{sample}_het_check_R.het",
		"{sample}_heterozygosity.pdf",
		"{sample}_excess_het_outliers.txt",
		"{sample}_excess_het_outliers_ids_only.txt"
	run:
		shell("plink --bfile {wildcards.sample}_P6.IBD --het --out {wildcards.sample}_het_check")
		shell("""awk '{{if(NR ==1) print $0, "Het_rate"; else if (NR > 1) print $0, ($5 - $3)/$5}}' {wildcards.sample}_het_check.het > {wildcards.sample}_het_check_R.het""")
		shell("Rscript bin/gwasQC_hetcheck_list_plot.R {wildcards.sample}_het_check_R.het {wildcards.sample}_")

rule excludeHeterozygosity:
	input:
		expand("{sample}_P6.IBD.bed", sample=SAMPLES),
		expand("{sample}_P6.IBD.bim", sample=SAMPLES),
		expand("{sample}_P6.IBD.fam", sample=SAMPLES),
		expand("{sample}_excess_het_outliers_ids_only.txt", sample=SAMPLES)
	output:
		"{sample}_P6.IBD.het.bed",
		"{sample}_P6.IBD.het.bim",
		"{sample}_P6.IBD.het.fam",
		"{sample}_P6.IBD.het.log"
		#"{sample}_P6_relateds.IBD.genome.gz"
	run:
		shell("plink --bfile {wildcards.sample}_P6.IBD --remove {wildcards.sample}_excess_het_outliers_ids_only.txt --make-bed --out {wildcards.sample}_P6.IBD.het")
		#shell("plink --bfile {wildcards.sample}_P6.IBD.het --Z-genome --out {wildcards.sample}_P6_relateds.IBD")

rule generateIBD:
	input:
		expand("{sample}_P6.IBD.het.bed", sample=SAMPLES),
		expand("{sample}_P6.IBD.het.bim", sample=SAMPLES),
		expand("{sample}_P6.IBD.het.fam", sample=SAMPLES)
	output:
		"{sample}_P6_relateds.IBD.genome.gz",
		"{sample}_IBD.drops"
	run:
		shell("plink --bfile {wildcards.sample}_P6.IBD.het --Z-genome --out {wildcards.sample}_P6_relateds.IBD")
		shell("zcat {wildcards.sample}_P6_relateds.IBD.genome.gz | awk '{{if($10>0.2) print }}' > {wildcards.sample}_IBD.drops")

rule excludeIBD:
	input:
		expand("{sample}_IBD.drops", sample=SAMPLES),
		expand("{sample}_P6.IBD.het.bed", sample=SAMPLES),
		expand("{sample}_P6.IBD.het.bim", sample=SAMPLES),
		expand("{sample}_P6.IBD.het.fam", sample=SAMPLES)
	output:
		"{sample}_MZ_twins",
		"{sample}_other_relatedness",
		"{sample}_MZ_twins_list1",
		"{sample}_MZ_twins_list2",
		"{sample}_MZ_twins_list",
		"{sample}_unique_MZ_twins_drop_list",
		"{sample}_other_drop_list",
		"{sample}_all_drop_list",
		"{sample}_all_drop_list_final",
		"{sample}_P7.bed",
		"{sample}_P7.bim",
		"{sample}_P7.fam",
		"{sample}_P7.log"
	run:
		shell("awk '{{if($10 >0.95) print }}' {wildcards.sample}_IBD.drops > {wildcards.sample}_MZ_twins")
		shell("awk '{{if($10 <0.95) print }}' {wildcards.sample}_IBD.drops > {wildcards.sample}_other_relatedness")
		shell("awk '{{if(NR != 1) print $1, $2}}' {wildcards.sample}_MZ_twins > {wildcards.sample}_MZ_twins_list1")
		shell("awk '{{if(NR != 1) print $3, $4}}' {wildcards.sample}_MZ_twins > {wildcards.sample}_MZ_twins_list2")
		shell("cat {wildcards.sample}_MZ_twins_list1 {wildcards.sample}_MZ_twins_list2 > {wildcards.sample}_MZ_twins_list")
		shell("sort {wildcards.sample}_MZ_twins_list | uniq > {wildcards.sample}_unique_MZ_twins_drop_list")
		shell("awk '{{print $1, $2}}' {wildcards.sample}_other_relatedness > {wildcards.sample}_other_drop_list")
		shell("cat {wildcards.sample}_other_drop_list {wildcards.sample}_unique_MZ_twins_drop_list > {wildcards.sample}_all_drop_list")
		shell("cat {wildcards.sample}_excess_het_outliers_ids_only.txt {wildcards.sample}_all_drop_list > {wildcards.sample}_all_drop_list_final")
		shell("plink --bfile {wildcards.sample}_P6.IBD.het --remove {wildcards.sample}_all_drop_list_final --make-bed --out {wildcards.sample}_P7")

# exercise 6
rule excludeHardyWeinbergEquilibrium:
	input:
		expand("{sample}_P7.bed", sample=SAMPLES),
		expand("{sample}_P7.bim", sample=SAMPLES),
		expand("{sample}_P7.fam", sample=SAMPLES)
	output:
		"{sample}_P8.bed",
		"{sample}_P8.bim",
		"{sample}_P8.fam",
		"{sample}_P8.log"
	shell:
		"plink --bfile {wildcards.sample}_P7 --hwe 0.000001 --make-bed --out {wildcards.sample}_P8"

# exercise 7
#rule excludeGenomicInflationFactors:
#	input:
#		expand("{sample}_P8.bed", sample=SAMPLES),
#		expand("{sample}_P8.bim", sample=SAMPLES),
#		expand("{sample}_P8.fam", sample=SAMPLES)
#	output:
#
#	shell:
#		"plink --bfile {wildcards.sample}_P3 --logistic --adjust --out {wildcards.sample}_P3"
#		"plink --bfile {wildcards.sample}_P8 --logistic --adjust --out {wildcards.sample}_P8"

rule summaryDropped:
	input:
		expand("{sample}_P2.log", sample=SAMPLES),
		expand("{sample}_P3.log", sample=SAMPLES),
        expand("{sample}_P3.1.log", sample=SAMPLES),
		expand("{sample}_P4.log", sample=SAMPLES),
		expand("{sample}_P5.log", sample=SAMPLES),
		expand("{sample}_P6.log", sample=SAMPLES),
		expand("{sample}_P6.IBD.log", sample=SAMPLES),
		expand("{sample}_P6.IBD.het.log", sample=SAMPLES),
		#expand("{sample}_P6_relateds.IBD.log", sample=SAMPLES),
		#expand("{sample}_het_check.log", sample=SAMPLES),
	output:
		"{sample}_dropSummary.txt"
	shell:
		"Rscript bin/dropSummary.R {wildcards.sample}_dropSummary.txt"












