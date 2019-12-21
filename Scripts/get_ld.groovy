bcftools view --regions I:12250000-12500000  /Users/Stefan/UCLA/Genomics_Data/VCFs/328_Soft_csq_annotated.vcf.gz -Oz -o propionate_finemap.vcf.gz

bcftools view --regions I:12250000-12500000 /Users/Stefan/UCLA/Genomics_Data/VCFs/328_Soft_csq_annotated.vcf.gz |\
awk '$0 !~ "#" {print $1":"$2}' > I_12250000-12500000.txt

plink --vcf propionate_finemap.vcf.gz \
	--snps-only \
	--maf 0.05 \
	--allow-extra-chr \
	--set-missing-var-ids @:# \
	--geno \
	--make-bed \
	--recode vcf-iid bgz \
	--extract I_12250000-12500000.txt \
	--out prop_markers

plink --r2 dprime with-freqs \
	--allow-extra-chr \
	--ld-window-r2 0 \
	--ld-snp I:12385811 \
	--chr I \
	--ld-window 28021 \
	--ld-window-kb 6000 \
	--out propionate_chrI_d \
	--set-missing-var-ids @:# \
	--vcf prop_markers.vcf.gz


plink --r2 with-freqs \
	--allow-extra-chr \
	--ld-window-r2 0 \
	--ld-snp I:12385811 \
	--chr I \
	--ld-window 28021 \
	--ld-window-kb 6000 \
	--out propionate_chrI \
	--set-missing-var-ids @:# \
	--vcf prop_markers.vcf.gz