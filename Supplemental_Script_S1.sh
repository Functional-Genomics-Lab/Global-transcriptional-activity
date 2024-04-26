# Define promoter regions of enhancers as [-400, 100] bp around TSS
bedtools flank -l 400 -r 0 -s -i eRNA_GM.bed -g chromInfo.txt |bedtools slop -l 0 -r 100 -s -i - -g chromInfo.txt >eRNA_GM.promoter.bed
bedtools flank -l 500 -r 0 -s -i eRNA_GM.bed -g chromInfo.txt |bedtools slop -l 0 -r 500 -s -i - -g chromInfo.txt >eRNA_GM.flank500.bed
bedtools flank -l 250 -r 0 -s -i eRNA_GM.bed -g chromInfo.txt |bedtools slop -l 0 -r 250 -s -i - -g chromInfo.txt >eRNA_GM.flank250.bed
bedtools flank -l 250 -r 0 -s -i refseq_groseq.bed -g chromInfo.txt |bedtools slop -l 0 -r 250 -s -i - -g chromInfo.txt >refseq_groseq.flank250.bed
bedtools flank -l 1 -r 0 -s -i eRNA_GM.bed -g chromInfo.txt >eRNA_GM.tss.bed
bedtools flank -l 1 -r 0 -s -i refseq_groseq.bed -g chromInfo.txt >refseq_groseq.tss.bed

# 1. Inducible enhancers
# Define eRNA promoters
grep -wf Inducible_eRNA.list eRNA_GM.promoter.bed >Inducible_eRNA.promoter.bed
grep -wf Inducible_eRNA.list eRNA_GM.flank500.bed >Inducible_eRNA.flank500.bed
# 1.1 Find motifs enriched in inducible enhancers
# Motif occurrence distributed in the promoter region
grep -wf Inducible_eRNA.list eRNA_GM.tss.bed >Inducible_eRNA.tss.bed
grep -wvf Inducible_eRNA.list eRNA_GM.tss.bed >nonInducible_eRNA.tss.bed
liftOver Inducible_eRNA.tss.bed ~/Documents/eRNA/link/UCSC/hg18ToHg19.over.chain Inducible_eRNA.tss.hg19.bed Unmapped
liftOver nonInducible_eRNA.tss.bed ~/Documents/eRNA/link/UCSC/hg18ToHg19.over.chain nonInducible_eRNA.tss.hg19.bed Unmapped
annotatePeaks.pl Inducible_eRNA.tss.bed hg18 -m HOCOMOCOv10_HUMAN_mono_homer_format_0.0001.motif -size 1000 -hist 10 > Motif_HOCO_Inducible_enhancer.profile
# 1.2 Conservation of enhancers
annotatePeaks.pl Inducible_eRNA.tss.bed hg18 -wig Conservation_hg18/phastCons44way/primates/all.phastCons44way.primates.wigFix -size 5000 -hist 10 -cpu 6 > phastCons44way_primate_Inducible_enhancer.profile
annotatePeaks.pl Inducible_eRNA.tss.bed hg18 -wig Conservation_hg18/phastCons44way/mammal/all.phastCons44way.placental.wigFix -size 5000 -hist 10 -cpu 6 > phastCons44way_mammal_Inducible_enhancer.profile
annotatePeaks.pl Inducible_eRNA.tss.bed hg18 -wig Conservation_hg18/phastCons44way/vertebrate/all.phastCons44way.wigFix -size 5000 -hist 10 -cpu 6 > phastCons44way_vertebrate_Inducible_enhancer.profile

# 2. Inducible promoters
# Define gene promoters
grep -wf Inducible_genes.list refseq_groseq.promoter.bed >Inducible_genes.promoter.bed
grep -wf Inducible_genes.list refseq_groseq.flank500.bed >Inducible_genes.flank500.bed
# 2.1 Find motifs enriched in inducible promoters
# Motif occurrence distributed in the promoter region
grep -wf Inducible_genes.list refseq_groseq.tss.bed >Inducible_genes.tss.bed
annotatePeaks.pl Inducible_genes.tss.bed hg18 -m HOCOMOCOv10_HUMAN_mono_homer_format_0.0001.motif -size 1000 -hist 10 > Motif_HOCO_Inducible_promoter.profile
# 2.2 Conservation of enhancers
annotatePeaks.pl Inducible_genes.tss.bed hg18 -wig Conservation_hg18/phastCons44way/primates/all.phastCons44way.primates.wigFix -size 5000 -hist 10 -cpu 6 > phastCons44way_primate_Inducible_promoter.profile
annotatePeaks.pl Inducible_genes.tss.bed hg18 -wig Conservation_hg18/phastCons44way/mammal/all.phastCons44way.placental.wigFix -size 5000 -hist 10 -cpu 6 > phastCons44way_mammal_Inducible_promoter.profile
annotatePeaks.pl Inducible_genes.tss.bed hg18 -wig Conservation_hg18/phastCons44way/vertebrate/all.phastCons44way.wigFix -size 5000 -hist 10 -cpu 6 > phastCons44way_vertebrate_Inducible_promoter.profile

