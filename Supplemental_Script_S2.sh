
##########################################################################
# 1. IRF binding sites in inducible enhancers
##########################################################################

computeMatrix reference-point \
--maxThreshold 500 \
-R Data/Inducible_eRNA.IRF.bed \
Data/Inducible_eRNA.tss.bed \
-S Data/bigwig/groseq_hg18/GM6h.bw \
Data/bigwig/UCSC/H3K4me1.bw \
Data/bigwig/UCSC/H3K27ac.bw \
--referencePoint TSS \
-b 1000 -a 2000 -p 6 --missingDataAsZero -q \
--outFileName Heatmap/Inducible_eRNA.IRF.gz

plotProfile --plotHeight 5 --plotWidth 10 --refPointLabel "IRF or TSS" --regionsLabel "IRF" "TSS" --samplesLabel "GROseq 6h" "H3K4me1" "H3K27ac" --numPlotsPerRow 1 -m Heatmap/Inducible_eRNA.IRF.gz -out Heatmap/Inducible_eRNA.IRF.profile.pdf


##########################################################################
# 2. H3K4me3
##########################################################################

computeMatrix reference-point \
--maxThreshold 500 \
-R Data/ISG.tss.hg19.bed \
-S Data/bigwig/chipseq_hg19/H3K4me3_GM0h_treat_pileup.sorted.bdg.bw \
Data/bigwig/chipseq_hg19/H3K4me3_GM6h_treat_pileup.sorted.bdg.bw \
Data/bigwig/chipseq_hg19/H3K4me3_GM12h_treat_pileup.sorted.bdg.bw \
Data/bigwig/chipseq_hg19/H3K4me3_GM18h_treat_pileup.sorted.bdg.bw \
Data/bigwig/chipseq_hg19/H3K4me3_GM24h_treat_pileup.sorted.bdg.bw \
--referencePoint TSS \
-b 2000 -a 2000 -p 6 --missingDataAsZero -q \
--outFileName Heatmap/Inducible_gene.H3K4me3.gz

plotProfile --plotHeight 5 --plotWidth 10 --refPointLabel "TSS" --regionsLabel "Inducible Genes" --samplesLabel "0h" "6h" "12h" "18h" "24h" --numPlotsPerRow 1 -m Heatmap/Inducible_gene.H3K4me3.gz -out Heatmap/Inducible_gene.H3K4me3.profile.pdf

computeMatrix reference-point \
--maxThreshold 500 \
-R Data/H3K4me3_6h-24h_peaks.narrowPeak \
Data/H3K4me3_12h-24h_peaks.narrowPeak \
Data/H3K4me3_18h-24h_peaks.narrowPeak \
-S Data/bigwig/chipseq_hg19/H3K4me3_GM0h_treat_pileup.sorted.bdg.bw \
Data/bigwig/chipseq_hg19/H3K4me3_GM6h_treat_pileup.sorted.bdg.bw \
Data/bigwig/chipseq_hg19/H3K4me3_GM12h_treat_pileup.sorted.bdg.bw \
Data/bigwig/chipseq_hg19/H3K4me3_GM18h_treat_pileup.sorted.bdg.bw \
Data/bigwig/chipseq_hg19/H3K4me3_GM24h_treat_pileup.sorted.bdg.bw \
--referencePoint TSS \
-b 2000 -a 2000 -p 6 --missingDataAsZero -q \
--outFileName Heatmap/Inducible_H3K4me3_peaks.gz

plotHeatmap --sortRegions no --colorList white,forestgreen --refPointLabel "Center" --regionsLabel "6-24h" "12-24h" "18-24h" --samplesLabel "0h" "6h" "12h" "18h" "24h" -m Heatmap/Inducible_H3K4me3_peaks.gz -out Heatmap/Inducible_H3K4me3_peaks.pdf

##########################################################################
# 3. H3K4me1
##########################################################################
computeMatrix reference-point \
--maxThreshold 500 \
-R Data/eRNA_GM.tss.hg19.bed \
-S /home/pengxie/Documents/eRNA/link/UCSC/wgEncodeBroadHistoneGm12878H3k04me1StdSigV2.bigWig \
--referencePoint TSS \
-b 2000 -a 2000 -p 6 --missingDataAsZero -q \
--outFileName Heatmap/eRNA_TSS.H3K4me1.gz

plotHeatmap --refPointLabel "TSS" --regionsLabel "" --samplesLabel "H3K4me1" --colorList white,blue -m Heatmap/eRNA_TSS.H3K4me1.gz -out Heatmap/eRNA_TSS.H3K4me1.heatmap.pdf

##########################################################################
# 4. H3K27ac
##########################################################################
computeMatrix reference-point \
--maxThreshold 500 \
-R Data/eRNA_GM.tss.hg19.bed \
-S /home/pengxie/Documents/eRNA/link/UCSC/wgEncodeBroadHistoneGm12878H3k27acStdSig.bigWig \
--referencePoint TSS \
-b 2000 -a 2000 -p 6 --missingDataAsZero -q \
--outFileName Heatmap/eRNA_TSS.H3K27ac.gz

plotHeatmap --refPointLabel "TSS" --regionsLabel "" --samplesLabel "H3K27ac" --colorList white,blue -m Heatmap/eRNA_TSS.H3K27ac.gz -out Heatmap/eRNA_TSS.H3K27ac.heatmap.pdf

##########################################################################
# 5. P300
##########################################################################
computeMatrix reference-point \
--maxThreshold 500 \
-R Data/eRNA_GM.tss.hg19.bed \
-S /home/pengxie/Documents/eRNA/link/UCSC/wgEncodeSydhTfbsGm12878P300bStdSig.bigWig \
--referencePoint TSS \
-b 2000 -a 2000 -p 6 --missingDataAsZero -q \
--outFileName Heatmap/eRNA_TSS.P300.gz

plotHeatmap --refPointLabel "TSS" --regionsLabel "" --samplesLabel "P300" --colorList white,blue -m Heatmap/eRNA_TSS.P300.gz -out Heatmap/eRNA_TSS.P300.heatmap.pdf

##########################################################################
# 6. DNase & Nucleosome
##########################################################################
computeMatrix reference-point \
--maxThreshold 500 \
-R Data/eRNA_GM.tss.hg19.bed \
-S /home/pengxie/Documents/eRNA/link/UCSC/wgEncodeOpenChromDnaseGm12878BaseOverlapSignal.bigWig \
--referencePoint TSS \
-b 2000 -a 2000 -p 6 --missingDataAsZero -q \
--outFileName Heatmap/eRNA_TSS.DNase.gz

plotHeatmap --refPointLabel "TSS" --regionsLabel "" --samplesLabel "DNase" --colorList white,blue -m Heatmap/eRNA_TSS.DNase.gz -out Heatmap/eRNA_TSS.DNase.heatmap.pdf

computeMatrix reference-point \
--maxThreshold 500 \
-R Data/eRNA_GM.tss.hg19.bed \
-S /home/pengxie/Documents/eRNA/link/UCSC/wgEncodeSydhNsomeGm12878Sig.bigWig \
--referencePoint TSS \
-b 2000 -a 2000 -p 6 --missingDataAsZero -q \
--outFileName Heatmap/eRNA_TSS.Nucleosome.gz

plotHeatmap --refPointLabel "TSS" --regionsLabel "" --samplesLabel "Nucleosome" --colorList white,blue -m Heatmap/eRNA_TSS.Nucleosome.gz -out Heatmap/eRNA_TSS.Nucleosome.heatmap.pdf


##########################################################################
# 7. GROseq
##########################################################################
computeMatrix reference-point \
--maxThreshold 500 \
-R Data/eRNA_GM.tss.plus.bed \
-S Data/bigwig/groseq_hg18/GM6h_plus.bw \
--referencePoint TSS \
-b 1000 -a 2000 -p 6 --missingDataAsZero -q \
--outFileName Heatmap/eRNA_TSS.plus.GROseq.gz

plotHeatmap --refPointLabel "TSS" --regionsLabel "" --samplesLabel "GROseq" --colorList white,orange -m Heatmap/eRNA_TSS.plus.GROseq.gz -out Heatmap/eRNA_TSS.plus.GROseq.heatmap.pdf

computeMatrix reference-point \
--maxThreshold 500 \
-R Data/eRNA_GM.tss.minus.bed \
-S Data/bigwig/groseq_hg18/GM6h_minus.bw \
--referencePoint TSS \
-b 1000 -a 2000 -p 6 --missingDataAsZero -q \
--outFileName Heatmap/eRNA_TSS.minus.GROseq.gz

plotHeatmap --refPointLabel "TSS" --regionsLabel "" --samplesLabel "GROseq" --colorList white,purple -m Heatmap/eRNA_TSS.minus.GROseq.gz -out Heatmap/eRNA_TSS.minus.GROseq.heatmap.pdf


