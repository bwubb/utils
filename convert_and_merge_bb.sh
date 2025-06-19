#!/bin/bash

# Required references
CHAIN="hg19ToHg38.over.chain"
GENE_BED="ref_genes_ucsc.bed"

# AGILENT S04380110
bigBedToBed S04380110_Covered.bb temp1.bed
cut -f1-3 temp1.bed > temp1.3col.bed
bedtools intersect -a temp1.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > temp1.annotated.bed
bedtools sort -i temp1.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > SureSelect-Exon_v5.S04380110.Covered.hg38.bed

bigBedToBed S04380110_Regions.bb temp2.bed
cut -f1-3 temp2.bed > temp2.3col.bed
bedtools intersect -a temp2.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > temp2.annotated.bed
bedtools sort -i temp2.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > SureSelect-Exon_v5.S04380110.Regions.hg38.bed

# AGILENT S07604514
bigBedToBed S07604514_Covered.bb temp3.bed
cut -f1-3 temp3.bed > temp3.3col.bed
bedtools intersect -a temp3.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > temp3.annotated.bed
bedtools sort -i temp3.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > SureSelect-Exon_v6.S07604514.Covered.hg38.bed

bigBedToBed S07604514_Regions.bb temp4.bed
cut -f1-3 temp4.bed > temp4.3col.bed
bedtools intersect -a temp4.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > temp4.annotated.bed
bedtools sort -i temp4.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > SureSelect-Exon_v6.S07604514.Regions.hg38.bed

# AGILENT S07604715
bigBedToBed S07604715_Covered.bb temp5.bed
cut -f1-3 temp5.bed > temp5.3col.bed
bedtools intersect -a temp5.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > temp5.annotated.bed
bedtools sort -i temp5.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > SureSelect-Exon_v6+COSMIC.S07604715.Covered.hg38.bed

bigBedToBed S07604715_Regions.bb temp6.bed
cut -f1-3 temp6.bed > temp6.3col.bed
bedtools intersect -a temp6.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > temp6.annotated.bed
bedtools sort -i temp6.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > SureSelect-Exon_v6+COSMIC.S07604715.Regions.hg38.bed

# AGILENT S31285117
bigBedToBed S31285117_Covered.bb temp7.bed
cut -f1-3 temp7.bed > temp7.3col.bed
bedtools intersect -a temp7.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > temp7.annotated.bed
bedtools sort -i temp7.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > SureSelect-Exon_v7.S31285117.Covered.hg38.bed

bigBedToBed S31285117_Regions.bb temp8.bed
cut -f1-3 temp8.bed > temp8.3col.bed
bedtools intersect -a temp8.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > temp8.annotated.bed
bedtools sort -i temp8.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > SureSelect-Exon_v7.S31285117.Regions.hg38.bed

# TWIST
bigBedToBed Twist_Exome_Target_hg38.bb temp9.bed
cut -f1-3 temp9.bed > temp9.3col.bed
bedtools intersect -a temp9.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > temp9.annotated.bed
bedtools sort -i temp9.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > Twist-Exome_v2.Regions.hg38.bed

# TWIST
bigBedToBed Twist_Exome_Target_hg38.bb temp9.bed
cut -f1-3 temp9.bed > temp9.3col.bed
bedtools intersect -a temp9.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > Twist_Exome_Target_hg38.annotated.hg38.bed
bedtools sort -i Twist_Exome_Target_hg38.annotated.hg38.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > Twist-Exome_v2.Regions.hg38.bed

# xgen-exome-hyb-panel-probes-hg38.bb
bigBedToBed xgen-exome-hyb-panel-probes-hg38.bb tmp10.bed
cut -f1-3 tmp10.bed > tmp10.3col.bed
bedtools intersect -a tmp10.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > xgen-exome-hyb-panel-probes-hg38.annotated.bed
bedtools sort -i xgen-exome-hyb-panel-probes-hg38.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > xGen-Exome-hyb.Covered.hg38.bed

# xgen-exome-hyb-panel-targets-hg38.bb
bigBedToBed xgen-exome-hyb-panel-targets-hg38.bb tmp11.bed
cut -f1-3 tmp11.bed > tmp11.3col.bed
bedtools intersect -a tmp11.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > xgen-exome-hyb-panel-targets-hg38.annotated.bed
bedtools sort -i xgen-exome-hyb-panel-targets-hg38.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > xGen-Exome-hyb.Regions.hg38.bed

# xgen-exome-hyb-panel-v2-probes-hg38.bb
bigBedToBed xgen-exome-hyb-panel-v2-probes-hg38.bb tmp12.bed
cut -f1-3 tmp12.bed > tmp12.3col.bed
bedtools intersect -a tmp12.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > xgen-exome-hyb-panel-v2-probes-hg38.annotated.bed
bedtools sort -i xgen-exome-hyb-panel-v2-probes-hg38.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > xGen-Exome-hyb_v2.Covered.hg38.bed

# xgen-exome-hyb-panel-v2-targets-hg38.bb
bigBedToBed xgen-exome-hyb-panel-v2-targets-hg38.bb tmp13.bed
cut -f1-3 tmp13.bed > tmp13.3col.bed
bedtools intersect -a tmp13.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > xgen-exome-hyb-panel-v2-targets-hg38.annotated.bed
bedtools sort -i xgen-exome-hyb-panel-v2-targets-hg38.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > xGen-Exome-hyb_v2.Regions.hg38.bed

# SeqCap capture targets (hg19) liftover
bigBedToBed sorted_SeqCap_EZ_Exome_v3_capture_targets_hg19.bb tmp14.bed
liftOver tmp14.bed $CHAIN tmp14.hg38.bed /dev/null
cut -f1-3 tmp14.hg38.bed > tmp14.3col.bed
bedtools intersect -a tmp14.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > sorted_SeqCap_EZ_Exome_v3_capture_targets_hg38.annotated.bed
bedtools sort -i sorted_SeqCap_EZ_Exome_v3_capture_targets_hg38.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > NimbleGen-SeqCap-Exome_v3.Covered.hg38.bed

# SeqCap primary targets (hg19) liftover
bigBedToBed sorted_SeqCap_EZ_Exome_v3_primary_targets_hg19.bb tmp15.bed
liftOver tmp15.bed $CHAIN tmp15.hg38.bed /dev/null
cut -f1-3 tmp15.hg38.bed > tmp15.3col.bed
bedtools intersect -a tmp15.3col.bed -b $GENE_BED -wa -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7,1000,"+"}' > sorted_SeqCap_EZ_Exome_v3_primary_targets_hg38.annotated.bed
bedtools sort -i sorted_SeqCap_EZ_Exome_v3_primary_targets_hg38.annotated.bed | bedtools merge -i - -c 4 -o distinct | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,1000,"+"}' > NimbleGen-SeqCap-Exome_v3.Regions.hg38.bed
