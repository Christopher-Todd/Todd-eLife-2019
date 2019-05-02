# Todd-eLife-2019

Repository for scripts used in the bioinformatic analyses of the following publication:

Functional evaluation of transposable elements as enhancers in mouse embryonic and trophoblast stem cells

Christopher D Todd, Özgen Deniz, Darren Taylor, Miguel R Branco

eLife 2019;8:e44344     doi: 10.7554/eLife.44344



Abstract:

Transposable elements (TEs) are thought to have helped establish gene regulatory networks. 
Both the embryonic and extraembryonic lineages of the early mouse embryo have seemingly co-opted 
TEs as enhancers, but there is little evidence that they play significant roles in gene regulation. 
Here we tested a set of long terminal repeat TE families for roles as enhancers in mouse embryonic 
and trophoblast stem cells. Epigenomic and transcriptomic data suggested that a large number of TEs 
helped to establish tissue-specific gene expression programmes. Genetic editing of individual TEs 
confirmed a subset of these regulatory relationships. However, a wider survey via CRISPR interference 
of RLTR13D6 elements in embryonic stem cells revealed that only a minority play significant roles 
in gene regulation. Our results suggest that a subset of TEs are important for gene regulation in 
early mouse development, and highlight the importance of functional experiments when evaluating 
gene regulatory roles of TEs.





Notes:
A few abbreviations used throughout these scripts are

REDE (RetroElement Derived Enhancer - referred to as "TE+ Enhancers" in publication)

NEDE (Non-RetroElement Derived Enhancer - referred to as "TE- Enhancers" in publication)

ROI (Regions of Interest)






Scripts used in order:

1.Annotation/plotting of ChIP-seq, DNase-seq & ATAC-seq signals

Rscript ~/todd/Homer_analysis/Homer_annotate_peaks.R 
Rscript ~/todd/Homer_analysis/Plot_complexheatmaps.R 
Rscript ~/todd/Homer_analysis/annotate_trend_lookups.R 


2.Plotting of DNase/ATAC-seq accessibility in various tissues

Rscript ~/todd/DHS/download_ENCODE_DHS.R
Rscript ~/todd/DHS/DNase_across_ENCODE.R
Rscript ~/todd/DHS/DNase_preimplantation.R


3.Motif analysis

Rscript ~/todd/Motif_analysis/Clustal_and_Motif_analysis.R
Rscript ~/todd/Motif_analysis/Plot_Motif_percentage.R
Rscript ~/todd/Motif_analysis/Motif_combination_analysis.R
Rscript ~/todd/Motif_analysis/Motif_analysis_of_FIREWACh.R

4.Analysis of Methylation levels

Rscript ~/todd/Methylation_analysis/TE_methylation_analysis.R

5.Pairing Genes with TEs using HiC

Rscript ~/todd/PCHiC/Bait_lookup_files.R
Rscript ~/todd/PCHiC/Gothic_lookup_files.R
Rscript ~/todd/PCHiC/Pairing_Gothic_Data_to_Expression.R

6.Analysis of Gene expression

Rscript ~/todd/Expression_analysis/Plotting_ESC_TSC_gene_expression.R
Rscript ~/todd/Expression_analysis/download_ENCODE_RNA.R
Rscript ~/todd/Expression_analysis/Multi_tissue_expression_Gingeras_data.R
Rscript ~/todd/Expression_analysis/ESC_differentiation.R
Rscript ~/todd/Expression_analysis/TSC_differentiation.R

7.Guide design for CRISPRi

Rscript ~/todd/CRISPRi/Guide_design/Step1_Repeat_bed_coord_to_fasta.R
Rscript ~/todd/CRISPRi/Guide_design/Step2_Get_All_potential_CRISPR_guides.R
Rscript ~/todd/CRISPRi/Guide_design/Step3_Prep_CasOFF_input.R
~/Cas_OFFinder/cas-offinder	   cas_off_input.txt     G      cas_off_output.txt
Rscript ~/todd/CRISPRi/Guide_design/Step5_CasOFF_Analysis_revised.R
Rscript ~/todd/CRISPRi/Guide_design/Step6_pick_guides.R



