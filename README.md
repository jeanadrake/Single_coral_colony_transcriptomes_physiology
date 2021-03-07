# Single_corol_colony_transcriptomes_physiology
Electronic Supplementary Material Table Legends, Figures, and Figure Legends for: Physiological and transcriptomic variability indicative of differences in key functions within a single coral colony

Jeana L. Drake1+, Assaf Malik1+, Yotam Popovits1+, Oshra Yosef1, Eli Shemesh1, Jarosław Stolarski2, Dan Tchernov1, Daniel Sher1* and Tali Mass1*

1Department of Marine Biology, University of Haifa, 199 Aba Koushy Avenue, Haifa, Israel, 3498838
2Institute of Paleobiology, Polish Academy of Sciences, Twarda 51/55, Warsaw, Poland, 00-818
+These authors contributed equally
*Corresponding authors

Tables
Table S1. Location information for branches used for transcriptional analysis. DESeq2 analysis used the additive model output so that the top-most branch or 9th branch, which points upward and thus does not have a defined cardinal direction, was not included.

Table S2. DE gene count statistics. The counts of over-expressed/under-expressed/unchanged/ NA genes using DESeq2 additive and interaction models. Also included are read count and filtered genes statistics.  Both host Stylophora pistillata and photosymbiont (reads mapped to the Symbiodinium microadriaticum genome, NCBI GCA_001939145.1) data are shown.

Table S3. S. pistillata additive model DE. DESeq2 additive model statistics of DE and gene annotations.

Table S4. S. pistillata interaction model DE. DESeq2 interaction model statistics of DE and gene annotations.

Table S5. Photosymbiont additive model DE. DESeq2 additive model statistics of DE and gene annotations.

Table S6. Dissolved oxygen concentration and total SYBR-stained or chlorophyll auto-fluorescence-based microbial cell counts and measured inside and outside five S. pistillata colonies in situ.

Table S7. Statistical analysis of environmental parameters measured between branches and at colony peripheries, and physiological parameters measured at two tips, a junction, and a base segment on each of five branches distributed throughout the S. pistillata colony.

Table S8. – Host and photosymbiont factor effect Adonis (non-parametric multivariate analysis of variance, known as NPMANOVA or PERMANOVA). The Adonis test, implemented in Vegan, was used to test the effect of factors branch axis, ring height, and cardinal direction, based on NMDS two-dimensional sample clustering.

Table S9. S. pistillata vs . Acropora spp. [1] DE tip vs base statistics for putative orthologous genes found using reciprocal BLAST.

Table S10. S. pistillata KEGG and GO (Trinotate) enrichment based on the DESeq2 additive-model. Enrichment analysis uses the GOSeq Wallenius (Wall) method.

Table S11. S. pistillata KEGG and GO (Trinotate) enrichment based on the DESeq2 interaction-model. Enrichment analysis uses the GOSeq Wallenius (Wall) method.

Table S12. - Photosymbiont GO enrichment (Uniprot database) of DE genes based on the DESeq2 additive-model. Enrichment analysis uses the Gene Set Enrichment Analysis (GSEA) method.

Table S13. Physiological parameters measured along five branches from the same single colony of S. pistillata from which differential gene expression was also examined.

Table S14. DE of putative toxins sequenced in this study obtained by reciprocal blast against known toxins from the literature.

Table S15. DE of orthologous putative toxins between S. pistillata and Acropora spp [1]. Most putative toxins from S. pistillata do not have orthologs in the Acropora spp. dataset.

Table S16. DE of known biomineralization-related genes. 60 genes code for proteins recently sequenced from S. pistillata skeleton [2]. CARPs 1, 2, and 3 [3] and a putative bicarbonate transporter [4] were appended to the list. Not all genes associated with known coral skeletal proteins were sequenced here.

Table S17. DE of orthologous biomineralization-related genes between S. pistillata and Acropora spp [1]. Just under half of such genes from S. pistillata have orthologs in the Acropora spp. dataset.

Supplementary Material Availability
The data underlying this article are available in the article, in its online supplementary material, and on GitHub in the following repository: https://github.com/jeanadrake/Single_corol_colony_transcriptomes_physiology.

Figures
 
Figure S1 – Single nucleotide polymorphism (SNP) identity by state. SNP identity by state heatmap showing the level of SNP identity between samples (0-1.0 indicates 0-100% identity), based on the Broad Institute’s GATK.

 
Figure S2. Interaction heatmap of host genes. The heatmap on the left shows the relative expression values of a subset of S. pistillata genes that demonstrate a significant interaction between base-tip (position factor) versus up-down (ring height factor).  The heatmap on the right represents the log2FC of the interactions. For clarity, only interactions for tips and bases are shown.
 

 
Figure S3. Gene expression patterns of putative toxins as log2FC (a, c), and log10FPKM  (b, d) between branch tips and junctions (a, b) and junctions and bases (c, d). Only significantly DE toxin genes are shown (p<0.05). An uncharacterized gene, for which reciprocal BLAST hits suggest it is a SCRiP (SI Table 14), is noted with a red star (b).
 

a
Query	Best Blast Hit	E-value	% Identity	BLAST Server
Stylophora_pistillata_9725	XP_022777760.1	8 e-122	87.4	NCBI
gene27671	XP_022777760.1	0	100	NCBI
XP_022777760.1	Stylophora_pistillata_9725	7.4 e-126	87.4	Comparative Reef Genomics
gene27671	Stylophora_pistillata_9725	7.4 e-126	87	Comparative Reef Genomics
b
 
Figure S4. Reciprocal protein blast information (a) and sequence alignment (b) for the D-Pocilopotoxin-Spi1 in [5] (Stylophora_pistillata_9725 from Bhattacharya, 2016]) and the gene27671 (rna35563) annotated here as DELTA-stichotoxin-She4a-like (NCBI accession number XP_022777760.1).



 
Figure S5. Gene expression patterns of known biomineralization-related genes as log2FC (a, c), and log10FPKM (b, d) between branch tips and junctions (a, b) and junctions and bases (c, d). Only significantly DE biomineralization genes are shown (p<0.05).
 
 
Figure S6. Percent fluorescence, relative to the injected fluoresceine dye aliquot’s fluorescence, of various volumes of water extracted by syringe either on the opposite side of the branch from the dye injection (a) or at colony peripheries (b). Fluorescence higher than the background water was not observed for extraction volumes of 1.6 ml or smaller.

References
1.	Hemond E.M., Kaluziak S.T., Vollmer S.V. 2014 The genetics of colony form and function in Caribbean Acropora corals. BMC Genomics 15(1), 1133.
2.	Peled Y., Drake J.L., Malik A., Almuly R., Lalzar M., Morgenstern D., Mass T. 2020 Optimization of skeletal protein preparation for LC–MS/MS sequencing yields additional coral skeletal proteins in Stylophora pistillata. BMC Materials 2(1), 8. (doi:10.1186/s42833-020-00014-x).
3.	Mass T., Drake Jeana L., Haramaty L., Kim J.D., Zelzion E., Bhattacharya D., Falkowski Paul G. 2013 Cloning and characterization of four novel coral acid-rich proteins that precipitate carbonates in vitro. Curr Biol 23(12), 1126-1131. (doi:http://dx.doi.org/10.1016/j.cub.2013.05.007).
4.	Zoccola D., Ganot P., Bertucci A., Caminiti-Segonds N., Techer N., Voolstra C.R., Aranda M., Tambutté E., Allemand D., Casey J.R. 2015 Bicarbonate transporters in corals point towards a key step in the evolution of cnidarian calcification. Scientific Reports 5.
5.	Ben-Ari H., Paz M., Sher D. 2018 The chemical armament of reef-building corals: inter-and intra-specific variation and the identification of an unusual actinoporin in Stylophora pistilata. Scientific Reports 8(1), 1-13.

 
Extended Methods for: Physiological and transcriptomic variability indicative of differences in key functions within a single coral colony

Jeana L. Drake1+, Assaf Malik1+, Yotam Popovits1+, Oshra Yosef1, Eli Shemesh1, Jarsoław Stolarski2, Dan Tchernov1, Daniel Sher1* and Tali Mass1*

1Department of Marine Biology, University of Haifa, 199 Aba Koushy Avenue, Haifa, Israel, 3498838
2Institute of Paleobiology, Polish Academy of Sciences, Twarda 51/55, Warsaw, Poland, 00-818
+These authors contributed equally
*Corresponding authors

Field measurements
Sample collection
Ex situ measurements were conducted on a 15 cm diameter colony of Stylophora pistillata collected during the morning in April 2016 at 8 m depth under a special permit from the Israeli Natural Parks Authority in the coral nursery of the Inter-University Institute for Marine Science (IUI), Eilat, Israel using SCUBA during April 2016. The colony was brought to the lab in a tank filled with seawater and divided into individual branches and branch segments within 20 minutes of collection (Figure 1: 3d printed model with all numbered branches labeled). Three branches each were taken from three horizontal concentric rings (Top, Middle, Bottom), for a total of nine branches, and classified according to their relative height within the coral colony. These branches were further divided into sections according to their position on the branch (Tip, Junction, Base) and flash-frozen in TRI-Reagent (Sigma) (Figure 1) for differential gene expression analysis. Five additional branches distributed around the colony were chosen for physiological examination, following the same study parameters and then kept at –80 ˚C. 
Oxygen measurements
Dissolved oxygen concentrations in seawater between branches of five S. pistillata colonies located in the IUI coral nursery were measured in-situ. Measurements were taken as five radial transects per colony, from inside to just outside each colony, using Unisense oxygen microsensors mounted on a Unisense UnderWater Meter System. The microsensors were moved through each colony at a controlled pace in order to keep the intra-colony water as undisturbed as possible and to prevent exchange of inter-branch water within the colony and with the surrounding water. Each measurement was taken for a duration of two minutes to allow the oxygen microsensor readings to stabilize. Dissolved oxygen concentrations were also measured for two minutes at 10 cm above each respective colony.  All measured oxygen concentrations were normalized to 1 ATM (sea surface pressure).  DO concentrations were found to be variable between colonies, with colonies growing at 8.1-9.4 m depth exhibiting significantly higher DO values than those at 6.6-7.7 m depth (electronic supplementary material, table S7); we therefore first standardized all values for each colony to the ‘inside’ value of that colony and then compared relative DO concentrations across colonies.
Bacterial counts
Bacterial counts outside and within the coral colony were performed on samples collected in-situ via radial transects on five different colonies growing at the same site. We conducted preliminary experiments using small aliquots of concentrated fluorescein dye injected to water between branches of S. pistillata skeletons, the same size as that used in our physiological experiments, in a controlled environment to test the efficacy of extracting water from between branches and at colony peripheries with minimal intra-colonial water mixing. Briefly, we injected the dye on one side of a branch and then used a syringe to carefully extract various volumes of water, from 0.4 to 12 ml in intervals of 0.4 ml, from the other side of the branch, with five replicate experiments per extract volume tested. We then measured the fluorescence of the extracted water, subtracted the fluorescence background of the surrounding water, and represented the results as a percentage of the fluorescence from the injected fluorescein dye aliquot. This showed that we could extract up to 1.6 ml of seawater from between coral branches and from colony peripheries without pulling in water from the immediate environment of the surrounding branches (electronic supplementary material, figure S6).  We therefore felt confident extracting 1 ml of seawater from between the branches (‘in’ position)  and the tip of the branch at the colony periphery (‘out’ position).of the live colonies. Upon completion of the dive, all water samples were immediately carried to the lab and fixed in cryovials in a final concentration of 0.0625% glutaraldehyde, incubated in the dark for 10 minutes, and then flash-frozen at − 80 °C until analysis by flow cytometry.
 	Just prior to analysis, bacterial count samples were thawed at room temperature and diluted 1:1 in ultra-pure water, and 2 μm-diameter fluorescent beads were added to each sample as an internal standard (Polysciences, Warminster, PA, USA). The samples were analyzed on a FACSCanto™ II Flow Cytometry Analyzer Systems (BD Biosciences). We first examined the natural fluorescence of the cells (chlorophyll and phycoerythrin pigments). We then stained the cells with SYBR Green I (Molecular Probes/ThermoFisher) according to the manufacturer’s instructions and counted the total microbial population as well as cyanobacterial sub-populations. Data were acquired and processed with FlowJo software. Flow rates were determined several times during each run and the average value for a sample-free test run was used for calculating cell per ml. 
Physiological analysis
Tissue removal 
Tissue was removed from the frozen coral fragment skeleton by airbrush with phosphate buffered saline (PBS) and homogenized using an electrical homogenizer (MRC, HOG-160-1/2) for 10 seconds. A sub-sample of each fragment’s homogenate was stored for cell counts, chlorophyll extraction, and total protein measurements and the rest of the homogenate was centrifuged for 5 minutes at 5000 rpm with a Lumintron Eppendorf Centrifuge.  The supernatant was retained for host protein and hemolysis assay. The pellet was observed under the microscope and found to consist primarily of discharged and non-discharged nematocysts and photosynthetic symbionts.  All quantities were adjusted to fragment extraction volume and normalized to cm2 of the fragment’s surface area.
Chlorophyll extraction 
Two ml of tissue homogenate was filtered on a Whatman GF/C filter and incubated with 1 ml 90% acetone for two hours at 4°C. After incubation, the filter was manually homogenized and the solution was filtered through a 0.22 µm syringe filter into a glass cuvette. Spectrophotometric measurements were conducted on a NanoDrop (Thermo-Fisher) and chlorophyll a concentrations were calculated from light absorbance results based on the following equation [1]: chl-a [mg/ml]= -0.3319(ABS630)-1.7485(ABS647)+11.9442(ABS664)-1,4306(ABS691).
Nematocysts and photosynthetic symbiont counts   
100 µl of tissue homogenate was used for nematocysts and zooxanthellae counts by hemocytometer (BOECO, Germany); cells were counted in five randomly selected fields per fragment (1 mm2 each) on a Nikon Eclipse Ti-S Inverted Microscope System. Nematocysts were counted under white light with a differential contrast filter whereas photosynthetic symbionts were viewed by fluorescence at 440 nm excitation and 590 nm emission. NIS ELEMENTS (Nikon) software was used for the cell counts with automatic counting settings limited to cell diameters smaller than 15 µm and circularity set to >0.5 but <1.
Total protein and host protein 
We quantified homogenate and host total protein; sonication (Ultrasonic Atomizer Probe Sonicator) was used for further extraction of symbiont protein. Protein concentration was measured by bichronoic acid (BCA) assay (Pierce) against a bovine serum albumin standard curve at 540 nm light. 
Hemolysis assay 
Hemolytic assays were performed against O Rh positive human blood cells obtained from the Rambam/Yoseftal Hospital Blood Bank, as described previously [2]. Briefly, 2 ml of whole blood was diluted to a final volume of 15 ml in PBS (pH 7.4) and centrifuged at 3000 g for 5 min. The supernatant was removed, was the process was repeated with the pelleted erythrocytes until the supernatant was clear. Washed erythrocytes were then resuspended in PBS to a final concentration of 20% (v/v). For each hemolysis assay, 160 μl of crude coral tissue extract was incubated with 40 μl of washed erythrocyte suspension (4% erythrocytes v/v) at 37°C for 30 minutes in a water bath. At the end of the incubation, 400 μl of PBS was added, and the assays were centrifuged at 3000 g for 3 min. The supernatant fluid containing the hemoglobin released from lysed erythrocytes was transferred to 96-well microplates and the absorbance at 540 nm was determined by using a spectrophotometric microplate reader (Perkin-Elmer). In addition, each experiment was normalized to a positive (100% hemolysis) and negative control (0% hemolysis) by incubating erythrocytes with DDW and PBS alone, respectively. HU50 was defined as the amount of coral protein required to cause 50% hemolysis [3] in a dilution series of the coral extract. 
Fragment surface area 
Surface areas were estimated using the aluminum foil method [4]. Briefly, aluminum foil was molded over the skeleton of each fragment (without measuring skeleton connective areas between the divided fragments), carefully removed, and weighed. This was repeated three times on every fragment and the surface area was estimated using a standard curve of the derived relationship between foil area and weight. 
Skeleton cleansing for assessments of surface area and corallite measurements 
After tissue extraction, all fragments were cleaned of organic residues overnight in 3% sodium hypochlorite; the skeletons were then washed in deionized water and dried at 55°C for three hours.
Micromorphological analysis
The distance between neighboring polyps and the surface area of the polyp was measured using a Nikon binocular microscope and analyzed using NIS-Elements (Nikon) software. For each coral fragment, at least four corallites located in a position parallel to the focal plane, were chosen. Major and minor axes and circularity were measured for each corallite and the ratio of the two axes (major divided by the minor) were calculated; a ratio of 1 indicates a circle and the larger the ratio, the more elliptical the shape. Surface area of each measured corallite was calculated from the two-measured axes, using the following equation: S = π∗ r (major axes)* r (minor axes). Finally, distances between neighboring corallites were also measured to evaluate polyp density in each area.
	
Differential Gene Expression Analysis
RNA extraction, processing, and sequencing
Total RNA was extracted from the holobiont in fragments stored in TRI-Reagent (Sigma) following the manufacturer's protocol, with some modification at the homogenization step. Briefly, samples frozen in TRI-Reagent were heated and centrifuged, and then bromochloropropane, at a ratio of 1:10, was added to the samples for separation. After incubation at room temperature for 10 minutes and centrifugation, RNA was purified from the clear phase using a Purelink RNA Mini Kit (Ambion) according to the manufacturer’s protocol. The RNA was washed in 70% ethanol and then on-column DNase digestion was performed using a Qiagen RNase-free DNaseI Kit according to the manufacturer's protocol. RNA-seq libraries were prepared using an in-house protocol at the Weizmann Institute of Science. Briefly, the polyA fraction (mRNA) was purified from 500 ng of total RNA per sample following by fragmentation and generation of double-stranded cDNA. Then, end repair, A base addition, adapter ligation and PCR amplification steps were performed. Libraries were evaluated by Qubit (Thermo Fisher) and TapeStation (Agilent). Sequencing libraries were constructed with barcodes to allow multiplexing of all samples to be run in each lane. Approximately 578 million total high-quality 125 bp paired end reads (15.45 ± 1.5 million paired reads per sample) were sequenced on an Illumina HiSeq 2500 across three different lanes (i.e., each sample run in triplicate to remove batch effects). 
Host and symbiont RNA-Seq differential expression analysis 
Standard RNA-Seq quality filtering were conducted as described previously [5]. Briefly, RNA-Seq reads of were adapter-trimmed using cutadapt 1.15 (https://cutadapt.readthedocs.io), and then low-quality regions were removed with Trimmomatic 0.3 [6]. In order to detect correct taxonomic classifications within our transcriptome data, high-quality reads were mapped to all available NCBI and Reefgenomics (reefgenomics.org; [7]) genome-based proteomes databases of Symbiodiniaceae  species Symbiodinium microadraticum, Cladocopium goreaui, and Fugacium kawagutii (formerly Symbiodinium spp. clades A, C1, and F, respectively [8]), and Cnidaria as well as selected stramenopiles/alveolates/Rhizaria and Metazoa databases, using Diamond [9]. Top hits almost exclusively belonged to robust corals (mainly S. pistillata), and S. microadriaticum (formerly “clade A” [8]). These results are expected since S. microadriaticum was found to be a predominant symbiont in shallow water in Eilat [5]. For symbiont RNA-Seq mapping, we further aligned the reads to the merged host genome assembly (NCBI GCA_002571385.1) and the symbiont genome assembly (NCBI GCA_001939145.1) using STAR [10-12]. Overall, 72-80% (~11-16 million reads) were concordantly-mapped to the host genome, and 2%-10% (~0.5-1.5 million reads) were concordantly-mapped to the symbiont genome. In total, 83-84% of the reads were aligned to either host or symbiont genomes. Differential expression (DE) analysis was conducted using Bioconductor DEseq2 [6], separately for the host and the symbiont genes, using a DEseq2 generalized linear model. Specifically, we tested: (1) the additive effect of branch position, ring height and compass direction; (2) the additive and interaction effect of branch position and ring height (excluding compass direction factor); (3) the effect of branch position and ring height as a single factor (SI Table ###, combined factor column, e.g.: “tip.down”, “tip.up”, “tip.top”, etc …). Note that the first model cannot include the four top branch samples, since top branch samples are represented by an incomparable compass direction (“center”); as such, these samples were excluded from the analysis of expression of putative toxins and known biomineralization-related genes.  For the symbiont analysis, only the first model was used. Although DESeq2 is designed to analyze samples of variable total read counts, in order to avoid read normalization problems we only selected symbiont genes whose 25% quantile of read count was at least 12 reads (with 12,735 symbiont genes remaining); further, only symbiont samples with >300k total reads were retained. After filtering, the correct symbiont reads normalization was indicated by a typical clost-to-symmetric distribution of fold-changes around the DESeq2 MA plot horizontal axis. NMDS analysis was conducted using metaDMA in the R Vegan package based on log10FPM values of all expressed genes. For the symbiont NMDS, since read counts were relatively low, pairs of tips from the same branch are represented together based on the sum of read counts for both tips.
Branch genetics and SNPs analysis
Single nucleotide polymorphisms (SNPs) analysis was conducted for each sample using the Broad Institute-recommended RNA-Seq SNPs practice (https://gatk.broadinstitute.org). The pipeline includes four steps: STAR reads mapping, a pre-processing step aimed at removing various alignment biases, variant calling using GATK version 3.5, and finally variant filtration and annotation [13]. In order to exclude variant-call biases due to changes in transcript abundance between samples, we considered only SNP loci with coverage greater than 30 RNA-Seq reads in all tested samples. Identity By Decent (IBD) analysis was conducted using SNPRealte in R.
Host and symbiont functional enrichment analysis 
Biological terms were assigned to genes based on Uniprot S. pistillata (www.uniprot.org), KEGG S. pistillata (www.kegg.jp), S. pistillata Trinotate annotations [14], and Uniprot S. microadriaticum databases. Enrichment analysis was conducted in Bioconductor GOSeq [15], which corrects for enrichment biases associated with correlations between gene size and DE significance. We also searched for functional enrichment using the score-based tool GSEA [16, 17], which may be more sensitive than p-value cutoff-based searches (such as GOSeq) when small non-significant changes in relatively large groups of genes are expected. In GSEA, we used log2 fold changes as scores. Because many enriched terms are functionally related, we hierarchically clustered the biological terms based on pairwise distances between groups of genes: D = |xa∩xb|/minimum(|a|,|b|), where a,b are two sets of gene ids, and xa and xa are DE genes ids from a and b, respectively. Terms trees were constructed using Bioconductor-ggtree [18].
Focused S. pistillata gene analyses
Differential expression of genes coding for known biomineralization-related proteins from S. pistillata [5, 19-22] skeleton  or the calicoblastic layer were examined (electronic supplementary material table S16). Putative toxin genes were identified using reciprocal best-BLAST hits as described in [23, 24]. 

Physiological Statistics 
Statistical tests of all skeletal, physiological, and environmental data were performed in RStudio [25]. If a parameter displayed both homogeneity of variance (Fligner-Killeen test) and normal distribution (Shapiro-Wilke test), either a t-test or a one-way ANOVA was performed, with TukeyHSD post hoc analysis for parameters with significant differences between branch locations. Parameters with non-homogenous variance and/or non-normal distribution were analyzed using Wilcoxon’s rank sum test or the Kruskal-Wallis rank sum test, with Dunn’s test for post hoc analysis for parameters with significant differences between branch locations. 

Electronic Supplementary Material Availability, Tables
All electronic supplementary tables are available in the Google Drive folder: https://drive.google.com/drive/folders/1jU-9fjZh5Qtr05w7tA1A47WayUGrs6dh?usp=sharing. Following receipt of a link to upload the data to the Dryad data repository, the dataset will be available under the dataset https://doi.org/10.5061/dryad.XXXX.

References
1.	Ritchie R. 2008 Universal chlorophyll equations for estimating chlorophylls a, b, c, and d and total chlorophylls in natural assemblages of photosynthetic organisms using acetone, methanol, or ethanol solvents. Photosynthetica 46(1), 115-126.
2.	Primor N., Zlotkin E. 1975 On the ichthyotoxic and hemolytic action of the skin secretion of the flatfish Pardachirus marmoratus (Soleidae). Toxicon 13(4), 227-231.
3.	Bartosz G., Finkelshtein A., Przygodzki T., Bsor T., Nesher N., Sher D., Zlotkin E. 2008 A pharmacological solution for a conspecific conflict: ROS-mediated territorial aggression in sea anemones. Toxicon 51(6), 1038-1050.
4.	Marsh Jr J.A. 1970 Primary productivity of reef‐building calcareous red algae. Ecology 51(2), 255-263.
5.	Peled Y., Drake J.L., Malik A., Almuly R., Lalzar M., Morgenstern D., Mass T. 2020 Optimization of skeletal protein preparation for LC–MS/MS sequencing yields additional coral skeletal proteins in Stylophora pistillata. BMC Materials 2(1), 8. (doi:10.1186/s42833-020-00014-x).
6.	Love M.I., Huber W., Anders S. 2014 Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15(12), 1-21.
7.	Liew Y.J., Aranda M., Voolstra C.R. 2016 Reefgenomics. Org-a repository for marine genomics data. Database 2016.
8.	LaJeunesse T.C., Parkinson J.E., Gabrielson P.W., Jeong H.J., Reimer J.D., Voolstra C.R., Santos S.R. 2018 Systematic revision of Symbiodiniaceae highlights the antiquity and diversity of coral endosymbionts. Curr Biol 28(16), 2570-2580. e2576.
9.	Buchfink B., Xie C., Huson D.H. 2015 Fast and sensitive protein alignment using DIAMOND. Nat Methods 12(1), 59-60.
10.	Dobin A., Davis C.A., Schlesinger F., Drenkow J., Zaleski C., Jha S., Batut P., Chaisson M., Gingeras T.R. 2013 STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29(1), 15-21.
11.	Voolstra C.R., Li Y., Liew Y.J., Baumgarten S., Zoccola D., Flot J.-F., Tambutté S., Allemand D., Aranda M. 2017 Comparative analysis of the genomes of Stylophora pistillata and Acropora digitifera provides evidence for extensive differences between species of corals. Scientific Reports 7(1), 17583.
12.	Aranda M., Li Y., Liew Y.J., Baumgarten S., Simakov O., Wilson M.C., Piel J., Ashoor H., Bougouffa S., Bajic V.B. 2016 Genomes of coral dinoflagellate symbionts highlight evolutionary adaptations conducive to a symbiotic lifestyle. Scientific Reports 6(1), 1-15.
13.	DePristo M.A., Banks E., Poplin R., Garimella K.V., Maguire J.R., Hartl C., Philippakis A.A., Del Angel G., Rivas M.A., Hanna M. 2011 A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nat Genet 43(5), 491.
14.	Bryant D.M., Johnson K., DiTommaso T., Tickle T., Couger M.B., Payzin-Dogru D., Lee T.J., Leigh N.D., Kuo T.-H., Davis F.G. 2017 A tissue-mapped axolotl de novo transcriptome enables identification of limb regeneration factors. Cell Reports 18(3), 762-776.
15.	Young M.D., Wakefield M.J., Smyth G.K., Oshlack A. 2012 goseq: Gene Ontology testing for RNA-seq datasets. R Bioconductor 8, 1-25.
16.	Sergushichev A.A. 2016 An algorithm for fast preranked gene set enrichment analysis using cumulative statistic calculation. BioRxiv, 060012.
17.	Subramanian A., Tamayo P., Mootha V.K., Mukherjee S., Ebert B.L., Gillette M.A., Paulovich A., Pomeroy S.L., Golub T.R., Lander E.S. 2005 Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences 102(43), 15545-15550.
18.	Yu G., Lam T.T.-Y., Zhu H., Guan Y. 2018 Two methods for mapping and visualizing associated data on phylogeny using ggtree. Mol Biol Evol 35(12), 3041-3043.
19.	Drake J.L., Mass T., Haramaty L., Zelzion E., Bhattacharya D., Falkowski P.G. 2013 Proteomic analysis of skeletal organic matrix from the stony coral Stylophora pistillata. Proceedings of the National Academy of Sciences 110(10), 3788-3793.
20.	Mass T., Drake Jeana L., Haramaty L., Kim J.D., Zelzion E., Bhattacharya D., Falkowski Paul G. 2013 Cloning and characterization of four novel coral acid-rich proteins that precipitate carbonates in vitro. Curr Biol 23(12), 1126-1131. (doi:http://dx.doi.org/10.1016/j.cub.2013.05.007).
21.	Puverel S., Tambutté E., Pereira-Mouries L., Zoccola D., Allemand D., Tambutté S. 2005 Soluble organic matrix of two Scleractinian corals: Partial and comparative analysis. Comparative Biochemistry and Physiology Part B 141, 480-487.
22.	Zoccola D., Ganot P., Bertucci A., Caminiti-Segonds N., Techer N., Voolstra C.R., Aranda M., Tambutté E., Allemand D., Casey J.R. 2015 Bicarbonate transporters in corals point towards a key step in the evolution of cnidarian calcification. Scientific Reports 5.
23.	Yosef O., Popovits Y., Malik A., Ofek-Lalzer M., Mass T., Sher D. 2020 A tentacle for every occasion: comparing the hunting tentacles and sweeper tentacles, used for territorial competition, in the coral Galaxea fascicularis. BMC Genomics 21(1), 1-16.
24.	Rachamim T., Morgenstern D., Aharonovich D., Brekhman V., Lotan T., Sher D. 2015 The dynamically evolving nematocyst content of an anthozoan, a scyphozoan, and a hydrozoan. Mol Biol Evol 32(3), 740-753.
25.	Team R. 2019 RStudio: Integrated Development for R.  (Boston, MA, RStudio, Inc.
