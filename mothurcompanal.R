#Kamimura, B. A. et al, 2020 

#Not paired end sequencing so no make.contigs command
Mothur> make.file(stab.file, typ=fastq, prefix=stab)

#filter sequences
##removesequences with ambig bases, shorter than 2.5%, and longer than 97.5%
Mothur>screen.seqs(fasta= stab.file.fasta, group=stab.groups, summary= stab.file.summary,  maxlength=305, maxambig=0)

#remove identical (grouped) sequences; representative sequences were picked and stored in fasta file, the other sequences are saved as sequence names to reduce the computational power/ work needed.

Mothur>unique.seqs(fasta= stab.good.file)

#count table created of the unique sequences
Mothur>count_seqs(name= stab.file.good.names, group= stab.good.groups)

#align the sequences to the silva reference database
Mothur>align.seqs(fasta= stab.good.unique.file, reference=silva.bacteria.fasta)

#summary statistics to get descriptive stats of alignment
Mothur>summary.seqs(fasta= stab.good.unique.align, count= stab.good.count_table)

#Remove sequences that were before and after where the sequences aligned for the alignment site alignment site for V3-V4 region of 16s
Mothur>screen.seqs(fasta= stab.good.unique.align, count= stab.good.count_table, summary= stab. good.unique.summary, start=6332, end=25316, maxahomop=8)
#remove overhangs, remove alignment characters that only have “-“, using vertical=T and remove the sequences that contain “.”
Mothur>filter.seqs(fasta= stab.good.unique.good.align, vertical=T, trump=.)

#run unique seqs again to check for new reducndant sequences
Mothur>unique.seqs(fasta= stab.good.unique.good.filter.fasta, count= stab.good.unique.good.filter.count_table)

#The sequences then need to be de-noised with a 2 threshold of mismatches in the sequence
Mothur>pre.cluster(fasta= stab.good.unique.good.filter.unique.fasta, count= stab.good.unique.good.filter.count_table, diffs=2)

#Read fasta and countfiles to chimera sequences and then remove chimera 

Mothur>chimer.vsearch(fasta= stab.good.unique.good.filter.unique.precluster.fasta, count= stab.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

Mothur>remove.seqs( fasta= stab. good.unique.good.filter.unique.precluster.fasta, accnos= stab.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

#Assign taxonomy
Mothur>classify.seqs(fasta= stab.good.unique.good.filter.unique.precluster.pick.fasta, count= stab.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=transient16_022016.pds.align, taxonomy= transient16_022016.pds.tax)

#Remove  Chloroplast, Mitochondria, Unknown, Archaea., and Eukaryotes 
Mothur>remove.lineage(fasta= stab.good.unique.good.filter.unique.precluster.pick.fasta, count= stab.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy= stab.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon= Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

#Calculate uncorrected pairwise distances between aligned DNA sequences. By default gaps are penalized; cutoff value indicate s that distances larger that 0.03 (>97%) won’t be saved
Mothur>dist.seqs(fasta=stab.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)

#Assign sequences to OTUs using default opticlsut method
Mothur>cluster(column= stab.good.unique.good.filter.unique.precluster.pick.pick.dist, count= stab.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

#Determine how many sequences are in each OTU at the 0.03 cutoff level.
Mothur>make.shared(list= stab.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count= stab.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

#Determine taxonomy for all OTUs
Mothur>classify.otus(list=stab.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stab.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy= stab.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy label=0.03)







#Zwirzitz, B. et al., 2020 

#make .files to include sample and paired end sequences  sequence associated
Mothur> make.file(bstability.file, typ=fastq, prefix=stability)

#combine paired end reads together
Mothur>make.contigs(file=stability.files)

#summary statistics to get descriptive stats of contigs 
Mothur>summary.seqs(fasta=stability.file.trim.contigs.fasta)

#filter sequences
##removesequences with ambig bases, shorter than 2.5%, and longer than 97.5%
Mothur>screen.seqs(fasta= stability.file.trim.contigs.fasta, group=stability.contigs.groups, summary= stability.trim.contigs.summary, minlength=391, maxlength=587, maxambig=0)

#remove identical (grouped) sequences; representative sequences were picked and stored in fasta , the other sequences are saved as sequence names to reduce the computational power/ work needed.

Mothur>unique.seqs(fasta= stability.trim.contigs.good.fasta)

#count table created of the unique sequences
Mothur>count_seqs(name= stability.trim.contigs.good.names, group= stability.contigs.good.groups)

#align the sequences to the silva reference database
Mothur>align.seqs(fasta= stability.trim.contigs.good.unique.fasta, reference=silva.bacteria.fasta)

#summary statistics to get descriptive stats of alignment
Mothur>summary.seqs(fasta= stability.trim.contigs.good.unique.align, count== stability.trim.contigs.good.count_table)

#Remove sequences that were before and after where the sequences aligned for the alignment site alignment site for V3-V4 region of 16s
Mothur>screen.seqs(fasta= stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary= stability.trim.contigs.good.unique.summary, start=6388, end=25316, maxahomop=8)

#remove overhangs, remove alignment characters that only have “-“, using vertical=T and remove the sequences that contain “.”
Mothur>filter.seqs(fasta= stabilityility.trim.contigs.good.unique.good.align, vertical=T, trump=.)

#run unique seqs again to check for new reducndant sequences
Mothur>unique.seqs(fasta= stability.trim.contigs.good.unique.good.filter.fasta, count= stability.trim.contigs.good.good.filter. count_table)

#The sequences then need to be de-noised with a 2 threshold of mismatches in the sequence
Mothur>pre.cluster(fasta= stability.trim.contigs.good.unique.good.filter.unique.fasta, count= stability.trim.contigs.good.unique.good.filter.count_table, diffs=2)

#Read fasta and counts to chimera sequences and then remove chimera 

Mothur>chimer.vsearch(fasta= stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count= stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

Mothur>remove.seqs( fasta= stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos= ab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

#Assign taxonomy
Mothur>classify.seqs(fasta= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=transient16_022016.pds.align, taxonomy= transient16_022016.pds.tax)

#Remove  Chloroplast, Mitochondria, Unknown, Archaea., and Eukaryotes 
Mothur>remove.lineage(fasta= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon= Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

#Calculate uncorrected pairwise distances between aligned DNA sequences. By default gaps are penalized; cutoff value indicate s that distances larger that 0.03 (>97%) won’t be saved
Mothur>dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)

#Assign sequences to OTUs using default opticlsut method
Mothur>cluster(column= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count= stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

#Determine how many sequences are in each OTU at the 0.03 cutoff level.
Mothur>make.shared(list= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count= stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

#Determine taxonomy for all OTUs
Mothur>classify.otus(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy label=0.03)



#Guzmán, A., et al., 2020
#make .files to include sample and paired end sequences  sequence associated
Mothur> make.file(gstabilityility.file, typ=fastq, prefix=gstability)

#combine paired end reads together
Mothur>make.contigs(file=gstability.files)

#summary statistics to get descriptive stats of contigs 
Mothur>summary.seqs(fasta=gstability.file.trim.contigs.fasta)

#filter sequences
##removesequences with ambig bases, shorter than 2.5%, and longer than 97.5%
Mothur>scree.seqs(fasta= gstability.file.trim.contigs.fasta, group=gstability.contigs.groups, summary= gstability.file.trim.contigs.summary, minlength=242, maxlength=310, maxambig=0)

#remove identical (grouped) sequences; representative sequences were picked and stored in fasta file, the other sequences are saved as sequence names to reduce the computational power/ work needed.

Mothur>unique.seqs(fasta= gstability.file.trim.contigs.good.fasta)

#count table created of the unique sequences
Mothur>count_seqs(name= gstability.file.trim.contigs.good.names, group= stability.contigs.good.groups)

#align the sequences to the silva reference database
Mothur>align.seqs(fasta= gstability.file.trim.contigs.good.unique.fasta, reference=gsilva.bacteria.fasta)

#summary statistics to get descriptive stats of alignment
Mothur>summary.seqs(fasta= gstability.file.trim.contigs.good.unique.align, count== gstability.file.trim.contigs.good.count_table)

#Remove sequences that were before and after where the sequences aligned for the alignment site alignment site for V3-V4 region of 16s
Mothur>screen.seqs(fasta= gstability.file.trim.contigs.good.unique.align, count=gstability.file.trim.contigs.good.count_table, summary= gstability.file.trim.contigs.good.unique.summary, start=10358, end=25318, maxahomop=8)

#remove overhangs, remove alignment characters that only have “-“, using vertical=T and remove the sequences that contain “.”
Mothur>filter.seqs(fasta= gstabilityility.file.trim.contigs.good.unique.good.align, vertical=T, trump=.)

#run unique seqs again to check for new reducndant sequences
Mothur>unique.seqs(fasta= gstability.file.trim.contigs.good.unique.good.filter.fasta, count= gstability.file.trim.contigs.good.good.filter. count_table)

#The sequences then need to be de-noised with a 2 threshold of mismatches in the sequence
Mothur>pre.cluster(fasta= gstability.file.trim.contigs.good.unique.good.filter.unique.fasta, count= gstability.file.trim.contigs.good.unique.good.filter.count_table, diffs=2)

#Read fasta and countfiles to chimera sequences and then remove chimera 

Mothur>chimer.vsearch(fasta= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

Mothur>remove.seqs( fasta= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

#Assign taxonomy
Mothur>classify.seqs(fasta= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=transient16_022016.pds.align, taxonomy= transient16_022016.pds.tax)

#Remove  Chloroplast, Mitochondria, Unknown, Archaea., and Eukaryotes 
Mothur>remove.lineage(fasta= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon= Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

#Calculate uncorrected pairwise distances between aligned DNA sequences. By default gaps are penalized; cutoff value indicate s that distances larger that 0.03 (>97%) won’t be saved
Mothur>dist.seqs(fasta=gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)

#Assign sequences to OTUs using default opticlsut method
Mothur>cluster(column= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

#Determine how many sequences are in each OTU at the 0.03 cutoff level.
Mothur>make.shared(list= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

#Determine taxonomy for all OTUs
Mothur>classify.otus(list=gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy= gstability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy label=0.03)


Kang, S., et al., 2019

Mothur> make.file(stab.file, typ=fastq, prefix=bstability)

#combine paired end reads together
Mothur>make.contigs(file=stab.files)

#summary statistics to get descriptive stats of contigs 
Mothur>summary.seqs(fasta=stab.file.trim.contigs.fasta)

#filter sequences
##removesequences with ambig bases, shorter than 2.5%, and longer than 97.5%
Mothur>screen.seqs(fasta= stab.file.trim.contigs.fasta, group=stab.contigs.groups, summary= stab.file.trim.contigs.summary, minlength=252, maxlength=254, maxambig=0)

#remove identical (grouped) sequences; representative sequences were picked and stored in fasta file, the other sequences are saved as sequence names to reduce the computational power/ work needed.

Mothur>unique.seqs(fasta= stab.file.trim.contigs.good.fasta)

#count table created of the unique sequences
Mothur>count_seqs(name= stab.file.trim.contigs.good.names, group= stability.contigs.good.groups)

#align the sequences to the silva reference database
Mothur>align.seqs(fasta= stab.file.trim.contigs.good.unique.fasta, reference=silva.bacteria.fasta)

#summary statistics to get descriptive stats of alignment
Mothur>summary.seqs(fasta= stab.file.trim.contigs.good.unique.align, count== stab.file.trim.contigs.good.count_table)

#Remove sequences that were before and after alignment site for V4 region of 16s
Mothur>screen.seqs(fasta= stab.file.trim.contigs.good.unique.align, count=stab.file.trim.contigs.good.count_table, summary= stab.file.trim.contigs.good.unique.summary, start=13862, end=23444, maxahomop=8)

#remove overhangs, remove alignment characters that only have “-“, using vertical=T and remove the sequences that contain “.”
Mothur>filter.seqs(fasta= stab.file.trim.contigs.good.unique.good.align, vertical=T, trump=.)

#run unique seqs again to check for new reducndant sequences
Mothur>unique.seqs(fasta= stab.file.trim.contigs.good.unique.good.filter.fasta, count= stab.file.trim.contigs.good.good.filter.count_table)

#The sequences then need to be de-noised with a 2 threshold of mismatches in the sequence
Mothur>pre.cluster(fasta= stab.file.trim.contigs.good.unique.good.filter.unique.fasta, count= stab.file.trim.contigs.good.unique.good.filter.count_table, diffs=2)

#Read fasta and countfiles to chimera sequences and then remove chimera 

Mothur>chimer.vsearch(fasta= stab.file.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count= stab.file.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

Mothur>remove.seqs( fasta= stab.file.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos= stab.file.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

#Assign taxonomy
Mothur>classify.seqs(fasta= stab.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= stab.file.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=transient16_022016.pds.align, taxonomy= transient16_022016.pds.tax)

#Remove  Chloroplast, Mitochondria, Unknown, Archaea., and Eukaryotes 
Mothur>remove.lineage(fasta= stability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= stability.file.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy= stability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon= Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

#Calculate uncorrected pairwise distances between aligned DNA sequences. By default gaps are penalized; cutoff value indicate s that distances larger that 0.03 (>97%) won’t be saved
Mothur>dist.seqs(fasta=stab.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)

#Assign sequences to OTUs using default opticlsut method
Mothur>cluster(column= stab.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count= stab.file.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

#Determine how many sequences are in each OTU at the 0.03 cutoff level.
Mothur>make.shared(list= stab.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count= stab.file.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

#Determine taxonomy for all OTUs
Mothur>classify.otus(list=stab.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stab.file.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy= stab.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy label=0.03)


#Wang, X., et al., 2018
Mothur> make.file(stab.file, typ=fastq, prefix=stab)

#combine paired end reads together
Mothur>make.contigs(file=stab.files)

#summary statistics to get descriptive stats of contigs 
Mothur>summary.seqs(fasta=stab.file.trim.contigs.fasta)

#filter sequences
##removesequences with ambig bases, shorter than 2.5%, and longer than 97.5%
Mothur>screen.seqs(fasta= stab.file.trim.contigs.fasta, group=stab.contigs.groups, summary= stab.file.trim.contigs.summary, minlength=358, maxlength=473, maxambig=0)

#remove identical (grouped) sequences; representative sequences were picked and stored in fasta file, the other sequences are saved as sequence names to reduce the computational power/ work needed.

Mothur>unique.seqs(fasta= stab.trim.contigs.good.fasta)

#count table created of the unique sequences
Mothur>count_seqs(name= stab.trim.contigs.good.names, group= stability.contigs.good.groups)

#align the sequences to the silva reference database
Mothur>align.seqs(fasta= stab.trim.contigs.good.unique.fasta, reference=silva.bacteria.fasta)

#summary statistics to get descriptive stats of alignment
Mothur>summary.seqs(fasta= stab.trim.contigs.good.unique.align, count== stab.trim.contigs.good.count_table)

#Remove sequences that were before and after where the sequences aligned for the alignment site  alignment site for V3-V4 region of 16s
Mothur>screen.seqs(fasta= stab.trim.contigs.good.unique.align, count=stab.trim.contigs.good.count_table, summary= stab.trim.contigs.good.unique.summary, start=6388, end=25316, maxahomop=8)

#remove overhangs, remove alignment characters that only have “-“, using vertical=T and remove the sequences that contain “.”
Mothur>filter.seqs(fasta= stab.file.trim.contigs.good.unique.good.align, vertical=T, trump=.)

#run unique seqs again to check for new reducndant sequences
Mothur>unique.seqs(fasta= stab.trim.contigs.good.unique.good.filter.fasta, count= stab.trim.contigs.good.good.filter.count_table)

#The sequences then need to be de-noised with a 2 threshold of mismatches in the sequence
Mothur>pre.cluster(fasta= stab.trim.contigs.good.unique.good.filter.unique.fasta, count= stab.trim.contigs.good.unique.good.filter.count_table, diffs=2)

#Read fasta and countfiles to chimera sequences and then remove chimera 

Mothur>chimer.vsearch(fasta= stab.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count= stab.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

Mothur>remove.seqs( fasta= stab.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos= stab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

#Assign taxonomy
Mothur>classify.seqs(fasta= stab.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= stab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=transient16_022016.pds.align, taxonomy= transient16_022016.pds.tax)

#Remove  Chloroplast, Mitochondria, Unknown, Archaea., and Eukaryotes 
Mothur>remove.lineage(fasta= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon= Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

#Calculate uncorrected pairwise distances between aligned DNA sequences. By default gaps are penalized; cutoff value indicate s that distances larger that 0.03 (>97%) won’t be saved
Mothur>dist.seqs(fasta=stab.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)

#Assign sequences to OTUs using default opticlsut method
Mothur>cluster(column= stab.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count= stab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

#Determine how many sequences are in each OTU at the 0.03 cutoff level.
Mothur>make.shared(list= stab.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count= stab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

#Determine taxonomy for all OTUs
Mothur>classify.otus(list=stab.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy= stab.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy label=0.03)


#Tan, X. et al., 2019
#combine paired end reads together
Mothur>make.contigs(file=stability.files)

#summary statistics to get descriptive stats of contigs 
Mothur>summary.seqs(fasta=stability.file.trim.contigs.fasta)

#filter sequences
##removesequences with ambig bases, shorter than 2.5%, and longer than 97.5%
Mothur>screen.seqs(fasta= stability.file.trim.contigs.fasta, group=stability.contigs.groups, summary= stability.trim.contigs.summary, minlength=246, maxlength=330, maxambig=0)

#remove identical (grouped) sequences; representative sequences were picked and stored in fasta , the other sequences are saved as sequence names to reduce the computational power/ work needed.

Mothur>unique.seqs(fasta= stability.trim.contigs.good.fasta)

#count table created of the unique sequences
Mothur>count_seqs(name= stability.trim.contigs.good.names, group= stability.contigs.good.groups)

#align the sequences to the silva reference database
Mothur>align.seqs(fasta= stability.trim.contigs.good.unique.fasta, reference=silva.bacteria.fasta)

#summary statistics to get descriptive stats of alignment
Mothur>summary.seqs(fasta= stability.trim.contigs.good.unique.align, count== stability.trim.contigs.good.count_table)

#Remove sequences that were before and after where the sequences aligned for the alignment site alignment site for V4 region of 16s
Mothur>screen.seqs(fasta= stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary= stability.trim.contigs.good.unique.summary, start=11895, end=25318, maxahomop=8)

#remove overhangs, remove alignment characters that only have “-“, using vertical=T and remove the sequences that contain “.”
Mothur>filter.seqs(fasta= stabilityility.trim.contigs.good.unique.good.align, vertical=T, trump=.)

#run unique seqs again to check for new reducndant sequences
Mothur>unique.seqs(fasta= stability.trim.contigs.good.unique.good.filter.fasta, count= stability.trim.contigs.good.good.filter. count_table)

#The sequences then need to be de-noised with a 2 threshold of mismatches in the sequence
Mothur>pre.cluster(fasta= stability.trim.contigs.good.unique.good.filter.unique.fasta, count= stability.trim.contigs.good.unique.good.filter.count_table, diffs=2)

#Read fasta and counts to chimera sequences and then remove chimera 

Mothur>chimer.vsearch(fasta= stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count= stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

Mothur>remove.seqs( fasta= stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos= ab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

#Assign taxonomy
Mothur>classify.seqs(fasta= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=transient16_022016.pds.align, taxonomy= transient16_022016.pds.tax)

#Remove  Chloroplast, Mitochondria, Unknown, Archaea., and Eukaryotes 
Mothur>remove.lineage(fasta= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon= Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

#Calculate uncorrected pairwise distances between aligned DNA sequences. By default gaps are penalized; cutoff value indicate s that distances larger that 0.03 (>97%) won’t be saved
Mothur>dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)

#Assign sequences to OTUs using default opticlsut method
Mothur>cluster(column= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count= stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

#Determine how many sequences are in each OTU at the 0.03 cutoff level.
Mothur>make.shared(list= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count= stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

#Determine taxonomy for all OTUs
Mothur>classify.otus(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy label=0.03)

#Einson, J. E. et al., 2018
#combine paired end reads together
Mothur>make.contigs(file=estability.files)

#summary statistics to get descriptive stats of contigs 
Mothur>summary.seqs(fasta=estability.file.trim.contigs.fasta)

#filter sequences
##removesequences with ambig bases, shorter than 2.5%, and longer than 97.5%
Mothur>screen.seqs(fasta= estability.file.trim.contigs.fasta, group=estability.contigs.groups, summary= estability.trim.contigs.summary, minlength=440 maxlength=492, maxambig=0)

#remove identical (grouped) sequences; representative sequences were picked and stored in fasta , the other sequences are saved as sequence names to reduce the computational power/ work needed.

Mothur>unique.seqs(fasta= estability.trim.contigs.good.fasta)

#count table created of the unique sequences
Mothur>count_seqs(name= estability.trim.contigs.good.names, group= estability.contigs.good.groups)

#align the sequences to the silva reference database
Mothur>align.seqs(fasta= estability.trim.contigs.good.unique.fasta, reference=silva.bacteria.fasta)

#summary statistics to get descriptive stats of alignment
Mothur>summary.seqs(fasta= estability.trim.contigs.good.unique.align, count== estability.trim.contigs.good.count_table)

#Remove sequences that were before and after where the sequences aligned for the alignment site  alignment site for V3-V4 region of 16s
Mothur>screen.seqs(fasta= estability.trim.contigs.good.unique.align, count=estability.trim.contigs.good.count_table, summary= estability.trim.contigs.good.unique.summary, start=6388, end=25316, maxahomop=8)

#remove overhangs, remove alignment characters that only have “-“, using vertical=T and remove the sequences that contain “.”
Mothur>filter.seqs(fasta= estabilityility.trim.contigs.good.unique.good.align, vertical=T, trump=.)

#run unique seqs again to check for new reducndant sequences
Mothur>unique.seqs(fasta= estability.trim.contigs.good.unique.good.filter.fasta, count= estability.trim.contigs.good.good.filter. count_table)

#The sequences then need to be de-noised with a 2 threshold of mismatches in the sequence
Mothur>pre.cluster(fasta= estability.trim.contigs.good.unique.good.filter.unique.fasta, count= estability.trim.contigs.good.unique.good.filter.count_table, diffs=2)

#Read fasta and counts to chimera sequences and then remove chimera 

Mothur>chimer.vsearch(fasta= estability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count= estability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

Mothur>remove.seqs( fasta= estability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos= ab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

#Assign taxonomy
Mothur>classify.seqs(fasta= estability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= estability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=transient16_022016.pds.align, taxonomy= transient16_022016.pds.tax)

#Remove  Chloroplast, Mitochondria, Unknown, Archaea., and Eukaryotes 
Mothur>remove.lineage(fasta= estability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= estability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy= estability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon= Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

#Calculate uncorrected pairwise distances between aligned DNA sequences. By default gaps are penalized; cutoff value indicate s that distances larger that 0.03 (>97%) won’t be saved
Mothur>dist.seqs(fasta=estability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)

#Assign sequences to OTUs using default opticlsut method
Mothur>cluster(column= estability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count= estability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

#Determine how many sequences are in each OTU at the 0.03 cutoff level.
Mothur>make.shared(list= estability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count= estability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

#Determine taxonomy for all OTUs
Mothur>classify.otus(list=estability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=estability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy= estability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy label=0.03)


#Falardeau, J., et al., 2019

Mothur> make.file(stab.file, typ=fastq, prefix=stab)

#combine paired end reads together
Mothur>make.contigs(file=stab.files)

#summary statistics to get descriptive stats of contigs 
Mothur>summary.seqs(fasta=stab.file.trim.contigs.fasta)

#filter sequences
##removesequences with ambig bases, shorter than 2.5%, and longer than 97.5%
Mothur>screen.seqs(fasta= stab.file.trim.contigs.fasta, group=stab.contigs.groups, summary= stab.file.trim.contigs.summary, minlength=185, maxlength=210, maxambig=0)

#remove identical (grouped) sequences; representative sequences were picked and stored in fasta file, the other sequences are saved as sequence names to reduce the computational power/ work needed.

Mothur>unique.seqs(fasta= stab.trim.contigs.good.fasta)

#count table created of the unique sequences
Mothur>count_seqs(name= stab.trim.contigs.good.names, group= stability.contigs.good.groups)

#align the sequences to the silva reference database
Mothur>align.seqs(fasta= stab.trim.contigs.good.unique.fasta, reference=silva.bacteria.fasta)

#summary statistics to get descriptive stats of alignment
Mothur>summary.seqs(fasta= stab.trim.contigs.good.unique.align, count== stab.trim.contigs.good.count_table)

#Remove sequences that were before and after where the sequences aligned for the alignment site for V3 region of 16s
Mothur>screen.seqs(fasta= stab.trim.contigs.good.unique.align, count=stab.trim.contigs.good.count_table, summary= stab.trim.contigs.good.unique.summary, start=6388, end=13862, maxahomop=8)

#remove overhangs, remove alignment characters that only have “-“, using vertical=T and remove the sequences that contain “.”
Mothur>filter.seqs(fasta= stab.file.trim.contigs.good.unique.good.align, vertical=T, trump=.)

#run unique seqs again to check for new reducndant sequences
Mothur>unique.seqs(fasta= stab.trim.contigs.good.unique.good.filter.fasta, count= stab.trim.contigs.good.good.filter.count_table)

#The sequences then need to be de-noised with a 2 threshold of mismatches in the sequence
Mothur>pre.cluster(fasta= stab.trim.contigs.good.unique.good.filter.unique.fasta, count= stab.trim.contigs.good.unique.good.filter.count_table, diffs=2)

#Read fasta and countfiles to chimera sequences and then remove chimera 

Mothur>chimer.vsearch(fasta= stab.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count= stab.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

Mothur>remove.seqs( fasta= stab.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos= stab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

#Assign taxonomy
Mothur>classify.seqs(fasta= stab.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= stab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=transient16_022016.pds.align, taxonomy= transient16_022016.pds.tax)

#Remove  Chloroplast, Mitochondria, Unknown, Archaea., and Eukaryotes 
Mothur>remove.lineage(fasta= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count= stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy= stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon= Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

#Calculate uncorrected pairwise distances between aligned DNA sequences. By default gaps are penalized; cutoff value indicate s that distances larger that 0.03 (>97%) won’t be saved
Mothur>dist.seqs(fasta=stab.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)

#Assign sequences to OTUs using default opticlsut method
Mothur>cluster(column= stab.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count= stab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

#Determine how many sequences are in each OTU at the 0.03 cutoff level.
Mothur>make.shared(list= stab.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count= stab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

#Determine taxonomy for all OTUs
Mothur>classify.otus(list=stab.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stab.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy= stab.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy label=0.03)

