1. Put all code on the private github repo so we both can access it and easier to comment etc.





~~2. Read and present the paintshop workflow  in detail (underlying code software used, and why) at one of our meetings.~~




~~3. regarding workflow using paintshop followed by quality control :~~

~~Check limits of how you can blast sequences for similarity. This works perfect for a few but will of course have limitations for many.~~

~~You can download blast and do it offline. How does it scale can one do that for so many probes ?~~

~~https://www.youtube.com/watch?v=A-DPU5xxyHs~~

~~Here is a list https://www.ebi.ac.uk/Tools/msa/  of massive sequence aligners~~

~~Other well known software are MAFFT or clustalW. Perhaps they perform better.~~ 



4. Write an intro (background)  section to multiplexed FISH and show at a meeting.

This will be parts of the intro section in your thesis and it feels good to know you have already started writing. I didn’t do that so learn from my mistakes : )



~~5. Regarding the other approach aligning the genomes and find unique regions or deletions that can be used to discriminate species. This is still a valid approach and we similarly to "3".  We struggled to find performant alternatives. I recently found tools such as FastANI and the newer SKANI:~~
~~https://github.com/bluenote-1577/skani They seem be able to aligning many large genomes in short time.~~
~~Perhaps you could have a go on Skani !? try just a few and see.~~



6. make probes for the human genome. Iso flattened to start with.

Once done please compare probes with published literature for a for a few genes.   

Find a couple of mRNA commonly used in some paper that we can try. Likely Malat1 (nuclear localised) and few 

seqFISH+ : https://www.nature.com/articles/s41586-019-1049-y and other more recent from long cai´s lab. 
MERFISH :  https://www.pnas.org/doi/abs/10.1073/pnas.1912459116     dataset s1 s2. has .csv files. 

7. Look into how to add readout seqences (the extended arms) what is needed. PCR primer sequence for sure but also landing sites that fits  imaging probe landing extended region.
 These are important regions and well functional sequences are likely reused. We need to find those from litterature. detective work! : )
 Initially just find and decide what to use - later work: How can we more formally have software that add these properly.
 check what paintshop or other software has for capabilities. 
 First look what paintshop can do. Q: can they do more than binary (color/signal vs no color/signal) or several colors (no signal vs color1 vs  
 color2 ) ? if so how is the error correction done?
 how do we select. overhang sequences? are they suggestions??  otherwise see previously mentioned detective work to find these in litterature.
 (perhaps also look at arcadia´s python script for overhang.)

8. Add document listing all our microbes with name and Importantly! the reference where we found info that this microbe was found in the lung.
8.1 Add folder with all the  genome files listed in document (point 8). 
 
~~9. Collect and compose a list of all bacterial species (with references) that are found in the lung. we need this list as complete as possible to know what potential "double-hits" we have (ie FISH probes designed for one bacterial sequence that targets another nearby species).~~

~~10. Generate unique probes using simlar technique as of the 4 in trial run. Aim for 50-80 species.~~

~~11.  Find out the position (what mRNA it codes for ) for each unique probe. Once we know what kind protein coding mRNA are targeted can we group them (replication, virulence etc.) In doing that can we compare different species more refined has one more upregulated replication than the other species~~

~~12. (related to 7)  how do we know which probes are targeting which mRNA given the output of paintshop ?~~
     In other words: Given the output of paintshop probe generation for a single species. We are expecting several probes per mRNA.
     is that so? can we double check that. and more generally what mRNAs are targeted with all probes?
      
13. Add more info in the introduction for (F)ISH and ISS. Next add a section for the overhang design   explain where the error robustness in MERFISH comes from. 

14. Try to find a general probe set with existing targets in all microbes. 

15.  Visualize how the number of unique probes change with chainging number of bases "probe matching criteria".

16.  gene description and gene type for all unique probes. We have redundancy for some of the species with 1000s of unique probes how do we choose the 100 best? helps by knowing what protein coding sequence that are targeted. etc. of course useful for other studies later on.  


17.   
Either generate new read out probe landing sequences from Random sequence  


OR check off-targets and change bps in existing read out probe landing sequences   

then of course blast No Cs etc. no “CCCTAA” or “TTAGGG” sequences to exclude the potential binding to telomeres etc. (see below for ref.).





https://ars.els-cdn.com/content/image/1-s2.0-S0006349517303430-mmc1.pdf

The mm10 mouse genomic sequence (UCSC Genome Bioinformatics) was used to design subtelomere oligonucleotide probe pools in this study. To selectively label subtelomeric genomic regions, 100 kb regions at the end of each chromosome were selected (Table S1). Across those regions, a set of non-overlapping 35-nt probes were designed which suffice several constraints including 40-60% GC content, no more than 5 contiguous identical nucleotides, no “CCCTAA” or “TTAGGG” sequences to exclude the potential binding to telomeres, and at least 2-nt spaces between adjacent probe

![image](https://github.com/npxhuy/thesis/assets/12096956/cbb284e6-7f2b-4458-b22e-dcefc292c29e)

18.
Please check that all probes for the microbes do not have off-targets in human. this will be the first selection criteria when we have many unique probes for a single species.
 
