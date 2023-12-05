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
 
