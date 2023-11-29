Possible Pipeline
### 1. Prepare list of microbes
Possible papers to look at:
   - Table 8 provide ~70 species of microbes through blast [Lower respiratory tract microbiome composition and community interactions in smokers | Microbiology Society (microbiologyresearch.org)]
   -  Table 3 abundance table in genus level [Analysis of the Lung Microbiome in the “Healthy” Smoker and in COPD | PLOS ONE]
### 2. Probe design with PaintSHOP
 Requirements: A GTF file and FNA/FASTA file of the microbes.
 As there is currently no module supports download GTF file, we might download it manually.
 The scripts to run PaintSHOP pipeline on multiple species (loop through species) is in scripts folder, either test_loop.sh or something similar.
### 3. Filter the probes and select final probes
 Probe filtering
 - Could be done by using bowtie2 tools, to align the designed probes to references genome of other species.
 - If the probes could align on that species, that probe is not a good probe as it could cause false positive.
 - The scripts will be updated *soon*
 Select final probes
 - After filter the probes, probably there are still a lot of optional probes (not tested yet), for each species we should select ??? number of probes
### 4. Adding "arms" to the probes
Designing of the "arms" is not well understood just yet, but the "arms" should also not match with the genome of any species also
### 5. Construct code table
Will have to look at Arcadia research [Google Colab](https://colab.research.google.com/gist/jasegehring/87feac0f48971598cf63e9ef50626b08/encodingprobesfilteringandqc_20220512.ipynb)
