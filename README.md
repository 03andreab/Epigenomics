Epigenomics Final Task
================
Andrea Bultó Sánchez
2024-02-14

- [Task 4](#task-4)
  - [EN‐TEx ATAC‐seq data: downstream
    analyses](#entex-atacseq-data-downstream-analyses)
    - [1. Move to folder ATAC-seq, and create folders to store bigBed
      data files and peaks analyses files. Make sure the files are
      organized in a consistent way as done for
      ChIP-seq.](#1-move-to-folder-atac-seq-and-create-folders-to-store-bigbed-data-files-and-peaks-analyses-files-make-sure-the-files-are-organized-in-a-consistent-way-as-done-for-chip-seq)
    - [2. Retrieve from a newly generated metadata file ATAC-seq peaks
      (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for
      stomach and sigmoid_colon for the same donor used in the previous
      sections. Make sure your md5sum values coincide with the ones
      provided by
      ENCODE.](#2-retrieve-from-a-newly-generated-metadata-file-atac-seq-peaks-bigbed-narrow-pseudoreplicated-peaks-assembly-grch38-for-stomach-and-sigmoid_colon-for-the-same-donor-used-in-the-previous-sections-make-sure-your-md5sum-values-coincide-with-the-ones-provided-by-encode)
    - [3. For each tissue, run an intersection analysis using
      BEDTools.](#3-for-each-tissue-run-an-intersection-analysis-using-bedtools)
- [Task 5](#task-5)
  - [Task 1](#task-1)
  - [Task 2](#task-2)
  - [Task 3](#task-3)
  - [Task 4](#task-4-1)
  - [Task 5](#task-5-1)
  - [Task 6](#task-6)
  - [Task 7](#task-7)

# Task 4

The following exercises come from
[here](https://github.com/bborsari/epigenomics_uvic/wiki/4.-EN%E2%80%90TEx-ATAC%E2%80%90seq-data:-downstream-analyses).

## EN‐TEx ATAC‐seq data: downstream analyses

### 1. Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files. Make sure the files are organized in a consistent way as done for ChIP-seq.

First we run the container inside our folder were we have already cloned
the [git repository](https://github.com/bborsari/epigenomics_uvic):

``` bash
docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
```

``` bash
cd ATAC-seq
mkdir analyses
mkdir data/bigBed.files data/bigWig.files
```

### 2. Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections. Make sure your md5sum values coincide with the ones provided by ENCODE.

We filter the data for only ATAC-seq and download those files
corresponding to donor ENCDO451RUA. The URL of the metadata file is
provided on the first line of the files.txt:

``` bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_title=ATAC-seq&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&type=Experiment&files.analyses.status=released&files.preferred_default=true"
```

We downloaded the metadata file.

Retreaving bigBed peak calling files from the metadeta file:

``` bash
grep -F "bigBed_narrowPeak" metadata.tsv|\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.txt

cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

ENCFF287UHP.bigBed ENCFF762IFP.bigBed

We make sure the md5sum values coincide with the ones provided by ENCODE

``` bash
# retrieve original MD5 hash from the metadata
../bin/selectRows.sh <(cut -f1 analyses/bigBed.*.ids.txt) metadata.tsv | cut -f1,46 > data/bigBed.files/md5sum.txt
# compute MD5 hash on the downloaded files
cat data/bigBed.files/md5sum.txt |\
while read filename original_md5sum; do 
  md5sum data/bigBed.files/"$filename".bigBed |\
  awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
done > tmp 
mv tmp data/bigBed.files/md5sum.txt
# make sure there are no files for which original and computed MD5 hashes differ
awk '$2!=$3' data/bigBed.files/md5sum.txt
```

No output means they coincide.

### 3. For each tissue, run an intersection analysis using BEDTools.

**1) the number of peaks that intersect promoter regions**

We download the promoter regions:

``` bash
wget -P ../ChIP-seq/annotation
 https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/gencode.v24.protein.coding.non.redundant.TSS.bed
```

Run bedtools intersect and retrieve the ATAC-seq peaks that intersect
with promoters:

``` bash
#filtered
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b ../ChIP-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -u | 
  sort -k 2,2 -k 3,3 | awk 'BEGIN{OFS="\t"} !seen[$1,$2,$3]++ {print $1, $2, $3, $4}' | sort -u > analyses/peaks.analysis/peaks."$tissue"_filtered.bed ; done  

#not filtered
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b ../ChIP-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -u | 
  sort -k 2,2 -k 3,3 | cut -f1,2,3,4 | sort -u > analyses/peaks.analysis/peaks."$tissue"_non.filtered.bed ; done 
```

We see sometimes the name of the peak is different but the coordinates
are the same, therefore, we will provide filtered and non-filtered peaks
for ATAC-seq with the promoter:

**Number of peaks that intersect promoter regions:** <br> Total: <br>
43711 filtered <br> 92620 not filtered <br>

Sigmoid colon *filtered* peaks: 21500 <br> Stomach *filtered* peaks:
22211 <br>

Sigmoid colon *non-filtered* peaks: 47871 <br> Stomach *non-filtered*
peaks: 44749 <br>

| Tissue        | Filtered | Non-filtered |
|---------------|----------|--------------|
| Sigmoid colon | 21500    | 47871        |
| Stomach       | 22211    | 44749        |
| Total         | 43711    | 92620        |

**2) the number of peaks that fall outside gene coordinates** <br>

Prepare a BED file with gene body coordinates of protein-coding genes:
*(we already did this in class, so we do not need to run the code
again)*

``` bash
awk '$3=="gene"' ../ChIP-seq/annotation/gencode.v24.primary_assembly.annotation.gtf |\
grep -F "protein_coding" |\
cut -d ";" -f1 |\
awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $10, 0, $7, $10}' |\
sed 's/\"//g' |\
awk 'BEGIN{FS=OFS="\t"}$1!="chrM"{$2=($2-1); print $0}' > ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed
```

ATAC-seq peaks that fall outside gene body:

``` bash
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed -v | 
  sort -k 2,2 -k 3,3 | awk 'BEGIN{OFS="\t"} !seen[$1,$2,$3]++ {print $1, $2, $3, $4}' | sort -u > analyses/peaks.analysis/peaks.body_"$tissue"_filtered.bed ; done

cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed -v | 
  sort -k 2,2 -k 3,3 | cut -f1,2,3,4 | sort -u > analyses/peaks.analysis/peaks.body_"$tissue"_non.filtered.bed ; done 
```

Sigmoid colon *filtered* peaks: 25635 <br> Stomach *filtered* peaks:
25665 <br>

Sigmoid colon *non-filtered* peaks: 37035 <br> Stomach *non-filtered*
peaks: 34537 <br>

**Number of peaks that fall outside of gene body coordinates:** <br>

| Tissue        | Filtered | Non-filtered |
|---------------|----------|--------------|
| Sigmoid colon | 25635    | 37035        |
| Stomach       | 25665    | 34537        |
| Total         | 60172    | 71572        |

# Task 5

The following exercises come from
[here](https://github.com/bborsari/epigenomics_uvic/wiki/5.-Distal-regulatory-activity).
\## Distal regulatory activity

### Task 1

Create a folder regulatory_elements inside epigenomics_uvic. This will
be the folder where you store all your subsequent results.

``` bash
mkdir ../regulatory_elements
cd ../regulatory_elements
mkdir data  peaks.analysis data/bigBed.files data/bed.files 
```

### Task 2

Distal regulatory regions are usually found to be flanked by both
H3K27ac and H3K4me1. From your starting catalogue of open regions in
each tissue, select those that overlap peaks of H3K27ac AND H3K4me1 in
the corresponding tissue. You will get a list of candidate distal
regulatory elements for each tissue. How many are they?

Get peaks of H3K27ac and H3K4me1 from the metadata stored in Chip-seq.

``` bash
for mark in H3K27ac H3K4me1; do grep -F $mark ../ChIP-seq/metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u;done > bigBed.peaks.ids.txt

cut -f1  bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done

cd regulatory_elements
# retrieve original MD5 hash from the metadata
../bin/selectRows.sh <(cut -f1 bigBed.*.ids.txt) ../ChIP-seq/metadata.tsv | cut -f1,46 > data/bigBed.files/md5sum.txt

# compute MD5 hash on the downloaded files
cat data/bigBed.files/md5sum.txt |\
while read filename original_md5sum; do 
  md5sum data/bigBed.files/"$filename".bigBed |\
  awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
done > tmp 
 mv tmp data/bigBed.files/md5sum.txt
# make sure there are no files for which original and computed MD5 hashes differ -- we dont expect an output
awk '$2!=$3' data/bigBed.files/md5sum.txt

cut -f1,2,3 bigBed.peaks.ids.txt |\
  while read filename tissue marker; \
  do bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$tissue"."$marker".bed;done
```

56661 sigmoid_colon.H3K27ac-human.bed <br> 97950
sigmoid_colon.H3K4me1-human.bed <br> 57121 stomach.H3K27ac-human.bed
<br> 68664 stomach.H3K4me1-human.bed <br>

**Number of peaks of each type of histone mark that fall outside the
gene body:** <br>

| Tissue        | H3K27ac | H3K4me1 |
|---------------|---------|---------|
| Sigmoid colon | 56661   | 97950   |
| Stomach       | 57121   | 68664   |

Now we intersect the ATAC peaks found outside the gene body with one
histone mark and we intersect the output with the other histone mark.
Per tissue.

``` bash
for tissue in stomach sigmoid_colon; do
  # Intersect ATAC outside peaks with each histone mark: 
    bedtools intersect -a ../ATAC-seq/analyses/peaks.analysis/peaks.body_"$tissue"_non.filtered.bed -b data/bed.files/"$tissue".H3K27ac-human.bed -u > peaks.analysis/ATAC_H3K27ac_overlap_"$tissue"_.bed \
    
    bedtools intersect -a peaks.analysis/ATAC_H3K27ac_overlap_"$tissue"_.bed -b data/bed.files/"$tissue".H3K4me1-human.bed -u > peaks.analysis/combined.peaks."$tissue".bed; done
```

14215 combined.peaks.sigmoid_colon.bed <br> 8022
combined.peaks.stomach.bed <br>

To check for non repeated coordinates:

``` bash
for file in peaks.analysis/combined.peaks.sigmoid_colon.bed peaks.analysis/combined.peaks.stomach.bed; \
  do sort -k 2,2 -k 3,3 "$file"| awk 'BEGIN{OFS="\t"} !seen[$1,$2,$3]++ {print $1, $2, $3, $4}' | sort -u | wc -l; done
```

8627 \# Filtered combined.peaks.sigmoid_colon.bed <br> 5148 \# Filtered
combined.peaks.stomach.bed <br>

**ATAC-seq, H3K27ac and H3K4me1 overlap peaks:** <br>

| Tissue        | Filtered | Non-filtered |
|---------------|----------|--------------|
| Sigmoid colon | 8627     | 14215        |
| Stomach       | 5148     | 8022         |
| Total         | 13775    | 22237        |

### Task 3

Focus on regulatory elements that are located on chromosome 1, and
generate a file regulatory.elements.starts.tsv that contains the name of
the regulatory region (i.e. the name of the original ATAC-seq peak) and
the start (5’) coordinate of the region.

regulatory.elements.starts.tsv:

``` bash
for tissue in sigmoid_colon stomach; do
  grep -w 'chr1' peaks.analysis/combined.peaks."$tissue".bed |sort -k 2,2 -k 3,3| awk '{print $4 "\t" $2}'  > regulatory.elements.starts_"$tissue".tsv;done
```

1521 regulatory.elements.starts_sigmoid_colon.tsv <br> 987
regulatory.elements.starts_stomach.tsv <br> 2508 total <br>

Filtered:

``` bash
for tissue in sigmoid_colon stomach; do
  grep -w 'chr1'  peaks.analysis/combined.peaks."$tissue".bed | sort -k 2,2 -k 3,3 | awk 'BEGIN{OFS="\t"} !seen[$1,$2,$3]++ {print $1, $2, $3, $4}' | awk '{print $4 "\t" $2}' | sort -u > regulatory.elements.starts_"$tissue"_filtered.tsv;done
```

939 regulatory.elements.starts_sigmoid_colon_filtered.tsv <br> 640
regulatory.elements.starts_stomach_filtered.tsv <br>

**Regulatory elements start position chr1:** <br>

| Tissue        | Filtered | Non-filtered |
|---------------|----------|--------------|
| Sigmoid colon | 939      | 1521         |
| Stomach       | 640      | 987          |
| Total         | 1579     | 2508         |

### Task 4

Focus on protein-coding genes located on chromosome 1. From the BED file
of gene body coordinates that you generated, prepare a tab-separated
file called gene.starts.tsv which will store the name of the gene in the
first column, and the start coordinate of the gene on the second column
(REMEMBER: for genes located on the minus strand, the start coordinate
will be at the 3’). Use the command below as a starting point:

``` bash
grep -w 'chr1' ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed | awk 'BEGIN{FS=OFS="\t"}{if ($6=="+"){start=$2} else {start=$3}; print $4, start}'> gene.starts.tsv
```

### Task 5

Download the python script inside the epigenomics_uvic/bin folder. Have
a look at the help page of this script to understand how it works:

``` bash
wget -P ../bin "https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/get.distance.py"
```

This script takes as input two distinct arguments: 1) –input corresponds
to the file gene.starts.tsv (i.e. the file you generated in Task \#4);
2) –start corresponds to the 5’ coordinate of a regulatory element.
Complete the python script so that for a given coordinate –start the
script returns the closest gene, the start of the gene and the distance
of the regulatory element.

``` python
#************

import sys
from optparse import OptionParser


#*****************
# OPTION PARSING *
#*****************

parser = OptionParser()
parser.add_option("-i", "--input", dest="input")
parser.add_option("-s", "--start", dest="start")
options, args = parser.parse_args()

open_input = open(options.input)
enhancer_start = int(options.start)


#********
# BEGIN *
#********

x=1000000 # set maximum distance to 1 Mb
selectedGene="" # initialize the gene as empty
selectedGeneStart=0 # initialize the start coordinate of the gene as empty

for line in open_input.readlines(): # for each line in the input file
    gene, y = line.strip().split('\t') # split the line into two columns based on a tab
    position= int(y)    # define a variable called position that correspond to the integer of the start of the gene
    absolute = abs(position - enhancer_start) # compute the absolute value of the difference between position and enhancer_start

    if absolute < x:    # if this absolute value is lower than x
        x = absolute # this value will now be your current x
        selectedGene = gene # save gene as selectedGene
        selectedGeneStart = position # save position as selectedGeneStart

print "\t".join([selectedGene, str(selectedGeneStart), str(x)])
```

### Task 6

For each regulatory element contained in the file
regulatory.elements.starts.tsv, retrieve the closest gene and the
distance to the closest gene using the python script you created above.
Use the command below as a starting point:

``` bash
for tissue in stomach sigmoid_colon; do
  cat regulatory.elements.starts_"$tissue".tsv |  while read element start; do 
   python ../bin/get.distance.py -i gene.starts.tsv --start "$start"; 
done > regulatoryElements.genes.distances_"$tissue".tsv; done
##### To retrieve non repeated coordinates:
for tissue in stomach sigmoid_colon; do
  cat regulatory.elements.starts_"$tissue"_filtered.tsv |  while read element start; do 
   python ../bin/get.distance.py -i gene.starts.tsv --start "$start"; 
done > regulatoryElements.genes.distances_"$tissue"_filtered.tsv; done
```

### Task 7

Use R to compute the mean and the median of the distances stored in
regulatoryElements.genes.distances.tsv.

``` bash
for tissue in stomach sigmoid_colon; do
  Rscript -e "x <- read.table('regulatoryElements.genes.distances_"$tissue".tsv', header = FALSE, sep='\t'); cat('${tissue} Mean:', mean(x[,3]), '\n'); cat('${tissue} Median:', median(x[,3]), '\n')";
done
```

stomach Mean: 45227.05 <br> stomach Median: 27735 <br> sigmoid_colon
Mean: 73635.89 <br> sigmoid_colon Median: 35802 <br>

| Tissue        | Mean     | Median |
|---------------|----------|--------|
| Sigmoid colon | 73635.89 | 35802  |
| Stomach       | 45227.05 | 27735  |

For filtered files:

``` bash
for tissue in stomach sigmoid_colon; do
  Rscript -e "x <- read.table('regulatoryElements.genes.distances_"$tissue"_filtered.tsv', header = FALSE, sep='\t'); cat('${tissue} Mean:', mean(x[,3]), '\n'); cat('${tissue} Median:', median(x[,3]), '\n')";
done
```

stomach Mean: 47013.77 <br> stomach Median: 27773.5 <br> sigmoid_colon
Mean: 73067.4 <br> sigmoid_colon Median: 36045 <br>
