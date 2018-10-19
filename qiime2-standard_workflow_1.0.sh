#!/bin/bash
#  This is a wrapper script for automation of QIIME2 workflow
#
#  Version 1.0 (October 18, 2018)
#
#  Copyright (c) 2018-2019 Deni Ribicic
#
#  This script is provided -as-is-, without any express or implied
#  warranty. In no event will the authors be held liable for any damages
#  arising from the use of this software.


#  load QIIME2
echo "activating QIIME2..."
source activate qiime2-2018.4
echo "qiime2-2018.4 active"


#  Get input files and locations
echo "looking for input files ..."
echo ""
if [ -f $1 ]; then
echo "using $1 file to import reads ..."
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path $1 --output-path 2_demux-paired-end.qza --source-format PairedEndFastqManifestPhred33
#  if no $1 has been found use directory 
elif [ -d $1 ]; then
echo "$1 is not a file"
echo "importing $1 as a folder ..."
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path $1 --source-format CasavaOneEightSingleLanePerSampleDirFmt --output-path 2_demux-paired-end.qza
fi

if [ -f "2_demux-paired-end.qza" ]; then
echo "2_demux-paired-end.qza generated"
echo ""
echo "preparing visualization file ..."
qiime demux summarize --i-data 2_demux-paired-end.qza --o-visualization 2_demux-paired-end.qzv
else
echo "2_demux-paired-end.qza FAILED TO GENERATE!!!"
fi

#  run QC, chimera check and OTU classification through DADA2
if [ -f "2_demux-paired-end.qza" ]; then
echo "running DADA2 step ... relax, this may take a while ..."
qiime dada2 denoise-paired --i-demultiplexed-seqs 2_demux-paired-end.qza --p-trim-left-f 0 --p-trunc-len-f 0 --p-trim-left-r 0 --p-trunc-len-r 0  --p-trunc-q $2 --o-representative-sequences 3_rep-seqs.qza --o-table 3_table.qza --p-n-threads $3 --output-dir dada2
fi

if [ -f "3_table.qza" ]; then
echo "preparing visualization file ..."
echo ""
qiime feature-table summarize --i-table 3_table.qza --o-visualization 3_table.qzv --m-sample-metadata-file $4
else
echo "3_table.qza FAILED TO GENERATE!!!"
fi

if [ -f "3_rep-seqs.qza" ]; then
echo "preparing visualization file ..."
qiime feature-table tabulate-seqs --i-data 3_rep-seqs.qza --o-visualization 3_rep-seqs.qzv
#  align  sequences
echo ""
echo "aligning sequences ..."
qiime alignment mafft --i-sequences 3_rep-seqs.qza --o-alignment 4_aligned-rep-seqs.qza
else
echo "3_rep-seqs.qza FAILED TO GENERATE!!!"
fi

if [ -f "4_aligned-rep-seqs.qza" ]; then
echo "sequences aligned :)"
echo "masking alignment to remove variable positions ..."
qiime alignment mask --i-alignment 4_aligned-rep-seqs.qza --o-masked-alignment 5_masked-aligned-rep-seqs.qza
else
echo "4_aligned-rep-seqs.qza FAILED TO GENERATE!!!"
fi

if [ -f "5_masked-aligned-rep-seqs.qza" ]; then
echo "variable positions masked or removed :)"
echo "generating phylogenetic tree ..."
qiime phylogeny fasttree --i-alignment 5_masked-aligned-rep-seqs.qza --o-tree 6_unrooted-tree.qza
else
echo "5_masked-aligned-rep-seqs.qza FAILED TO GENERATE!!!"
fi

if [ -f "6_unrooted-tree.qza" ]; then
echo "unrooted phylogenetic tree generated :)"
echo "applying midpoint rooting ..."
qiime phylogeny midpoint-root --i-tree 6_unrooted-tree.qza --o-rooted-tree 7_rooted-tree.qza
else
echo "6_unrooted-tree.qza FAILED TO GENERATE!!!"
fi

#  assigning taxonomy
echo "asigning taxonomy ... relax this may take a while :)"
qiime feature-classifier classify-sklearn --i-classifier $5 --i-reads 3_rep-seqs.qza --o-classification 9_taxonomy.qza
if [ -f "9_taxonomy.qza" ]; then
echo "taxonomy assigned :)"
else
echo "9_taxonomy.qza FAILED TO GENERATE!!!"
fi

#  converting all necessary files for Phyloseq input
echo "export 3_table.qza to biom file (feature-table.biom)"
qiime tools export 3_table.qza --output-dir exported
if [ -f "exported/feature-table.biom" ]; then
echo "biom file successfully exported :)"
else
echo "biom file FAILED TO EXPORT!!!"
fi

echo "export 9_taxonomy.qza to text file (taxonomy.tsv)"
qiime tools export 9_taxonomy.qza --output-dir exported
if [ -f "exported/taxonomy.tsv" ]; then
echo "taxonomy file successfully exported :)"
else
echo "taxonomy file FAILED TO EXPORT!!!"
fi

echo "export 7_rooted-tree.qza to newick file format (tree.nwk)"
qiime tools export 7_rooted-tree.qza --output-dir exported
if [ -f "exported/tree.nwk" ]; then
echo "tree file successfully exported :)"
else
echo "tree file FAILED TO EXPORT!!!"
fi
echo ""
echo "export 3_rep-seqs.qza to fasta file (dna-sequences.fasta)"
qiime tools export 3_rep-seqs.qza --output-dir exported
if [ -f "exported/dna-sequences.fasta" ]; then
echo ""
echo "fasta file successfully exported :)"
else
echo "fasta file FAILED TO EXPORT!!!"
fi

echo ""
echo ""
echo ""
echo "wrapper END"
