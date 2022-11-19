#!/bin/bash

#sudo mount /mnt/neurogen-rw/egfr-viii-determination
#rsync -avP --exclude ".*" --no-links 'output/' 'neuro-genomic-1.erasmusmc.nl:/volume1/gsam/output/'
rsync -avP --exclude ".*" --exclude "output/alignments" --exclude "alignments/" --exclude "tables/rseqc/" --exclude "rseqc/" --no-links 'output/' 'neuro-genomic-1.erasmusmc.nl:/volume1/gsam/output/'
rsync -avP --exclude ".*" --exclude "__pycache__"  --no-links 'scripts/' 'neuro-genomic-1.erasmusmc.nl:/volume1/gsam/scripts/'
rsync -avP --exclude ".*" --exclude "__pycache__"  --no-links 'documents/' 'neuro-genomic-1.erasmusmc.nl:/volume1/gsam/documents/'

rsync -avP --exclude ".*" --exclude "*.fa" --exclude "SAindex" --exclude "__pycache__" --exclude "SA" --exclude "Genome" --exclude "rsem-hg19/" --no-links 'ref/' 'neuro-genomic-1.erasmusmc.nl:/volume1/gsam/ref/'


rsync -avP --exclude ".*" --exclude "__pycache__" --exclude "tin.py" --exclude "*.bam" --exclude "*.bai" --no-links 'output/tables/' 'neuro-genomic-1.erasmusmc.nl:/volume1/gsam/output/tables/'

# use -n for dry run


