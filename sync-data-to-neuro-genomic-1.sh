#!/bin/bash

#sudo mount /mnt/neurogen-rw/egfr-viii-determination
rsync -avP --exclude ".*" --no-links 'output/' 'gsamadmin@neuro-genomic-1.erasmusmc.nl:/volume1/gsam/output/'
rsync -avP --exclude ".*" --exclude "__pycache__"  --no-links 'scripts/' 'gsamadmin@neuro-genomic-1.erasmusmc.nl:/volume1/gsam/scripts/'
rsync -avP --exclude ".*" --exclude "__pycache__"  --no-links 'documents/' 'gsamadmin@neuro-genomic-1.erasmusmc.nl:/volume1/gsam/documents/'

rsync -avP --exclude ".*" --exclude "*.fa" --exclude "SAindex" --exclude "__pycache__" --exclude "SA" --exclude "Genome" --exclude "rsem-hg19/" --no-links 'ref/' 'gsamadmin@neuro-genomic-1.erasmusmc.nl:/volume1/gsam/ref/'

# use -n for dry run


