#!/bin/bash

#sudo mount /mnt/neurogen-rw/egfr-viii-determination
rsync -avP --exclude ".*" --no-links 'output/' 'gsamadmin@neuro-genomic-1.erasmusmc.nl:/volume1/gsam/output/'
rsync -avP --exclude ".*" --exclude "__pycache__"  --no-links 'scripts/' 'gsamadmin@neuro-genomic-1.erasmusmc.nl:/volume1/gsam/scripts/'
rsync -avP --exclude ".*" --exclude "__pycache__"  --no-links 'docs/' 'gsamadmin@neuro-genomic-1.erasmusmc.nl:/volume1/gsam/docs/'



