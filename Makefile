
data:
	ln -s /mnt/neurogen-ro/gsam data

misc/h.all.v6.2.entrez.gmt:
	mkdir -p misc ; cd misc ; wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/6.2/h.all.v6.2.entrez.gmt

misc/msigdb.v6.2.entrez.gmt:
	mkdir -p misc; cd misc ; wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/6.2/msigdb.v6.2.entrez.gmt

