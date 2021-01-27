 SKraken: Fast and Sensitive Classification of Short Metagenomic Reads based on filtering uninformative k-mers

The study of microbial communities is an emerging field that is revolutionizing many disciplines from ecology to medicine. The major problem when analyzing a metagenomic sample is to taxonomic annotate its reads in order to identify the \begin{flushright} \end{flushright} species in the sample and their relative abundance. Many tools have been developed in the recent years, however the performance in terms of precision and speed are not always adequate for these very large datasets. In this work we present SKraken an efficient approach to accurately classify metagenomic reads against a set of reference genomes, e.g. the NCBI/RefSeq database. SKraken is based on $k$-mers statistics combined with the taxonomic tree. Given a set of target genomes
SKraken is able to detect the most representative $k$-mers for each species, filtering out uninformative $k$-mers. The classification performance on several synthetic and real metagenomics datasets shows that SKraken achieves in most cases the best performances in terms of precision and recall w.r.t. Kraken. In particular, at species level classification, the estimation of the abundance ratios improves by 6\% and the precision by $8\%$. This behavior is confirmed also on a real stool metagenomic sample where SKraken is able to detect species with high precision. Because of the efficient filtering of uninformative $k$-mers, SKraken requires less RAM than Kraken.

Licence

The software is freely available for academic use.

For questions about the tool, please contact Matteo Comin.

Reference

Please cite the following paper:

J. Qian, D. Marchiori, M. Comin
"Fast and Sensitive Classification of Short Metagenomic Reads with SKraken"
Communications in Computer and Information Science 2018, vol 881, pp. 212-226.

D. Marchiori, M. Comin
"SKraken: Fast and Sensitive Classification of Short Metagenomic Reads based on Filtering Uninformative k-mers".
In Proceedings of the 10th International Conference on Bioinformatics Models, Methods and Algorithms (Bioinformatics 2017), pp. 59-67.




--------------------------


# **SKRAKEN** #

# **Installation** #

## **Download Skraken:** ##

### Download Skraken at: ###
* [Skraken](https://bitbucket.org/marchiori_dev/skraken/downloads/skraken.tar.gz)

## **Compilation:** ##

Open terminal and go to skraken-master/ and then:



```
#!

./install_skraken.sh $SKRAKEN_DIR
```




(Replace "$SKRAKEN_DIR" above with the directory where you want to install Skraken's programs/directories.)

## **Standard Skraken Database:** ##

To create the standard Skraken database, you can use the following command:



```
#!


skraken-build --standard --db $DBNAME
```




(Replace "$DBNAME" above with your preferred database name/location. Please note that the database will use approximately 160 GB of disk space during creation.)
This will download NCBI taxonomic information, as well as the complete genomes in RefSeq for the bacterial, archaeal, and viral domains.

If you have multiple processing cores, you can run this process with multiple threads, e.g.:


```
#!


skraken-build --standard --threads 16 --db $DBNAME
```






## **Custom Databases** ##

Skraken also allows creation of customized databases.
### To build a custom database: ###
1. Install a taxonomy.
2. Install a genomic library.
3. Build the database with :

```
#!


skraken-build --build --db $DBNAME
```


### Parameter added: ###


```
#!

--score-filter NUM
```


Remove k-mers stored in DB with a SCORE less than NUM. 
Paramenter NUM is a float representing a percentage in [0,100] range.
If at least one k-mer is removed this mode creates four DB file:

* **database.kdb**
Filtered Skraken DB.
* **database.idx**
Filtered Skraken DB Index
* **database_full.kdb**
 Full Kraken DB
* **database_full.idx**
Full Kraken DB Index

If parameter is not set, output DB is equal to Kraken.
A full list of options for skraken-build can be obtained using skraken-build --help.

## **Classification** ##

To classify a set of sequences (reads), use the skraken command:


```
#!

skraken --db $DBNAME seqs.fa
```


## **References:** ##

* [Kraken webpage](http://ccb.jhu.edu/software/kraken/)
* [Kraken manual](http://ccb.jhu.edu/software/kraken/MANUAL.html)


## **Dataset:** ##

* [SRR1804065](https://bitbucket.org/marchiori_dev/skraken/downloads/SRR1804065_blast_filtered_paired.tar.gz)
