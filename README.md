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