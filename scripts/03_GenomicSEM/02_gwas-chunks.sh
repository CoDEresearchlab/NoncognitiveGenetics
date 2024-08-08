#create chunks for GWAS

#capture header in varaible header
header=`head -1 ../out/formatted_sumstats_07102020.txt`. 

#split in numbered shunks starting name with DEP_SMK_CAN_ALC_
split -d -l 50000 ../out/formatted_sumstats_07102020.txt ../out/COG_nonCOG_exp_

#find chunks and add header 
find .|grep COG_nonCOG_exp_ | xargs sed -i "1s/^/$header\n/"

#filter out redundant header in chunk _00
sed -i '1d' COG_nonCOG_exp_00

#make a list of sumstats chunks for submission as array job
ls COG_nonCOG_exp_* > listGwasChunks

# N of chunks to be submitted as array job:
wc -l listGwasChunks 

#135 listGwasChunks