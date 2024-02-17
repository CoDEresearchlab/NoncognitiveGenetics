# format GWAS output

## concatenate chunks
for i in C_chol NC_chol

do
{
awk '
    FNR==1 && NR!=1 { while (/^SNP/) getline; }
    1 {print}
' ../out/$i.chunk.*.txt > ../out/$i.sumstat.txt
}&
done
wait


## format sumstats
echo "SNP CHR BP MAF A1 A2 est SE Z_Estimate Pval_Estimate chisq chisq_df chisq_pval AIC" > header

for i in C_chol.sumstat NC_chol.sumstat
do 
awk '/^rs/' $i.txt >  $i.tmp #exclude all lines that do not begin with rs (SNP names)
awk '{print $1, $2, $3, $4, $5, $6, $11, $12, $13, $14, $15, $16, $17, $18}'  $i.tmp >  $i.tmp2
cat header $i.tmp2  > $i.format
done &
wait

# remove redundant files
rm *.tmp  *.tmp2

# compress 
gzip *.format &
