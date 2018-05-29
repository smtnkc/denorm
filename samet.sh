#!/bin/bash

for f in GPR/*.gpr
do
	awk 'BEGIN {FS="\t"} {OFS=";"} NR > 31 {print $1,$2,$3,$4,$5,$9}' < "$f" > "CSV/$(basename ${f%.*}).csv"
	echo "$f is done!"
done

f="GSM577967_CTRL_05.csv"
sed 's/\"//g' < "CSV/$f" > "CSV/v2.csv"
awk 'BEGIN {FS=";"} {OFS=";"} {print $1,$2,$3,"\""$4"\"","\""$5"\"",$6}' < "CSV/v2.csv" > "CSV/v3.csv"

var="\"Block\";\"Column\";\"Row\";\"Name\";\"ID\";\"F635 Median\""
sed -i "1s/.*/$var/" "CSV/v3.csv"
rm "CSV/v2.csv"
mv "CSV/v3.csv" "CSV/$(basename ${f%.*}).csv"
