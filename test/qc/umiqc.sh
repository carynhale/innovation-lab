#!/usr/bin/env sh

echo -e 'FamilySize\tFamilyType\tFrequency\tSAMPLE' > family-sizes.txt
echo -e 'SAMPLE\tType\tCount' > family-types-A.txt
echo -e 'SAMPLE\tType\tCount' > family-types-B.txt
ln -s `readlink -e $1` A-positions.txt
ln -s `readlink -e $2` B-positions.txt

# family sizes
grep "@Marianas" collapsed_R1_.fastq | awk -v sample=$3 'BEGIN{FS=OFS="\t"; while(getline < "A-positions.txt"){pos=$1"\t"$2; positions[pos]=pos} FS=":"}{if(positions[$3"\t"($4+20)]==null) next; if($5!=0 && $6!=0 && $9!=0 && $10!=0) type="Duplex"; else if($5+$6<=2 || $9+$10<=2) type="SingletonOrSubSimplex"; else type="Simplex"; key=$5+$6"\t"type; freqs[key]++; allKey=$5+$6"\tAll"; allFreqs[allKey]++;}END{ for(key in freqs) print key, freqs[key], sample; for(key in allFreqs) print key, allFreqs[key], sample}' >> family-sizes.txt

grep "@Marianas" collapsed_R1_.fastq | awk -v sample=$3 'BEGIN{FS=OFS="\t"; while(getline < "A-positions.txt"){pos=$1"\t"$2; positions[pos]=pos} FS=":"}{if(positions[$3"\t"($4+20)]==null) next; if($5!=0 && $6!=0 && $9!=0 && $10!=0) duplex++; else if($5+$6==1 || $9+$10==1) singletons++; else if($5+$6==2 || $9+$10==2) subSimplex++; else simplex++}END{print sample, "Singletons", singletons; print sample, "Sub-Simplex", subSimplex; print sample, "Simplex", simplex; print sample, "Duplex", duplex}' >> family-types-A.txt

grep "@Marianas" collapsed_R1_.fastq | awk -v sample=$3 'BEGIN{FS=OFS="\t"; while(getline < "B-positions.txt"){pos=$1"\t"$2; positions[pos]=pos} FS=":"}{if(positions[$3"\t"($4+20)]==null) next; if($5!=0 && $6!=0 && $9!=0 && $10!=0) duplex++; else if($5+$6==1 || $9+$10==1) singletons++; else if($5+$6==2 || $9+$10==2) subSimplex++; else simplex++}END{print sample, "Singletons", singletons; print sample, "Sub-Simplex", subSimplex; print sample, "Simplex", simplex; print sample, "Duplex", duplex}' >> family-types-B.txt
