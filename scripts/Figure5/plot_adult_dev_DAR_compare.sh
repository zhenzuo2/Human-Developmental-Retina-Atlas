cd /storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set;
celltypes=("AC" "BC" "Cone" "HC" "RGC" "Rod" "PRPC" "NRPC" "MG")
for cy in "${celltypes[@]}"; do
    echo ${cy};
    bedtools intersect -a Dev_DAR_peak_${cy}.bed -b Adult_DAR_peak_AC.bed Adult_DAR_peak_BC.bed Adult_DAR_peak_Cone.bed Adult_DAR_peak_HC.bed Adult_DAR_peak_RGC.bed Adult_DAR_peak_Rod.bed Adult_DAR_peak_MG.bed -f 0.20 -wa -wb -names AC BC Cone HC RGC Rod MG > ${cy}_Dev_Adult_Compare.bed;
done
