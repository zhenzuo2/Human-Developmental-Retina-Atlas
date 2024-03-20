cd /storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set;

bedtools intersect -a all_peaks.bed -b H3K4me2.bed H3K27ac.bed -f 0.20 -wa -wb -names H3K4me2 H3K27ac > compare_histone_all.bed;
bedtools intersect -a p2g_unique.bed -b H3K4me2.bed H3K27ac.bed -f 0.20 -wa -wb -names H3K4me2 H3K27ac > compare_histone_p2g.bed;




