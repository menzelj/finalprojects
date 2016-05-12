#! /usr/bin/env bash


####
# Investigating if copy number variations (CNVs) are
# more frequently associated with early or late replicating regions.
#
####

# library(bedtools)
# moldule load bedops 

###
# the following uses bedops tool called closest feature
# will measure distance to closest replication region. 
#
# this will add to each entry in our sorted outcnv.bed file the position
# and distance to the nearest region in the early replication
# file outrepstarts.bed, negative distances are downstream positive upstream and zeros are overlapping  

# sort bed files 
sort-bed cnv.bed > sortedcnv.bed
sort-bed repstarts.bed > sortedrepstarts.bed
sort-bed repvalleys.bed > sortedrepvalleys.bed


# closest early replication signal
closest-features --dist sortedcnv.bed sortedrepstarts.bed > early.bed

# replaces '|' with 'tabs' then prints cnv regions and their distance to nearest replication signal
tr '\|' '\t' < early.bed | awk ' {print $1, $2, $3, $4, $19, $29} ' > earlydist.bed 

# closest late replication signal 
closest-features --dist sortedcnv.bed sortedrepvalleys.bed > late.bed
tr '\|' '\t' < late.bed | awk ' {print $1, $2, $3, $4, $19, $29} ' > latedist.bed

# finds closest feature and makes new file for distances only 
cat earlydist.bed | awk '{ if ($5*-1 < $6) print $5*-1; else print $6; }' > dist_to_earlyreps.bed
cat latedist.bed | awk '{ if ($5*-1 < $6) print $5*-1; else print $6 }' > dist_to_latereps.bed

cat sortedcnv.bed | awk '{print $1, $2, $3, $4}' > cnvs.bed

# concatenates early and late distances into one file and ranks each CNV as early or late replicating or overlaping with both
paste cnvs.bed dist_to_earlyreps.bed dist_to_latereps.bed > cnvearly_late.bed
cat cnvearly_late.bed | awk '{ if ($5 < $6) print $0, "early"; if ( $5 == $6) print $0, "overlap"; if ($5 > $6) print $0, "late";}' > cnv_replication_timing.bed

# sort files into CNV type 
cat cnv_replication_timing.bed | grep -v "normal" > all.bed
cat cnv_replication_timing.bed | grep "homo.del" > homo.bed
cat cnv_replication_timing.bed | grep "het.del" > het.bed
cat cnv_replication_timing.bed | grep "amp" > amp.bed


# counts how many cnvs are early and late for each type

cat all.bed | grep "early" | wc -l | awk '{print $1 " total CNVs replicate early"}' >> counts.txt
cat all.bed | grep "late" | wc -l | awk '{print $1 " total CNVs replicate late"}' >> counts.txt
cat all.bed | grep "overlap" | wc -l | awk '{print $1 " total CNVs overlap with both early and late replication regions"}' >> counts.txt
 
cat homo.bed | grep "early" | wc -l | awk '{print $1 " homozygous deletion CNVs replicated early"}' >> counts.txt
cat homo.bed | grep "late" | wc -l | awk '{print $1 " homozygous deletion CNVs replicated late"}' >> counts.txt
cat homo.bed | grep "overlap" | wc -l | awk '{print $1 " homozygous deletion CNVs overlap with both early and late replication regions"}' >> counts.txt

cat het.bed | grep "early" | wc -l | awk '{print $1 " heterozygous deletion CNVs replicated early"}' >> counts.txt
cat het.bed | grep "late" | wc -l | awk '{print $1 " heterozygous deletion CNVs replicated late"}' >> counts.txt
cat het.bed | grep "overlap" | wc -l | awk '{print $1 " heterozygous deletion CNVs overlap with both early and late replication regions"}' >> counts.txt

cat amp.bed | grep "early" | wc -l | awk '{print $1 " amplification CNVs replicated early"}' >> counts.txt
cat amp.bed | grep "late" | wc -l | awk '{print $1 " amplification CNVs replicated late"}' >> counts.txt
cat amp.bed | grep "overlap" | wc -l | awk '{print $1 " amplification CNVs overlap with both early and late replication regions"}' >> counts.txt




