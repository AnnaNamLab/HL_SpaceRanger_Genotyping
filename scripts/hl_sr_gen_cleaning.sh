#!/bin/bash
SAMPLES=('HL8' 'HL8R' 'HL21' 'HL24')


############################################ STEP 1 ############################################################
# filter Space Ranger BAM for records that contribute to UMI counting ("xf:i:25" tag)

for sample in "${SAMPLES[@]}"; do

        input_bam="${sample}_sr_gen_bam.bam"
        output_sam="${sample}_sr_gen_umis.sam"
        samtools view "$input_bam" | grep -wF "xf:i:25" > "$output_sam"
        echo "Processed $sample"
done


############################################ STEP 2 ############################################################
# For each target, extract records labeled with WT/MUT probes
# 'probe_tags.txt' is a space-separated text file with one row per target containing the target name and probe tags
# example entry: 'ABCF1 DEPRECATED_ENSG00000204574_MUT DEPRECATED_ENSG00000204574_WT' 

while read -r target mut_id wt_id; do

    for sample in "${SAMPLES[@]}"; do

        input_sam="${sample}_sr_gen_umis.sam"
        mut_sam="${sample}/mut_umis/${target}_umis.sam"
        wt_sam="${sample}/wt_umis/${target}_umis.sam"

        grep -wF "$mut_id" "$input_sam" > "$mut_sam"
        grep -wF "$wt_id" "$input_sam" > "$wt_sam"

        mut_seqs="${sample}/mut_umis/seqs/${target}_seqs.txt"
        wt_seqs="${sample}/wt_umis/seqs/${target}_seqs.txt"

        awk -F'\t' '{print $10}' "$mut_sam" > "$mut_seqs"
        awk -F'\t' '{print $10}' "$wt_sam" > "$wt_seqs"

    done

done < probe_tags.txt


############################################ STEP 3 ############################################################
# determine orientation of R2 read (forward vs reverse complement) for each target and genotype (WT/MUT)
# extract any records that have a discrepant probe label and allele (i.e. MUT probe label but no MUT allele)
# 'std_rev.txt' is a space-separated text file with probe sequence + allele + position information
# example entry: 'ABCF1 TCTCTCCGTCCCCGATCCACTCGG CCGAGTGGATCGGGGACGGAGAGA C T 26 G A 25'

for sample in "${SAMPLES[@]}"; do

    while read -r target forward reverse_comp for_wt for_mut for_pos rc_wt rc_mut rc_pos; do

        mut_seqs="${sample}/mut_umis/seqs/${target}_seqs.txt"
        wt_seqs="${sample}/wt_umis/seqs/${target}_seqs.txt"
        mut_sam="${sample}/mut_umis/${target}_umis.sam"
        wt_sam="${sample}/wt_umis/${target}_umis.sam"

        mut_f_seqs="${sample}/mut_umis/forward_seqs.txt"
        mut_rc_seqs="${sample}/mut_umis/rc_seqs.txt"
        wt_f_seqs="${sample}/wt_umis/forward_seqs.txt"
        wt_rc_seqs="${sample}/wt_umis/rc_seqs.txt"

        good_mut_out="${sample}/mut_umis/${target}_umis_with_mut.sam"
        good_wt_out="${sample}/mut_umis/${target}_umis_no_mut.sam"
        bad_mut_out="${sample}/wt_umis/${target}_umis_with_wt.sam"
        bad_wt_out="${sample}/wt_umis/${target}_umis_no_wt.sam"

        grep -F "$forward" "$mut_seqs" > "$mut_f_seqs"
        grep -F "$reverse_comp" "$mut_seqs" > "$mut_rc_seqs"

        mutf_count=$(wc -l < "$mut_f_seqs")
        mutrc_count=$(wc -l < "$mut_rc_seqs")

        total=$((mutf_count + mutrc_count)) 
        mut_count=$(wc -l < "$mut_sam")
        echo "'$target' total mut_count vs f + rc: '$mut_count'/'$total' ('$mutf_count' and '$mutrc_count' for '$sample')" 

        echo "f: $mutf_count rc: $mutrc_count ($target MUT $sample)"

        # just get the ones that have the right allele first before subsetting sam
        if [ "$mutf_count" -gt "$mutrc_count" ]; then
            awk -v char="$for_mut" -v pos="$for_pos" 'substr($0, pos, 1) == char' "$mut_seqs" > intermediate_v2.txt
        else
            awk -v char="$rc_mut" -v pos="$rc_pos" 'substr($0, pos, 1) == char' "$mut_seqs" > intermediate_v2.txt
        fi

        grep -wFf intermediate_v2.txt "$mut_sam" > "$good_mut_out"
        grep -wvFf intermediate_v2.txt "$mut_sam" > "$bad_mut_out"

        yes_mut=$(wc -l < "$good_mut_out")
        no_mut=$(wc -l < "$bad_mut_out")
        echo "contains mut: $yes_mut doesn't contain: $no_mut"

        grep -F "$forward" "$wt_seqs" > "$wt_f_seqs"
        grep -F "$reverse_comp" "$wt_seqs" > "$wt_rc_seqs"

        wtf_count=$(wc -l < "$wt_f_seqs")
        wtrc_count=$(wc -l < "$wt_rc_seqs")

        total=$((wtf_count + wtrc_count))
        wt_count=$(wc -l < "$wt_sam")
        echo "'$target' total wt_count vs f + rc: '$wt_count'/'$total' ('$wtf_count' and '$wtrc_count' for '$sample')" 

        echo "f: $wtf_count rc: $wtrc_count ($target WT $sample)"

        if [ "$wtf_count" -gt "$wtrc_count" ]; then
            awk -v char="$for_wt" -v pos="$for_pos" 'substr($0, pos, 1) == char' "$wt_seqs" > intermediate_v2.txt
        else
            awk -v char="$rc_wt" -v pos="$rc_pos" 'substr($0, pos, 1) == char' "$wt_seqs" > intermediate_v2.txt
        fi

        grep -wFf intermediate_v2.txt "$wt_sam" > "$good_wt_out"
        grep -wvFf intermediate_v2.txt "$wt_sam" > "$bad_wt_out"

        yes_wt=$(wc -l < "$good_wt_out")
        no_wt=$(wc -l < "$bad_wt_out")
        echo "contains wt: $yes_wt doesn't contain: $no_wt"

    done < std_rev.txt

done


############################################# STEP 4 ############################################################
# reassign any cross-matched records (i.e. a record with MUT probe label but WT allele gets reassigned to WT probe)

for sample in "${SAMPLES[@]}"; do

    while read -r target forward reverse_comp for_wt for_mut for_pos rc_wt rc_mut rc_pos; do

        good_mut_out="${sample}/mut_umis/${target}_umis_with_mut.sam"
        good_wt_out="${sample}/mut_umis/${target}_umis_no_mut.sam"
        bad_mut_out="${sample}/wt_umis/${target}_umis_with_wt.sam"
        bad_wt_out="${sample}/wt_umis/${target}_umis_no_wt.sam"

        recovered_mut="${sample}/mut_umis/${target}_final.sam"
        recovered_wt="${sample}/wt_umis/${target}_final.sam"

        awk -F'\t' '{print $10}' "$bad_mut_out" > intermediate.txt
        grep -F "$forward" intermediate.txt > intermediate_v2.txt
        grep -F "$reverse" intermediate.txt > intermediate_v3.txt

        mutf_count=$(wc -l < intermediate_v2.txt)
        mutrc_count=$(wc -l < intermediate_v3.txt)

        if [ "$mutf_count" -gt "$mutrc_count" ]; then
            awk -v char="$for_wt" -v pos="$for_pos" 'substr($0, pos, 1) == char' intermediate.txt > intermediate_v2.txt
        else
            awk -v char="$rc_wt" -v pos="$rc_pos" 'substr($0, pos, 1) == char' intermediate.txt > intermediate_v2.txt
        fi

        grep -wFf intermediate_v2.txt "$bad_mut_out" > saved_wt_reads.sam
        cat "$good_wt_out" saved_wt_reads.sam > "$recovered_wt"

        awk -F'\t' '{print $10}' "$bad_wt_out" > intermediate.txt
        grep -F "$forward" intermediate.txt > intermediate_v2.txt
        grep -F "$reverse" intermediate.txt > intermediate_v3.txt

        wtf_count=$(wc -l < intermediate_v2.txt)
        wtrc_count=$(wc -l < intermediate_v3.txt)

        if [ "$wtf_count" -gt "$wtrc_count" ]; then
            awk -v char="$for_mut" -v pos="$for_pos" 'substr($0, pos, 1) == char' intermediate.txt > intermediate_v2.txt
        else
            awk -v char="$rc_mut" -v pos="$rc_pos" 'substr($0, pos, 1) == char' intermediate.txt > intermediate_v2.txt
        fi

        grep -wFf intermediate_v2.txt "$bad_wt_out" > saved_mut_reads.sam
        cat "$good_mut_out" saved_mut_reads.sam > "$recovered_mut"


    done < std_rev.txt

done


############################################ STEP 5 (OPTIONAL) ######################################################
# for each target, look at allele distribution for discrepant records (if they don't have the correct allele, determine what allele do they have)


while read -r target forward reverse_comp for_wt for_mut for_pos rc_wt rc_mut rc_pos; do

    for sample in "${SAMPLES[@]}"; do
        
        bad_mut_out="${sample}/wt_umis/${target}_umis_with_wt.sam"
        bad_wt_out="${sample}/wt_umis/${target}_umis_no_wt.sam"

        awk -F'\t' '{print $10}' "$bad_mut_out" > intermediate.txt
        grep -F "$forward" intermediate.txt > intermediate_v2.txt
        grep -F "$reverse" intermediate.txt > intermediate_v3.txt

        mutf_count=$(wc -l < intermediate_v2.txt)
        mutrc_count=$(wc -l < intermediate_v3.txt)

        if [ "$mutf_count" -gt "$mutrc_count" ]; then
            echo "$target ($sample) wt allele: $for_wt"
            cut -c "$for_pos" intermediate.txt | sort | uniq -c
        else
            echo "$target ($sample) wt allele: $rc_wt"
            cut -c "$rc_pos" intermediate.txt | sort | uniq -c
        fi

        awk -F'\t' '{print $10}' "$bad_wt_out" > intermediate.txt
        grep -F "$forward" intermediate.txt > intermediate_v2.txt
        grep -F "$reverse" intermediate.txt > intermediate_v3.txt

        wtf_count=$(wc -l < intermediate_v2.txt)
        wtrc_count=$(wc -l < intermediate_v3.txt)

        if [ "$wtf_count" -gt "$wtrc_count" ]; then
            echo "$target ($sample) mut allele: $for_mut"
            cut -c "$for_pos" intermediate.txt | sort | uniq -c
        else
            echo "$target ($sample) mut allele: $rc_mut"
            cut -c "$rc_pos" intermediate.txt | sort | uniq -c
        fi

    done

done < std_rev.txt
