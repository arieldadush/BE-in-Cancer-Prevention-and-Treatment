#!/bin/bash
#this script add three columns to the table first col is the bases before the mutation -5' end,
# second col  is the bases after the mutation 3' end 
# third col is the three bases including the mutation according to the frame 
#note: the mutation is sometimes in the first col sometimes in the second it is dosen't matter because for get the blast i merge the two col
#the reference genome
genome="/private/dropbox/Genomes/Human/hg19/hg19.fa";
#"/private/dropbox/Genomes/Human/hg38/hg38.fa"; 
#A temporary BED file used to find the sequence 
out_temp="$PWD/TEMP.bed";
#create the tmp file
touch $out_temp;
#open the of the bedtools intersect and loop runs on lines
# old ver: "/private5/Projects/dadush/clinvar_base_editing/Erez_15_03_23/Erez_merged.tsv"
# old ver: "/private5/Projects/dadush/clinvar_base_editing/Erez_15_03_23/Erez_merged_ver_2.tsv"
#old ver: "/private5/Projects/dadush/clinvar_base_editing/Erez_15_03_23/Erez_without_merged_ver_1.tsv"
# old ver: "/private6/Projects/dadush/cancer_mutation/sample_of_100_of_tables3_all_with_strand_ready_to_add_dna_seq.tsv"
# old ver: "/private6/Projects/dadush/cancer_mutation/tables3_all_with_strand_ready_to_add_dna_seq.tsv"
# ols ver : "/private6/Projects/dadush/cancer_mutation/table_MSK_for_ariel_with_strand_ready_to_add_dna_seq.tsv"
#"/private6/Projects/dadush/cancer_mutation/msk_entire_info_for_ariel_with_strand_ready_to_add_dna_seq.tsv"
cat "/private6/Projects/dadush/cancer_mutation/tables3_all_with_strand_ready_to_add_dna_seq_fix.tsv" | while read -r line; do 


#take the base of mutation as as a bed6 format
	reg_base=$(echo -e $line |cut -f1-6 | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" $6}' ); 
	#create variable contain the strand
	strand=$(echo -e $line | awk '{ print $6}');
	#number of basese that have before and after the mutation 
	bases_before_and_after=20;
	if [[ $strand != "+" ]] && [[ $strand != "-" ]]; then
	strand=$(echo -e $line | rev |cut -d ' ' -f9 | rev)
	fi
	#change evrey  space delimiter in tab delimiter
	echo -e $reg_base| sed 's/ /\t/g' > $out_temp;
	#if the mutation on the plus strand
	if [[ $strand == "+" ]];  then 
	#get the range before the mutation and create it as a bed6 format
		long_range=$(echo -e $reg_base | awk '{print $1 "\t" ($2 - '$bases_before_and_after') "\t" ($3 + '$bases_before_and_after')"\t" "." "\t" "." "\t" $6}')
		##change evrey  space delimiter in tab delimiter and save it to the temp file
		echo -e $long_range | sed 's/ /\t/g' > $out_temp;
		#get the sequence on the right strand and print it as upper case
		long_seq=$(bedtools getfasta -fi $genome -bed $out_temp| awk 'NR==2' | awk '{print toupper($0)}')
		#add the sequence before to the original line
		echo -e "${line}\t${long_seq}";

		#if the mutation is on minus strand
	elif [[ $strand == "-" ]]; then 
	#find the sequence before the mutaion 5' end in the minus strand meaning that i take the oppsite in the plus strand
		long_range=$(echo -e $reg_base | awk -v strand=$strand -v bases_before_and_after=$bases_before_and_after '{print $1 "\t" ($2 - '$bases_before_and_after') "\t" ($3 + '$bases_before_and_after') "\t" "." "\t" "."}')
		#had a bag with the strand so i add it after the bed6 format
		long_range=${long_range}" ${strand}"
		echo -e $long_range | sed 's/ /\t/g' > $out_temp;
		long_seq=$(bedtools getfasta -fi $genome -bed $out_temp -s | awk 'NR==2' | awk '{print toupper($0)}')
		echo -e "${line}\t${long_seq}";
	else
	echo -e "${line}\tnAn"
	fi;
done
