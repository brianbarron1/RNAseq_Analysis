#/bin/sh

#####
#	GenerateStarRsem_PairedEnd.sh - Generates .PBS files to align paired-end reads
#####

#####Initialize Variables.
#Folder for .PBS Files
PBS_folder="/gpfs/scratch/fas/gerstein/ky26/ky26_brian/Chupp_RNA-seq_Controls/PBS_Files"
#Folder for output alignment files
aligned_output="/gpfs/scratch/fas/gerstein/ky26/ky26_brian/Chupp_RNA-seq_Controls/Aligned_Reads/Paired_End"
#Folder for trimmed reads
data_folder="/gpfs/scratch/fas/gerstein/ky26/ky26_brian/Chupp_RNA-seq_Controls/Trimmed_Reads/Paired_End"

#Create folders for .PBS files and output alignments
mkdir $PBS_folder
mkdir $aligned_output

#Change WD to that directory containing the trimmed data
cd $data_folder 

#Initiate for loop "F" to detect samples and then write corresponding .pbs file

for F in *trimmed; do
	echo "#PBS -q gerstein -Wgroup_list=gerstein" > $PBS_folder/${F}.pbs
	echo "#PBS -m abe -M brian.barron@yale.edu" >> $PBS_folder/${F}.pbs
	echo "#PBS -l nodes=1:ppn=8 -W PARTITION:m610" >> $PBS_folder/${F}.pbs
	echo "mkdir $aligned_output/${F}_aligned" >> $PBS_folder/${F}.pbs
	
	#Enter sample folder to extract read names
	cd $F
	read1="$(ls *1.fq.gz -1)";
        read2="$(ls *2.fq.gz -1)";
	
	#Write line of code to execute STAR_RSEM
	echo "cd $aligned_output/${F}_aligned"  >> $PBS_folder/${F}.pbs
	echo "/home2/bab99/Scripts/STAR_RSEM.sh" "$data_folder/${F}/"$read1 "$data_folder/${F}/"$read2 "/net/gerstein/as2665/Data/STAR2 /net/gerstein/ENCODE3_RNA-seq/RSEM/RSEM str_PE 4 8" >> $PBS_folder/${F}.pbs;
	#Return to directory with trimmed reads
	cd $data_folder

done
