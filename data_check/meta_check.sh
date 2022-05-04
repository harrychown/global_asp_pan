# Function to check for presence
checkPresence () {
	# Store output prefix
	PREFIX=$(echo $1 | rev | cut -d"." -f2 | rev)
	# Remove pre-exisitng files
	if [[ -e "${PREFIX}.check" ]]
	then
		rm $PREFIX.*.check
	fi

	# Check whether names in $1 are found in $2
	# Presumes that the name in $1 is found to be preceeding a comma in $2
	while read isolate
	do
		BOOL=$(grep -c "${isolate}," $2)
		echo "${isolate},${BOOL}" >> "${PREFIX}.check"
		if [[ $BOOL -eq 0 ]]
		then
			echo $isolate >> "${PREFIX}.no_meta.check"
		fi
	done <$1
}

# Extract isolate names from SRA groups
sed 1d sra_meta.csv | awk -F',' '{print $2}' > sra_names.txt

 # Identify what metadata is missing from the data on the SSD and SRA
checkPresence non_ref_ssd_files.txt dropbox_meta.csv
checkPresence sra_names.txt dropbox_meta.csv

# Identify what data is missing from the metadata
cat sra_names.txt non_ref_ssd_files.txt > available_data.names
sed 's/$/,/g' available_data.names > available_data.search
sed 1d dropbox_meta.csv | grep -Fv -f available_data.search | awk -F',' '{print $1}' > dropbox.no_data.check

mkdir output
mv available_data* *.check sra_names.txt ./output/


