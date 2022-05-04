# Simple script to update any filenames incase of naming errors causing data to be labelled as "missing"
# Make a copy of the file Jo sent over, change the formatting and remove any BritishIsles genomes
tr -d $'\r' < dropbox.no_data.altnames | grep "#N/A" > dropbox.no_data.altnames.copy
SSDNAME=$(sed 's/RBH_/RBH/g' dropbox.no_data.altnames.copy | sed 's/MRL_/MRL/g' | \
sed 's/ARAF03 /ARAF03/g' | sed 's/CXH_/CXH/g' | sed 's/CBS /CBS/g' | awk '{print $3}' | sed 's/Og_code/SSD_code/')
echo "$SSDNAME" > updated.names 
paste dropbox.no_data.altnames.copy updated.names | column -s $'\t' -t > dropbox.no_data.updatednames

# Check to see if these files are now comparable
sed 1d sra_meta.csv | awk -F',' '{print $2}' > sra_names.txt
cat sra_names.txt non_ref_ssd_files.txt > available_data.names
sed 1d updated.names | grep -Fv -f available_data.names  > dropbox.updated.no_data.check2
