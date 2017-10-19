library(phyloseq)

#Import data from QIIME. It's a .biom file that contains the OTU, taxonomy and metadata for the Argentina samples for the forward and reverse read separately.

fdata<-import_biom('/Volumes/USB DISK/MGM_OTUpick_forward/otu_table_mc2_metadata_wtax.biom')
rdata<-import_biom('/Volumes/USB DISK/MGM_OTUpick_reverse/otu_table_mc2_metadata_wtax.biom')

#Check that those 3 objects are in the file and what they look like

sample_data(fdata)
head(otu_table(fdata))
head(tax_table(fdata))

#Check column names of taxonomy table

colnames(tax_table(fdata))

#Change column names of taxonomy table

colnames(tax_table(fdata))<-c("Kingdom", "Phylum", "Class", 
"Order", "Family", "Genus","Species")
colnames(tax_table(rdata))<-c("Kingdom", "Phylum", "Class", 
"Order", "Family", "Genus","Species")

library(ggplot2)

#Ordinate forward and reverse data

fdata_pcoa<-ordinate(physeq=fdata, method="PCoA", distance="bray")
rdata_pcoa<-ordinate(physeq=rdata, method="PCoA", distance="bray")

#Make PCoA for forward and reverse reads

plot_ordination(physeq=fdata, ordination=fdata_pcoa, color="Area", shape="Plant.type")

#Organize data to make bar plot by phylum of Forward reads

f_phylum <- fdata %>%
#agglomerate at phylum level 
tax_glom(taxrank = "Phylum") %>%      
#Transform to relative abundance              
transform_sample_counts(function(x) {x/sum(x)} ) %>% 
#Melt to long format
psmelt() %>%      
#Filter low abundance phyla                                   
filter(Abundance > 0.02) %>%     
#Sort alphabetically                    
arrange(Phylum)	

#Repeat organization of Phyla for Reverse read

r_phylum <- fdata %>%
#agglomerate at phylum level 
tax_glom(taxrank = "Phylum") %>%      
#Transform to relative abundance              
transform_sample_counts(function(x) {x/sum(x)} ) %>% 
#Melt to long format
psmelt() %>%      
#Filter low abundance phyla                                   
filter(Abundance > 0.02) %>%     
#Sort alphabetically                    
arrange(Phylum)	


#######try different bar plots for phyla

#Phyla for F and R by sample
ggplot(f_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + geom_bar(stat = "identity")
ggplot(r_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + geom_bar(stat = "identity")

#Phyla for F and R by plant type
ggplot(r_phylum, aes(x = Plant.type, y = Abundance, fill = Phylum)) + geom_bar(stat = "identity")
ggplot(f_phylum, aes(x = Plant.type, y = Abundance, fill = Phylum)) + geom_bar(stat = "identity")

#Phyla for F and R by area							
ggplot(r_phylum, aes(x = Area, y = Abundance, fill = Phylum)) + geom_bar(stat = "identity")
ggplot(f_phylum, aes(x = Area, y = Abundance, fill = Phylum)) + geom_bar(stat = "identity")