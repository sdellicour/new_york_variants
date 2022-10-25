# 1. Selection of the background genomic sequences to include
# 2. Preparing the inputs for the preliminary discrete phylogeographic analyses 
# 3. Analysing the number of distinct intropduction events (not performed)
# 4. Getting the MCC trees from the preliminary discrete phylogeographic analyses
# 5. Identifying the different clusters (clades following introduction events)
# 6. Preparing the discrete phylogeographic analyses among counties
# 7. Preparing the continuous phylogeographic analyses (RRW, Cauchy model)
# 8. Building the maximum clade consensus (MCC) trees for both analyses
# 9. Extracting spatio-temporal information embedded in MCC and posterior trees
# 10. Visualisations of the discrete and continuous phylogeographic reconstructions
# 11. Dispersal statistics based on the continuous phylogeographic reconstructions
# 12. Dispersal statistics based on the discrete phylogeographic reconstructions

library(diagram)
library(lubridate)
library(seraphim)
library(treeio)

writingFiles = TRUE; showingPlots = TRUE
writingFiles = FALSE; showingPlots = FALSE

nberOfExtractionFiles = 900

NYC_counties = c("NewYork","Bronx","Kings","Queens","Richmond","Nassau","Suffolk","Westchester")
variants = c("Alpha","Beta","Delta","Epsilon","Eta","Gamma","Kappa","Lambda","Mu","O-BA1","O-BA2","O-BA3","Theta","Zeta")
variants = c("Iota","Alpha","Delta","O-BA1")

# 1. Selection of the background genomic sequences to include

	# 1.1. Selection of the Nextstrain genomic sequences to include

files = list.files("All_NYU_alignment_data/All_NYU_alignment_data_26032022")
files = files[which(grepl(".tsv",files))] # to do: convert the ".tsv" into ".csv" files (";" --> "," and "\t" --> ";")
for (i in 1:length(files))
	{
		newName = gsub(".tsv",".csv",files[i])
		tab = read.csv(paste0("All_NYU_alignment_data/All_NYU_alignment_data_26032022/",newName), sep=";")
		dates = mdy(tab[,"date"]); missingDates = tab[which(is.na(dates)),"date"]
		cat(files[i],"  -  ",as.character(max(dates, na.rm=T)),"\n",sep="")
	}
system("`curl http://data.nextstrain.org/ncov_gisaid_north-america_2022-02-08.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2022-02-08_Alpha.json`") # 2022-02-08
system("`curl http://data.nextstrain.org/ncov_gisaid_north-america_2021-11-11.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2021-11-11_Beta.json`") # 2021-11-11
system("`curl http://data.nextstrain.org/ncov_gisaid_north-america_2022-02-21.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2022-02-21_Delta.json`") # 2022-02-21
system("`curl http://data.nextstrain.org/ncov_north-america_2021-05-26.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2021-05-26_Epsilon.json`") # 2021-05-26
system("`curl http://data.nextstrain.org/ncov_north-america_2021-06-02.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2021-06-02_Eta.json`") # 2021-05-29
system("`curl http://data.nextstrain.org/ncov_gisaid_north-america_2021-10-25.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2021-10-25_Gamma.json`") # 2021-10-23
system("`curl http://data.nextstrain.org/ncov_gisaid_north-america_2021-08-04.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2021-08-04_Iota.json`") # 2021-08-03
system("`curl http://data.nextstrain.org/ncov_north-america_2021-05-13.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2021-05-13_Kappa.json`") # 2021-05-13
system("`curl http://data.nextstrain.org/ncov_gisaid_north-america_2021-09-03.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2021-09-03_Lambda.json`") # 2021-09-02
system("`curl http://data.nextstrain.org/ncov_gisaid_north-america_2021-12-24.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2021-12-24_Mu.json`") # 2021-12-24
system("`curl http://data.nextstrain.org/ncov_gisaid_north-america_2022-02-27.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2022-02-27_O-BA1.json`") # 2022-02-27
system("`curl http://data.nextstrain.org/ncov_gisaid_north-america_2022-02-27.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2022-02-27_O-BA2.json`") # 2022-02-27
system("`curl http://data.nextstrain.org/ncov_gisaid_north-america_2021-12-06.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2021-12-06_O-BA3.json`") # 2021-12-06
system("`curl http://data.nextstrain.org/ncov_gisaid_north-america_2022-02-10.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2022-02-10_Other.json`") # 2022-02-10
system("`curl http://data.nextstrain.org/ncov_north-america_2021-06-22.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2021-06-22_Theta.json`") # 2021-06-22
system("`curl http://data.nextstrain.org/ncov_north-america_2021-05-11.json --compressed -o All_Nextstrain_extractions/Nextstrain_NA_2021-05-11_Zeta.json`") # 2021-05-09
files = list.files("All_Nextstrain_extractions/"); files = files[which(grepl(".json",files))]
for (i in 1:length(files))
	{
		newName = gsub(".json",".txt",files[i])
		txt = scan(paste0("All_Nextstrain_extractions/",files[i]), what="", sep="\n", quiet=T)
		IDs = gsub(" ","",gsub("\",","",gsub("\"name\": \"","",txt[which(grepl("\"name\":",txt))])))
		IDs = IDs[which((IDs!="GISAID\"")&(IDs!="theNextstrainteam")&(!grepl("NODE_",IDs)))]
		write(IDs, paste0("All_Nextstrain_extractions/",newName))
	}

	# 1.2. Filtering the New York sequences to decrease the size of some alignments

tab = read.csv("All_NYU_alignment_data/All_NYU_alignment_data_26032022/Delta_metadata.csv", head=T, sep=";")
indices = which((tab[,"division"]=="New York")&((tab[,"location"]=="New York City")|(tab[,"location"]=="")|(is.na(tab[,"location"]))))
sequencesToDiscard = tab[indices,"strain"]; sequencesToDiscard = paste0(">",sequencesToDiscard); tab = tab[-indices,]
write.csv(tab, "All_NYU_alignment_data/All_NYU_alignment_data_26032022/Delta_metadata_light.csv")
seqs = scan("All_NYU_alignment_data/All_NYU_alignment_data_26032022/Delta_sequences.fasta", what="", sep="\n", quiet=T)
indices = which(grepl(">",seqs)); sink("All_NYU_alignment_data/All_NYU_alignment_data_26032022/Delta_sequences_light.fasta")
for (i in 1:length(indices))
	{
		if (!seqs[indices[i]]%in%sequencesToDiscard)
			{
				cat(seqs[indices[i]],"\n",sep="")
				if (i < length(indices)) cat(seqs[(indices[i]+1):(indices[i+1]-1)],"\n",sep="")
				if (i == length(indices)) cat(seqs[(indices[i]+1):length(seqs)],"\n",sep="")
			}
	}
sink(NULL)

tab = read.csv("All_NYU_alignment_data/All_NYU_alignment_data_26032022/Iota_metadata.csv", head=T, sep=";")
indices = which((tab[,"division"]=="New York")&((tab[,"location"]=="New York City")|(tab[,"location"]=="")|(is.na(tab[,"location"]))))
sequencesToDiscard = tab[indices,"strain"]; sequencesToDiscard = paste0(">",sequencesToDiscard); tab = tab[-indices,]
write.csv(tab, "All_NYU_alignment_data/All_NYU_alignment_data_26032022/Iota_metadata_light.csv")
seqs = scan("All_NYU_alignment_data/All_NYU_alignment_data_26032022/Iota_sequences.fasta", what="", sep="\n", quiet=T)
indices = which(grepl(">",seqs)); sink("All_NYU_alignment_data/All_NYU_alignment_data_26032022/Iota_sequences_light.fasta")
for (i in 1:length(indices))
	{
		if (!seqs[indices[i]]%in%sequencesToDiscard)
			{
				cat(seqs[indices[i]],"\n",sep="")
				if (i < length(indices)) cat(seqs[(indices[i]+1):(indices[i+1]-1)],"\n",sep="")
				if (i == length(indices)) cat(seqs[(indices[i]+1):length(seqs)],"\n",sep="")
			}
	}
sink(NULL)

tab = read.csv("All_NYU_alignment_data/All_NYU_alignment_data_26032022/O-BA1_metadata.csv", head=T, sep=";")
indices = which((tab[,"division"]=="New York")&((tab[,"location"]=="New York City")|(tab[,"location"]=="")|(is.na(tab[,"location"]))))
sequencesToDiscard = tab[indices,"strain"]; sequencesToDiscard = paste0(">",sequencesToDiscard); tab = tab[-indices,]
write.csv(tab, "All_NYU_alignment_data/All_NYU_alignment_data_26032022/O-BA1_metadata_light.csv")
seqs = scan("All_NYU_alignment_data/All_NYU_alignment_data_26032022/O-BA1_sequences.fasta", what="", sep="\n", quiet=T)
indices = which(grepl(">",seqs)); sink("All_NYU_alignment_data/All_NYU_alignment_data_26032022/O-BA1_sequences_light.fasta")
for (i in 1:length(indices))
	{
		if (!seqs[indices[i]]%in%sequencesToDiscard)
			{
				cat(seqs[indices[i]],"\n",sep="")
				if (i < length(indices)) cat(seqs[(indices[i]+1):(indices[i+1]-1)],"\n",sep="")
				if (i == length(indices)) cat(seqs[(indices[i]+1):length(seqs)],"\n",sep="")
			}
	}
sink(NULL)

# 2. Preparing the inputs for the preliminary discrete phylogeographic analyses 

for (h in 1:length(variants))
	{
		# tree = read.nexus(paste0("IQ-TREE_TreeTime_runs/",variants[h],"_TreeTime.tre"))
		tree = read.nexus(paste0("Thorney_BEAST_analyses/",variants[h],"_Thorney.tree"))
		data1 = read.csv(paste0("IQ-TREE_TreeTime_runs/",variants[h],"_TreeTime.csv"), head=T)
		data2 = read.csv(paste0("All_NYU_alignment_data/All_NYU_alignment_data_26032022/",variants[h],"_metadata.csv"), head=T, sep=";")
		seqIDs1 = tree$tip.label; seqIDs2 = data1[,"name"]; NYtipsToExclude = c()
		locations = rep(NA, length(seqIDs1)); collectionDates = rep(NA, length(seqIDs1))
		data2[,"location"] = gsub(" ","",data2[,"location"])
		for (i in 1:length(seqIDs1))
			{
				if (grepl("USA",seqIDs1[i]))
					{
						state = unlist(strsplit(gsub("hCoV-19/USA/","",seqIDs1[i]),"\\/"))[1]
						state = unlist(strsplit(state,"-"))[1]
						if (state == "NY")
							{
								seqID = unlist(strsplit(seqIDs1[i],"\\|"))[1]
								index = which(data2[,"strain"]==seqID)
								if (length(index) == 1)
									{
										if (data2[index,"division"] == "New York")
											{
												county = data2[index,"location"]
												if (grepl("County",county))
													{
														county = gsub("County","",gsub(" ","",county))
														county = gsub("County","",gsub(" ","",county))
														if (county%in%NYC_counties)
															{
																locations[i] = "NYC_counties"
															}	else	{
																locations[i] = "other"
															}
													}	else	{
														if ((county == "New York")&(!grepl("NYULH",seqIDs1[i])))
															{
																county = "Unknown NYC county"
															}
														county = gsub(" ","",county)
														if (county%in%NYC_counties)
															{
																locations[i] = "NYC_counties"
															}	else	{
																NYtipsToExclude = rbind(NYtipsToExclude, cbind(seqIDs1[i],data2[index,"division"],data2[index,"location"]))
															}	
													}
											}	else	{
												locations[i] = "other"
											}
									}	else	{
										NYtipsToExclude = rbind(NYtipsToExclude, cbind(seqIDs1[i],NA,NA))
									}
							}	else	{
								locations[i] = "other"
							}
					}	else	{
						locations[i] = "other"
					}
				if (length(which(seqIDs2==seqIDs1[i])) != 1) print(c(h,i))
				collectionDates[i] = data1[which(seqIDs2==seqIDs1[i])[1],"date"]
			}
		tab = cbind(seqIDs1,locations,collectionDates)
		colnames(tab) = c("trait","location","collection_date")
		tab = tab[-which(tab[,"trait"]%in%NYtipsToExclude),]; txt = c()
		for (i in 1:dim(tab)[1])
			{
				if (!grepl("NY",tab[i,"location"])) tab[i,"location"] = "other"
				txt = c(txt, paste0(">",tab[i,"trait"]),"NNNN")
			}
		tree = drop.tip(tree, NYtipsToExclude[,1])
		if (writingFiles) write.table(tab, paste0("Preliminary_discrete_runs/",variants[h],"_DTA.txt"), row.names=F, quote=F, sep="\t")
		if (writingFiles) write(txt, paste0("Preliminary_discrete_runs/",variants[h],"_DTA.fasta"))
		if (writingFiles) write.tree(tree, paste0("Preliminary_discrete_runs/",variants[h],"_DTA.tre"))
	}
for (h in 1:length(variants))
	{
		tab = read.table(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.txt"), head=T)
		tre = read.tree(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.tre"))
		sequencesToDiscard = c()
		for (i in 1:dim(tab)[1])
			{
				if (length(unlist(strsplit(tab[i,"collection_date"],"-"))) != 3)
					{
						sequencesToDiscard = c(sequencesToDiscard, tab[i,"trait"])
					}
				if (tab[i,"collection_date"] == "XXXX-XX-XX") # at least one case for Zeta
					{
						date = unlist(strsplit(tab[i,"trait"],"\\|"))
						tab[i,"collection_date"] = date[length(date)]
					}
			}
		tab = tab[which(!tab[,"trait"]%in%sequencesToDiscard),]
		tre = drop.tip(tre, sequencesToDiscard); # print(c(h,sequencesToDiscard))
		write.tree(tre, paste0("Preliminary_discrete_runs/",variants[h],"_DTA.tre"))
		tre = scan(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.tre"), what="", sep="\n", quiet=T)
		xml = scan("Template_DTA_runs.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
		xml = gsub("TEMPLATE",paste0(variants[h],"_DTA"),xml)
		sink(file=paste0("Preliminary_discrete_runs/",variants[h],"_DTA.xml"))
		for (i in 1:length(xml))
			{
				cat(xml[i],"\n")
				if (grepl("<taxa id=\"taxa\">",xml[i]))
					{
						for (j in 1:dim(tab)[1])
							{
								cat("\t\t<taxon id=\"",tab[j,"trait"],"\">\n",sep="")
								cat("\t\t\t<date value=\"",decimal_date(ymd(tab[j,"collection_date"])),"\" direction=\"forwards\" units=\"years\"/>\n",sep="")
								cat("\t\t\t<attr name=\"location\">\n",sep="")
								cat("\t\t\t\t",tab[j,"location"],"\n",sep="")
								cat("\t\t\t</attr>\n",sep="")
								cat("\t\t</taxon>\n",sep="")
							}
					}
				if (grepl("<alignment id=\"alignment\" dataType=\"nucleotide\">",xml[i]))
					{
						for (j in 1:dim(tab)[1])
							{
								cat("\t\t<sequence>\n",sep="")
								cat("\t\t\t<taxon idref=\"",tab[j,"trait"],"\"/>\n",sep="")
								cat("\t\t\t\tNNNN\n",sep="")
								cat("\t\t</sequence>\n",sep="")
							}
					}
				if (grepl("<newick id=\"startingTree\">",xml[i]))
					{
						cat("\t\t",tre,"\n",sep="")
					}
			}
		sink(NULL)
	}

# 3. Analysing the number of distinct intropduction events (not performed)

burnIns = rep(31, length(variants))
for (h in 1:length(variants))
	{
		trees = scan(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		data = read.table(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.txt"), head=T)
		indices1 = which(!grepl("tree STATE_",trees)); indices2 = which(grepl("tree STATE_",trees))
		mostRecentSamplingDatum = max(decimal_date(ymd(data[,"collection_date"])), na.rm=T)
		NYC_branches_list = rep(NA,length(trees)); NYC_introductions_list = rep(NA,length(trees))
		NYC_tipBranches_list = rep(NA,length(trees)); NYC_tMRCAs_list = list()
		clusters_of_one_sequence_list = rep(NA,length(trees))
		for (i in (burnIns[h]+1):length(indices2))
			{
				tree1 = trees[c(indices1[1:(length(indices1)-1)],indices2[i],indices1[length(indices1)])]
				write(tree1, paste0(variants[h],"_",i,".tree"))
				tree2 = readAnnotatedNexus(paste0(variants[h],"_",i,".tree"))
				NYC_branches = 0; NYC_introductions = 0; NYC_tipBranches = 0
				NYC_tMRCAs = c(); clusters_of_one_sequence = 0
				for (j in 1:dim(tree2$edge)[1])
					{
						if ((!is.null(tree2$annotations[[j]]$location))&&(tree2$annotations[[j]]$location == "NYC_counties"))
							{
								NYC_branches = NYC_branches + 1
								index = which(tree2$edge[,2]==tree2$edge[j,1])
								if ((!is.null(tree2$annotations[[index]]$location))&&(tree2$annotations[[index]]$location != "NYC_counties"))
									{
										NYC_introductions = NYC_introductions + 1
										tMRCA = mostRecentSamplingDatum-nodeheight(tree2,tree2$edge[j,1])
										NYC_tMRCAs = c(NYC_tMRCAs, tMRCA)
										if (!tree2$edge[j,2]%in%tree2$edge[,1])
											{
												clusters_of_one_sequence = clusters_of_one_sequence + 1
											}
									}
								if (!tree2$edge[j,2]%in%tree2$edge[,1])
									{
										NYC_tipBranches = NYC_tipBranches + 1
									}
							}
					}
				NYC_branches_list[i] = NYC_branches
				NYC_introductions_list[i] = NYC_introductions
				NYC_tipBranches_list[i] = NYC_tipBranches
				NYC_tMRCAs_list[[i]] = NYC_tMRCAs
				clusters_of_one_sequence_list[i] = clusters_of_one_sequence
				file.remove(paste0(variants[h],"_",i,".tree"))
			}
		quantiles = quantile(NYC_introductions_list[!is.na(NYC_introductions_list)],probs=c(0.025,0.975))
		cat("A minimum number of ",median(NYC_introductions_list[!is.na(NYC_introductions_list)])," lineage introductions (95% HPD interval = [",
			quantiles[1],"-",quantiles[2],"])"," identified from the global phylogenetic analysis of ",NYC_tipBranches," ",variants[h]," sequences sampled in the NYC area","\n",sep="")
		proportions_of_cluster_1 = clusters_of_one_sequence_list/NYC_introductions_list
		quantiles = quantile(proportions_of_cluster_1[!is.na(proportions_of_cluster_1)],probs=c(0.025,0.975))
		cat("Proportion of clusters of n = 1: ",median(proportions_of_cluster_1[!is.na(proportions_of_cluster_1)])," (95% HPD interval = [",quantiles[1],"-",quantiles[2],"])","\n",sep="")
	}
registerDoMC(cores=10)
for (h in 1:length(variants))
	{
		trees = scan(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		data = read.table(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.txt"), head=T)
		indices1 = which(!grepl("tree STATE_",trees)); indices2 = which(grepl("tree STATE_",trees))
		mostRecentSamplingDatum = max(decimal_date(ymd(data[,"collection_date"])), na.rm=T)
		NYC_introductions_list = rep(NA,length(trees)); clusters_of_one_sequence_list = rep(NA,length(trees))
		buffer = foreach(i = (burnIns[h]+1):length(indices2)) %dopar% {
		# for (i in (burnIns[h]+1):length(indices2)) {
				tree1 = trees[c(indices1[1:(length(indices1)-1)],indices2[i],indices1[length(indices1)])]
				write(tree1, paste0(variants[h],"_",i,".tree"))
				tree2 = readAnnotatedNexus(paste0(variants[h],"_",i,".tree"))
				NYC_branches = 0; NYC_introductions = 0
				NYC_tipBranches = 0; clusters_of_one_sequence = 0
				for (j in 1:dim(tree2$edge)[1])
					{
						if ((!is.null(tree2$annotations[[j]]$location))&&(tree2$annotations[[j]]$location == "NYC_counties"))
							{
								NYC_branches = NYC_branches + 1
								index = which(tree2$edge[,2]==tree2$edge[j,1])
								if ((!is.null(tree2$annotations[[index]]$location))&&(tree2$annotations[[index]]$location != "NYC_counties"))
									{
										NYC_introductions = NYC_introductions + 1
										if (!tree2$edge[j,2]%in%tree2$edge[,1])
											{
												clusters_of_one_sequence = clusters_of_one_sequence + 1
											}
									}
								if (!tree2$edge[j,2]%in%tree2$edge[,1])
									{
										NYC_tipBranches = NYC_tipBranches + 1
									}
							}
					}
				file.remove(paste0(variants[h],"_",i,".tree"))
				cbind(NYC_tipBranches, NYC_introductions, clusters_of_one_sequence)
				cbind(runif(1,1,10),runif(1,1,10),runif(1,1,10))
			}
		NYC_tipBranches = buffer[[length(buffer)]][1]
		for (i in (burnIns[h]+1):length(indices2))
			{
				NYC_introductions_list[[i]] = buffer[[i]][2]
				clusters_of_one_sequence_list[[i]] = buffer[[i]][3]
			}
		quantiles = quantile(NYC_introductions_list[!is.na(NYC_introductions_list)],probs=c(0.025,0.975))
		cat("A minimum number of ",median(NYC_introductions_list[!is.na(NYC_introductions_list)])," lineage introductions (95% HPD interval = [",
			quantiles[1],"-",quantiles[2],"])"," identified from the global phylogenetic analysis of ",NYC_tipBranches," ",variants[h]," sequences sampled in the NYC area",sep="")
		proportions_of_cluster_1 = clusters_of_one_sequence_list/NYC_introductions_list
		quantiles = quantile(proportions_of_cluster_1[!is.na(proportions_of_cluster_1)],probs=c(0.025,0.975))
		cat("Proportion of clusters of n = 1: ",median(proportions_of_cluster_1[!is.na(proportions_of_cluster_1)])," (95% HPD interval = [",quantiles[1],"-",quantiles[2],"])")
	}
		# A minimum number of 297 lineage introductions (95% HPD interval = [285-315]) identified from the global phylogenetic analysis of 2519 Iota sequences sampled in the NYC area.
			# Proportion of clusters of n = 1: 0.723 (95% HPD interval = [0.697-0.745])
		# A minimum number of 787 lineage introductions (95% HPD interval = [754-835]) identified from the global phylogenetic analysis of 1655 Alpha sequences sampled in the NYC area.
			# Proportion of clusters of n = 1: 0.688 (95% HPD interval = [0.675-0.702])
		# A minimum number of 3258 lineage introductions (95% HPD interval = [3227-3292]) identified from the global phylogenetic analysis of 4796 O-BA1 sequences sampled in the NYC area.
			# Proportion of clusters of n = 1: 0.763 (95% HPD interval = [0.757-0.769])
		# A minimum number of 2368 lineage introductions (95% HPD interval = [2333-2397]) identified from the global phylogenetic analysis of 3123 O-BA1 sequences sampled in the NYC area.
			# Proportion of clusters of n = 1: 0.833 (95% HPD interval = [0.827-0.840])

# 4. Getting the MCC trees from the preliminary discrete phylogeographic analyses

	# Note: issue due to RAM limitation (Java heap space) for the larger trees. Solution: editing the "Info.plist" files (TreeAnnotator app --> right click --> show content)
	# 		to increase the RAM that can be allocated to TreeAnnotator (e.g. --> Xmx20480M)

burnIns = rep(31, length(variants))
wd = getwd(); setwd(paste0(wd,"/Preliminary_discrete_runs/"))
for (h in 1:length(variants))
	{
		system(paste0("BEAST_1104/bin/treeannotator -burninTrees ",burnIns[h]," -heights keep ",variants[h],"_DTA.trees ",variants[h],"_DTA.tree"), ignore.stdout=F, ignore.stderr=F)
	}
setwd(wd)

# 5. Identifying the different clusters (clades following introduction events)

clusters1_list = list(); clusters2_list = list(); NYC_introductions_list = list()
zipCodes = shapefile("NY_state_all_shapefiles/ZipCodes_US.shp")
for (h in 1:length(variants))
	{
		if (file.exists(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.tree")))
			{
				tree = readAnnotatedNexus(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.tree")); indices = c()
			}	else	{
				tree = readAnnotatedNexus(paste0("Preliminary_discrete_runs/",variants[h],"_last.tree")); indices = c()
			}
		metadata1 = read.csv(paste0("All_NYU_alignment_data/All_NYU_alignment_data_26032022/",variants[h],"_metadata.csv"), head=T, sep=";")
		metadata2 = read.csv(paste0("IQ-TREE_TreeTime_runs/",variants[h],"_TreeTime.csv"), head=T, sep=",")
		if (file.exists(paste0("IQ-TREE_TreeTime_runs/",variants[h],"_lineages.csv")))
			{
				metadata3 = read.csv(paste0("IQ-TREE_TreeTime_runs/",variants[h],"_lineages.csv"), head=T, sep=",")
			}	else	{
				metadata3 = metadata2 # for those variants the "lineage" column is included in the main metadat file
				colnames(metadata3)[1] = "taxon"
			}
		metadata4 = read.table(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.txt"), head=T)
		mostRecentSamplingDatum = max(decimal_date(ymd(metadata4[,"collection_date"])), na.rm=T)
		NYC_branches = c(); NYC_introductions = c(); NYC_TipBranches = c(); sampledSequences = c()
		for (i in 1:dim(tree$edge)[1])
			{
				if (is.null(tree$annotations[[i]]))
					{
						print(c(h,i))
					}	else		{
						if (tree$annotations[[i]]$location == "NYC_counties")
							{
								NYC_branches = c(NYC_branches,i)
								index = which(tree$edge[,2]==tree$edge[i,1])
								if (tree$annotations[[index]]$location != "NYC_counties")
									{
										NYC_introductions = c(NYC_introductions, i)
									}
								if (!tree$edge[i,2]%in%tree$edge[,1])
									{
										NYC_TipBranches = c(NYC_TipBranches, i)
										sampledSequences = c(sampledSequences, tree$tip.label[tree$edge[i,2]])
									}
							}
					}
			}
		NYC_introductions_list[[h]] = NYC_introductions
		for (i in 1:length(NYC_introductions))
			{
				if (i == 1) clusters1 = list()
				if (tree$edge[NYC_introductions[i],2]%in%tree$edge[,1])
					{
						subtree = tree_subset(tree, tree$edge[NYC_introductions[i],2], levels_back=0)
						clusters1[[i]] = gsub("'","",subtree$tip.label)
					}	else		{
						clusters1[[i]] = gsub("'","",tree$tip.label[tree$edge[NYC_introductions[i],2]])
					}
			}
		for (i in 2:length(clusters1))
			{
				for (j in 1:(i-1))
					{
						if (sum(clusters1[[i]]%in%clusters1[[j]]) == length(clusters1[[i]]))
							{
								clusters1[[j]] = clusters1[[j]][which(!clusters1[[j]]%in%clusters1[[i]])]
							}
						if (sum(clusters1[[j]]%in%clusters1[[i]]) == length(clusters1[[j]]))
							{
								clusters1[[i]] = clusters1[[i]][which(!clusters1[[i]]%in%clusters1[[j]])]
							}
					}
			}
		sampledSequences = gsub("'","",sampledSequences)
		if (!file.exists(paste0("Preliminary_discrete_runs/",variants[h],"_data.csv")))
			{
				samplingData = matrix(nrow=length(sampledSequences), ncol=9)
				colnames(samplingData) = c("sequence_ID","collection_date","lineage","state","county","location","zip_code","latitude","longitude")
				samplingData[,"sequence_ID"] = sampledSequences
				for (i in 1:dim(samplingData)[1])
					{
						if (sum(samplingData[i,"sequence_ID"]==metadata3[,"taxon"]) == 1)
							{
								index = which(metadata3[,"taxon"]==samplingData[i,"sequence_ID"])
								samplingData[i,"lineage"] = metadata3[index,"lineage"]
							}
						if (sum(samplingData[i,"sequence_ID"]==metadata2[,"name"]) == 1)
							{
								index = which(metadata2[,"name"]==samplingData[i,"sequence_ID"])
								samplingData[i,"collection_date"] = decimal_date(ymd(metadata2[index,"date"]))
							}
						if (sum(unlist(strsplit(samplingData[i,"sequence_ID"],"\\|"))[1]==metadata1[,"strain"]) == 1)
							{
								index = which(metadata1[,"strain"]==unlist(strsplit(samplingData[i,"sequence_ID"],"\\|"))[1])
								samplingData[i,"state"] = metadata1[index,"division"]
								samplingData[i,"county"] = gsub("County","",gsub(" ","",metadata1[index,"location"]))
								samplingData[i,"zip_code"] = metadata1[index,"zip"]
								if (!is.na(samplingData[i,"county"]))
									{
										if ((samplingData[i,"county"] == "")|(samplingData[i,"county"] == "New York City")) samplingData[i,"county"] = NA
									}
								if (!is.na(samplingData[i,"zip_code"]))
									{
										if (samplingData[i,"zip_code"] == "unknown") samplingData[i,"zip_code"] = NA
									}
							}
						if (!is.na(samplingData[i,"county"]))
							{
								if ((samplingData[i,"state"]=="New York")&(samplingData[i,"county"]%in%NYC_counties))
									{
										samplingData[i,"location"] = samplingData[i,"county"]
									}
							}
						if ((!is.na(samplingData[i,"zip_code"]))&(samplingData[i,"county"]%in%NYC_counties))
							{
								zipCode = gsub(" ","",unlist(strsplit(samplingData[i,"zip_code"],"-"))[1])
								if (nchar(zipCode) == 4) zipCode = paste0("0",zipCode)
								index = which(zipCodes@data[,"ZCTA5CE10"]==zipCode); maxArea = 0; pol = NULL
								if (length(index) == 1)
									{
										for (j in 1:length(zipCodes@polygons[[index]]@Polygons))
											{
												if (maxArea < zipCodes@polygons[[index]]@Polygons[[j]]@area)
													{
														maxArea = zipCodes@polygons[[index]]@Polygons[[j]]@area
														pol = zipCodes@polygons[[index]]@Polygons[[j]]
													}
											}
										coords = spsample(pol, 1, type="random")@coords
										samplingData[i,"longitude"] = coords[,"x"]
										samplingData[i,"latitude"] = coords[,"y"]
									}	else	{
										cat("Unmatched zip code: ",zipCode," (for ",variants[h],")\n",sep="")
									}
							}
					}
				write.csv(samplingData, paste0("Preliminary_discrete_runs/",variants[h],"_data.csv"), quote=F, row.names=F)
			}	
		samplingData = read.csv(paste0("Preliminary_discrete_runs/",variants[h],"_data.csv"), head=T)
		for (i in 1:length(NYC_introductions))
			{
				tab = c()
				if (i == 1)
					{
						clusters2 = list(); centroids = list()
					}
				for (j in 1:length(clusters1[[i]]))
					{
						index = which(samplingData[,"sequence_ID"]==clusters1[[i]][j])
						if (length(index) == 1)
							{
								collection_date = samplingData[index,"collection_date"]
								lineage = samplingData[index,"lineage"]
								state = samplingData[index,"state"]
								county = samplingData[index,"county"]
								location = gsub(" ","",samplingData[index,"location"])
								zip_code = samplingData[index,"zip_code"]
								latitude = samplingData[index,"latitude"]
								longitude = samplingData[index,"longitude"]
								line = cbind(collection_date, lineage, state, county, location, zip_code, latitude, longitude)
								row.names(line) = clusters1[[i]][j]; tab = rbind(tab, line)
							}
					}
				colnames(tab) = c("collection_date","lineage","state","county","location","zip_code","latitude","longitude")
				clusters2[[i]] = tab
			}
		clusters1_list[[h]] = clusters1; clusters2_list[[h]] = clusters2
	}
saveRDS(clusters1_list, "Step5_Clusters1_list.rds"); saveRDS(clusters2_list, "Step5_Clusters2_list.rds")
clusters1_list = readRDS("Step5_Clusters1_list.rds"); clusters2_list = readRDS("Step5_Clusters2_list.rds")
cluster_sizes_list = list()
for (i in 1:length(clusters2_list))
	{
		cluster_sizes = c()
		for (j in 1:length(clusters2_list[[i]])) cluster_sizes = c(cluster_sizes, dim(clusters2_list[[i]][[j]])[1])
		cluster_sizes_list[[i]] = cluster_sizes
	}
for (i in 1:length(clusters2_list))
	{
		if (i == 1) cat("\tNumber of distinct clusters:\n",sep="")
		cat("\t\t",variants[i],":\t",length(cluster_sizes_list[[i]]),"\n",sep="")
	}	# 287 (Iota), 786 (Alpha), 3246 (Delta), 2364 (O-BA1)
for (i in 1:length(clusters2_list))
	{
		if (i == 1) cat("\tProportion of clusters of size = 1:\n",sep="")
		cat("\t\t",variants[i],":\t",round(sum(cluster_sizes_list[[i]]==1)/length(cluster_sizes_list[[i]]),3),"\n",sep="")
	}	# 0.74 (Iota), 0.70 (Alpha), 0.77 (Delta), 0.84 (O-BA1)
source("DTA_tree_extraction2.r")
for (h in 1:length(variants))
	{
		metadata = read.table(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.txt"), head=T)
		mostRecentSamplingDatum = max(decimal_date(ymd(metadata[,"collection_date"])), na.rm=T)
		if (file.exists(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.rds")))
			{
				tree = readRDS(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.rds"))
			}	else	{
				tree = readAnnotatedNexus(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.tree"))
			}
		tab = DTA_tree_extraction2(tree, mostRecentSamplingDatum) # alternative DTA tree extraction to avoid memory issue
		write.csv(tab, paste0("Preliminary_discrete_runs/",variants[h],"_DTA.csv"), row.names=F, quote=F)
	}
NYC_introduction_dates_list = list()
for (h in 1:length(variants))
	{
		tab = read.csv(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.csv"), head=T)
		indices = which((tab[,"startLoc"]=="other")&(tab[,"endLoc"]=="NYC_counties"))
		# NYC_introduction_dates_list[[h]] = (tab[indices,"endYear"]+tab[indices,"startYear"])/2
		NYC_introduction_dates_list[[h]] = tab[indices,"endYear"]
	}

# 6. Preparing the discrete phylogeographic analyses among counties

tipSwapping = TRUE; tipSwapping = FALSE
for (h in 1:length(variants))
	{
		if (file.exists(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.tree")))
			{
				tree = readAnnotatedNexus(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.tree"))
			}	else	{
				tree = readAnnotatedNexus(paste0("Preliminary_discrete_runs/",variants[h],"_last.tree"))
			}
		clusters2 = clusters2_list[[h]]; NYC_introductions = NYC_introductions_list[[h]]
		template = scan("Template_for_BBSVS.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
		if (tipSwapping == FALSE)
			{
				dir.create(file.path(paste0("DTA_boroughs_analyses/",variants[h],"_DTA")), showWarnings=F)
				sink(file=paste0("DTA_boroughs_analyses/",variants[h],"_DTA/All_clades.xml"))
			}	else	{
				dir.create(file.path(paste0("DTA_boroughs_analyses/",variants[h],"_TSW")), showWarnings=F)
				sink(file=paste0("DTA_boroughs_analyses/",variants[h],"_TSW/All_clades.xml"))				
			}
		for (i in 1:length(template))
			{
				if ((grepl("</operators>",template[i]))&(tipSwapping == TRUE))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t\t<tipStateSwapOperator weight=\"2\" uniformRandomization=\"true\">","\n"))
										cat(paste0("\t\t\t<ancestralTreeLikelihood idref=\"location.treeLikelihood_",j,"\"/>","\n"))
										cat(paste0("\t\t</tipStateSwapOperator>","\n"))
									}
							}
					}
				cat(template[i],"\n")
				if (grepl("Insert taxa blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t<taxa id=\"taxa_",j,"\">","\n"))
										for (k in 1:dim(clusters2[[j]])[1])
											{
												if (!is.na(clusters2[[j]][k,"location"]))
													{
														cat(paste0("\t\t<taxon id=\"",row.names(clusters2[[j]])[k],"\">","\n"))
														cat(paste0("\t\t\t<date value=\"",clusters2[[j]][k,"collection_date"],"\" direction=\"forwards\" units=\"years\"/>","\n"))
														cat("\t\t\t<attr name=\"location\">\n")
														cat(paste0("\t\t\t\t",clusters2[[j]][k,"location"],"\n"))
														cat("\t\t\t</attr>\n")
														cat("\t\t</taxon>\n")
													}
											}
										cat("\t</taxa>","\n")
									}
							}
					}
				if (grepl("Insert alignment blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t<alignment id=\"alignment_",j,"\" dataType=\"nucleotide\">","\n"))
										for (k in 1:dim(clusters2[[j]])[1])
											{
												if (!is.na(clusters2[[j]][k,"location"]))
													{
														cat("\t\t<sequence>\n")
														cat(paste0("\t\t\t<taxon idref=\"",row.names(clusters2[[j]])[k],"\"/>","\n"))
														cat("\t\t\tNNNN\n")
														cat("\t\t</sequence>\n")
													}
											}
										cat("\t</alignment>","\n")
									}
							}
					}
				if (grepl("Insert pattern blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t<patterns id=\"patterns_",j,"\" from=\"1\" strip=\"false\">","\n"))
										cat(paste0("\t\t<alignment idref=\"alignment_",j,"\"/>","\n"))
										cat("\t</patterns>","\n")
									}
							}
					}
				if (grepl("Insert starting tree blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										tre = tree_subset(tree, tree$edge[NYC_introductions[j],2], levels_back=0)
										tips = row.names(clusters2[[j]]); tips = tips[which(!is.na(clusters2[[j]][,"location"]))]
										tips_to_drop = tre$tip.label[which(!gsub("'","",tre$tip.label)%in%tips)]
										if (length(tips_to_drop) > 0) tre = ape::drop.tip(tre, tips_to_drop)
										if (tipSwapping == FALSE)
											{
												write.tree(tre, paste0("DTA_boroughs_analyses/",variants[h],"_DTA/Clade_",j,".tre"))
												tre = scan(paste0("DTA_boroughs_analyses/",variants[h],"_DTA/Clade_",j,".tre"), what="", sep="\n", quiet=T)
												txt = c("#NEXUS","begin trees;",paste0("\ttree tree_1 = [&R] ",tre),"end;")
												write(txt, paste0("DTA_boroughs_analyses/",variants[h],"_DTA/Clade_",j,".tre"))
											}	else	{
												write.tree(tre, paste0("DTA_boroughs_analyses/",variants[h],"_TSW/Clade_",j,".tre"))
												tre = scan(paste0("DTA_boroughs_analyses/",variants[h],"_TSW/Clade_",j,".tre"), what="", sep="\n", quiet=T)
												txt = c("#NEXUS","begin trees;",paste0("\ttree tree_1 = [&R] ",tre),"end;")
												write(txt, paste0("DTA_boroughs_analyses/",variants[h],"_TSW/Clade_",j,".tre"))												
											}
										cat(paste0("\t<empiricalTreeDistributionModel id=\"treeModel_",j,"\" fileName=\"Clade_",j,".tre\">","\n"))
										cat(paste0("\t\t<taxa idref=\"taxa_",j,"\"/>","\n"))
										cat("\t</empiricalTreeDistributionModel>","\n")
									}
							}
					}
				if (grepl("Insert location.pattern blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t<attributePatterns id=\"location.pattern_",j,"\" attribute=\"location\">","\n"))
										cat(paste0("\t\t<taxa idref=\"taxa_",j,"\"/>","\n"))
										cat(paste0("\t\t<generalDataType idref=\"location.dataType\"/>","\n"))
										cat("\t</attributePatterns>","\n")
									}
							}
					}
				if (grepl("Insert rateStatistic blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t<rateStatistic id=\"location.meanRate_",j,"\" name=\"location.meanRate_",j,"\" mode=\"mean\" internal=\"true\" external=\"true\">","\n"))
										cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
										cat(paste0("\t\t<strictClockBranchRates idref=\"location.branchRates\"/>","\n"))
										cat("\t</rateStatistic>","\n")
									}
							}
					}
				if (grepl("Insert ancestralTreeLikelihood blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t<ancestralTreeLikelihood id=\"location.treeLikelihood_",j,"\" stateTagName=\"location.states\" useUniformization=\"true\" saveCompleteHistory=\"false\" logCompleteHistory=\"false\">","\n"))
										cat(paste0("\t\t<attributePatterns idref=\"location.pattern_",j,"\"/>","\n"))
										cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
										cat(paste0("\t\t<siteModel idref=\"location.siteModel\"/>","\n"))
										cat(paste0("\t\t<generalSubstitutionModel idref=\"location.model\"/>","\n"))
										cat(paste0("\t\t<strictClockBranchRates idref=\"location.branchRates\"/>","\n"))
										cat(paste0("\t\t<frequencyModel id=\"location.root.frequencyModel_",j,"\" normalize=\"true\">","\n"))
										cat(paste0("\t\t\t<generalDataType idref=\"location.dataType\"/>","\n"))
										cat(paste0("\t\t\t<frequencies>","\n"))
										cat(paste0("\t\t\t\t<parameter id=\"location.root.frequencies_",j,"\" dimension=\"",length(NYC_counties),"\"/>","\n"))
										cat(paste0("\t\t\t</frequencies>","\n"))
										cat(paste0("\t\t</frequencyModel>","\n"))
										cat(paste0("\t</ancestralTreeLikelihood>","\n"))
									}
							}
					}
				if (grepl("Insert deltaExchange blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t\t<deltaExchange delta=\"0.75\" weight=\"1\">","\n"))
										cat(paste0("\t\t\t<parameter idref=\"location.root.frequencies_",j,"\"/>","\n"))
										cat("\t\t</deltaExchange>","\n")
									}
							}
					}
				if (grepl("Insert uniformPrior blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t\t\t\t<ctmcScalePrior>","\n"))
										cat(paste0("\t\t\t\t\t<ctmcScale>","\n"))
										cat(paste0("\t\t\t\t\t\t<parameter idref=\"location.clock.rate\"/>","\n"))
										cat(paste0("\t\t\t\t\t</ctmcScale>","\n"))
										cat(paste0("\t\t\t\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
										cat(paste0("\t\t\t\t</ctmcScalePrior>","\n"))
									}
							}
					}
				if (grepl("Insert uniformPrior blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1.0\">","\n"))
										cat(paste0("\t\t\t\t\t<parameter idref=\"location.root.frequencies_",j,"\"/>","\n"))
										cat(paste0("\t\t\t\t</uniformPrior>","\n"))
									}
							}
					}
				if (grepl("Insert ancestralTreeLikelihood lines 1",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t\t\t\t<ancestralTreeLikelihood idref=\"location.treeLikelihood_",j,"\"/>","\n"))
									}
							}
					}
				if (grepl("Insert rateStatistic lines",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t\t\t<rateStatistic idref=\"location.meanRate_",j,"\"/>","\n"))
									}
							}
					}
				if (grepl("Insert ancestralTreeLikelihood lines 2",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t\t\t<ancestralTreeLikelihood idref=\"location.treeLikelihood_",j,"\"/>","\n"))
									}
							}
					}
				if (grepl("Insert treeFileLog blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3))
									{
										cat(paste0("\t\t<logTree id=\"treeFileLog_",j,"\" logEvery=\"50000\" nexusFormat=\"true\" fileName=\"Clade_",j,".trees\" sortTranslationTable=\"true\">","\n"))
										cat(paste0("\t\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
										cat(paste0("\t\t\t<trait name=\"rate\" tag=\"location.rate\">","\n"))
										cat(paste0("\t\t\t\t<strictClockBranchRates idref=\"location.branchRates\"/>","\n"))
										cat(paste0("\t\t\t</trait>","\n"))
										cat(paste0("\t\t\t<joint idref=\"joint\"/>","\n"))
										cat(paste0("\t\t\t<trait name=\"location.states\" tag=\"location\">","\n"))
										cat(paste0("\t\t\t\t<ancestralTreeLikelihood idref=\"location.treeLikelihood_",j,"\"/>","\n"))
										cat(paste0("\t\t\t</trait>","\n"))
										cat(paste0("\t\t</logTree>","\n"))
									}
							}
					}
			}
		sink(NULL)
	}

# 7. Preparing the continuous phylogeographic analyses (RRW, Cauchy model)

for (h in 1:length(variants))
	{
		if (file.exists(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.tree")))
			{
				tree = readAnnotatedNexus(paste0("Preliminary_discrete_runs/",variants[h],"_DTA.tree")); indices = c()
			}	else	{
				tree = readAnnotatedNexus(paste0("Preliminary_discrete_runs/",variants[h],"_last.tree")); indices = c()
			}
		clusters2 = clusters2_list[[h]]; NYC_introductions = NYC_introductions_list[[h]]
		template = scan("RRW_template_file2.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
		dir.create(file.path(paste0("RRW_diffusion_analyses/",variants[h],"_RRW")), showWarnings=F)
		sink(file=paste0("RRW_diffusion_analyses/",variants[h],"_RRW/All_clades.xml"))
		for (i in 1:length(template))
			{
				cat(template[i],"\n")
				if (grepl("Insert taxa blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<taxa id=\"taxa_",j,"\">","\n"))
										for (k in 1:dim(clusters2[[j]])[1])
											{
												if (!is.na(clusters2[[j]][k,"longitude"]))
													{
														cat(paste0("\t\t<taxon id=\"",row.names(clusters2[[j]])[k],"\">","\n"))
														cat(paste0("\t\t\t<date value=\"",clusters2[[j]][k,"collection_date"],"\" direction=\"forwards\" units=\"years\"/>","\n"))
														cat("\t\t\t<attr name=\"latitude\">\n")
														cat(paste0("\t\t\t\t",clusters2[[j]][k,"latitude"],"\n"))
														cat("\t\t\t</attr>\n")
														cat("\t\t\t<attr name=\"longitude\">\n")
														cat(paste0("\t\t\t\t",clusters2[[j]][k,"longitude"],"\n"))
														cat("\t\t\t</attr>\n")
														cat("\t\t\t<attr name=\"coordinates\">\n")
														cat(paste0("\t\t\t\t",clusters2[[j]][k,"latitude"]," ",clusters2[[j]][k,"longitude"],"\n"))
														cat("\t\t\t</attr>\n")
														cat("\t\t</taxon>\n")
													}
											}
										cat("\t</taxa>","\n")
									}
							}
					}
				if (grepl("Insert alignment blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<alignment id=\"alignment_",j,"\" dataType=\"nucleotide\">","\n"))
										for (k in 1:dim(clusters2[[j]])[1])
											{
												if (!is.na(clusters2[[j]][k,"longitude"]))
													{
														cat("\t\t<sequence>\n")
														cat(paste0("\t\t\t<taxon idref=\"",row.names(clusters2[[j]])[k],"\"/>","\n"))
														cat("\t\t\tNNNN\n")
														cat("\t\t</sequence>\n")
													}
											}
										cat("\t</alignment>","\n")
									}
							}
					}
				if (grepl("Insert pattern blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<patterns id=\"patterns_",j,"\" from=\"1\" strip=\"false\">","\n"))
										cat(paste0("\t\t<alignment idref=\"alignment_",j,"\"/>","\n"))
										cat("\t</patterns>","\n")
									}
							}
					}
				if (grepl("Insert starting tree blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										tre = tree_subset(tree, tree$edge[NYC_introductions[j],2], levels_back=0)
										tips = row.names(clusters2[[j]]); tips = tips[which(!is.na(clusters2[[j]][,"longitude"]))]
										tips_to_drop = tre$tip.label[which(!gsub("'","",tre$tip.label)%in%tips)]
										if (length(tips_to_drop) > 0) tre = ape::drop.tip(tre, tips_to_drop)
										write.tree(tre, paste0("RRW_diffusion_analyses/",variants[h],"_RRW//Clade_",j,".tre"))
										tre = scan(paste0("RRW_diffusion_analyses/",variants[h],"_RRW//Clade_",j,".tre"), what="", sep="\n", quiet=T)
										txt = c("#NEXUS","begin trees;",paste0("\ttree tree_1 = [&R] ",tre),"end;")
										write(txt, paste0("RRW_diffusion_analyses/",variants[h],"_RRW//Clade_",j,".tre"))
										cat(paste0("\t<empiricalTreeDistributionModel id=\"treeModel_",j,"\" fileName=\"Clade_",j,".tre\">","\n"))
										cat(paste0("\t\t<taxa idref=\"taxa_",j,"\"/>","\n"))
										cat("\t</empiricalTreeDistributionModel>","\n")
									}
							}
					}
				if (grepl("Insert tree model blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<treeModel id=\"treeModel_",j,"\">","\n"))
										cat(paste0("\t\t<coalescentTree idref=\"startingTree_",j,"\"/>","\n"))
										cat("\t\t<rootHeight>","\n")
										cat(paste0("\t\t\t<parameter id=\"treeModel.rootHeight_",j,"\"/>","\n"))
										cat("\t\t</rootHeight>","\n")
										cat("\t\t<nodeHeights internalNodes=\"true\">","\n")
										cat(paste0("\t\t\t<parameter id=\"treeModel.internalNodeHeights_",j,"\"/>","\n"))
										cat("\t\t</nodeHeights>","\n")
										cat("\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">","\n")
										cat(paste0("\t\t\t<parameter id=\"treeModel.allInternalNodeHeights_",j,"\"/>","\n"))
										cat("\t\t</nodeHeights>","\n")
										cat("\t</treeModel>","\n")
									}
							}
					}
				if (grepl("Insert arbitraryBranchRates blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<arbitraryBranchRates id=\"coordinates.diffusion.branchRates",j,"\">","\n"))
										cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
										cat("\t\t<rates>","\n")
										cat(paste0("\t\t\t<parameter id=\"coordinates.diffusion.rates",j,"\" lower=\"0.0\"/>","\n"))
										cat("\t\t</rates>","\n")
										cat("\t</arbitraryBranchRates>","\n")
									}
							}
					}
				if (grepl("Insert distributionLikelihood blocks 1",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<distributionLikelihood id=\"coordinates.diffusion.prior",j,"\">","\n"))
										cat("\t\t<data>","\n")
										cat(paste0("\t\t\t<parameter idref=\"coordinates.diffusion.rates",j,"\"/>","\n"))
										cat("\t\t</data>","\n")
										cat("\t\t<distribution>","\n")
										cat(paste0("\t\t\t<onePGammaDistributionModel>","\n"))
										cat("\t\t\t\t<shape>","\n")
										cat("\t\t\t\t\t<parameter value=\"0.5\"/>","\n")
										cat("\t\t\t\t</shape>","\n")
										cat("\t\t\t</onePGammaDistributionModel>","\n")
										cat("\t\t</distribution>","\n")
										cat("\t</distributionLikelihood>","\n")
									}
							}
					}
				if (grepl("Insert coordinates.traitLikelihood blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<multivariateTraitLikelihood id=\"coordinates.traitLikelihood",j,"\" traitName=\"coordinates\" useTreeLength=\"true\" scaleByTime=\"true\" reportAsMultivariate=\"true\" reciprocalRates=\"true\" integrateInternalTraits=\"true\">","\n"))
										cat("\t\t<multivariateDiffusionModel idref=\"coordinates.diffusionModel\"/>","\n")
										cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>"))
										cat("\t\t<traitParameter>","\n")
										cat(paste0("\t\t\t<parameter id=\"leaf.coordinates",j,"\"/>","\n"))
										cat("\t\t</traitParameter>","\n")
										cat("\t\t<conjugateRootPrior>","\n")
										cat("\t\t\t<meanParameter>","\n")
										cat("\t\t\t\t<parameter value=\"0.0 0.0\"/>","\n")
										cat("\t\t\t</meanParameter>","\n")
										cat("\t\t\t<priorSampleSize>","\n")
										cat("\t\t\t\t<parameter value=\"0.000001\"/>","\n")
										cat("\t\t\t</priorSampleSize>","\n")
										cat("\t\t</conjugateRootPrior>","\n")
										cat(paste0("\t\t<arbitraryBranchRates idref=\"coordinates.diffusion.branchRates",j,"\"/>","\n"))
										cat("\t</multivariateTraitLikelihood>","\n")
									}
							}
					}
				if (grepl("Insert continuousDiffusionStatistic blocks 1",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t<continuousDiffusionStatistic id=\"coordinates.diffusionRate",j,"\" greatCircleDistance=\"true\">","\n"))
										cat(paste0("\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
										cat("\t</continuousDiffusionStatistic>","\n")
									}
							}
					}
				if (grepl("Insert scaleOperator blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"30\">","\n"))
										cat(paste0("\t\t\t<parameter idref=\"coordinates.diffusion.rates",j,"\"/>","\n"))
										cat("\t\t</scaleOperator>","\n")
									}
							}
					}
				if (grepl("Insert precisionGibbsOperator blocks",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t<precisionGibbsOperator weight=\"2\">","\n"))
										cat(paste0("\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
										cat("\t\t\t<multivariateWishartPrior idref=\"coordinates.precisionPrior\"/>","\n")
										cat("\t\t</precisionGibbsOperator>","\n")
									}
							}
					}
				if (grepl("Insert distributionLikelihood blocks 2",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t\t\t<distributionLikelihood idref=\"coordinates.diffusion.prior",j,"\"/>","\n"))
									}
							}
					}
				if (grepl("Insert multivariateTraitLikelihood blocks 1",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
									}
							}
					}
				if (grepl("Insert continuousDiffusionStatistic blocks 2",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t\t\t<continuousDiffusionStatistic idref=\"coordinates.diffusionRate",j,"\"/>","\n"))
									}
							}
					}
				if (grepl("Insert multivariateTraitLikelihood blocks 2",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
									}
							}
					}
				if (grepl("<!-- Insert logTree blocks -->",template[i]))
					{
						for (j in 1:length(clusters2))
							{
								if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
									{
										cat(paste0("\t\t<logTree id=\"treeFileLog",j,"\" logEvery=\"100000\" nexusFormat=\"true\" fileName=\"Clade_",j,".trees\" sortTranslationTable=\"true\">","\n"))
										cat(paste0("\t\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
										cat("\t\t\t<joint idref=\"joint\"/>","\n")
										cat("\t\t\t<trait name=\"coordinates\" tag=\"coordinates\">","\n")
										cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
										cat("\t\t\t</trait>","\n")
										cat("\t\t\t<multivariateDiffusionModel idref=\"coordinates.diffusionModel\"/>","\n")
										cat("\t\t\t<trait name=\"rate\" tag=\"coordinates.rate\">","\n")
										cat(paste0("\t\t\t\t<arbitraryBranchRates idref=\"coordinates.diffusion.branchRates",j,"\"/>","\n"))
										cat("\t\t\t</trait>","\n")
										cat("\t\t</logTree>","\n")
									}
							}
					}
			}
		sink(NULL)
	}
	
# 8. Building the maximum clade consensus (MCC) trees for both analyses

wd = getwd()
for (h in 1:length(variants))
	{
		setwd(paste0(wd,"/DTA_boroughs_analyses/",variants[h],"_DTA/"))
		treeFiles = list.files(); treeFiles = gsub(".trees","",treeFiles[which(grepl(".trees",treeFiles))])
		for (j in 1:length(treeFiles))
			{
				system(paste0("BEAST_1104/bin/treeannotator -burninTrees 101 -heights keep ",treeFiles[j],".trees ",treeFiles[j],".tree"), ignore.stdout=F, ignore.stderr=F)
			}
		setwd(paste0(wd,"/RRW_diffusion_analyses/",variants[h],"_RRW/"))
		treeFiles = list.files(); treeFiles = gsub(".trees","",treeFiles[which(grepl(".trees",treeFiles))])
		for (j in 1:length(treeFiles))
			{
				system(paste0("BEAST_1104/bin/treeannotator -burninTrees 101 -heights keep ",treeFiles[j],".trees ",treeFiles[j],".tree"), ignore.stdout=F, ignore.stderr=F)
			}
	}
setwd(wd)

# 9. Extracting spatio-temporal information embedded in MCC and posterior trees

source("DTA_tree_extraction1.r"); source("MCC_tree_extractions.r")
wd = getwd(); previousVersion = F; registerDoMC(cores=5)
for (h in 1:length(variants))
	{
		clusters2 = clusters2_list[[h]]
		setwd(paste0(wd,"/DTA_boroughs_analyses/",variants[h],"_DTA/"))
		treeFiles = list.files(); treeFiles = treeFiles[which(grepl(".trees",treeFiles))]
		for (i in 1:length(treeFiles))
			{
				dir.create(file.path(gsub(".trees","_ext",treeFiles[i])), showWarnings=F)
				trees = readAnnotatedNexus(treeFiles[i]); trees = trees[102:1001]
				# buffer = foreach(j = 1:length(trees)) %dopar% {
				for (j in 1:length(trees)) {
						tree = trees[[j]]
						if (previousVersion)
							{
								tab = matrix(nrow=dim(tree$edge)[1], ncol=4)
								colnames(tab) = c("node1","node2","startLoc","endLoc")
								tab[,"node1"] = tree$edge[,1]; tab[,"node2"] = tree$edge[,2]
								for (k in 1:dim(tree$edge)[1])	
									{
										tab[k,"endLoc"] = tree$annotations[[k]]$location
										index = which(tree$edge[,2]==tree$edge[k,1])
										if (length(index) == 1)
											{
												tab[k,"startLoc"] = tree$annotations[[index]]$location
											}	else		{
												if (!tree$edge[k,1]%in%tree$edge[,2])
													{
														tab[k,"startLoc"] = tree$root.annotation$location
													}
											}
									}
							}	else	{
								index = as.numeric(unlist(strsplit(gsub(".trees","",treeFiles[i]),"_"))[2])
								mostRecentSamplingDatum = max(as.numeric(clusters2[[index]][which(!is.na(clusters2[[index]][,"location"])),"collection_date"]))
								tab = DTA_tree_extraction1(tree, mostRecentSamplingDatum)
								tab$cladeID = rep(index, dim(tab)[1])
							}
						write.csv(tab, paste0(gsub(".trees","_ext",treeFiles[i]),"/TreeExtractions_",j,".csv"), row.names=F, quote=F)
						# j
					}
			}
		dir.create(file.path("All_clades_ext"), showWarnings=F); tab = NULL
		for (i in 1:nberOfExtractionFiles)
			{
				for (j in 1:length(treeFiles))
					{
						if (j == 1)
							{
								tab = read.csv(paste0(gsub(".trees","_ext",treeFiles[j]),"/TreeExtractions_",i,".csv"))
							}	else	{
								tab = rbind(tab, read.csv(paste0(gsub(".trees","_ext",treeFiles[j]),"/TreeExtractions_",i,".csv")))
							}
					}
				write.csv(tab, paste0("All_clades_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
			}
		matrices = list()
		for (i in 1:nberOfExtractionFiles)
			{
				mat = matrix(0, nrow=length(NYC_counties), ncol=length(NYC_counties))
				row.names(mat) = NYC_counties; colnames(mat) = NYC_counties
				tab = read.csv(paste0("All_clades_ext/TreeExtractions_",i,".csv"), head=T)
				for (j in 1:dim(tab)[1])
					{
						index1 = which(NYC_counties==tab[j,"startLoc"])
						index2 = which(NYC_counties==tab[j,"endLoc"])
						mat[index1,index2] = mat[index1,index2]+1
					}
				matrices[[i]] = mat
			}
		saveRDS(matrices, "Matrices.rds")
		log1 = scan(paste0("All_clades1.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		write(log1[which(!grepl("# ",log1))], paste0("All_clades3.log"))
		log1 = read.table(paste0("All_clades3.log"), header=T, sep="\t"); log1 = log1[102:1001,]
		setwd(paste0(wd,"/DTA_boroughs_analyses/",variants[h],"_TSW/"))
		log2 = scan(paste0("All_clades1.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		write(log2[which(!grepl("# ",log2))], paste0("All_clades3.log"))
		log2 = read.table(paste0("All_clades3.log"), header=T, sep="\t"); log2 = log2[102:1001,]
		BFs1 = matrix(nrow=length(NYC_counties), ncol=length(NYC_counties))
		BFs2 = matrix(nrow=length(NYC_counties), ncol=length(NYC_counties))
		row.names(BFs1) = NYC_counties; colnames(BFs1) = NYC_counties
		row.names(BFs2) = NYC_counties; colnames(BFs2) = NYC_counties
		for (i in 1:length(NYC_counties))
			{
				for (j in 1:length(NYC_counties))
					{
						if (i != j)
							{
								colName = paste0("location.indicators.",gsub(" ",".",NYC_counties[i]),".",gsub(" ",".",NYC_counties[j]))
								index1 = which(colnames(log1)==colName); index2 = which(colnames(log2)==colName)
								p = sum(log1[,index1]==1)/dim(log1)[1]
								K = 56 # length(locations)*(length(locations)-1) # K shoulf be divided by 2 if "symetric" case
								q = (log(2)+K-1)/(K*(K-1))
								BFs1[i,j] = (p/(1-p))/(q/(1-q))
								p1 = sum(log1[,index1]==1)/dim(log1)[1]
								p2 = sum(log2[,index2]==1)/dim(log2)[1]
								BFs2[i,j] = (p1/(1-p1))/(p2/(1-p2))
							}
					}
			}
		setwd(paste0(wd,"/DTA_boroughs_analyses/",variants[h],"_DTA/"))
		write.table(round(BFs1,1), paste0("BF_values.csv"), sep=",", quote=F)
		setwd(paste0(wd,"/DTA_boroughs_analyses/",variants[h],"_TSW/"))
		write.table(round(BFs2,1), paste0("BF_values.csv"), sep=",", quote=F)
			setwd(paste0(wd,"/RRW_diffusion_analyses/",variants[h],"_RRW/"))
		treeFiles = list.files(); treeFiles = gsub(".trees","",treeFiles[which(grepl(".trees",treeFiles))])
		for (i in 1:length(treeFiles))
			{
				index = as.numeric(unlist(strsplit(treeFiles[i],"_"))[2])
				mostRecentSamplingDatum = max(as.numeric(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collection_date"]))
				mcc_tre = readAnnotatedNexus(paste0(treeFiles[i],".tree")); dates = c()
				mcc_tab = MCC_tree_extractions(mcc_tre, mostRecentSamplingDatum)
				write.csv(mcc_tab, paste0(treeFiles[i],".csv"), row.names=F, quote=F)
			}
		nberOfTreesToSample = nberOfExtractionFiles; burnIn = 101; randomSampling = FALSE; coordinateAttributeName = "coordinates"; nberOfCores = 5
		for (i in 1:length(treeFiles))
			{
				localTreesDirectory = paste0(treeFiles[i],"_ext")
				index = as.numeric(unlist(strsplit(treeFiles[i],"_"))[2])
				mostRecentSamplingDatum = max(as.numeric(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collection_date"]))
				allTrees = scan(file=paste0(treeFiles[i],".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
			}
		for (i in 1:length(treeFiles))
			{
				tab = read.csv(paste0(treeFiles[i],".csv"), head=T)
				if (i == 1)
					{
						all = tab
					}	else	{
						maxNodeID = max(all[,c("node1","node2")])
						tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
						all = rbind(all, tab)
					}
			}
		write.csv(all, "All_clades.csv", row.names=F, quote=F)
		dir.create(file.path("All_clades_ext"), showWarnings=F)
		nberOfExtractionFiles = nberOfTreesToSample
		for (i in 1:nberOfExtractionFiles)
			{
				if (!file.exists(paste0("All_clades_ext/TreeExtractions_",i,".csv")))
					{
						for (j in 1:length(treeFiles))
							{
								tab = read.csv(paste0(treeFiles[j],"_ext/TreeExtractions_",i,".csv"), head=T)
								if (j == 1)
									{
										all = tab
									}	else	{
										maxNodeID = max(all[,c("node1","node2")])
										tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
										all = rbind(all, tab)
									}
							}
						write.csv(all, paste0("All_clades_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
					}
			}
	}
setwd(wd)

# 10. Visualisations of the discrete and continuous phylogeographic reconstructions

variants = c("Iota","Alpha","Delta","O-BA1")
variant_names = c("Iota","Alpha","Delta","Omicron (BA1)")
all_counties = shapefile("NY_state_all_shapefiles/GADM_USA_2.shp")
NY_state_counties = subset(all_counties, all_counties@data$NAME_1=="New York")
selected_counties = subset(NY_state_counties, gsub(" ","",NY_state_counties@data$NAME_2)%in%NYC_counties)
centroids = coordinates(selected_counties); row.names(centroids) = gsub(" ","",selected_counties@data$NAME_2)
centroids = centroids[NYC_counties,]; matrix_mean_list = list()
for (h in 1:length(variants))
	{
		if (h == 1) nberOfExtractionFiles = 400 # TO BE REMOVED !!
		if (h != 1) nberOfExtractionFiles = 900 # TO BE REMOVED !!
		matrices = readRDS(paste0("DTA_boroughs_analyses/",variants[h],"_DTA/Matrices.rds"))
		matrix_mean = matrix(0, nrow=length(NYC_counties), ncol=length(NYC_counties))
		for (i in 1:nberOfExtractionFiles) matrix_mean = matrix_mean+matrices[[i]]
		matrix_mean = matrix_mean/nberOfExtractionFiles; matrix_mean_list[[h]] = matrix_mean
	}
minVals1 = min(diag(matrix_mean_list[[1]])); maxVals1 = max(diag(matrix_mean_list[[1]]))
mat = matrix_mean_list[[1]]; diag(mat) = NA; minVals2 = min(mat, na.rm=T); maxVals2 = max(mat, na.rm=T)
if (length(variants) > 1)
	{
		for (h in 2:length(variants))
			{
				mat1 = matrix_mean_list[[h]]; mat2 = mat1; diag(mat2) = NA
				if (minVals1 > min(diag(mat1))) minVals1 = min(diag(mat1))
				if (maxVals1 < max(diag(mat1))) maxVals1 = max(diag(mat1))
				if (minVals2 > min(mat2,na.rm=T)) minVals2 = min(mat2,na.rm=T)
				if (maxVals2 < max(mat2,na.rm=T)) maxVals2 = max(mat2,na.rm=T)
			}
	}
mccs = list(); polygons_list = list()
for (h in 1:length(variants))
	{
		if (h == 1) nberOfExtractionFiles = 400 # TO BE REMOVED !!
		if (h != 1) nberOfExtractionFiles = 900 # TO BE REMOVED !!
		mccs[[h]] = read.csv(paste0("RRW_diffusion_analyses/",variants[h],"_RRW/All_clades.csv"), head=T)
		localTreesDirectory = paste0("RRW_diffusion_analyses/",variants[h],"_RRW/All_clades_ext")
		percentage = 80; prob = percentage/100; precision = 1/(365/7); startDatum = min(mccs[[h]][,"startYear"])
		polygons_list[[h]] = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
	}
minYear = min(mccs[[1]][,"startYear"]); maxYear = max(mccs[[1]][,"endYear"])
if (length(variants) > 1)
	{
		for (h in 2:length(variants))
			{
				if (minYear > min(mccs[[h]][,"startYear"])) minYear = min(mccs[[h]][,"startYear"])
				if (maxYear < max(mccs[[h]][,"endYear"])) maxYear = max(mccs[[h]][,"endYear"])
			}
	}
colourScale = rev(colorRampPalette(brewer.pal(11,"PuOr"))(141)[26:126])

pdf(paste0("Figure_A1_NEW.pdf"), width=11.0, height=11) # dev.new(width=11.0, height=11)
par(mfrow=c(4,2), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30"); cexNode = 0.7
plottingLegend = TRUE; croppingPolygons = TRUE; adjustedBFs = TRUE; onlyInternalNodesOfTipBranches = FALSE
for (h in 1:length(variants))
	{
		# mat = matrix_mean_list[[h]]; multiplier1 = 500; multiplier2 = 15; multiplier3 = 0.4
		# if (adjustedBFs == TRUE) BFs = read.csv(paste0("DTA_boroughs_analyses/",variants[h],"_TSW/BF_values.csv"), head=T)
		# if (adjustedBFs == FALSE) BFs = read.csv(paste0("DTA_boroughs_analyses/",variants[h],"_DTA/BF_values.csv"), head=T)
		# plot(selected_counties, col="gray90", border="gray50", lwd=0.2)
		# points(centroids, cex=sqrt((multiplier1*((diag(mat)-minVals1)/(maxVals1-minVals1)))/pi), pch=16, col="#DE432750")
		# for (i in 1:dim(selected_counties)[1])
			# {
				# for (j in 1:dim(selected_counties)[1])
					# {
						# if ((i!=j)&(mat[i,j]>=1)&(!is.na(BFs[i,j]))&&(BFs[i,j]>3))
							# {
								# LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(mat[i,j]/maxVals2))+0.04
								# curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
											# lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
							# }
					# }
			# }
		# mtext(variant_names[h], side=3, line=-5, cex=0.8)
		# if ((h == 1)&(plottingLegend))
			# {
				# vS = 5; LWD = (((vS-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(vS/maxVals2))+0.04
				# curvedarrow(cbind(-74.20,41.000), cbind(-74.12,41.000), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
							# lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
				# vS = 20; LWD = (((vS-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(vS/maxVals2))+0.04
				# curvedarrow(cbind(-74.20,40.975), cbind(-74.12,40.975), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
							# lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
				# vS = 50; LWD = (((vS-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(vS/maxVals2))+0.04
				# curvedarrow(cbind(-74.20,40.950), cbind(-74.12,40.950), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
							# lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
				# points(cbind(rep(-74.15,4),rep(41.12,4)), cex=sqrt((multiplier1*((c(100,200,500)-minVals1)/(maxVals1-minVals1)))/pi), pch=1, col="#DE432750", lwd=0.3)
			# }
		polygons = polygons_list[[h]]; mcc = mccs[[h]]; selectedBranches = 1:dim(mcc)[1]
		startYears_indices = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		startYears_colours = colourScale[startYears_indices]
		endYears_colours = colourScale[endYears_indices]
		polygons_colours = rep(NA, length(polygons_list[[]]))
		for (i in 1:length(polygons))
			{
				date = as.numeric(names(polygons[[i]])); polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				polygons_colours[i] = paste0(colourScale[polygon_index],"40")
			}
		plot(selected_counties, col="gray90", border=NA, lwd=0.01)
		for (i in length(polygons):1)
			{
				for (j in 1:length(polygons[[i]]@polygons))
					{
						polygons[[i]]@polygons[[j]] = maptools::checkPolygonsHoles(polygons[[i]]@polygons[[j]])
					}
				pol = polygons[[i]]
				if (croppingPolygons == TRUE)
					{
						pol = crop(pol, selected_counties)
					}
				plot(pol, axes=F, col=polygons_colours[i], add=T, border=NA)
			}
		plot(selected_counties, col=NA, border="gray50", lwd=0.2, add=T)
		for (i in selectedBranches)
			{
				curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
			}
		for (i in rev(selectedBranches))
			{
				if (onlyInternalNodesOfTipBranches != TRUE)
					{
						if (!mcc[i,"node1"]%in%mcc[selectedBranches,"node2"])
							{
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[i], cex=cexNode)
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", cex=cexNode, lwd=0.2)
							}
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=cexNode)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", cex=cexNode, lwd=0.2)
					}	else	{
						if (!mcc[i,"node2"]%in%mcc[selectedBranches,"node1"])
							{
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[i], cex=cexNode)
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", cex=cexNode, lwd=0.2)
							}
					}			
			}
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minYear; rast[2] = maxYear
		selectedDates = NYC_introduction_dates_list[[h]]; selectedLabels = rep("", length(selectedDates))	
		plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.47,0.96,0.12,0.14),
			 legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
			 axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.25,0),
			 at=selectedDates, labels=selectedLabels))
		if (h == length(variants))
			{
				selectedDates = decimal_date(ymd(c("2021-01-01","2021-05-01","2021-09-01","2022-01-01")))
				selectedLabels = c("01-01-21","01-04-21","01-09-22","01-01-22")	
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.47,0.96,0.09,0.11),
					 legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
		  			 axis.args=list(cex.axis=0.9, lwd=0, lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.25,0),
		   			 at=selectedDates, labels=selectedLabels))
		     }
	}
dev.off()

pdf(paste0("Figure_A2_NEW.pdf"), width=11.0, height=11) # dev.new(width=11.0, height=11)
par(mfrow=c(4,2), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30"); cexNode = 0.7
plottingLegend = TRUE; croppingPolygons = TRUE; onlyInternalNodesOfTipBranches = FALSE
for (h in 1:length(variants))
	{
		mat = matrix_mean_list[[h]]; multiplier1 = 500; multiplier2 = 15; multiplier3 = 0.4
		BFs = read.csv(paste0("DTA_boroughs_analyses/",variants[h],"_DTA/BF_values.csv"), head=T)
		plot(selected_counties, col="gray90", border="gray50", lwd=0.2)
		points(centroids, cex=sqrt((multiplier1*((diag(mat)-minVals1)/(maxVals1-minVals1)))/pi), pch=16, col="#DE432750")
		for (i in 1:dim(selected_counties)[1])
			{
				for (j in 1:dim(selected_counties)[1])
					{
						if ((i!=j)&(mat[i,j]>=1)&(!is.na(BFs[i,j]))&&(BFs[i,j]>3))
							{
								LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(mat[i,j]/maxVals2))+0.04
								curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
											lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
							}
					}
			}
		mtext(variant_names[h], side=3, line=-5, cex=0.8)
		BFs = read.csv(paste0("DTA_boroughs_analyses/",variants[h],"_TSW/BF_values.csv"), head=T)
		plot(selected_counties, col="gray90", border="gray50", lwd=0.2)
		points(centroids, cex=sqrt((multiplier1*((diag(mat)-minVals1)/(maxVals1-minVals1)))/pi), pch=16, col="#DE432750")
		for (i in 1:dim(selected_counties)[1])
			{
				for (j in 1:dim(selected_counties)[1])
					{
						if ((i!=j)&(mat[i,j]>=1)&(!is.na(BFs[i,j]))&&(BFs[i,j]>3))
							{
								LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(mat[i,j]/maxVals2))+0.04
								curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
											lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
							}
					}
			}
	}
dev.off()

# 11. Dispersal statistics based on the continuous phylogeographic reconstructions

for (h in 1:length(variants))
	{
		localTreesDirectory = paste0("RRW_diffusion_analyses/",variants[h],"_RRW/All_clades_ext")
		timeSlices = 100; onlyTipBranches = FALSE; showingPlots = FALSE; nberOfCores = 1; slidingWindow = 1
		outputName = paste0("RRW_dispersal_statistics/",variants[h])
		spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
	}	
for (h in 1:length(variants))
	{
		stats = read.table(paste0("RRW_dispersal_statistics/",variants[h],"_estimated_dispersal_statistics.txt"), head=T)
		wldv = stats[,"weighted_branch_dispersal_velocity"]/365; wldv_md = round(median(wldv),2); wldv_qs = round(quantile(wldv,c(0.025,0.975)),2)
		stats = read.table(paste0("RRW_dispersal_statistics/",variants[h],"_estimated_dispersal_statistics.txt"), head=T)
		wdc = stats[,"weighted_diffusion_coefficient"]/365; wdc_md = round(median(wdc),1); wdc_qs = round(quantile(wdc,c(0.025,0.975)),1)
		cat(variants[h],":\tWLDV = ",wldv_md," km/day [",wldv_qs[1],"-",wldv_qs[2],"]\t\tWDC = ",wdc_md," km2/year [",wdc_qs[1],"-",wdc_qs[2],"]","\n",sep="")	
	}
		# Iota:		WLDV = 0.25 km/day [0.22-0.34]		WDC = 1.5 km2/year [1.3-2.9]
		# Alpha:	WLDV = 0.50 km/day [0.39-0.87]		WDC = 3.0 km2/year [1.8-9.4]
		# Delta:	WLDV = 1.13 km/day [0.82-2.21]		WDC = 8.8 km2/year [4.9-38.1]
		# O-BA1:	WLDV = 1.71 km/day [1.23-3.32]		WDC = 12.0 km2/year [6.2-55.0]

# 12. Dispersal statistics based on the discrete phylogeographic reconstructions

variants = c("Iota","Alpha","Delta","O-BA1")
variantColours1 = c(); variantColours2 = c()
variantColours1[1] = rgb(150,150,150,255,maxColorValue=255); variantColours2[1] = rgb(150,150,150,100,maxColorValue=255) # light grey (Iota)
variantColours1[2] = rgb(250,165,26,255,maxColorValue=255); variantColours2[2] = rgb(250,165,26,100,maxColorValue=255) # orange (Alpha)
variantColours1[3] = rgb(222,67,39,255,maxColorValue=255); variantColours2[3] = rgb(222,67,39,100,maxColorValue=255) # red (Delta)
variantColours1[4] = rgb(70,118,187,255,maxColorValue=255); variantColours2[4] = rgb(70,118,187,100,maxColorValue=255) # blue (Omicron-BA.1)
variantColours1[5] = rgb(129,81,161,255,maxColorValue=255); variantColours2[5] = rgb(129,81,161,100,maxColorValue=255) # purple (not used)
variantColours1[6] = rgb(139,101,8,255,maxColorValue=255); variantColours2[6] = rgb(139,101,8,100,maxColorValue=255) # brown (not used)
variantColours1[7] = rgb(76,76,76,255,maxColorValue=255); variantColours2[7] = rgb(60,60,60,100,maxColorValue=255) # dark grey (not used)

	# 12.1. Evolution of the number of clusters, averaged cluster size, and averaged duration since cluster TMRCA

nberOfDays = as.numeric(ymd("2022-03-31")-ymd("2020-08-01"))
minYear = decimal_date(ymd("2020-08-01")); maxYear = decimal_date(ymd("2022-03-31"))
timeSlices = nberOfDays; timePoints = seq(minYear,maxYear,(maxYear-minYear)/timeSlices)
mat1s = list(); mat2s = list(); mat3s = list(); mat4s = list()
for (h in 1:length(variants))
	{
		if (h == 1) nberOfExtractionFiles = 400 # TO BE REMOVED !!
		if (h != 1) nberOfExtractionFiles = 900 # TO BE REMOVED !!
		localTreesDirectory = paste0("DTA_boroughs_analyses/",variants[h],"_DTA/All_clades_ext")
		mat1a = matrix(nrow=timeSlices, ncol=nberOfExtractionFiles) # p2: evolution of the ratio between the number of circulating clusters and lineages (branches)
		mat2a = matrix(nrow=timeSlices, ncol=nberOfExtractionFiles) # p1a: evolution of the averaged proportion of circulating lineages belonging to the same cluster
		mat3a = matrix(nrow=timeSlices, ncol=nberOfExtractionFiles) # p1b: evolution of the probability that two circulating lineages drawn at random belong to the same cluster
		mat4a = matrix(nrow=timeSlices, ncol=nberOfExtractionFiles) # evolution of the averaged duration since cluster TMRCA
		for (i in 1:nberOfExtractionFiles)
			{
				tab = read.csv(paste0("DTA_boroughs_analyses/",variants[h],"_DTA/All_clades_ext/TreeExtractions_",i,".csv"), head=T)
				for (j in 2:length(timePoints))
					{
						sub = tab[which((tab[,"startYear"]<timePoints[j-1])&(tab[,"endYear"]>timePoints[j])),]
						circulatingClades = sub[,"cladeID"]; mat1a[j-1,i] = length(unique(circulatingClades))/dim(sub)[1]
						if (length(circulatingClades) != 0)
							{
								mat2a[j-1,i] = mean(table(circulatingClades))/dim(sub)[1]
								mat3a[j-1,i] = sum((table(circulatingClades)/dim(sub)[1])^2)
								tMRCAs = rep(NA, length(circulatingClades))
								for (k in 1:length(tMRCAs))
									{
										tMRCAs[k] = min(tab[which(tab[,"cladeID"]==circulatingClades[k]),"startYear"])
									}
								mat4a[j-1,i] = mean(c(timePoints[j-1],timePoints[j]))-mean(tMRCAs)
							}	else	{
								mat2a[j-1,i] = 0; mat4a[j-1,i] = NA
							}
					}
			}
		mat1b = matrix(nrow=timeSlices, ncol=4); colnames(mat1b) = c("time","median","lower95hpd","higher95hpd")
		mat2b = matrix(nrow=timeSlices, ncol=4); colnames(mat2b) = c("time","median","lower95hpd","higher95hpd")
		mat3b = matrix(nrow=timeSlices, ncol=4); colnames(mat3b) = c("time","median","lower95hpd","higher95hpd")
		mat4b = matrix(nrow=timeSlices, ncol=4); colnames(mat4b) = c("time","median","lower95hpd","higher95hpd")
		for (i in 1:dim(mat1b)[1])
			{
				mat1b[i,"time"] = mean(c(timePoints[i],timePoints[i+1])); mat1b[i,"median"] = median(mat1a[i,],na.rm=T)
				mat1b[i,c("lower95hpd","higher95hpd")] = quantile(mat1a[i,],c(0.025,0.975),na.rm=T)
				mat2b[i,"time"] = mean(c(timePoints[i],timePoints[i+1])); mat2b[i,"median"] = median(mat2a[i,],na.rm=T)
				mat2b[i,c("lower95hpd","higher95hpd")] = quantile(mat2a[i,],c(0.025,0.975),na.rm=T)
				mat3b[i,"time"] = mean(c(timePoints[i],timePoints[i+1])); mat3b[i,"median"] = median(mat3a[i,],na.rm=T)
				mat3b[i,c("lower95hpd","higher95hpd")] = quantile(mat3a[i,],c(0.025,0.975),na.rm=T)
				mat4b[i,"time"] = mean(c(timePoints[i],timePoints[i+1])); mat4b[i,"median"] = median(mat4a[i,],na.rm=T)
				mat4b[i,c("lower95hpd","higher95hpd")] = quantile(mat4a[i,],c(0.025,0.975),na.rm=T)
			}
		# N.B.: logical to do not observe variations among trees as the BSSVS was performed on a fixed tree topology
		mat1c = matrix(nrow=timeSlices, ncol=2); mat2c = matrix(nrow=timeSlices, ncol=2); mat3c = matrix(nrow=timeSlices, ncol=2); mat4c = matrix(nrow=timeSlices, ncol=2)
		colnames(mat1c) = c("time","sliddingW7d"); colnames(mat2c) = c("time","sliddingW7d"); colnames(mat3c) = c("time","sliddingW7d"); colnames(mat4c) = c("time","sliddingW7d")
		for (i in 8:(dim(mat1c)[1]-7))
			{
				indices = seq(i-7,i+7)
				mat1c[i,"time"] = mat1b[i,"time"]; mat1c[i,"sliddingW7d"] = mean(mat1b[indices,"median"],na.rm=T)
				mat2c[i,"time"] = mat2b[i,"time"]; mat2c[i,"sliddingW7d"] = mean(mat2b[indices,"median"],na.rm=T)
				mat3c[i,"time"] = mat3b[i,"time"]; mat3c[i,"sliddingW7d"] = mean(mat3b[indices,"median"],na.rm=T)
				mat4c[i,"time"] = mat4b[i,"time"]; mat4c[i,"sliddingW7d"] = mean(mat4b[indices,"median"],na.rm=T)
			}
		mat1c = mat1c[-which(is.na(mat1c[,"time"])),]; mat2c = mat2c[-which(is.na(mat2c[,"time"])),]; mat4c = mat4c[-which(is.na(mat4c[,"time"])),]
		mat1s[[h]] = mat1c; mat2s[[h]] = mat2c; mat3s[[h]] = mat3c; mat4s[[h]] = mat4c
	}

pdf(paste0("Figure_B2_NEW.pdf"), width=11, height=2.2) # dev.new(width=11, height=2.2)
par(mfrow=c(1,2), oma=c(0,0,0,0), mar=c(1.5,3.0,0,1), mgp=c(1.2,0.2,0), lwd=0.2, col="gray30")
plot(mat3s[[1]], col=NA, axes=F, ann=F, ylim=c(0,1))
for (h in 1:length(variants)) # reporting p1b
	{
		vS = 1-mat3s[[h]][,2]; xx_l = c(mat3s[[h]][,1],rev(mat3s[[h]][,1])); yy_l = c(rep(1,dim(mat3s[[h]])[1]),rev(vS))
		getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=variantColours2[h], border=0)
		lines(mat3s[[h]][,1], 1-mat3s[[h]][,2], lwd=1, col=variantColours1[h])
	}
for (h in 1:length(variants)) # reporting p2
	{
		lines(mat1s[[h]][,1], 1-mat1s[[h]][,2], lwd=1, col=variantColours1[h], lty=2)
	}
dates = c("2020-09-01","2021-01-01","2021-05-01","2021-09-01","2022-01-01","2022-05-01"); ats = decimal_date(ymd(dates))
axis(side=1, lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.03, col="gray30", col.axis="gray30", mgp=c(0,0.11,0), at=ats, label=dates)
axis(side=2, lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.03, col="gray30", col.axis="gray30", mgp=c(1,0.30,0), at=seq(0,1,0.2),
	 label=c("1","0.8","0.6","0.4","0.2","0"))
title(ylab="proportions p1 and p2", cex.lab=0.80, mgp=c(1.5,0,0), col.lab="gray30")
for (h in 2:length(variants))
	{
		disparsalStatistics = read.table(paste0("RRW_dispersal_statistics/",variants[h],"_estimated_dispersal_statistics.txt"), head=T)
		WDC = disparsalStatistics[,"weighted_diffusion_coefficient"]/365 # to get km2/day
		if (h == 2) plot(density(WDC), col=NA, axes=F, ann=F, xlim=c(0,30), ylim=c(0,0.43))
		polygon(density(WDC), col=variantColours2[h], border=NA); lines(density(WDC), lwd=1, col=variantColours1[h])
	}
axis(side=1, lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.03, col="gray30", col.axis="gray30", mgp=c(0,0.11,0))
axis(side=2, lwd.tick=0.2, cex.axis=0.65, lwd=0.2, tck=-0.03, col="gray30", col.axis="gray30", mgp=c(1,0.30,0))
title(ylab="weighted diffusion coefficient", cex.lab=0.80, mgp=c(1.5,0,0), col.lab="gray30")
dev.off()

	# 12.2. Averaged ratio between the number of county transition events and the cluster size
		  # (somehow a measure of a diffusion coefficient in the discrete phylogeographic framework)

wd = getwd()
for (h in 1:length(variants))
	{
		if (h == 1) nberOfExtractionFiles = 400 # TO BE REMOVED !!
		if (h != 1) nberOfExtractionFiles = 900 # TO BE REMOVED !!
		setwd(paste0(wd,"/DTA_boroughs_analyses/",variants[h],"_DTA/"))
		treeFiles = list.files(); treeFiles = treeFiles[which(grepl(".trees",treeFiles))]
		AR_TE_CS = matrix(nrow=nberOfExtractionFiles, ncol=1)
		for (i in 1:nberOfExtractionFiles)
			{				
				R_TE_CS = rep(NA, length(treeFiles))
				for (j in 1:length(treeFiles))
					{
						tab = read.csv(paste0(gsub(".trees","_ext",treeFiles[j]),"/TreeExtractions_",i,".csv"), head=T)
						R_TE_CS[j] = length(which(tab[,"startLoc"]!=tab[,"endLoc"]))/length(which(!tab[,"node2"]%in%tab[,"node1"]))
					}
				AR_TE_CS[i,1] = mean(R_TE_CS)
			}
		median = round(median(AR_TE_CS),2); quantiles = round(quantile(AR_TE_CS,c(0.025,0.975)),2)
		cat(variants[h],": ",median,", 95% HPD = [",quantiles[1],"-",quantiles[2],"]\n",sep="")
			# Iota:  0.37, 95% HPD = [0.35-0.39]
			# Alpha: 0.34, 95% HPD = [0.33-0.36]
			# Delta: 0.30, 95% HPD = [0.29-0.30]
			# O-BA1: 0.44, 95% HPD = [0.42-0.46]
	}
setwd(wd)

	# 12.3. Averaged number of different counties invaded by a distinct cluster (introduction event)

wd = getwd()
for (h in 1:length(variants))
	{
		if (h == 1) nberOfExtractionFiles = 400 # TO BE REMOVED !!
		if (h != 1) nberOfExtractionFiles = 900 # TO BE REMOVED !!
		setwd(paste0(wd,"/DTA_boroughs_analyses/",variants[h],"_DTA/"))
		treeFiles = list.files(); treeFiles = treeFiles[which(grepl(".trees",treeFiles))]
		AN_DC_DC = matrix(nrow=nberOfExtractionFiles, ncol=1)
		for (i in 1:nberOfExtractionFiles)
			{				
				N_DC_DC = rep(NA, length(treeFiles))
				for (j in 1:length(treeFiles))
					{
						tab = read.csv(paste0(gsub(".trees","_ext",treeFiles[j]),"/TreeExtractions_",i,".csv"), head=T)
						N_DC_DC[j] = length(unique(c(tab[,"startLoc"],tab[,"endLoc"])))
					}
				AN_DC_DC[i,1] = mean(N_DC_DC)
			}
		median = round(median(AN_DC_DC),2); quantiles = round(quantile(AN_DC_DC,c(0.025,0.975)),2)
		cat(variants[h],": ",median,", 95% HPD = [",quantiles[1],"-",quantiles[2],"]\n",sep="")
			# Iota:  2.93, 95% HPD = [2.90-3.03]
			# Alpha: 2.70, 95% HPD = [2.68-2.75]
			# Delta: 2.21, 95% HPD = [2.20-2.23]
			# O-BA1: 2.71, 95% HPD = [2.67-2.74]
	}
setwd(wd)

	# 12.4. Averaged proportion of phylogeny branches associated with a transition event between counties

wd = getwd()
for (h in 1:length(variants))
	{
		if (h == 1) nberOfExtractionFiles = 400 # TO BE REMOVED !!
		if (h != 1) nberOfExtractionFiles = 900 # TO BE REMOVED !!
		setwd(paste0(wd,"/DTA_boroughs_analyses/",variants[h],"_DTA/"))
		treeFiles = list.files(); treeFiles = treeFiles[which(grepl(".trees",treeFiles))]
		avgPropTransitionEvents = matrix(nrow=nberOfExtractionFiles, ncol=1)
		for (i in 1:nberOfExtractionFiles)
			{				
				propTransitionEvents = rep(NA, length(treeFiles))
				for (j in 1:length(treeFiles))
					{
						tab = read.csv(paste0(gsub(".trees","_ext",treeFiles[j]),"/TreeExtractions_",i,".csv"), head=T)
						propTransitionEvents[j] = length(which(tab[,"startLoc"]!=tab[,"endLoc"]))/dim(tab)[1]
					}
				avgPropTransitionEvents[i,1] = mean(propTransitionEvents)
			}
		median = round(median(avgPropTransitionEvents),2); quantiles = round(quantile(avgPropTransitionEvents,c(0.025,0.975)),2)
		cat(variants[h],": ",median,", 95% HPD = [",quantiles[1],"-",quantiles[2],"]\n",sep="")
			# Iota:  0.23, 95% HPD = [0.22-0.25]
			# Alpha: 0.22, 95% HPD = [0.21-0.24]
			# Delta: 0.20, 95% HPD = [0.20-0.21]
			# O-BA1: 0.30, 95% HPD = [0.29-0.31]
	}
setwd(wd)

