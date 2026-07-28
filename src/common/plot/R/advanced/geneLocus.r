######################
# Locus Zoom, Make LD Data
# September 2016/February 2017
# Tanya Flynn
# Uni of Otago

### Important Running Notes:
## You need to have acces to the data that the associations were originally done on OR data to calculate ethnicity with
## Genes.Data expects a data.frame with headers Gene, Chrom, Start, & End - this will be used to annotate the graph.

## LD Colours
# 1.0 - 0.8 = "#D53E4F"
# 0.8 - 0.6 = "#FDAE61"
# 0.6 - 0.4 = "#E6F598"
# 0.4 - 0.2 = "#66C2A5"
# 0.2 - 0.0 = "#3288BD"
# NA = "#5E4FA2"
# top-hit = "#9E0142"
library(RColorBrewer)


locus.zoom = function(data = NULL, snp = NA, gene = NA, region = NA, ld.file = NULL, offset = NA, genes.data = NULL, plot.title = NULL, nominal = 6, significant = 7.3, file.name = NULL, secondarysnp = NA, population = "EUR")
{
	# Load Data and define constants:
	data$SNP = as.character(data$SNP)
	LD.colours = data.frame(LD = c(seq(from = 0, to = 1, by = 0.1), "NaN"), Colour = c("#3288BD",rep(c("#3288BD", "#66C2A5", "#E6F598", "#FDAE61", "#D53E4F"), each = 2), "#5E4FA2"))

	# Stop if CHR:POS ID, but with no alias:
	if (!grepl('rs', snp) & is.null(alias)) {
		stop("You didn't provide an alias for your CHR:POS ID")
	}

	# Stop if there are duplicate SNPs:
	if (length(which(duplicated(data$SNP))) > 0) {
		stop('There are duplicates in your results file - Please remove them before running again')
	}

	# Get start and end regions for plotting and for pulling out data:
	if (all(is.na(region))) {
		# If SNP is given:
		if (!is.na(snp)) {
			snp.ind = which(data$SNP == snp)
			snp.chr = data$CHR[snp.ind]
			snp.pos = data$BP[snp.ind]
			region = c(snp.chr, snp.pos, snp.pos)
		} else if (!is.na(gene)) {
			# If Gene is given
			gene.ind = which(genes.data$Gene == gene)
			gene.chr = genes.data$Chrom[gene.ind]
			gene.start = genes.data$Start[gene.ind]
			gene.end = genes.data$End[gene.ind]
			region = c(gene.chr, gene.start, gene.end)
		} else {
			# If nothing was given
			stop("You must specify a SNP, Gene or Region to plot")
		}
	} else {
		offset = ifelse(is.na(offset), 0, offset)
	}

	# Now re-define region to work with:
	region[2] = region[2] - offset # start position
	region[3] = region[3] + offset # end position

	# Pull out the relevant information from the UCSC gene data.
	# Any gene that overlaps/intersect with the defined region is
	# included:
	genes.data = genes.data[genes.data$Chrom == region[1], ]
	genes.data = genes.data[genes.data$End > region[2], ]
	genes.data = genes.data[genes.data$Start < region[3], ]

	# Likewise, pull out the relevant data from the results file:
	data = data[data$CHR == region[1], ]
	data = data[data$BP >= region[2], ]
	data = data[data$BP <= region[3], ]

	# Identify SNP with most significant association:
	lead.ind = which(data$P %in% min(data$P, na.rm = T))
	lead.snp = data$SNP[lead.ind]
	lead.chr = data$CHR[lead.ind]
	lead.pos = data$BP[lead.ind]
	lead.p = data$P[lead.ind]
	print(lead.snp)
	print(lead.pos)
	print(lead.p)

	# If LD information is not supplied, calculate it from the 1000 genomes
	# data:
	if (is.null(ld.file)) {
		ld.file = get.ld(region, lead.snp, population)
	}

	# Check if results file is rsID-based or CHR:POS-based
	if (length(which(grepl('rs', data$SNP))) < 1) {
		ld.file$SNP_B = paste(paste0("chr",ld.file$CHR_B), ld.file$BP_B, sep = ':')
	}

	# Add LD to Results
	new.ld.file = ld.file[,c("SNP_B", "R2")]

	round.up = function(x, decimals = 1){
		round(x + (5 * 10 ^ (-decimals - 1)), digits = decimals)
	}

	data = merge(data, new.ld.file, by.x = "SNP", by.y = "SNP_B", all.x = TRUE)

	data$plot.ld = round.up(data$R2, decimals = 1)
	data$plot.ld[data$plot.ld > 1 & !is.na(data$plot.ld)] = 1

	data.plot = merge(data, LD.colours, by.x = "plot.ld", by.y = "LD", all.x = TRUE)

	# Make Plotting Variables
	y.max = max(-log10(lead.p), 8)
	x.min = region[2]
	x.max = region[3]
	y_offset = y.max / 20

	# Make Plot
	jpeg(width = 150, height = 160, units = "mm", res = 300, file = file.name)
	layout(matrix(c(1, 2, 3), byrow = TRUE), heights = c(2, 6, 3))

	## Plot 1 - SNP presence
	par(mar = c(0.4, 4, 2, 2), mgp = c(2, 1, 0))
	plot(x = data.plot$BP[-log10(data.plot$P) > significant], y = rep(1, times = nrow(data.plot[-log10(data.plot$P) > significant,])), axes = FALSE, pch = "|", xlab = "", ylab = "SNPs", cex.lab=1.3, las = 1, main = plot.title, xlim = c(x.min, x.max))

	## Plot 2 - Manhattan/LocusZoom
	par(mar = c(0.5, 4, 0.4, 2), mgp = c(2, 1, 0))
	plot(x = data.plot$BP, y = -log10(data.plot$P), ylim = c(0, y.max*1.1), pch = 20, col = as.character(data.plot$Colour), xlab = "", ylab = expression(-log[10](italic(P))), cex = 2, cex.lab = 1.3, xaxt = "n",yaxt="n", xlim = c(x.min, x.max))
	# Plot the lead SNP
	axis(side=2,las=1,cex.axis=1.3)
	points(x = lead.pos, y = -log10(lead.p), pch = 18, cex = 3, col = "#9E0142")
	abline(h = nominal, col = "black", lty = 2)
	abline(h = significant, col = "black", lty = 1)
	text(x = lead.pos[1], y = (-log10(lead.p)[1] + y_offset), labels = "chr8:2130270", cex = 1.5)

	# TODO: perhaps put the code block below into a separate function
	## Extra - plotting the label/text for the secondary SNP
	# Calculate the coordinates to put the labels for the secondary SNP
	if(!is.na(secondarysnp)){
		x_offset = abs(x.max - x.min) / 10
		ind = which(data$SNP %in% secondarysnp)
		if (length(ind) != 1) {
			stop('The secondary SNP was not in the region you wanted to plot:\n\tPlease check the location of your secondary SNP')
		}
		snp.pos = data$BP[ind]
		snp.p = -log10(data$P[ind])
		label.x.offset = snp.pos + x_offset
		line.x.offset = snp.pos + (x_offset / 3)
		data = data[data$BP > (snp.pos - (abs(x.max - x.min) / 6)) & data$BP < (snp.pos + (abs(x.max - x.min) / 6)), ]
		label.y.offset = snp.p * 1.03
		if(abs(nominal - label.y.offset) <= 1){
			label.y.offset = max(label.y.offset, nominal) * 1.03
		}
		# Label the secondary SNP
		text(x = label.x.offset, y = label.y.offset, labels = secondarysnp)
		segments(x0 = snp.pos, x1 = line.x.offset, y0 = snp.p, y1 = label.y.offset)
	}

	# Make a legend
	legend(x = "topright", legend = c("> 0.8", "> 0.6", "> 0.4", "> 0.2","< 0.2"), col = c("#D53E4F","#FDAE61","#E6F598","#66C2A5","#3288BD"), border = c("#D53E4F","#FDAE61","#E6F598","#66C2A5","#3288BD"), pch = 16, pt.cex = 1.5, cex = 1.2, bg = "white", bty = "n", title = expression(italic("r")^2))

	## Plot 3 - Gene tracks
	par(mar = c(3, 4, 0.5, 2), mgp = c(2, 1, 0))
	plot(1, type = "n", yaxt = "n", xlab = paste("Position on Chromosome", "8"), ylab="", xlim = c(x.min, x.max), ylim = c(0,2), cex.lab = 1.4, cex.axis = 1.25)

	# Function to plot the gene tracks and labels properly:
	gene.position = function(data) {
		odd = 1
		for (i in 1:length(data$Gene)) {
			lines(x = c(data$Start[i], data$End[i]), y = c(data$Y[i], data$Y[i]), lwd = 3, col = "#000080")
			if(length(data$Gene) >= 10){
				length = abs(data$Start[i] - data$End[i])
				if(length > 8000) {
					if (odd%%2 == 0) {
						text(x = (data$Start[i] + data$End[i])/2, y = data$Y[i] - 0.2, labels = as.character(data$Gene[i]), cex = 1.1)
					} else {
						text(x = (data$Start[i] + data$End[i])/2, y = data$Y[i] + 0.2, labels = as.character(data$Gene[i]), cex = 1.1)
					}
				}
			} else {
				text(x = (data$Start[i] + data$End[i])/2, y = data$Y[i] + 0.2, labels = as.character(data$Gene[i]),cex = 1.1)
			}
			odd = odd + 1
		}
	}

	# Stagger the genes
	y = rep(c(1.5, 0.5), times = length(genes.data$Gene))
	genes.data$Y = y[1:length(genes.data$Gene)]
	genes.top = genes.data[genes.data$Y == 1.5,]
	genes.bot = genes.data[genes.data$Y == 0.5,]

	# Plot the gene tracks:
	if (nrow(genes.top) > 0) {
		gene.position(genes.top)
	}
	if (nrow(genes.bot) > 0) {
		gene.position(genes.bot)
	}

	dev.off()
}



# Function to get the LD information of specified population from the 1000
# Genomes data (March 2017 release).
#
# Note that the input SNP MUST be in rsID format, not CHR:POS-based.
#
# This function will leave/save the LD information in the working directory for
# future reference (e.g. if the user wanted to use the same LD information)
get.ld = function(region, snp, population) {
	ld.snp = snp
	# If the SNP is in CHR:POS ID, then find the rsID:
	if (!grepl('rs', snp)) {
		command = "awk '{print $1\":\"$2, $3}' /Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/1kgp_chrZZ_biallelic_snps.txt | grep -w SNP"
		command = gsub('ZZ', region[1], command)
		command = gsub('SNP', snp, command)
		rsid = system(command, intern = T)
		rsid = unlist(strsplit(rsid, ' '))[2]
		ld.snp = rsid
	}

	# gsub the command and filename for chr, start/end positions and the
	# population:
	base.command = "bcftools view -r ZZ:Y1-Y2 -S /Volumes/archive/merrimanlab/riku/1kgp_sample_files/POP.sample /Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/ALL.chrZZ.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -Oz -o tmp.vcf.gz && plink2 --vcf tmp.vcf.gz --allow-no-sex --snps-only --r2 --inter-chr --ld-snp SNP --ld-window-r2 0 --out POP_region_ZZ_Y1_Y2 && rm tmp.vcf.gz POP_region_ZZ_Y1_Y2.nosex"
	base.command = gsub('ZZ', region[1], base.command)
	base.command = gsub('Y1', region[2], base.command)
	base.command = gsub('Y2', region[3], base.command)
	base.command = gsub('POP', population, base.command)
	base.command = gsub('SNP', ld.snp, base.command)

	# Make a system call to run the bcftools/plink command.
	# I'm only assigning it to a variable to suppress any possible form of
	# messages/outputs from the command, just in case
	messages = system(base.command, ignore.stdout = T, intern = T)

	# Import the LD data:
	ld.file = "POP_region_ZZ_Y1_Y2.ld"
	ld.file = gsub('ZZ', region[1], ld.file)
	ld.file = gsub('Y1', region[2], ld.file)
	ld.file = gsub('Y2', region[3], ld.file)
	ld.file = gsub('POP', population, ld.file)

	# Import the ld file and return the data
	ld = read.table(ld.file, stringsAsFactors = F, header = T)
	return(ld)
}
