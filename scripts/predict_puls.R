#!/usr/bin/env Rscript

#################################
#
# Functions
#
#################################

process_pairs <- function(idx, td) {

        mindbc <- dbcan[dbcan$protein_id %in% td$protein_id,]

	mindbc <- unique(mindbc[,c("protein_id","hmm")])
	
	gh <- data.frame(protein_id="",hmm="", stringsAsFactors=FALSE)
	if (nrow(mindbc) > 0) {
        	gh <- aggregate(mindbc$hmm, by=list(protein_id=mindbc$protein_id), paste, collapse=";")
        }
	colnames(gh) <- c("protein_id","hmm")

        cmin = idx -1
        while(cmin > 0) {

                contains = 0
                for (j in cmin:(cmin-4)) {

                        if (j>0) {

				dist <- td$start[j+1] - td$end[j]
				if (dist > 500) {
	                                break
        	                }

                                # check if it is susC/D
                                if (td$sus[j] == "susC" || td$sus[j] == "susD") {
                                        contains = 1
                                }

                                # check if it is a GH
                                prot_id = td$protein_id[j]
                                ghvec <- gh$hmm[gh$protein_id==prot_id]
                                if (length(ghvec) > 0) {
                                        contains = 1
                                }

                        }

                }

                if (contains > 0) {
                        cmin <- cmin -1
                } else {
                        break
                }
        }

        cmin <- cmin + 1

        cmax = idx + 1
        while(1) {

                contains = 0
                for (j in cmax:(cmax+4)) {

                        if (j <= nrow(td)) {

				dist <- td$start[j] - td$end[j-1]
				if (dist > 500) {
					break
				}

                                # check if it is susC/D
                                if (td$sus[j] == "susC" || td$sus[j] == "susD") {
                                        contains = 1
                                }

                                # check if it is a GH
                                prot_id = td$protein_id[j]
                                ghvec <- gh$hmm[gh$protein_id==prot_id]
                                if (length(ghvec) > 0) {
                                        contains = 1
                                }
                        }

                }

                if (contains > 0) {
                        cmax <- cmax + 1
                } else {
                        break
                }
        }

        cmax <- cmax - 1

	if (cmax > cmin) {

        	pultdf <- merge(td[cmin:cmax,], gh, all.x=TRUE, sort=FALSE)
        	pultdf <- pultdf[order(pultdf$order),]
        	pultdf[is.na(pultdf)] <- ""
        	pultdf$pulid <- rep(paste("PUL",PULCOUNTER, sep=""), nrow(pultdf))

		print(paste("created PUL", PULCOUNTER))
		print(paste(cmin,"-",cmax))

        	ALLPULS <<- rbind(ALLPULS, pultdf)

        	PULCOUNTER <<- PULCOUNTER + 1
	}
        return(cmax)

}

#################################
#
# code
#
#################################

# get command line arguments as an array
args <- commandArgs(trailingOnly = TRUE)

# arguments
ft.name <- args[1]
pf.name <- args[2]
db.name <- args[3]

# file base
file.out <- args[4]
file.sum <- args[5]

# genome
genome <- args[6]

cn <- c("contig","start","end","strand","protein_id","protein_name")
ft <- read.table(ft.name, 
		header=FALSE, 
		sep="\t", 
		stringsAsFactors=FALSE, 
		col.names=cn)

ft$order = 1:nrow(ft)

cn <- c("protein_id","as","ae","es","ee","ha","hn","type","hs","he","hl","bit","e","sig","clan","active")
pfam <- read.table(pf.name, 
		skip=29, header=FALSE, 
		stringsAsFactors=FALSE, 
		fill=TRUE, 
		col.names=cn)

cn <- c("hmm","hit_len","protein_id","query_len","evalue","hit_start","hit_end","query_start","query_end","cov")
dbcan <- read.table(db.name, 
		sep="\t", 
		header=FALSE, 
		stringsAsFactors=FALSE, 
		col.names=cn)

dbcan$hmm <- gsub(".hmm","",dbcan$hmm)

# GLOBAL VARIABLES
PULCOUNTER <- 1
ALLPULS <- data.frame(protein_id=character(),
		      contig=character(),
		      start=numeric(),
		      end=numeric(),
		      strand=character(),
		      protein_name=character(),
		      order=numeric(), 
		      active=character(), 
		      sus=character(), 
		      hmm=character(), 
		      pulid=character(),
		      stringsAsFactors=FALSE)


# filter pfam
halign_prop = (pfam$he - pfam$hs + 1) / pfam$hl
pfam <- pfam[halign_prop>=0.6,]

# find and annotate susC/susD
pfam$sus <- rep("none", nrow(pfam))
pfam$sus[grep("PF00593",pfam$ha)] <- "susC"
pfam$sus[grep("PF13715",pfam$ha)] <- "susC"
pfam$sus[grep("PF07715",pfam$ha)] <- "susC"
pfam$sus[grep("PF07980",pfam$ha)] <- "susD"
pfam$sus[grep("PF12741",pfam$ha)] <- "susD"
pfam$sus[grep("PF12771",pfam$ha)] <- "susD"
pfam$sus[grep("PF14322",pfam$ha)] <- "susD"

# limit pfam to relevant rows
pfam <- pfam[pfam$sus != "none",]

# merge with feature table
ftp <- merge(ft, pfam, by.x="protein_id", by.y="protein_id", all.x=TRUE, sort=FALSE)[,c(1:7,22,23)]
ftp <- unique(ftp)
ftp <- ftp[order(ftp$order),]

# unique list of contigs
cons <- unique(ftp$contig)

# go through each contig
for (c in cons[order(cons)]) {

	tdf <- ftp[ftp$contig==c,]
	tdf <- tdf
	tdf[is.na(tdf)] <- ""

	i <- 1
	while(i <= nrow(tdf)) {
	#for (i in 1:nrow(tdf)) {

		#print(c)
		#print(nrow(tdf))
		#print(i)
		if (i>=nrow(tdf)) {
			break
		}
		if (tdf$sus[i] == "susC" && tdf$sus[i+1] == "susD") {
			# we have a pair - do something
			print(paste("pair at",c,i))
			i <- process_pairs(i,tdf)
			print(i)
		}

		if (i>=nrow(tdf)) {
                        break
                } 

		if (tdf$sus[i] == "susD" && tdf$sus[i+1] == "susC") {
			# we have a pair - do something
			print(paste("pair at",c, i))
			i <- process_pairs(i,tdf)
			print(i)
		}
		
		i <- i+1
	}
	
}


if (nrow(ALLPULS) >= 2) {
	two <- ALLPULS$end[1:(nrow(ALLPULS)-1)]
	one <- ALLPULS$start[2:nrow(ALLPULS)]

	ALLPULS$dist <- c(0,one - two)

	ALLPULS$genome <- rep(genome, nrow(ALLPULS))

	ALLPULS <- ALLPULS[, c("genome","pulid","protein_id","contig","start","end","strand","dist","protein_name","sus","hmm","active")]
	
	ALLPULS <- ALLPULS[order(ALLPULS$contig, ALLPULS$start),]

	write.table(ALLPULS, file.out, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

} else {
	ALLPULS$dist <- numeric()
	file.create(file.out)
}


if (nrow(ALLPULS) >= 2) {

	ALLPULS$GENE <- ALLPULS$hmm
	ALLPULS$GENE[ALLPULS$sus!=""] <- ALLPULS$sus[ALLPULS$sus!=""]
	ALLPULS$GENE[ALLPULS$GENE==""] <- "unk"

	agg <- aggregate(ALLPULS$GENE, by=list(pulid=ALLPULS$pulid), paste, sep="-", collapse="-")
	uni <- unique(ALLPULS[,c("genome","pulid","contig")])
	start <- aggregate(ALLPULS$start,  by=list(pulid=ALLPULS$pulid), function(x) return(x[1]))
	end   <- aggregate(ALLPULS$end,  by=list(pulid=ALLPULS$pulid), function(x) return(x[length(x)]))

	pos <- merge(start, end, by="pulid")
	deets <- merge(uni, pos, by="pulid")

	out <- merge(deets,agg,by="pulid")

	colnames(out) <- c("pulid","genome","contigid","start","end","pattern")

	out <- out[,c("genome","pulid","contigid","start","end","pattern")]

	out <- out[order(out$contigid, out$start),]

	write.table(out, file.sum, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
} else {
	file.create(file.sum)
}


