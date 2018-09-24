## required R codes
source("dapp-helper-functions.R)
## list of cells to analyze
triplet.meta <- read.csv("363triplets.csv")
## local path to pull data from (set to current directory)
local.path <- "./"
## if pulling data from the web
web.url <- "http://www2.stat.duke.edu/~st118/Jenni/STCodes/"

require(parallel)
trip.set <- 1:363

BW <- 50
data.path <- paste0(local.path, "Data")
save.path <- paste0(local.path, "ResultsBW", BW)
if(!dir.exists(save.path)) dir.create(save.path)
if(!dir.exists(paste0(save.path, "/Figures"))) dir.create(paste0(save.path, "/Figures"))
if(!dir.exists(paste0(save.path, "/Summaries"))) dir.create(paste0(save.path, "/Summaries"))

all.set <- mclapply(trip.set, function(jj) try(fitter.fn(jj, triplet.meta, plot = FALSE, verbose = FALSE, bw = BW, save.figure = TRUE, save.out = TRUE, data.path = data.path, local.pull = TRUE, save.path = save.path)), mc.cores = 4)
trip.err <- (sapply(all.set, class) == "try-error")


### Post processing for summary and plotting

file.tag <- rep(NA, length(trip.set))
for(jj in 1:length(trip.set)){
    triplet <- trip.set[jj]
    freq <- triplet.meta[triplet, "AltFreq"]
    angl <- triplet.meta[triplet, "AltPos"]
    cell <- triplet.meta[triplet, "PairId"]
    fname <- triplet.meta[triplet, "CellId"]
    file.tag[jj] <- paste0(fname, "_cell", cell, "_freq", freq, "_pos", angl)
}

triplet.combo <- cbind(triplet.meta, tag = file.tag)

summary.files <- list.files(paste0(save.path, "/Summaries"))
nfiles <- length(summary.files)
ncuts <- 4
summary.meta <- as.data.frame(matrix(nrow = nfiles, ncol = 1+4*ncuts))
names(summary.meta) <- c("tag", paste0("mmPriorQ", 1:ncuts), paste0("mmPostQ", 1:ncuts), paste0("avPriorQ", 1:ncuts), paste0("avPostQ", 1:ncuts))
for(jj in 1:nfiles){
    ff <- summary.files[jj]
    load(paste0(save.path, "/Summaries/", ff))
    ff.tag <- substring(ff, 1, nchar(ff)-3)
    sm <- summary(fit.post, tilt = TRUE, cut = 1/ncuts, mesh = 1/ncuts, nprior = 1e4)
    summary.meta[jj,1] <- ff.tag
    summary.meta[jj,-1] <- c(sm)
    rm(fit.post)
}

select.tags <- c("YKIC092311_1Loc_DoubleSound_cellNA_freq609_pos24",
"YKIC140121Loc_DoubleSound_cell1_freq500_pos24",
"YKIC141016Loc_DoubleSound_cell1_freq500_pos-24",
"YKIC080311Loc_DoubleSound_cellNA_freq1100_pos24",
"YKIC140501Loc_DoubleSound_cell1_freq500_pos24",
"YKIC140226Loc_DoubleSound_cell1_freq1988_pos-6",
"YKIC102711Loc_DoubleSound_cellNA_freq1100_pos24"
)

triplet.final <- merge(triplet.combo, summary.meta, by = "tag", all = TRUE)
triplet.final$wavieness <- with(triplet.final, (mmPostQ3 + mmPostQ4)/(mmPostQ1 + mmPostQ2))
triplet.final$centrality <- with(triplet.final, (avPostQ2 + avPostQ3)/(avPostQ1 + avPostQ4))
triplet.final$skewness <- with(triplet.final, pmax((avPostQ3 + avPostQ4)/(1e-3 + avPostQ1 + avPostQ2), (avPostQ1 + avPostQ2)/(1e-3 + avPostQ3 + avPostQ4)))
triplet.final$scaledskew <- log(pmax(1, pmin(triplet.final$skewness, 100))) / log(1e2)
triplet.final$chosen.tag <- ""
select.triplets <- sapply(select.tags, match, table = triplet.final$tag)
triplet.final$chosen.tag[select.triplets] <- as.character(triplet.final$tag[select.triplets])
ggplot(triplet.final, aes(centrality, wavieness)) + geom_point(aes(color = rank(skewness), shape = WinModels), cex = 4, alpha = 0.7) + scale_colour_gradientn(colours=rainbow(2)) + geom_hline(yintercept = 1, lty = 2, col = "gray") + geom_vline(xintercept = 7/3, lty = 2, col = "gray") + geom_text(aes(label = chosen.tag), size = 2, vjust = -.01, hjust = -.01)


three.tags <- list(wavieness = c("flat", "", "wavy"),
centrality = c("extreme", "", "central"),
skewness = c("symmetric", "", "skewed"))
pr.cuts <- list(wavieness = c(0, .78, 1.3, Inf), centrality = c(0,1.68,3.24,Inf), skewness = c(0,2,4,Inf))
#Prior wavieness = 1, centrality = 2.33 and skewness = 1.
# The cut-offs used above are resonable extremeties -- but could be altered.

three.summaries <- data.frame(matrix(NA, nrow = nrow(triplet.final), ncol = 3))
names(three.summaries) <- names(three.tags)

for(sm.name in c("centrality", "wavieness", "skewness")){
    groups <- cut(triplet.final[,sm.name], pr.cuts[[sm.name]], labels = FALSE, include.lowest = TRUE)
    three.summaries[,sm.name] <- three.tags[[sm.name]][groups]
}
triplet.final$dynamic.analysis.tag <- apply(three.summaries, 1, paste, collapse = "-")

whole3summ <- with(triplet.final, table(WinModels, dynamic.analysis.tag))
whole3summ <- whole3summ / rowSums(whole3summ)
whole3summ <- whole3summ[,apply(whole3summ, 2, max) >= 0.1]
whole3summ <- cbind(whole3summ, other = 1 - rowSums(whole3summ))
class(whole3summ) <- "table"
whole3summ.df <- as.data.frame(whole3summ)
names(whole3summ.df)[1:2] <- c("WinModels", "dynamic.analysis.tag")
ggplot(subset(whole3summ.df, Freq > 0.1), aes(WinModels, dynamic.analysis.tag)) + geom_tile(aes(fill = Freq)) + scale_fill_gradientn(colors = rainbow(4)) +  labs(title = "All 363 triplets, BW = 50")

top.cut <- 0.95
whole3summ.top <- with(subset(triplet.final, WinPr > top.cut), table(WinModels, dynamic.analysis.tag))
whole3summ.top <- whole3summ.top / rowSums(whole3summ.top)
whole3summ.top <- whole3summ.top[,apply(whole3summ.top, 2, max, na.rm = TRUE) >= 0.1]
whole3summ.top <- cbind(whole3summ.top, other = 1 - rowSums(whole3summ.top))
class(whole3summ.top) <- "table"
whole3summ.top.df <- as.data.frame(whole3summ.top)
names(whole3summ.top.df)[1:2] <- c("WinModels", "dynamic.analysis.tag")
ggplot(subset(whole3summ.top.df, Freq > 0.1), aes(WinModels, dynamic.analysis.tag)) + geom_tile(aes(fill = Freq)) + scale_fill_gradientn(colors = rainbow(4)) +  labs(title = paste0("WinPr > ", top.cut, ", BW = 50"))

## whole trial summary
qplot(factor(WinModels), data = triplet.final, fill = factor(WinPr < 0.95))
qplot(factor(WinModels), data = triplet.final, fill = factor(SepBF < 7))


## dynamic analysis
ggplot(triplet.final, aes(dynamic.analysis.tag, fill=(WinPr < 0.95))) + geom_bar() + coord_flip()
ggplot(triplet.final, aes(dynamic.analysis.tag, fill=(SepBF > 7.0))) + geom_bar() + coord_flip()


sep.cut <- 14.0
whole3summ.top <- with(subset(triplet.final, SepBF > sep.cut), table(WinModels, dynamic.analysis.tag))
whole3summ.top <- whole3summ.top / rowSums(whole3summ.top)
whole3summ.top <- whole3summ.top[,apply(whole3summ.top, 2, max, na.rm = TRUE) >= 0.1]
whole3summ.top <- cbind(whole3summ.top, other = 1 - rowSums(whole3summ.top))
class(whole3summ.top) <- "table"
whole3summ.top.df <- as.data.frame(whole3summ.top)
names(whole3summ.top.df)[1:2] <- c("WinModels", "dynamic.analysis.tag")
ggplot(subset(whole3summ.top.df, Freq > 0.1), aes(WinModels, dynamic.analysis.tag)) + geom_tile(aes(fill = Freq)) + scale_fill_gradientn(colors = rainbow(4)) +  labs(title = paste0("WinPr > ", top.cut, ", BW = 50"))
