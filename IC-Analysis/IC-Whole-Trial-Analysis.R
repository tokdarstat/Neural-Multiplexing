## Code for carrying out whole trial Poisson analysis
## of IC data. Intendedn to be run in the batch mode


## Filenames from JA's and VC's experiments

fnames.ja230.int <- scan("JA230_list_randomly_interleaved.txt", "a")
fnames.vc <- scan("VC_extended_list.txt", "a")

## Codes needed
source("poisson-helper-functions.R")
require(parallel)

pdf(height = 8, width = 5, file = "poi-ja-600-zero-removed.pdf")
for(fname in fnames.ja230.int){
  try(poi.JA(fname, outfile="bfs-ja-int-end600-zero-removed.txt", start=0, end=600))
}
dev.off()

pdf(height = 8, width = 5, file = "poi-vc-sof-zero-removed.pdf")
for(fname in fnames.vc){
  try(poi.VC(fname, cell=1, outfile="bfs-vc-sof-zero-removed.txt", start=0, end=1000, go.by.soff=TRUE))
  try(poi.VC(fname, cell=2, outfile="bfs-vc-sof-zero-removed.txt", start=0, end=1000, go.by.soff=TRUE))
}
dev.off()


q("no")

