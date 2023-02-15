###TUTORIAL: HOW TO VISUALIZE PHYLOGENIES AND METADATA USING R
###IF you do not have R locally, you can use Rstudio Cloud https://rstudio.cloud/ 
## If you just want a simple way to view a phylogeny, you can using the program FigTree - link below
#https://github.com/rambaut/figtree/releases/tag/v1.4.4
#However, R will give you more flexibility and customization to make publication ready figures

##First you will need to install the R libraries used for this work 

#This only has to be run the first time you want to run this code in your local instance of R
#This code will install the requisite R packages
#We will mainly be using a program called GGTREE - you can learn a lot more about this tool here:
#https://bioconductor.org/packages/release/bioc/html/ggtree.html
#https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html 

install.packages("ape")
install.packages("ggtree")
install.packages("ggplot2")
install.packages("phytools")
#If installation of ggtree failes, you can try this
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
aBiocManager::install("ggtree")

#Then you can load the libraries
library(ape)
library(ggtree)
library(ggplot2)
library(phytools)

#Setting the working directory where your files are located
##If you don't do this you will need to specify the entire path when you import/export a file
setwd("/Users/tajazarian/Documents/Tree_tutorial/")

#import metadata
metadata <- read.delim("metadata.txt")

#import treefile
tree <- read.newick(file.choose())
#plot the tree as currently imported.  This will be un-ladderized
plot(tree)
#Now we can "ladderize" (i.e., order the tree by branch length) and midpoint root the tree
tree_ladderized <- ladderize(midpoint.root(tree))
#Now we can see the changes this made.  We can also turn off branch labels to make it easier to 
#view the tree topology
plot(tree_ladderized, show.tip.label = FALSE)

#Now we can use ggtree to view the tree. We can also show the node labels, which will allow us to place annotations on clades
ggtree(tree_ladderized) + 
  geom_text2(aes(label=node), hjust=-.3) + 
  geom_tiplab(cex=2.5)

#Using ggtree to view the tree and annotate with pertinent information
#We can label clades using the node information from above. We can also add annotations to the plot
annotated_tree <- ggtree(tree_ladderized, layout = "rectangular", ladderize=F) +
  theme(legend.position = "left") + 
  geom_treescale() +
  geom_cladelabel(node=94, label="ST121", align=FALSE, offset = .0, barsize = .1, fontsize = 3) +
  annotate(geom = 'text', label = '2174 core genes: 1.94 Mbps - SNP diversity: 538.00 SNPs', x = 0.0001, y = Inf, hjust = 0, vjust = 1, size=2)


#create metadata for heatmap by subsampling
#The row label needs to be the taxa IDs that are found in the tree file
#Each column will be a column in the heatmap. So you need to subsample the data accordingly
metadata_sub <- data.frame(Region = metadata$Region, Phylogroup = metadata$Phage_Phylo_Group)
  rownames(metadata_sub) <- metadata$Unique_ID #This changes row labels to the unique taxa ID

#Adding the heatmap to the tree 
annotatedtree_metadata <- gheatmap(annotated_tree, metadata_sub, offset = 0.0001, width=0.1, font.size=3, colnames_position= "top", colnames_angle = 45, colnames_offset_y = 0, hjust = 0) +
  scale_fill_manual(breaks=c("Africa","Asia","North America","UK","Europe","South America","Australia", "Unknown", "P05", "P07", "P08", "P09","P10", "P11", "P12", "P13", "P14", "P15", "P16", "P17", "P28", "P52", "P61"), values=c("Africa"="#899DA4","Asia"="#C93312","North America"="#FAEFD1","UK"="#DC863B","Europe"="#9A8822","South America"="#F5CDB4","Australia"="#74A089", "Unknown"="#FEFDFB", "P05"="#0F2E0F", "P07"="#2D4D19","P08"="#5C7326", "P09"="#63C600","P10"="#8BD000", "P11"="#B6DB00", "P12"="#E6E600", "P13"="#E7CE1D", "P14"="#E9BD3A", "P15"="#C76857", "P16"="#ECB176", "P17"="#D59E81", "P28"="#F1D6D3", "P52"="#D5D7F1", "P61"="#F2F2F2"))

#save heatmap+tree
