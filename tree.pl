#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;

use constant TEMP_LOCATION             => '/scratch';
use constant TEMP_DIR                  => tempdir(DIR => TEMP_LOCATION, CLEANUP => 1);

# Handle input options
my ($Help, $Last);
my $Method   = "flat";
my $Title    = "Pan-genomic Dendrogram";
my $ProgName = "tree.pl";
my $Version  = "0.2";

my %Methods = ("flat"  => "All weights set to 1",   "zero"  => "ORFans weight 0, all other weights 1",
               "shell" => "Shell genes emphasized", "cloud" => "Cloud genes emphasized");

# Main options
&GetOptions (
	"l|last=i"          => \$Last,   # Number of last organism to include - defaults to the last in the coreplot
	"m|method=s"        => \$Method, # Select the method used, should be in keys %Methods
	"t|title=s"         => \$Title,  # Title for plot
	"h|help"            => \$Help    # Prints the help text
);

# Want help?
print_help() if (defined $Help);

my $Dir = $ARGV[0];
die("ERROR! You must specify a valid directory. Run coregenome version 2.3 or later with the '-keep' option.\nRun $ProgName with the '-h' option for more help.\n")
  unless ((-r $Dir) && (-e "$Dir/tbl"));

$Last  = (defined $Last ? $Last - 1 : `tail -1 $Dir/tbl | cut -f 1`);
chomp $Last;





# --------------------------------------------------------------------
# %% Main Program %%
#

open TBL, "$Dir/tbl" or die;
my %md2id;
while (<TBL>) {
  next unless (/^([0-9]+)/ && $1 <= $Last);
  my ($id, $desc) = split /\t/, $_;
  my $md5 = `head -1 $Dir/$id.fsa | perl -pe 's/^>|[.][0-9]+|\\s.*//g'`;
  $md2id{$md5} = $desc;
}
close TBL;

open GRP, "$Dir/group_".$Last.".dat" or die;
my (%Out, $NrFam);
while (<GRP>) {
  my ($group_id, $count, @members) = split /\t/;
  my %repr;
  for (@members) {
    next unless /^(.*)(\..*)/;
    $repr{$1} = 1;
  }
  while (my ($k, $v) = each %md2id) {
    $Out{$v}->[$group_id-1] = (exists $repr{$k} ? 1 : 0);
  }
  $NrFam = $group_id;
}
close GRP;


my ($fh, $file) = tempfile(DIR => TEMP_DIR);
print $fh join("\t", 1..$NrFam), "\n";
my @str;
for (keys %Out) {
  push @str, join("\t", $_, @{ $Out{$_} });
}
print $fh join("\n", @str);
close $fh;


open R , "| /usr/cbs/tools/bin/linux64/R-2.8.1 --slave --vanilla";
my $title = (defined $ARGV[0] ? $ARGV[0] : "R Graphics Output");
source_lars();
print R "
data <- as.matrix(read.table('$file', header=TRUE, row.names=1, sep='\t', as.is=TRUE));
postscript('| cat', title = '$title');
w <- computeWeights(data, method='$Method');
panGenomeTree(data, weights=w, nboot=100, main='$Title');

#data.dist <- dist(data, method='manhattan');
#data.hclust <- hclust(data.dist);
#data.dend <- as.dendrogram(data.hclust);
#op <- par(mar = c(4,0,3,8));
#nP <- list(cex=c(1.0, 0.9), lab.cex = 0.9);
#plot(data.dend, nodePar = nP, horiz = T);

graphics.off();
q();
";

close R;




# --------------------------------------------------------------------
# %% Help Page/Documentation %%
#

sub print_help {
  my @method;
  while (my @a = each %Methods) {
    push @method, join("\t# ", @a);
  }
  my $method = join("\n\t\t", @method), "\n";
  print <<EOH;

NAME
	$ProgName - Draws whole genome based dendrograms

SYNOPSIS
	$ProgName [Options] [Directory]

DESCRIPTION
	Outputs a postscript with a dendrogram representation of the data in
	[Directory]. The directory argument must be of the type	created as
	temporary directory by the coregenome script version 2.3 or later.

OPTIONS
     -l or --last [integer]
	The index of the last organism from the core genome plot to include in
	the calculation of the core genome. Defaults to the last organism in the
	plot.

     -m or --method [string]
	Select the method used for clustering and constructing the dendrogram.
	The default is "$Method" and possible values are:
		$method

     -t or --title [string]
	Set the title for the plot.

     -h or --help
	Prints this text.

VERSION
    Current: $Version

    This is the development version.


AUTHOR
    Carsten Friis, carsten\@cbs.dtu.dk,
    Integrates code written by Lars Snipen and users should consider citing:

	"Microbial comparative pan-genomics using binomial mixture models"
	Snipen L, Almøy T, Ussery DW.
	http://www.ncbi.nlm.nih.gov/pubmed/19691844

EOH
  exit;
}




#--------------------------------------------------------#
#---- Functions for creating pan-genome family trees ----#
#--------------------------------------------------------#
    
sub source_lars {
  print R '
panGenomeTree <- function( X, weights=NA, main=NA, linkage="average", nboot=0 ){
#   Creating a pan-genome tree from a pan-matrix.
#   
#   mapping.frame is used to map the GID-code in the pan-matrix to the Name specified for the
#   corresponding GID in mapping.frame. mapping.frame is a data.frame with columns GID, File and Name.
#   weights is a vector with a numeric for each gene family, see computeWeights().
#   linkage could be "single", "average" or "complete", see hclust().
#   nboot is the number of bootstrap samples used to produce bootstrap-values.
#   
#   Originally by Lars Snipen,
#   Biostatistics group, Norwegian University of Life Sciences.
    dd <- dim( X )
    if( !is.na( weights[1] ) ){
        X <- X * matrix( weights, nrow=dd[1], ncol=dd[2], byrow=T )
        N <- sum( weights )
    } else {
        N <- dd[2]
    }
    dtab <- dist( X, method="manhattan" )/N
    clst <- hclust( dtab, method=linkage )
    maxl <- max( nchar( rownames( X ) ) )    
    mar <- c(5,1,1,10)
    mar[4] <- 2 + maxl*0.8*par("cin")[1]*2.54
    par( mar=mar )
    plot( as.dendrogram( clst ), horiz=T, nodePar=list(lab.cex=0.8, pch=16, cex=c(0.1,0.5)), xlab="Relative manhattan distance", main=main )
    if( nboot > 0 ){
        cat( "...bootstrapping", file=stderr() )
        signatur <- clusterSignature( clst$merge )
        nb <- dim( clst$merge )[1]
        nbranch <- rep( 0, times=nb )
        for( i in 1:nboot ){
            cat( ".", file=stderr() )
            idx <- sample( (1:dd[2]), dd[2], replace=T )
            Xb <- X[,idx]
            dtab <- dist( Xb, method="manhattan" )
            cls <- hclust( dtab, method=linkage )
            sig <- clusterSignature( cls$merge )
            nbranch <- nbranch + as.numeric( signatur %in% sig )
        }
        cat( "\n", file=stderr() )
        lab.nod <- as.character( round( 100*nbranch/nboot ) )
        bp <- branchPos( clst$merge, clst$order )
        for( i in 1:length( clst$height ) ){
            text( clst$height[i], bp[i], lab.nod[i], cex=0.75, col="red", pos=4, offset=0.1 )
        }
    }
}
computeWeights <- function( X, method="flat", ..., pan.mat=t(X) ){
#   
#   Computes weights for the weighted distance metric in panGenomeTree().
#   
#   method is either "flat" (all weights 1), "zero" (ORFans weight 0, all other weights 1),
#   "shell" (shell genes emphasized) and "cloud" (cloud genes emphasized).
#   
#   Originally by Lars Snipen,
#   Biostatistics, Norwegian University of Life Sciences.
    meths <- c( "flat", "zero", "shell", "cloud" )
    midx <- grep( method, meths )
    dd <- dim( pan.mat )
    rs <- rowSums( pan.mat )
    if( midx == 1 ){
        W <- rep( 1, dd[1] )
    } else if( midx == 2 ) {
        W <- rep( 1, dd[1] )
        orfans.idx <- which( rs == 1 )
        W[orfans.idx] <- 0
    } else if( midx == 3 ) {
        W <- rep( 1, dd[1] )
        x <- 1:dd[2]
        ww <- 1 - 1/(1+exp( ((x-1)-(max(x)-1)/2)/((max(x)-1)/10) ))
        for( i in 1:dd[2] ){
            idx <- which( rs == i )
            W[idx] <- ww[i]
        }
    } else if( midx == 4 ){
        W <- rep( 1, dd[1] )
        x <- 1:dd[2]
        ww <- 1/(1+exp( ((x-1)-(max(x)-1)/2)/((max(x)-1)/10) ))
        for( i in 1:dd[2] ){
            idx <- which( rs == i )
            W[idx] <- ww[i]
        }
    } else {
        cat( "Unknown method:", method, "specified, returning flat weights\n" )
        W <- rep( 1, dd[1] )
    }
    return( W )
}
clusterSignature <- function( mergeMatrix ){
    dd <- dim( mergeMatrix )
    signature <- character( dd[1] )
    for( i in 1:dd[1] ){
        if( mergeMatrix[i,1] < 0 ){
            left <- as.character( -1*mergeMatrix[i,1] )
        } else {
            left <- gsub( ";", ",", signature[mergeMatrix[i,1]] )
        }
        if( mergeMatrix[i,2] < 0 ){
            right <- as.character( -1*mergeMatrix[i,2] )
        } else {
            right <- gsub( ";", ",", signature[mergeMatrix[i,2]] )
        }
        left <- paste( sort( unlist( strsplit( left, split="," ) ) ), collapse="," )
        right <- paste( sort( unlist( strsplit( right, split="," ) ) ), collapse="," )
        zig <- sort( c( left, right ) )
        signature[i] <- paste( zig[1], zig[2], sep=";" )
    }
    return( signature )
}
branchPos <- function( mergeMatrix, ordering ){
    dd <- dim( mergeMatrix )
    branchp <- numeric( dd[1] )
    for( i in 1:dd[1] ){
        if( mergeMatrix[i,1] < 0 ){
            left <- which( ordering == -1*mergeMatrix[i,1] )
        } else {
            left <- branchp[mergeMatrix[i,1]]
        }
        if( mergeMatrix[i,2] < 0 ){
            right <- which( ordering == -1*mergeMatrix[i,2] )
        } else {
            right <- branchp[mergeMatrix[i,2]]
        }
        branchp[i] <- mean( c( left, right ) )
    }
    return( branchp )
}
'
}
