#!/usr/bin/perl

use strict;
use Parallel::ChildManager;
use Digest::MD5;
use Getopt::Long;

my $pid = $$;

my $scratch      = "/tmp/cache-$pid-blastmatrix";
my $FORMATDB     = "/usr/biotools/blast/bin/formatdb";
my $TIGRCUT      = "/usr/biotools/indirect/tigrcut"; # puts 50% ALR and 50% cutoff identity on an XML output of blastall

my $BLASTALL     = "/usr/biotools/blast/bin/blastall"; # static version of blastall
#my $BLASTALL     = "/home/people/carsten/scripts/blastall.pl"; # static version of blastall

my $cache_source = "blastp-core+matrix";
my $ProgramName  = "coregenome.pl";
my $Version      = "3.0c";

#$ENV{BLASTMAT} = '/usr/cbs/bio/lib/ncbi/blast/matrix';

my ($keep,$slack,$cpu,$clean,$old,$wanthelp);
my $size   = 1;
my $RSpace = 10;
my $A      = 50;
my $M      = 50;
my $P      = 'blastp';

&GetOptions (
	# main options
 "a=i"        =>  \$A,          # The alignment identity used for tigr cutoff
 "clean"      =>  \$clean,      # wipe all existing results for this run
 "cpu:s"      =>  \$cpu,        # number of CPUs - 5 is default
 "keep:s"     =>  \$keep,       # keep temporary files in this dir ( TRY AND CREATED UNLESS EXIST )
 "m=i"        =>  \$M,          # The match length in percont of longest sequence for tigr cutoff
 "old"        =>  \$old,        # make old style plot
 "p=s"        =>  \$P,          # The blast program to use, default: BLASTP
 "r|rspace=i" =>  \$RSpace,     # Add extra space on the right side of the plot used for large names
 "size!"      =>  \$size,       # Enables the printing of genome size columns (default)
 "slack=i"    =>  \$slack,      # no. genomes a gene is allowed to be missing from and still be core
 "help"       =>  \$wanthelp    # print the help and exit
);

print_help() if ($wanthelp);

$cpu = 500 unless defined $cpu;
$slack = 0 unless defined $slack;

mkdir $scratch unless -d $scratch;

# create config reading from stdin
my @config = parse_config ();

@config = prepare_fasta(@config);

# make blast databases
&make_blastdbs(@config);

# make blast reports
my $cm = new ChildManager($cpu);
make_blastreports(@config);

&group(@config);

@config = core_genome ( @config ) ;
@config = pan_genome ( @config ) ;
#@config = new_genes ( @config ) ;
@config = total_families ( @config ) ;
@config = new_families ( @config ) ;

# make_table ( @config ) ;

if (defined $old) {
  plot_old(@config);
} else {
  plot(@config);
}

if ( defined ( $keep ) ) {
	system "mv $keep $keep.old_from_$$" if -e "$keep.old_from_$$";
	system "mv $scratch $keep";
	print STDERR "# kept temporary files in $keep\n";
}

cleanup ();

warn "#\n\n# ADVERTISEMENT :-)\n#   -> Try out the newer /home/people/carsten/scripts/coregenome4.pl using MCL. \n#\tRuns on the CGE cluster with faster alignments, better clustering, and many more options for customization.# -Have fun!\n\n";

sub cleanup {
	system "rm -rf $scratch";
}

sub make_blastreports {
	my @config = @_;
	my %tab_h;
	foreach my $a ( 0 .. $#config ) {
		foreach my $b ( 0 .. $#config ) {
			$tab_h{"$a-$b"} = $cm->start();
			unless ( $tab_h{"$a-$b"} ) {
				my $blastreport_scratch = "$scratch/$a-$b.blastout.gz";
				print STDERR "# make blast report $blastreport_scratch\n";
				my $jobid = md5 ( "$scratch/$a.fsa" , "$scratch/$b.fsa" ) ;
                                $jobid .= "a".$A.$M;               # Need $A and $M to ensure that the cutoff is respected in the cache...
				system "$perl /usr/biotools/indirect/cacher --id='$jobid' --source='$cache_source' -action get > $blastreport_scratch";
				if ( $? != 0 or $clean ) {
					print STDERR "# jobid $jobid not in cache - redoing\n";
					system "$BLASTALL -F 0 -i $scratch/$a.fsa -p $P -e 1e-5 -d $scratch/$b.fsa | $TIGRCUT -a $A -m $M | gawk '{print \$1\"\\t\"\$2}' | gzip > $blastreport_scratch  \n";
					system "$perl /usr/biotools/indirect/cacher --id='$jobid' --source='$cache_source' -action put -expire 100 < $blastreport_scratch";
				} else {
					print STDERR "# fetched jobid $jobid from cache\n";
				}
				my $waiting = ($#config+1)**2;
				foreach my $a ( 0 .. $#config ) {
					foreach my $b ( 0 .. $#config ) {
						$waiting-- if -e "$scratch/$a-$b.blastout.gz";
					}
				}
				print STDERR "# jobid $jobid finished - $waiting left\n";
				exit;
			}
		}
	}
	$cm->wait_all_children; 
}

sub make_blastdbs {
	my @config = @_;
	foreach my $id ( 0 .. $#config ) {
		my $target = $config[$id]->{target};
		my $p = ($P =~ /blastn|tblastn|tblastx/i ? "F" : "T");
		print STDERR "# building blast database '$target', Protein=$p\n";
		system "$FORMATDB -i $target -p $p -t $id";
	}
}

sub prepare_fasta {
	my @config = @_;
	foreach my $id ( 0 .. $#config ) {
		my ($nprot,$skipped,$checksum) = fasta2fasta( $config[$id]->{source} , "$scratch/$id.fsa");
		warn "# prepare_fasta(): nprot=$nprot, skipped=$skipped\n";
		$config[$id]->{target} = "$scratch/$id.fsa";
		$config[$id]->{'total genes'} = $nprot;
		$config[$id]->{skipped} = $skipped;
		$config[$id]->{checksum} = $checksum;
		die "nprot == 0 for $config[$id]->{source}\n" if $nprot == 0;
	}
	return @config;
}

sub fasta2fasta {
	my ($input,$output) = @_;
	my $temp = "$output.temp.raw";
	my $nprot = 0;
	my $skipped = 0;
	open RAW , "| saco_convert -I fasta -O raw | sort > $temp" or die $!;
	open FASTA , $input;
	while (<FASTA>) {
		print RAW;
	}
	close FASTA;
	close RAW;
	my $checksum = md5 ( $temp );
	my %ids;
	for (`saco_convert -I fasta -O tab $input`) {
		my ($id, $seq) = split /\t|\n/, $_;
		push @{ $ids{$seq} }, $id;
	}
	open FASTA , "| saco_convert -I tab -O fasta > $output";
	open TAB , $temp or die $!;
	while (<TAB>) {
		chomp;
		if ( ! /^([A-Za-z]+)/ ) {
			$skipped++;
			next;
		}
		$nprot++;
		my $id = shift @{ $ids{$1} };
		print FASTA "$checksum.$nprot\t$1\t\t$id /sourcefile=\"$input\"\n";
	}
	close TAB;
	close FASTA;
	unlink $temp;
	return ($nprot,$skipped,$checksum);
}

sub md5 {
	my @files = @_;
	my $md5 = Digest::MD5->new;
	foreach my $file (@files) {
		open(FILE, $file) or die "Can't open '$file': $!";
		binmode(FILE);
		while (<FILE>) {
			$md5->add($_);
		}
		close(FILE);
	}
	return $md5->hexdigest;
}

sub parse_config {
	my @config;
	while ( <> ) {
		next if /^#/;
		chomp;
		my ($description,$source,$tag) = split /\t/;
		warn "# description=$description, source=$source\n";
		my $id = $#config + 1;
		( $config[$id]->{description} , $config[$id]->{source} , 
                  $config[$id]->{id} , $config[$id]->{tag}) = ( $description , $source , $id , $tag);
	}
	return @config;
}

sub group {
	my @config = @_;
	foreach my $id (0 .. $#config) {
	print STDERR "# Building group $id of $#config\n";
		open GRP , "| perl /home/people/pfh/scripts/group/group > $scratch/group_$id.dat";
		foreach my $x (0 .. $id) {
			foreach my $y (0 .. $id) {
				foreach my $z ($x , $y) {
					foreach my $n ( 1 .. $config[$z]->{'total genes'} ) {
						print GRP $config[$z]->{checksum}.".$n\n";
					}
				}
				open BLASTOUT , "gunzip -c $scratch/$x-$y.blastout.gz |";
				while (<BLASTOUT>) {
					print GRP $_;
				}
				close BLASTOUT;
			}
		}
		close GRP;
	}
}

sub pan_genome {
	my @config = @_;
	foreach my $id (0 .. $#config) {
		$config[$id]->{'pan genome'} = 0;
		open GRP , "$scratch/group_$id.dat";
		while (<GRP>) {
			$config[$id]->{'pan genome'}++;
		}
		close GRP;
	}
	return @config;
}

sub core_genome {
	my @config = @_;
	foreach my $id (0 .. $#config) {
		$config[$id]->{'core genome'} = 0;
		open GRP , "$scratch/group_$id.dat";
		while (<GRP>) {
			my ($group_id,$count,@members) = split /\t/;
			my %repr;
			foreach (@members) {
				next unless /^(.*)\.(.*)/;
				$repr{$1} = 1;
			}
			$config[$id]->{'core genome'}++ if (scalar ( keys %repr ) <= $id + 1 and
                                                            scalar ( keys %repr ) >= $id + 1 - $slack);
		}
		close GRP;
	}
	return @config;
}

sub total_families {
	my @config = @_;
	foreach my $id (0 .. $#config) {
		$config[$id]->{'total families'} = 0;
		open GRP , "$scratch/group_$id.dat";
		while (<GRP>) {
			my ($group_id,$count,@members) = split /\t/;
			my $is_member = 0;
			foreach (@members) {
				next unless /^(.*)\.(.*)/;
				$is_member = 1 if $1 eq $config[$id]->{checksum};
			}
			$config[$id]->{'total families'} += $is_member;
		}
		close GRP;
	}
	return @config;
}

sub new_genes {
	my @config = @_;
	foreach my $a ( 0 .. $#config ) {
		my %HITS;
		$config[$a]->{'new genes'} = 0;
		print STDERR "# parsing blast reports (id $a)\n";
		foreach my $b ( 0 .. $a ) {
			open GUNZIP , "gunzip -c $scratch/$a-$b.blastout.gz $scratch/$b-$a.blastout.gz |";
			while (<GUNZIP>) {
				chomp;
				my ($q,$s) = split /\t/;
				$HITS{$q}{$s} = 1;
			}
			close GUNZIP;
		}
		open FSA , "$scratch/$a.fsa" or die $!;
		while (<FSA>) {
			chomp;
			next unless /^>(.*)/;
			my $q = $1;
			my $hits_in_other_samples;
			foreach (keys %{$HITS{$q}}) {
				next unless /^(.*)\./; 
				if ( $config[$a]->{checksum} ne $1) {
					$hits_in_other_samples = 1;
					last;
				}
			}
			$config[$a]->{'new genes'}++ if ! defined $hits_in_other_samples;
		}
	}
	return @config;
}

sub new_families {
	my @config = @_;
	foreach my $id (0 .. $#config) {
		$config[$id]->{'new families'} = 0;
		open GRP , "$scratch/group_$id.dat";
		while (<GRP>) {
			my ($group_id,$count,@members) = split /\t/;
			my $is_new_family = 1;
			foreach (@members) {
				next unless /^(.*)\.(.*)/;
				$is_new_family = 0 if $1 ne $config[$id]->{checksum};
			}
			$config[$id]->{'new families'} += $is_new_family;
		}
		close GRP;
	}
	return @config;
}

sub plot {
	my @config = @_;
	my $tbl = "$scratch/tbl";
	open TBL , ">$tbl";
	my @keys = ("id", "description", "total genes", "total families", "new families",
                    "pan genome", "core genome", "tag");
	print TBL join ("\t",@keys)."\n";
	foreach my $id ( 0 .. $#config ) {
		foreach my $key (@keys) {
			if ($id >= $slack or $key ne "core genome") {
				print TBL $config[$id]->{$key},"\t";
			} else {
				print TBL "NA","\t";
			}
		}
	print TBL "\n";
	}
	close TBL;
	open R , "| /usr/bin/R --vanilla > /dev/null";
	print R "
postscript('$scratch/ps', title='$ProgramName v$Version');
data <- read.table('$tbl',skip=1,sep='\t',dec='.',header=FALSE);
layout( matrix(c(1,2), nrow= 2), heights=c(10,3))
op <- par(mar = c(0,2,2,1))
x<-rbind(data[,3], data[,4], data[,5])
data[which.max(data[,7]),7] <- data[1,5]
ymax <- max(1.2*max(data[,6]), 2.0 * max(data[,3]))
col <- c(gray(0.8), gray(0.5), gray(0.2))
";
	if ($size) {
	print R "
rspace <- $RSpace * ncol(x) / 100
r <- barplot(x, beside=TRUE, plot=F)
xlim <- c(min(r)-0.5, max(r)+0.5+rspace)
r <- barplot(x, beside=TRUE, ylim=c(0,ymax), xlim=xlim, col=col)
lines(r[3,1:(ncol(r)-sum(is.na(data[,7])))], type='b', data[!is.na(data[,7]),7], col='red', lwd=4)
lines(r[3,], type='b', data[,6], col='blue', lwd=4)
legend(1,ymax,c('Total Genes', 'Total Gene Families', 'New Gene Families','Core Genome','Pan Genome'),col=c(col,'red','blue'), lwd=c(4,4,4,4))
	"} else {
	print R "
rspace <- $RSpace * 1 / 100
r <- barplot(x[3,], beside=TRUE, plot=F)
xlim <- c(min(r)-0.5, max(r)+0.5+rspace)
r <- barplot(x[3,], beside=TRUE, ylim=c(0,ymax), xlim=xlim, col=col[3])
lines(r[1:(length(r)-sum(is.na(data[,7])))], type='b', data[!is.na(data[,7]),7], col='red', lwd=4)
lines(r, type='b', data[,6], col='blue', lwd=4)
legend(1,ymax,c('New Gene Families','Core Genome','Pan Genome'),col=c(col[3],'red','blue'), lwd=c(4,4,4))
r <- t(r)
"}
	print R "
op <- par(mar = c(1,2,0,1))
plot.new()
plot.window(ylim=c(0,10),xlim=xlim)
for (i in 1:nrow(data)) {
  text(mean(r[,i]), 9, adj=0, data[i,2], cex=0.7, srt=-45)
}
dev.off();
";
	close R;
	system "cat $scratch/ps";
}

sub plot_old {
	my @config = @_;
	my $tbl = "$scratch/tbl";
	open TBL , ">$tbl";
	my @keys = ("id","description","total genes","new genes","new families","pan genome","core genome");
	print TBL join ("\t",@keys)."\n";
	foreach my $id ( 0 .. $#config ) {
		foreach my $key (@keys) {
			print TBL $config[$id]->{$key},"\t";
		}
	print TBL "\n";
	}
	close TBL;
	open R , "| /usr/bin/R --vanilla > /dev/null";
	print R "
postscript('$scratch/ps', title='$ProgramName v$Version');
layout( matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(0.6,0.4))
data <- read.table('$tbl',skip=1,sep='\t',dec='.',header=FALSE);
tN <- data[,3]
x<-rbind(data[,4],data[,5])
r <- barplot(x,beside=TRUE,cex.names=0.7,names.arg=as.vector(data[,1])+1,ylim=c(0,1.3*max(data[,6])))
lines(r[2,(1:nrow(data))[!is.na(data[,7])]], type='b', data[!is.na(data[,7]),7], col='red', lwd=4)
lines(r[2,], type='b', data[,6], col='blue', lwd=4)
legend(1,1.3*max(data[,6]),c('New genes','New gene families','Core genome','Pan genome'),col=c('darkgray','lightgray','red','blue'), lwd=c(4,4,4,4))
plot.new()
plot.window(xlim=c(0,2),ylim=c(-nrow(data),1))
for (i in 1:nrow(data)) {
 text(0,-i, adj=0,paste(i,': ',data[i,2]),cex=0.7)
}
dev.off();
";
	close R;
	system "cat $scratch/ps";
}


sub print_help {
  print <<EOH;

NAME
    $ProgramName - derive core-/pan genome as well as count new genes/families
    in a list of genomes/samples.


SYNOPSIS
    $ProgramName [options] < list > output.ps

DESCRIPTION
    coregenome reads a list from STDIN or a file, defining the proteomes to be
    compared. The order of the list indicate in what order the genomes will
    appear in the plot. Each line must contain a description and data source,
    separated by tab. Example:

    Campylobacter jejuni subsp. jejuni NCTC 11168	saco_extract -t -I genbank -O fasta /home/projects/pfh/genomes/data/AL111168/AL111168.gbk |
    Campylobacter jejuni RM1221	saco_extract -t -I genbank -O fasta /home/projects/pfh/genomes/data/CP000025/CP000025.gbk |
    Campylobacter fetus subsp. fetus 82-40	saco_extract -t -I genbank -O fasta /home/projects/pfh/genomes/data/CP000487/CP000487.gbk |
    Campylobacter jejuni subsp. jejuni 81-176	saco_extract -t -I genbank -O fasta /home/projects/pfh/genomes/data/CP000538/CP000538.gbk |
    Campylobacter jejuni subsp. doylei 269.97	saco_extract -t -I genbank -O fasta /home/projects/pfh/genomes/data/CP000768/CP000768.gbk |
    Campylobacter curvus 525.92	saco_extract -t -I genbank -O fasta /home/projects/pfh/genomes/data/CP000767/CP000767.gbk |
    Campylobacter hominis ATCC BAA-381	saco_extract -t -I genbank -O fasta /home/projects/pfh/genomes/data/CP000776/CP000776.gbk |
    Campylobacter concisus 13826	saco_extract -t -I genbank -O fasta /home/projects/pfh/genomes/data/CP000792/CP000792.gbk |
    Campylobacter jejuni subsp. jejuni 81116	saco_extract -t -I genbank -O fasta /home/projects/pfh/genomes/data/CP000814/CP000814.gbk |

    Note that the above example does not specify files as the data source,
    rather a command to obtain the data (tailing pipe is important). It is also
    possible to just specify the files if these already exist.

    Each proteome is blasted against all previous proteoms in the list. All
    blast results are kept in the directory './cache', and the naming is
    derived based in MD5 checksums of the raw proteins sequences. Example:

    GenomeA	A.fsa
    GenomeB	B.fsa
    GenomeC	C.fsa

		Providing the list as above will cause all proteoms to be BLAST againast all.

    The program will for each proteome X calculate the follwing:

    - Number of proteins observed in X
    - Number of new genes observed in genome in X compared to 0 .. (X-1) 
    - Number of new gene families observed in genome in X compared to 0 .. (X-1) 
    - Size of pan genome at genome X
    - Size of core genome at genome X

    The program assumes two proteins being similar if the aligment identified by
    BLASTP (evalue <= 1e-5) spans 50% or more of the longest of the sequence AND
    that 50% of the residues are conserved within the alignment. Only the reciprocal match
    is considered.

OPTIONS
    -a <N>
	Set the required alignment identity in percent for two sequences to be
	joined into one family. Also see '-m'.
	Default is: $A

    -clean
	Ignores existing and cached results and regenerates the blast reports.

    -cpu <N>
	Specify how many BLAST searches run in parallel.

    -keep <S>
	Keep temporary files in this dir ( TRY AND CREATED UNLESS EXIST )

    -m <N>
	Set the required match length in percent of the longest participating
	sequence in each alignment for the two to be joined into one family.
	Also see '-a'.
	Default is: $M

    -slack <N>
	Number of genomes a gene is allowed to be missing from and still be
	considered part of the core genome. Notice that the core genome curve
	will not reach the rightmost genomes when using this option. This is
	natural, because the core genome of species 'n' now requires data from
	genome 'n+slack', which is unavailable for the rightmost species.

    -nosize
	Suppreses the printing of the Genome Size columns.

    -r or --rspace [Integer]
	Add empty space to the right of the plot. Increase if the column names
	don't fit in the plotting window. Given as a relative value roughly
	equal to a percentage of the (original) plotting window.
	Default is $RSpace.

    -old
	Make old-style plot.
 
VERSION
    Current Version: $ProgramName $Version

    Version Changes:

    Version 3.0c has been updated to use the sbiology cueing system. This also
    changes the behaviour of the '-cpu' option which now de-facto controls the
    amount of simultaneous submissions to the cue. Thus this number should now
    generally be higher than for previous versions else "submission-lag" will be
    a factor. New '-cpu' default has been set.
    New options: '-a' and '-m' options to finetune the cutoff

    Version 2.8c includes smaller bugfixes and adjusts composition of temporary
    files to better integrate with the coregenes script.

    Version 2.7c added the 'rspace' option as well as the 'total gene families'
    column. Also replaced 'saco_convert' with 'tab2fasta' for preparing files
    for blast. 'saco_convert' has problems with long identifiers which can be
    disasterous.

    Version 2.6c mainly adjusts the composition of temporary files to better
    integrate with the coregenes script for extracting individual genes.

    Version 2.5c introduced the '-nosize' option.

    Version 2.4c introduced new options and fixed the bug which prevented the
    script from using newer versions of BLAST. As a consequence it should run
    faster.

EXAMPLE
    mysql -B -N -e "select genbank,concat('saco_extract -t -I genbank -O fasta \
      /home/projects/pfh/genomes/data/',genbank,'/',genbank,'.gbk |') from pfh_public.genbank_complete_seq as s , \
      pfh_public.genbank_complete_prj as p where s.pid = p.pid and p.organism_name like 'campy%' \
      and genbank not like 'genome%' and segment_name like 'chromo%' order by released" | \
      perl coregenome.pl -cpu 10 > output.ps


AUTHOR
    Carsten Friis, carsten\@cbs.dtu.dk,
    Based on script made by Peter Fischer Hallin, pfh\@cbs.dtu.dk

EOH
  exit;
}
