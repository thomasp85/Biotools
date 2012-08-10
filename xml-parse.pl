#!/usr/bin/perl

use Getopt::Long;

&GetOptions 	(
				"alrmethod:s", \$alrmethod,
				);
###############################################################################################################################
# This script reads an XML blast report from STDIN. The report must be formated using 'blastall -m7'!
# You can extract whatever information that is stored within the following structure: m/<.*><\/<.*>/
# When a </Hit> is observed, all retrieved data are printet to STDOUT as tab separated file.
# ARGV requires a string, containing the tags that must be printed on each line
#
# EXAMPLE: "Hsp_bit-score,Hsp_query-to,Hsp_positive,H"
# 
# Alternativley, standard outputs can be specified, 1,2,3... 
###############################################################################################################################


# we check to see, if the user specified a standard output format - otherwise,
# use the custom model, specified by the user.

%stdModels = ('std' => "BlastOutput_query-def,Hit_def,Hsp_identity,Hsp_align-len,Hsp_gaps,Hsp_query-from,Hsp_query-to,Hsp_hit-from,Hsp_hit-to,Hsp_evalue,Hsp_score",
              'std-with-length' => "BlastOutput_query-def,Hit_def,Hsp_identity,Hsp_align-len,BlastOutput_query-len,Hit_len,Hsp_gaps,Hsp_query-from,Hsp_query-to,Hsp_hit-from,Hsp_hit-to,Hsp_evalue,Hsp_score",
              'std-with-length-alr' => "BlastOutput_query-def,Hit_def,Hsp_identity,Hsp_align-len,BlastOutput_query-len,Hit_len,ALR,Hsp_gaps,Hsp_query-from,Hsp_query-to,Hsp_hit-from,Hsp_hit-to,Hsp_evalue,Hsp_score",
              'std-with-length-alr-seq' => "BlastOutput_query-def,Hit_def,Hsp_identity,Hsp_align-len,BlastOutput_query-len,Hit_len,ALR,Hsp_gaps,Hsp_query-from,Hsp_query-to,Hsp_hit-from,Hsp_hit-to,Hsp_evalue,Hsp_score,Hsp_hseq");

$model = $stdModels{$ARGV[0]}  if (defined($stdModels{$ARGV[0]}));
$model = $ARGV[0] unless (defined($stdModels{$ARGV[0]}));

if (! defined $ARGV[0] ) 
    {
    warn ("No output model specified\nYou can choose between one of the following standard output formats:\n");
    foreach $key (keys %stdModels)
        {
        print "\n$key:\n\t$stdModels{$key}\n";
        }
    die("\n");
    }

@print = split(/,/,$model);

print "#".join("\t",@print)."\n";

$line = 1;
my $ret = 1;
while ( defined ( $line = <STDIN> ) ) {
		$ret = 0;
		# capture END
    if ($line =~ m/<\/Hsp>/) {
        if ($model =~ m/ALR/)
            {
			$qlen = $hash{'BlastOutput_query-len'}+0;
			$slen = $hash{'Hit_len'}+0;
            if ($qlen > $slen || $alrmethod eq 'query')
	        	{
	        	$hash{"ALR"} = sprintf("%0.2f",$hash{'Hsp_align-len'} / $qlen );
	            }
	        elsif($qlen < $slen || $alrmethod eq 'subject')
	        	{
	        	$hash{"ALR"} = sprintf("%0.2f",$hash{'Hsp_align-len'} / $slen );
	            }
	        elsif($alrmethod eq '')
	        	{
	        	$hash{"ALR"} = sprintf("%0.2f",$hash{'Hsp_align-len'} / $qlen ) if ($qlen > $slen );
	        	$hash{"ALR"} = sprintf("%0.2f",$hash{'Hsp_align-len'} / $slen ) if ($qlen <= $slen );
	            }
            }
        foreach $key (@print) {
            $hash{$key} = 0 unless defined ($hash{$key});
            print $hash{$key}."\t"; }
        print "\n";
        undef($hash); }
    if ($line =~ m/<(.*)>(.*)<\/.*>/) {
       $hash{$1} = $2; } }
exit $ret;