#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ( $input_file, $help, $output_dir, $fasta_flag, @selected_profiles );

#-------------------------------------------------------------------------
#Get input options

&GetOptions(
    "i=s"    => \$input_file,
    "h"      => \$help,
    "o=s"    => \$output_dir,
    "fasta!" => \$fasta_flag,         #"fasta:o"=> \$fasta_flag,	#boolean
    "p:s"    => \@selected_profiles
);

#-------------------------------------------------------------------------
#Timestamp
my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst )
    = localtime();

my $timestamp = sprintf("%d%02d%02d_%02d%02d%02d", $year + 1900, ++$mon, $mday, $hour, $min, $sec);

#-------------------------------------------------------------------------
#Deal with input options

if ($help) {
    help();
}

if ( ( !$input_file ) or ( !$output_dir ) ) {
    print "Unknown option, try running $0 -h for more information\n";
}

unless ( -e $input_file ) {
    print "$input_file not found\n";
    help();
}

validate_fasta();

unless ( defined $fasta_flag ) {
    $fasta_flag = 1;
}

unless (@selected_profiles) {
    @selected_profiles
        = ( "Galpha", "Gs", "Gio", "Gq11", "G1213", "Gbeta", "Ggamma" );
}

system("mkdir -p $output_dir/$timestamp");
print STDOUT "Running hmmsearch...\n";

@selected_profiles = split( /,/, join( ',', @selected_profiles ) );

foreach my $profile (@selected_profiles) {

    system(
        "hmmsearch --cut_tc profiles/$profile.hmm $input_file >$output_dir/$timestamp/$profile.res"
    );
    print STDOUT "Storing hmmsearch output for $profile profile...\n";
    create_summary( "$output_dir/$timestamp/$profile.res", $profile );

}

print STDOUT "Summary file created...\n";

#print $fasta_flag."\n";
if ($fasta_flag) {
    print STDOUT "Creating fasta output files...\n";
    system("mkdir -p $output_dir/$timestamp/fasta_output");

    foreach my $profile (@selected_profiles) {

        #print $profile."\n";
        create_fasta(
            $profile, $input_file,
            "$output_dir/$timestamp/summary.txt",
            "$output_dir/$timestamp/fasta_output/$profile.fa"
        );
    }
}

print STDOUT "Process completed successfully!\n";

#-------------------------------------------------------------------------
# Subroutines

sub create_summary {
    my ( $res_file, $profile ) = @_;

    open my $summary_fh, ">>", dirname($res_file) . "/summary.txt";

    print $summary_fh
        "List of predicted $profile proteins (by Pfam database model):\n\n"
        . "E-value  Score  Sequence  Model Start  Model End  Alignment Start  Alignment End\n"
        . "-------- ------  --------  -----------  ---------  ---------------  -------------\n";

    open my $res_fh, "<", $res_file;

    my $counter = 0;
    my @general_info;

    while (<$res_fh>) {
        if ( $_ =~ /^\s\s\s(.{7}\d)\s([\s\.\d]{6})[e\-\s\.\d]{42}(.*)/ ) {
            chomp $3;
            push( @general_info, "$1 $2  $3  " );
        }

        if ( $_ =~ /No hits detected that satisfy reporting thresholds/ ) {
            print $summary_fh "NONE\n";
        }

    }
    close $res_fh;

    open $res_fh, "<", $res_file;

    local $/ = "\n>>";
    my @align_info;

    while (<$res_fh>) {
        my $align_info = "";
        my $flag       = 0;

        if ( $_
            =~ /No individual domains that satisfy reporting thresholds (although complete target did)/
            )
        {
            push( @align_info, "no align info" );
            next;
        }
        while ( $_
            =~ /\s\s\s\d\s!.*\s+(\d+)\s+(\d+)\s[\.\[\]]{2}\s+(\d+)\s+(\d+)\s[\.\[\]]{2}\s+(\d+)\s+\d+\s[\.\[\]]{2}/g
            )
        {
            $align_info .= "$1	$2		$3	$4;";
            $flag = 1;
        }
        if ( $flag == 1 ) {
            $align_info .= " $profile";
            push( @align_info, $align_info );
            $counter++;
        }
    }

    close $res_fh;

    if ( (@general_info) && (@align_info) ) {
        for ( my $i = 0; $i <= $#general_info; $i++ ) {
            if ( $align_info[$i] ne "no align info" ) {
                print $summary_fh "$general_info[$i] --- $align_info[$i]\n";
            }
        }
    }

    print $summary_fh "\n$counter hit(s).\n\n";
    print $summary_fh
        "********************************************************************************\n";
    print $summary_fh
        "********************************************************************************\n";

}

sub create_fasta {
    my ( $profile, $input_file, $summary_file, $output_fasta ) = @_;

    open my $summary_fh, "<", $summary_file;
    open my $fasta_fh,   ">", $output_fasta;

    my @fam;

    while (<$summary_fh>) {
        while ( $_ =~ /^.{14}\d\s\s(\S+)\s.*\s$profile$/mg ) {

            push( @fam, $1 );
        }
    }

    close $summary_fh;

    my $pattern;
    for ( my $i = 0; $i <= $#fam; $i++ ) {
        $pattern = $fam[$i];
        open my $input_fh, "<", $input_file;

        $/ = ">";
        while (<$input_fh>) {
            if ( $_ =~ /\Q$pattern\E/ ) {
                chomp $_;
                print $fasta_fh ">$_";
                last;
            }

        }
        close $input_fh;

    }
    close $fasta_fh;
}

sub validate_fasta {
    my $valid_protein_seq_characters = "A-Z";

    open my $checkfile, '<', $input_file;

    my $counter = 0;
    my $error   = 0;

    while (<$checkfile>) {
        $counter++;
        if ( $_ =~ /^>\S.*$/ ) {
            $error = 0;
        }
        elsif ($_ =~ /^[$valid_protein_seq_characters]+\n$/g
            || $_ =~ /^\n/ )
        {
            $error = 0;

        }
        else {
            $error = 1;
            last;
        }

    }

    close $checkfile;

    if ($error) {
        print "ERROR: Your file is not in fasta format. Check line $counter"
            . "\n";
        exit;
    }

}

sub help {
    print STDERR << "EOF";
    
	GprotPRED is a tool for accurate detection and classification of G-proteins. 
	It uses profile Hidden Markov Models (pHMMs) for the four known heterotrimeric Galpha protein families, 
	the Gbeta and the Ggamma subunit in order to classify a set of protein sequences 
	into the appropriate G-protein family.
	This is a standalone version of GprotPRED tool for local execution.


	Usage:

	$0 -i <input_fasta_file> -o <output_directory> -nofasta <no_fasta_output> -p <selected_profiles> 

	i | input_fasta_file		Input file that contains the query protein sequences in fasta format. Required

	o | output_directory		The directory where the result output files will be stored. Required

	nofasta | no_fasta_output	Use this option if you don't want fasta output. Optional, default: fasta files of the predicted proteins are produced

	p | selected_profiles		Use selected profiles only, i.e. -p Gs,Gio. Optional, default: All profiles are used (Galpha,Gs,Gio,Gq11,G1213,Gbeta,Ggamma)
				

EOF

    exit(0);
}

