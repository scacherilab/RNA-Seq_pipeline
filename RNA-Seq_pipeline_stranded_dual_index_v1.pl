#!/usr/bin/perl -w
BEGIN {$^W=0}
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
=head1 NAME

    Using RNA-Seq-pipeline.pl

    
=head1 SYNOPSIS

     Example: ./RNA-Seq_pipeline.pl -i sample_adapter_mapping_file.txt -bwindex /home1/genomes/hg18/hg18 < other options ...>
     
     Required Parameters:
    
        -infile,-i          =   (REQUIRED) File listing all input FASTQs and corresponding adapters (if applicable), one sample per line. Paired end samples should have 2 files
                                listed per line, delimited by a comma.
        
            Example file if samples need to be adapter-trimmed:
            
            full/path/to/untrimmed_sample1_single_end.fastq <TAB(s) or SPACE(s)> ADAPTER1_int
            full/path/to/untrimmed_sample2_paired_end_1.fastq,full/path/to/untrimmed_sample2_paired_end_2.fastq <TAB(s) or SPACE(s)> ADAPTER2_int,ADAPTER3_int
            full/path/to/untrimmed_sample3_single_end.fastq <TAB(s) or SPACE(s)> ADAPTER3_int
            ...
            
            Example file if samples are already adapter-trimmed (must also enable -noatrim):
            
            full/path/to/trimmed_sample1_single_end.fastq <additional columns that will be ignored>
            full/path/to/trimmed_sample2_paired_end_1.fastq,full/path/to/trimmed_sample2_paired_end_2.fastq  <additional columns that will be ignored>
            full/path/to/trimmed_sample3_single_end.fastq <additional columns that will be ignored>
            
            
        -bwindex         =   (REQUIRED) Path to bowtie index (example: /home1/genomes/hg18/hg18 where the second hg18 is the basename for .ebwt files hg18.1.ebwt, hg18.2.ebwt, etc. )
        
        -path               =  (REQUIRED only if running outside evolution/genomecruncher) path to directory containing merge.plx and
                                the following working executables: cufflinks,tophat,cuffcompare [ default on evolution/genomecruncher: /homeG/ScacheriLab/ars51/rna_seq_tools ]
                                To run the pipeline from HPCC, set the path to /home/ars51/programs/rna_seq_tools_for_pipeline
                                Note: fastx_clipper and fastq_quality_trimmer must also be installed and accessible through the PATH variable on your system.
        
        
        Other Options:
    
        -h, -help, ?        =   Prints this message
        -noatrim            =  Do not perform adapter trimming
        -noqtrim            =  Do not perform quality trimming
        -minlen, -l         =  (will be ignored if -noatrim and -noqtrim are both set) Length threshold - sequences shorter than this (after trimming) will be discarded. [default: 25]
        -minqual, -q        =  (will be ignored if -noqtrim is set) Quality threshold - nucleotides with lower quality will be trimmed (from the end of the sequence). [default: 20]
        -hits               =   Max multiple hits (tophat) [default: 10]
        -ilen               =   Max intron length (cufflinks)- ignore alignments with gaps longer than this [default: 20000]
	-seglen		    =	Segment length [default: 25]
        -r                  =   Mate inner distance (paired end samples only for tophat) [default: 60]
        -thread, -t         =   Number of cores to run tophat and cufflinks on [default: 1]
        -cuffref            =   Reference transcript annotations file (gtf reference file for running cufflinks, cuffcompare or cuffidff) [default: none]
        -genome, -g         =   Path to reference genome file for cufflinks to use bias correction(/home1/genomes/hg18/hg18.fa)
        
        To run cuffcompare , set this flag (must also set -cuffref):
        
        -cuffcomp           =   Run cuffcompare on each sample
       
=head1 DESCRIPTION
        
        This Perl script is designed to run the RNA-Seq pipeline  on at least 1 sample. It requires a path to the reference genome (ex. /home1/genomes/hg18/hg18 where the second hg18 is a prefix) 
        and file listing all input FASTQs and corresponding adapters (if applicable), one sample per line. Paired end samples should have 2 files listed per line, delimited by a comma:
        
            full/path/to/untrimmed_sample1_single_end.fastq <TAB(s) or SPACE(s)> ADAPTER1_int
            full/path/to/untrimmed_sample2_paired_end_1.fastq,full/path/to/untrimmed_sample2_paired_end_2.fastq <TAB(s) or SPACE(s)> ADAPTER2_int,ADAPTER3_int
            full/path/to/untrimmed_sample3_single_end.fastq <TAB(s) or SPACE(s)> ADAPTER3_int
        
        The adapter integers correspond to standard Illumina adapters.
        
                    '2'=> 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG'
                   '4'=> 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG'
                   '5'=> 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG'
                   '6'=> 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG'
                   '7'=> 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG'
                   '12'=> 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG'

        
        The steps performed during the execution of the pipeline are:
        1) Adapter trimming using standard Illumina adapters (requires fastx-clipper of the Fastx Toolkit)
        2) Quality trimming (requires fastq_quality_trimmer of the Fastx Toolkit)
        3) Running Tophat to align RNA-Seq reads to a reference genome and analyze the mapping results to identify splice junctions between exons (requires Tophat)
        4) Running Cufflinks to assemble transcripts and estimate their abundances (requires Cufflinks)
        5) (Optional step) Running Cuffcompare on each sample to compare your cufflinks-assembled transcripts to a reference annotation
        6) Rounding, flooring and merging (if applicable) FPKM expression values (requires merge.plx script)
        
        Final processed FPKM values will be stored in  fpkm_files/*fpkm_tomerge and fpkm_files/fpkm_all_samples_rounded.txt (if applicable)
        Cuffcompare output files will be stored in cuffcompare_files.
        The tophat and cufflinks output files will be stored in the individual sample fastq output directories.
        For a paired end sample, the tophat and cufflinks output files will be stored in the first of the output fastq directories for that sample (Ex. untrimmed_sample2_paired_end_1).
        
        NOTE: To test for differential expression in RNA-Seq samples, run cuffdiff. See http://cufflinks.cbcb.umd.edu/manual.html#cuffdiff for more information
     in RNA-Seq samples, run cuffdiff. See http://cufflinks.cbcb.umd.edu/manual.html#cuffdiff for more information
    
=cut

### internal variables
my $test_mode=0; ### set 1 if testing only, otherwise 0. 

###Standard Illumina Adapters map###


my %adapters_map=('704' => 'AATGATACGGCGACCACCGAGATCTACACGAGATTCCTCGTCGGCAGCGTC',
'705' => 'AATGATACGGCGACCACCGAGATCTACACATTCAGAATCGTCGGCAGCGTC',
'706' => 'AATGATACGGCGACCACCGAGATCTACACGAATTCGTTCGTCGGCAGCGTC',
'501' => 'CAAGCAGAAGACGGCATACGAGATTATAGCCTGTCTCGTGGGCTCGG',
'502' => 'CAAGCAGAAGACGGCATACGAGATATAGAGGCGTCTCGTGGGCTCGG',
'503' => 'CAAGCAGAAGACGGCATACGAGATCCTATCCTGTCTCGTGGGCTCGG',
'504' => 'CAAGCAGAAGACGGCATACGAGATGGCTCTGAGTCTCGTGGGCTCGG',
'505' => 'CAAGCAGAAGACGGCATACGAGATAGGCGAAGGTCTCGTGGGCTCGG',
'506' => 'CAAGCAGAAGACGGCATACGAGATTAATCTTAGTCTCGTGGGCTCGG',
'507' => 'CAAGCAGAAGACGGCATACGAGATCAGGACGTGTCTCGTGGGCTCGG',
'508'=> 'CAAGCAGAAGACGGCATACGAGATGTACTGACGTCTCGTGGGCTCGG'
);

###Parameters
my %h_options = ('minlen'=>25,'minqual'=>20,'hits'=>10,'ilen'=>20000,'r'=>60,'thread'=>1,'c'=>10,'path'=>"/homeG/ScacheriLab/ars51/rna_seq_tools",'seglen'=>25);
GetOptions (\%h_options,'infile|i=s', 'minlen|l=i','minqual|q=i' ,'hits=i','genome|g=s','ilen=i', 'r=i','help|h|?','thread|t=i','cuffref=s','c=i','cuffcomp','cuffdiff=s','noatrim','noqtrim','path=s','bwindex=s','seglen=i');


##Other stuff
my $path_to_rna_seq_tools=$h_options{'path'}; ##this is where the most recent cufflinks and tophat files are
my $sample_adapter_map_file=$h_options{'infile'};
my $bwindex=$h_options{'bwindex'};
my $genome=$h_options{'genome'};
my $min_seq_len=$h_options{'minlen'};
my $num_cores=$h_options{'thread'};
my $valid_input=$input_line;
my %file_adapter_prefixes_map;

if($h_options{'help'}){pod2usage( -verbose => 2);}
elsif($sample_adapter_map_file eq "" || $bwindex eq "" ){ pod2usage()}
else{ ##do various initial pre-pipeline checks
   
print "Validating input file\n";
($valid_input, $input_line, %file_adapter_prefixes_map)=validateSampleAdapterFile($sample_adapter_map_file);
if(!$valid_input){  die "Invalid $sample_adapter_map_file line: $input_line.\n";} else{}; ##success!!!

} ##passed various checks




####MAIN PIPELINE START

foreach my $sample_prefix_key ( keys %file_adapter_prefixes_map){
  
    if($sample_prefix_key =~m/,/){ ##paired end sample
       
    my $success=processPairedEnd($sample_prefix_key, $file_adapter_prefixes_map{$sample_prefix_key}[0], $file_adapter_prefixes_map{$sample_prefix_key}[1] );
    
    delete $file_adapter_prefixes_map{$sample_prefix_key} if $success == 0; ##sample processing failed, exclude it from the list       
    
    }else{ ##single end sample
    
    my $success=processSingleEnd($sample_prefix_key, $file_adapter_prefixes_map{$sample_prefix_key}[0], $file_adapter_prefixes_map{$sample_prefix_key}[1] );
    delete $file_adapter_prefixes_map{$sample_prefix_key} if $success == 0; ##sample processing failed, exclude it from the list
 
     } 



} ##finish looping through  keys




###PART 2###merging FPKM files, only those files that were  properly processed during all previous step
print "Rounding and flooring FPKM data.\n";
mkdir "fpkm_files" if (! -e "fpkm_files"); ##create this if it doesn't exist
my $num_valid_samples=keys %file_adapter_prefixes_map;
my $string_for_merge_script=joinpath($path_to_rna_seq_tools, "merge.plx")." -a ";

for my $sample_prefix_key  ( keys %file_adapter_prefixes_map){
    
    my $base_sample_name=getBaseNameForMerge($sample_prefix_key);
    my $sample_file_for_merging=getPrepSampleForMerge($base_sample_name, $sample_prefix_key);
    $string_for_merge_script=$string_for_merge_script. $sample_file_for_merging .":1 " if $sample_file_for_merging ne 0 and $num_valid_samples > 1;
    
}
if($sample_file_for_merging ne 0 and $num_valid_samples > 1){
    
    print "Merging FPKM data\n";
    my $return_code=system("$string_for_merge_script  > fpkm_files/temp.txt");

    if($return_code != 0){
    
    errorHandler($return_code, "Merging of FPKM data failed. Exiting pipeline.");
    
}
    else {
            getFinalMergedFile("fpkm_files/temp.txt");

        print "Pipeline finished\n";}
}







###SUB ROUTINESS###
sub processSingleEnd{
    my $return_code=0; 
    my $prefix=$_[0];
    my $fastq=$_[1];
    my $adapter=$adapters_map{$_[2]};
        
    mkdir $prefix;
    my $fastq_for_tophat=$fastq;
    
    if(!$h_options{'noatrim'}){
    print "Adapter trimming $prefix (adapter $adapter)\n";

   $return_code=system("fastx_clipper -a $adapter -l $min_seq_len -i $fastq_for_tophat -o $prefix/$prefix"."_clipped.fastq -Q33" ) if !$test_mode;
       $fastq_for_tophat="$prefix/$prefix"."_clipped.fastq";

    if($return_code != 0 ){errorHandler($return_code, "Adapter trimming failed on $prefix. Sample will be skipped. "); return 0};
    }
    
    if(! $h_options{'noqtrim'}){
    print "Quality trimming $prefix\n";
   $return_code= system("fastq_quality_trimmer -t $h_options{'minqual'} -l $min_seq_len -i $fastq_for_tophat -o $prefix/$prefix"."_qual.fastq -Q33") if !$test_mode;
       $fastq_for_tophat="$prefix/$prefix"."_qual.fastq";

    if($return_code!=0 ){ errorHandler($return_code,"Quality trimming failed on $prefix. Sample will be skipped. "); return 0};
    }
    print "Running tophat on $prefix using $num_cores cores\n";
    $return_code=system(joinpath($path_to_rna_seq_tools,"tophat")." --segment-length $h_options{'seglen'} --output-dir $prefix/$prefix"."_tophat_out  -g $h_options{'hits'} -p $num_cores $h_options{'bwindex'} $fastq_for_tophat ") if !$test_mode;
    if ($return_code != 0 ){ errorHandler($return_code, "Tophat failed on $prefix. Sample will be skipped. "); return 0;}
    
   `mv $prefix/$prefix"_tophat_out/accepted_hits.bam" $prefix/$prefix"_tophat_out/"$prefix"_accepted_hits.bam"`;
       
       ###START CUFFLINKS
    print "Running cufflinks on $prefix/$prefix"."_tophat_out/$prefix"."_accepted_hits.bam\n";
    
    my $cufflinks_string=joinpath($path_to_rna_seq_tools,"cufflinks")." -I $h_options{'ilen'} -p $num_cores -o $prefix/cufflinks_output "; ##initial part
    
    if(defined($h_options{'cuffref'})){ #run cufflinks with reference
            print "with a reference file\n";
            $cufflinks_string=$cufflinks_string."-G  $h_options{'cuffref'} ";}
    if(defined($h_options{'genome'})){ ##run cufflinks
            print "use bias correction\n";
            $cufflinks_string=$cufflinks_string."-b $genome ";
        
    }
    $cufflinks_string=$cufflinks_string."$prefix/$prefix"."_tophat_out/$prefix"."_accepted_hits.bam";
    
    $return_code=system("$cufflinks_string") if !$test_mode;
    if($return_code !=0 ){errorHandler($return_code, "Cufflinks failed on $prefix. Sample will be skipped. "); return 0;}
     ####END CUFFLINKS  
       
       ##(Optional)  cuffcompare analysis
      
       if($h_options{'cuffcomp'}){
                print "Running cuffcompare on $prefix/cufflinks_output/transcripts.gtf\n";

        if(!$h_options{'cuffref'}){print "Can't run cuffcompare on $prefix because -cuffref is not set\n";}
        else{$return_code=system(joinpath($path_to_rna_seq_tools,"cuffcompare")." -r $h_options{'cuffref'} -o $prefix"."_cuffcompare_output.gtf $prefix/cufflinks_output/transcripts.gtf");}
        if($return_code != 0) {errorHandler($return_code, "Cuffcompare failed on $prefix")}
        else{
            mkdir "cuffcompare_files";
            `mv *cuffcompare_output* cuffcompare_files`;
        }
            
       }
       
       
  
        return 1 ##successfully  processed single-end sample

    
}
sub processPairedEnd{
    
    
    my $return_code=0;
    
    my $sample_prefix_key=$_[0];
    my $fastqs_string=$_[1];
    my $adapters_string=$_[2];
    
        my @prefixes=split(/,/,$sample_prefix_key);
        my @fastqs=split(/,/,$fastqs_string);
        my @adapters=split(/,/, $adapters_string);
        my @fastqs_for_tophat=@fastqs;
        
    for(my $i=0; $i<scalar @prefixes;$i++){
    mkdir $prefixes[$i];
    
    if(! $h_options{'noatrim'}){
    print "Adapter trimming $prefixes[$i] (adapter: $adapters_map{$adapters[$i]})\n";
    $return_code=system("fastx_clipper -a $adapters_map{$adapters[$i]} -l $min_seq_len -i $fastqs_for_tophat[$i] -o  $prefixes[$i]/$prefixes[$i]"."_clipped.fastq -Q33 ") if !$test_mode;
        $fastqs_for_tophat[$i]="$prefixes[$i]/$prefixes[$i]"."_clipped.fastq";

    if($return_code != 0){errorHandler($return_code, "Adapter trimming failed on $prefixes[$i]. Sample will be skipped. "); return 0;}
    
    }
   
   if (! $h_options{'noqtrim'}){
    print "Quality trimming $prefixes[$i]\n";
    $return_code= system("fastq_quality_trimmer -t $h_options{'minqual'} -l $min_seq_len  -i $fastqs_for_tophat[$i] -o $prefixes[$i]/$prefixes[$i]"."_qual.fastq -Q33") if !$test_mode;
         $fastqs_for_tophat[$i]="$prefixes[$i]/$prefixes[$i]"."_qual.fastq";

    if($return_code!=0 ){errorHandler($return_code,"Quality trimming failed on $prefixes[$i]. Sample will be skipped. ");  return 0;}

   }
} ##finished  adapter and quality trimming both fastq's for paired end sample
    print "Running tophat on $prefixes[0] and $prefixes[1] using $num_cores cores\n";
#print joinpath($path_to_rna_seq_tools,"tophat")." --segment-length $h_options{'seglen'} --output-dir $prefixes[0]/$prefixes[0]"."_tophat_out  -g $h_options{'hits'} -p $num_cores -r $h_options{'r'} $bwindex $fastqs_for_tophat[0] $fastqs_for_tophat[1] \n";
    $return_code=system(joinpath($path_to_rna_seq_tools,"tophat")." --segment-length $h_options{'seglen'} --output-dir $prefixes[0]/$prefixes[0]"."_tophat_out  -g $h_options{'hits'} -p $num_cores -r $h_options{'r'} $bwindex $fastqs_for_tophat[0] $fastqs_for_tophat[1] ") if !$test_mode;
    if ($return_code != 0 ){ errorHandler($return_code, "Tophat failed on $prefixes[0] and $prefixes[1]. Sample will be skipped. ");return 0;}

 `mv $prefixes[0]/$prefixes[0]_tophat_out/accepted_hits.bam $prefixes[0]/$prefixes[0]_tophat_out/$prefixes[0]_accepted_hits.bam`;
  
      
       ###START CUFFLINKS
    print "Running cufflinks on $prefixes[0]/$prefixes[0]"."_tophat_out/$prefixes[0]"."_accepted_hits.bam\n";
    
    my $cufflinks_string=joinpath($path_to_rna_seq_tools,"cufflinks")." -I $h_options{'ilen'} -p $num_cores -o $prefixes[0]/cufflinks_output "; ##initial part
    
    if(defined($h_options{'cuffref'})){ #run cufflinks with reference
            print "with a reference file\n";
            $cufflinks_string=$cufflinks_string."-G  $h_options{'cuffref'} ";}
    if(defined($h_options{'genome'})){ ##run cufflinks
            print "use bias correction\n";
            $cufflinks_string=$cufflinks_string."-b $genome ";
        
    }
    $cufflinks_string=$cufflinks_string."$prefixes[0]/$prefixes[0]"."_tophat_out/$prefixes[0]"."_accepted_hits.bam";
    
    $return_code=system("$cufflinks_string") if !$test_mode;
   
    if($return_code !=0 ){errorHandler($return_code, "Cufflinks failed on $prefixes[0]. Sample will be skipped. "); return 0;}
     ####END CUFFLINKS  
        
        ##(Optional)  cuffcompare analysis
      
       if($h_options{'cuffcomp'}){
        print "Running cuffcompare on $prefixes[0]/cufflinks_output/transcripts.gtf\n";
        if(!$h_options{'cuffref'}){print "Can't run cuffcompare on $prefixes[0] because -cuffref is not set\n";}
        else{$return_code=system(joinpath($path_to_rna_seq_tools,"cuffcompare")." -r $h_options{'cuffref'} -o $prefixes[0]"."_cuffcompare_output.gtf $prefixes[0]/cufflinks_output/transcripts.gtf");}
        if($return_code != 0) {errorHandler($return_code, "Cuffcompare failed on $prefixes[0]");} else{
            mkdir "cuffcompare_files";
            `mv *cuffcompare_output* cuffcompare_files`;
        }
            
       }
    
    return 1 ##successfully  processed paired-end sample
    
}



sub getFinalMergedFPKM{ ##NOT USED
    my $file=$_[0];
    my $prefix_key=$_[1];
    
    open(INPUT, "<$file") || die "Can't open $file. Exiting pipeline.\n";
    open(OUTPUT, ">fpkm_files/fpkm_all_samples_rounded.txt") || die "Can't create output file: fpkm_all_samples_rounded. Exiting pipeline. ";
    while(my $line=<INPUT>){
        chomp $line;
        my @liner=split(/\t|\s+/, $line);
        my @col_1=split(/;/, $liner[0]);
       my $num_samples=scalar @liner/3;
            my $sample1=$liner[1];
          if($sample1 < 0.3  and $sample1!~m/FPKM/){ $sample1=0.3;}
          elsif ($sample1 > 0.3 and $sample1!~m/FPKM/){ $sample1=$sample1+0.3;}
          print OUTPUT "$liner[0]\t$col_1[0]\t$col_1[1]\t$sample1\t";
          
            for (my $i=3; $i<scalar @liner;$i++){
                if($i % 3 ==1){
                    my $sample_i=$liner[$i];
                    if($sample_i < 0.3 and $sample_i!~m/FPKM/){ $sample_i=0.3;}
          elsif ($sample_i > 0.3 and $sample_i!~m/FPKM/) { $sample_i=$sample_i+0.3;}
                    
                     print OUTPUT $sample_i."\t" 

                }
    
        
    }
            print OUTPUT "\n";
    
}
}
sub prepSampleForMerging{ ##NOT USED
    my $sample_prefix_key=shift @_;
    my $prefix="";

#process cufflinks gene.fpkm file
if($sample_prefix_key=~m/,/){ #paired end sample
 my @prefixes=split(/,/, $sample_prefix_key);
 $prefix=$prefixes[0];
 
}else{ #single end sample

    $prefix=$sample_prefix_key;
    
}

   
  #  prep FPKM file for merging with other samples (merge command not included)
    open(INPUT, "<$prefix/cufflinks_output/genes.fpkm_tracking") || errorHandler("none", "Cannot open $prefix/cufflinks_output/genes.fpkm_tracking. The sample will be skipped." ) & delete $file_adapter_prefixes_map{$sample_prefix_key} & return 0;
    open(OUTPUT, ">fpkm_files/$prefix"."_fpkm_tomerge");
 print OUTPUT "tracking_id;locus\t$prefix"."_FPKM\tstatus\n";
   while(my $line=<INPUT>){
    chomp $line;
    my @liner=split(/\t+|\s+/, $line);
    print OUTPUT "$liner[0]".";$liner[6]\t$liner[9]\t$liner[12]\n";
   }
   
   return "fpkm_files/$prefix"."_fpkm_tomerge";
   
   }

    



sub validateSampleAdapterFile{
    
    my $sample_adapter_map_file=$_[0];
    my %sample_adapter_map=();
    
    open(INPUT, "<$sample_adapter_map_file") || die "Cannot open $sample_adapter_map_file\n";
    
    while(my $line=<INPUT>){
        chomp $line;
        if($line ne ""){
            my @liner=split(/\t+|\s+/, $line);
        if($h_options{'noatrim'}){ ##only 1 column file required
            
          if($line=~m/^([^\s\t]+,[^\s\t]+).*/){ ##paired end
             my @paired_end_sample_files=split(/,/,$1);
            my $output_prefix1=getOutputPrefix($paired_end_sample_files[0]);
           my $output_prefix2=getOutputPrefix($paired_end_sample_files[1]);
            $sample_adapter_map{$output_prefix1.",".$output_prefix2}[0]=$1;
            $sample_adapter_map{$output_prefix1.",".$output_prefix2}[1]="";
            
          }elsif($line=~m/^([^\s\t\,]+).*/) ##single end
          {
            
            my $output_prefix=getOutputPrefix($1);
            $sample_adapter_map{$output_prefix}[0]=$1; ##store fastq file info
            $sample_adapter_map{$output_prefix}[1]="";
            
          }
          else {
            
            die ("Badly formatted line $line\n");
          }
      
    }else{ ##Expect at least 2 columns
        
        if($line=~m/^([^\s\t]+,[^\s\t]+)[\t\s]+([\d]+,[\d]+).*/){ ##paired end   
        
           my @paired_end_sample_files=split(/,/, $1);
           
    die ("Badly formatted line $line\n") if (scalar @liner <2 or ($paired_end_sample_files[0] eq "" or $paired_end_sample_files[1] eq "" ));       
         my $output_prefix1=getOutputPrefix($paired_end_sample_files[0]);
           my $output_prefix2=getOutputPrefix($paired_end_sample_files[1]);
            $sample_adapter_map{$output_prefix1.",".$output_prefix2}[0]=$1; ##store fastq file information
            
            my @adapters=split(/,/,$2);
            die ("Adapter $adapters[0] in line: $line is invalid. Exiting.\n") if !$adapters_map{$adapters[0]};
            die ("Adapter $adapters[1] in line: $line is invalid. Exiting.\n") if !$adapters_map{$adapters[1]};
             $sample_adapter_map{$output_prefix1.",".$output_prefix2}[1]=$2; ##store adapter file informatio     

            
        }
        elsif($line=~m/^([^\s\t\,]+)[\t\s]+([\d]+).*/){ ## single end
            
            my $output_prefix=getOutputPrefix($1);
        $sample_adapter_map{$output_prefix}[0]=$1; ##store fastq file info
            die ("Adapter $2 in: $line is invalid. Exiting pipeline.\n") if !$adapters_map{$2};
        $sample_adapter_map{$output_prefix}[1]=$2; ##store adapter file info 
            
        }
        else{
            die ("Badly formatted line $line\n");
            
        }
            
        }
        
        
        } else {#empty line
            }
        
    
    }
    
    return (1,"",%sample_adapter_map);
}
sub getNumCoresForPipeline{ ##Unused
    my $total_cores=`grep -ic ^ processor /proc/cpuinfo`;
    my $cores_to_run=1;
    
    $cores_to_run = floor($total_cores/4) if $total_cores>=8;
    
    
    return $cores_to_run;
    
    }  
sub errorHandler{
my $standard_error_string=$_[0];
my $custom_error_string=$_[1];
print "$custom_error_string Error code: $standard_error_string\n";
	
}
sub getOutputPrefix{
    my $input_fastq=$_[0];
    my $output_prefix;
    
    if($input_fastq=~ m/(.*\/)*([^\.\/]+)(\.[^\.]+){0,1}/){ ##if has an extension
        
       $output_prefix=$2;
       
    }
    else {
        die ("Unrecognized file name format of $input_fastq.\n");
    }
    return $output_prefix;
    
}
    sub joinpath{
        
        my $path1=$_[0];
        my $path2=$_[1];
        
        my $last_char= substr($path1,length($path1)-1,1);
        if($last_char ne "/"){
            
            return $path1."/".$path2;
        }
        else{
            return $path1.$path2;
        }
    }
    
sub getBaseNameForMerge{
    
    my $prefix=$_[0];
    my @arr=();
    if($prefix=~m/,/) {
        @arr=split(/,/, $prefix);
        
    return $arr[0];}
    else{
        
        return $prefix;
    }
}
sub getPrepSampleForMerge{
    my $base_sample_name=$_[0];
    my $sample_prefix_key=$_[1];
    
    open(INPUT, "<$base_sample_name/cufflinks_output/genes.fpkm_tracking") || errorHandler("none", "Cannot open $base_sample_name/cufflinks_output/genes.fpkm_tracking. The sample will be skipped." ) & delete $file_adapter_prefixes_map{$sample_prefix_key} & return 0;
    open(OUTPUT1, ">fpkm_files/$base_sample_name"."_fpkm_tomerge");
    open(OUTPUT2, ">fpkm_files/$base_sample_name"."_select_fields_unfloored");
    open(OUTPUT3, ">fpkm_files/$base_sample_name"."_select_fields_floored_and_rounded");

 print OUTPUT1 "tracking_id;locus\t$base_sample_name"."_FPKM\n";
 print OUTPUT2 "tracking_id;locus\ttracking_id\tlocus\t$base_sample_name"."_FPKM\tstatus\n";
  print OUTPUT3 "tracking_id;locus\ttracking_id\tlocus\t$base_sample_name"."_FPKM\tstatus\n";

 
   while(my $line=<INPUT>){
    chomp $line;
    my @liner=split(/\t+|\s+/, $line);
    if($line!~/FPKM/){
        my $sample_FPKM=$liner[9];
        my $new_sample_FPKM=$sample_FPKM;
        $new_sample_FPKM=0 if($sample_FPKM < 0.3);
        
        $new_sample_FPKM+=0.3;
        
    print OUTPUT1 "$liner[0]".";$liner[6]\t$new_sample_FPKM\n";
    print OUTPUT2 "$liner[0]".";$liner[6]\t$liner[0]\t$liner[6]\t$sample_FPKM\t$liner[12]\n";
    print OUTPUT3 "$liner[0]".";$liner[6]\t$liner[0]\t$liner[6]\t$new_sample_FPKM\t$liner[12]\n";

    }
   }
   
   return "fpkm_files/$base_sample_name"."_fpkm_tomerge";
    
    
    
}

sub getFinalMergedFile{
    
    my $file=$_[0];
    open(INPUT, "<$file") || die ("Can't open $file for final merging\n");
    open(OUTPUT, ">fpkm_files/all_FPKM_merged") || die ("Can't open fpkm_files/all_FPKM_merged") ;
    
    while(my $line=<INPUT>){
        chomp $line;
        
        my @liner=split(/\t|\s+/, $line);
        my @id=split(/;/, $liner[0]);
        print OUTPUT "$liner[0]\t$id[0]\t$id[1]\t$liner[1]\t";
        
        for(my $i=2; $i<scalar @liner;$i++){
            
         print  OUTPUT "$liner[$i]\t"   if ($i % 2 == 1);
        }
        
        print OUTPUT "\n";
       
    }
    
    close INPUT;
    `rm $file`;
    
    
}


