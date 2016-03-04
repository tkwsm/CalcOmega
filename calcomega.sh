#!/usr/bin/bash
#
# = calcomega.sh  (v1.0) - Calculate Omega (dN/dS) for two set of Proteomes
#
# Copyright::   Copyright (C) 2016-
#               Takeshi Kawashima <kawashima38@gmail.com>
#
# Usage: calcomega.sh -o -o <Ortholog Table>    \
#                        -a <species-A-prot.fa> \
#                        -b <species-B-prot.fa> \
#                        -x <species-A-transcript.fa> \
#                        -y <species-B-transcript.fa> \
#                        [options]
#
#     <Ortholog Table>:          Tab Separated Gene IDs
#                                example. 
#                                SPU_000001	HPU_000100
#                                SPU_000002	HPU_001301
#                                SPU_000003	HPU_000023
#                                ...
#
#     <species-A-prot.fa>:       Protein Sequences (multi-fasta)
#     <species-B-prot.fa>:       Protein Sequences (multi-fasta)
#     <species-A-transcript.fa>: Protein Sequences (multi-fasta)
#     <species-B-transcript.fa>: Protein Sequences (multi-fasta)
#
#     Options:           -h <<show usage>>"
#
#     Requirement
#
#     (1) screen_list2.pl by Jarrod Chapman in JGI
#
#     (2) yn00 in PAML from the following
#         http://abacus.gene.ucl.ac.uk/software/paml.html
#
#     (3) pal2nal from the following
#         http://www.bork.embl.de/pal2nal/
#
#     (4) clustalw2
#         http://www.ebi.ac.uk/Tools/msa/clustalw2/
#

usage="USAGE: ./program -o <Ortholog Table> -a <species-A-prot.fa> -b <species-B-prot.fa> -x <species-A-transcript.fa> -y <species-B-transcript.fa> -h <<show usage>>"

# PATH to pal2nal
pal2nal="/home/takeshik/src/pal2nal/pal2nal.v14/pal2nal.pl"

# Template to create a control file for yn00 in the PAML
ctl1="      seqfile = "
ctl2="      outfile = " 
ctl3="      verbose = 0  "
ctl4="        icode = 0  "
ctl5="    weighting = 0  "
ctl6="   commonf3x4 = 0  "

# Option Arguments Definition
aprot="a"
bprot="b"
amrna="x"
bmrna="y"
orthologtable="o";
helpopt="h";
while getopts "o:a:b:x:y:h" opt; do
  case $opt in
    a ) aprot=$OPTARG;;
    b ) bprot=$OPTARG;;
    x ) amrna=$OPTARG;;
    y ) bmrna=$OPTARG;;
    o ) orthologtable=$OPTARG;;
    h ) echo $usage
        exit 1 ;;
    \? ) echo $usage 
         exit 1 ;; 
  esac
done

# Variables Definition
agid=""
bgid=""
lwl85=""
lwl85m=""
lpb93=""

# Start Main Loop
cat $orthologtable | while read line
do 

  # Create Temporal Files
  tmppep_fa=$(mktemp)
  tmprna_fa=$(mktemp)
  tmpyn_ctl=$(mktemp)
  tmpout_paml=$(mktemp)
  tmpaln_4paml=$(mktemp)
  tmpyn00_ctl=$(mktemp)
  tmppep_aln=$(mktemp)

  # Variables re-Definition
  agid=""
  bgid=""
  agid=` echo $line | awk '{print $1}' `
  bgid=` echo $line | awk '{print $2}' `

  # Create Multifasta both for Proteins and Nucleotides
  echo $agid | screen_list2.pl -l - -f $aprot -k  >$tmppep_fa
  echo $bgid | screen_list2.pl -l - -f $bprot -k >>$tmppep_fa
  echo $agid | screen_list2.pl -l - -f $amrna -k  >$tmprna_fa
  echo $bgid | screen_list2.pl -l - -f $bmrna -k >>$tmprna_fa

  # Create a Temporal File in Local then do ClustalW2 for the Proteins
  newtmppep_fa=`echo $tmppep_fa | sed 's/\/tmp\///' `
  cp $tmppep_fa $newtmppep_fa.fa
  clustalw2 -infile=./$newtmppep_fa.fa 1>dnds.err 2>dnds.err
  
  # Do PAL2NAL
  perl $pal2nal $newtmppep_fa.aln $tmprna_fa -nomismatch -nogap -output paml 1> $tmpaln_4paml 2>>dnds.err

  # Create Control Files for the yn00 in PAML
  echo $ctl1 $tmpaln_4paml  >$tmpyn00_ctl 
  echo $ctl2 $tmpout_paml  >>$tmpyn00_ctl
  echo $ctl3 >>$tmpyn00_ctl
  echo $ctl4 >>$tmpyn00_ctl
  echo $ctl5 >>$tmpyn00_ctl
  echo $ctl6 >>$tmpyn00_ctl

  # Do yn00 in PAML
  yn00 $tmpyn00_ctl 1>>dnds.err 2>>dnds.err

  # Variables re-Definition
  lwl85=""
  lwl85m=""
  lpb93=""
#  lwl85=`cat $tmpout_paml | grep LWL85:  | grep -v LWL85m | awk '{print $10}' 2>>dnds.err`
#  lwl85m=`cat $tmpout_paml | grep LWL85m:| awk '{print $10}' 2>>dnds.err`
#  lpb93=`cat $tmpout_paml | grep LPB93:  | awk '{print $10}' 2>>dnds.err`
  lwl85=`cat $tmpout_paml | grep LWL85:   | sed -e 's/.*w[ \f\n\r\t]\?\=[ \f\n\r\t]\?\([^ \f\n\r\t]\+\).*/\1/' 2>>dnds.err`
  lwl85m=`cat $tmpout_paml | grep LWL85m: | sed -e 's/.*w[ \f\n\r\t]\?\=[ \f\n\r\t]\?\([^ \f\n\r\t]\+\).*/\1/' 2>>dnds.err`
  lpb93=`cat $tmpout_paml | grep LPB93:   | sed -e 's/.*w[ \f\n\r\t]\?\=[ \f\n\r\t]\?\([^ \f\n\r\t]\+\).*/\1/' 2>>dnds.err`

  lwl85=`echo $lwl85`
  lwl85m=`echo $lwl85m`
  lpb93=`echo $lpb93`

  # Output the Result
  echo "$agid $bgid $lwl85 $lwl85m $lpb93 "

  # Remove Temporal Files
  rm $tmppep_fa
  rm $tmprna_fa
  rm $tmpyn_ctl
  rm $tmpout_paml
  rm $tmpaln_4paml
  rm $tmpyn00_ctl
  rm $tmppep_aln
  rm ./$newtmppep_fa.fa
  rm ./$newtmppep_fa.aln
  rm ./$newtmppep_fa.dnd

done

