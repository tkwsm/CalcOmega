#!/usr/bin/bash
#
# = calcomega.sh  (v1.0) - Calculate Omega (dN/dS) for two set of Proteomes
#
# Copyright::   Copyright (C) 2016-
#               Takeshi Kawashima <kawashima38@gmail.com>
#
# Usage: calcomega.sh -o <Ortholog Table>    \
#                     -a <species-A-prot.fa> \
#                     -b <species-B-prot.fa> \
#                     -x <species-A-transcript.fa> \
#                     -y <species-B-transcript.fa> \
#                     [options]
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

usage="USAGE: ./program -o <Ortholog Table> -a <species-A-prot.fa> -b <species-B-prot.fa> -x <species-A-transcript.fa> -y <species-B-transcript.fa> -h <<show usage>> \n \n Please edit the following variable \"pal2nal\" to your installed path. \n"

# PATH to pal2nal
 pal2nal="/home/takeshik/bin/pal2nal.pl"
# pal2nal="/Users/takeshik/bin/pal2nal.pl"
# pal2nal="/home/takeshik/src/pal2nal/pal2nal.v14/pal2nal.pl"

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
echo "AGID	BGID	CLUSTALSCORE" >$orthologtable.clw
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
  newtmppep_fa=${tmppep_fa##*/}
  cp $tmppep_fa $newtmppep_fa.fa
  gset=$agid$'\t'$bgid
  gset=`echo $gset`
#  clwscore=`clustalw2 -infile=./$newtmppep_fa.fa 2>dnds.err | grep Alignment |grep Score | awk '{print $3}' `
  clwscore=`clustalw2 -infile=./$newtmppep_fa.fa 2>dnds.err | grep Sequences | grep Aligned |grep Score | awk '{print $5}' `
  echo -e "$gset\t$clwscore" >>$orthologtable.clw

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
#################################################################
#  Get dN/dS Only
#################################################################
  lwl85=`cat $tmpout_paml | grep LWL85:  | grep -v LWL85m | awk '{print $10}' 2>>dnds.err`
  lwl85m=`cat $tmpout_paml | grep LWL85m:| awk '{print $10}' 2>>dnds.err`
  lpb93=`cat $tmpout_paml | grep LPB93:  | awk '{print $10}' 2>>dnds.err`
#################################################################
#  Get more information
#################################################################
#  lwl85=`cat $tmpout_paml | grep LWL85:   | sed -e 's/.*w[ \f\n\r\t]\?\=[ \f\n\r\t]\?\([^ \f\n\r\t]\+\).*/\1/' 2>>dnds.err`
#  lwl85m=`cat $tmpout_paml | grep LWL85m: | sed -e 's/.*w[ \f\n\r\t]\?\=[ \f\n\r\t]\?\([^ \f\n\r\t]\+\).*/\1/' 2>>dnds.err`
#  lpb93=`cat $tmpout_paml | grep LPB93:   | sed -e 's/.*w[ \f\n\r\t]\?\=[ \f\n\r\t]\?\([^ \f\n\r\t]\+\).*/\1/' 2>>dnds.err`
#################################################################
#  Get more and more information
#################################################################
#  lwl85=`cat $tmpout_paml | grep LWL85:   2>>dnds.err`
#  lwl85m=`cat $tmpout_paml | grep LWL85m: 2>>dnds.err`
#  lpb93=`cat $tmpout_paml | grep LPB93:   2>>dnds.err`
#################################################################
  lwl85=`echo $lwl85`
  lwl85m=`echo $lwl85m`
  lpb93=`echo $lpb93`

  # Output the Result
  echo "$agid $bgid $lwl85 $lwl85m $lpb93 "

  # Remove Temporal Files
  if [ -e $tmppep_fa ]; then 
    rm $tmppep_fa
  fi 
  if [ -e $tmprna_fa ]; then
    rm $tmprna_fa
  fi
  if [ -e $tmpyn_ctl ]; then
    rm $tmpyn_ctl
  fi
  if [ -e $tmpout_paml ]; then
    rm $tmpout_paml
  fi
  if [ -e $tmpaln_4paml ]; then
    rm $tmpaln_4paml
  fi
  if [ -e $tmpyn00_ctl ]; then
    rm $tmpyn00_ctl
  fi
  if [ -e $tmppep_aln ]; then
    rm $tmppep_aln
  fi
  if [ -e ./$newtmppep_fa.fa ]; then 
    rm ./$newtmppep_fa.fa
  fi
  if [ -e ./$newtmppep_fa.aln ]; then
    rm ./$newtmppep_fa.aln
  fi
  if [ -e ./$newtmppep_fa.dnd ]; then
    rm ./$newtmppep_fa.dnd
  fi
  if [ -e ./2YN.dN ]; then
    rm ./2YN.dN
  fi
  if [ -e ./2YN.dS ]; then
    rm ./2YN.dS
  fi
  if [ -e ./2YN.t ]; then
    rm ./2YN.t
  fi
  if [ -e ./rst ]; then
    rm ./rst
  fi
  if [ -e ./rst1 ]; then
    rm ./rst1
  fi
  if [ -e ./rub ]; then
    rm ./rub
  fi

done

