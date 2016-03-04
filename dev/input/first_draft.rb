#!/usr/bin/ruby
#
require 'bio'
include Bio

ortholog_table_f = open( ARGV.shift )

speA_prot_f = FlatFile.new( FastaFormat, open( ARGV.shift ) )
speB_prot_f = FlatFile.new( FastaFormat, open( ARGV.shift ) )

orthologs = {}
a = []
ortholog_table_f.each do |x|
  a = x.chomp.split("\t")
  speAgid = a[0]
  speBgid = a[1]
  orthologs[ speAgid ] = speBgid
end

speA_prots = {}
speB_prots = {}
speA_prot_f.each{ |e| speA_prots[ e.definition ] = e.aaseq.to_fasta( e.definition, 60 ) }
speB_prot_f.each{ |e| speB_prots[ e.definition ] = e.aaseq.to_fasta( e.definition, 60 ) }

orthologs.each_key do |speAgid|
  speBgid = orthologs[ speAgid ]
  print speA_prots[ speAgid ]
  print speB_prots[ speBgid ]
end


