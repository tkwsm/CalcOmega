#!/usr/bin/ruby
#
require 'rubygems'
require 'bio'
include Bio

original_table_f = open( ARGV.shift )
domain_fasta_f   = open( ARGV.shift )

seqh = {}
ff =FlatFile.new( FastaFormat, domain_fasta_f )
origina_seq_name = ""
ff.each do |e|
  dom_seq_name  = e.definition.split("\s")[0]
  base_seq_name = dom_seq_name.slice(/(^\S+)_\d+_\d+$/, 1)
  seqh[ base_seq_name ] = [] if seqh[ base_seq_name ] == nil
  seqh[ base_seq_name ] << [ dom_seq_name, e.aaseq,  ]
end

a = []
original_table_f.each do |x|
  a = x.chomp.split("\s")
  speAgi = a[0] 
  speBgi = a[1] 
  next unless seqh[ speAgi ]
  seqh[ speAgi ].each do |speAdom_info|
    speAdomgi  = speAdom_info[0]
    speAdomseq = speAdom_info[1]
    print "#{speAdomgi}\t#{speBgi}\n"
  end
end


