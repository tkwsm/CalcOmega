#!/usr/bin/ruby

require 'rubygems'
require 'bio'
include Bio

if ARGV.size == 0
  print "create_new_table.rb <speAprot.fa> <speBprot.fa> <original-table>\n"
  exit
end

fileA=FlatFile.new( FastaFormat, open( ARGV.shift ))
fileB=FlatFile.new( FastaFormat, open( ARGV.shift ))
originaltable=open( ARGV.shift )

oh = {}
originaltable.each_with_index do |x, i|
  a = x.chomp.split("\t")
  speAid = a[0]
  speBid = a[1]
  oh[i] = [ speAid, speBid ]
end

fah = {}
fbh = {}

def create_hash( afile )
  h = {}
  doms = ""
  gid  = ""
  afile.each do |e|
    gid   = e.definition.split("_")[0..-3].join("_")
    subid = e.definition
    doms  = e.definition.slice(/^\S+\s(.+)/, 1)
    h[ gid ] = {} if h[ gid  ] == nil
    h[ gid ][ doms ] = [] if h[ gid ][ doms ] == nil
    h[ gid ][ doms ] << subid
  end
  return h
end

fah = create_hash( fileA )
fbh = create_hash( fileB )

oh.each_key do |i|
  speAsubid = ""
  speBsubid = ""
  speAid, speBid = oh[i]
  next if fah[ speAid ] == nil
  fah[ speAid ].each_key do |doms|
    next if fah[ speAid ][ doms ].size > 1
    next if fbh[ speBid ] == nil
    next if fbh[ speBid ][ doms ] == nil
    next if fbh[ speBid ][ doms ].size > 1

    speAsubid = fah[ speAid ][ doms ][0]
    speBsubid = fbh[ speBid ][ doms ][0] 

    print "#{speAsubid.slice(/^(\S+)/, 1)}\t#{speBsubid.slice(/^(\S+)/, 1)}\n"

  end
end
