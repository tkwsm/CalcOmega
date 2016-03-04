#!/usr/bin/ruby
#

a2b_f = open( ARGV.shift )
b2a_f = open( ARGV.shift )

threshold_evalue = 1.0e-100


def x2y_hash( x2y_blastf, threshold_evalue )
  x2y_h = {}
  a = []
  x2y_blastf.each do |x|
    a = x.chomp.split("\t")
    xgid = a[0]
    ygid = a[1]
    eval = a[10].to_f
    next if eval > threshold_evalue
    x2y_h[ xgid ] = ygid
  end
  return x2y_h
end

a2b_h = {}
b2a_h = {}

a2b_h = x2y_hash( a2b_f, threshold_evalue )
b2a_h = x2y_hash( b2a_f, threshold_evalue )

a2b_h.each_key do |agid|
  bgid = a2b_h[ agid ]
  next unless b2a_h[ bgid ]
  zigzagagid = b2a_h[ bgid ]
  next unless zigzagagid == agid
  print "#{agid}\t#{bgid}\t#{agid}\n"
end
