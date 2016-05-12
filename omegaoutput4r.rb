#!/usr/bin/ruby

##################################################
if ARGV.size == 0
  print "omegaoutput4r.rb <output of calcomega>\n"
  exit
end
##################################################

h = {}
speAgid = ""
speBgid = ""
lwl85  = ""
lwl85m = ""
lpb93  = ""
a = []
ARGF.each do |x|
  a = x.chomp.split("\s")
  next if x =~ /nan/
  next if a.size != 5
  speAgid = a[0]
  speBgid = a[1]
  lwl85   = a[2]
  lwl85m  = a[3]
  lpb93   = a[4]
  akey = [speAgid, speBgid]
  aval = [ lwl85, lwl85m, lpb93 ]
  h[ akey ] = aval
end

sorted_keys = h.keys.sort{|x, y| h[x][0] <=> h[y][0] }

sorted_keys.each_with_index do |akey, i|
  h[ akey ].each_with_index do |dnds, j|
    if    j == 0
      print "#{(i+1)}\t#{akey.join("\t")}\t#{dnds}\tlwl85\n" 
    elsif j == 1
      print "#{(i+1)}\t#{akey.join("\t")}\t#{dnds}\tlwl85m\n" 
    elsif j == 2
      print "#{(i+1)}\t#{akey.join("\t")}\t#{dnds}\tlpb93\n" 
    end
  end
end
