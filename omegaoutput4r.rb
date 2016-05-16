#!/usr/bin/ruby

##################################################
if ARGV.size == 0
  print "omegaoutput4r.rb <output of calcomega> <clustalw score table>\n"
  exit
end
##################################################

omegatable = open( ARGV.shift )
clwscoref  = open( ARGV.shift )

omegah = {}
clwsch = {}

speAgid = ""
speBgid = ""
lwl85  = ""
lwl85m = ""
lpb93  = ""
a = []
omegatable.each do |x|
  next if x.downcase =~ /nan/
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
  omegah[ akey ] = aval
end

clwscoref.each do |x|
  a = x.chomp.split("\s")
  next if a.size != 3
  speAgid = a[0]
  speBgid = a[1]
  clwscor = a[2].to_f
  akey = [speAgid, speBgid]
  clwsch[ akey ] = clwscor 
end

sorted_keys = omegah.keys.sort{|x, y| omegah[x][0] <=> omegah[y][0] }

clwscore = 0.0
sorted_keys.each_with_index do |akey, i|
  clwscore = 0.0
  clwscore = clwsch[ akey ] if clwsch[ akey ] 
  omegah[ akey ].each_with_index do |dnds, j|
    if    j == 0
      print "#{(i+1)}\t#{akey.join("\t")}\t#{dnds}\tlwl85\t#{clwscore}\n" 
    elsif j == 1
      print "#{(i+1)}\t#{akey.join("\t")}\t#{dnds}\tlwl85m\t#{clwscore}\n" 
    elsif j == 2
      print "#{(i+1)}\t#{akey.join("\t")}\t#{dnds}\tlpb93\t#{clwscore}\n" 
    end
  end
end


