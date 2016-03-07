#!/usr/bin/ruby
#

require 'bio'
include Bio

if ARGV.size == 0 
  STDERR.puts "USAGE: get_dom.rb <hmmscan.out> <Proteins.fa> <mRNAs.fa> "
  exit
end

domtbloutf = open(ARGV.shift)
protfname = ARGV.shift
mrnafname = ARGV.shift
pf = FlatFile.new(FastaFormat, open( protfname ) )
rf = FlatFile.new(FastaFormat, open( mrnafname ) )
outpf = File.new( "./#{protfname.split("/")[-1]}.dom.fa", "w+")
outrf = File.new( "./#{mrnafname.split("/")[-1]}.dom.fa", "w+")

seqh = {}
rnah = {}
pf.each{ |e| seqh[ e.definition ] = e.aaseq }
rf.each{ |e| rnah[ e.definition ] = e.naseq }

STDERR.puts "STEP 1/4 Done"

dh = {}
a = []
domtbloutf.each do |x|
  next if x =~ /^\#/
  a = x.chomp.split("\s")
  pfamid  = a[0]
  pfamacc = a[1]
  tid     = a[3]
  dfrom   = a[19].to_i
  dto     = a[20].to_i
  dh[ tid ] = {} if dh[ tid ] == nil
  dh[ tid ][ pfamacc ] = [] if dh[ tid ][ pfamacc ] == nil
  dh[ tid ][ pfamacc ] << [ dfrom, dto ]
end

STDERR.puts "STEP 2/4 Done"

def merge_region( region_array )


  region_array.sort!{|x, y| x[0].to_i <=> y[0].to_i }
  merged_region_array = []
  original_size       = 0

  while region_array.size != original_size 

    break if region_array.size == 1
    region_array[0..-2].each_with_index do |r, i|
      if region_array[(i+1)][0] <= ( r[1] + 1 )
        merged_region_array << [ r[0], region_array[(i+1)][1] ]
        merged_region_array.concat(region_array[(i+2)..-1]) if region_array[(i+2)]
        break
      else
        if i + 2 == region_array.size
          merged_region_array << r
          merged_region_array << region_array[-1]
        else
          merged_region_array << r
        end
      end
    end
    merged_region_array.uniq!
    merged_region_array.sort!{|x, y| x[0].to_i <=> y[0].to_i }
    region_array.sort!{|x, y| x[0].to_i <=> y[0].to_i }
    original_size        = region_array.size
    region_array        = merged_region_array
    merged_region_array = []
    
  end
  return region_array 

end

# a = [[35, 112], [59, 137], [116, 183], [155, 231]]
# a  = [[1, 10], [20, 40], [30, 50], [60, 70], [71, 80], [90, 100] ]
# ma = [[1, 10], [20, 50], [60, 80], [90, 100] ]
#  p a
#  p merge_region( a )
# "expected"
# p ma

STDERR.puts "STEP 3/4 Done"

def create_domain_cluster_hash( region_array, dh, tid )
  domain_cluster_hash = {}
  region_array.each do |dc|
    domain_cluster_hash[ dc ] = []
    dcfrom = dc[0]
    dcto   = dc[1]
    dh[tid].each_key do |pfamacc|
      dh[ tid ][ pfamacc ].each do |dfrom, dto|
        if     dcto < dfrom
        elsif  dfrom <= dcto 
          domain_cluster_hash[ dc ] << pfamacc
          break
        elsif dcfrom <= dfrom and  dfrom <= dcto 
          domain_cluster_hash[ dc ] << pfamacc
          break
        elsif dfrom <= dcfrom and dcfrom <=  dto
          domain_cluster_hash[ dc ] << pfamacc
          break
        elsif dto   <  dcfrom 
        else
          STDERR.puts "Something Strange!"
        end
        domain_cluster_hash[ dc ].uniq!
      end
    end
  end
  return domain_cluster_hash
end

domain_cluster_h = {}
region_array = []
dh.each_key do |tid|
  domain_cluster_h[ tid ] = {}
  region_array = []
  dh[ tid ].each_key do |pfamacc|
    dh[ tid ][ pfamacc ].each do |dfrom, dto|
      region_array << [ dfrom, dto ]
    end
  end
  region_array.sort!{|x, y| x[0].to_i <=> y[0].to_i}
  region_array = merge_region( region_array )
  domain_cluster_h[ tid ] = create_domain_cluster_hash( region_array, dh, tid )
end

domain_cluster_h.each_key do |tid|
  domain_cluster_h[ tid ].each_key do |r|
    domains = domain_cluster_h[ tid ][r]
    if    seqh[ tid ] == nil
      STDERR.puts " No Sequence for #{tid}!!"
    elsif rnah[ tid ] == nil
      STDERR.puts " No Sequence for #{tid}!!"
    else
      outpf.print seqh[ tid ].subseq( r[0], r[1] ).to_fasta( "#{tid}_#{r[0]}_#{r[1]} #{domains.join(", ")}", 60 )
      outrf.print rnah[ tid ].to_fasta( "#{tid}_#{r[0]}_#{r[1]} #{domains.join(", ")}", 60 )
    end
  end
end

STDERR.puts "STEP 4/4 Done, All Steps are Completed!"

