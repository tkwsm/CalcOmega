#!/usr/bin/ruby

require 'tempfile'
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

peph = {}
rnah = {}
cdsh = {}
pf.each{ |e| peph[ e.definition ] = e.aaseq }
rf.each{ |e| rnah[ e.definition ] = e.naseq }

STDERR.puts "STEP 1/6 Done"

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

STDERR.puts "STEP 2/6 Done"

def merge_region( region_array )

  region_array.sort!{|x, y| x[0].to_i <=> y[0].to_i }
  merged_region_array = []
  original_size       = 0

  while region_array.size != original_size 

    break if region_array.size == 1
    region_array[0..-2].each_with_index do |r, i|
      if region_array[(i+1)][0] <= ( r[1] + 1 )
        if    ( r[1] <= region_array[(i+1)][1] )
          merged_region_array << [ r[0], region_array[(i+1)][1] ]
        elsif ( region_array[(i+1)][1] < r[1] )
          merged_region_array << [ r[0], r[1] ]
        end
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

STDERR.puts "STEP 3/6 Done"

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

STDERR.puts "STEP 4/6 Done Domain-Cluster-Hash Created."

unless File.exist?("#{protfname}.phr")
  `makeblastdb -in #{protfname} -dbtype prot`
end

domain_cluster_h.keys.each_with_index do |tid, i|
  STDERR.puts " within STEP 5/6, #{i} / #{domain_cluster_h.keys.size} done." if i % 1000 == 0
  outdat=`echo "#{tid}" | screen_list2.pl -l - -f #{mrnafname} -k | blastx -query - -db #{protfname} -outfmt 6 -max_target_seqs 1 | head -1 `
  tstart = outdat.chomp.split("\t")[6].to_i
  tend   = outdat.chomp.split("\t")[7].to_i
  pstart = outdat.chomp.split("\t")[8].to_i
  pend   = outdat.chomp.split("\t")[9].to_i
  cdsh[tid] = [ tstart, tend, pstart, pend ]
end

STDERR.puts "STEP 5/6 BlastX executed. CDS-Hash Created."

domain_cluster_h.each_key do |tid|
  ts, te, ps, pe = cdsh[ tid ]
  domain_cluster_h[ tid ].each_key do |r|
    rs = r[0]
    re = r[1]
    domains = domain_cluster_h[ tid ][r]
    if    peph[ tid ] == nil
      STDERR.puts " No Sequence for #{tid}!!"
    elsif rnah[ tid ] == nil
      STDERR.puts " No Sequence for #{tid}!!"
    else
      cs = ( ts - 1 ) + (( rs - ps )*3 ) + 1
      ce = ( ts - 1 ) + (( re - ps )*3 ) + 1 + 2
      next if cs < 0
      next if ce < 0
      STDERR.puts "Tid: #{tid}, Rs: #{rs}, Re : #{re}, Ts: #{ts}, Te: #{te}, Ps: #{ps}, Pe: #{pe}, Cs: #{cs}, Ce: #{ce}\n"
      outpf.print peph[ tid ].subseq( rs, re ).to_fasta( "#{tid}_#{rs}_#{re} #{domains.join(", ")}", 60 )
      outrf.print rnah[ tid ].subseq( cs, ce ).to_fasta( "#{tid}_#{rs}_#{re} #{domains.join(", ")}", 60 )
    end
  end
end

STDERR.puts "STEP 6/6 Done, All Steps are Completed!"

