# CalcOmega

# How to use the calcomega
# For calculate gene level omega, please try step (1) Only.

(1)
bash calcomega.sh -o <Ortholog Table>        \
                  -a <spe-A-prot.fa>         \
                  -b <spe-B-prot.fa>         \
                  -x <spe-A-transcript.fa>   \
                  -y <spe-B-transcript.fa>   > <omega.out>

# For calculate domain level omega, please try step (1), then the following steps (2) - (6) .

(2)
hmmscan --domeblout <spe-A.domtblout> <path2Pfam><spe-A-prot.fa> 
hmmscan --domeblout <spe-B.domtblout> <path2Pfam><spe-B-prot.fa> 

(3)
ruby get_dom.rb <spe-A.domtblout> <spe-A-prot.fa> <spe-A-transcript.fa> 
# OUTPUT: spe-A-prot.fa.dom.fa 
# OUTPUT: spe-A-transcript.fa.dom.fa 

(4)
ruby get_dom.rb <spe-B.domtblout> <spe-B-prot.fa> <spe-B-transcript.fa> 
# OUTPUT: spe-B-prot.fa.dom.fa 
# OUTPUT: spe-B-transcript.fa.dom.fa 

(5)
ruby create_new_table.rb <spe-A-prot.fa> <spe-B-prot.fa>  \
                         <Ortholog table> > <new.table>

(6)
bash calcomega.sh -o <new.table>        \
                  -a <spe-A-prot.fa.dom.fa>         \
                  -b <spe-B-prot.fa.dom.fa>         \
                  -x <spe-A-transcript.fa.dom.fa>   \
                  -y <spe-B-transcript.fa.dom.fa>   > <dom.omega.out>
