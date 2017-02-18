#phyml to find initial tree using BioNJ
phyml --quiet -i $1.phy -q -t e -a e -o lr -b 0
mv $1.phy_phyml_stats.txt initial_phyml_stats.txt
mv $1.phy_phyml_tree.txt initial_tree.newick
#raxml to find ml topology
raxml -f d -T 6 -j -s $1.fasta -n topology -m GTRGAMMA -t initial_tree.newick 
#phyml to optimise branch lengths on ml topology
phyml --quiet -i $1.phy -q -t e -a e -o lr -u RAxML_result.topology -b 0
cp $1.phy_phyml_tree.txt $1.ml.tree
#use raxml to produce boostrap trees
raxml -f a -s $1.fasta -#200 -m GTRGAMMA -n bootstrap -T 16 ­x 12345 ­p 98765 
#use raxml to project boostrap values on to ML tree
raxml -f b -t $1.ml.tree -z RAxML_bootstrap.bootstrap -m GTRGAMMA -n final
cp RAxML_result.final $1.ml_bootstrap.tree
