
#2.5 未満、dep 2020 Jan 以前、残基数 50-600、非天然アミノ酸 5 %未満 のものを収集することを目指す
# 25% identity, short 80% coverage でフィルタリング

#PDB は 20210126 に RSYNC
#rsync -avz --delete ftp.pdbj.org::ftp_data/structures/divided/pdb/ ~/vshare/bioo/database/pdb/
#MMCIF しかないものについては無視
#pdb_seqres.txt.gz は同日に  https://www.rcsb.org/downloads/fasta から DL したと思う。

find pdb |grep ent.gz|xargs -I {} echo "zcat {} | grep -v -E '^(ATOM|HETATM|CONECT)'"|bash > meta_info_lines.dat 
python get_resolution_deposition.py > resotable.dat 
perl check_stat_pdb.pl
cat res*.dat > results/pdb_prop.dat

get_filtered_fasta.py 内のファイルパスを変更し、
python .\get_filtered_fasta.py > samples/filtered_aa.fas 2> err.dat

で filtered_aa.fas としてフィルタリングしたものを生成。



perl iter_blastclust.pl -in samples/filtered_aa.fas -out samples/clustered.dat -outdir samples/clusout -ident 30 -cov_long 0.5 -cov_short 0.5
