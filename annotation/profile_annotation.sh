out_dir=../processed_data/annotation/profile_annotation/
mkdir -p $out_dir
python3 -m cProfile -o $out_dir/stat.prof annotation/annotation.py -f data/example.fasta -o results/
python3 -m gprof2dot -f pstats -o $out_dir/stat.dot $out_dir/stat.prof
dot -o $out_dir/stat.pdf -Tpdf $out_dir/stat.dot