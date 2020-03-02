/mnt/scratch/boxiang/projects/open_reading_frame/src/glimmer3.02/scripts/g3-from-training.csh \
    /mnt/scratch/boxiang/projects/open_reading_frame/data/COVID19_dna.txt \
    /mnt/scratch/boxiang/projects/open_reading_frame/data/COVID19.nh \
    ../processed_data/glimmer/COVID19.from-training

/mnt/scratch/boxiang/projects/open_reading_frame/src/glimmer3.02/scripts/g3-iterated.csh \
    /mnt/scratch/boxiang/projects/open_reading_frame/data/COVID19_dna.txt \
    ../processed_data/glimmer/COVID19.iterated


/mnt/scratch/boxiang/projects/open_reading_frame/src/glimmer3.02/scripts/g3-iterated.csh \
    ../processed_data/prepare_g3_input/g3_input.fasta \
    ../processed_data/glimmer/gisaid/iterated


python3 alignment.py ../data/blast_aligned_wuhan_vs_90seqs ../processed_data/alignment/

python3 coord2ref.py

