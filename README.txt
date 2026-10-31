Midterm Assignment: Oxford Nanopore Direct RNA Sequencing Analysis
Vignesh Gutta (vgutta@iu.edu)
System: IU Big Red 200, Slurm account c02374
Working folder: /N/scratch/vgutta/midterm
Data: SRA/Zenodo record 15853348, rna004_hsapiens_r10k.pod5
      (10,000 human direct RNA reads, RNA004 chemistry, UHR polyA-selected)

Compute environment
  Cluster        IU Big Red 200 (HPE Cray EX), SUSE Linux
  Scheduler      Slurm (account c02374)
  GPU node       NVIDIA A100-SXM4-40GB
  GPU driver     595.71.05, CUDA 13.2
  Package mgr    conda 26.3.2
  Python         3.13.13 (conda env info636)
  Shell          bash

Software versions and language
  Dorado         2.1.2+8b8fc5d     C++ / CUDA (downloaded binary)
  pod5           0.3.x (lib_pod5 0.3.49)   Python (pip)
  samtools       1.25              C (bioconda)
  NanoPlot       1.48.0            Python (pip)
  minimap2       2.31-r1302        C (bioconda)
  blue-crab      0.6.0             Python (pip)
  Squigualiser   0.7.0             Python (pip)

  NanoPlot pulled in: pandas 3.1.0, numpy 2.5.4, scipy 1.18.1,
  plotly 6.9.0, pysam 0.24.1, biopython 1.88

  The pipeline is driven entirely from bash; the tools themselves are
  written in C/C++ (Dorado, minimap2, samtools) or Python (pod5,
  NanoPlot, blue-crab, Squigualiser). No custom code was written.

Version check commands
  dorado --version
  pod5 --version
  samtools --version | head -1
  NanoPlot --version
  minimap2 --version
  blue-crab --version
  squigualiser --version
  python --version
  nvidia-smi    # on a GPU node: driver, CUDA, GPU model


Step 0. Environment setup (login node)

  module load conda
  conda create -n info636 -c conda-forge -c bioconda sra-tools velvet oases samtools minimap2 -y
  conda activate info636
  pip install pod5 NanoPlot squigualiser blue-crab

  Dorado (downloaded binary; not in conda):
    mkdir -p ~/tools && cd ~/tools
    wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-2.1.2-linux-x64.tar.gz
    tar xf dorado-2.1.2-linux-x64.tar.gz
    export PATH="$HOME/tools/dorado-2.1.2-linux-x64/bin:$PATH"

  mkdir -p /N/scratch/vgutta/midterm && cd /N/scratch/vgutta/midterm


Step 1. Download data (login node)

  wget -O rna004_hsapiens_r10k.pod5 \
    "https://zenodo.org/records/15853348/files/rna004_hsapiens_r10k.pod5?download=1"
  pod5 inspect summary rna004_hsapiens_r10k.pod5 | head
  pod5 view rna004_hsapiens_r10k.pod5 | head -3
    # Confirms RNA004 chemistry, sample_rate 4000, UHR_polyA_RNA sample


Step 2. Basecalling (interactive GPU session)

  srun -p gpu-debug -A c02374 --gres=gpu:1 --cpus-per-task=8 \
       --mem=32G --time=01:00:00 --pty bash
  # on compute node:
  module load conda; conda activate info636
  export PATH="$HOME/tools/dorado-2.1.2-linux-x64/bin:$PATH"
  nvidia-smi    # confirm A100

  mkdir -p basecall
  dorado basecaller sup rna004_hsapiens_r10k.pod5 \
    --emit-moves --device cuda:0 > basecall/basecalls.bam

  # Model auto-selected: rna004_130bps_sup@v6.0.0
  # --emit-moves is required for the move table (needed in Step 5)
  # Runtime: ~35 seconds on NVIDIA A100-SXM4-40GB
  # Output: 10,050 reads (49 filtered, 53 chimera splits)

  exit    # release the GPU


Step 3. Quality assessment with NanoPlot (login node)

  mkdir -p nanoplot
  NanoPlot --ubam basecall/basecalls.bam -o nanoplot --tsv_stats --raw

  # Note: --ubam (not --bam) because the BAM is unaligned at this stage.
  # Outputs: nanoplot/NanoPlot-report.html, NanoStats.txt, PNG plots.


Step 4. Alignment with minimap2 (login node)

  # Reference: Ensembl GRCh38 release 112, cDNA + ncRNA combined
  mkdir -p reference && cd reference
  wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
  wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
  gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.ncrna.fa.gz
  cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.ncrna.fa \
      > Homo_sapiens.GRCh38.transcriptome.fa
  cd ..

  # minimap2 needs FASTQ input, not BAM
  samtools fastq -T "*" basecall/basecalls.bam > basecall/basecalls.fastq
  # -T "*" preserves all tags (including the move table) in the FASTQ comment line

  mkdir -p align
  minimap2 -ax map-ont -uf -k 14 -t 8 -y --MD \
    reference/Homo_sapiens.GRCh38.transcriptome.fa \
    basecall/basecalls.fastq | samtools sort -o align/aligned.bam -
  samtools index align/aligned.bam

  # Flags:
  #   -ax map-ont   Nanopore preset, SAM output
  #   -uf           forward strand only (direct RNA is strand-specific)
  #   -k 14         k-mer length (ONT recommendation for direct RNA)
  #   -y            carry input SAM tags into output (preserves move table)
  #   --MD          emit MD tag for mismatches
  #   -t 8          8 CPU threads

  # Results check
  samtools flagstat align/aligned.bam
    # Primary mapped: 9562 / 10050 = 95.14%
  samtools view -c align/aligned.bam ENST00000387347.2
    # 183 reads mapped to MT-RNR2

  # Top transcripts
  samtools view -F 260 align/aligned.bam | awk '{print $3}' \
    | sort | uniq -c | sort -rn | head -10 > align/top10_transcripts.txt


Step 5. Raw signal visualization with Squigualiser (login node)

  # POD5 → BLOW5 (Squigualiser requires BLOW5/SLOW5, not POD5)
  blue-crab p2s rna004_hsapiens_r10k.pod5 -o rna004_hsapiens_r10k.blow5

  mkdir -p squig

  # Sub-step 1: reform — convert Dorado move table to PAF
  squigualiser reform \
    --sig_move_offset 0 --kmer_length 5 \
    -b basecall/basecalls.bam \
    -c -o squig/reform.paf

  # Sub-step 2: realign — combine reform PAF with aligned BAM
  squigualiser realign \
    --bam align/aligned.bam \
    --paf squig/reform.paf \
    -o squig/realign.bam
  samtools index squig/realign.bam

  # Sub-step 3: plot — generate HTML signal plot for MT-RNR2 positions 1360-1400
  squigualiser plot \
    --file reference/Homo_sapiens.GRCh38.transcriptome.fa \
    --slow5 rna004_hsapiens_r10k.blow5 \
    --alignment squig/realign.bam \
    --region ENST00000387347.2:1360-1400 \
    --output_dir squig/plot \
    --tag_name "MT-RNR2_1360-1400" \
    --plot_limit 10 --rna

  # Outputs: 10 interactive Bokeh HTML files (one per read covering the region)


Reproducing from scratch
  1. Steps 0-1 on the login node
  2. Step 2 on an interactive GPU node (gpu or gpu-debug partition)
  3. Steps 3-5 on the login node
  Everything written to /N/scratch/vgutta/midterm; copy small outputs
  (HTMLs, text summaries) to ~/precision_medicine/midterm before the
  30-day scratch purge.
