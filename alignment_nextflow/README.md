# alignment_nextflow

Nextflow port of `alignment.smk` for **parity**: same per-lane steps, same `bam_input/work/...` and `bam_input/final/...` paths, same `metrics/...` and `bam.table`, same FASTQ YAML parsing as the Snakemake preamble.

## Snakemake ↔ Nextflow

| Snakemake | Nextflow |
|-----------|----------|
| `bwa_mem` | `BWA_MEM` |
| `samtools_readgroup` | `SAMTOOLS_READGROUP` |
| `samtools_markdup` | `SAMTOOLS_MARKDUP` |
| `input_ready` | `INPUT_READY` |
| `ValidateSamFile` | `VALIDATE_SAMFILE` |
| `ready_bam` | `READY_BAM` |
| `write_bam_table` | `WRITE_BAM_TABLE` |

Wildcards `{sample}`, `{reference}`, `{run}`, `{lane}`, `{index}` map to the same directory tree under `bam_input/work/`.

## Run

From the project directory (where `fastq.yml` and `sample.list` live), on a node where `bsub` works and paths are visible to compute nodes.

### LSF (cluster jobs via `bsub`)

1. Load Nextflow on the login/submit host (`module load nextflow`, conda env, etc.).
2. Edit `alignment_nextflow/nextflow.config` → profile `lsf` → set **`process.queue`** to your queue (default placeholder is `general`).
3. Run:

```bash
cd /path/to/your/project

nextflow run /path/to/utils/alignment_nextflow/main.nf \
  -profile lsf \
  -c /path/to/utils/alignment_nextflow/nextflow.config \
  --fastq_config fastq.yml \
  --sample_list sample.list \
  --reference_fasta /path/to/GRCh38.fa \
  --reference_key GRCh38 \
  --library_key YOUR_LIB
```

Optional LSF extras (queue, host group, account) without editing the config file:

```bash
nextflow run ... -profile lsf \
  -process.queue=your_queue_name \
  -process.clusterOptions='-P your_account -R "rusage[mem=32000]"'
```

Preview the job graph (no submissions):

```bash
nextflow run ... -profile lsf -preview
```

Run in background and write a report:

```bash
nextflow run ... -profile lsf -bg -with-report run_report.html
```

Watch LSF: `bjobs` / `bpeek <jobid>`. Nextflow logs: `.nextflow.log`; per-task logs under `work/*/.command.log`.

**Requirements:** shared filesystem for project dir, `fastq.yml` paths, reference FASTA, and `work/`; `bwa`, `samtools`, `java`, and Picard available on compute nodes (modules in `nextflow.config` or `~/.bashrc` on workers).

### Slurm

```bash
nextflow run /path/to/utils/alignment_nextflow/main.nf \
  -profile slurm \
  --fastq_config fastq.yml \
  ...
```

### Local (test, no batch system)

```bash
nextflow run ... 
```

(default executor in `nextflow.config` is `local`)

## Outputs (same as Snakemake)

- `bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/1.mapped.bam` … `5.markdup.bam`, `5.stats.txt`
- `bam_input/work/{sample}/{reference}/input.bam` (+ `.bai`)
- `metrics/{reference}/{sample}/validation_data.table`
- `bam_input/final/{sample}/{sample}.{reference}.bam` (+ `.bai`)
- `bam.table`

## Notes

- FASTQ pairing logic is ported from the Python block at the top of `alignment.smk` (`lib/inputs.nf`).
- Lane merge order uses **sorted** markdup paths, matching `map_input()` in Snakemake.
- Picard jar path defaults to `$HOME/software/picard/2.20.7/picard.jar` (override with `--picard_jar`).
- Nextflow keeps its own `work/` cache; final artifacts are written to the paths above relative to the launch directory, as in Snakemake.
