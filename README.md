exom sequencing workflow for use with qproject.

TODO: Adapt readme to exomseq.

Add a config file "params.json" in `etc`:

```json
{
    "gtf": "path/to/gtf",
    "indexed_genome": "path/to/genome/basename",
    "stranded": "yes",
    "overlap_mode": "union",
    "feature_type": "exon",
    "gff_attribute": "gene_id",
    "normalize_counts": "deseq2"
}
```

where `indexed_genome` and `gtf` are paths relative to `ref`.

`indexed_genome` is the basename of a bowtie2 index.

The parameters `stranded`, `overlap_mode`, `feature_type` and `gff_attribute`
are explained in the htseq documentation.
