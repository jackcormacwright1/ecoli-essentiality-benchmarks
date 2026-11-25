# Source tables

These three files are redistributed here so that the analysis is reproducible without a live internet
connection, and so that CI can run the tests that parse them. All three are from open-access sources.
Each is unmodified from the publisher's copy; parsing happens in `src/eeb/essentiality.py`.

| File | Source | Licence |
|---|---|---|
| `keio_baba2006_table3.xls` | Baba et al. (2006), *Construction of Escherichia coli K-12 in-frame, single-gene knockout mutants: the Keio collection*, Molecular Systems Biology 2:2006.0008. Supplementary Table 3. | open access (EMBO Press / MSB) |
| `tradis_goodall2018_tableS4.xlsx` | Goodall et al. (2018), *The Essential Genome of Escherichia coli K-12*, mBio 9:e02096-17. Supplementary Table S4. | open access (ASM / mBio) |
| `pec_genes.txt.rtf` | [PEC: Profiling of the E. coli Chromosome](https://shigen.nig.ac.jp/ecoli/pec/), gene table export. | database export, factual content |

If you would rather fetch them yourself, `scripts/01_fetch.py` prints the exact source for each one.

The SBML models are **not** committed: they are about 50MB and freely available from
[BiGG](http://bigg.ucsd.edu/). `scripts/01_fetch.py --download` retrieves the four BiGG models;
`iML1515_model_adjusted_all_corrections.xml` is a supplementary file of Bernstein et al. (2023) and the
script prints where to get it.
