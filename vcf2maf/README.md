### vcf2maf: a VCF to MAF conversion utility

[VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) is a standard text file format for storing variation data. It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome.

[Mutation Annotation Format](https://software.broadinstitute.org/software/igv/MutationAnnotationFormat) (MAF) is a tab-delimited text file with aggregated mutation information from VCF Files.

vcf2maf is a tool for converting files in Variant Call Format (VCF) to MAF format.

#### Usage

The vcf2maf converter takes vcf files as input and produces the corresponding maf files:

```
python vcf2maf.py [--help] --input-data <path/to/vcf/files> --output-directory <path/to/output/data/directory>
```

#### Options
```
-i | --input-data: comma seperated paths to .vcf files or data directories with .vcf files [Required].
-o | output-directory: path to output data directory [Optional]. If not specified, the output will be saved to vcf2maf_output folder in the current working directory.
```

### Annotating the variants

#### Merging the output MAF's

The vcf2maf tool produces one MAF per VCF file. The per-VCF MAF files can be merged to generate a single MAF with one of the existing tools provided: [merge_mafs.py](https://github.com/genome-nexus/annotation-tools/blob/master/merge_mafs.py)

The MAF merging script supports two different input arguments for merging MAFs: (1) a comma-delimited list of MAF filenames and (2) an input directory containing all MAFs to be merged. 

#### MAF Annotation

The [Genome Nexus Annotation](https://github.com/genome-nexus/genome-nexus-annotation-pipeline) tools allow for the annotation of genomic variants from a MAF using [Genome Nexus](https://www.genomenexus.org/). Please follow the pre-build and usage instructions [here](https://github.com/genome-nexus/genome-nexus-annotation-pipeline#maf-annotation)