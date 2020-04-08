# scmap-prod-workflow

Workflow for generating predicted labels on a library of pre-built scmap indices. Run as a part of [control workflow](https://github.com/ebi-gene-expression-group/cell-types-prod-control-workflow) that generates predictions using a variety of tools.  

This workflow relies on [dropletutils-scripts](https://github.com/ebi-gene-expression-group/dropletutils-scripts) to initialise an SCE object from input expression data. This object is then used as an input to the classifiers. A list of prediction tables in standardised format is created as output.

### Running the workflow 
To run the workflow, you will need to have [nextflow](https://www.nextflow.io/)installed. It is recommended to run the workflow in a clean conda environment. Specify the input parameters in `nextflow.config`.

Expression data are supplied as stadnard 10X-type directories that must contain three files: `matrix.mtx`, `genes.tsv`, and `barcodes.tsv`. Specify the path to this directory in the `query_10x_dir` parameter in config file.

Then run the following commands:

```
conda install nextflow 
nextflow run main.nf -profile <profile> 
```

The `-profile` flag can be set to 'local' or 'cluster' depending on where you're running the workflow.
