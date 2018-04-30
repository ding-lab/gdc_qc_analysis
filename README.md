Conda environment setup:

    conda create -n gdc_qc python=3.6 notebook matplotlib sqlalchmey pandas snakemake

And a different env for crossmap (since it uses Python 2.7):

    conda create -y -n crossmap crossmap

And update `config.yaml` to point to the correct path of all input data.

Now the pipeline is managed by Snakemake:

    snakemake -l            # List all the possible rules
    snakemake make_db       # Generate the variant database
    snakemake all           # Generate all output
