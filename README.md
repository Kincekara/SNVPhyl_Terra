# SNVPhyl_Terra

This is a WDL translation of [SNVPhyl_Nextflow](https://github.com/DHQP/SNVPhyl_Nextflow) based on [SNVPhyl](https://github.com/phac-nml/snvphyl-galaxy/blob/development/docs/workflows/SNVPhyl/1.0.1/snvphyl-workflow-1.0.1.ga) pipeline. The workflow is optimized for [Terra](https://terra.bio/) with a few changes:

1. You can enter a full taxon name (e.g. "Escherichia coli") or a NCBI accession number (e.g. "GCF_000005845.2") instead of giving a fasta file for the reference genome. The workflow will download the reference genome automatically. The pipeline hierarchically selects one of them:
```fasta file > accession number > taxon name```. Although taxon name is a required input, it will be ignored when you give accession number or a reference fasta file. 

2. This pipeline creates an html summary report besides the standard outputs of SNVPhyl. You can find a heatmap styled SNV matrix and a phylogenetic tree inside the report. 
<br>[Example report](files/snvphyl_report.html)

The original SNVPhyl pipeline was written by Aaron Petkau. You can find more information in SNVPhyl [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5628696/) and [documentation](https://snvphyl.readthedocs.io/en/latest/). Please keep in mind that this is indirect adaptation of SNVPhyl from Jill Hagey's SNVPhyl_Nextflow pipeline.

![snvphyl](https://snvphyl.readthedocs.io/en/latest/images/snvphyl-overview.png)

## Terra
### Installation
You can install SNVPhyl_Terra from [Dockstore](https://dockstore.org/workflows/github.com/Kincekara/SNVPhyl_Terra/SNVPhyl:main?tab=info) as usual.

### Running pipeline
SNVPhyl_Terra requires a sample set to run. 

| Variable | Attribute | Required? |
| --- | --- | --- |
| read1 | this.*{sampleset name}s*.*{read1*} | required |  
| read2 | this.*{sampleset name}s*.*{read2*} | required |
| samplename | this.*{sampleset name}s*.*{id}* | required|
| taxon | *"{reference taxon}"* | required |
| accession | *{reference accession number}* | optional |
|reference | *{reference.fasta}* | optional |
|window_size | *{integer} (default: 11)* | optional |
|density_threshold | *{integer} (default: 2 )* | optional |


## Local Run
If you want to use the workflow in your local computer, you can use **wf_snvphyl_local.wdl** which is prepared for that purpose. You will need a workflow manager (miniwdl or cromwell) and a container runtime (docker, singularity, etc.) in your path.

Install workflow from github:

```git clone https://github.com/Kincekara/SNVPhyl_Terra.git```

Prepare inputs:

Create a tab separated text file (e.g. samples.tsv) including your sample ids and file paths like below:
```
sample1 /path_to_sample1_read1.fastq.gz /path_to_sample1_read2.fastq.gz
sample2 /path_to_sample2_read1.fastq.gz /path_to_sample2_read2.fastq.gz
...
```

Create a json file for inputs (e.g. inputs.json)
```
{
    "snvphyl.reference": "/path/to/reference.fasta",
    "snvphyl.inputSamplesFile": "/path/to/samples.tsv"
}
```
Run the workflow by using **wf_snvphyl_local.wdl**
```
miniwdl run ~/SNVPhlyl_Terra/workflows/wf_snvphyl_local.wdl -i inputs.json
```