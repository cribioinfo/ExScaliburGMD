# ExScaliburGMD: Scalable and Accurate Detection of Germline Mutations from Whole Exome Sequencing in the Cloud #


ExScaliburGMD is a germline mutation pipeline for Whole Exome Sequencing (WES) data developed by the [Center for Research Informatics (CRI) bioinformatics core](http://cri.uchicago.edu/?page_id=1185) at the University of Chicago. 

The pipeline is designed to be

* **Scalable** *Manages analysis of tens to hundreds of samples*
* **Accurate** *Over 99% sensitivity and accuracy in variant detection*
* **Automate** *Full analysis pipeline from QC to annotated germline variants with submission of one master script*
* **Comprehensive** *Implements multiple aligners and callers allowing for comparison and integration of different variant sets*
* **Flexible** *Fine control of analysis modules, parameters, and environmental settings*
* **Robust** *Extensive checkpoints and error detection features*
* **Reproducible** *Quick restart of analysis and sharing of results*
* **Real-time** *Easy access and monitor of the pipeline progress*
* **Transferable** *Runs on a desktop, work station, HPC and cloud with one single switch*

The pipeline is implemented in 

* [BigDataScript](http://pcingola.github.io/BigDataScript/) *Scripting language for data pipelines*
* [Perl](https://www.perl.org/) *Utilities*

The GMD pipeline is part of the ExScalibur suite (CRI, University of Chicago): https://exscalibur.cri.uchicago.edu.

# Quick Demo#

To gain the first look into the pipelines, we provide sample data and scripts to build a small project with 3 samples.

In directory example/

* **NA12878trio.metadata.txt** *Sample metadata table*
* **NA12878trio.pipeline.yaml** *Pipeline configuration file*

The project directory **myProject** has been prepared for you. To regenerate the project directory and all config files in the directory, try to run 
 
```
#!bash
./Build_ExScaliburGMD.NA12878trio.sh

```

Then the project directory will be regenerated, as well as the master job submission script **Submit_ExScaliburGMD.NA12878trio.sh**, which can be used to initiate the pipeline run. The full project directory structure can be viewed in **myProject/NA12878trio.tree.txt**.

This is a quick demo of the initialization of a project. To run real tests, please follow the instructions on ExScalibur website (https://exscalibur.cri.uchicago.edu).

# Documentation #

Please see the [Wiki](https://bitbucket.org/cribioinformatics/exscaliburgmd/wiki) for full documentation.

# Communication #

INSERT CONTACT INFORMATION HERE

# Release #

* **Version 0.5.0** *2015-03-31*

# Development #

* Add checkpoints
* Add VCMM variant caller
* Add hg38 genome

# More #

* [ExScaliburSMD](https://bitbucket.org/cribioinformatics/exscalibursmd) *somatic mutation detection*
* [ExScaliburViz](https://bitbucket.org/cribioinformatics/exscaliburviz) *visualization of data analysis report*