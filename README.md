# PathogenID

PathogenID is a bioinformatics pipeline that takes in raw sequencing data from Illumina machines and aligns the data to the Mycobacterium Tuberculosis reference genome. The alignment is the analyzed for variants and the variants are fed out to the GUI. Various graphs to visualize the data are generated and saved into subdirectories as well.

## Getting Started
in the PathogenID directory

```
go build
./PathogenID
```

Select
```
Next
```
then
```
Run Pipeline
```
in the PathogenID directory the is a subdiretory called SampleData. Here you will find 2 fastq files you can run your pipeline on.

Resulting analysis files can be found in .../PathogenID/Analyses

after running an analysis, you can visualize your data by clicking on
```
Next
```
The analysis menu will generate plots in .../PathogenID/resources/
in subdirectories: Fastq_Plots, SamPlots, VCFPlots

### Prerequisites

In the go/src folder download
go get gonum.org/v1/plot/...
go get -u github.com/murlokswarm/mac

hg mercurial must be installed as well. You can determine this by typing
```
which hg
```
if hg turns up, and the two packages are in the go src along side this package and you're good to go!

## Authors

* **Brendan Wee** - [Bweestyle](https://github.com/Bweestyle)

* **Jinke Liu**

* **Xiaodi Pan**

* **Chaoying Wang**

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
