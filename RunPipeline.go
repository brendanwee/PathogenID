package main

import (
	"bufio"
	"fmt"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"
)

func CheckError(err error) {
	if err != nil {
		fmt.Println("ERROR:")
		fmt.Println(err)
		os.Exit(1)
	}
}

//Takes in a sam filename and a genome length
//returns true if @SN is found and the value of LN is equal to LN on the first line of the sam file
func CheckSamFile(samFile string, LN int) {
	file, err := os.Open(samFile)
	CheckError(err)
	scanner := bufio.NewScanner(file)
	scanner.Scan()
	words := strings.Split(scanner.Text(), ":")
	if len(words) < 3 { //check that there are at least 3 items in the first line
		fmt.Println("Less than three items in the first line. BWA mem Failed")
		os.Exit(1)
	}
	genomeLength, err := strconv.Atoi(words[2]) //the third word should be LN:num, get num
	if genomeLength != LN {
		fmt.Println("Genome length:", genomeLength, "!= ", LN, "BWA mem Failed")
		os.Exit(1)
	}
}

//indexes a reference genome
func IndexReference(reference string) {
	//index genome
	bwaIndex := CreateCommand(cwd + "/bin/bwa index " + reference)
	RunCommand(bwaIndex)
}

//checks the bin to ensure the correct binaries are stored there.
func CheckBin() {
	cmd := "ls " + cwd + "/bin"
	lsBin := CreateCommand(cmd)
	binContents := WriteOutputToString(lsBin)
	if !strings.Contains(binContents, "bwa") || !strings.Contains(binContents, "samtools") || !strings.Contains(binContents, "bcftools") {
		fmt.Println("ERROR: Repository did not download correctly. Please redownload from https://github.com/bweestyle/PathogenID")
		os.Exit(1)
	}
}

//uses chmod to change the permissions of the necessary executives in this program
func MakeBinExecutable() {
	executables := []string{"bwa", "samtools", "bcftools"}
	for i := 0; i < len(executables); i++ {
		fmt.Println("chmod 755 " + cwd + "/" + executables[i])
		cmd := CreateCommand("chmod 755 " + cwd + "/bin/" + executables[i])
		cmd.Run()
	}
}

//Checks to be sure the right files are there.Makes the bin executable
func PrepareBin() {
	CheckBin()
	MakeBinExecutable()
}

//given a reference, the length of the reference and a set of reads,
//MakeSamFile aligns the reads to the reference genome and returns a path to
//the samFile
func MakeSamFile(reference string, LN int, readFiles ...string) string {
	folders := strings.Split(strings.TrimSuffix(readFiles[0], ".fastq"), "/")
	filename := folders[len(folders)-1]
	samFile := outputPath + filename + ".sam"
	var bwaMem *exec.Cmd
	if len(readFiles) == 2 {
		bwaMem = CreateCommand(cwd + "/bin/bwa mem -t " + strconv.Itoa(numProcs) + " " + reference + " " + readFiles[0] + " " + readFiles[1])
	} else {
		bwaMem = CreateCommand(cwd + "/bin/bwa mem -t " + strconv.Itoa(numProcs) + " " + reference + " " + readFiles[0])
	}
	OutputCommandToFile(bwaMem, samFile)
	CheckSamFile(samFile, LN)
	sam = append(sam, samFile)
	return samFile
}

// calls samtools view to turn a sm file into a more compressed bam file
func MakeBamFile(samFile string) string {
	samtoolsView := CreateCommand(cwd + "/bin/samtools view -@ " + strconv.Itoa(numProcs) + " -bS " + samFile)
	bamFile := strings.Split(samFile, ".")[0] + ".bam"
	OutputCommandToFile(samtoolsView, bamFile)
	return bamFile
}

/*
calls samtools sort to sort the binary alignment file
*/
func SortBamFile(bamFile string) string {
	samtoolsSort := CreateCommand(cwd + "/bin/samtools sort -@ " + strconv.Itoa(numProcs) + " " + bamFile)
	sortedBam := strings.Split(bamFile, ".")[0] + ".sorted.bam"
	OutputCommandToFile(samtoolsSort, sortedBam)
	return sortedBam
}


/*
Aligns a single or a set of fastq files to a reference gneome.
*/
func AlignReads(reference string, readFiles []string, LN int, pairedEnd bool) string {
	var samFile string
	if pairedEnd {
		samFile = MakeSamFile(reference, LN, readFiles[0], readFiles[1])
	} else {
		samFile = MakeSamFile(reference, LN, readFiles[0])
	}

	bamFile := MakeBamFile(samFile)
	sortedBam := SortBamFile(bamFile)
	return sortedBam
}

//Takes a sorted bam file and returns a VCF file of it to the reference.
func MakeVCF(reference, sortedBam string) string {
	samtoolsMpileup := CreateCommand(cwd + "/bin/samtools mpileup -v --reference " + reference + " " + sortedBam)
	vcfFile := strings.Split(sortedBam, ".")[0] + ".vcf"
	OutputCommandToFile(samtoolsMpileup, vcfFile)
	vcf = append(vcf, vcfFile)
	return vcfFile
}

//Tkaes a VCF files and returns a called VCF file by processing it with bcfTool
func CallVCF(vcfFile string) string {
	bcfToolsCall := CreateCommand(cwd + "/bin/bcftools call --threads " + strconv.Itoa(numProcs) + " --ploidy 1 -c " + vcfFile)
	calledVCFFile := strings.Split(vcfFile, ".")[0] + ".called.vcf"
	OutputCommandToFile(bcfToolsCall, calledVCFFile)
	return calledVCFFile
}

//Given a reference file and a bam file, it call the variants of the data
// by calling samtools and bcftool
func CallVariants(reference, sortedBam string) string {
	vcfFile := MakeVCF(reference, sortedBam)
	calledVCFFile := CallVCF(vcfFile)
	return calledVCFFile
}

//Makes a folder
func MakeFolder(folder string) string {
	analysisFolder := cwd + "/" + folder + "/"
	RunCommand(CreateCommand("mkdir " + analysisFolder))
	return analysisFolder
}

/*
in older versions of the code, it would automaticall select the sample data
for you if you did not import files before clicking.
*/
func GetSampleData() []string {
	path := cwd + "/SampleData/"
	var readFiles []string
	forward := UnzipFile(path + "test_data.fastq.gz")
	reverse := UnzipFile(path + "test_rev.fastq.gz")
	readFiles = append(readFiles, forward)
	readFiles = append(readFiles, reverse)
	return readFiles
}

//returns the current working directory
func pwd() string {
	pwd := CreateCommand("pwd")
	s := WriteOutputToString(pwd)
	return s
}

/*
Only Align takes in fasta or fastq files and aligns the reads to the MTB reference
genome
*/
func OnlyAlign(readFiles []string) string {
	//identify oraganism
	reference, LN := HandleReference()
	IndexReference(reference)
	numProcs = runtime.NumCPU() //weird bug where global variable resets to 0

	var samFile string
	if len(readFiles) == 2 {
		samFile = MakeSamFile(reference, LN, readFiles[0], readFiles[1])
	} else {
		samFile = MakeSamFile(reference, LN, readFiles[0])
	}
	return samFile
}
/*
Given a bamfile, ONLYVCF calls the variants and predicts resistances. then
it returns the path to the VCFFile
*/
func OnlyVCF(bamFiles []string) string {
	reference, _ := HandleReference()
	IndexReference(reference)

	calledVCFFile := CallVariants(reference, bamFiles[0])
	PredictResistance(reference, calledVCFFile)
	return calledVCFFile
}

/*
Takes a reference filename and a VCF file name and predicts the drug resistances
based on the amino acid changes.
*/

func PredictResistance(reference, calledVCFFile string) string {
	ReadReference(reference)
	ShowMutation(calledVCFFile)
	DrugRecom()
	resistanceFile := WriteResult(allDrug)
	return resistanceFile
}

/*
RunPipeline takes two fasta or fastq files and aligns the reads to the MTB
reference genome. Then it calls its variants
*/
func RunPipeline(readFiles []string) {
	if len(readFiles) == 0 {
		fmt.Println("No files recieved")
		return
	}
	numProcs = runtime.NumCPU() //weird bug where global variable resets to 0

	pairedEnd := (len(readFiles) == 2)
	//get Filename reads()

	reference, LN := HandleReference()
	IndexReference(reference)
	sortedBam := AlignReads(reference, readFiles, LN, pairedEnd)

	calledVCFFile := CallVariants(reference, sortedBam)

	PredictResistance(reference, calledVCFFile)

}
