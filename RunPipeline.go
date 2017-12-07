package main

import (
	"bufio"
	"fmt"
	"os"
	"os/exec"
	"strconv"
	"strings"
	"runtime"
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

func IndexReference(reference string) {
	//index genome
	bwaIndex := CreateCommand(cwd + "/bin/bwa index " + reference)
	RunCommand(bwaIndex)
}

func CheckBin() {
	cmd := "ls " + cwd + "/bin"
	lsBin := CreateCommand(cmd)
	binContents := WriteOutputToString(lsBin)
	if !strings.Contains(binContents, "bwa") || !strings.Contains(binContents, "samtools") || !strings.Contains(binContents, "bcftools") {
		fmt.Println("ERROR: Repository did not download correctly. Please redownload from https://github.com/bweestyle/PathogenID")
		os.Exit(1)
	}
}

func MakeBinExecutable() {
	executables := []string{"bwa", "samtools", "bcftools"}
	for i := 0; i < len(executables); i++ {
		fmt.Println("chmod 755 " + cwd + "/" + executables[i])
		cmd := CreateCommand("chmod 755 " + cwd + "/bin/" + executables[i])
		cmd.Run()
	}
}

func PrepareBin() {
	CheckBin()
	MakeBinExecutable()
}

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
	sam = append(sam,samFile)
	return samFile
}

func MakeBamFile(samFile string) string {
	samtoolsView := CreateCommand(cwd + "/bin/samtools view -@ " + strconv.Itoa(numProcs) + " -bS " + samFile)
	bamFile := strings.Split(samFile, ".")[0] + ".bam"
	OutputCommandToFile(samtoolsView, bamFile)
	return bamFile
}

func SortBamFile(bamFile string) string {
	samtoolsSort := CreateCommand(cwd + "/bin/samtools sort -@ " + strconv.Itoa(numProcs) + " " + bamFile)
	sortedBam := strings.Split(bamFile, ".")[0] + ".sorted.bam"
	OutputCommandToFile(samtoolsSort, sortedBam)
	return sortedBam
}

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

func MakeVCF(reference, sortedBam string) string {
	samtoolsMpileup := CreateCommand(cwd + "/bin/samtools mpileup -v --reference " + reference + " " + sortedBam)
	vcfFile := strings.Split(sortedBam, ".")[0] + ".vcf"
	OutputCommandToFile(samtoolsMpileup, vcfFile)
	vcf = append(vcf,vcfFile)
	return vcfFile
}

func CallVCF(vcfFile string) string {
	bcfToolsCall := CreateCommand(cwd + "/bin/bcftools call --threads " + strconv.Itoa(numProcs) + " --ploidy 1 -c " + vcfFile)
	calledVCFFile := strings.Split(vcfFile, ".")[0] + ".called.vcf"
	OutputCommandToFile(bcfToolsCall, calledVCFFile)
	return calledVCFFile
}

func CallVariants(reference, sortedBam string) string {
	vcfFile := MakeVCF(reference, sortedBam)
	calledVCFFile := CallVCF(vcfFile)
	return calledVCFFile
}

func MakeFolder(folder string) string {
	analysisFolder := cwd + "/"+folder+"/"
	RunCommand(CreateCommand("mkdir " + analysisFolder))
	return analysisFolder
}

func DownloadDriver() {
  download := CreateCommand("go get -u github.com/murlokswarm/mac")
  RunCommand(download)
}

func GetSampleData() []string {
	path := cwd + "/SampleData/"
	var readFiles []string
	forward := UnzipFile(path + "test_data.fastq.gz")
	reverse := UnzipFile(path + "test_rev.fastq.gz")
	readFiles = append(readFiles, forward)
	readFiles = append(readFiles, reverse)
	return readFiles
}

func pwd() string {
	pwd := CreateCommand("pwd")
	s := WriteOutputToString(pwd)
	return s
}

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


func OnlyVCF(bamFiles []string) string{
  reference,_ := HandleReference()
	IndexReference(reference)

  calledVCFFile := CallVariants(reference, bamFiles[0])
	PredictResistance(reference,calledVCFFile)
	return calledVCFFile
}

func PredictResistance(reference, calledVCFFile string) string{
	ReadReference(reference)
	ShowMutation(calledVCFFile)
	DrugRecom()
	resistanceFile := WriteResult(allDrug)
	return resistanceFile
}

func RunPipeline(readFiles []string) {
	if len(readFiles) ==0{
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
