package main

import(
  "strings"
  "fmt"
  "os"
  "os/exec"
  "bufio"
  "strconv"
  "runtime"
)

func CheckError(err error) {
  if err!=nil{
    fmt.Println(err)
    os.Exit(1)
  }
}

//Takes in a sam filename and a genome length
//returns true if @SN is found and the value of LN is equal to LN on the first line of the sam file
func CheckSamFile(samFile string , LN int) {
  file,err := os.Open(samFile)
  CheckError(err)
  scanner := bufio.NewScanner(file)
  scanner.Scan()
  words := strings.Split(scanner.Text(),":")
  if len(words)<3{ //check that there are at least 3 items in the first line
    fmt.Println("Less than three items in the first line. BWA mem Failed")
    os.Exit(1)
  }
  genomeLength, err := strconv.Atoi(words[2]) //the third word should be LN:num, get num
  if genomeLength!=LN{
    fmt.Println("Genome length:", genomeLength, "!= ", LN, "BWA mem Failed")
    os.Exit(1)
  }
}

func IndexReference(cwd,reference string) {
  //index genome
  bwaIndex := CreateCommand(cwd+"/bin/bwa index "+reference)
  RunCommand(bwaIndex)
}

func CheckBin(cwd string){
  cmd := "ls "+cwd+"/bin"
  lsBin:=CreateCommand(cmd)
  binContents := WriteOutputToString(lsBin)
  if !strings.Contains(binContents, "bwa") || !strings.Contains(binContents, "samtools") || !strings.Contains(binContents, "bcftools") {
    fmt.Println("ERROR: Repository did not download correctly. Please redownload from https://github.com/bweestyle/PathogenID")
    os.Exit(1)
  }
}

func MakeBinExecutable(cwd string) {
  executables := []string{"bwa", "samtools", "bcftools"}
  for i:=0;i<len(executables);i++{
    fmt.Println("chmod 755 " + cwd+"/"+executables[i])
    cmd := CreateCommand("chmod 755 " + cwd+"/bin/"+executables[i])
    cmd.Run()
  }
}

func PrepareBin(cwd string){
  CheckBin(cwd)
  MakeBinExecutable(cwd)
}



func MakeSamFile(cwd,reference string, numProcs, LN int, analysisFolder string, readFiles... string) string{
  folders := strings.Split(strings.TrimSuffix(readFiles[0],".fastq"),"/")
  filename := folders[len(folders)-1]
  samFile := analysisFolder + filename+".sam"
  var bwaMem *exec.Cmd
  if len(readFiles)==2{
    bwaMem = CreateCommand(cwd+"/bin/bwa mem -t " +strconv.Itoa(numProcs) + " " + reference + " " + readFiles[0] + " " + readFiles[1])
  } else {
    bwaMem = CreateCommand(cwd+"/bin/bwa mem -t " +strconv.Itoa(numProcs) + " " + reference + " " + readFiles[0])
  }
  OutputCommandToFile(bwaMem, samFile)
  CheckSamFile(samFile, LN)
  return samFile
}

func MakeBamFile(cwd, samFile string, numProcs int) string{
  samtoolsView := CreateCommand(cwd+"/bin/samtools view -@ " + strconv.Itoa(numProcs-1) + " -bS " + samFile)
  bamFile := strings.Split(samFile, ".")[0]+".bam"
  OutputCommandToFile(samtoolsView,bamFile)
  return bamFile
}

func SortBamFile(cwd, bamFile string, numProcs int) string {
  samtoolsSort:= CreateCommand(cwd+"/bin/samtools sort -@ " + strconv.Itoa(numProcs-1) + " " + bamFile)
  sortedBam := strings.Split(bamFile,".")[0]+".sorted.bam"
  OutputCommandToFile(samtoolsSort, sortedBam)
  return sortedBam
}

func AlignReads(cwd, reference string, readFiles []string, numProcs, LN int, pairedEnd bool, analysisFolder string) string{
  var samFile string
  if pairedEnd{
    samFile = MakeSamFile(cwd,reference, numProcs, LN, analysisFolder, readFiles[0], readFiles[1])
  } else {
    samFile = MakeSamFile(cwd,reference, numProcs, LN, analysisFolder, readFiles[0])
  }

  bamFile := MakeBamFile(cwd,samFile, numProcs)
  sortedBam := SortBamFile(cwd, bamFile, numProcs)
  return sortedBam
}

func MakeVCF(cwd, reference, sortedBam string)string {
  samtoolsMpileup:= CreateCommand(cwd+"/bin/samtools mpileup -v --reference "+reference+" "+sortedBam)
  vcfFile := strings.Split(sortedBam,".")[0]+".vcf"
  OutputCommandToFile(samtoolsMpileup, vcfFile)
  return vcfFile
}

func CallVCF(cwd, vcfFile string, numProcs int)string{
  bcfToolsCall := CreateCommand(cwd+"/bin/bcftools call --threads " + strconv.Itoa(numProcs-1) + " --ploidy 1 -c " + vcfFile)
  calledVCFFile := strings.Split(vcfFile,".")[0]+".called.vcf"
  OutputCommandToFile(bcfToolsCall, calledVCFFile)
  return calledVCFFile
}

func CallVariants(cwd, reference, sortedBam string, numProcs int) string{
  vcfFile:= MakeVCF(cwd, reference, sortedBam)
  calledVCFFile := CallVCF(cwd, vcfFile, numProcs)
  return calledVCFFile
}

func MakeAnalysisFolder(readFiles []string) string{
  folders := strings.Split(readFiles[0],"/")
  dataLocation :=""
  for i:=1;i<len(folders)-1;i++{
    dataLocation = dataLocation+"/"+folders[i]
  }
  analysisFolder:=dataLocation+"/Analysis/"
  RunCommand(CreateCommand("mkdir "+analysisFolder))
  return analysisFolder
}

func GetSampleData(cwd string) []string{
  path := cwd+"/SampleData/"
  var readFiles []string
  forward := UnzipFile(path+"test_data.fastq.gz")
  reverse := UnzipFile(path+"test_rev.fastq.gz")
  readFiles = append(readFiles, forward)
  readFiles = append(readFiles, reverse)
  return readFiles
}

func pwd()string{
  pwd:= CreateCommand("pwd")
  cwd := WriteOutputToString(pwd)
  return cwd
}

func OnlyAlign(readFiles []string) string{
  cwd := pwd()
  PrepareBin(cwd)

  numProcs := runtime.NumCPU()
  if numProcs >1 { //use all but one Processor just in case
    numProcs -= 1
  }

  analysisFolder := MakeAnalysisFolder(readFiles)

  //identify oraganism
  reference,LN := HandleReference(cwd)
  IndexReference(cwd, reference)

  var samFile string
  if len(readFiles)==2{
    samFile = MakeSamFile(cwd,reference, numProcs, LN, analysisFolder, readFiles[0], readFiles[1])
  } else {
    samFile = MakeSamFile(cwd,reference, numProcs, LN, analysisFolder, readFiles[0])
  }
  return samFile
}

func RunPipeline(readFiles []string) {
  pwd:= CreateCommand("pwd")
  cwd := WriteOutputToString(pwd)
  PrepareBin(cwd)

  numProcs := runtime.NumCPU()
  if numProcs >1 { //use all but one Processor just in case
    numProcs -= 1
  }
  pairedEnd := (len(readFiles)==2)
  //get Filename reads()

  if len(readFiles)==0{
    //Use sample data
    readFiles = GetSampleData(cwd)
  }

  analysisFolder := MakeAnalysisFolder(readFiles)

  //identify oraganism
  reference,LN := HandleReference(cwd)
  IndexReference(cwd, reference)
  sortedBam:= AlignReads(cwd, reference, readFiles, numProcs, LN, pairedEnd, analysisFolder)

  calledVCFFile := CallVariants(cwd, reference, sortedBam, numProcs)
  fmt.Println(calledVCFFile, "Created")

  //ProcessVCF()
}

func main(){
  var s []string
  RunPipeline(s)
}
