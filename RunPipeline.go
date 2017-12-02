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


func CreateCommand(input string) *exec.Cmd{

  items := strings.Fields(input)
  command := items[0]
  args := items[1:]
  if len(args)>0{
    fmt.Printf("creating %s command \n", input)
  } else {
    fmt.Printf("creating %s command \n", input)
  }


  cmd := exec.Command(command, args...)
  return cmd
}

func RunCommand(cmd *exec.Cmd){
  cmd.Run()
  cmd.Wait()
}

func OutputCommandToFile(cmd *exec.Cmd, filename string){
  file, err := os.Create(filename)
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }
  cmd.Stdout = file
  cmd.Run()
  cmd.Wait()
}

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
  words := strings.Split(scanner.Text(),"	")
  if len(words)<3{ //check that there are at least 3 items in the first line
    fmt.Println("Less than three items in the first line. BWA mem Failed")
    os.Exit(1)
  }
  genomeLength, err := strconv.Atoi((strings.Split(words[2],":")[1])) //the third word should be LN:num, get num
  if genomeLength!=LN{
    fmt.Println("Genome length:", genomeLength, "!= ", LN, "BWA mem Failed")
    os.Exit(1)
  }
}

//takes a genomeFile and returns the amount of nucleotides in it.
func GetGenomeLength(genomeFile string) int {
  file,err := os.Open(genomeFile)
  CheckError(err)
  scanner := bufio.NewScanner(file)
  var LN int
  for scanner.Scan(){
    if scanner.Text()[0]=='>' || scanner.Text()[0]=='@' {
      continue
    } else if scanner.Text()[0]=='+' { //skip two lines
      scanner.Scan()// quality line
      continue
    } else {
      LN += len(scanner.Text())
    }
  }
  return LN
}

func CheckReferenceFile(reference string) {
  md5 := MD5(reference)
  errorMessage := `ERROR: Reference file did not download properly. Please check
                   your internet connection and try again`
  if strings.HasSuffix(reference, ".gz"){
    if md5 != "c34fb6593a6cbdbcfc0ac8d0c7db58ee" {
      fmt.Println(errorMessage)
      os.Exit(1)
    }
  } else {
    if md5 != "8c6a53ab340a9429c0db9a30801235c4" {
      fmt.Println(errorMessage)
      os.Exit(1)
    }
  }
}

func MD5(filename string) string {
  md5 := CreateCommand("md5 "+filename)
  OutputCommandToFile(md5, filename+".md5")
  file,err:=os.Open(filename+".md5")
  CheckError(err)
  scanner := bufio.NewScanner(file)
  scanner.Scan()
  fmt.Println(scanner.Text())
  md5Text := strings.Split(scanner.Text(), " ")[3]
  return md5Text
}

func UnzipFile(file string) string {
  if strings.HasSuffix(file, ".gz") {
    gunzip := CreateCommand("gunzip " + file)
    RunCommand(gunzip)
  }
  return strings.TrimSuffix(file, ".gz")
}

func DownloadReference() string{
  downloadReference:= CreateCommand("curl ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/dna/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa.gz")
  reference := "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa.gz"
  OutputCommandToFile(downloadReference, reference)
  CheckReferenceFile(reference)
  reference = UnzipFile(reference)
  return reference
}

func IndexReference(cwd,reference string) {
  //index genome
  bwaIndex := CreateCommand(cwd+"/bin/bwa index "+reference)
  RunCommand(bwaIndex)
}

func WriteOutputToString(cmd *exec.Cmd) string{
  out,err := cmd.Output()
  CheckError(err)
  return string(out)[:len(out)-1]
}

func CheckBin(cwd string){
  fmt.Println(cwd)
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

func ReferenceExists() (bool, string){
  ls:=exec.Command("ls")
  contents := WriteOutputToString(ls)
  if !strings.Contains(contents, "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa"){
    return false, ""
  }
  return true, contents
}

func RetrieveReference(contents string) string{
  var reference string
  if strings.Contains(contents, "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa"){
    reference = "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa"
  } else if strings.Contains(contents, "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa.gz"){
    reference = UnzipFile("Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa.gz")
  } else {
    fmt.Println("ERROR: reference not detected! Downloading Reference Genome")
    reference = DownloadReference()
  }
  return reference
}

func GetReference() (string, int) {
  var reference string
  var LN int
  if exists,contents:= ReferenceExists(); exists {
    reference = RetrieveReference(contents)
  } else {
    reference = DownloadReference()
  }
  LN = GetGenomeLength(reference)
  return reference,LN
}

func PrepareBin(cwd string){
  CheckBin(cwd)
  MakeBinExecutable(cwd)
}

func MakeSamFile(cwd,reference string, numProcs, LN int, readFiles... string) string{
  samFile := strings.Split(readFiles[0],".")[0]+".sam"
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

func AlignReads(cwd, reference string, readFiles []string, numProcs, LN int, pairedEnd bool) string{
  var samFile string
  if pairedEnd{
    samFile = MakeSamFile(cwd,reference, numProcs, LN, readFiles[0], readFiles[1])
  } else {
    samFile = MakeSamFile(cwd,reference, numProcs, LN, readFiles[0])
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

func HandleReference(cwd string)  (string, int){
  reference, LN := GetReference()
  IndexReference(cwd,reference)
  return reference, LN
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

  //identify oraganism
  reference,LN := HandleReference(cwd)

  sortedBam:= AlignReads(cwd, reference, readFiles, numProcs, LN, pairedEnd)

  calledVCFFile := CallVariants(cwd, reference, sortedBam, numProcs)
  fmt.Println(calledVCFFile, "Created")

  //ProcessVCF()
}

func main(){
  var s []string
  RunPipeline(s)
}
