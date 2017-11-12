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
  fmt.Printf("creating %s %s command \n", command, args[0])

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
    fmt.Println("BWA mem Failed")
    os.Exit(1)
  }
  genomeLength, err := strconv.Atoi((strings.Split(words[2],":")[1])) //the third word should be LN:num, get num
  if genomeLength!=LN{
    fmt.Println("BWA mem Failed")
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

//downlaods the Tuberculosis reference genome and returns the filename as well as its genome length
func DownloadAndIndexReference() (string, int) {
  downloadReference:= CreateCommand("curl ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/dna/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa.gz")
  reference := "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa.gz"
  OutputCommandToFile(downloadReference, reference)
  CheckReferenceFile(reference)
  if strings.HasSuffix(reference, ".gz") {
    gunzip := CreateCommand("gunzip " + reference)
    RunCommand(gunzip)
    reference = "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa"
  }
  LN:=GetGenomeLength(reference)
  //index genome
  bwaIndex := CreateCommand("bwa index "+reference)
  RunCommand(bwaIndex)
  return reference, LN
}


func main() {
  numProcs := runtime.NumCPU()
  if numProcs >1 { //use all but one Processor just in case
    numProcs -= 1
  }
  reference, LN := DownloadAndIndexReference()

  //get Filename reads()
  forwardReads := "A70376.fastq"
  reverseReads := "A70376_2.fastq"
  samFile := strings.Split(forwardReads,".")[0]+".sam"
  bwaMem := CreateCommand("bwa mem -t " +strconv.Itoa(numProcs) + " " + reference + " " + forwardReads + " " + reverseReads)
  OutputCommandToFile(bwaMem, samFile)
  CheckSamFile(samFile, LN)

  samtoolsView := CreateCommand("samtools view -@ " + strconv.Itoa(numProcs-1) + " -bS " + samFile)
  bamFile := strings.Split(samFile, ".")[0]+".bam"
  OutputCommandToFile(samtoolsView,bamFile)

  samtoolsSort:= CreateCommand("samtools sort -@ " + strconv.Itoa(numProcs-1) + " " + bamFile)
  sortedBam := strings.Split(bamFile,".")[0]+".sorted.bam"
  OutputCommandToFile(samtoolsSort, sortedBam)


  samtoolsMpileup:= CreateCommand("samtools mpileup -v --reference "+reference+" "+sortedBam)
  vcfFile := strings.Split(sortedBam,".")[0]+".vcf"
  OutputCommandToFile(samtoolsMpileup, vcfFile)
  bcfToolsCall := CreateCommand("bcftools call --threads " + strconv.Itoa(numProcs-1) + " --ploidy 1 -c " + vcfFile)
  calledVcfFile := strings.Split(sortedBam,".")[0]+".called.vcf"
  OutputCommandToFile(bcfToolsCall, calledVcfFile)
  //ProcessVCF()
}
