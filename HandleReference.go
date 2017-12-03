package main

import(
  "strings"
  "fmt"
  "os"
  "bufio"
  "os/exec"
)

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


func DownloadReference() string{
  downloadReference:= CreateCommand("curl ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/dna/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa.gz")
  reference := "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa.gz"
  OutputCommandToFile(downloadReference, reference)
  CheckReferenceFile(reference)
  reference = UnzipFile(reference)
  return reference
}

func ReferenceExists() (bool, string){
  ls:=exec.Command("ls")
  contents := WriteOutputToString(ls)
  if !strings.Contains(contents, "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa"){
    return false, ""
  }
  return true, contents
}

func HandleReference()  (string, int){
  var reference string
  var LN int
  if exists,contents:= ReferenceExists(); exists {
    reference = PrepareReference(contents)
  } else {
    reference = DownloadReference()
  }
  LN = GetGenomeLength(reference)
  return reference, LN
}

func PrepareReference(contents string) string{
  var reference string
  if strings.Contains(contents, "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa.gz"){
    reference = UnzipFile("Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa")
  } else if strings.Contains(contents, "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa"){
    reference = "Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa"
  } else {
    fmt.Println("ERROR: reference not detected! Downloading Reference Genome")
    reference = DownloadReference()
  }
  return reference
}
