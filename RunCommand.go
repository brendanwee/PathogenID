package main

import(
  "strings"
  "fmt"
  "os"
  "os/exec"
)


func CreateCommand(input string) *exec.Cmd{
  items := strings.Fields(input)
  command := items[0]
  args := items[1:]
  fmt.Printf("creating %s command \n", command)

  cmd := exec.Command(command, args...)
  return cmd
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

func main() {
  CreateCommand("ls")
  CreateCommand("pwd")
  OutputCommandToFile(CreateCommand("bwa mem Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa ERR027083_1.fastq"), "test.sam")
}
