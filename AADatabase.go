package main

type AA struct{
  name string
  polarity string
  charge int
  codon []string
}

var aminoAcid []AA
var ala, arg, asn, asp, cys, glu, gln, gly, his, ile, leu, lys, met, phe, pro,
ser, thr, trp, tyr,val, stp AA
var CodonTable map[string]AA

// Initiate the codon table
func InitializeCodonTable() {
  CodonTable = make(map[string]AA)
  for i := range aminoAcid {
    for j := range aminoAcid[i].codon {
      CodonTable[aminoAcid[i].codon[j]] = aminoAcid[i]
    }
  }
}

// Initiate amino acid database
func InitializeAA() {
  ala.name = "ala"
  ala.polarity = "nonpolar"
  ala.charge = 0
  ala.codon = append(ala.codon, "CGA")
  ala.codon = append(ala.codon, "CGG")
  ala.codon = append(ala.codon, "CGT")
  ala.codon = append(ala.codon, "CGC")
  aminoAcid = append(aminoAcid, ala)

  arg.name = "arg"
  arg.polarity = "basic polar"
  arg.charge = 1
  arg.codon = append(arg.codon, "GCA")
  arg.codon = append(arg.codon, "GCG")
  arg.codon = append(arg.codon, "GCT")
  arg.codon = append(arg.codon, "GCC")
  arg.codon = append(arg.codon, "TCT")
  arg.codon = append(arg.codon, "TCC")
  aminoAcid = append(aminoAcid, arg)

  asn.name = "asn"
  asn.polarity = "polar"
  asn.charge = 0
  asn.codon = append(asn.codon, "TTA")
  asn.codon = append(asn.codon, "TTG")
  aminoAcid = append(aminoAcid, asn)

  asp.name = "asp"
  asp.polarity = "acidic polar"
  asp.charge = -1
  asp.codon = append(asp.codon, "CTA")
  asp.codon = append(asp.codon, "CTG")
  aminoAcid = append(aminoAcid, asp)

  cys.name = "cys"
  cys.polarity = "nonpolar"
  cys.charge = 0
  cys.codon = append(cys.codon, "ACA")
  cys.codon = append(cys.codon, "ACG")
  aminoAcid = append(aminoAcid, cys)

  glu.name = "glu"
  glu.polarity = "acidic polar"
  glu.charge = -1
  glu.codon = append(glu.codon, "CTT")
  glu.codon = append(glu.codon, "CTC")
  aminoAcid = append(aminoAcid, glu)

  gln.name = "gln"
  gln.polarity = "polar"
  gln.charge = 0
  gln.codon = append(gln.codon, "GTT")
  gln.codon = append(gln.codon, "GTC")
  aminoAcid = append(aminoAcid, gln)

  gly.name = "gly"
  gly.polarity = "nonpolar"
  gly.charge = 0
  gly.codon = append(gly.codon, "CCA")
  gly.codon = append(gly.codon, "CCG")
  gly.codon = append(gly.codon, "CCT")
  gly.codon = append(gly.codon, "CCC")
  aminoAcid = append(aminoAcid, gly)

  his.name = "his"
  his.polarity = "basic polar"
  his.charge = 0
  his.codon = append(his.codon, "GTA")
  his.codon = append(his.codon, "GTG")
  aminoAcid = append(aminoAcid, his)

  ile.name = "ile"
  ile.polarity = "nonpolar"
  ile.charge = 0
  ile.codon = append(ile.codon, "TAA")
  ile.codon = append(ile.codon, "TAG")
  ile.codon = append(ile.codon, "TAT")
  aminoAcid = append(aminoAcid, ile)

  leu.name = "leu"
  leu.polarity = "nonpolar"
  leu.charge = 0
  leu.codon = append(leu.codon, "GAA")
  leu.codon = append(leu.codon, "GAG")
  leu.codon = append(leu.codon, "GAT")
  leu.codon = append(leu.codon, "GAC")
  leu.codon = append(leu.codon, "AAT")
  leu.codon = append(leu.codon, "AAC")
  aminoAcid = append(aminoAcid, leu)

  lys.name = "lys"
  lys.polarity = "basic polar"
  lys.charge = 1
  lys.codon = append(lys.codon, "TTT")
  lys.codon = append(lys.codon, "TTC")
  aminoAcid = append(aminoAcid, lys)

  met.name = "met"
  met.polarity = "nonpolar"
  met.charge = 0
  met.codon = append(met.codon, "TAC")
  aminoAcid = append(aminoAcid, met)

  phe.name = "phe"
  phe.polarity = "nonpolar"
  phe.charge = 0
  phe.codon = append(phe.codon, "AAA")
  phe.codon = append(phe.codon, "AAG")
  aminoAcid = append(aminoAcid, phe)

  pro.name = "pro"
  pro.polarity = "nonpolar"
  pro.charge = 0
  pro.codon = append(pro.codon, "GGA")
  pro.codon = append(pro.codon, "GGG")
  pro.codon = append(pro.codon, "GGT")
  pro.codon = append(pro.codon, "GGC")
  aminoAcid = append(aminoAcid, pro)

  ser.name = "ser"
  ser.polarity = "polar"
  ser.charge = 0
  ser.codon = append(ser.codon, "AGA")
  ser.codon = append(ser.codon, "AGG")
  ser.codon = append(ser.codon, "AGT")
  ser.codon = append(ser.codon, "AGC")
  ser.codon = append(ser.codon, "TCA")
  ser.codon = append(ser.codon, "TCG")
  aminoAcid = append(aminoAcid, ser)

  thr.name = "thr"
  thr.polarity = "polar"
  thr.charge = 0
  thr.codon = append(thr.codon, "TGA")
  thr.codon = append(thr.codon, "TGG")
  thr.codon = append(thr.codon, "TGT")
  thr.codon = append(thr.codon, "TGC")
  aminoAcid = append(aminoAcid, thr)

  trp.name = "trp"
  trp.polarity = "nonpolar"
  trp.charge = 0
  trp.codon = append(trp.codon, "ACC")
  aminoAcid = append(aminoAcid, trp)

  tyr.name = "tyr"
  tyr.polarity = "polar"
  tyr.charge = 0
  tyr.codon = append(tyr.codon, "ATA")
  tyr.codon = append(tyr.codon, "ATG")
  aminoAcid = append(aminoAcid, tyr)

  val.name = "val"
  val.polarity = "nonpolar"
  val.charge = 0
  val.codon = append(val.codon, "CAA")
  val.codon = append(val.codon, "CAG")
  val.codon = append(val.codon, "CAT")
  val.codon = append(val.codon, "CAC")
  aminoAcid = append(aminoAcid, val)

  stp.name = "stop"
  stp.codon = append(stp.codon, "ATT")
  stp.codon = append(stp.codon, "ATC")
  stp.codon = append(stp.codon, "ACT")
  aminoAcid = append(aminoAcid, stp)
}
