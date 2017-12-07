package main

func InitializeGene() {
	rpoB.boundary[0] = 759807
	rpoB.boundary[1] = 763325
	rpoB.name = "rpoB"
	rpoB.orientation = true
	rpoB.drug = append(rpoB.drug, "RIF")
	allGenes = append(allGenes, rpoB)

	embB.boundary[0] = 4246514
	embB.boundary[1] = 4249810
	embB.name = "embB"
	embB.orientation = true
	embB.drug = append(embB.drug, "RIF")
	embB.drug = append(embB.drug, "INH")
	embB.drug = append(embB.drug, "EMB")
	allGenes = append(allGenes, embB)

	fbpC.boundary[0] = 156578
	fbpC.boundary[1] = 157600
	fbpC.name = "fbpC"
	fbpC.orientation = false
	fbpC.drug = append(fbpC.drug, "INH")
	allGenes = append(allGenes, fbpC)

	Rv0340.boundary[0] = 408634
	Rv0340.boundary[1] = 409173
	Rv0340.name = "Rv0340"
	Rv0340.orientation = true
	Rv0340.drug = append(Rv0340.drug, "INH")
	Rv0340.drug = append(Rv0340.drug, "EMB")
	allGenes = append(allGenes, Rv0340)

	iniB.boundary[0] = 409362
	iniB.boundary[1] = 410801
	iniB.name = "iniB"
	iniB.orientation = true
	iniB.drug = append(iniB.drug, "INH")
	Rv0340.drug = append(Rv0340.drug, "EMB")
	allGenes = append(allGenes, iniB)

	iniA.boundary[0] = 410838
	iniA.boundary[1] = 412760
	iniA.name = "iniA"
	iniA.orientation = true
	iniA.drug = append(iniA.drug, "INH")
	iniA.drug = append(iniA.drug, "EMB")
	allGenes = append(allGenes, iniA)

	iniC.boundary[0] = 412757
	iniC.boundary[1] = 414238
	iniC.name = "iniC"
	iniC.orientation = true
	iniC.drug = append(iniC.drug, "INH")
	iniC.drug = append(iniC.drug, "EMB")
	allGenes = append(allGenes, iniC)

	mabA.boundary[0] = 1673440
	mabA.boundary[1] = 1674183
	mabA.name = "mabA"
	mabA.orientation = true
	mabA.drug = append(mabA.drug, "INH")
	mabA.drug = append(mabA.drug, "ETH")
	allGenes = append(allGenes, mabA)

	inhA.boundary[0] = 1674202
	inhA.boundary[1] = 1675011
	inhA.name = "inhA"
	inhA.orientation = true
	inhA.drug = append(inhA.drug, "INH")
	inhA.drug = append(inhA.drug, "ETH")
	allGenes = append(allGenes, inhA)

	Rv1592c.boundary[0] = 1792400
	Rv1592c.boundary[1] = 1793740
	Rv1592c.name = "Rv1592c"
	Rv1592c.orientation = false
	Rv1592c.drug = append(Rv1592c.drug, "INH")
	allGenes = append(allGenes, Rv1592c)

	Rv1772.boundary[0] = 2006636
	Rv1772.boundary[1] = 2006947
	Rv1772.name = "Rv1772"
	Rv1772.orientation = true
	Rv1772.drug = append(Rv1772.drug, "INH")
	allGenes = append(allGenes, Rv1772)

	ndh.boundary[0] = 2101651
	ndh.boundary[1] = 2103042
	ndh.name = "ndh"
	ndh.orientation = false
	ndh.drug = append(ndh.drug, "INH")
	allGenes = append(allGenes, ndh)

	katG.boundary[0] = 2153889
	katG.boundary[1] = 2156111
	katG.name = "katG"
	katG.orientation = false
	katG.drug = append(katG.drug, "INH")
	allGenes = append(allGenes, katG)

	furA.boundary[0] = 2156149
	furA.boundary[1] = 2156592
	furA.name = "furA"
	furA.orientation = false
	furA.drug = append(furA.drug, "INH")
	allGenes = append(allGenes, furA)

	srmR.boundary[0] = 2515304
	srmR.boundary[1] = 2516548
	srmR.name = "srmR"
	srmR.orientation = true
	srmR.drug = append(srmR.drug, "INH")
	allGenes = append(allGenes, srmR)

	fabD.boundary[0] = 2516787
	fabD.boundary[1] = 2517695
	fabD.name = "fabD"
	fabD.orientation = true
	fabD.drug = append(fabD.drug, "INH")
	allGenes = append(allGenes, fabD)

	kasA.boundary[0] = 2518115
	kasA.boundary[1] = 2519365
	kasA.name = "kasA"
	kasA.orientation = true
	kasA.drug = append(kasA.drug, "INH")
	allGenes = append(allGenes, kasA)

	accD6.boundary[0] = 2520743
	accD6.boundary[1] = 2522164
	accD6.name = "accD6"
	accD6.orientation = true
	accD6.drug = append(accD6.drug, "INH")
	allGenes = append(allGenes, accD6)

	oxyR.boundary[0] = 2725571
	oxyR.boundary[1] = 2726087
	oxyR.name = "oxyR"
	oxyR.orientation = false
	oxyR.drug = append(oxyR.drug, "INH")
	allGenes = append(allGenes, oxyR)

	aphC.boundary[0] = 2726193
	aphC.boundary[1] = 2726780
	aphC.name = "aphC"
	aphC.orientation = true
	aphC.drug = append(aphC.drug, "INH")
	allGenes = append(allGenes, aphC)

	efpA.boundary[0] = 3153039
	efpA.boundary[1] = 3154631
	efpA.name = "efpA"
	efpA.orientation = false
	efpA.drug = append(efpA.drug, "INH")
	allGenes = append(allGenes, efpA)

	fadE24.boundary[0] = 3505363
	fadE24.boundary[1] = 3506769
	fadE24.name = "fadE24"
	fadE24.orientation = false
	fadE24.drug = append(fadE24.drug, "INH")
	allGenes = append(allGenes, fadE24)

	nhoA.boundary[0] = 4007331
	nhoA.boundary[1] = 4008182
	nhoA.name = "nhoA"
	nhoA.orientation = true
	nhoA.drug = append(nhoA.drug, "INH")
	allGenes = append(allGenes, nhoA)

	pncA.boundary[0] = 2288681
	pncA.boundary[1] = 2289241
	pncA.name = "pncA"
	pncA.orientation = false
	pncA.drug = append(pncA.drug, "PZA")
	allGenes = append(allGenes, pncA)

	rpsL.boundary[0] = 781560
	rpsL.boundary[1] = 781934
	rpsL.name = "rpsL"
	rpsL.orientation = true
	rpsL.drug = append(rpsL.drug, "SM")
	allGenes = append(allGenes, rpsL)

	rrs.boundary[0] = 1471846
	rrs.boundary[1] = 1473382
	rrs.name = "rrs"
	rrs.orientation = false
	rrs.drug = append(rrs.drug, "SM")
	rrs.drug = append(rrs.drug, "AMI")
	allGenes = append(allGenes, rrs)

	gidB.boundary[0] = 2288681
	gidB.boundary[1] = 2289241
	gidB.name = "gidB"
	gidB.orientation = false
	gidB.drug = append(gidB.drug, "SM")
	allGenes = append(allGenes, gidB)

	tlyA.boundary[0] = 1917940
	tlyA.boundary[1] = 1918746
	tlyA.name = "tlyA"
	tlyA.orientation = true
	tlyA.drug = append(tlyA.drug, "AMI")
	allGenes = append(allGenes, tlyA)

	embR.boundary[0] = 1416181
	embR.boundary[1] = 1417347
	embR.name = "embR"
	embR.orientation = false
	embR.drug = append(embR.drug, "EMB")
	allGenes = append(allGenes, embR)

	Rv3124.boundary[0] = 3489506
	Rv3124.boundary[1] = 3490375
	Rv3124.name = "Rv3124"
	Rv3124.orientation = true
	Rv3124.drug = append(Rv3124.drug, "EMB")
	allGenes = append(allGenes, Rv3124)

	Rv3125c.boundary[0] = 3490476
	Rv3125c.boundary[1] = 3491651
	Rv3125c.name = "Rv3125c"
	Rv3125c.orientation = false
	Rv3125c.drug = append(Rv3125c.drug, "EMB")
	allGenes = append(allGenes, Rv3125c)

	Rv3264c.boundary[0] = 3644898
	Rv3264c.boundary[1] = 3645977
	Rv3264c.name = "Rv3264c"
	Rv3264c.orientation = false
	Rv3264c.drug = append(Rv3264c.drug, "EMB")
	allGenes = append(allGenes, Rv3264c)

	Rv3266c.boundary[0] = 3646895
	Rv3266c.boundary[1] = 3647809
	Rv3266c.name = "Rv3266c"
	Rv3266c.orientation = false
	Rv3266c.drug = append(Rv3266c.drug, "EMB")
	allGenes = append(allGenes, Rv3266c)

	embC.boundary[0] = 4239863
	embC.boundary[1] = 4243147
	embC.name = "embC"
	embC.orientation = true
	embC.drug = append(embC.drug, "EMB")
	allGenes = append(allGenes, embC)

	embA.boundary[0] = 4243233
	embA.boundary[0] = 4246517
	embA.name = "embA"
	embA.orientation = true
	embA.drug = append(embA.drug, "EMB")
	allGenes = append(allGenes, embA)

	// ethA,etaA,Rv3854c
	ethA.boundary[0] = 4326004
	ethA.boundary[0] = 4327473
	ethA.name = "ethA"
	ethA.orientation = false
	ethA.drug = append(ethA.drug, "ETH")
	allGenes = append(allGenes, ethA)

	gyrB.boundary[0] = 5240
	gyrB.boundary[0] = 7267
	gyrB.name = "gyrB"
	gyrB.orientation = true
	gyrB.drug = append(gyrB.drug, "FLQ")
	allGenes = append(allGenes, gyrB)

	gyrA.boundary[0] = 7302
	gyrA.boundary[0] = 9818
	gyrA.name = "gyrA"
	gyrA.orientation = true
	gyrA.drug = append(gyrA.drug, "FLQ")
	allGenes = append(allGenes, gyrA)

	thyA.boundary[0] = 3073680
	thyA.boundary[0] = 3074471
	thyA.name = "thyA"
	thyA.orientation = false
	thyA.drug = append(thyA.drug, "PAS")
	allGenes = append(allGenes, thyA)

}

func InitializeDrug() {
	RIF.name = "RIF"
	allDrug = append(allDrug, RIF)

	INH.name = "INH"
	allDrug = append(allDrug, INH)

	PZA.name = "PZA"
	allDrug = append(allDrug, PZA)

	SM.name = "SM"
	allDrug = append(allDrug, SM)

	AMI.name = "AMI"
	allDrug = append(allDrug, AMI)

	EMB.name = "EMB"
	allDrug = append(allDrug, EMB)

	ETH.name = "ETH"
	allDrug = append(allDrug, ETH)

	FLQ.name = "FLQ"
	allDrug = append(allDrug, FLQ)

	PAS.name = "PAS"
	allDrug = append(allDrug, PAS)
}
