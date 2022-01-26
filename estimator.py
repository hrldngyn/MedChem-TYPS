from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Fragments
from rdkit.Chem import Descriptors

# import matplotlib.pyplot as pyplot
# import matplotlib.image as mpimg

import requests


from collections.abc import Iterable

import os

os.chdir(os.path.join(__file__, '..'))

SP3 = Chem.rdchem.HybridizationType.SP3

groups = {
        "amide": "C(=O)-N",
        "amine_primary": "[NH2,nH2]",
        "amine_secondary": "[NH1,nH1]",
        "amine_tertiary": "[NH,nH]",
        "aniline": "c-[NX3;!$(N=*)]",
        #"beta_lactam": "N1C(=O)CC1",
        "oxygen": "[#8]",
        "COOH": "C(=O)[O;H1,-]",
        #"COOH_al": "C-C(=O)[O;H1,-]",
        #"COOH_ar": "c-C(=O)[O;H1,-]",
        "carbonyl": "[CX3]=[OX1]" ,
        "C_aromatic": "[$([cX3](:*):*),$([cX2+](:*):*)]",
        "ester": "[#6][CX3](=O)[OX2H0][#6]",
        "phenol": "[OX2H]-c1ccccc1",
        "pyridine": "n1ccccc1",
        "sulfur": "[S]",
        "sulfonamide": "N-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]",
        "tetrazole": "c1nnnn1",
        "phosphorous": "[P]",
        "halogen": "[F,Cl,Br,I]",
}

#chem.fragments docs
acids = {
    "COOH_al": "C-C(=O)[O;H1,-]",
    "COOH_ar": "c-C(=O)[O;H1,-]",

}

bases = {

}

def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el

def getEWGs(mol):
    EWGindices = []
    for key in groups:
        patt = Chem.MolFromSmarts(groups[key])
        EWGindices.append(mol.GetSubstructMatches(patt))
    print(EWGindices)
    return set(flatten(EWGindices))

def draw(m):    
    d = rdMolDraw2D.MolDraw2DCairo(1000, 1000) # or MolDraw2DSVG to get SVGs

    d.drawOptions().addStereoAnnotation = True
    d.drawOptions().addAtomIndices = True
    d.DrawMolecule(m)
    d.FinishDrawing()
    d.WriteDrawingText('atom_annotation_1.png')

#ALGORITHM CODE

def setMol(smiles):
    m = Chem.MolFromSmiles(smiles)

#implement hierarchical acid/base identification, i.e. benzoic acid recognized before normal acid, amide recognized before amine?
    #id groups, exclude
        #check if acidic or basic based on above dictionaries
            #estimate(index, true)
    #also create array of groups with output estimation data

#make below nicer with enums and generator expressions
#basic: bool

def isAliphaticCarbon(atom):
    return (atom.GetHybridization() == SP3 and atom.GetSymbol() == 'C')

#deprecated
def estimate(index, ewgids, parent, basic=False): #index is the 0 atom to draw distance from
    searched = [index]
    b = int(basic) + 1
    subtracted = []
    subtraction = 0

    #find carbonyl locations
        #check if aliphatic carbon
    #if in searched, check that path is shortest

    for c in m.GetAtomWithIdx(index).GetNeighbors(): #for the neighbors of the 0 atom (c is 0 away)
        print("Searching " + str(c.GetIdx()))
        searched.append(c.GetIdx())
        # print(c.GetHybridization() == SP3)
        # print(c.GetSymbol() == 'C')
        if(c.GetHybridization() == SP3 and c.GetSymbol() == 'C'):
            #for d in (d for d in c.GetNeighbors() if d.GetIdx() not in searched): #d is 1 away
            for d in c.GetNeighbors():
                if d.GetIdx() not in searched:   
                    print("Searching " + str(d.GetIdx()))
                    searched.append(d.GetIdx())
                    if(d.GetHybridization() == SP3 and d.GetSymbol() == 'C'):
                        for e in (e for e in d.GetNeighbors() if e.GetIdx() not in searched): #e is 2 away
                            print("Searching " + str(e.GetIdx()))
                            searched.append(e.GetIdx())
                            if(e.GetHybridization() == SP3 and e.GetSymbol() == 'C'):
                                for f in (f for f in e.GetNeighbors() if f.GetIdx() not in searched): #f is 3 away
                                    print("Searching " + str(f.GetIdx()))
                                    searched.append(f.GetIdx())
                                    if(f.GetHybridization() == SP3 and f.GetSymbol() == 'C'):
                                        for g in (g for g in f.GetNeighbors() if g.GetIdx() not in searched): #g is 4 away
                                            print("Searching " + str(g.GetIdx()))
                                            searched.append(g.GetIdx())
                                            if g.GetIdx() in ewgids:
                                                if (g.GetIsAromatic()):
                                                    subtraction += 0.075*b
                                                else:
                                                    subtraction += 0.15*b
                                                subtracted.append(g.GetIdx())
                                    elif(f.GetIdx() in ewgids):
                                        if(f.GetIsAromatic()):
                                            subtraction += 0.2*b
                                        else:
                                            subtraction += 0.4*b
                                        subtracted.append(f.GetIdx())
                            elif e.GetIdx() in ewgids:
                                if (e.GetIsAromatic()):
                                    subtraction += 0.45*b
                                else:
                                    subtraction += 0.9*b
                                subtracted.append(e.GetIdx())
                    elif d.GetIdx() in ewgids:
                        if (d.GetIsAromatic()):
                            subtraction += 0.45*b
                        else:
                            subtraction += 0.9*b
                        subtracted.append(d.GetIdx())
    #add in check for amine in ring
    if m.GetAtomWithIdx(index).GetSymbol() == 'N' and m.GetAtomWithIdx(index).IsInRing():
        subtraction += 1


    print("Setting prop for index " + str(index) + " as " + str(round(parent - subtraction, 4)))
    print(m.GetAtomWithIdx(index).SetProp("atomNote", "\nEst. pKa = " + str(round(parent - subtraction, 4))))
    print("Subtracted indices from " + str(index) + ": " + str(subtracted))
    return round(parent - subtraction, 4)

def BFestimate(m, index, ewgids, parent, basic=False): #breadth-first algorithm
    searched = [index]
    bval = int(basic) + 1
    subtracted = []
    subtraction = 0
    zeroaway, oneaway, twoaway, threeaway, fouraway = [], [], [], [], []

    for a in m.GetAtomWithIdx(index).GetNeighbors():
        searched.append(a.GetIdx())
        if(isAliphaticCarbon(a)):
            zeroaway.append(a.GetIdx())

    for a in zeroaway:
        for b in m.GetAtomWithIdx(a).GetNeighbors():
            if b.GetIdx() not in searched:
                searched.append(b.GetIdx())
                if(isAliphaticCarbon(b)):
                    oneaway.append(b.GetIdx())
                elif b.GetIdx() in ewgids:
                    if (b.GetIsAromatic()):
                        subtraction += 0.45*bval
                    else:
                        subtraction += 0.9*bval
                    subtracted.append(b.GetIdx())
    print("oneaway")
    print(oneaway)
    for b in oneaway:
        for c in m.GetAtomWithIdx(b).GetNeighbors():
            if c.GetIdx() not in searched:
                searched.append(c.GetIdx())
                if(isAliphaticCarbon(c)):
                    twoaway.append(c.GetIdx())
                elif c.GetIdx() in ewgids:
                    if (c.GetIsAromatic()):
                        subtraction += 0.45*bval
                    else:
                        subtraction += 0.9*bval
                    subtracted.append(c.GetIdx())
    print("twoaway")
    print(twoaway)
    for c in twoaway:
        for d in m.GetAtomWithIdx(c).GetNeighbors():
            if d.GetIdx() not in searched:
                searched.append(d.GetIdx())
                if(isAliphaticCarbon(d)):
                    threeaway.append(d.GetIdx())
                elif d.GetIdx() in ewgids:
                    if (d.GetIsAromatic()):
                        subtraction += 0.2*bval
                    else:
                        subtraction += 0.4*bval
                    subtracted.append(d.GetIdx())
    for d in threeaway:
        for e in m.GetAtomWithIdx(d).GetNeighbors():
            if e.GetIdx() not in searched:
                searched.append(e.GetIdx())
                if(isAliphaticCarbon(e)):
                    fouraway.append(e.GetIdx())
                elif e.GetIdx() in ewgids:
                    if (e.GetIsAromatic()):
                        subtraction += 0.075*bval
                    else:
                        subtraction += 0.15*bval
                    subtracted.append(e.GetIdx())
        
    if m.GetAtomWithIdx(index).GetSymbol() == 'N' and m.GetAtomWithIdx(index).IsInRing():
        subtraction += 1

    print("searched:")
    print(searched)
    print("Setting prop for index " + str(index) + " as " + str(round(parent - subtraction, 4)))
    print(m.GetAtomWithIdx(index).SetProp("atomNote", "\nEst. pKa = " + str(round(parent - subtraction, 4))))
    print("Subtracted indices from " + str(index) + ": " + str(subtracted))
    return round(parent - subtraction, 4)

#SMILES INPUT
#molecule_SMILES = input("Enter a SMILES string: ")
#molecule_SMILES = "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C" #penicillin [3.0] BF CORRECT
#molecule_SMILES = "CC(C1=CC2=C(C=C1)C=C(C=C2)OC)C(=O)O" #naproxen [4.35] BF CORRECT
#molecule_SMILES = "C(CC(=O)O)C(C(=O)O)N" #glutamic acid [4,  3.5, 8] BF CORRECT
#molecule_SMILES = "CC(C)CC(CC(=O)O)CN" #pregabalin [4.4, 9.8] BF CORRECT
#molecule_SMILES = "c1cc(ccc1)C2CCNCC2" #4-phenylpiperidine [9.2] BF CORRECT
#molecule_SMILES =   "C1CC(N(C1)C(=O)C(CCCCN)NC(CCC2=CC=CC=C2)C(=O)O)C(=O)O" #lisinopril [3.9, 3.7, 10.6, 6.6] BF CORRECT
##pubchem requests - remember to print out actual p ka


# molecule_NAME = input("Enter a drug name: ")

# properties = 'Title,CanonicalSMILES,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount'
# titleindex = properties.split(",").index("Title") - 1 #subtract 1 because CID always shows
# smilesindex = properties.split(",").index("CanonicalSMILES") + 1
# url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + molecule_NAME + "/property/" + properties + "/csv"
# #url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/" + molecule_SMILES + "/property/" + properties + "/csv"
# try:
#     u_a = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0) Gecko/20100101 Firefox/39.0"
#     response = requests.get(url, headers={"USER-AGENT":u_a, "Accept": "application/xml"})
#     response = requests.get(url)
#     response.raise_for_status()
#     title = response.text.splitlines()[0].split(",")[titleindex]
#     smiles = response.text.splitlines()[1].split(",")[smilesindex].strip('"')
#     print(smiles)
#     m = Chem.MolFromSmiles(smiles)
#     EWGindices = getEWGs(m)
    
#     COOHindices = []
#     for x in m.GetSubstructMatches(Chem.MolFromSmarts(acids["COOH_al"])):
#         COOHindices.append(x[1])
#     for i in COOHindices:
#         print(BFestimate(i, EWGindices, 4.8))
    
#     amideNindices = []
#     for x in m.GetSubstructMatches(Chem.MolFromSmarts(groups["amide"])):
#         amideNindices.append(x[2])
#     print("amide Ns")
#     print(amideNindices)
#     Nindices = []
#     for x in m.GetSubstructMatches(Chem.MolFromSmarts('N')):
#         if x[0] not in amideNindices:
#             Nindices.append(x[0])
#     for i in Nindices:
#         print(BFestimate(i, EWGindices, 10.6, True))
#     draw(m)


#     for line in response.text.splitlines():
#         print(line)
    
# except requests.exceptions.HTTPError as errh:
#     print(errh)
# except requests.exceptions.ConnectionError as errc:
#     print(errc)
# except requests.exceptions.Timeout as errt:
#     print(errt)
# except requests.exceptions.RequestException as err:
#     print(err)

#for external

def estimateMolecule(m):
    EWGindices = getEWGs(m)
    
    COOHindices = []
    for x in m.GetSubstructMatches(Chem.MolFromSmarts(acids["COOH_al"])):
        COOHindices.append(x[1])
    for i in COOHindices:
        print(BFestimate(m, i, EWGindices, 4.8))
    
    amideNindices = []
    for x in m.GetSubstructMatches(Chem.MolFromSmarts(groups["amide"])):
        amideNindices.append(x[2])
    print("amide Ns")
    print(amideNindices)
    Nindices = []
    for x in m.GetSubstructMatches(Chem.MolFromSmarts('N')):
        if x[0] not in amideNindices:
            Nindices.append(x[0])
    for i in Nindices:
        print(BFestimate(m, i, EWGindices, 10.6, True))
    return m

def nameToSMILES(name):
    properties = 'Title,CanonicalSMILES,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount'
    titleindex = properties.split(",").index("Title") - 1 #subtract 1 because CID always shows
    smilesindex = properties.split(",").index("CanonicalSMILES") + 1
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + name + "/property/" + properties + "/csv"
    
    if name == "":
        return Chem.MolFromSmiles("")
    try:
        u_a = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0) Gecko/20100101 Firefox/39.0"
        response = requests.get(url, headers={"USER-AGENT":u_a, "Accept": "application/xml"})
        response = requests.get(url)
        response.raise_for_status()
        
        data = response.text.splitlines()
        title = response.text.splitlines()[0].split(",")[titleindex]
        
        smiles = response.text.splitlines()[1].split(",")[smilesindex].strip('"')
        for l in data:
            if name.capitalize() in l:
                smiles = l.split(",")[smilesindex].strip('"')
        print(smiles)
        return smiles

    except requests.exceptions.HTTPError as errh:
        print(errh)
    except requests.exceptions.ConnectionError as errc:
        print(errc)
    except requests.exceptions.Timeout as errt:
        print(errt)
    except requests.exceptions.RequestException as err:
        print(err)


def getProperties(m):
    #MW, LogP, TPSA, HBD, HBA, RotB
    properties = {
        "molw": Chem.Descriptors.ExactMolWt(m),
        "logp": Chem.Crippen.MolLogP(m),
        "tpsa": Chem.MolSurf.TPSA(m),
        "hbd": Chem.Lipinski.NumHDonors(m),
        "hba": Chem.Lipinski.NumHAcceptors(m),
        "rotb": Chem.Lipinski.NumRotatableBonds(m)
    }

    print(properties)
    return properties