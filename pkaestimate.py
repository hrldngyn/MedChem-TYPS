from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Fragments

import matplotlib.pyplot as pyplot
import matplotlib.image as mpimg

import requests


from collections.abc import Iterable
#%matplotlib inline

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
    return set(flatten(EWGindices))

def atomArrayToIndex(array):
    indexArray = []
    for a in array:
        indexArray.append(a.GetIdx())
    return indexArray

def draw():    
    d = rdMolDraw2D.MolDraw2DCairo(1000, 1000) # or MolDraw2DSVG to get SVGs

    d.drawOptions().addStereoAnnotation = True
    d.drawOptions().addAtomIndices = True
    d.DrawMolecule(m)
    d.FinishDrawing()
    d.WriteDrawingText('atom_annotation_1.png')   

def getAcidEstimate(COOHindices):
    for n in COOHindices: #for each carbon of a COOH
        oneCAway = []
        twoCAway = []
        threeCAway = []
        fourCAway = []
        searched = []
        neighbors = m.GetAtomWithIdx(n).GetNeighbors() #find neighbors, returns list of Atoms
        for p in neighbors:
            searched.append(p.GetIdx())   
            if(p.GetIdx() not in m.GetSubstructMatches(Chem.MolFromSmarts(groups["COOH"]))): #check if neighbor is not part of the COOH
                for q in p.GetNeighbors():
                    searched.append(q.GetIdx()) 
                    if(q.GetIdx() not in COOHindices):
                            oneCAway.append(q)
        for q in oneCAway:
            for r in q.GetNeighbors():
                if(r.GetIdx() not in searched): #if r is not p
                    twoCAway.append(r)
                searched.append(r.GetIdx())
        for r in twoCAway:
            for s in r.GetNeighbors():
                #print(s.GetIdx())
                if(s.GetIdx() not in searched): #if s is not q
                    threeCAway.append(s)
                searched.append(s.GetIdx())
        #print(atomArrayToIndex(threeCAway))
        for s in threeCAway:
            for t in s.GetNeighbors():
                if(t.GetIdx() not in searched): #if t is not 
                    fourCAway.append(t)
                searched.append(t.GetIdx())
        deduction = getCOOHDeduction(n, oneCAway, twoCAway, threeCAway, fourCAway)
        print(str(n) + ": " + str(round(4.8-deduction,3)) + " Acid")
        m.GetAtomWithIdx(n).SetProp("atomNote", "\nEst. pKa = " + str(round(4.8-deduction,3)) + ", Acid")     
       
def getCOOHDeduction(n, oneCAway, twoCAway, threeCAway, fourCAway):
    total = 0
    for q in oneCAway:
        if q.GetIdx() in EWGindices: #if q is an ewg
            if q.GetIsAromatic(): #if q aromatic add 0.45
                total += 0.45
            else: #else q not aromatic add full 0.9
                total += 0.9
    if total == 0:
        for r in twoCAway:
            if r.GetIdx() in EWGindices:
                if r.GetIsAromatic(): #if q aromatic add 0.45
                    total += 0.45
                else: #else q not aromatic add full 0.9
                    total += 0.9
    if total == 0:
        for s in threeCAway:
            if s.GetIdx() in EWGindices:
                if s.GetIsAromatic(): #if q aromatic add 0.45
                    total += 0.2
                else: #else q not aromatic add full 0.9
                    total += 0.4
    if total == 0:
        for t in fourCAway:
            if t.GetIdx() in EWGindices:
                if t.GetIsAromatic(): #if q aromatic add 0.45
                    total += (0.15/2)
                else: #else q not aromatic add full 0.9
                    total += 0.15
    return (total)

def getAmineEstimate(amineindices):
    deduction = 0
    
    for n in amineindices: #for each amine
        deductedindices = [] #keep track of duplicates like in cyclic compounds
        neighbors = m.GetAtomWithIdx(n).GetNeighbors() #find neighbors, returns list of Atoms
        for p in neighbors: # for each leg of the amine
            oneCAway = []
            twoCAway = []
            threeCAway = []
            fourCAway = []
            searched = []
            searched.append(p.GetIdx())   
            #if(p.GetIdx() not in m.GetSubstructMatches(Chem.MolFromSmarts("[[N;H2;X3][CX4],[N;H;X3]([CX4])[CX4],[NX3]([CX4])([CX4])[CX4]]"))): #check if neighbor is not part of the amines
            for q in p.GetNeighbors():
                searched.append(q.GetIdx()) 
                if(q.GetIdx() not in amineindices):
                        oneCAway.append(q)
            for q in oneCAway:
                for r in q.GetNeighbors():
                    if(r.GetIdx() not in searched): #if r is not p
                        twoCAway.append(r)
                    searched.append(r.GetIdx())
            for r in twoCAway:
                for s in r.GetNeighbors():
                    #print(s.GetIdx())
                    if(s.GetIdx() not in searched): #if s is not q
                        threeCAway.append(s)
                    searched.append(s.GetIdx())
            #print(atomArrayToIndex(threeCAway))
            for s in threeCAway:
                for t in s.GetNeighbors():
                    if(t.GetIdx() not in searched): #if t is not s
                        fourCAway.append(t)
                    searched.append(t.GetIdx())
            deduction += getAmineDeduction(oneCAway, twoCAway, threeCAway, fourCAway, deductedindices)
        if(n in cyclicamineindices):
            deduction += 1
        print(str(n) + ": " + str(round(10.6-deduction,2)) + " Amine" + str(deductedindices))
        m.GetAtomWithIdx(n).SetProp("atomNote", "\nEst. pKa = " + str(round(10.6-deduction,2)) + ", Amine")

    #return n, parent pka, est pka

def getAmineDeduction(oneCAway, twoCAway, threeCAway, fourCAway, deductedindices):
    total = 0
    dedind = []
    for q in oneCAway:
        if (q.GetIdx() in EWGindices): #if q is an ewg
            if (q.GetIdx() not in deductedindices):
                if q.GetIsAromatic(): #if q aromatic add 0.45
                    total += 0.9
                else: #else q not aromatic add full 0.9
                    total += 1.8
                dedind.append(q.GetIdx())
            else:
                for n in dedind:
                    deductedindices.append(n)
                return total
    if len(dedind) == 0:
        for r in twoCAway:
            if (r.GetIdx() in EWGindices): #if r is an ewg
                if (r.GetIdx() not in deductedindices):
                    if r.GetIsAromatic(): #if q aromatic add 0.45
                        total += 0.9
                    else: #else q not aromatic add full 0.9
                        total += 1.8
                    dedind.append(r.GetIdx())
                else:
                    for n in dedind:
                        deductedindices.append(n)
                    return total
    if len(dedind) == 0:
        for s in threeCAway:
            if (s.GetIdx() in EWGindices): #if s is an ewg
                if (s.GetIdx() not in deductedindices):
                    if s.GetIsAromatic(): #if q aromatic add 0.4
                        total += 0.4
                    else: #else q not aromatic add full 0.8
                        total += 0.8
                    dedind.append(s.GetIdx())
                else:
                    for n in dedind:
                        deductedindices.append(n)
                    return total
    if len(dedind) == 0:
        for t in fourCAway:
            if (t.GetIdx() in EWGindices): #if t is an ewg
                if (t.GetIdx() not in deductedindices):
                    if t.GetIsAromatic(): #if t aromatic add 0.15
                        total += 0.15
                    else: #else t not aromatic add full 0.3
                        total += 0.3
                    dedind.append(t.GetIdx())
                else:
                    for n in dedind:
                        deductedindices.append(n)
                    return total
    for n in dedind:
        deductedindices.append(n)
    return total
groups = {"amine": "[#7]",
          "oxygen": "[#8]",
          "COOH": "[CX3](=O)[OX2H1]",
          "carbonyl": "[CX3]=[OX1]" ,
          "aromatic C": "[$([cX3](:*):*),$([cX2+](:*):*)]",
          "sulfur": "[S]",
          "phosphorous": "[P]",
          "halogen": "[F,Cl,Br,I]",
    }

#SMILES INPUT
    #molecule_SMILES = input("Enter a SMILES string: ")
    #molecule_SMILES = "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C" #penicillin
    #molecule_SMILES = "CC(C1=CC2=C(C=C1)C=C(C=C2)OC)C(=O)O" #naproxen  
    #molecule_SMILES = "C(CC(=O)O)C(C(=O)O)N" #glutamic acid
#molecule_SMILES = "CC(C)CC(CC(=O)O)CN" #pregabalin
molecule_SMILES = "c1cc(ccc1)C2CCNCC2" #4-phenylpiperidine

##pubchem requests - remember to print out actual p ka
drugsmiles = molecule_SMILES #input("Enter a SMILES string")
url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/" + drugsmiles + "/property/Title,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount/csv" #,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount/"
try:
    u_a = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0) Gecko/20100101 Firefox/39.0"
    response = requests.get(url, headers={"USER-AGENT":u_a, "Accept": "application/xml"})
    response = requests.get(url)
    response.raise_for_status()
    for line in response.text.splitlines():
        print(line)
except requests.exceptions.HTTPError as errh:
    print(errh)
except requests.exceptions.ConnectionError as errc:
    print(errc)
except requests.exceptions.Timeout as errt:
    print(errt)
except requests.exceptions.RequestException as err:
    print(err)


#convert SMILES to rdkit molecule
m = Chem.MolFromSmiles(molecule_SMILES)

EWGindices = getEWGs(m)

#COOH estimates
numCOOH = Fragments.fr_Al_COO(m)
COOHindices = []
for x in m.GetSubstructMatches(Chem.MolFromSmarts(groups["COOH"])):
    COOHindices.append(x[0])

#amine estimates
#[NX3;H2,H1;!$(NC=O)] = primary or secondary amine, no amide
#"[[N;H2;X3][CX4],[N;H;X3]([CX4])[CX4],[NX3]([CX4])([CX4])[CX4]]" all amines
amineindices = []
cyclicamineindices = []
aromaticAmineindices = [] #[N;R]
#for x in m.GetSubstructMatches(Chem.MolFromSmarts("[[N;H2;X3][CX4],[N;H;X3]([CX4])[CX4],[NX3]([CX4])([CX4])[CX4]]")):
for x in m.GetSubstructMatches(Chem.MolFromSmarts("[$([NX4+]),$([NX3]);!$(*=*)&!$(*:*)]")):
    amineindices.append(x[0])
for x in m.GetSubstructMatches(Chem.MolFromSmarts("[N&R]")):
    cyclicamineindices.append(x[0])
#print pka estimates
getAcidEstimate(COOHindices)
getAmineEstimate(amineindices)
                     
#Info Output
draw()
print("Number of COOHs " + str(numCOOH))
print("Indices of COOH carbons " + str(COOHindices))

print("Number of Amines " + str(len(amineindices)))
print("Indices of Amine Nitrogens" + str(amineindices))
print("Number of cyclic Amines " + str(len(cyclicamineindices)))
print("Indices of cyclic Amine Nitrogens" + str(cyclicamineindices))


print("EWG indices" + str(getEWGs(m)))