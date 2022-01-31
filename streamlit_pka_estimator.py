import streamlit as st
from stmol import showmol
import py3Dmol

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Fragments

import pickle
import os

import estimator as e


def strToMol(input, inputmode):
    if(inputmode == 'SMILES'):
        mol = Chem.MolFromSmiles(input)
        if mol is None:
            mol = Chem.MolFromSmiles("")
            st.error("Not a valid SMILES, probably")
    elif(inputmode == 'Drug Name'):
        try:
            mol = Chem.MolFromSmiles(e.nameToSMILES(input))
        except TypeError:
            mol = Chem.MolFromSmiles("")
            st.error("PubChem search failed, probably")
    return mol

def draw(c, mol):
    try:

        

        if view == "RDKit 2D":
            d = rdMolDraw2D.MolDraw2DSVG(800,800) # or MolDraw2DSVG to get SVGs
            d.drawOptions().addStereoAnnotation = True
            d.drawOptions().addAtomIndices = True
            d.drawOptions().fixedBondLength = 100
            d.drawOptions().comicMode = True
            d.DrawMolecule(mol)
            d.TagAtoms(mol)
            d.FinishDrawing()
            im = d.GetDrawingText()
            c.image(im)
        elif view == "PyMol3D":
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            mblock = Chem.MolToMolBlock(mol)
            v = py3Dmol.view()#(width=400,height=400)
            v.addModel(mblock,'mol')
            v.setStyle({'stick':{}})
            v.setBackgroundColor('white')
            v.zoomTo()
            showmol(v,height=600,width=800)

        #c.write(e.getProperties(mol))
    
    except AttributeError:
        st.error("Invalid string, probably")
    

def getdeltas(old, props):
    deltas = {}
    if old['molw'] > 0: #if oldprops existed
        for label in props:
            deltas[label] = props[label] - old[label]
    else:
        for label in props:        
            deltas[label] = None
    return deltas
def displayproperties(props, deltas, showanswers):
    if showanswers:
        col1, col2, col3, col4, col5, col6 = st.columns(6)
        col1.metric(label="Molecular Weight", value=props["molw"], delta = deltas["molw"])
        col2.metric(label="Crippen LogP", value=props["logp"], delta = deltas["logp"])
        col3.metric(label="TPSA", value=props["tpsa"], delta = deltas["tpsa"])
        col4.metric(label="HB Donors", value=props["hbd"], delta = deltas["hbd"])
        col5.metric(label="HB Acceptors", value=props["hba"], delta = deltas["hba"])
        col6.metric(label="Rotatable Bonds", value=props["rotb"], delta = deltas["rotb"])

        col7, col8, col9, col10, col11, col12 = st.columns(6)
        col7.metric(label="Carbonyls", value = props["carbonyls"], delta = deltas["carbonyls"])


def pkamodule():
    def updatepka():
        try:
            mol = strToMol(molecule_input, inputmode)
            if mol.GetNumAtoms() != 0:
                if showanswers:
                    e.estimateMolecule(mol)
                    draw(c, mol)
                else:
                    draw(c, mol)

            props = e.getProperties(mol)
            oldprops = {'molw': 0}
            filename = 'oldprops.pk'
            if os.path.exists(filename) and os.stat(filename).st_size != 0:
                with open(filename, 'rb') as fi:
                    oldprops = pickle.load(fi)
                    print(oldprops)
            if props['molw'] != 0:
                with open(filename, 'wb') as fi:
                    pickle.dump(props, fi)
            
            deltas = getdeltas(oldprops, props)
            if deltas['molw'] == 0 or props['molw'] == 0:
                for label in deltas:        
                    deltas[label] = None
            
            name = e.SMILESToName(Chem.MolToSmiles(mol))
            st.header(name)
            displayproperties(props, deltas, showanswers)

            history.append((mol, molecule_input))
        
        except AttributeError:
            st.error("Not a valid molecule")

    ph = {
        'SMILES': 'SMILES for the camera, please',
        'Drug Name': 'Enter a drug name to search PubChem'
    }

    st.title('pKa Estimation')
    inputmode = st.selectbox('Input Format', options = ['SMILES', 'Drug Name'])
    molecule_input = st.text_input('Input', placeholder=ph[inputmode], key='input')
    showanswers = st.checkbox("Show Answers")

    c = st.container()

    history = []

    if molecule_input != "":
        updatepka()






def substratecrystalsmodule():
    sclass = optscol.radio(
        label = "Substrate Class",
        options = ["SSRIs", "NSAIDs", "custom"])
    style = st.sidebar.selectbox('style',['line','cross','stick','sphere','cartoon','clicksphere'])
    spin = st.sidebar.checkbox('Spin', value = False)
    if sclass != "custom":
        protviews = loadproteins(drugclass = sclass)
        protpicks = []
        for v in protviews:
            protpicks.append(v[1])
    

        pick = st.sidebar.selectbox('Select Protein', protpicks)
        activeview = protviews[protpicks.index(pick)][0]
        activeview.setStyle({style:{'color':'spectrum'}})
        activeview.setBackgroundColor('white')
        activeview.zoomTo()
        activeview.spin(spin)
        showmol(protviews[protpicks.index(pick)][0], height=2000, width=2000)

    elif sclass == "custom":
        uploaded_files = st.sidebar.file_uploader("Choose .pdb files", accept_multiple_files=True)
        for f in uploaded_files:
            pdb = f.getvalue().decode("utf-8")
            pdbview = py3Dmol.view(width=1000, height=1000)
            pdbview.addModel(pdb, 'pdb')
            pdbview.setStyle({style:{'color':'spectrum'}})
            pdbview.setBackgroundColor('white')
            pdbview.zoomTo()
            pdbview.spin(spin)
            showmol(pdbview, height=2000, width=2000)
            


def loadproteins(pdblist = None, drugclass = None):

    pdbviewlist = []
    drugclasses = {
        "SSRIs": '5I71,5I6X,6AWO',
        "NSAIDs": ''
    }

    drugs = {
        "5I71": "escitalopram central site", 
        "5I6X": "paroxetine central site",
        "6AWO": "sertraline central site" 
    }

    for p in drugclasses[drugclass].split(","):
        pdbviewlist.append((py3Dmol.view(query='pdb:'+ p), drugs[p]))
    
    return pdbviewlist






modecol, optscol = st.sidebar.columns(2)
mode = modecol.radio(
    label = "Mode", 
    options = ["pKa Estimation", "Protein Exploration"])

if mode == "pKa Estimation":
    view = optscol.radio(
    label = "View",
    options = ["RDKit 2D", "PyMol3D"]
    )

    pkamodule()

elif mode == "Protein Exploration":
    
    substratecrystalsmodule()
#add history
#add tags (default as pubchem data)
#add resize canvas

#target mode: log p and pka


#log p visualizer
