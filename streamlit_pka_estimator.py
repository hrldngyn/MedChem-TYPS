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

history = []
oldprops = {}

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

def estanddraw(c, mol):
    try:
        print(mol)
        e.estimateMolecule(mol)

        if view == "RDKit 2D":
            d = rdMolDraw2D.MolDraw2DSVG(800,800) # or MolDraw2DSVG to get SVGs
            d.drawOptions().addStereoAnnotation = True
            d.drawOptions().addAtomIndices = True
            d.drawOptions().fixedBondLength = 100
            d.DrawMolecule(mol)
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
def displayproperties(props, deltas):
    col1, col2, col3, col4, col5, col6 = st.columns(6)
    col1.metric(label="Molecular Weight", value=props["molw"], delta = deltas["molw"])
    col2.metric(label="Crippen LogP", value=props["logp"], delta = deltas["logp"])
    col3.metric(label="TPSA", value=props["tpsa"], delta = deltas["tpsa"])
    col4.metric(label="HB Donors", value=props["hbd"], delta = deltas["hbd"])
    col5.metric(label="HB Acceptors", value=props["hba"], delta = deltas["hba"])
    col6.metric(label="Rotatable Bonds", value=props["rotb"], delta = deltas["rotb"])


def pkamodule():
    def updatepka():
        try:
            mol = strToMol(molecule_input, inputmode)
            if mol.GetNumAtoms() != 0:
                estanddraw(c, mol)

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
            displayproperties(props, deltas)

            history.append(molecule_input)
        
        except AttributeError:
            st.error("Not a valid molecule")

    ph = {
        'SMILES': 'SMILES for the camera, please',
        'Drug Name': 'Enter a drug name to search PubChem'
    }

    st.title('pKa Estimation')
    inputmode = st.selectbox('Input Format', options = ['SMILES', 'Drug Name'])
    molecule_input = st.text_input('Input', placeholder=ph[inputmode], key='input')

    c = st.container()
    if molecule_input != "":
        updatepka()
    
modecol, viewcol = st.sidebar.columns(2)
mode = modecol.radio(
    label = "Mode", 
    options = ["pKa Estimation"])
view = viewcol.radio(
    label = "View",
    options = ["RDKit 2D", "PyMol3D"]
)
if mode == "pKa Estimation":
    pkamodule()
#add history
#add tags (default as pubchem data)
#add resize canvas

#target mode: log p and pka


#log p visualizer
