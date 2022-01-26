import streamlit as st
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Fragments

import estimator as e

history = []

def estanddraw(input):
    mol = Chem.MolFromSmiles("")
    try:
        if(choice == 'SMILES'):
            mol = e.estimatesmiles(input)
        elif(choice == 'Drug Name'):
            mol = e.estimatename(input)
        d = rdMolDraw2D.MolDraw2DSVG(800,800) # or MolDraw2DSVG to get SVGs
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().addAtomIndices = True
        d.drawOptions().fixedBondLength = 100
        d.DrawMolecule(mol)
        d.FinishDrawing()
        im = d.GetDrawingText()

        c.image(im)
    
    except AttributeError:
        st.header("Invalid string, probably")


def update():
    # if(choice == 'SMILES'):
    #     st.session_state.input.placeholder = 'SMILES for the camera, please'
    # elif(choice == 'Drug Name'):
    #     st.session_state.input.placeholder = 'Name, please'
    estanddraw(molecule_input)
    history.append(molecule_input)
    


st.title('pKa Estimation')
molecule_input = st.text_input('Input', placeholder='SMILES for the camera, please', key='input')
choice = st.selectbox('Input Format', options = ['SMILES', 'Drug Name'])
c = st.container()
update()

#add history
#add tags (default as pubchem data)
#add resize canvas