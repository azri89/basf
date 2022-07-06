from tkinter import E
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from rest_framework.decorators import api_view
from rest_framework.response import Response

import requests
import json
import openbabel

from .models import MoleculeData
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.rdMolDescriptors import CalcTPSA

@api_view()
def generic_endpoint(request):
    """
    Generic endpoint to check if the service is running.
    """
    return Response()

@api_view(['POST'])
def retrieve_add_molecule_data(request):
    """
    Receive smiles string and retrieve molecule data from  
    ChemBL Database via its API (https://www.ebi.ac.uk/chembl/api/data/docs).

    Create local database entry consisting of molecule_chembl_id as primary key,
    smiles, molecule_type, and alogp if it does not yet exist.

    Return status as reponse.
    """
    try:
        smiles = request.data['smiles']
        url = f'https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__canonical_smiles__flexmatch={smiles}'
        cheml_response = requests.get(url=url)

        all_molecule_data = json.loads(cheml_response.text)['molecules']
        all_molecule_data_stored = []
        for mol_data in all_molecule_data:
            molecule_chembl_id = mol_data['molecule_chembl_id']
            molecule_type = mol_data['molecule_type']
            alogp = mol_data['molecule_properties']['alogp']

            if len(MoleculeData.objects.filter(molecule_chembl_id=molecule_chembl_id)) == 0:
                molecule_data_entry = MoleculeData()
                molecule_data_entry.molecule_chembl_id = molecule_chembl_id
                molecule_data_entry.smiles = smiles
                molecule_data_entry.molecule_type = molecule_type
                molecule_data_entry.alogp = alogp
                molecule_data_entry.save()
                molecule_data = {'molecule_chembl_id': molecule_chembl_id,
                                'smiles': smiles,
                                'molecule_type': molecule_type,
                                'alogp': alogp}
                all_molecule_data_stored.append(molecule_data)
        
        if len(all_molecule_data_stored) > 0:
            response = {'status': f'{len(all_molecule_data_stored)} new molecule data added to the database.',
                        'molecule_data': all_molecule_data_stored}

        else:
            response = {'status': f'All smiles({smiles}) molecule data already exist in the database.'}

        return Response(response)

    except Exception:
        response = {'status': 'Fail to retrieve and add molecule data entry in the database.'}
        return Response(response)

@api_view(['GET'])
def retrieve_molecule_data(request):
    """
    GET request to retrieve molecule information from local database
    based on the molecule_chembl_id submitted. 
    """
    try:
        molecule_chembl_id = request.query_params['molecule_chembl_id']
        molecule_db_entry = MoleculeData.objects.filter(molecule_chembl_id=molecule_chembl_id)

        if len(molecule_db_entry) > 0:
            smiles = molecule_db_entry[0].smiles
            molecule_type = molecule_db_entry[0].molecule_type
            alogp = molecule_db_entry[0].alogp

            response = {'status': 'Molecule data queried from the database.',
                        'molecule_data': {'molecule_chembl_id': molecule_chembl_id,
                                        'smiles': smiles,
                                        'molecule_type': molecule_type,
                                        'alogp': alogp}}

        else:
            response = {'status': 'Molecule data does not exist in the database.'}

        return Response(response)

    except Exception:
        response = {'status': 'Fail to retrieve molecule data from the database.'}
        return Response(response)

@api_view(['POST'])
def add_molecule_data(request):
    """
    POST request to update local database with molecule information
    submitted by the user.
    """
    try:
        molecule_chembl_id = request.data['molecule_chembl_id']
        molecule_db_entry = MoleculeData.objects.filter(molecule_chembl_id=molecule_chembl_id)

        if len(molecule_db_entry) == 0:
            new_molecule_entry = MoleculeData()
            new_molecule_entry.molecule_chembl_id = molecule_chembl_id
            new_molecule_entry.smiles = request.data['smiles']
            if 'molecule_type' in request.data:
                new_molecule_entry.molecule_type = request.data['molecule_type']
            if 'alogp' in request.data:
                new_molecule_entry.alogp = request.data['alogp']
            new_molecule_entry.save()
            response = {'status': 'New molecule data added to the database.'}
        else:
            response = {'status': 'Molecule data already exists in the database.'}

        return Response(response)

    except Exception:
        response = {'status': 'Fail to add molecule data entry in the database.'}

        return Response(response)

@api_view(['GET'])
def convert_smiles_to_inchi(request):
    """
    GET request to convert SMILES string submitted by the user 
    to InchI string as a response using openbabel module.

    Reference: https://open-babel.readthedocs.io/en/latest/UseTheLibrary/Python.html
    """
    try:
        smiles = request.query_params['smiles']
        molecule = openbabel.OBMol()
        molecule_conversion = openbabel.OBConversion()
        molecule_conversion.SetInAndOutFormats("smi", "inchi")

        if molecule_conversion.ReadString(molecule, smiles):
            inchi = molecule_conversion.WriteString(molecule).strip()
            response = {'status': 'SMILES to InChI conversion successful.',
                        'SMILES': smiles,
                        'InChI' : inchi}
            return Response(response)

        else:
            response = {'status': f'{smiles} is not a valid SMILES.'}
            return Response(response)

    except Exception:
        response = {'status': 'SMILES to InChI conversion fail.'}

        return Response(response)

@api_view(['POST'])
def substructure_smiles_tpsa_hits(request):
    """
    Perform substructure SMILES search on ChemBL, 
    https://www.ebi.ac.uk/chembl/api/data/substructure/{smiles}

    Calculate the Topological Polar Surface Area (TPSA) using RDKit
    Reference: https://www.rdkit.org/docs/index.html

    Store molecule with the top 5 highest TPSA in the database
    """
    try:
        substructure_smiles = request.data['substructure_smiles']
        url = f'https://www.ebi.ac.uk/chembl/api/data/substructure/{substructure_smiles}?format=json&limit=0'
        substructure_smiles_response = requests.get(url=url)
        substructure_smiles_json = json.loads(substructure_smiles_response.text)

        all_molecules = substructure_smiles_json['molecules']
        all_molecule_data = []
        for molecule in all_molecules:
            smiles = molecule['molecule_structures']['canonical_smiles']
            mol = MolFromSmiles(smiles)
            tpsa = CalcTPSA(mol)
            molecule_data = {'molecule_chembl_id': molecule['molecule_chembl_id'],
                            'smiles': smiles,
                            'molecule_type': molecule['molecule_type'],
                            'alogp': molecule['molecule_properties']['alogp'],
                            'tpsa': tpsa}
            all_molecule_data.append(molecule_data)
        top_molecule_data = sorted(all_molecule_data, key=lambda d: d['tpsa'], reverse=True)[:5]

        for molecule in top_molecule_data:
            if len(MoleculeData.objects.filter(molecule_chembl_id=molecule['molecule_chembl_id'])) == 0:
                molecule_data_entry = MoleculeData()
                molecule_data_entry.molecule_chembl_id = molecule['molecule_chembl_id']
                molecule_data_entry.smiles = molecule['smiles']
                molecule_data_entry.molecule_type = molecule['molecule_type']
                molecule_data_entry.alogp = molecule['alogp']
                molecule_data_entry.save()

        response = {'status': 'Top 5 TPSA hits stored in the database.',
                    'Top 5 TPSA': top_molecule_data}
        return Response(response)

    except Exception:
        response = {'status': 'Substructure SMILES search and TPSA calculation fail.'}

        return Response(response)