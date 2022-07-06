from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from rest_framework.decorators import api_view
from rest_framework.response import Response

import requests
import json
import openbabel

from .models import MoleculeData

@api_view()
def generic_endpoint(request):
    """
    Generic endpoint to check if the service is running.
    """
    return Response({'status': 'Service is running'})

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

        molecule_data = json.loads(cheml_response.text)['molecules'][0]
        molecule_chembl_id = molecule_data['molecule_chembl_id']
        molecule_type = molecule_data['molecule_type']
        alogp = molecule_data['molecule_properties']['alogp']

        if len(MoleculeData.objects.filter(molecule_chembl_id=molecule_chembl_id)) == 0:
            molecule_data_entry = MoleculeData()
            molecule_data_entry.molecule_chembl_id = molecule_chembl_id
            molecule_data_entry.smiles = smiles
            molecule_data_entry.molecule_type = molecule_type
            molecule_data_entry.alogp = alogp
            molecule_data_entry.save()
            status_message = 'New molecule data added to the database.'
        else:
            status_message = 'Molecule data already exists in the database.'
            
        response = {'status': status_message,
                    'molecule_data': {'molecule_chembl_id': molecule_chembl_id,
                                    'smiles': smiles,
                                    'molecule_type': molecule_type,
                                    'alogp': alogp}}

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

    except Exception as e:
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
    to InchI string as a response.

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