from django.urls import path
from . import views

urlpatterns = [   
    # api
    path('', views.generic_endpoint),
    path('retrieve_add_molecule_data/', views.retrieve_add_molecule_data),
    path('retrieve_molecule_data/', views.retrieve_molecule_data),
    path('add_molecule_data/', views.add_molecule_data),
    path('convert_smiles_to_inchi/', views.convert_smiles_to_inchi),
    path('substructure_smiles_tpsa_hits/', views.substructure_smiles_tpsa_hits)
]