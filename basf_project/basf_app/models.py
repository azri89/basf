from django.db import models

# Create your models here.
class MoleculeData(models.Model):
    molecule_chembl_id = models.CharField(db_column='molecule_chembl_id', primary_key=True, max_length=255)
    smiles = models.CharField(db_column='smiles', blank=True, null=True, max_length=255)
    molecule_type = models.CharField(db_column='molecule_type', blank=True, null=True, max_length=255)
    alogp = models.CharField(db_column='alogp', blank=True, null=True, max_length=255)