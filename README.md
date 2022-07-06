# Prerequisites
You need to have docker (https://docs.docker.com/) install on your machine to run the app.

# Installation
Clone the repository in your local directory.

Run `docker compose build` to start building the docker containers.

# Starting the app
a) Run `docker compose up` in your terminal from the app directory.

b) Run `docker compose up -d` to run the app in detached mode.

# Available Endpoints
**1.** Generic endpoint just returning a heartbeat.
* **URL:**
 /basf_app/
 
 * **Method:**
  `GET`
  
  * **URL Params:**
  None
  
  * **Data Params:**
  None
  
  * **Success Response:**
    **Code:** 200 <br />



**2.** Provide an endpoint to retrieve data from the ChEMBL Database via its API. The user sends in
a molecule in a string format (called SMILES), the API then creates a local database entry and returns a
status for the database entry generation.
* **URL:**
    /basf_app/retrieve_add_molecule_data/
 
 * **Method:**
      `POST`
  
  * **URL Params:**
      None
  
  * **Data Params:**
      ```
      {"smiles": "CCCC"}
      ```
  
  * **Success Response:**
    **Content:** 
    ```
    { "status": "1 new molecule added to the database.", 
      "molecule_data": [{ "molecule_chembl_id": "CHEMBL123",
                         "smiles": "CCCC",
                         "molecule_type": "Small molecule",
                         "alogp": "1.00"}, ....]}`
    ```
                                       
 
 
**3.** An endpoint to check and retrieve data from the local database via the molecule_chembl_id (JSON format).
  * **URL:**
   /basf_app/retrieve_molecule_data/?molecule_chembl_id=<molecule_chembl_id>

  * **Method:**
    `GET`

  * **URL Params:**
  molecule_chembl_id

  * **Data Params:**
  None

  * **Success Response:**
    **Content:** 
    ```
    { "status": "Molecule data queried from the database.", 
      "molecule_data": { "molecule_chembl_id": "CHEMBL123",
                         "smiles": "CCCC",
                         "molecule_type": "Small molecule",
                         "alogp": "1.00"}}`
    ```


**4.** An endpoint to create a new record in the local database.
* **URL:**
    /basf_app/add_molecule_data/
 
 * **Method:**
      `POST`
  
  * **URL Params:**
      None
  
  * **Data Params:**
      ```
      { "molecule_chembl_id": "CHEMBL123",
         "smiles": "CCCC",
         "molecule_type": "Small molecule",
         "alogp": "1.00"}
      ```
  
  * **Success Response:**
    **Content:** 
    ```
    {"status": "New molecule data added to the database."}
    ```
    

**5.** An endpoint accepting the SMILES of a molecule as input, returning the InChi (another molecular string representation).
* **URL:**
 /basf_app/convert_smiles_to_inchi/?smiles=<smiles>

* **Method:**
  `GET`

* **URL Params:**
smiles

* **Data Params:**
None

* **Success Response:**
  **Content:** 
    ```
    {"status": "SMILES to InChI conversion successful.",
    "SMILES": "CC(=O)OCC",
    "InChI": "InChI=1S/C4H8O2/c1-3-6-4(2)5/h3H2,1-2H3"}
    ```
                
                
                
**6.** An endpoint accepting a substructure SMILES, performing a substructure search on ChEMBL,
calculating the topological polar surface area (TPSA) of all hits via RDKit and storing the 5 hits with the
highest TPSA in the local database.
* **URL:**
    /basf_app/substructure_smiles_tpsa_hits/
 
 * **Method:**
      `POST`
  
  * **URL Params:**
      None
  
  * **Data Params:**
    ```
    { "substructure_smiles": "Cn(c(=O)c1n(C)c(c2ccccc2)nc1n3C)c3=O" }
    ```
  
  * **Success Response:**
    **Content:** 
      ```
      { "status": "Top 5 TPSA hits stored in the database.",
        "Top 5 TPSA": [{ "molecule_chembl_id": "CHEMBL3359729",
                          "smiles": "CN(CCOc1ccc(-c2nc3c(c(=O)n(CC(=O)O)c(=O)n3C)n2CC(=O)O)cc1)c1ccccn1",
                          "molecule_type": "Small molecule",
                          "alogp": "0.64",
                          "tpsa": 161.78000000000003}, ....]}
      ```
    
                
