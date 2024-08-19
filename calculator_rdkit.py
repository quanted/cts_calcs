from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

from .calculator import Calculator



class RdkitCalc(Calculator):

    def __init__(self):
        # Establish smarts objects, NOTE these strings will not change!
        self.CAE_smarts = '[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]'
        self.CAE = Chem.MolFromSmarts(self.CAE_smarts)
        self.CarbAnhydride_smarts = '[CX3;$([H0][#6]),$([H1])](=[OX1])[#8X2][CX3;$([H0][#6]),$([H1])](=[OX1])'
        self.CarbAnhydride = Chem.MolFromSmarts(self.CarbAnhydride_smarts)

        self.meta_info = {
            'metaInfo': {
                'model': "rdkit",
                'collection': "qed",
                'modelVersion': "2024.3.5",
                'description': "A collection of chemoinformatics and machine-learning software written in C++ and Python.",
                'status': '',
                'timestamp': self.gen_jid(),
                'url': "https://www.rdkit.org/",
                'props': []
            }
        }

    def increment_atom_number(self, atom_list):
        atoms = []
        for site in atom_list:
            new=[i+1 for i in site]
            atoms.append(new)
        logging.warning("Incremented atoms: {}".format(atoms))
        return(atoms)

    def get_functional_groups_anhydride(self, smiles):
        """
        List of atoms in anhydride functional group.
        """
        mol = Chem.MolFromSmiles(smiles)
        anhydride_atom = list(Chem.Mol.GetSubstructMatches(mol, self.CarbAnhydride, uniquify=True))
        logging.info("Anhydride group: {}".format(anhydride_atom))
        # anhydride_atom = self.increment_atom_number(anhydride_atom)
        # logging.info("Updated Anhydride group: {}".format(anhydride_atom))
        return anhydride_atom

    def get_functional_groups_cae(self, smiles):
        """
        List of atoms in carboxylic acid ester functional group.
        """
        mol = Chem.MolFromSmiles(smiles)
        CAE_atom = list(Chem.Mol.GetSubstructMatches(mol, self.CAE, uniquify=True))
        logging.info("CAE group: {}".format(CAE_atom))
        # CAE_atom = self.increment_atom_number(CAE_atom)
        # logging.info("Updated CAE group: {}".format(CAE_atom))
        return CAE_atom

    def get_functional_groups(self, route, smiles):
        """
        Calls rdkit to get functional groups.
        """
        func_group = []
        if "anhydride" in route.lower():
            func_group = self.get_functional_groups_anhydride(smiles)
        elif route.lower() == "carboxylic acid ester hydrolysis":
            func_group = self.get_functional_groups_cae(smiles)

        # increments atom numbers:
        func_group = self.increment_atom_number(func_group)

        logging.warning("Incremented groups: {}".format(func_group))

        return func_group

    def get_diffusivity(self, request_dict):
        """
        Returns diffusivity in air and water.
        Diffusivity in Water returns W-C (Wilke-Chang equation) and H-L (Hayduk and Luadie equation) methods.
        Diffusivity in Air returns FSG (Fuller, Schettler, and Giddings equation) method.
        """

        smiles = request_dict.get("chemical")

        #get molecular weight
        mw=Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles(smiles))
        
        #get molar mass for air
        a_mass=(1.55/mw**0.65)
        #get molar mass for water
        w_mass=(2.7E-4)/(mw**0.71)
            
        #get van der waal volume (molecular volume)    
        mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol)
        vol=Chem.AllChem.ComputeMolVolume(mol)


        #constants
        n=8.90E-04 #dynamic viscosity of water, 25C
        T=298 #temp in K
        k=1.38E-23 #boltzman constant in kg-m2/s2-K
        pi=3.1415926
        X=2.6 #constant that depends on solvent, this is the constant for water
        
         #calculate FSG (air) diffusivity coefficinent
        FSG=10**-3*((T**1.75*((1/28.97)+(1/mw))**0.5)/(20.1**(1/3)+(vol**(1/3)))**2)
        
        #calculate Wilke-Chang (water) diffusivity coefficinent
        WC=(7.4E-8)*((((X*w_mass)**0.5)*T)/(n*((vol)**0.6)))
        
         #calculate Hayduk-Laudie (water) diffusivity coefficinent
        HL=13.26E-5/((n**1.14)*((vol)**0.589))

        diff_vals = {
            'FSG': FSG,
            'W-C': WC,
            'H-L': HL
        }

        response_obj = dict(request_dict)
        response_obj["data"] = diff_vals
        
        return response_obj