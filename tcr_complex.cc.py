import math
import os

class TCRComplex:
    PI = math.pi

    def __init__(self, tcrfile, pepmhcfile, mhc_type):
        # Initialization
        self.my_tcra_chain = "D"
        self.my_tcrb_chain = "E"
        self.my_mhc_chain = "A"
        self.my_pep_chain = "C"
        
        self.my_mhc_type = mhc_type
        self.IsClassII = self.my_mhc_type == 1

        self.my_z_offset = 25.0
        self.my_tcrz_shift = 0.0

        # Alternative TCR docking start site (start2)
        self.my_tcrx_shift = 20.0
        self.my_tcr_tilt_shift = 25.0
        self.my_tcr_cross_shift = -45.0

        # Uncomment for start3
        # self.my_tcrx_shift = 10.0
        # self.my_tcr_tilt_shift = 12.5
        # self.my_tcr_cross_shift = -22.5

        # Uncomment for start4
        # self.my_tcrx_shift = 25.0
        # self.my_tcr_tilt_shift = 40.0
        # self.my_tcr_cross_shift = -60.0

        if not pepmhcfile:
            pepmhcfile = tcrfile

        if self.IsClassII:
            self.my_mhcb_chain = "B"

        self.my_tcr_file = tcrfile
        self.my_pepmhc_file = pepmhcfile

        self.verbose = True

        if self.my_mhc_type == 6:  # Antibody
            self.my_tcra_chain = "L"
            self.my_tcrb_chain = "H"

        # Load TCR and PepMHC data
        self.my_tcra = None  # Placeholder for TCR Alpha chain data
        self.my_tcrb = None  # Placeholder for TCR Beta chain data
        self.my_mhc = None  # Placeholder for MHC data
        self.my_pep = None  # Placeholder for peptide data

        self.load_tcr(self.my_tcr_file)
        if self.my_mhc_type != 6:
            self.load_pepmhc(self.my_pepmhc_file)

    def __del__(self):
        # Destructor to free resources
        if hasattr(self, 'my_tcra') and self.my_tcra and hasattr(self.my_tcra, 'Atoms'):
            del self.my_tcra.Atoms
        if hasattr(self, 'my_tcrb') and self.my_tcrb and hasattr(self.my_tcrb, 'Atoms'):
            del self.my_tcrb.Atoms
        if self.my_mhc_type < 6:
            if hasattr(self, 'my_mhc') and self.my_mhc and hasattr(self.my_mhc, 'Atoms'):
                del self.my_mhc.Atoms
            if hasattr(self, 'my_pep') and self.my_pep and hasattr(self.my_pep, 'Atoms'):
                del self.my_pep.Atoms
            if self.IsClassII and hasattr(self, 'my_mhcb') and self.my_mhcb and hasattr(self.my_mhcb, 'Atoms'):
                del self.my_mhcb.Atoms

    def load_tcr(self, tcrfile):
        # Placeholder for TCR loading logic
        print(f"Loading TCR data from {tcrfile}...")

    def load_pepmhc(self, pepmhcfile):
        # Placeholder for PepMHC loading logic
        print(f"Loading PepMHC data from {pepmhcfile}...")

class Atom:
        def __init__(self, x=0.0, y=0.0, z=0.0, res=0, ins_code="", atom="", res_name="", line=""):
            self.x = x
            self.y = y
            self.z = z
            self.res = res
            self.ins_code = ins_code
            self.atom = atom
            self.res_name = res_name
            self.line = line

class ProteinChain:
        def __init__(self):
            self.Atoms = []
            self.chain_id = ""
            self.num_atoms = 0

def __init__(self):
        self.IsClassII = False

        # Initialize placeholders for chains
        self.my_mhc = self.ProteinChain()
        self.my_pep = self.ProteinChain()
        self.my_mhcb = self.ProteinChain()
        self.my_tcra = self.ProteinChain()
        self.my_tcrb = self.ProteinChain()

def load_pepmhc(self, filename):
        self.load_chain_from_file(filename, self.my_mhc_chain, self.my_mhc)
        self.load_chain_from_file(filename, self.my_pep_chain, self.my_pep)

        if self.my_mhc.num_atoms == 0:
            raise ValueError("No atoms loaded for MHC!")

        if self.IsClassII:
            self.load_chain_from_file(filename, self.my_mhcb_chain, self.my_mhcb)
            if self.my_mhcb.num_atoms == 0:
                raise ValueError("No atoms loaded for MHC B chain!")

def load_tcr(self, filename):
        self.load_chain_from_file(filename, self.my_tcra_chain, self.my_tcra)
        self.load_chain_from_file(filename, self.my_tcrb_chain, self.my_tcrb)

        if self.my_tcra.num_atoms == 0 or self.my_tcrb.num_atoms == 0:
            raise ValueError("No atoms loaded for TCR alpha or beta!")

def load_chain_from_file(self, filename, chain_id, prot):
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"Error opening PDB file: {filename}")

        # First pass to count atoms
        with open(filename, "r") as pdbfile:
            num_atoms = sum(1 for line in pdbfile if line.startswith(("ATOM", "HETATM")) and line[21] == chain_id)

        prot.Atoms = [self.Atom() for _ in range(num_atoms)]
        prot.chain_id = chain_id
        prot.num_atoms = 0

        # Second pass to load atoms
        with open(filename, "r") as pdbfile:
            for line in pdbfile:
                if not line.startswith(("ATOM", "HETATM")):
                    continue

                ch_id = line[21]
                if ch_id == chain_id:
                    xcoord = float(line[30:38].strip())
                    ycoord = float(line[38:46].strip())
                    zcoord = float(line[46:54].strip())
                    res = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    ins_code = line[26:27].strip()
                    atom = line[12:16].strip()
                    res_name = line[17:20].strip()

                    if self.is_already_loaded(prot, res_num, ins_code, atom):
                        print(f"Warning: duplicate entries found in chain {ch_id} "
                              f"res {res_num}{ins_code} atom {atom}")

                    prot.Atoms[prot.num_atoms] = self.Atom(
                        x=xcoord, y=ycoord, z=zcoord, res=res_num,
                        ins_code=ins_code, atom=atom, res_name=res_name, line=line.strip()
                    )
                    prot.num_atoms += 1

def is_already_loaded(self, prot, res_num, ins_code, atom):
        # Check if an atom with the same properties is already loaded
        for existing_atom in prot.Atoms:
            if (existing_atom.res == res_num and
                    existing_atom.ins_code == ins_code and
                    existing_atom.atom == atom):
                return True
        return False

def __init__(self):
        self.my_z_offset = 7.5  # Example offset, replace as needed

def is_already_loaded(self, prot, res_num, ins_code, atom):
        for existing_atom in prot.Atoms:
            if (existing_atom.res == res_num and
                existing_atom.ins_code == ins_code and
                existing_atom.atom == atom):
                return True
        return False

def get_coords(self, prot, res_num, ins_code, atom):
        for existing_atom in prot.Atoms:
            if (existing_atom.res == res_num and
                existing_atom.ins_code == ins_code and
                existing_atom.atom == atom):
                return existing_atom.x, existing_atom.y, existing_atom.z
        return None

def get_coords_res(self, prot, res_num, ins_code, atom, res_name):
        for existing_atom in prot.Atoms:
            if res_name == "AAA":  # Wild-card residue
                if (existing_atom.res == res_num and
                    existing_atom.ins_code == ins_code and
                    existing_atom.atom == atom):
                    return existing_atom.x, existing_atom.y, existing_atom.z
            elif (existing_atom.res == res_num and
                  existing_atom.ins_code == ins_code and
                  existing_atom.atom == atom and
                  existing_atom.res_name == res_name):
                return existing_atom.x, existing_atom.y, existing_atom.z
        return None

def calc_docking_angle(self):
        # Placeholder for variables
        tcrx, tcry, tcrz = 0.0, 0.0, 0.0
        tcrcentx, tcrcenty, tcrcentz = 0.0, 0.0, 0.0
        tcrrotx, tcrroty, tcrrotz = 0.0, 0.0, 0.0
        mhcx, mhcy, mhcz = 0.0, 0.0, 0.0
        mhccentx, mhccenty, mhccentz = 0.0, 0.0, 0.0
        mhcnormx, mhcnormy, mhcnormz = 0.0, 0.0, 0.0

        # Calculations for vectors (need to replace)
        self.calc_tcr_sg_vector(tcrx, tcry, tcrz, tcrcentx, tcrcenty, tcrcentz)
        self.calc_tcr_rot_vector(tcrrotx, tcrroty, tcrrotz)
        self.calc_mhc_vectors(mhcx, mhcy, mhcz, mhccentx, mhccenty, mhccentz, mhcnormx, mhcnormy, mhcnormz)
        self.align_tcr_rot_vector(tcrrotx, tcrroty, tcrrotz, tcrcentx, tcrcenty, tcrcentz)

        # Calculate docking angle
        angle = math.acos(tcrx * mhcx + tcry * mhcy + tcrz * mhcz) * 180.0 / self.PI

        # Offset calculations
        line_dist = (self.my_z_offset - 7.5 - tcrcentz) / tcrrotz
        x_off = tcrrotx * line_dist + tcrcentx
        y_off = tcrroty * line_dist + tcrcenty

        # Cross product calculation for normal vector
        tmpx, tmpy, tmpz = self.calc_cross_prod(
            tcrcentx - mhccentx, tcrcenty - mhccenty, tcrcentz - mhccentz,
            tcrcentx - mhccentx - mhcnormx, tcrcenty - mhccenty - mhcnormy, tcrcentz - mhccentz - mhcnormz
        )

        # Normalize cross product
        mag = math.sqrt(tmpx**2 + tmpy**2 + tmpz**2)
        tmpx /= mag
        tmpy /= mag
        tmpz /= mag

        # Output docking angle
        print(f"Docking angle: {angle}")
        print(f"Offsets: x: {x_off}, y: {y_off}")
        return angle