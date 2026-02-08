#!/usr/bin/env python3

import re

# Van der Waals radii in Angstroms (standard values)
VDW_RADII = {
    'H': 1.20,
    'C': 1.70,
    'N': 1.55,
    'O': 1.52,
    'S': 1.80,
    'P': 1.80,
    'F': 1.47,
    'CL': 1.75,
    'BR': 1.85,
    'I': 1.98,
}

def extract_element(atom_name):
    """Extract element symbol from atom name (first 1-2 characters)"""
    # Remove leading spaces and numbers
    elem = atom_name.lstrip().lstrip('0123456789')
    # Take first 1-2 characters that form the element
    if len(elem) > 0:
        if len(elem) > 1 and elem[1].islower():
            return elem[:2].upper()
        else:
            return elem[0].upper()
    return 'C'  # default to carbon

def parse_pdb(pdb_file, output_file):
    """Parse PDB file and extract coordinates with VdW radii"""
    atoms = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Parse PDB line format
                atom_num = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                element = line[76:78].strip() if len(line) > 77 else extract_element(atom_name)
                
                # If element column is empty, extract from atom name
                if not element:
                    element = extract_element(atom_name)
                
                # Get VdW radius
                vdw = VDW_RADII.get(element.upper(), 1.70)
                
                # Convert from Angstroms to nm (1 nm = 10 Angstroms)
                x_nm = x / 10.0
                y_nm = y / 10.0
                z_nm = z / 10.0
                vdw_nm = vdw / 10.0
                
                atoms.append((x_nm, y_nm, z_nm, vdw_nm, element, atom_name))
    
    # Write output file
    with open(output_file, 'w') as f:
        f.write(f"# Extracted from {pdb_file}\n")
        f.write("# Format: x(nm) y(nm) z(nm) vdw_radius(nm) element atom_name\n")
        for x, y, z, vdw, element, atom_name in atoms:
            f.write(f"{x:.6f} {y:.6f} {z:.6f} {vdw:.6f} {element} {atom_name}\n")
    
    print(f"Extracted {len(atoms)} atoms from {pdb_file}")
    print(f"Coordinates saved to {output_file}")

if __name__ == '__main__':
    parse_pdb('protein.pdb', 'protein_coords.txt')
