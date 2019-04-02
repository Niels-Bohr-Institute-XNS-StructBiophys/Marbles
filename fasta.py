import re
import sys

print("# Sequence from PDB generator")

if( len(sys.argv) < 2 ):
    print("# Usage (python version 3.6.x): python fasta.py my_pdb.pdb")

else:
    try:
        import mdtraj as md
    except ImportError:
        print("Install MDTraj: http://mdtraj.org")
        sys.exit()

    three2one_code = { "MET":"M", "ALA":"A", "LEU":"L", "VAL":"V", "PHE":"F", \
                        "TYR":"Y", "GLY":"G", "THR":"T", "HIS":"H", "SER":"S", \
                        "LYS":"K", "ILE":"I", "PRO":"P", "ARG":"R", "ASP":"D", \
                        "GLU":"E", "ASN":"N", "GLN":"Q", "CYS":"C", "TRP":"W"}

    pdb = sys.argv[1]

    struct = md.load(pdb)
    topology = struct.topology
    table, bonds = topology.to_dataframe()
    outname = pdb.rstrip('.pdb') + '.fasta.txt'

    with open( outname, "w" ) as file:
        fstr = f"> FASTA sequence {pdb} | {topology.n_residues} residues | Generated with MDTraj\n"
        file.write( fstr )
        for i in range( topology.n_residues ):
            residue = re.sub(r'[^a-zA-Z]', "", str( topology.residue(i) ) )
            single_letter = three2one_code[ residue ]
            file.write( single_letter )

    print( "# Sequence saved in ", outname )
