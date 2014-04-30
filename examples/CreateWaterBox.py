from pdbfixer.pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile
from simtk.openmm import *
from simtk.unit import *
import optparse

def parse_cmdln():
        import os
        parser=optparse.OptionParser()
        parser.add_option('-p','--pdb', dest='pdb', type='string')
        parser.add_option('-o','--out-file', dest='outfile', type='string', default='2.5_2.5_2.5_WaterBox.pdb')
        (options, args) = parser.parse_args()
        return (options, args)

def run(options):
        fixer = PDBFixer(options['pdb'])
        fixer.addMissingHydrogens(7.0)
        fixer.addSolvent(boxSize=Vec3(2.5,2.5,2.5)*nanometers, padding=None, 
                    positiveIon='Na+', negativeIon='Cl-', ionicStrength=0.0*molar)
        PDBFile.writeFile(fixer.topology, fixer.positions, open(options['outfile'], 'w'))

if __name__=="__main__":
        (options, args) = parse_cmdln()
        run(vars(options))