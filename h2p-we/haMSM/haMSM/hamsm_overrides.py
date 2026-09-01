import logging

from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])
log = logging.getLogger("msm_we")
log.setLevel(logging.DEBUG)

import numpy as np
import mdtraj as md
import io

def processCoordinates(self, coords):
    log.debug("Processing coordinates")

    if self.dimReduceMethod == "none":
        nC = np.shape(coords)
        nC = nC[0]
        new_coords = coords[:,:935, :] # Specific to BdpA, Only looking at the 935 atoms of the protein
        data = new_coords.reshape(nC, 3 * 935) # Specific for BdpA, Only looking at the 935 atoms of protein.
        self.nAtoms = 935
        return data

    if self.dimReduceMethod == "pca" or self.dimReduceMethod == "vamp":
        # Still work in Progress, don't use. 
        xyz = io.StringIO()
        nS = np.shape(coords)[0] # segments
        

        for i in range(0,nC): #Iterate through all segs

            for i in np.shape(coord):
        

        ### RMSD dimensionality reduction
        log.warning("Hardcoded selection: Doing dim reduction for Na, Cl. This is only for testing!")
        indTraj = self.reference_structure.topology.select("resid 0 to 58")
        g = md.load('/ocean/projects/mcb180038p/jml230/bdpa_wsh2045_p2_1_rebin/reference/wsh2045_eq3.pdb')
        g.atom_slice(indTraj, inplace=True)
        md.rmsd(g, self.reference_structure)
        
        diff = np.subtract(coords[:, indNA], coords[:, indCL])

        dist = np.array(np.sqrt(np.mean(np.power(diff, 2), axis=-1)))

        return dist


log.info("Loading user-override functions for modelWE")
