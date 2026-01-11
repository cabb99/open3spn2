from ..parser.xyz import parse_xyz
class DNAForce(object):
    """ Wrapper for the openMM force. """

    def __init__(self, dna, OpenCLPatch=True):
        self.periodic = dna.periodic
        self.force = None
        self.dna = dna
        # The patch allows the crosstacking force to run in OpenCL
        # introducing a small difference in the crosstacking energy
        self.OpenCLPatch = OpenCLPatch

        # Define the dna force
        self.reset()

        # Define the interaction pairs
        self.defineInteraction()

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        elif 'force' in self.__dict__:
            return getattr(self.force, attr)
        else:
            if '__repr__' in self.__dict__:
                raise AttributeError(f"type object {str(self)} has no attribute {str(attr)}")
            else:
                raise AttributeError()

    def computeEnergy(self, system, trajectory):
        # Parse trajectory
        traj = parse_xyz('Tests/adna/traj.xyz')

        # clear all forces on the system
        system.clearForces()
        # setup the force
        self.setUpInteraction()
        # for each item of the table:
        # add the force item

        # compute the energy for every frame

        # return a Table with the energy

    def computeSingleEnergy(self, system, trajectory):
        # Parse trajectory
        traj = parse_xyz('Tests/adna/traj.xyz')
        # for each item of the table:
        # clear all forces on the system
        # system.clearForces()
        # setup the force
        # self.setUpInteraction()
        # add the force item

        # compute the energy for every frame

        # return a table with the energy

class ProteinDNAForce(DNAForce):
    def __init__(self, dna, protein):
        self.protein = protein
        super().__init__(dna)

