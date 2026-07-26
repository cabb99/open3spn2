import numpy as np
import pandas
import itertools


import openmm
import openmm.unit as unit
from .template import DNAForce
_af = 1 * unit.degree / unit.radian  # angle scaling factor
_dnaResidues = ['DA', 'DC', 'DT', 'DG']
_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
# Base-type codes for the in-expression base-pair mask in Exclusion. Complementary pairs are exactly
# the pairs whose codes sum to 5 (A+T = 1+4, G+C = 2+3); no other pair of codes (0 for P, S, non-DNA)
# sums to 5, so a single window on the sum masks exactly the base-pair exclusions.
_base_type = {'A': 1, 'G': 2, 'C': 3, 'T': 4}


class Bond(DNAForce, openmm.CustomBondForce):
    def __init__(self, dna, k=1, k_name=None, force_group=6, OpenCLPatch=True):
        self.k = k
        self.k_name = k_name or 'k_bond'
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def getParameterNames(self):
        self.perInteractionParameters = []
        self.GlobalParameters = []
        for i in range(self.force.getNumPerBondParameters()):
            self.perInteractionParameters += [self.force.getPerBondParameterName(i)]
        for i in range(self.force.getNumGlobalParameters()):
            self.GlobalParameters += [self.force.getGlobalParameterName(i)]
        return [self.perInteractionParameters, self.GlobalParameters]

    def reset(self):
        bondForce = openmm.CustomBondForce(f"{self.k_name}*(Kb2*(r-r0)^2+Kb3*(r-r0)^3+Kb4*(r-r0)^4)")
        bondForce.addPerBondParameter('r0')
        bondForce.addPerBondParameter('Kb2')
        bondForce.addPerBondParameter('Kb3')
        bondForce.addPerBondParameter('Kb4')
        bondForce.addGlobalParameter(self.k_name, self.k)
        bondForce.setUsesPeriodicBoundaryConditions(self.periodic)
        bondForce.setForceGroup(self.force_group)
        self.force = bondForce

    def defineInteraction(self):
        for i, b in self.dna.bonds.iterrows():
            parameters = [b['r0'],
                          b['Kb2'],
                          b['Kb3'],
                          b['Kb4']]
            self.force.addBond(int(b['aai']), int(b['aaj']), parameters)


class Angle(DNAForce, openmm.CustomAngleForce):
    def __init__(self, dna, k=1, k_name=None, force_group=7, OpenCLPatch=True):
        self.k = k
        self.k_name = k_name or 'k_angle'
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        # A CustomAngleForce (rather than HarmonicAngleForce) so it can carry the k multiplier;
        # with the per-angle stiffness k = 2*epsilon this is 0.5*k*(theta-t0)^2 = epsilon*(theta-t0)^2.
        angleForce = openmm.CustomAngleForce(f"{self.k_name}*0.5*k*(theta-t0)^2")
        angleForce.addPerAngleParameter('t0')
        angleForce.addPerAngleParameter('k')
        angleForce.addGlobalParameter(self.k_name, self.k)
        angleForce.setUsesPeriodicBoundaryConditions(self.periodic)
        angleForce.setForceGroup(self.force_group)
        self.force = angleForce

    def defineInteraction(self):
        for i, a in self.dna.angles.iterrows():
            parameters = [a['t0'] * _af,
                          a['epsilon'] * 2]
            self.force.addAngle(int(a['aai']), int(a['aaj']), int(a['aak']), parameters)


class Stacking(DNAForce, openmm.CustomCompoundBondForce):
    def __init__(self, dna, k=1, k_name=None, force_group=8, OpenCLPatch=True):
        self.k = k
        self.k_name = k_name or 'k_stacking'
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        stackingForce = openmm.CustomCompoundBondForce(3, f"""{self.k_name}*energy;
                        energy=rep+f2*attr;
                        rep=epsilon*(1-exp(-alpha*(dr)))^2*step(-dr);
                        attr=epsilon*(1-exp(-alpha*(dr)))^2*step(dr)-epsilon;
                        dr=distance(p2,p3)-sigma;
                        f2=max(f*pair2,pair1);
                        pair1=step(dt+pi/2)*step(pi/2-dt);
                        pair2=step(dt+pi)*step(pi-dt);
                        f=1-cos(dt)^2;
                        dt=rng*(angle(p1,p2,p3)-t0);""")
        stackingForce.setUsesPeriodicBoundaryConditions(self.periodic)
        stackingForce.addPerBondParameter('epsilon')
        stackingForce.addPerBondParameter('sigma')
        stackingForce.addPerBondParameter('t0')
        stackingForce.addPerBondParameter('alpha')
        stackingForce.addPerBondParameter('rng')
        stackingForce.addGlobalParameter('pi', np.pi)
        stackingForce.addGlobalParameter(self.k_name, self.k)
        stackingForce.setForceGroup(self.force_group)
        self.force = stackingForce

    def defineInteraction(self):
        for i, a in self.dna.stackings.iterrows():
            parameters = [a['epsilon'],
                          a['sigma'],
                          a['t0'] * _af,
                          a['alpha'],
                          a['rng']]
            self.force.addBond([a['aai'], a['aaj'], a['aak']], parameters)


class Dihedral(DNAForce, openmm.CustomTorsionForce):
    def __init__(self, dna, k=1, k_name=None, force_group=9, OpenCLPatch=True):
        self.k = k
        self.k_name = k_name or 'k_dihedral'
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        dihedralForce = openmm.CustomTorsionForce(f"""{self.k_name}*energy;
                        energy = K_periodic*(1-cs)-K_gaussian*exp(-dt_periodic^2/2/sigma^2);
                        cs = cos(dt);
                        dt_periodic = dt-floor((dt+pi)/(2*pi))*(2*pi);
                        dt = theta-t0""")
        dihedralForce.setUsesPeriodicBoundaryConditions(self.periodic)
        dihedralForce.addPerTorsionParameter('K_periodic')
        dihedralForce.addPerTorsionParameter('K_gaussian')
        dihedralForce.addPerTorsionParameter('sigma')
        dihedralForce.addPerTorsionParameter('t0')
        dihedralForce.addGlobalParameter('pi', np.pi)
        dihedralForce.addGlobalParameter(self.k_name, self.k)
        dihedralForce.setForceGroup(self.force_group)
        self.force = dihedralForce

    def defineInteraction(self):
        for i, a in self.dna.dihedrals.iterrows():
            parameters = [a['K_dihedral'],
                          a['K_gaussian'],
                          a['sigma'],
                          (180 + a['t0']) * _af]
            particles = [a['aai'], a['aaj'], a['aak'], a['aal']]
            self.force.addTorsion(*particles, parameters)


class BasePair(DNAForce, openmm.CustomHbondForce):
    def __init__(self, dna, k=1, k_name=None, force_group=10, OpenCLPatch=True):
        self.k = k
        self.k_name = k_name or 'k_basepair'
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        def basePairForce():
            pairForce = openmm.CustomHbondForce(f'''{self.k_name}*energy;
                        energy=rep+1/2*(1+cos(dphi))*fdt1*fdt2*attr;
                        rep  = epsilon*(1-exp(-alpha*dr))^2*(1-step(dr));
                        attr = epsilon*(1-exp(-alpha*dr))^2*step(dr)-epsilon;
                        fdt1 = max(f1*pair0t1,pair1t1);
                        fdt2 = max(f2*pair0t2,pair1t2);
                        pair1t1 = step(pi/2+dt1)*step(pi/2-dt1);
                        pair1t2 = step(pi/2+dt2)*step(pi/2-dt2);
                        pair0t1 = step(pi+dt1)*step(pi-dt1);
                        pair0t2 = step(pi+dt2)*step(pi-dt2);
                        f1 = 1-cos(dt1)^2;
                        f2 = 1-cos(dt2)^2;
                        dphi = dihedral(d2,d1,a1,a2)-phi0;
                        dr    = distance(d1,a1)-sigma;
                        dt1   = rng*(angle(d2,d1,a1)-t01);
                        dt2   = rng*(angle(a2,a1,d1)-t02);''')
            if self.periodic:
                pairForce.setNonbondedMethod(pairForce.CutoffPeriodic)
            else:
                pairForce.setNonbondedMethod(pairForce.CutoffNonPeriodic)
            pairForce.setCutoffDistance(1.8)  # Paper
            pairForce.addPerDonorParameter('phi0')
            pairForce.addPerDonorParameter('sigma')
            pairForce.addPerDonorParameter('t01')
            pairForce.addPerDonorParameter('t02')
            pairForce.addPerDonorParameter('rng')
            pairForce.addPerDonorParameter('epsilon')
            pairForce.addPerDonorParameter('alpha')
            pairForce.addGlobalParameter('pi', np.pi)
            pairForce.addGlobalParameter(self.k_name, self.k)
            self.force = pairForce
            pairForce.setForceGroup(self.force_group)
            return pairForce

        basePairForces = {}
        pair_definition = self.dna.pair_definition[self.dna.pair_definition['DNA'] == self.dna.DNAtype]
        for i, pair in pair_definition.iterrows():
            basePairForces.update({i: basePairForce()})
        self.forces = basePairForces

    def defineInteraction(self):
        pair_definition = self.dna.pair_definition[self.dna.pair_definition['DNA'] == self.dna.DNAtype]
        atoms = self.dna.atoms.copy()
        atoms['index'] = atoms.index
        atoms.index = zip(atoms['chainID'], atoms['resSeq'], atoms['name'])
        is_dna = atoms['resname'].isin(_dnaResidues)

        for i, pair in pair_definition.iterrows():
            D1 = atoms[(atoms['name'] == pair['Base1']) & is_dna].copy()
            A1 = atoms[(atoms['name'] == pair['Base2']) & is_dna].copy()

            try:
                D2 = atoms.loc[[(c, r, 'S') for c, r, n in D1.index]]
            except KeyError:
                for c, r, n in D1.index:
                    if (c, r, 'S') not in atoms.index:
                        print(f'Residue {c}:{r} does not have a Sugar atom (S)')
                raise KeyError

            try:
                A2 = atoms.loc[[(c, r, 'S') for c, r, n in A1.index]]
            except KeyError:
                for c, r, n in A1.index:
                    if (c, r, 'S') not in atoms.index:
                        print(f'Residue {c}:{r} does not have a Sugar atom (S)')
                raise KeyError

            D1_list = list(D1['index'])
            A1_list = list(A1['index'])
            D2_list = list(D2['index'])
            A2_list = list(A2['index'])

            # Define parameters
            parameters = [pair.torsion * _af,
                          pair.sigma,
                          pair.t1 * _af,
                          pair.t2 * _af,
                          pair.rang,
                          pair.epsilon,
                          pair.alpha]

            # Add donors and acceptors
            # Here I am including the same atom twice,
            # it doesn't seem to break things
            for d1, d2 in zip(D1_list, D2_list):
                self.forces[i].addDonor(d1, d2, -1, parameters)
            for a1, a2 in zip(A1_list, A2_list):
                self.forces[i].addAcceptor(a1, a2, -1)
            # Exclude interactions
            D1['donor_id'] = [i for i in range(len(D1))]
            A1['aceptor_id'] = [i for i in range(len(A1))]

            for (_i, atom_a), (_j, atom_b) in itertools.product(D1.iterrows(), A1.iterrows()):
                # Neighboring residues
                # The sequence exclusion was reduced to two residues
                # since the maximum number of exclusions in OpenCL is 4.
                # In the original 3SPN2 it was 3 residues (6 to 9)
                # This change has no noticeable effect
                if (atom_a.chainID == atom_b.chainID) and (abs(atom_a.resSeq - atom_b.resSeq) <= 2):
                    self.forces[i].addExclusion(atom_a['donor_id'], atom_b['aceptor_id'])

    def addForce(self, system):
        for f in self.forces:
            system.addForce(self.forces[f])


class BasePair2(DNAForce, openmm.CustomHbondForce):
    def __init__(self, dna, k=1, k_name=None, force_group=10, OpenCLPatch=True):
        self.k = k
        self.k_name = k_name or 'k_basepair'
        self.force_group = force_group
        # moves pair exclusions into the hamiltonian
        # so we don't need to use the openmm exclusion list;
        # also applies a minimum sequence separation of 5,
        # consistent with the lammps 3SPN2 code
        # (actually, we'll use 3 to be consistent with the OpenCLPatch version)
        self.min_seq_sep = 3
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        def basePairForce():
            pairForce = openmm.CustomHbondForce(f'''{self.k_name}*seqsep*energy;
                        seqsep=max(1-delta(chainID_a-chainID_d),step(abs(resSeq_a-resSeq_d)-{self.min_seq_sep}));
                        energy=rep+1/2*(1+cos(dphi))*fdt1*fdt2*attr;
                        rep  = epsilon*(1-exp(-alpha*dr))^2*(1-step(dr));
                        attr = epsilon*(1-exp(-alpha*dr))^2*step(dr)-epsilon;
                        fdt1 = max(f1*pair0t1,pair1t1);
                        fdt2 = max(f2*pair0t2,pair1t2);
                        pair1t1 = step(pi/2+dt1)*step(pi/2-dt1);
                        pair1t2 = step(pi/2+dt2)*step(pi/2-dt2);
                        pair0t1 = step(pi+dt1)*step(pi-dt1);
                        pair0t2 = step(pi+dt2)*step(pi-dt2);
                        f1 = 1-cos(dt1)^2;
                        f2 = 1-cos(dt2)^2;
                        dphi = dihedral(d2,d1,a1,a2)-phi0;
                        dr    = distance(d1,a1)-sigma;
                        dt1   = rng*(angle(d2,d1,a1)-t01);
                        dt2   = rng*(angle(a2,a1,d1)-t02);''')
            if self.periodic:
                pairForce.setNonbondedMethod(pairForce.CutoffPeriodic)
            else:
                pairForce.setNonbondedMethod(pairForce.CutoffNonPeriodic)
            pairForce.setCutoffDistance(1.8)  # Paper
            pairForce.addPerDonorParameter('phi0')
            pairForce.addPerDonorParameter('sigma')
            pairForce.addPerDonorParameter('t01')
            pairForce.addPerDonorParameter('t02')
            pairForce.addPerDonorParameter('rng')
            pairForce.addPerDonorParameter('epsilon')
            pairForce.addPerDonorParameter('alpha')
            pairForce.addPerDonorParameter('chainID_d')
            pairForce.addPerDonorParameter('resSeq_d')
            pairForce.addPerAcceptorParameter('chainID_a')
            pairForce.addPerAcceptorParameter('resSeq_a')
            pairForce.addGlobalParameter('pi', np.pi)
            pairForce.addGlobalParameter(self.k_name, self.k)
            self.force = pairForce
            pairForce.setForceGroup(self.force_group)
            return pairForce

        basePairForces = {}
        pair_definition = self.dna.pair_definition[self.dna.pair_definition['DNA'] == self.dna.DNAtype]
        for i, pair in pair_definition.iterrows():
            basePairForces.update({i: basePairForce()})
        self.forces = basePairForces

    def defineInteraction(self):
        pair_definition = self.dna.pair_definition[self.dna.pair_definition['DNA'] == self.dna.DNAtype]
        atoms = self.dna.atoms.copy()
        atoms['index'] = atoms.index
        atoms.index = zip(atoms['chainID'], atoms['resSeq'], atoms['name'])
        is_dna = atoms['resname'].isin(_dnaResidues)

        for i, pair in pair_definition.iterrows():
            D1 = atoms[(atoms['name'] == pair['Base1']) & is_dna].copy()
            A1 = atoms[(atoms['name'] == pair['Base2']) & is_dna].copy()

            try:
                D2 = atoms.loc[[(c, r, 'S') for c, r, n in D1.index]]
            except KeyError:
                for c, r, n in D1.index:
                    if (c, r, 'S') not in atoms.index:
                        print(f'Residue {c}:{r} does not have a Sugar atom (S)')
                raise KeyError

            try:
                A2 = atoms.loc[[(c, r, 'S') for c, r, n in A1.index]]
            except KeyError:
                for c, r, n in A1.index:
                    if (c, r, 'S') not in atoms.index:
                        print(f'Residue {c}:{r} does not have a Sugar atom (S)')
                raise KeyError

            D1_list = list(D1['index'])
            A1_list = list(A1['index'])
            D2_list = list(D2['index'])
            A2_list = list(A2['index'])

            # Define parameters
            parameters = [pair.torsion * _af,
                          pair.sigma,
                          pair.t1 * _af,
                          pair.t2 * _af,
                          pair.rang,
                          pair.epsilon,
                          pair.alpha]

            # Add donors and acceptors
            # Here I am including the same atom twice,
            # it doesn't seem to break things
            for d1, d2 in zip(D1_list, D2_list):
                info = D1[D1['index']==d1] # returns a row of a DataFrame, whose attributes are numpy object arrays
                #raise ValueError(f'info.chainID.shape: {info.chainID.shape}, info.resSeq.shape:{info.resSeq.shape})')
                self.forces[i].addDonor(d1, d2, -1, parameters + [ord(info.chainID[0]), int(info.resSeq[0])])
            for a1, a2 in zip(A1_list, A2_list):
                info = A1[A1['index']==a1] # returns a row of a DataFrame, whose attributes are numpy object arrays
                self.forces[i].addAcceptor(a1, a2, -1, [ord(info.chainID[0]), int(info.resSeq[0])])

    def addForce(self, system):
        for f in self.forces:
            system.addForce(self.forces[f])


class CrossStacking(DNAForce):
    def __init__(self, dna, k=1, k_name=None, force_group=11, OpenCLPatch=True):
        self.k = k
        self.k_name = k_name or 'k_crossstacking'
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        def crossStackingForce(parametersOnDonor=False):
            crossForce = openmm.CustomHbondForce(f'''{self.k_name}*energy;
                         energy   = fdt3*fdtCS*attr/2;
                         attr     = epsilon*(1-exp(-alpha*dr))^2*step(dr)-epsilon;
                         fdt3     = max(f1*pair0t3,pair1t3);
                         fdtCS    = max(f2*pair0tCS,pair1tCS);
                         pair0t3  = step(pi+dt3)*step(pi-dt3);
                         pair0tCS = step(pi+dtCS)*step(pi-dtCS);
                         pair1t3  = step(pi/2+dt3)*step(pi/2-dt3);
                         pair1tCS = step(pi/2+dtCS)*step(pi/2-dtCS);
                         f1       = 1-cos(dt3)^2;
                         f2       = 1-cos(dtCS)^2;
                         dr       = distance(d1,a3)-sigma;
                         dt3      = rng_BP*(t3-t03);
                         dtCS     = rng_CS*(tCS-t0CS);
                         tCS      = angle(d2,d1,a3);
                         t3       = acos(cost3lim);
                         cost3lim = min(max(cost3,-0.99),0.99);
                         cost3    = sin(t1)*sin(t2)*cos(phi)-cos(t1)*cos(t2);
                         t1       = angle(d2,d1,a1);
                         t2       = angle(d1,a1,a2);
                         phi      = dihedral(d2,d1,a1,a2);''')
            if self.periodic:
                crossForce.setNonbondedMethod(crossForce.CutoffPeriodic)
            else:
                crossForce.setNonbondedMethod(crossForce.CutoffNonPeriodic)
            crossForce.setCutoffDistance(1.8)  # Paper
            parameters = ['t03', 't0CS', 'rng_CS', 'rng_BP', 'epsilon', 'alpha', 'sigma']
            for p in parameters:
                if parametersOnDonor:
                    crossForce.addPerDonorParameter(p)
                else:
                    crossForce.addPerAcceptorParameter(p)
            crossForce.addGlobalParameter('pi', np.pi)
            crossForce.addGlobalParameter(self.k_name, self.k)
            crossForce.setForceGroup(self.force_group)
            return crossForce

        crossStackingForces = {}
        for base in ['A', 'T', 'G', 'C']:
            crossStackingForces.update({base: (crossStackingForce(), crossStackingForce())})
        self.crossStackingForces = crossStackingForces

    def defineInteraction(self):
        atoms = self.dna.atoms.copy()
        atoms['index'] = atoms.index
        atoms.index = zip(atoms['chainID'], atoms['resSeq'], atoms['name'].replace(['A', 'C', 'T', 'G'], 'B'))
        is_dna = atoms['resname'].isin(_dnaResidues)
        bases = atoms[atoms['name'].isin(['A', 'T', 'G', 'C']) & is_dna]
        D1 = bases
        D2 = atoms.reindex([(c, r, 'S') for c, r, n in bases.index])
        D3 = atoms.reindex([(c, r + 1, 'B') for c, r, n in bases.index])
        A1 = D1
        A2 = D2
        A3 = atoms.reindex([(c, r - 1, 'B') for c, r, n in bases.index])

        # Select only bases where the other atoms exist
        D2.index = D1.index
        D3.index = D1.index
        temp = pandas.concat([D1, D2, D3], axis=1, keys=['D1', 'D2', 'D3'])
        sel = temp[temp['D3', 'name'].isin(['A', 'T', 'G', 'C']) &  # D3 must be a base
                   temp['D2', 'name'].isin(['S']) &  # D2 must be a sugar
                   (temp['D3', 'chainID'] == temp['D1', 'chainID']) &  # D3 must be in the same chain
                   (temp['D2', 'chainID'] == temp['D1', 'chainID'])].index  # D2 must be in the same chain
        D1 = atoms.reindex(sel)
        D2 = atoms.reindex([(c, r, 'S') for c, r, n in sel])
        D3 = atoms.reindex([(c, r + 1, 'B') for c, r, n in sel])

        # Aceptors
        A2.index = A1.index
        A3.index = A1.index
        temp = pandas.concat([A1, A2, A3], axis=1, keys=['A1', 'A2', 'A3'])
        sel = temp[temp['A3', 'name'].isin(['A', 'T', 'G', 'C']) &  # A3 must be a base
                   temp['A2', 'name'].isin(['S']) &  # A2 must be a sugar
                   (temp['A3', 'chainID'] == temp['A1', 'chainID']) &  # A3 must be in the same chain
                   (temp['A2', 'chainID'] == temp['A1', 'chainID'])].index  # A2 must be in the same chain
        A1 = atoms.reindex(sel)
        A2 = atoms.reindex([(c, r, 'S') for c, r, n in sel])
        A3 = atoms.reindex([(c, r - 1, 'B') for c, r, n in sel])

        # Parameters
        cross_definition = self.dna.cross_definition[self.dna.cross_definition['DNA'] == self.dna.DNAtype].copy()
        i = [a for a in zip(cross_definition['Base_d1'], cross_definition['Base_a1'], cross_definition['Base_a3'])]
        cross_definition.index = i

        donors = {i: [] for i in ['A', 'T', 'G', 'C']}
        for donator, donator2, d1, d2, d3 in zip(D1.itertuples(), D3.itertuples(), D1['index'], D2['index'],
                                                 D3['index']):
            d1t = donator.name
            d3t = donator2.name
            c1, c2 = self.crossStackingForces[d1t]
            a1t = _complement[d1t]
            param = cross_definition.loc[[(a1t, d1t, d3t)]].squeeze()
            parameters = [param['t03'] * _af,
                          param['T0CS_2'] * _af,
                          param['rng_cs2'],
                          param['rng_bp'],
                          param['eps_cs2'],
                          param['alpha_cs2'],
                          param['Sigma_2']]
            c1.addDonor(d1, d2, d3)
            c2.addAcceptor(d1, d2, d3, parameters)
            donors[d1t] += [d1]

        aceptors = {i: [] for i in ['A', 'T', 'G', 'C']}
        for aceptor, aceptor2, a1, a2, a3 in zip(A1.itertuples(), A3.itertuples(), A1['index'], A2['index'],
                                                 A3['index']):
            a1t = aceptor.name
            a3t = aceptor2.name
            c1, c2 = self.crossStackingForces[_complement[a1t]]
            d1t = _complement[a1t]
            param = cross_definition.loc[[(d1t, a1t, a3t)]].squeeze()
            parameters = [param['t03'] * _af,
                          param['T0CS_1'] * _af,
                          param['rng_cs1'],
                          param['rng_bp'],
                          param['eps_cs1'],
                          param['alpha_cs1'],
                          param['Sigma_1']]
            c1.addAcceptor(a1, a2, a3, parameters)
            c2.addDonor(a1, a2, a3)
            aceptors[_complement[a1t]] += [a1]

        # Exclusions
        for base in ['A', 'T', 'G', 'C']:
            c1, c2 = self.crossStackingForces[base]
            for ii, i in enumerate(donors[base]):
                for jj, j in enumerate(aceptors[base]):
                    # The sequence exclusion was reduced to two residues
                    # since the maximum number of exclusions in OpenCL is 4.
                    # In the original 3SPN2 it was 3 residues (6 to 9)
                    # This change has a small effect in B-DNA and curved B-DNA
                    # The second change is to make the interaction symetric and dividing the energy over 2
                    # This also reduces the number of exclusions in the force
                    maxn = 6 if self.OpenCLPatch else 9
                    if (self.dna.atoms.at[i, 'chainID'] == self.dna.atoms.at[j, 'chainID'] and abs(i - j) <= maxn) or \
                            (not self.OpenCLPatch and i > j):
                        c1.addExclusion(ii, jj)
                        c2.addExclusion(jj, ii)

    def addForce(self, system):
        for c1, c2 in self.crossStackingForces.values():
            system.addForce(c1)
            system.addForce(c2)

    def getForceGroup(self):
        fg = 0
        for c1, c2 in self.crossStackingForces.values():
            fg = c1.getForceGroup()
            break
        for c1, c2 in self.crossStackingForces.values():
            assert fg == c1.getForceGroup()
            assert fg == c2.getForceGroup()
        return fg


class CrossStacking2(DNAForce):
    def __init__(self, dna, k=1, k_name=None, force_group=11, OpenCLPatch=True):
        self.k = k
        self.k_name = k_name or 'k_crossstacking'
        self.force_group = force_group
        # moves pair exclusions into the hamiltonian
        # so we don't need to use the openmm exclusion list;
        # also applies a minimum sequence separation of 5,
        # consistent with the lammps 3SPN2 code
        # (actually, we'll use 3 to be consistent with the OpenCLPatch version)
        self.min_seq_sep = 3
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        def crossStackingForce(parametersOnDonor=False):
            crossForce = openmm.CustomHbondForce(f'''{self.k_name}*seqsep*energy;
                         seqsep   = max(1-delta(chainID_a-chainID_d),step(abs(resSeq_a-resSeq_d)-{self.min_seq_sep}));
                         energy   = fdt3*fdtCS*attr/2;
                         attr     = epsilon*(1-exp(-alpha*dr))^2*step(dr)-epsilon;
                         fdt3     = max(f1*pair0t3,pair1t3);
                         fdtCS    = max(f2*pair0tCS,pair1tCS);
                         pair0t3  = step(pi+dt3)*step(pi-dt3);
                         pair0tCS = step(pi+dtCS)*step(pi-dtCS);
                         pair1t3  = step(pi/2+dt3)*step(pi/2-dt3);
                         pair1tCS = step(pi/2+dtCS)*step(pi/2-dtCS);
                         f1       = 1-cos(dt3)^2;
                         f2       = 1-cos(dtCS)^2;
                         dr       = distance(d1,a3)-sigma;
                         dt3      = rng_BP*(t3-t03);
                         dtCS     = rng_CS*(tCS-t0CS);
                         tCS      = angle(d2,d1,a3);
                         t3       = acos(cost3lim);
                         cost3lim = min(max(cost3,-0.99),0.99);
                         cost3    = sin(t1)*sin(t2)*cos(phi)-cos(t1)*cos(t2);
                         t1       = angle(d2,d1,a1);
                         t2       = angle(d1,a1,a2);
                         phi      = dihedral(d2,d1,a1,a2);''')
            if self.periodic:
                crossForce.setNonbondedMethod(crossForce.CutoffPeriodic)
            else:
                crossForce.setNonbondedMethod(crossForce.CutoffNonPeriodic)
            crossForce.setCutoffDistance(1.8)  # Paper
            parameters = ['t03', 't0CS', 'rng_CS', 'rng_BP', 'epsilon', 'alpha', 'sigma']
            for p in parameters:
                if parametersOnDonor:
                    crossForce.addPerDonorParameter(p)
                else:
                    crossForce.addPerAcceptorParameter(p)
            crossForce.addPerDonorParameter('chainID_d')
            crossForce.addPerDonorParameter('resSeq_d')
            crossForce.addPerAcceptorParameter('chainID_a')
            crossForce.addPerAcceptorParameter('resSeq_a')
            crossForce.addGlobalParameter('pi', np.pi)
            crossForce.addGlobalParameter(self.k_name, self.k)
            crossForce.setForceGroup(self.force_group)
            return crossForce

        crossStackingForces = {}
        for base in ['A', 'T', 'G', 'C']:
            crossStackingForces.update({base: (crossStackingForce(), crossStackingForce())})
        self.crossStackingForces = crossStackingForces

    def defineInteraction(self):
        atoms = self.dna.atoms.copy()
        atoms['index'] = atoms.index
        atoms.index = zip(atoms['chainID'], atoms['resSeq'], atoms['name'].replace(['A', 'C', 'T', 'G'], 'B'))
        is_dna = atoms['resname'].isin(_dnaResidues)
        bases = atoms[atoms['name'].isin(['A', 'T', 'G', 'C']) & is_dna]
        D1 = bases
        D2 = atoms.reindex([(c, r, 'S') for c, r, n in bases.index])
        D3 = atoms.reindex([(c, r + 1, 'B') for c, r, n in bases.index])
        A1 = D1
        A2 = D2
        A3 = atoms.reindex([(c, r - 1, 'B') for c, r, n in bases.index])

        # Select only bases where the other atoms exist
        D2.index = D1.index
        D3.index = D1.index
        temp = pandas.concat([D1, D2, D3], axis=1, keys=['D1', 'D2', 'D3'])
        sel = temp[temp['D3', 'name'].isin(['A', 'T', 'G', 'C']) &  # D3 must be a base
                   temp['D2', 'name'].isin(['S']) &  # D2 must be a sugar
                   (temp['D3', 'chainID'] == temp['D1', 'chainID']) &  # D3 must be in the same chain
                   (temp['D2', 'chainID'] == temp['D1', 'chainID'])].index  # D2 must be in the same chain
        D1 = atoms.reindex(sel)
        D2 = atoms.reindex([(c, r, 'S') for c, r, n in sel])
        D3 = atoms.reindex([(c, r + 1, 'B') for c, r, n in sel])

        # Aceptors
        A2.index = A1.index
        A3.index = A1.index
        temp = pandas.concat([A1, A2, A3], axis=1, keys=['A1', 'A2', 'A3'])
        sel = temp[temp['A3', 'name'].isin(['A', 'T', 'G', 'C']) &  # A3 must be a base
                   temp['A2', 'name'].isin(['S']) &  # A2 must be a sugar
                   (temp['A3', 'chainID'] == temp['A1', 'chainID']) &  # A3 must be in the same chain
                   (temp['A2', 'chainID'] == temp['A1', 'chainID'])].index  # A2 must be in the same chain
        A1 = atoms.reindex(sel)
        A2 = atoms.reindex([(c, r, 'S') for c, r, n in sel])
        A3 = atoms.reindex([(c, r - 1, 'B') for c, r, n in sel])

        # Parameters
        cross_definition = self.dna.cross_definition[self.dna.cross_definition['DNA'] == self.dna.DNAtype].copy()
        i = [a for a in zip(cross_definition['Base_d1'], cross_definition['Base_a1'], cross_definition['Base_a3'])]
        cross_definition.index = i
        atom_index_to_group_index = {key: [{'donors':{}, 'acceptors':{}}, {'donors':{}, 'acceptors':{}}] for key in self.crossStackingForces}
        # this dictionary looks like {'A': [{'donors':{}, 'acceptors':{}}, {'donors':{}, 'acceptors':{}}],
        #                             'T': [{'donors':{}, 'acceptors':{}}, {'donors':{}, 'acceptors':{}}],
        #                             'G': [{'donors':{}, 'acceptors':{}}, {'donors':{}, 'acceptors':{}}],
        #                             'C': [{'donors':{}, 'acceptors':{}}, {'donors':{}, 'acceptors':{}}]}
        # there are two Forces in each type A/T/G/C, each with donors and acceptors,
        # whose dictionary values (the innermost dictionaries) map tuples of atom indices onto donor/acceptor indices

        donors = {i: [] for i in ['A', 'T', 'G', 'C']}
        for donator, donator2, d1, d2, d3 in zip(D1.itertuples(), D3.itertuples(), D1['index'], D2['index'],
                                                 D3['index']):
            d1t = donator.name
            d3t = donator2.name
            c1, c2 = self.crossStackingForces[d1t]
            a1t = _complement[d1t]
            param = cross_definition.loc[[(a1t, d1t, d3t)]].squeeze()
            #raise ValueError(donator)
            parameters = [param['t03'] * _af,
                          param['T0CS_2'] * _af,
                          param['rng_cs2'],
                          param['rng_bp'],
                          param['eps_cs2'],
                          param['alpha_cs2'],
                          param['Sigma_2']]
            top_info = [ord(donator.chainID), int(donator.resSeq)]
            assert isinstance(top_info[0], int)
            assert isinstance(top_info[1], int)
            atom_index_to_group_index[d1t][0]['donors'].update({(d1,d2,d3): c1.addDonor(d1, d2, d3, top_info)})
            atom_index_to_group_index[d1t][1]['acceptors'].update({(d1,d2,d3): c2.addAcceptor(d1, d2, d3, parameters + top_info)})
            donors[d1t] += [d1]

        aceptors = {i: [] for i in ['A', 'T', 'G', 'C']}
        for aceptor, aceptor2, a1, a2, a3 in zip(A1.itertuples(), A3.itertuples(), A1['index'], A2['index'],
                                                 A3['index']):
            a1t = aceptor.name
            a3t = aceptor2.name
            c1, c2 = self.crossStackingForces[_complement[a1t]]
            d1t = _complement[a1t]
            param = cross_definition.loc[[(d1t, a1t, a3t)]].squeeze()
            parameters = [param['t03'] * _af,
                          param['T0CS_1'] * _af,
                          param['rng_cs1'],
                          param['rng_bp'],
                          param['eps_cs1'],
                          param['alpha_cs1'],
                          param['Sigma_1']]
            top_info = [ord(aceptor.name), int(aceptor.resSeq)]
            assert isinstance(top_info[0], int)
            assert isinstance(top_info[1], int)
            atom_index_to_group_index[_complement[a1t]][0]['acceptors'].update({(a1,a2,a3): c1.addAcceptor(a1, a2, a3, parameters + top_info)})
            atom_index_to_group_index[_complement[a1t]][1]['donors'].update({(a1,a2,a3): c2.addDonor(a1, a2, a3, top_info)})
            aceptors[_complement[a1t]] += [a1]

        # Close-in-sequence donors/acceptors mathematically have 0 energy
        # (due to the seqsep condition),
        # but there could be numerical issues with evaluating the finite
        # but huge morse potential at short distances before multiplying by 0.
        # Also, donor-acceptor pairs sharing (a) common atom(s) 
        # could have issues differentiating the potential, as described in
        # https://github.com/cabb99/openawsem/issues/94.
        # So we will build a small pair exclusion list only for 
        # donor-acceptor pairs with (a) common atom(s).
        for force_type in self.crossStackingForces:
            forces = self.crossStackingForces[force_type]
            for force_index in [0,1]: # forces list should only have two elements
                donor_dict = atom_index_to_group_index[force_type][force_index]['donors']
                acceptor_dict = atom_index_to_group_index[force_type][force_index]['acceptors']
                for donor in donor_dict:
                    for acceptor in acceptor_dict:
                        # donor and acceptor are tuples of openmm particle indices
                        # donor: (d1,d2,d3)
                        # acceptor: (a1,a2,a3)
                        if not set(donor).isdisjoint(set(acceptor)): # checks that set intersection is not empty
                            forces[force_index].addExclusion(donor_dict[donor], acceptor_dict[acceptor]) # get donor/acceptor indices
            self.crossStackingForces[force_type] = forces
        
    def addForce(self, system):
        for c1, c2 in self.crossStackingForces.values():
            system.addForce(c1)
            system.addForce(c2)

    def getForceGroup(self):
        fg = 0
        for c1, c2 in self.crossStackingForces.values():
            fg = c1.getForceGroup()
            break
        for c1, c2 in self.crossStackingForces.values():
            assert fg == c1.getForceGroup()
            assert fg == c2.getForceGroup()
        return fg


def addNonBondedExclusions(dna, force, OpenCLPatch=True):
    """Adds the identity-independent intra-residue and neighboring-residue exclusions that every
    CustomNonbondedForce shares. Complementary base-pair exclusions are NOT added here -- Exclusion
    masks them inside its energy expression -- so all CustomNonbondedForce keep one identical
    exclusion list (required on the CPU/GPU platforms) and mix cleanly with other force fields."""
    is_dna = dna.atoms['resname'].isin(_dnaResidues)
    atoms = dna.atoms.copy()
    selection = atoms[is_dna].sort_index()
    selection['index'] = selection.index
    selection['neighbor'] = selection['chainID'].astype(str) + '_' + (selection['resSeq'] - 1).astype(str)
    selection.index = selection['chainID'].astype(str) + '_' + (selection['resSeq']).astype(str)

    exclusions = []
    for i, neighbor_res, self_res in zip(selection['index'], selection['neighbor'], selection.index):
        # Add exclusions for the same residue
        for j in selection.loc[self_res, 'index']:
            if i > j:
                exclusions += [(j, i)]
        # Add exclusions with the neighboring residue on the same chain
        try:
            for j in selection.loc[neighbor_res, 'index']:
                exclusions += [(j, i)]
        except KeyError:
            continue

    for i, j in set(exclusions):
        force.addExclusion(i, j)


class Exclusion(DNAForce, openmm.CustomNonbondedForce):
    def __init__(self, dna, k=1, k_name=None, force_group=12, OpenCLPatch=True):
        self.k = k
        self.k_name = k_name or 'k_exclusion'
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        # Complementary base pairs are excluded in-expression via `mask` (0 when bt1+bt2 == 5), so a
        # masked pair contributes 0 energy and 0 force -- identical to an addExclusion, but without
        # putting the base-pair exclusions into the shared exclusion list.
        exclusionForce = openmm.CustomNonbondedForce(f"""{self.k_name}*mask*core;
                         core=(epsilon*((sigma/r)^12-2*(sigma/r)^6)+epsilon)*step(sigma-r);
                         mask=1-step(bt1+bt2-4.5)*step(5.5-bt1-bt2);
                         sigma=0.5*(sigma1+sigma2);
                         epsilon=sqrt(epsilon1*epsilon2)""")
        exclusionForce.addPerParticleParameter('epsilon')
        exclusionForce.addPerParticleParameter('sigma')
        exclusionForce.addPerParticleParameter('bt')
        exclusionForce.addGlobalParameter(self.k_name, self.k)
        exclusionForce.setCutoffDistance(1.8)
        exclusionForce.setForceGroup(self.force_group)  # There can not be multiple cutoff distance on the same force group
        if self.periodic:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffPeriodic)
        else:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffNonPeriodic)
        self.force = exclusionForce

    def defineInteraction(self):
        # addParticles
        particle_definition = self.dna.particle_definition[self.dna.particle_definition['DNA'] == self.dna.DNAtype]
        particle_definition.index = particle_definition.name

        # Reduces or increases the cutoff to the maximum particle radius
        self.force.setCutoffDistance(particle_definition.radius.max())

        # Select only dna atoms
        is_dna = self.dna.atoms['resname'].isin(_dnaResidues)
        atoms = self.dna.atoms.copy()
        atoms['is_dna'] = is_dna
        for i, atom in atoms.iterrows():
            if atom.is_dna:
                param = particle_definition.loc[atom['name']]
                parameters = [param.epsilon, param.radius, _base_type.get(atom['name'], 0)]
            else:
                parameters = [0, .1, 0]  # Null energy, some radius, non-base
            self.force.addParticle(parameters)

        # addExclusions
        addNonBondedExclusions(self.dna, self.force)

class Exclusion2(DNAForce, openmm.CustomNonbondedForce):
    def __init__(self, dna, k=1, k_name=None, force_group=12, OpenCLPatch=True):
        self.k = k
        self.k_name = k_name or 'k_exclusion'
        self.force_group = force_group
        # moves pair exclusions into the hamiltonian
        # so we don't need to use the openmm exclusion list;
        # also applies a minimum sequence separation of 5,
        # consistent with the lammps 3SPN2 code
        # (actually, we'll use 2 to be consistent with the regular Exclusion term)
        self.min_seq_sep = 2
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        # Complementary base pairs are excluded in-expression via `mask` (0 when bt1+bt2 == 5), so a
        # masked pair contributes 0 energy and 0 force -- identical to an addExclusion, but without
        # putting the base-pair exclusions into the shared exclusion list.
        exclusionForce = openmm.CustomNonbondedForce(f"""{self.k_name}*mask*seqsep*core;
                         core=(epsilon*((sigma/r)^12-2*(sigma/r)^6)+epsilon)*step(sigma-r);
                         seqsep=max(1-delta(chainID1-chainID2),step(abs(resSeq2-resSeq1)-{self.min_seq_sep}));
                         mask=1-step(bt1+bt2-4.5)*step(5.5-bt1-bt2);
                         sigma=0.5*(sigma1+sigma2);
                         epsilon=sqrt(epsilon1*epsilon2)""")
        exclusionForce.addPerParticleParameter('epsilon')
        exclusionForce.addPerParticleParameter('sigma')
        exclusionForce.addPerParticleParameter('bt')
        exclusionForce.addPerParticleParameter('chainID')
        exclusionForce.addPerParticleParameter('resSeq')
        exclusionForce.addGlobalParameter(self.k_name, self.k)
        exclusionForce.setCutoffDistance(1.8)
        exclusionForce.setForceGroup(self.force_group)  # There can not be multiple cutoff distance on the same force group
        if self.periodic:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffPeriodic)
        else:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffNonPeriodic)
        self.force = exclusionForce

    def defineInteraction(self):
        # addParticles
        particle_definition = self.dna.particle_definition[self.dna.particle_definition['DNA'] == self.dna.DNAtype]
        particle_definition.index = particle_definition.name

        # Reduces or increases the cutoff to the maximum particle radius
        self.force.setCutoffDistance(particle_definition.radius.max())

        # Select only dna atoms
        is_dna = self.dna.atoms['resname'].isin(_dnaResidues)
        atoms = self.dna.atoms.copy()
        atoms['is_dna'] = is_dna
        for i, atom in atoms.iterrows():
            if atom.is_dna:
                param = particle_definition.loc[atom['name']]
                parameters = [param.epsilon, param.radius, _base_type.get(atom['name'], 0), ord(atom.chainID), int(atom.resSeq)]
            else:
                parameters = [0, .1, 0, 0, 0]  # Null energy, some radius, non-base, some chain, some residue
            self.force.addParticle(parameters)

        # addExclusions
        #addNonBondedExclusions(self.dna, self.force)

class Electrostatics(DNAForce, openmm.CustomNonbondedForce):
    def __init__(self, dna, k=1, k_name=None, force_group=13, temperature=300*unit.kelvin, salt_concentration=100*unit.millimolar, OpenCLPatch=True):
        self.k = k
        self.k_name = k_name or 'k_electrostatics'
        self.force_group = force_group
        self.T = temperature
        self.C = salt_concentration
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        T = self.T
        C = self.C
        e = 249.4 - 0.788 * (T / unit.kelvin) + 7.2E-4 * (T / unit.kelvin) ** 2
        a = 1 - 0.2551 * (C / unit.molar) + 5.151E-2 * (C / unit.molar) ** 2 - 6.889E-3 * (C / unit.molar) ** 3
        #print(e, a)
        dielectric = e * a
        # Debye length
        kb = unit.BOLTZMANN_CONSTANT_kB  # Bolztmann constant
        Na = unit.AVOGADRO_CONSTANT_NA  # Avogadro number
        ec = 1.60217653E-19 * unit.coulomb  # proton charge
        pv = 8.8541878176E-12 * unit.farad / unit.meter  # dielectric permittivity of vacuum

        ldby = np.sqrt(dielectric * pv * kb * T / (2.0 * Na * ec ** 2 * C))
        ldby = ldby.in_units_of(unit.nanometer)
        denominator = 4 * np.pi * pv * dielectric / (Na * ec ** 2)
        denominator = denominator.in_units_of(unit.kilocalorie_per_mole**-1 * unit.nanometer**-1)
        #print(ldby, denominator)

        electrostaticForce = openmm.CustomNonbondedForce(f"""{self.k_name}*energy;
                                                                energy=q1*q2*exp(-r/dh_length)/denominator/r;""")
        electrostaticForce.addPerParticleParameter('q')
        electrostaticForce.addGlobalParameter('dh_length', ldby)
        electrostaticForce.addGlobalParameter('denominator', denominator)
        electrostaticForce.addGlobalParameter(self.k_name, self.k)

        electrostaticForce.setCutoffDistance(5)
        if self.periodic:
            electrostaticForce.setNonbondedMethod(electrostaticForce.CutoffPeriodic)
        else:
            electrostaticForce.setNonbondedMethod(electrostaticForce.CutoffNonPeriodic)
        electrostaticForce.setForceGroup(self.force_group)
        self.force = electrostaticForce

    def defineInteraction(self):
        # addParticles
        particle_definition = self.dna.particle_definition[self.dna.particle_definition['DNA'] == self.dna.DNAtype]
        particle_definition.index = particle_definition.name

        # Select only dna atoms
        is_dna = self.dna.atoms['resname'].isin(_dnaResidues)
        atoms = self.dna.atoms.copy()
        atoms['is_dna'] = is_dna

        for i, atom in atoms.iterrows():
            if atom.is_dna:
                param = particle_definition.loc[atom['name']]
                parameters = [param.charge]
            else:
                parameters = [0]  # No charge if it is not DNA
            # print (i,parameters)
            self.force.addParticle(parameters)

        # add neighbor exclusion
        addNonBondedExclusions(self.dna, self.force, self.OpenCLPatch)

class Electrostatics2(DNAForce, openmm.CustomNonbondedForce):
    def __init__(self, dna, k=1, k_name=None, force_group=13, temperature=300*unit.kelvin, salt_concentration=100*unit.millimolar, OpenCLPatch=True):
        self.k = k
        self.k_name = k_name or 'k_electrostatics'
        self.force_group = force_group
        self.T = temperature
        self.C = salt_concentration
        # moves pair exclusions into the hamiltonian
        # so we don't need to use the openmm exclusion list;
        # also applies a minimum sequence separation of 5,
        # consistent with the lammps 3SPN2 code
        # (actually, we'll use 2 to be consistent with the regular Exclusion term)
        self.min_seq_sep = 2
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        T = self.T
        C = self.C
        e = 249.4 - 0.788 * (T / unit.kelvin) + 7.2E-4 * (T / unit.kelvin) ** 2
        a = 1 - 0.2551 * (C / unit.molar) + 5.151E-2 * (C / unit.molar) ** 2 - 6.889E-3 * (C / unit.molar) ** 3
        #print(e, a)
        dielectric = e * a
        # Debye length
        kb = unit.BOLTZMANN_CONSTANT_kB  # Bolztmann constant
        Na = unit.AVOGADRO_CONSTANT_NA  # Avogadro number
        ec = 1.60217653E-19 * unit.coulomb  # proton charge
        pv = 8.8541878176E-12 * unit.farad / unit.meter  # dielectric permittivity of vacuum

        ldby = np.sqrt(dielectric * pv * kb * T / (2.0 * Na * ec ** 2 * C))
        ldby = ldby.in_units_of(unit.nanometer)
        denominator = 4 * np.pi * pv * dielectric / (Na * ec ** 2)
        denominator = denominator.in_units_of(unit.kilocalorie_per_mole**-1 * unit.nanometer**-1)
        #print(ldby, denominator)

        electrostaticForce = openmm.CustomNonbondedForce(f"""{self.k_name}*seqsep*energy;
                                                                energy=q1*q2*exp(-r/dh_length)/denominator/r;
                                                                seqsep=max(1-delta(chainID1-chainID2),step(abs(resSeq2-resSeq1)-{self.min_seq_sep}));""")
        electrostaticForce.addPerParticleParameter('q')
        electrostaticForce.addPerParticleParameter('chainID')
        electrostaticForce.addPerParticleParameter('resSeq')
        electrostaticForce.addGlobalParameter('dh_length', ldby)
        electrostaticForce.addGlobalParameter('denominator', denominator)
        electrostaticForce.addGlobalParameter(self.k_name, self.k)

        electrostaticForce.setCutoffDistance(5)
        if self.periodic:
            electrostaticForce.setNonbondedMethod(electrostaticForce.CutoffPeriodic)
        else:
            electrostaticForce.setNonbondedMethod(electrostaticForce.CutoffNonPeriodic)
        electrostaticForce.setForceGroup(self.force_group)
        self.force = electrostaticForce

    def defineInteraction(self):
        # addParticles
        particle_definition = self.dna.particle_definition[self.dna.particle_definition['DNA'] == self.dna.DNAtype]
        particle_definition.index = particle_definition.name

        # Select only dna atoms
        is_dna = self.dna.atoms['resname'].isin(_dnaResidues)
        atoms = self.dna.atoms.copy()
        atoms['is_dna'] = is_dna

        for i, atom in atoms.iterrows():
            if atom.is_dna:
                param = particle_definition.loc[atom['name']]
                parameters = [param.charge, ord(atom.chainID), int(atom.resSeq)]
            else:
                parameters = [0, 0, 0]  # No charge if it is not DNA, some chain, some residue
            # print (i,parameters)
            self.force.addParticle(parameters)

        # add neighbor exclusion
        #addNonBondedExclusions(self.dna, self.force, self.OpenCLPatch)