import numpy as np
import pandas
import itertools


import openmm
import openmm.unit as unit
from .template import DNAForce
_af = 1 * unit.degree / unit.radian  # angle scaling factor
_dnaResidues = ['DA', 'DC', 'DT', 'DG']
_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

class Bond(DNAForce, openmm.CustomBondForce):
    def __init__(self, dna, force_group=6, OpenCLPatch=True):
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
        bondForce = openmm.CustomBondForce("Kb2*(r-r0)^2+Kb3*(r-r0)^3+Kb4*(r-r0)^4")
        bondForce.addPerBondParameter('r0')
        bondForce.addPerBondParameter('Kb2')
        bondForce.addPerBondParameter('Kb3')
        bondForce.addPerBondParameter('Kb4')
        bondForce.setUsesPeriodicBoundaryConditions(self.periodic)
        bondForce.setForceGroup(self.force_group)
        self.force = bondForce

    def defineInteraction(self):
        for i, b in self.dna.bonds.iterrows():
            # Units converted from
            parameters = [b['r0'],
                          b['Kb2'],
                          b['Kb3'],
                          b['Kb4']]
            self.force.addBond(int(b['aai']), int(b['aaj']), parameters)

class Angle(DNAForce, openmm.HarmonicAngleForce):

    def __init__(self, dna, force_group=7, OpenCLPatch=True):
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        angleForce = openmm.HarmonicAngleForce()
        angleForce.setUsesPeriodicBoundaryConditions(self.periodic)
        angleForce.setForceGroup(self.force_group)
        self.force = angleForce

    def defineInteraction(self):
        for i, a in self.dna.angles.iterrows():
            parameters = [a['t0'] * _af,
                          a['epsilon'] * 2]
            self.force.addAngle(int(a['aai']), int(a['aaj']), int(a['aak']), *parameters)

class Stacking(DNAForce, openmm.CustomCompoundBondForce):
    def __init__(self, dna, force_group=8, OpenCLPatch=True):
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        stackingForce = openmm.CustomCompoundBondForce(3, """energy;
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
    def __init__(self, dna, force_group=9, OpenCLPatch=True):
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        dihedralForce = openmm.CustomTorsionForce("""energy;
                        energy = K_periodic*(1-cs)-K_gaussian*exp(-dt_periodic^2/2/sigma^2);
                        cs = cos(dt);
                        dt_periodic = dt-floor((dt+pi)/(2*pi))*(2*pi);
                        dt = theta-t0""")
        # dihedralForce=simtk.openmm.CustomTorsionForce("theta/60.")
        dihedralForce.setUsesPeriodicBoundaryConditions(self.periodic)
        dihedralForce.addPerTorsionParameter('K_periodic')
        dihedralForce.addPerTorsionParameter('K_gaussian')
        dihedralForce.addPerTorsionParameter('sigma')
        dihedralForce.addPerTorsionParameter('t0')
        dihedralForce.addGlobalParameter('pi', np.pi)
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
    def __init__(self, dna, force_group=10, OpenCLPatch=True):
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        def basePairForce():
            pairForce = openmm.CustomHbondForce('''energy;
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
            self.force = pairForce
            pairForce.setForceGroup(self.force_group)
            return pairForce

        basePairForces = {}
        pair_definition = self.dna.pair_definition[self.dna.pair_definition['DNA'] == self.dna.DNAtype]
        for i, pair in pair_definition.iterrows():
            basePairForces[i] = basePairForce()
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
                #print(d1, d2, d2, parameters)
            for a1, a2 in zip(A1_list, A2_list):
                self.forces[i].addAcceptor(a1, a2, -1)
                #print(a1, a2, a2)
            # Exclude interactions
            D1['donor_id'] = list(range(len(D1)))
            A1['aceptor_id'] = list(range(len(A1)))

            for (_i, atom_a), (_j, atom_b) in itertools.product(D1.iterrows(), A1.iterrows()):
                # Neighboring residues
                # The sequence exclusion was reduced to two residues
                # since the maximum number of exclusions in OpenCL is 4.
                # In the original 3SPN2 it was 3 residues (6 to 9)
                # This change has no noticeable effect
                if (atom_a.chainID == atom_b.chainID) and (abs(atom_a.resSeq - atom_b.resSeq) <= 2):
                    self.forces[i].addExclusion(atom_a['donor_id'], atom_b['aceptor_id'])
                    #print(_i, _j)

    def addForce(self, system):
        for f in self.forces:
            system.addForce(self.forces[f])


class CrossStacking(DNAForce):
    def __init__(self, dna, force_group=11, OpenCLPatch=True):
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        def crossStackingForce(parametersOnDonor=False):
            crossForce = openmm.CustomHbondForce(
                '''energy;
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
                         phi      = dihedral(d2,d1,a1,a2);'''
            )
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
            crossForce.setForceGroup(self.force_group)
            return crossForce

        crossStackingForces = {}
        for base in ['A', 'T', 'G', 'C']:
            crossStackingForces[base] = (crossStackingForce(), crossStackingForce())
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
        i = list(
            zip(
                cross_definition['Base_d1'],
                cross_definition['Base_a1'],
                cross_definition['Base_a3'],
            )
        )
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
            # print("Donor", d1t, d1, d2, d3)
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
            # print("Aceptor", a1t, a1, a2, a3)
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

def addNonBondedExclusions(dna, force, OpenCLPatch=True):
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

    exclusions_bp = []
    if OpenCLPatch:
        # Add basepair exclusions
        for key in _complement:
            selection_N = atoms[is_dna & (atoms['name'] == key)]
            selection_C = atoms[is_dna & (atoms['name'] == _complement[key])]
            for i, j in itertools.product(selection_N.index, selection_C.index):
                if i > j:
                    exclusions_bp += [(j, i)]

    exclusions = list(set(exclusions+exclusions_bp))
    for i, j in exclusions:
        force.addExclusion(i, j)


class Exclusion(DNAForce, openmm.CustomNonbondedForce):
    def __init__(self, dna, force_group = 12, OpenCLPatch=True):
        self.force_group = force_group
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        exclusionForce = openmm.CustomNonbondedForce("""energy;
                         energy=(epsilon*((sigma/r)^12-2*(sigma/r)^6)+epsilon)*step(sigma-r);
                         sigma=0.5*(sigma1+sigma2);
                         epsilon=sqrt(epsilon1*epsilon2)""")
        exclusionForce.addPerParticleParameter('epsilon')
        exclusionForce.addPerParticleParameter('sigma')
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
                parameters = [param.epsilon,
                              param.radius]
            else:
                parameters = [0, .1]  # Null energy and some radius)
            # print(i, parameters)
            self.force.addParticle(parameters)

        # addExclusions
        addNonBondedExclusions(self.dna, self.force)


class Electrostatics(DNAForce, openmm.CustomNonbondedForce):
    def __init__(self, dna, force_group=13, temperature=300*unit.kelvin, salt_concentration=100*unit.millimolar, OpenCLPatch=True):
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

        electrostaticForce = openmm.CustomNonbondedForce("""energy;
                                                                energy=q1*q2*exp(-r/dh_length)/denominator/r;""")
        electrostaticForce.addPerParticleParameter('q')
        electrostaticForce.addGlobalParameter('dh_length', ldby)
        electrostaticForce.addGlobalParameter('denominator', denominator)

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