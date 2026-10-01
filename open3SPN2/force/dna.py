import numpy as np
import pandas
import itertools


import openmm
import openmm.unit as unit
from .template import DNAForce
from ..circular import circular_chain_info, residue_wrapper
_af = 1 * unit.degree / unit.radian  # angle scaling factor
_dnaResidues = ['DA', 'DC', 'DT', 'DG']
_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
# Base-type codes for the in-expression base-pair mask in Exclusion. Complementary pairs are exactly
# the pairs whose codes sum to 5 (A+T = 1+4, G+C = 2+3); no other pair of codes (0 for P, S, non-DNA)
# sums to 5, so a single window on the sum masks exactly the base-pair exclusions.
_base_type = {'A': 1, 'G': 2, 'C': 3, 'T': 4}

def chained_sequence_index(atoms, chain_gap=100, circular_chains=()):
    """Residue index used by the sequence separation masks of the forces.

    Residues keep their relative separation within a chain and consecutive chains are offset by
    `chain_gap`, so a single ``abs(seqid_i - seqid_j)`` comparison masks pairs that are close in
    sequence without also having to compare chain identifiers on every step.

    The gap is widened to the length of the longest circular chain, so pairs on different chains are
    always further apart than the pairs around the end of a circular chain (see `sequence_separation`).

    OpenMM requires every CustomNonbondedForce of a system to carry the same exclusions, so the
    masks are evaluated inside the energy expressions instead of through the exclusion list."""
    lengths = circular_sequence_length(atoms, circular_chains)
    chain_gap = max(chain_gap, int(lengths.max()) if len(lengths) else 0)
    resSeq = atoms['resSeq'].astype(int)
    chains = resSeq.groupby(atoms['chainID'])
    chain_start = chains.min()
    chain_offset = (chains.max() - chain_start + chain_gap).cumsum().shift(fill_value=0)
    return resSeq - atoms['chainID'].map(chain_start) + atoms['chainID'].map(chain_offset)


def circular_sequence_length(atoms, circular_chains=()):
    """Number of residues of the circular chain of each particle, 0 for particles in linear chains."""
    info = circular_chain_info(atoms, circular_chains)
    return atoms['chainID'].map({c: length for c, (start, length) in info.items()}).fillna(0).astype(int)


def sequence_separation(min_seq_sep, seqid1, seqid2, seqlen1=None, seqlen2=None):
    """Energy expression that is 0 for pairs closer than `min_seq_sep` residues in a chain and 1 otherwise.

    If the per-particle lengths of the circular chains are given (`seqlen1`, `seqlen2`, 0 for linear chains),
    the separation of a pair in a circular chain of length L is also measured around the end of the
    chain, ``L - abs(seqid1 - seqid2)``."""
    separation = f'abs({seqid2}-{seqid1})'
    mask = f'step({separation}-{min_seq_sep})'
    if seqlen1 is not None:
        length = f'max({seqlen1},{seqlen2})'
        # Masked if L - min_seq_sep < separation < L, i.e. closer than min_seq_sep around the circle
        mask += f'*(1-step({separation}-{length}+{min_seq_sep}-0.5)*step({length}-0.5-{separation}))'
    return mask


class Bond(DNAForce, openmm.CustomBondForce):
    def __init__(self, dna, k=1, k_name=None, force_group=6):
        self.k = k
        self.k_name = k_name or 'k_bond'
        self.force_group = force_group
        super().__init__(dna)

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
    def __init__(self, dna, k=1, k_name=None, force_group=7):
        self.k = k
        self.k_name = k_name or 'k_angle'
        self.force_group = force_group
        super().__init__(dna)

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
    def __init__(self, dna, k=1, k_name=None, force_group=8):
        self.k = k
        self.k_name = k_name or 'k_stacking'
        self.force_group = force_group
        super().__init__(dna)

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
    def __init__(self, dna, k=1, k_name=None, force_group=9):
        self.k = k
        self.k_name = k_name or 'k_dihedral'
        self.force_group = force_group
        super().__init__(dna)

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
    def __init__(self, dna, k=1, k_name=None, force_group=10):
        self.k = k
        self.k_name = k_name or 'k_basepair'
        self.force_group = force_group
        self.min_seq_sep = 3
        self.circular = bool(getattr(dna, 'circular_chains', ()))
        self.seqlen_names = ('seqlen_d', 'seqlen_a') if self.circular else ()
        super().__init__(dna)

    def reset(self):
        def basePairForce():
            pairForce = openmm.CustomHbondForce(f'''{self.k_name}*seqsep*energy;
                        seqsep={sequence_separation(self.min_seq_sep, 'seqid_d', 'seqid_a', *self.seqlen_names)};
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
            pairForce.addPerDonorParameter('seqid_d')
            pairForce.addPerAcceptorParameter('seqid_a')
            if self.circular:
                pairForce.addPerDonorParameter('seqlen_d')
                pairForce.addPerAcceptorParameter('seqlen_a')
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
        circular_chains = getattr(self.dna, 'circular_chains', ())
        atoms['seqid'] = chained_sequence_index(atoms, circular_chains=circular_chains)
        atoms['seqlen'] = circular_sequence_length(atoms, circular_chains)
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
            for d1, d2, seqid, seqlen in zip(D1_list, D2_list, D1['seqid'], D1['seqlen']):
                self.forces[i].addDonor(d1, d2, -1, parameters + [int(seqid)] + ([int(seqlen)] if self.circular else []))
            for a1, a2, seqid, seqlen in zip(A1_list, A2_list, A1['seqid'], A1['seqlen']):
                self.forces[i].addAcceptor(a1, a2, -1, [int(seqid)] + ([int(seqlen)] if self.circular else []))

    def addForce(self, system):
        for f in self.forces:
            system.addForce(self.forces[f])


class CrossStacking(DNAForce):
    def __init__(self, dna, k=1, k_name=None, force_group=11):
        self.k = k
        self.k_name = k_name or 'k_crossstacking'
        self.force_group = force_group
        # Donor-acceptor pairs closer than min_seq_sep residues along a chain are masked in the
        # energy expression. 3SPN2 masks 3 residues; open3SPN2 masks 2 (min_seq_sep=3), which has
        # a small effect on B-DNA and curved B-DNA energies.
        self.min_seq_sep = 3
        self.circular = bool(getattr(dna, 'circular_chains', ()))
        self.seqlen_names = ('seqlen_d', 'seqlen_a') if self.circular else ()
        super().__init__(dna)

    def reset(self):
        def crossStackingForce(parametersOnDonor=False):
            crossForce = openmm.CustomHbondForce(f'''{self.k_name}*seqsep*energy;
                         seqsep   = {sequence_separation(self.min_seq_sep, 'seqid_d', 'seqid_a', *self.seqlen_names)};
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
            crossForce.addPerDonorParameter('seqid_d')
            crossForce.addPerAcceptorParameter('seqid_a')
            if self.circular:
                crossForce.addPerDonorParameter('seqlen_d')
                crossForce.addPerAcceptorParameter('seqlen_a')
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
        circular_chains = getattr(self.dna, 'circular_chains', ())
        wrap = residue_wrapper(atoms, circular_chains)
        atoms['seqid'] = chained_sequence_index(atoms, circular_chains=circular_chains)
        atoms['seqlen'] = circular_sequence_length(atoms, circular_chains)
        atoms.index = zip(atoms['chainID'], atoms['resSeq'], atoms['name'].replace(['A', 'C', 'T', 'G'], 'B'))
        is_dna = atoms['resname'].isin(_dnaResidues)
        bases = atoms[atoms['name'].isin(['A', 'T', 'G', 'C']) & is_dna]
        D1 = bases
        D2 = atoms.reindex([(c, r, 'S') for c, r, n in bases.index])
        D3 = atoms.reindex([(c, wrap(c, r + 1), 'B') for c, r, n in bases.index])
        A1 = D1
        A2 = D2
        A3 = atoms.reindex([(c, wrap(c, r - 1), 'B') for c, r, n in bases.index])

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
        D3 = atoms.reindex([(c, wrap(c, r + 1), 'B') for c, r, n in sel])

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
        A3 = atoms.reindex([(c, wrap(c, r - 1), 'B') for c, r, n in sel])

        # Parameters
        cross_definition = self.dna.cross_definition[self.dna.cross_definition['DNA'] == self.dna.DNAtype].copy()
        i = [a for a in zip(cross_definition['Base_d1'], cross_definition['Base_a1'], cross_definition['Base_a3'])]
        cross_definition.index = i
        # Particle triplet of each donor and acceptor mapped onto the index its force assigned to it
        groups = {base: [{'donors': {}, 'acceptors': {}} for force in forces]
                  for base, forces in self.crossStackingForces.items()}

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
            seqid = [int(donator.seqid)] + ([int(donator.seqlen)] if self.circular else [])
            groups[d1t][0]['donors'][(d1, d2, d3)] = c1.addDonor(d1, d2, d3, seqid)
            groups[d1t][1]['acceptors'][(d1, d2, d3)] = c2.addAcceptor(d1, d2, d3, parameters + seqid)
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
            seqid = [int(aceptor.seqid)] + ([int(aceptor.seqlen)] if self.circular else [])
            groups[d1t][0]['acceptors'][(a1, a2, a3)] = c1.addAcceptor(a1, a2, a3, parameters + seqid)
            groups[d1t][1]['donors'][(a1, a2, a3)] = c2.addDonor(a1, a2, a3, seqid)
            aceptors[d1t] += [a1]

        # The sequence separation mask already sets the energy of close-in-sequence pairs to zero,
        # but a donor and an acceptor that share a particle define a degenerate geometry whose
        # angles and derivatives are undefined (https://github.com/cabb99/openawsem/issues/94), so
        # those few pairs still need an exclusion.
        for base, forces in self.crossStackingForces.items():
            for force, force_groups in zip(forces, groups[base]):
                acceptors_by_particle = {}
                for particles, acceptor_index in force_groups['acceptors'].items():
                    for particle in particles:
                        acceptors_by_particle.setdefault(particle, set()).add(acceptor_index)
                for particles, donor_index in force_groups['donors'].items():
                    shared = set().union(*[acceptors_by_particle.get(p, set()) for p in particles])
                    for acceptor_index in sorted(shared):
                        force.addExclusion(donor_index, acceptor_index)

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


class Exclusion(DNAForce, openmm.CustomNonbondedForce):
    def __init__(self, dna, k=1, k_name=None, force_group=12):
        self.k = k
        self.k_name = k_name or 'k_exclusion'
        self.force_group = force_group
        # Pairs closer than min_seq_sep residues along a chain are masked in the energy expression;
        # masking one residue (min_seq_sep=2) reproduces the intra-residue and neighboring-residue
        # exclusions of 3SPN2.
        self.min_seq_sep = 2
        self.circular = bool(getattr(dna, 'circular_chains', ()))
        super().__init__(dna)

    def reset(self):
        # `mask` is 0 for complementary base pairs (bt1+bt2 == 5) and `seqsep` is 0 for pairs close
        # in sequence, so those pairs contribute no energy and no force.
        exclusionForce = openmm.CustomNonbondedForce(f"""{self.k_name}*mask*seqsep*core;
                         core=(epsilon*((sigma/r)^12-2*(sigma/r)^6)+epsilon)*step(sigma-r);
                         seqsep={sequence_separation(self.min_seq_sep, 'seqid1', 'seqid2', *(('seqlen1', 'seqlen2') if self.circular else ()))};
                         mask=1-step(bt1+bt2-4.5)*step(5.5-bt1-bt2);
                         sigma=0.5*(sigma1+sigma2);
                         epsilon=sqrt(epsilon1*epsilon2)""")
        exclusionForce.addPerParticleParameter('epsilon')
        exclusionForce.addPerParticleParameter('sigma')
        exclusionForce.addPerParticleParameter('bt')
        exclusionForce.addPerParticleParameter('seqid')
        if self.circular:
            exclusionForce.addPerParticleParameter('seqlen')
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
        circular_chains = getattr(self.dna, 'circular_chains', ())
        atoms['seqid'] = chained_sequence_index(atoms, circular_chains=circular_chains)
        atoms['seqlen'] = circular_sequence_length(atoms, circular_chains)
        for i, atom in atoms.iterrows():
            if atom.is_dna:
                param = particle_definition.loc[atom['name']]
                parameters = [param.epsilon, param.radius, _base_type.get(atom['name'], 0)]
            else:
                parameters = [0, .1, 0]  # Null energy, some radius, non-base
            self.force.addParticle(parameters + [int(atom.seqid)] + ([int(atom.seqlen)] if self.circular else []))


class Electrostatics(DNAForce, openmm.CustomNonbondedForce):
    def __init__(self, dna, k=1, k_name=None, force_group=13, temperature=300*unit.kelvin, salt_concentration=100*unit.millimolar):
        self.k = k
        self.k_name = k_name or 'k_electrostatics'
        self.force_group = force_group
        self.T = temperature
        self.C = salt_concentration
        # Pairs closer than min_seq_sep residues along a chain are masked in the energy expression;
        # masking one residue (min_seq_sep=2) reproduces the intra-residue and neighboring-residue
        # exclusions of 3SPN2.
        self.min_seq_sep = 2
        self.circular = bool(getattr(dna, 'circular_chains', ()))
        super().__init__(dna)

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
                                                                seqsep={sequence_separation(self.min_seq_sep, 'seqid1', 'seqid2', *(('seqlen1', 'seqlen2') if self.circular else ()))};""")
        electrostaticForce.addPerParticleParameter('q')
        electrostaticForce.addPerParticleParameter('seqid')
        if self.circular:
            electrostaticForce.addPerParticleParameter('seqlen')
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
        circular_chains = getattr(self.dna, 'circular_chains', ())
        atoms['seqid'] = chained_sequence_index(atoms, circular_chains=circular_chains)
        atoms['seqlen'] = circular_sequence_length(atoms, circular_chains)

        for i, atom in atoms.iterrows():
            if atom.is_dna:
                param = particle_definition.loc[atom['name']]
                parameters = [param.charge]
            else:
                parameters = [0]  # No charge if it is not DNA
            # print (i,parameters)
            self.force.addParticle(parameters + [int(atom.seqid)] + ([int(atom.seqlen)] if self.circular else []))
