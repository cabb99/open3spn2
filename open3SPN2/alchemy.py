"""Alchemical base transformations for Open3SPN2.

An :class:`AlchemicalTransformation` turns one nucleic base into another through a coupling parameter
lambda, building the single-topology hybrid Hamiltonian

    U(lambda) = (1 - lambda) * U_initial + lambda * U_target.

It builds every identity-dependent DNA force twice -- from the initial DNA and from the target DNA
(``DNA.create_mutant``) -- weighting the two copies with the per-force ``k`` multiplier
(``k_bond``-style global) so their contributions add up to the linear interpolation. Electrostatics is
identity-independent (every DNA base carries charge 0) and is added once. Because ``Exclusion`` encodes
its base-pair exclusions in its energy expression, the two lambda-scaled copies share one identical
exclusion list and coexist on the CPU/GPU platforms -- no special case is needed.
"""
import openmm.unit as unit

from .ff3SPN2 import forces


class AlchemicalTransformation:
    """A single-topology transformation of one or more DNA bases into new identities.

    Parameters
    ----------
    dna : DNA
        The initial-state DNA (already built, with its topology computed).
    mutations : iterable
        ``(chainID, resSeq, target)`` tuples (dicts or a DataFrame also work) with ``target`` in
        {A, T, G, C}. The target-state DNA is ``dna.create_mutant(mutations)``.
    lambda_name : str
        Base name for the coupling; the two state weights are ``f'{lambda_name}_initial'`` and
        ``f'{lambda_name}_target'``.
    """

    def __init__(self, dna, mutations, lambda_name='lambda'):
        self.initial = dna
        self.target = dna.create_mutant(mutations)
        self.mutations = mutations
        self.lambda_name = lambda_name

    @property
    def initial_weight(self):
        return f'{self.lambda_name}_initial'

    @property
    def target_weight(self):
        return f'{self.lambda_name}_target'

    def add_forces(self, system, verbose=False):
        """Build the hybrid on `system` (an ``openmm.System`` or ``open3SPN2.System``).

        Adds two lambda-scaled copies of every identity-dependent DNA force -- the initial copy
        (from ``self.initial``, weight ``initial_weight``) and the target copy (from ``self.target``,
        weight ``target_weight``) -- plus a single electrostatics force. Returns
        ``{force_name: initial-state force}`` for per-force-group energy readout."""
        added = {}
        for force_name, ForceClass in forces.items():
            if verbose:
                print(force_name)
            if force_name == 'Electrostatics':                 # identity-independent -> one copy
                force = ForceClass(self.initial)
                force.addForce(system)
                added[force_name] = force
            else:
                initial = ForceClass(self.initial, k=1.0, k_name=self.initial_weight)
                initial.addForce(system)
                target = ForceClass(self.target, k=0.0, k_name=self.target_weight)
                target.addForce(system)
                added[force_name] = initial
        return added

    def set_lambda(self, context, value):
        """Set the coupling on an OpenMM context (0 -> initial identity, 1 -> target identity)."""
        context.setParameter(self.initial_weight, 1.0 - value)
        context.setParameter(self.target_weight, value)

    def get_lambda(self, context):
        """Return the current coupling value from an OpenMM context."""
        return context.getParameter(self.target_weight)

    def energy_derivative(self, context, energy_unit=unit.kilojoule_per_mole):
        """Return dE/dlambda (default kJ/mol) for thermodynamic integration.

        The coupling is linear, so at fixed coordinates dE/dlambda = U(1) - U(0) exactly (independent
        of lambda); evaluating it this way also captures the base-pair / cross-stacking
        ``CustomHbondForce`` contributions that OpenMM cannot expose through energy parameter
        derivatives."""
        saved = (context.getParameter(self.initial_weight), context.getParameter(self.target_weight))
        context.setParameter(self.initial_weight, 1.0)
        context.setParameter(self.target_weight, 0.0)
        u0 = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(energy_unit)
        context.setParameter(self.initial_weight, 0.0)
        context.setParameter(self.target_weight, 1.0)
        u1 = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(energy_unit)
        context.setParameter(self.initial_weight, saved[0])
        context.setParameter(self.target_weight, saved[1])
        return u1 - u0
