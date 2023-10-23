from open3SPN2 import DNA, System, parse_log, forces
import pandas
import numpy as np
import openmm.app
import openmm
import openmm.unit as unit
import pytest
import pandas as pd
from pathlib import Path

# Define a fixture to load the test sets
# Define a module-level variable to store the test sets
_test_sets = pd.read_csv('Tests/test_cases.csv', comment='#')
_ef = 1 * unit.kilocalorie / unit.kilojoule  # energy scaling factor

@pytest.fixture
def skip_platform(request):
    return request.config.getoption("--skip-platform")

id_list = [f"{row['Name']}_{row['Platform']}" for _, row in _test_sets.iterrows()]
class TestEnergies:
    """Tests that the energies are the same as the example outputs from lammps"""
    @pytest.mark.parametrize(
        "name, dna_type, folder, trajectory, log, energy_term, periodic_size, platform",
        [(row['Name'], row['DNA type'], row['Folder'], row['Trajectory'], row['Log'],
          row['Energy term'], row['periodic size'], row['Platform']) for _, row in _test_sets.iterrows()],
        ids=id_list  # Use the function here
    )
    def test_energies(self, name, dna_type, folder, trajectory, log, energy_term, periodic_size, platform, skip_platform):
        if skip_platform and skip_platform.lower() == platform.lower():
            pytest.skip(f"Skipping tests for platform: {platform}")
        folder = Path(folder)
        dna = DNA.fromXYZ(f'{folder}/in00_conf.xyz', dna_type, template_from_X3DNA=False)
        system = System(dna)
        self._test_energy(energy_term, folder/log, folder/trajectory, name, periodic_size, platform, dna, system)

    @pytest.mark.parametrize(
        "name, dna_type, folder, trajectory, log, energy_term, periodic_size, platform",
        [(row['Name'], row['DNA type'], row['Folder'], row['Trajectory'], row['Log'],
          row['Energy term'], row['periodic size'], row['Platform']) for _, row in _test_sets.iterrows()],
        ids=id_list  # Use the function here
    )
    def test_forces(self, name, dna_type, folder, trajectory, log, energy_term, periodic_size, platform, skip_platform):
        if skip_platform and skip_platform.lower() == platform.lower():
            pytest.skip(f"Skipping tests for platform: {platform}")
        folder = Path(folder)
        dna = DNA.fromXYZ(f'{folder}/in00_conf.xyz', dna_type, template_from_X3DNA=False)
        system = System(dna)
        self._test_force(energy_term, folder/log, folder/trajectory, name, periodic_size, platform, dna, system)

  
    def _test_energy(self,
                     log_energy='E_bond',
                     log_file='Tests/adna/sim.log',
                     traj_file='Tests/adna/traj.xyz',
                     force='Bond', periodic_size=94.2,
                     platform_name='Reference', dna=None, system=None):
        self.dna = dna
        self.system = system
        self.system.clearForces()
        self.system.setDefaultPeriodicBoxVectors(*np.diag([2 * periodic_size / 10] * 3))

        log = parse_log(log_file)

        self.system.clearForces()
        f = forces[force]
        tempforce = f(self.dna)
        try:
            tempforce.addForce(self.system)
        except AttributeError:
            self.system.addForce(tempforce)
        energies = self.system.recomputeEnergy(traj_file, platform_name=platform_name)
        d = (energies / _ef - log[log_energy])
        diff = np.sqrt((d ** 2).sum() / len(energies))
        print(f'The difference in the energy term {log_energy} is {diff} Kcal/mol')
        print(f'The DNA type of the system is {self.dna.DNAtype}')
        data = np.array([energies / _ef, np.array(log[log_energy]), np.array(d)])
        results = pandas.DataFrame(data.T, columns=['Openmm energy', 'Lammps energy', 'Difference'])
        # print(data)
        print(pandas.DataFrame(data.T, columns=['Openmm energy', 'Lammps energy', 'Difference']))
        # The maximum error seems to be bonds in curved BDNA (0.002)
        assert diff < 2E-3, diff
        assert len(results.dropna()) == len(results), results

    def _test_force(self,
                    log_energy='E_bond',
                    log_file='Tests/adna/sim.log',
                    traj_file='Tests/adna/traj.xyz',
                    force='Bond', periodic_size=94.2,
                    platform_name='Reference', dna=None, system=None):
        self.dna = dna
        self.system = system
        self.system.clearForces()
        self.system.setDefaultPeriodicBoxVectors(*np.diag([2 * periodic_size / 10] * 3))
        log = parse_log(log_file)
        self.system.clearForces()
        f = forces[force]
        tempforce = f(self.dna)
        try:
            tempforce.addForce(self.system)
        except AttributeError:
            self.system.addForce(tempforce)
        temperature = 300 * openmm.unit.kelvin
        integrator = openmm.LangevinIntegrator(temperature, 1 / openmm.unit.picosecond,
                                                     2 * openmm.unit.femtoseconds)
        platform = openmm.Platform.getPlatformByName(platform_name)
        simulation = openmm.app.Simulation(self.system.top, self.system, integrator, platform)
        simulation.context.setPositions(self.system.coord.getPositions())
        energy_unit = openmm.unit.kilojoule_per_mole
        nan_force_particles = 0
        for i in range(10):
            state = simulation.context.getState(getForces=True)
            ff = state.getForces()
            sf = (np.array(ff) ** 2).sum(axis=1) ** .5
            for j, f in enumerate(sf):
                if np.isnan(f.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)):
                    print(f"Particle {j + 1}/{len(sf)} has force {f} at step {i}")
                    nan_force_particles += 1
            simulation.step(1)
        assert nan_force_particles == 0, "At least one particle has undefined force"

    def _test_energies_slow(self):
        test_sets = pandas.read_csv('Tests/test_cases.csv', comment='#')
        for i, test in test_sets.iterrows():
            dna_type = test['DNA type']
            folder = test['Folder']
            yield self._test_energy, test['Energy term'], f'{folder}/{test.Log}', f'{folder}/{test.Trajectory}', \
                  test['Name'], folder, dna_type, test['periodic size'], test['Platform']

if __name__=='__main__':
    skip_platform=False
    for _, (name, dna_type, folder, trajectory, log, energy_term, periodic_size, platform) in _test_sets.iterrows():
        test=TestEnergies()
        test.test_energies(name, dna_type, folder, trajectory, log, energy_term, periodic_size, platform, skip_platform)
        test.test_forces(name, dna_type, folder, trajectory, log, energy_term, periodic_size, platform, skip_platform)