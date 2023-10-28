import pandas

def parse_xyz(filename=''):
    columns = ['N', 'timestep', 'id', 'name', 'x', 'y', 'z']
    data = []
    with open(filename, 'r') as traj_file:
        atom = pandas.Series(index=columns)
        atom['id'] = None
        for line in traj_file:
            s = line.split()
            if len(s) == 1:
                atom['N'] = int(s[0])
                if atom['id'] > -1:
                    assert atom['id'] == atoms
                atoms = int(s[0])
            elif len(s) == 3:
                atom['timestep'] = int(s[2])
                atom['id'] = 0
            elif len(s) > 3:
                atom['name'] = int(s[0])
                atom['x'], atom['y'], atom['z'] = [float(a) for a in s[1:4]]
                data += [atom.copy()]
                atom['id'] += 1
    xyz_data = pandas.concat(data, axis=1).T
    for i in ['N', 'timestep', 'id', 'name']:
        xyz_data[i] = xyz_data[i].astype(int)
    return xyz_data
