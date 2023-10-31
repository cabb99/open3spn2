import pandas

def parse_log(filename=''):
    columns = ''
    log_data = []
    with open(filename, 'r') as log_file:
        start = False
        for line in log_file:
            if line[:4] == 'Step':
                columns = line.split()
                start = True
                continue
            if start:
                try:
                    log_data += [[float(a) for a in line.split()]]
                except ValueError:
                    break
    log_data = pandas.DataFrame(log_data, columns=columns)

    try:
        for i in ['Step', 'nbp']:
            log_data[i] = log_data[i].astype(int)
    except KeyError:
        for i in ['Step', 'v_nbp']:
            log_data[i] = log_data[i].astype(int)
    return log_data