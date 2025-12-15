# scope_utils.py
import pickle
from collections import defaultdict

def parse_scope_file(filename, min_fields=5, domain_filter=None):
    """Parse SCOPe des file to dict: domain -> sccs.
    Args:
        filename (str): path to SCOPe des file
        min_fields (int): minimum number of whitespace-separated fields expected
        domain_filter (callable or None): function(domain_name)->bool to include domain. 
            Default: keep domain != '-' and startswith 'd'
    Returns:
        dict: {domain_name: sccs}
    """
    if domain_filter is None:
        domain_filter = lambda dn: dn != '-' and dn.startswith('d')
    out = {}
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < min_fields:
                continue
            sccs = parts[2]
            domain_name = parts[3]
            if domain_filter(domain_name):
                out[domain_name] = sccs
    return out

def save_pickle(data, out_path):
    """Save object to pickle."""
    with open(out_path, 'wb') as f:
        pickle.dump(data, f)
    return out_path

def load_pickle(path):
    """Load pickle and return object."""
    with open(path, 'rb') as f:
        return pickle.load(f)

def invert_mapping(domain_to_sccs):
    """Return sccs -> list(domains)."""
    r = defaultdict(list)
    for d, s in domain_to_sccs.items():
        r[s].append(d)
    return dict(r)

def group_by_level(domain_to_sccs, level=2):
    """Group domains by first `level` parts of sccs string.
    Example: level=2 -> 'a.1' (fold), level=3 -> 'a.1.1' (superfamily).
    """
    r = defaultdict(list)
    for d, s in domain_to_sccs.items():
        key = '.'.join(s.split('.')[:level])
        r[key].append(d)
    return dict(r)