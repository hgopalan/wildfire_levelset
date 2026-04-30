#!/usr/bin/env python3
"""
Parser for prodfunc.cpp - generates one YAML file per reaction block.

Usage:
    python parse_prodfunc.py
    python parse_prodfunc.py --src /path/to/prodfunc.cpp --out /path/to/reactions/
"""

import re
import os
import sys
import argparse
import yaml

SPECIES_NAMES = [
    'E',      # 0  - electron
    'H',      # 1
    'Hn',     # 2  - H-
    'H2',     # 3
    'H2v',    # 4  - H2 vibrational
    'H2p',    # 5  - H2+
    'H3p',    # 6  - H3+
    'O',      # 7
    'On',     # 8  - O-
    'O2',     # 9
    'O2n',    # 10 - O2-
    'CO',     # 11
    'CO2',    # 12
    'CO2v1',  # 13
    'CO2v2',  # 14
    'CO2v3',  # 15
    'CO2v4',  # 16
    'CO2p',   # 17 - CO2+
    'C2O3p',  # 18 - C2O3+
    'C2O4p',  # 19 - C2O4+
    'CO3n',   # 20 - CO3-
    'C',      # 21
    'OH',     # 22
    'OHp',    # 23 - OH+
    'H2O',    # 24
    'H2Op',   # 25 - H2O+
    'H3Op',   # 26 - H3O+
    'HCO',    # 27
    'CHp',    # 28 - CH+
    'CH2',    # 29
    'CH2p',   # 30 - CH2+
    'CH3',    # 31
    'CH3p',   # 32 - CH3+
    'CH4',    # 33
    'CH2O',   # 34
    'CH3O',   # 35
    'CH3OH',  # 36
    'CH2OH',  # 37
    'AR',     # 38
    'ARe',    # 39 - Ar*
    'ARp',    # 40 - Ar+
    'ARHp',   # 41 - ArH+
    'AR2e',   # 42 - Ar2*
    'AR2p',   # 43 - Ar2+
]
NUM_SPECIES = 44
NA = 6.02214085774e23

# Map from C++ identifier names (e.g. 'E_ID') to integer indices
ID_MAP = {name + '_ID': i for i, name in enumerate(SPECIES_NAMES)}


def resolve_species(s):
    """Convert a C++ species identifier or integer string to an index."""
    s = s.strip()
    if s in ID_MAP:
        return ID_MAP[s]
    try:
        return int(s)
    except ValueError:
        raise ValueError(f"Unknown species identifier: '{s}'")


# ---------------------------------------------------------------------------
# Block extraction
# ---------------------------------------------------------------------------

def extract_reaction_blocks(content):
    """
    Extract all depth-2 {} blocks within productionRate() that contain
    '// reaction'.  Depth-1 is the function body; depth-2 are the individual
    reaction scopes.
    """
    func_start = content.find('void productionRate')
    if func_start == -1:
        raise ValueError("Could not find productionRate function in source")
    func_body_start = content.find('{', func_start)

    blocks = []
    depth = 0
    i = func_body_start
    block_start = -1

    while i < len(content):
        c = content[i]
        if c == '{':
            depth += 1
            if depth == 2:
                block_start = i
        elif c == '}':
            if depth == 2 and block_start >= 0:
                block_text = content[block_start: i + 1]
                if '// reaction' in block_text:
                    blocks.append(block_text)
                block_start = -1
            depth -= 1
        i += 1

    return blocks


# ---------------------------------------------------------------------------
# Comment extraction
# ---------------------------------------------------------------------------

def parse_comment(block):
    """Collect all // reaction ... comment lines from a block."""
    parts = []
    for line in block.split('\n'):
        line = line.strip()
        if line.startswith('//') and 'reaction' in line and 'amrex::Print' not in line:
            parts.append(line[2:].strip())
    return ' | '.join(parts) if parts else ''


# ---------------------------------------------------------------------------
# Helpers for finding fit coefficients
# ---------------------------------------------------------------------------

def find_coefs_in(text, fit_type):
    """
    Find and parse a Jfit_coefs or Ffit_coefs brace-initialiser list.
    Returns a list of floats, or None if not found.
    """
    pattern = rf'{fit_type}_coefs\s*=\s*\{{([^}}]+)\}}'
    m = re.search(pattern, text, re.DOTALL)
    if not m:
        return None
    nums = re.findall(r'[+-]?(?:\d+\.\d*|\d+|\.\d+)(?:[eE][+-]?\d+)?', m.group(1))
    return [float(v) for v in nums]


def find_A_in(text, fit_type):
    """Find Jfit_A or Ffit_A numeric value."""
    pattern = rf'double\s+{fit_type}_A\s*=\s*([+-]?(?:\d+\.\d*|\d+|\.\d+)(?:[eE][+-]?\d+)?)'
    m = re.search(pattern, text)
    return float(m.group(1)) if m else None


# ---------------------------------------------------------------------------
# Brace-matched sub-block extraction (for if/else)
# ---------------------------------------------------------------------------

def extract_brace_body(text, start):
    """
    Starting from position `start` (which should be at '{'), return the
    content between the outermost matching { }.
    """
    assert text[start] == '{'
    depth = 0
    for i in range(start, len(text)):
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
            if depth == 0:
                return text[start + 1: i], i
    raise ValueError("Unmatched brace")


def split_if_else(flat):
    """
    Split 'if (TeeV < X) { LOW_BODY } else { HIGH_BODY }' from a flattened
    block string.  Returns (threshold, low_body, high_body) or None.
    """
    m = re.search(r'if\s*\(\s*TeeV\s*<\s*([\d.]+)\s*\)', flat)
    if not m:
        return None
    threshold = float(m.group(1))

    # Find opening { of if body
    open1 = flat.find('{', m.end())
    if open1 == -1:
        return None
    low_body, close1 = extract_brace_body(flat, open1)

    # Find 'else' then its opening {
    else_pos = flat.find('else', close1)
    if else_pos == -1:
        return None
    open2 = flat.find('{', else_pos)
    if open2 == -1:
        return None
    high_body, _ = extract_brace_body(flat, open2)

    return threshold, low_body, high_body


# ---------------------------------------------------------------------------
# Arrhenius exp() argument parser
# ---------------------------------------------------------------------------

def parse_exp_arg(arg):
    """
    Parse exp() argument to extract beta (logT coefficient) and invT_coef.

    The expression is of the form:
        (beta) * logT - (Ea) * invT
    where logT = log(T/300) and tc[0] is equivalent.

    Returns (beta, invT_coef) such that:
        exp_value = exp(beta * log(T/300) + invT_coef / T)
    and therefore  Ea = -invT_coef  in the standard form
        k = A * exp(beta * log(T/300) - Ea / T).
    """
    beta = 0.0
    invT_coef = 0.0
    num_pat = r'[+-]?(?:\d+\.\d*|\d+|\.\d+)(?:[eE][+-]?\d+)?'

    # Coefficient of logT (or tc[0])
    for m in re.finditer(
        r'([+-]?)\s*\((' + num_pat + r')\)\s*\*\s*(?:logT|tc\[0\])', arg
    ):
        outer = -1.0 if m.group(1) == '-' else 1.0
        beta += outer * float(m.group(2))

    # Coefficient of invT
    for m in re.finditer(
        r'([+-]?)\s*\((' + num_pat + r')\)\s*\*\s*invT', arg
    ):
        outer = -1.0 if m.group(1) == '-' else 1.0
        invT_coef += outer * float(m.group(2))

    return beta, invT_coef


# ---------------------------------------------------------------------------
# k_max extraction from std::min
# ---------------------------------------------------------------------------

_NUM_PAT = r'[+-]?(?:\d+\.\d*|\d+|\.\d+)(?:[eE][+-]?\d+)?'


def _parse_kmax_expr(expr):
    """
    Evaluate a k_max expression of the form 'A' or 'A * B' where A and B are
    numeric literals (possibly in scientific notation).  Returns a float or
    None when the expression does not match.

    This avoids the use of eval() on arbitrary strings.
    """
    expr = expr.strip()
    # Single number
    m = re.fullmatch(_NUM_PAT, expr)
    if m:
        return float(expr)
    # Product of two numbers: A * B
    m = re.fullmatch(r'(' + _NUM_PAT + r')\s*\*\s*(' + _NUM_PAT + r')', expr)
    if m:
        return float(m.group(1)) * float(m.group(2))
    return None


def extract_kmax(flat):
    """
    Find std::min(k_expr, k_max_expr) and parse k_max_expr.
    Returns float or None.
    """
    m = re.search(r'std::min\(\s*[^,]+,\s*([^)]+)\)', flat)
    if not m:
        return None
    return _parse_kmax_expr(m.group(1))


# ---------------------------------------------------------------------------
# Rate type parser
# ---------------------------------------------------------------------------

def parse_rate(block):
    """
    Determine the rate type and extract parameters.
    Returns a dict with 'rate_type' and associated keys.
    """
    # Flatten whitespace for simpler regex matching
    flat = re.sub(r'\s+', ' ', block)

    # --- Conditional (if TeeV < threshold) ---
    cond = split_if_else(flat)
    if cond is not None:
        threshold, low_body, high_body = cond
        result = {'rate_type': 'conditional', 'TeeV_threshold': threshold}

        low_coefs = find_coefs_in(low_body, 'Jfit')
        low_A = find_A_in(low_body, 'Jfit')
        if low_coefs and low_A is not None:
            result['low_type'] = 'Jfit'
            result['low_Jfit_A'] = low_A
            result['low_Jfit_coefs'] = low_coefs

        high_coefs = find_coefs_in(high_body, 'Ffit')
        high_A = find_A_in(high_body, 'Ffit')
        if high_coefs and high_A is not None:
            result['high_type'] = 'Ffit'
            result['high_Ffit_A'] = high_A
            result['high_Ffit_coefs'] = high_coefs

        return result

    # --- Jfit ---
    if 'Jfit_coefs' in block:
        coefs = find_coefs_in(flat, 'Jfit')
        A = find_A_in(flat, 'Jfit')
        k_max = extract_kmax(flat)
        return {
            'rate_type': 'Jfit',
            'Jfit_A': A,
            'Jfit_coefs': coefs,
            'k_max': k_max,
        }

    # --- Ffit ---
    if 'Ffit_coefs' in block:
        coefs = find_coefs_in(flat, 'Ffit')
        A = find_A_in(flat, 'Ffit')
        return {
            'rate_type': 'Ffit',
            'Ffit_A': A,
            'Ffit_coefs': coefs,
        }

    # --- Arrhenius or constant: inspect the k_f assignment line ---
    # Match 'const double k_f = EXPR;' (may span lines in original, but flat here)
    kf_m = re.search(r'(?:const\s+)?double\s+k_f\s*=\s*(.+?);', flat)
    if not kf_m:
        return {'rate_type': 'constant', 'k_f': 0.0}

    kf_expr = kf_m.group(1).strip()

    if 'exp(' not in kf_expr:
        try:
            val = float(kf_expr)
        except ValueError:
            val = 0.0
        return {'rate_type': 'constant', 'k_f': val}

    # Arrhenius: A * exp(...)
    exp_idx = kf_expr.find('* exp(')
    if exp_idx == -1:
        exp_idx = kf_expr.find('*exp(')
    if exp_idx == -1:
        return {'rate_type': 'constant', 'k_f': 0.0}

    A_str = kf_expr[:exp_idx].strip().rstrip('*').strip()
    try:
        A_val = float(A_str)
    except ValueError:
        A_val = 0.0

    exp_arg_m = re.search(r'exp\((.+)\)', kf_expr)
    if not exp_arg_m:
        return {'rate_type': 'Arrhenius', 'A': A_val, 'beta': 0.0, 'Ea': 0.0}

    beta, invT_coef = parse_exp_arg(exp_arg_m.group(1))
    # Use + 0.0 to avoid negative-zero in YAML output
    Ea   = -invT_coef + 0.0
    beta = beta + 0.0

    return {'rate_type': 'Arrhenius', 'A': A_val, 'beta': beta, 'Ea': Ea}


# ---------------------------------------------------------------------------
# Stoichiometry parser
# ---------------------------------------------------------------------------

def parse_stoichiometry(block):
    """
    Parse wdot[X] += coeff * qdot and wdot[X] -= coeff * qdot lines.
    Returns a list of [species_index, net_delta] pairs (zeros filtered out).
    """
    flat = re.sub(r'\s+', ' ', block)
    stoich = {}

    num_pat = r'(?:\d+\.\d*|\d+|\.\d+)(?:[eE][+-]?\d+)?'
    for m in re.finditer(
        r'wdot\[([^\]]+)\]\s*([+-]=)\s*(?:(' + num_pat + r')\s*\*\s*)?qdot',
        flat
    ):
        sp_idx = resolve_species(m.group(1))
        sign = 1.0 if m.group(2) == '+=' else -1.0
        coef = float(m.group(3)) if m.group(3) else 1.0
        stoich[sp_idx] = stoich.get(sp_idx, 0.0) + sign * coef

    return [[idx, delta] for idx, delta in sorted(stoich.items())
            if abs(delta) > 1e-10]


# ---------------------------------------------------------------------------
# Reactants parser
# ---------------------------------------------------------------------------

def parse_reactants(block):
    """
    Extract species indices from the qf = k_f * (...) expression.
    Returns a list of indices (with repetition for quadratic terms).
    """
    flat = re.sub(r'\s+', ' ', block)
    m = re.search(r'const double qf\s*=\s*(.+?);', flat)
    if not m:
        return []
    qf_expr = m.group(1)
    return [resolve_species(s) for s in re.findall(r'sc\[([^\]]+)\]', qf_expr)]


# ---------------------------------------------------------------------------
# Gibbs reverse term parser
# ---------------------------------------------------------------------------

def parse_gibbs(block):
    """
    Parse the Gibbs exponent from the qr expression.

    The C++ pattern is:
        exp(-(EXPR)) * (refC|refCinv|1.0) * (0.0)
    where EXPR contains g_RT[i] terms.

    Returns dict with 'gibbs_exponent' (list of [idx, coef]) and 'refC_factor'
    (1 = refC, -1 = refCinv, 0 = neither).
    """
    flat = re.sub(r'\s+', ' ', block)
    m = re.search(r'const double qr\s*=\s*(.+?);', flat)
    if not m or 'g_RT' not in m.group(1):
        return {'gibbs_exponent': [], 'refC_factor': 0}

    qr_expr = m.group(1)

    # Extract argument of exp(-( ... ))
    exp_m = re.search(r'exp\(-\(([^()]+)\)\)', qr_expr)
    if not exp_m:
        return {'gibbs_exponent': [], 'refC_factor': 0}

    gibbs_arg = exp_m.group(1)
    num_pat = r'(?:\d+\.\d*|\d+|\.\d+)(?:[eE][+-]?\d+)?'
    pairs = []
    for gm in re.finditer(
        r'([+-]?)\s*(?:(' + num_pat + r')\s*\*\s*)?g_RT\[(\d+)\]',
        gibbs_arg
    ):
        sign = -1.0 if gm.group(1) == '-' else 1.0
        coef = float(gm.group(2)) if gm.group(2) else 1.0
        idx = int(gm.group(3))
        pairs.append([idx, sign * coef])

    refC_factor = 0
    if '(refCinv)' in qr_expr:
        refC_factor = -1
    elif '(refC)' in qr_expr:
        refC_factor = 1

    return {'gibbs_exponent': pairs, 'refC_factor': refC_factor}


# ---------------------------------------------------------------------------
# Energy exchange parser
# ---------------------------------------------------------------------------

def parse_energy_exchange(block):
    """
    Extract comp_ener_exch parameters (rxntype, eexci, elidx).
    Returns a dict or None.
    """
    flat = re.sub(r'\s+', ' ', block)
    if 'comp_ener_exch' not in flat:
        return None

    rt_m = re.search(r'int\s+rxntype\s*=\s*(\d+)', flat)
    ee_m = re.search(r'double\s+eexci\s*=\s*([+-]?(?:\d+\.\d*|\d+|\.\d+)(?:[eE][+-]?\d+)?)', flat)
    el_m = re.search(r'int\s+elidx\s*=\s*([^\s;,)]+)', flat)

    if not (rt_m and ee_m and el_m):
        return None

    elidx_str = el_m.group(1).strip()
    try:
        elidx = int(elidx_str)
    except ValueError:
        elidx = resolve_species(elidx_str)

    return {
        'rxntype': int(rt_m.group(1)),
        'eexci': float(ee_m.group(1)),
        'elidx': elidx,
    }


# ---------------------------------------------------------------------------
# Main block parser
# ---------------------------------------------------------------------------

def parse_block(rxn_id, block):
    """Parse a complete reaction block into a YAML-ready dict."""
    result = {'id': rxn_id}
    result['comment'] = parse_comment(block)

    # Rate parameters
    rate_info = parse_rate(block)
    result.update(rate_info)

    # Stoichiometry, reactants, third body
    result['has_third_body'] = 'const double Corr = mixture' in block
    result['stoichiometry'] = parse_stoichiometry(block)
    result['reactants'] = parse_reactants(block)

    # Gibbs reverse term
    result.update(parse_gibbs(block))

    # Energy exchange (optional)
    ee = parse_energy_exchange(block)
    if ee is not None:
        result['energy_exchange'] = ee

    return result


# ---------------------------------------------------------------------------
# Custom YAML dumper for compact float lists
# ---------------------------------------------------------------------------

class _FloatList(list):
    pass


def _float_list_representer(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)


def _float_representer(dumper, data):
    # Avoid scientific-notation confusion: use repr for compactness
    return dumper.represent_scalar('tag:yaml.org,2002:float', repr(data))


def build_dumper():
    dumper = yaml.Dumper
    dumper.add_representer(_FloatList, _float_list_representer)
    return dumper


def to_yaml_friendly(rxn):
    """Convert a reaction dict so that coefficient lists render inline."""
    out = {}
    for k, v in rxn.items():
        if k in ('Jfit_coefs', 'Ffit_coefs',
                 'low_Jfit_coefs', 'high_Ffit_coefs'):
            out[k] = _FloatList(v) if v is not None else v
        elif k == 'stoichiometry':
            out[k] = [_FloatList(pair) for pair in v]
        elif k == 'gibbs_exponent':
            out[k] = [_FloatList(pair) for pair in v]
        elif k == 'reactants':
            out[k] = _FloatList(v) if v else v
        else:
            out[k] = v
    return out


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(here)

    parser = argparse.ArgumentParser(description='Parse prodfunc.cpp into YAML reaction files.')
    parser.add_argument('--src', default=os.path.join(repo, 'src', 'prodfunc.cpp'),
                        help='Path to prodfunc.cpp')
    parser.add_argument('--out', default=os.path.join(here, 'reactions'),
                        help='Output directory for YAML files')
    args = parser.parse_args()

    with open(args.src) as fh:
        content = fh.read()

    blocks = extract_reaction_blocks(content)
    print(f"Found {len(blocks)} reaction blocks")

    os.makedirs(args.out, exist_ok=True)

    for i, block in enumerate(blocks):
        try:
            rxn = parse_block(i, block)
        except Exception as exc:
            comment = parse_comment(block)
            print(f"  WARNING: block {i} ({comment!r}) failed: {exc}", file=sys.stderr)
            rxn = {'id': i, 'comment': comment, 'parse_error': str(exc)}

        out_path = os.path.join(args.out, f'reaction_{i:03d}.yaml')
        friendly = to_yaml_friendly(rxn)
        with open(out_path, 'w') as fh:
            yaml.dump(friendly, fh, Dumper=build_dumper(),
                      default_flow_style=False, sort_keys=False, allow_unicode=True)

    print(f"Generated {len(blocks)} YAML files in {args.out}")


if __name__ == '__main__':
    main()
