"""The unit-conversion spellings of the code before 9460f8d..34f0ce1, rewritten as the units:: names they became.

The token checks compare the code with a text from before a move (an old SisProd.cpp, or a module before a
split). The unit cleanup (b067fe6, 9460f8d, 04e77e2, c895d5c, 34f0ce1) replaced every unit-conversion spelling
by a constant of UnitConversions.h, some with a new value: 98600 and 98066.52 became kPascalPerKgfPerCm2, 273.14
became kZeroCelsiusInKelvin, and so on. A reference rewritten by the same rules reads like the code after the
cleanup, so the comparison stays about the move. The rules are the ones that rewrote the branch
(estagio-13/ferramentas/unidades.py, branch mode); on code already rewritten they change nothing.

The base patches 20 to 23 apply the same rules to the corrected base, with the exact values as literals.
"""
import re

_NUMBER = r'(?<![\w.])(?:{})(?![\w.])'


def _lit(*spellings):
    return _NUMBER.format('|'.join(re.escape(s) for s in spellings))


def _name(constant):
    return 'units::' + constant


def _additive(text, start, end):
    before = text[:start].rstrip()
    after = text[end:].lstrip()
    if before.endswith(('+', '-', '+=', '-=')):
        return True
    return after.startswith('+') or (after.startswith('-') and not after.startswith('->'))


def _near_cubic_foot(text, start, end):
    window = text[max(0, start - 60):end + 60]
    return '35.31467' in window or 'kCubicFootPerCubicMetre' in window


def _qualify(*names):
    return (r'(?<!::)\b(' + '|'.join(names) + r')\b', lambda m: _name(m.group(1)), None)


def _using(*names):
    return (r'(?m)^using sisprod::(?:' + '|'.join(names) + r');[ \t]*\n', lambda m: '', None)


# In the order the commits applied them.
_RULES = [
    # b067fe6: the exact spellings, by name
    _using('kPascalPerKgfPerCm2', 'kSecondsPerDay', 'kPascalSecondPerCentipoise'),
    _qualify('kPascalPerKgfPerCm2', 'kSecondsPerDay', 'kPascalSecondPerCentipoise', 'celsiusToFahrenheit'),
    (_lit('98066.5'), lambda m: _name('kPascalPerKgfPerCm2'), None),
    (_lit('101325.0', '101325.', '101325'), lambda m: _name('kPascalPerAtmosphere'), None),
    (_lit('459.67'), lambda m: _name('kZeroFahrenheitInRankine'), None),
    (_lit('273.15'), lambda m: _name('kZeroCelsiusInKelvin'), _additive),
    # 9460f8d: kgf/cm2 in pascals
    _using('kPascalPerKgfPerCm2Variant', 'kPascalPerKgfPerCm2Coarse', 'kPascalPerKgfPerCm2PvtSim', 'kKgfPerCm2PerPascal'),
    (r'(?<!::)\bkPascalPerKgfPerCm2(?:Variant|Coarse|PvtSim)\b', lambda m: _name('kPascalPerKgfPerCm2'), None),
    (r'\*\s*kKgfPerCm2PerPascal\b', lambda m: '/ ' + _name('kPascalPerKgfPerCm2'), None),
    (_lit('98066.52', '98066.22', '98068.059233', '98600.', '98600'), lambda m: _name('kPascalPerKgfPerCm2'), None),
    (r'\*\s*' + _lit('1.01971621e-5'), lambda m: '/ ' + _name('kPascalPerKgfPerCm2'), None),
    # 04e77e2: psi, atmosphere and bar
    _using('kAtmospherePerKgfPerCm2', 'kPsiPerAtmosphere', 'kPsiPerPascal', 'kPsiPerKgfPerCm2',
           'kAtmosphereInKgfPerCm2', 'kAtmosphereInPsi'),
    (r'\(\s*([^()]+?)\s*\*\s*kAtmospherePerKgfPerCm2\s*\)\s*\*\s*kPsiPerAtmosphere\b',
     lambda m: m.group(1) + ' * ' + _name('kPsiPerKgfPerCm2'), None),
    (r'\*\s*kPsiPerPascal\b', lambda m: '/ ' + _name('kPascalPerPsi'), None),
    _qualify('kPsiPerKgfPerCm2', 'kAtmosphereInKgfPerCm2', 'kAtmosphereInPsi'),
    (r'\(\s*([^()]+?)\s*\*\s*0\.9678411\s*\)\s*\*\s*' + _lit('14.69595'),
     lambda m: m.group(1) + ' * ' + _name('kPsiPerKgfPerCm2'), None),
    (r'\(\s*0\.9678411\s*\)\s*\*\s*' + _lit('14.69595'), lambda m: _name('kPsiPerKgfPerCm2'), None),
    (_lit('0.9678411') + r'\s*\*\s*' + _lit('14.69595'), lambda m: _name('kPsiPerKgfPerCm2'), None),
    (_lit('14.69595') + r'\s*\*\s*' + _lit('0.9678411'), lambda m: _name('kPsiPerKgfPerCm2'), None),
    (_lit('14.223595'), lambda m: _name('kPsiPerKgfPerCm2'), None),
    (_lit('14.6959488', '14.69595', '14.69'), lambda m: _name('kAtmosphereInPsi'), None),
    (_lit('1.033211', '1.033210485', '1.03322745279996', '1.0332274497825', '1.03322745', '1.03322', '1.033'),
     lambda m: _name('kAtmosphereInKgfPerCm2'), None),
    (_lit('1.01971621', '1.0197'), lambda m: _name('kKgfPerCm2PerBar'), None),
    (_lit('0.980665'), lambda m: _name('kBarPerKgfPerCm2'), None),
    (_lit('6894.8'), lambda m: _name('kPascalPerPsi'), None),
    (r'\*\s*' + _lit('0.00014503773800722'), lambda m: '/ ' + _name('kPascalPerPsi'), None),
    # c895d5c: the Celsius and Fahrenheit offsets
    (_lit('273.16', '273.14', '273.23', '273.0', '273.1', '273.', '273', '272.15'),
     lambda m: _name('kZeroCelsiusInKelvin'), _additive),
    (_lit('460.67'), lambda m: _name('kZeroFahrenheitInRankine'), None),
    # 34f0ce1: the barrel and the cubic foot
    _using('kBarrelPerCubicMetre', 'kCubicFootPerCubicMetre'),
    _qualify('kBarrelPerCubicMetre', 'kCubicFootPerCubicMetre'),
    (_lit('6.29'), lambda m: _name('kBarrelPerCubicMetre'), _near_cubic_foot),
    (_lit('35.31467'), lambda m: _name('kCubicFootPerCubicMetre'), None),
    (_lit('0.1589876', '0.158987'), lambda m: _name('kCubicMetrePerBarrel'), None),
]
_COMPILED = [(re.compile(pattern), replace, check) for pattern, replace, check in _RULES]
_PARENTHESISED = re.compile(r'\(\s*(units::k\w+)\s*\)')


def _segments(text):
    """(is_code, piece) pieces of C++ source: comments and string or character literals are not code."""
    out, i, n, start = [], 0, len(text), 0
    while i < n:
        if text.startswith('//', i) or text.startswith('/*', i):
            if i > start:
                out.append((True, text[start:i]))
            if text.startswith('//', i):
                j = text.find('\n', i)
                j = n if j < 0 else j
            else:
                j = text.find('*/', i + 2)
                j = n if j < 0 else j + 2
            out.append((False, text[i:j]))
            i = start = j
            continue
        if text[i] in '"\'':
            if i > start:
                out.append((True, text[start:i]))
            quote, j = text[i], i + 1
            while j < n and text[j] != quote:
                j += 2 if text[j] == '\\' else 1
            j = min(j + 1, n)
            out.append((False, text[i:j]))
            i = start = j
            continue
        i += 1
    if start < n:
        out.append((True, text[start:]))
    return out


def to_unit_names(text):
    """`text` with every unit-conversion spelling of before the cleanup replaced by its units:: name."""
    pieces = []
    for is_code, piece in _segments(text):
        if is_code:
            for pattern, replace, check in _COMPILED:
                parts, position = [], 0
                for match in pattern.finditer(piece):
                    if check is not None and not check(piece, match.start(), match.end()):
                        continue
                    parts.append(piece[position:match.start()])
                    parts.append(replace(match))
                    position = match.end()
                parts.append(piece[position:])
                piece = ''.join(parts)
            piece = _PARENTHESISED.sub(r'\1', piece)
        pieces.append(piece)
    return ''.join(pieces)
