import io
import os
import re
import sys
import tokenize
from glob import glob

from sidekick import fn

CTRL_TOKENS = {
    tokenize.ENCODING,
    tokenize.NEWLINE,
    tokenize.INDENT,
    tokenize.DEDENT,
    tokenize.ENDMARKER,
}

IS_FILE = re.compile(r'''
(?P<indent>\s*)
if\s+os.path.isfile\(\s*
    (?P<string>"[^"]+")
    (?P<params>[^):]+)
\s*\):
''', re.MULTILINE | re.VERBOSE)

PATH = re.compile(r'''
(?P<indent>\s*)
if\s+os.path.isfile\(\s*
    (?P<string>f"[^"]+")
\s*\):
''', re.MULTILINE | re.VERBOSE)

OPEN_FILE = re.compile(r'''
open\(\s*
    (?P<string>"[^"]+")
    (?P<params>.+?)
    \s*,\s*(?:"r"|'r')\s*
\s*\)\s*
''', re.MULTILINE | re.VERBOSE)


@fn.curry(3)
def transform_with_regex(regex, transform, data, debug=False):
    fragments = []

    if debug:
        _transform = transform

        def transform(inpt, **kwargs):
            out = _transform(inpt, **kwargs)
            print('FROM:')
            print(inpt)
            print('TO:')
            print(out)
            print()
            return out

    idx = 0
    for m in regex.finditer(data):
        i, j = m.span()
        fragments.append(data[idx:i])
        source = data[i:j]
        fragments.append(transform(source, **m.groupdict()))
        idx = j
    fragments.append(data[idx:])

    return "".join(fragments)


def to_fstring(string, params):
    last_kind = tokenize.STRING
    parts = [string[1:-1]]
    token_list = list(tokens(params))
    for i, tk in enumerate(token_list):
        if tk.type in CTRL_TOKENS:
            continue
        if tk.type == tokenize.STRING:
            parts.append(tk.string[1:-1])
        elif tk.string == '+' and last_kind == tokenize.STRING:
            parts.append('{')
        elif tk.string == '+':
            parts.append('}')
            if token_list[i + 1].type != tokenize.STRING:
                parts.append('{')
        else:
            parts.append(tk.string)
        last_kind = tk.type

    out = ''.join(parts)
    return f'f"{out}"'


@transform_with_regex(OPEN_FILE)
def replace_open_file(fragment, string, params):
    out = to_fstring(string, params)
    return f'open({out}, "r")'


@transform_with_regex(IS_FILE)
def replace_is_file(fragment, indent, string, params):
    out = to_fstring(string, params)
    return f'{indent}if os.path.isfile({out}):'


@transform_with_regex(PATH)
def replace_path(fragment, indent, string):
    return f"{indent}path = {string}\n{indent}if os.path.isfile(path):"


def tokens(st):
    """Return an iterator of Python tokens"""
    return tokenize.tokenize(io.BytesIO(st.encode('utf8')).readline)


@fn.curry(2)
def apply_replacement(replacement, path):
    """
    Apply replacement to data
    """
    debug = os.environ.get('DEBUG', 'false').lower() in ('1', 'true')

    if not callable(replacement):
        replacement = globals()['replace_' + replacement]

    for path in glob(path):
        print('Processing', path)
        with open(path) as fd:
            source = fd.read()
            transformed = replacement(source, debug=debug)

        if source != transformed:
            # Check if it is valid Python source
            compile(transformed, 'string', 'exec')
            print('writing to', path)
            if not debug:
                with open(path + '.bak', 'w') as fd:
                    fd.write(source)
                with open(path, 'w') as fd:
                    fd.write(transformed)
        else:
            print('no changes detected')
        print()


def main():
    *_, transform, path = sys.argv
    apply_replacement(transform, path)


if __name__ == '__main__':
    main()
