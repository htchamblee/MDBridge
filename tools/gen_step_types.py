"""Write the step-type schemas of MDATK's Project Planner as C++ data.

    python3 -I -B tools/gen_step_types.py MDATK_CHECKOUT

Imports MDATK (a checkout at v0.0.4.1; nothing is written into it) and dumps every registered
step type (mdatk.workflows.registry.TOOLS, in registration order) with its static ports, sections
and options, plus CATEGORIES and FILE_KINDS:

* src/core/step_types.inc                 the C++ data, included by src/core/step_types.cpp
* tests/data/step-types/step_types.json   the same as JSON, for tests/core/test_step_types.cpp

Values keep their Python type (tagged JSON, see tests/parity/make_step_type_cases.py); tuples
become arrays. Both outputs are deterministic; the build needs no Python.
"""

from __future__ import annotations

import datetime
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INC = os.path.join(ROOT, "src", "core", "step_types.inc")
JSON_OUT = os.path.join(ROOT, "tests", "data", "step-types", "step_types.json")
WIDTH = 100


def enc(v):
    """Tagged JSON of a plain value (see tests/parity/make_step_type_cases.py)."""
    if v is None:
        return {"t": "none"}
    if isinstance(v, bool):
        return {"t": "bool", "v": v}
    if isinstance(v, int):
        return {"t": "int", "v": str(v)}
    if isinstance(v, float):
        return {"t": "float", "hex": v.hex()}
    if isinstance(v, str):
        return {"t": "str", "v": v}
    if isinstance(v, (list, tuple)):
        return {"t": "array", "v": [enc(x) for x in v]}
    if isinstance(v, dict):
        return {"t": "table", "v": [[k, enc(x)] for k, x in v.items()]}
    for kind in (datetime.datetime, datetime.date, datetime.time):  # datetime before date
        if isinstance(v, kind):
            return {"t": kind.__name__, "iso": v.isoformat()}
    raise TypeError(f"unsupported value {v!r}")


def number(v):
    if v is None:
        return None
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise TypeError(f"not a number: {v!r}")
    return float(v)


def port(p):
    return {"name": p.name, "label": p.label, "kinds": list(p.kinds), "optional": bool(p.optional),
            "multiple": bool(p.multiple), "tooltip": p.tooltip}


def option(o):
    for c in o.choices:
        if not isinstance(c[0], str) or not isinstance(c[1], str):
            raise TypeError(f"{o.key}: choice {c!r} is not (str, str)")
    return {"key": o.key, "label": o.label, "kind": o.kind, "default": enc(o.default),
            "tooltip": o.tooltip, "choices": [[c[0], c[1]] for c in o.choices],
            "minimum": number(o.minimum), "maximum": number(o.maximum), "step": number(o.step),
            "digits": int(o.digits), "optional": bool(o.optional), "placeholder": o.placeholder,
            "unit": o.unit, "requires": [[k, enc(v)] for k, v in o.requires.items()]}


def section(s):
    return {"key": s.key, "title": s.title, "tooltip": s.tooltip, "repeat": bool(s.repeat),
            "item_label": s.item_label, "options": [option(o) for o in s.options],
            "templates": [{"id": t[0], "label": t[1], "values": enc(dict(t[2]))}
                          for t in s.templates],
            "default_items": list(s.default_items),
            "requires": [[k, enc(v)] for k, v in s.requires.items()]}


def tool(t):
    return {"name": t.name, "label": t.label, "category": t.category, "tag": t.tag,
            "description": t.description, "generates": bool(t.generates),
            "path_settings": list(t.path_settings), "inputs": [port(p) for p in t.inputs],
            "outputs": [port(p) for p in t.outputs], "sections": [section(s) for s in t.sections]}


def dump(registry, schema):
    return {"categories": list(registry.CATEGORIES),
            "file_kinds": [[k, label, list(exts)]
                           for k, (label, exts) in schema.FILE_KINDS.items()],
            "tools": [tool(t) for t in registry.TOOLS.values()]}


# ---------------------------------------------------------------------------------------- C++

def cpp_pieces(text: str) -> list[str]:
    """The escaped characters of a C++ string literal (non-ASCII as octal UTF-8 bytes)."""
    out = []
    for ch in text:
        if ch == "\\":
            out.append("\\\\")
        elif ch == '"':
            out.append('\\"')
        elif ch == "\n":
            out.append("\\n")
        elif ch == "\t":
            out.append("\\t")
        elif " " <= ch <= "~":
            out.append(ch)
        else:
            out.append("".join(f"\\{b:03o}" for b in ch.encode("utf-8")))
    return out


def cpp_string(text: str, column: int, indent: int) -> str:
    """A string literal that starts at `column`, split into adjacent literals to stay in WIDTH."""
    pieces = cpp_pieces(text)
    room = WIDTH - column - 3  # quotes and a comma
    if sum(len(p) for p in pieces) <= room:
        return '"' + "".join(pieces) + '"'
    words, word = [], []
    for p in pieces:  # words end after a space
        word.append(p)
        if p == " ":
            words.append(word)
            word = []
    if word:
        words.append(word)
    lines, current = [], ""
    width = WIDTH - indent - 3
    for w in words:
        for p in (w if sum(len(p) for p in w) > width else ["".join(w)]):
            if len(current) + len(p) > width and current:
                lines.append(current)
                current = ""
            current += p
    lines.append(current)
    pad = " " * indent
    return "\n" + "\n".join(f'{pad}"{line}"' for line in lines)


def cpp_double(x: float) -> str:
    text = repr(float(x))
    if text in ("inf", "-inf", "nan"):
        raise ValueError(f"no literal for {text}")
    return text


def cpp_value(v, column: int = 0, indent: int = 0) -> str:
    if v is None:
        return "Value()"
    if isinstance(v, bool):
        return "Value(true)" if v else "Value(false)"
    if isinstance(v, int):
        if not -2**63 <= v < 2**63:
            raise ValueError(f"integer beyond 64 bits: {v}")
        return f"Value(std::int64_t{{{v}}})"
    if isinstance(v, float):
        return f"Value({cpp_double(v)})"
    if isinstance(v, str):
        return f"Value({cpp_string(v, column + 6, indent)})"
    if isinstance(v, (list, tuple)):
        items = [cpp_value(x, 0, indent + 4) for x in v]
        one = "Value(Array{" + ", ".join(items) + "})"
        if column + len(one) + 2 <= WIDTH:
            return one
        pad = " " * (indent + 4)
        return "Value(Array{\n" + "".join(f"{pad}{x},\n" for x in items) + " " * indent + "})"
    raise TypeError(f"unsupported value {v!r}")


class Emitter:
    def __init__(self):
        self.lines: list[str] = []

    def add(self, indent: int, text: str):
        for k, line in enumerate(text.split("\n")):
            self.lines.append(((" " * indent + line) if k == 0 else line).rstrip())

    def field(self, indent: int, name: str, text: str):
        """`.name = "text",` with the literal wrapped below when long."""
        head = f".{name} = "
        literal = cpp_string(text, indent + len(head), indent + 4)
        self.add(indent, head + literal + ",")

    def value_field(self, indent: int, name: str, v):
        head = f".{name} = "
        self.add(indent, head + cpp_value(v, indent + len(head), indent + 4) + ",")


def kind_enum(kind: str) -> str:
    return "OptionKind::" + kind.capitalize()


def emit_requirements(e: Emitter, indent: int, requires):
    if not requires:
        return
    items = [f"{{{cpp_string(k, 0, 0)}, {cpp_value(v)}}}" for k, v in requires.items()]
    one = ".requirements = {" + ", ".join(items) + "},"
    if indent + len(one) <= WIDTH:
        e.add(indent, one)
        return
    e.add(indent, ".requirements = {")
    for k, v in requires.items():
        key = cpp_string(k, 0, 0)
        e.add(indent + 4, f"{{{key}, {cpp_value(v, indent + 4 + len(key) + 5, indent + 4)}}},")
    e.add(indent, "},")


def emit_option(e: Emitter, indent: int, o):
    e.add(indent, "Option{")
    i = indent + 4
    e.field(i, "key", o.key)
    e.field(i, "label", o.label)
    e.add(i, f".kind = {kind_enum(o.kind)},")
    e.value_field(i, "defaultValue", o.default)
    if o.tooltip:
        e.field(i, "tooltip", o.tooltip)
    if o.choices:
        e.add(i, ".choices = {")
        for value, label in o.choices:
            one = f"{{{cpp_string(value, 0, 0)}, {cpp_string(label, 0, 0)}}},"
            if i + 4 + len(one) > WIDTH:
                label = cpp_string(label, i + 5, i + 5).lstrip("\n").lstrip()
                one = f"{{{cpp_string(value, 0, 0)},\n{' ' * (i + 5)}{label}}},"
            e.add(i + 4, one)
        e.add(i, "},")
    for name in ("minimum", "maximum", "step"):
        if getattr(o, name) is not None:
            e.add(i, f".{name} = {cpp_double(number(getattr(o, name)))},")
    if o.digits:
        e.add(i, f".digits = {int(o.digits)},")
    if o.optional:
        e.add(i, ".optional = true,")
    if o.placeholder:
        e.field(i, "placeholder", o.placeholder)
    if o.unit:
        e.field(i, "unit", o.unit)
    emit_requirements(e, i, o.requires)
    e.add(indent, "},")


def emit_strings(e: Emitter, indent: int, name: str, items):
    if not items:
        return
    one = f".{name} = {{" + ", ".join(cpp_string(x, 0, 0) for x in items) + "},"
    if indent + len(one) <= WIDTH:
        e.add(indent, one)
        return
    e.add(indent, f".{name} = {{")
    for x in items:
        e.add(indent + 4, cpp_string(x, 0, 0) + ",")
    e.add(indent, "},")


def emit_section(e: Emitter, indent: int, s):
    e.add(indent, "Section{")
    i = indent + 4
    e.field(i, "key", s.key)
    e.field(i, "title", s.title)
    if s.tooltip:
        e.field(i, "tooltip", s.tooltip)
    if s.repeat:
        e.add(i, ".repeat = true,")
    if s.item_label != "name":
        e.field(i, "itemLabel", s.item_label)
    e.add(i, ".options = {")
    for o in s.options:
        emit_option(e, i + 4, o)
    e.add(i, "},")
    if s.templates:
        e.add(i, ".templates = {")
        for tid, label, values in s.templates:
            head = f"RecordTemplate{{{cpp_string(tid, 0, 0)}, {cpp_string(label, 0, 0)}, Table{{"
            e.add(i + 4, head)
            for k, v in values.items():
                key = cpp_string(k, 0, 0)
                e.add(i + 8, f"{{{key}, {cpp_value(v, i + 8 + len(key) + 3, i + 12)}}},")
            e.add(i + 4, "}},")
        e.add(i, "},")
    emit_strings(e, i, "defaultItems", list(s.default_items))
    emit_requirements(e, i, s.requires)
    e.add(indent, "},")


def emit_port(e: Emitter, indent: int, p):
    e.add(indent, "Port{")
    i = indent + 4
    e.field(i, "name", p.name)
    e.field(i, "label", p.label)
    emit_strings(e, i, "kinds", list(p.kinds))
    if p.optional:
        e.add(i, ".optional = true,")
    if p.multiple:
        e.add(i, ".multiple = true,")
    if p.tooltip:
        e.field(i, "tooltip", p.tooltip)
    e.add(indent, "},")


def function_name(name: str) -> str:
    parts = [p for p in name.replace("-", " ").replace("_", " ").split() if p]
    return "stepType" + "".join(p[:1].upper() + p[1:] for p in parts)


def emit_tool(e: Emitter, t):
    e.add(0, f"StepType {function_name(t.name)}() {{")
    e.add(4, "return StepType{")
    i = 8
    for name in ("name", "label", "category", "tag", "description"):
        e.field(i, name, getattr(t, name))
    if not t.generates:
        e.add(i, ".generates = false,")
    emit_strings(e, i, "pathSettings", list(t.path_settings))
    for name in ("inputs", "outputs"):
        ports = getattr(t, name)
        if ports:
            e.add(i, f".{name} = {{")
            for p in ports:
                emit_port(e, i + 4, p)
            e.add(i, "},")
    e.add(i, ".sections = {")
    for s in t.sections:
        emit_section(e, i + 4, s)
    e.add(i, "},")
    e.add(4, "};")
    e.add(0, "}")
    e.add(0, "")


def cpp(registry, schema) -> str:
    e = Emitter()
    e.add(0, "// Step types of MDATK v0.0.4.1's Project Planner (mdatk/workflows): registry.TOOLS,")
    e.add(0, "// CATEGORIES and schema.FILE_KINDS. Generated by tools/gen_step_types.py; do not edit.")
    e.add(0, "// Included by src/core/step_types.cpp inside its anonymous namespace.")
    e.add(0, "")
    for t in registry.TOOLS.values():
        emit_tool(e, t)
    e.add(0, "std::vector<StepType> builtinStepTypes() {")
    e.add(4, "std::vector<StepType> out;")
    e.add(4, f"out.reserve({len(registry.TOOLS)});")
    for t in registry.TOOLS.values():
        e.add(4, f"out.push_back({function_name(t.name)}());")
    e.add(4, "return out;")
    e.add(0, "}")
    e.add(0, "")
    e.add(0, "std::vector<std::string> builtinCategories() {")
    e.add(4, "return {")
    for c in registry.CATEGORIES:
        e.add(8, cpp_string(c, 0, 0) + ",")
    e.add(4, "};")
    e.add(0, "}")
    e.add(0, "")
    e.add(0, "std::vector<FileKindInfo> builtinFileKinds() {")
    e.add(4, "return {")
    for k, (label, exts) in schema.FILE_KINDS.items():
        e.add(8, f"{{{cpp_string(k, 0, 0)}, {cpp_string(label, 0, 0)}, "
                 f"{{{', '.join(cpp_string(x, 0, 0) for x in exts)}}}}},")
    e.add(4, "};")
    e.add(0, "}")
    return "\n".join(e.lines) + "\n"


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: gen_step_types.py MDATK_CHECKOUT", file=sys.stderr)
        return 2
    sys.dont_write_bytecode = True
    sys.path.insert(0, os.path.abspath(sys.argv[1]))
    from mdatk.workflows import registry, schema

    os.makedirs(os.path.dirname(JSON_OUT), exist_ok=True)
    with open(JSON_OUT, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(dump(registry, schema), fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    with open(INC, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(cpp(registry, schema))
    for path in (INC, JSON_OUT):
        print(f"wrote {os.path.relpath(path, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
