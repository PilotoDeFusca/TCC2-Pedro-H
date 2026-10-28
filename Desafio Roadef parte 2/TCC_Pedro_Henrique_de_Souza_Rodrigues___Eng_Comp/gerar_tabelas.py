#!/usr/bin/env python3
"""
gerar_tabelas.py

Gera as tabelas de resultados do TCC a partir de ../figures/summary.csv
(produzido por compare_solutions.py):

    capitulos/tab_resultados.tex   uma linha por instancia
    capitulos/tab_conjuntos.tex    media/mediana por conjunto

Uso (a partir desta pasta):  python gerar_tabelas.py
"""
import csv
import os
import statistics as st
from pathlib import Path

here = Path(__file__).resolve().parent
rows = [r for r in csv.DictReader(open(here.parent / "figures" / "summary.csv"))
        if r["trivial_max"]]          # so instancias com solucao trivial


def num(x, d=3):
    return f"{float(x):.{d}f}".replace(".", "{,}")


def pct(x):
    v = float(x)
    return ("$-$" if v < 0 else "") + f"{abs(v):.1f}".replace(".", "{,}") + "\\%"


def write(path, text):
    tmp = str(path) + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    os.replace(tmp, path)


# ---- tabela por instancia ----------------------------------------------------
lines = [
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{MLU do pior período (avaliador oficial) e tempo de execução por instância.}",
    "\\label{tab:resultados}",
    "\\footnotesize",
    "\\begin{tabular}{lrrrrrrr}",
    "\\toprule",
    " & \\multicolumn{3}{c}{MLU do pior período} & \\multicolumn{2}{c}{Variação vs. trivial}"
    " & \\multicolumn{2}{c}{Tempo (s)} \\\\",
    "\\cmidrule(lr){2-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}",
    "Instância & Trivial & AMNS & VND & AMNS & VND & AMNS & VND \\\\",
    "\\midrule",
]
prev = None
for r in rows:
    s = r["instance"].split("-")[0]
    if prev and s != prev:
        lines.append("\\midrule")
    prev = s
    lines.append(" & ".join([
        r["instance"], num(r["trivial_max"]), num(r["amns_max"]), num(r["vnd_max"]),
        pct(r["amns_vs_trivial_pct"]), pct(r["vnd_vs_trivial_pct"]),
        num(r["amns_s"], 1), num(r["vnd_s"], 1)]) + " \\\\")
lines += ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
write(here / "capitulos" / "tab_resultados.tex", "\n".join(lines))

# ---- tabela por conjunto -----------------------------------------------------
lines = [
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Variação da MLU do pior período em relação à solução trivial e tempo total, por conjunto.}",
    "\\label{tab:conjuntos}",
    "\\begin{tabular}{lrrrrrrr}",
    "\\toprule",
    " & & \\multicolumn{2}{c}{Média} & \\multicolumn{2}{c}{Mediana} & \\multicolumn{2}{c}{Tempo total (s)} \\\\",
    "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}",
    "Conjunto & Instâncias & AMNS & VND & AMNS & VND & AMNS & VND \\\\",
    "\\midrule",
]
for g in ["setA", "setB", "setC", "setX", "todos"]:
    sel = [r for r in rows if g == "todos" or r["instance"].startswith(g)]
    a = [float(r["amns_vs_trivial_pct"]) for r in sel]
    v = [float(r["vnd_vs_trivial_pct"]) for r in sel]
    ta = sum(float(r["amns_s"]) for r in sel)
    tv = sum(float(r["vnd_s"]) for r in sel)
    name = "Todos" if g == "todos" else g.replace("set", "")
    if g == "todos":
        lines.append("\\midrule")
    lines.append(" & ".join([name, str(len(sel)), pct(st.mean(a)), pct(st.mean(v)),
                             pct(st.median(a)), pct(st.median(v)),
                             num(ta, 0), num(tv, 0)]) + " \\\\")
lines += ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
write(here / "capitulos" / "tab_conjuntos.tex", "\n".join(lines))
print("ok:", len(rows), "instancias")
