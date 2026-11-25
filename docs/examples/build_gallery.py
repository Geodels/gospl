"""
Build the examples gallery of the documentation from the goSPL-examples
repository (https://github.com/Geodels/goSPL-examples).

Called by ``docs/conf.py`` at the start of every build. It finds a copy of
goSPL-examples (``GOSPL_EXAMPLES_DIR``, or a shallow clone of the branch
that matches the docs being built), and for every example folder of
``Local-examples`` and ``Global-examples`` writes, under
``docs/examples/generated/`` (not tracked by git):

- ``<name>/index.md``: the example's ``README.md`` (or a short page when it
  has none) followed by links to its notebooks;
- ``<name>/view_Results.ipynb`` (and the other notebooks): copied as they
  are, rendered by myst-nb without being executed, so the figures are the
  ones stored in the notebook;
- ``<name>/thumb.png``: the card image of the gallery, the example's figure
  of the goSPL model-description paper when it has one
  (``docs/_static/examples/paper_figures.json``, shown at the top of its page
  with the paper's caption), else the first figure of ``view_Results.ipynb``;

and ``docs/examples/generated/gallery.rst``: the grid of cards included by
``docs/examples/index.rst``. Nothing is written into goSPL-examples.
"""

from __future__ import annotations

import base64
import json
import os
import re
import shutil
import subprocess

REPO = "https://github.com/Geodels/goSPL-examples.git"
HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(HERE, "generated")
SRC = os.path.join(os.path.dirname(HERE), "_examples_src")
GROUPS = [("Local-examples", "Regional examples"), ("Global-examples", "Global examples")]
NOTEBOOKS = ["view_Results.ipynb", "model_setup.ipynb", "build_inputs.ipynb"]
PAPER = os.path.join(os.path.dirname(HERE), "_static", "examples")
PAPER_NOTE = ("*Figure of the goSPL model-description paper "
              "(Geoscientific Model Development, in review).*")


def _paper_figure(page, fig):
    """Insert the paper figure and its caption below the page title."""
    block = ("\n```{figure} /_static/examples/%s\n:width: 100%%\n:alt: Figure of the example in the goSPL paper\n\n%s\n\n%s\n```\n"
             % (fig["image"], fig["caption"], PAPER_NOTE))
    m = re.search(r"^#\s+.+$", page, re.M)
    if m is None:
        return block + page
    return page[:m.end()] + "\n" + block + page[m.end():]


def _branch():
    """goSPL-examples branch matching the docs being built (GOSPL_EXAMPLES_BRANCH
    overrides): the feature branch for feature docs, dev for the development
    docs, master for a release."""
    if os.environ.get("GOSPL_EXAMPLES_BRANCH"):
        return os.environ["GOSPL_EXAMPLES_BRANCH"]
    rtd = os.environ.get("READTHEDOCS_VERSION", "")
    if rtd.startswith("feat"):
        return "feat/gravity-carbonates-karst"
    if rtd in ("dev", "latest"):
        return "dev"
    return "master"


def _source():
    """Path of a goSPL-examples checkout, or None (no gallery)."""
    local = os.environ.get("GOSPL_EXAMPLES_DIR")
    if local and os.path.isdir(local):
        return local
    if os.path.isdir(os.path.join(SRC, "Local-examples")):
        return SRC
    try:
        subprocess.run(["git", "clone", "--depth", "1", "--branch", _branch(), REPO, SRC],
                       check=True, capture_output=True, timeout=600)
        return SRC
    except Exception as exc:                       # offline, or branch missing
        print("[examples] goSPL-examples not available (%s): no gallery" % exc)
        return None


def _title(readme, name):
    if readme:
        m = re.search(r"^#\s+(.+)$", readme, re.M)
        if m:
            return m.group(1).strip()
    return name.replace("_", " ").capitalize()


def _summary(readme):
    """First paragraph of the README (plain text, for the card)."""
    if not readme:
        return ""
    body = re.split(r"^#\s+.+$", readme, maxsplit=1, flags=re.M)[-1]
    for para in re.split(r"\n\s*\n", body):
        p = para.strip()
        if p and not p.startswith(("#", "!", "|", "```", "<", "[![")):
            p = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", p)
            p = re.sub(r"[`*_]", "", " ".join(p.split()))
            return p if len(p) < 170 else p[:167].rsplit(" ", 1)[0] + "..."
    return ""


def _thumb(nb_path, out_png):
    """Save the first PNG figure of a notebook; False when it has none."""
    try:
        nb = json.load(open(nb_path))
    except Exception:
        return False
    for cell in nb.get("cells", []):
        for o in cell.get("outputs", []):
            data = o.get("data", {})
            if "image/png" in data:
                png = data["image/png"]
                png = "".join(png) if isinstance(png, list) else png
                open(out_png, "wb").write(base64.b64decode(png))
                return True
    return False


def _executed(nb_path):
    try:
        nb = json.load(open(nb_path))
    except Exception:
        return False
    return any(c.get("outputs") for c in nb.get("cells", []) if c.get("cell_type") == "code")


def build():
    src = _source()
    if os.path.isdir(GEN):
        shutil.rmtree(GEN)
    os.makedirs(GEN)
    cards = []
    if src is None:
        open(os.path.join(GEN, "gallery.rst"), "w").write(
            ".. note::\n\n    The examples gallery is built from the goSPL-examples "
            "repository, which was not available for this build.\n")
        return
    branch = _branch()
    pfile = os.path.join(PAPER, "paper_figures.json")
    paper = json.load(open(pfile)) if os.path.exists(pfile) else {}
    # one-line descriptions of the examples table of the repository README:
    # the card text and the page of an example that has no README of its own
    index = {}
    top = os.path.join(src, "README.md")
    if os.path.exists(top):
        for m in re.finditer(r"^\|\s*\[`([^`]+)`\]\([^)]*\)\s*\|\s*(.+?)\s*\|\s*$",
                             open(top).read(), re.M):
            index[m.group(1)] = re.sub(r"\s*Needs goSPL with .*$", "", m.group(2)).strip()
    for group, gtitle in GROUPS:
        gdir = os.path.join(src, group)
        if not os.path.isdir(gdir):
            continue
        items = []
        for name in sorted(os.listdir(gdir)):
            edir = os.path.join(gdir, name)
            if not os.path.isdir(edir) or name.startswith((".", "_")):
                continue
            readme_p = os.path.join(edir, "README.md")
            readme = open(readme_p).read() if os.path.exists(readme_p) else None
            out = os.path.join(GEN, name)
            os.makedirs(out)
            nbs = [n for n in NOTEBOOKS if os.path.exists(os.path.join(edir, n))]
            for n in nbs:
                shutil.copy(os.path.join(edir, n), os.path.join(out, n))
            title = _title(readme, name)
            url = "https://github.com/Geodels/goSPL-examples/tree/%s/%s/%s" % (branch, group, name)
            page = readme if readme else "# %s\n\n%s\n" % (
                title, index.get(name, "This example does not have a description yet."))
            # links in the README to other files of the example point to GitHub
            page = re.sub(r"\]\((?!https?://|#)([^)]+)\)", lambda m: "](%s/%s)" % (url, m.group(1)), page)
            if name in paper:
                page = _paper_figure(page, paper[name])
            page += "\n\n---\n\nSource, inputs and scripts: [%s on GitHub](%s).\n" % (name, url)
            done = [n for n in nbs if _executed(os.path.join(out, n))]
            if nbs:
                page += "\n```{toctree}\n:maxdepth: 1\n:caption: Notebooks\n\n"
                page += "".join("%s\n" % n[:-6] for n in nbs) + "```\n"
            open(os.path.join(out, "index.md"), "w").write(page)
            if name in paper:
                shutil.copy(os.path.join(PAPER, paper[name]["image"]), os.path.join(out, "thumb.png"))
                has_thumb = True
            else:
                has_thumb = "view_Results.ipynb" in done and _thumb(
                    os.path.join(out, "view_Results.ipynb"), os.path.join(out, "thumb.png"))
            summ = _summary(readme) or index.get(name, "")
            if len(summ) > 170:
                summ = summ[:167].rsplit(" ", 1)[0] + "..."
            summ = re.sub(r"[`*_]", "", re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", summ))
            items.append((name, title, summ, has_thumb))
        cards.append((gtitle, items))

    L = []
    for gtitle, items in cards:
        L += [gtitle, "-" * len(gtitle), "", ".. grid:: 1 2 3 3", "    :gutter: 3", ""]
        for name, title, summ, thumb in items:
            L += ["    .. grid-item-card:: %s" % title, "        :link: generated/%s/index" % name,
                  "        :link-type: doc", "        :class-card: gospl-example-card"]
            if thumb:
                L += ["        :img-top: generated/%s/thumb.png" % name]
            L += ["", "        %s" % (summ or name.replace("_", " ")), ""]
        L += [""]
    L += [".. toctree::", "    :hidden:", ""]
    L += ["    generated/%s/index" % n for _, items in cards for n, *_ in items]
    open(os.path.join(GEN, "gallery.rst"), "w").write("\n".join(L) + "\n")
    print("[examples] gallery: %d examples from goSPL-examples (%s)"
          % (sum(len(i) for _, i in cards), src))


if __name__ == "__main__":
    build()
