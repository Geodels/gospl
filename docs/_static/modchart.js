/* goSPL module chooser: pick the processes you want, see what they need.

   Data: window.GOSPL_MODULES = {modules: [{id, label, block, page, summary,
   group, always, requires: [{id, why, enforced}], requires_any: [[{id, why}],
   ...], excludes: [{id, why, enforced}], recommends: [{id, why}], mesh}]}
   written into the page by docs/user_guide/build_reference.py.

   The core is selected as soon as anything is. A selected module is green
   when everything it needs is selected too, red while something is missing;
   a missing requirement is drawn as a dashed "+ add" block (click it, or
   its checkbox, to add it). The graph is laid out in columns (what is
   needed on the left), ordered to limit crossings, and the links are curves
   that pass between the blocks. */
(function () {
  "use strict";

  function el(tag, attrs, children) {
    var e = document.createElement(tag);
    if (attrs) {
      Object.keys(attrs).forEach(function (k) {
        if (k === "text") e.textContent = attrs[k];
        else if (k !== "html") e.setAttribute(k, attrs[k]);
      });
    }
    (children || []).forEach(function (c) { if (c) e.appendChild(c); });
    return e;
  }

  function svg(tag, attrs) {
    var e = document.createElementNS("http://www.w3.org/2000/svg", tag);
    Object.keys(attrs || {}).forEach(function (k) { e.setAttribute(k, attrs[k]); });
    return e;
  }

  function init(root) {
    var data = window.GOSPL_MODULES;
    if (!data) return;
    var mods = {}, coreId = null;
    data.modules.forEach(function (m) {
      mods[m.id] = m;
      if (m.always || m.id === "core") coreId = coreId || m.id;
    });
    var picked = {};
    var domain = "regional";            // regional (flat mesh) or global (sphere)

    // the core is selected as soon as anything is, with the climate rain
    // unless the orographic rain replaces it
    var rainId = mods.climate ? "climate" : null, oroId = mods.orography ? "orography" : null;
    function selected() {
      var s = {};
      Object.keys(picked).forEach(function (id) { s[id] = true; });
      if (Object.keys(s).length && coreId) {
        s[coreId] = true;
        if (rainId && !(oroId && s[oroId])) s[rainId] = true;
      }
      return s;
    }
    function automatic(id) {
      return id === coreId || (id === rainId && !picked[rainId] && Object.keys(picked).length > 0);
    }

    function fits(m) {
      return !m.domain || m.domain === "both" ||
        (m.domain === "flat" && domain === "regional") || (m.domain === "sphere" && domain === "global");
    }

    // ---- state of the selection: which requirements are missing, conflicts
    function analyse() {
      var sel = selected(), missing = {}, status = {}, conflicts = [], anyMissing = [];
      Object.keys(sel).forEach(function (id) {
        var m = mods[id], ok = true;
        (m.requires || []).forEach(function (r) {
          if (!mods[r.id] || sel[r.id]) return;
          ok = false;
          (missing[r.id] = missing[r.id] || []).push({by: id, why: r.why, enforced: r.enforced});
        });
        (m.requires_any || []).forEach(function (alts) {
          var alive = alts.filter(function (r) { return mods[r.id]; });
          if (!alive.length || alive.some(function (r) { return sel[r.id]; })) return;
          ok = false;
          anyMissing.push({by: id, alts: alive});
          alive.forEach(function (r) {
            (missing[r.id] = missing[r.id] || []).push({by: id, why: r.why, alt: true});
          });
        });
        if (!fits(m)) {
          ok = false;
          conflicts.push({a: id, b: null, why: "not available on a " + (domain === "global" ?
            "global (spherical)" : "regional (flat)") + " mesh"});
        }
        (m.excludes || []).forEach(function (x) {
          if (sel[x.id]) {
            ok = false;
            if (id < x.id || !(mods[x.id].excludes || []).some(function (y) { return y.id === id; }))
              conflicts.push({a: id, b: x.id, why: x.why});
          }
        });
        status[id] = ok ? "ok" : "bad";
      });
      return {sel: sel, missing: missing, status: status, conflicts: conflicts, anyMissing: anyMissing};
    }

    // ---- controls: one checkbox per module, by group (core is automatic)
    var groups = {}, order = [];
    data.modules.forEach(function (m) {
      if (m.id === coreId) return;
      var g = m.group || "Other";
      if (!groups[g]) { groups[g] = []; order.push(g); }
      groups[g].push(m);
    });
    var controls = el("div", {"class": "modchart-controls"});
    var dom = el("div", {"class": "modchart-domain", role: "radiogroup", "aria-label": "Kind of model"});
    [["regional", "Regional model", "flat mesh, metres to hundreds of km"],
     ["global", "Global model", "spherical mesh of the Earth"]].forEach(function (d) {
      var b = el("button", {type: "button", "class": "modchart-dom" + (d[0] === domain ? " on" : ""),
                            role: "radio", "aria-checked": d[0] === domain ? "true" : "false"},
        [el("strong", {text: d[1]}), el("span", {text: d[2]})]);
      b.addEventListener("click", function () {
        domain = d[0];
        dom.querySelectorAll(".modchart-dom").forEach(function (x) {
          x.classList.remove("on"); x.setAttribute("aria-checked", "false"); });
        b.classList.add("on"); b.setAttribute("aria-checked", "true");
        resources.setDomain(domain);
        render();
      });
      dom.appendChild(b);
    });
    order.forEach(function (g) {
      var box = el("fieldset", {"class": "modchart-group"}, [el("legend", {text: g})]);
      groups[g].forEach(function (m) {
        var cb = el("input", {type: "checkbox", id: "mc-" + m.id});
        cb.addEventListener("change", function () { toggle(m.id, cb.checked); });
        var lab = el("label", {"for": "mc-" + m.id, title: m.summary || ""}, [cb,
          el("span", {text: " " + m.label})]);
        box.appendChild(el("div", {"class": "modchart-item", id: "mci-" + m.id}, [lab]));
      });
      controls.appendChild(box);
    });
    var reset = el("button", {type: "button", "class": "modchart-reset", text: "Clear"});
    reset.addEventListener("click", function () {
      picked = {};
      controls.querySelectorAll("input").forEach(function (c) { c.checked = false; });
      render();
    });
    controls.appendChild(reset);

    var graph = el("div", {"class": "modchart-graph"});
    var report = el("div", {"class": "modchart-report"});
    var extras = el("div", {"class": "modchart-extras"});
    root.appendChild(dom);
    root.appendChild(controls);
    root.appendChild(graph);
    root.appendChild(report);
    root.appendChild(extras);
    var resources = makeResources();
    root.appendChild(resources.node);

    function toggle(id, on) {
      if (id === coreId) return;
      if (id === rainId && !on && !picked[rainId] && Object.keys(picked).length) {
        // the rain is automatic: tick the orographic rain to replace it
        var cb0 = root.querySelector("#mc-" + id);
        if (cb0) cb0.checked = true;
        return;
      }
      if (on) picked[id] = true; else delete picked[id];
      var cb = root.querySelector("#mc-" + id);
      if (cb) cb.checked = !!on;
      render();
    }

    function link(m) {
      if (!m.page) return el("strong", {text: m.label});
      var a = el("a", {href: (window.GOSPL_DOCROOT || "../") + m.page + ".html"});
      a.appendChild(el("strong", {text: m.label}));
      return a;
    }

    function render() {
      var st = analyse();
      data.modules.forEach(function (m) {
        var item = root.querySelector("#mci-" + m.id);
        if (!item) return;
        var box = item.querySelector("input");
        if (box) box.checked = !!picked[m.id] || (automatic(m.id) && !!st.sel[m.id]);
        item.classList.toggle("is-auto", automatic(m.id) && !!st.sel[m.id]);
        item.classList.toggle("is-missing", m.id in st.missing);
        item.classList.toggle("is-bad", st.status[m.id] === "bad");
        item.classList.toggle("is-ok", st.status[m.id] === "ok");
        item.classList.toggle("is-off", !fits(m));
        item.title = fits(m) ? (m.summary || "") : "Not available on a " + domain + " model";
      });
      drawGraph(st);
      writeReport(st);
      extras.innerHTML = "";
      writeExamples(st);
      writeParams(st);
      resources.update(st);
    }

    // ---- report under the graph
    function writeReport(st) {
      report.innerHTML = "";
      var ids = Object.keys(st.sel);
      if (!ids.length) {
        report.appendChild(el("p", {"class": "modchart-hint",
          text: "Tick the processes you want to use. A selected module turns green " +
                "when everything it needs is selected, red while something is missing; " +
                "click a dashed block to add it."}));
        return;
      }
      var miss = Object.keys(st.missing);
      if (miss.length) {
        report.appendChild(el("h4", {"class": "modchart-bad", text: "Still to add"}));
        var u = el("ul");
        miss.forEach(function (id) {
          var li = el("li", {}, [link(mods[id])]);
          var add = el("button", {type: "button", "class": "modchart-add", text: "+ add"});
          add.addEventListener("click", function () { toggle(id, true); });
          li.appendChild(add);
          var why = el("ul", {"class": "modchart-why"});
          st.missing[id].forEach(function (r) {
            why.appendChild(el("li", {text: (r.alt ? "one of the options for " : "needed by ") +
              mods[r.by].label + ": " + r.why}));
          });
          li.appendChild(why);
          u.appendChild(li);
        });
        report.appendChild(u);
      }
      if (st.conflicts.length) {
        report.appendChild(el("h4", {"class": "modchart-bad", text: "Not compatible"}));
        var c = el("ul");
        st.conflicts.forEach(function (x) {
          c.appendChild(el("li", {text: x.b ? mods[x.a].label + " and " + mods[x.b].label + ": " + x.why
                                             : mods[x.a].label + ": " + x.why}));
        });
        report.appendChild(c);
      }
      report.appendChild(el("h4", {text: "Blocks to set in the input file"}));
      var ul = el("ul");
      ids.sort(function (a, b) { return a === coreId ? -1 : b === coreId ? 1 : 0; });
      ids.forEach(function (id) {
        var m = mods[id], li = el("li", {}, [link(m)]);
        if (m.block) li.appendChild(el("code", {"class": "modchart-block", text: m.block}));
        li.appendChild(el("span", {"class": "modchart-tag " + (st.status[id] === "ok" ? "ok" : "bad"),
          text: st.status[id] === "ok" ? "complete" : "incomplete"}));
        if (m.mesh && m.mesh !== "both")
          li.appendChild(el("div", {"class": "modchart-why", text: "Mesh: " + m.mesh}));
        ul.appendChild(li);
      });
      report.appendChild(ul);
      var recs = [];
      ids.forEach(function (id) {
        (mods[id].recommends || []).forEach(function (r) {
          if (!st.sel[r.id] && mods[r.id] && !(r.id in st.missing))
            recs.push(mods[r.id].label + " (with " + mods[id].label + "): " + r.why);
        });
      });
      if (recs.length) {
        report.appendChild(el("h4", {text: "Often used with"}));
        var r = el("ul");
        recs.forEach(function (t) { r.appendChild(el("li", {text: t})); });
        report.appendChild(r);
      }
    }

    // ---- layered layout (Sugiyama style): columns by depth of requirement,
    // dummy points where a link crosses a column, barycentre ordering
    function drawGraph(st) {
      graph.innerHTML = "";
      var nodes = Object.keys(st.sel).concat(Object.keys(st.missing).filter(function (id) {
        return !st.sel[id];
      }));
      if (!nodes.length) return;
      var shown = {};
      nodes.forEach(function (id) { shown[id] = true; });

      // links: from a module to what it needs (only from selected modules)
      var edges = [];
      Object.keys(st.sel).forEach(function (id) {
        var m = mods[id];
        (m.requires || []).forEach(function (r) {
          if (shown[r.id]) edges.push({from: id, to: r.id, ok: !!st.sel[r.id], alt: false});
        });
        (m.requires_any || []).forEach(function (alts) {
          alts.forEach(function (r) {
            if (shown[r.id]) edges.push({from: id, to: r.id, ok: !!st.sel[r.id], alt: true});
          });
        });
      });

      // depth: longest chain of requirements below a node
      var depth = {};
      function d(id, stack) {
        if (id in depth) return depth[id];
        stack = stack || {};
        if (stack[id]) return 0;
        stack[id] = true;
        var k = 0;
        edges.forEach(function (e) { if (e.from === id) k = Math.max(k, d(e.to, stack) + 1); });
        delete stack[id];
        depth[id] = k;
        return k;
      }
      nodes.forEach(function (id) { d(id); });
      var ncol = 1 + Math.max.apply(null, nodes.map(function (id) { return depth[id]; }));
      var cols = [];
      for (var c = 0; c < ncol; c++) cols.push([]);
      nodes.forEach(function (id) { cols[depth[id]].push({id: id}); });

      // dummy points on links that skip columns
      var chains = [];
      edges.forEach(function (e, i) {
        var a = depth[e.from], b = depth[e.to], chain = [e.from];
        for (var k = a - 1; k > b; k--) {
          var dummy = {id: "_d" + i + "_" + k, dummy: true};
          cols[k].push(dummy);
          chain.push(dummy.id);
        }
        chain.push(e.to);
        chains.push({edge: e, ids: chain});
      });

      // neighbours of every point, for the barycentre ordering
      var nbr = {};
      chains.forEach(function (ch) {
        for (var k = 0; k + 1 < ch.ids.length; k++) {
          (nbr[ch.ids[k]] = nbr[ch.ids[k]] || []).push(ch.ids[k + 1]);
          (nbr[ch.ids[k + 1]] = nbr[ch.ids[k + 1]] || []).push(ch.ids[k]);
        }
      });
      function index() {
        var idx = {};
        cols.forEach(function (col) { col.forEach(function (n, i) { idx[n.id] = i / Math.max(1, col.length - 1); }); });
        return idx;
      }
      for (var it = 0; it < 6; it++) {
        var idx = index();
        cols.forEach(function (col) {
          col.forEach(function (n) {
            var ns = nbr[n.id] || [];
            n.key = ns.length ? ns.reduce(function (s, x) { return s + idx[x]; }, 0) / ns.length : idx[n.id];
          });
          col.sort(function (p, q) { return p.key - q.key; });
        });
      }

      // geometry: required modules on the left, selected ones on the right
      var W = 200, H = 36, GX = 90, ROW = 54, DH = 10;
      var nrow = Math.max.apply(null, cols.map(function (col) { return col.length; }));
      var width = ncol * W + (ncol - 1) * GX + 24, height = Math.max(nrow, 1) * ROW + 24;
      var pos = {};
      cols.forEach(function (col, k) {
        var x = 12 + k * (W + GX);
        var top = 12 + (nrow - col.length) * ROW / 2;
        col.forEach(function (n, i) {
          pos[n.id] = {x: x, y: top + i * ROW + (ROW - H) / 2, dummy: !!n.dummy};
        });
      });

      var s = svg("svg", {viewBox: "0 0 " + width + " " + height, width: "100%",
                          "class": "modchart-svg", role: "img",
                          "aria-label": "Selected modules and what they need"});
      var defs = svg("defs");
      ["ok", "bad"].forEach(function (k) {
        var mk = svg("marker", {id: "mc-arrow-" + k, viewBox: "0 0 10 10", refX: 9, refY: 5,
                                markerWidth: 9, markerHeight: 9, markerUnits: "userSpaceOnUse",
                                orient: "auto"});
        mk.appendChild(svg("path", {d: "M0 0 L10 5 L0 10 z", fill: k === "ok" ? "#6b8f71" : "#c0392b"}));
        defs.appendChild(mk);
      });
      s.appendChild(defs);

      // links: horizontal-tangent cubic curves through the dummy points, so
      // they leave and enter the blocks sideways and run between them
      chains.forEach(function (ch) {
        // from the left side of the module to the right side of what it
        // needs, through the centre of each dummy slot (free rows between
        // the blocks of the columns it crosses), with horizontal tangents
        var pts = ch.ids.map(function (id, k) {
          var p = pos[id];
          if (p.dummy) return [p.x + W / 2, p.y + H / 2];
          return [k === 0 ? p.x : p.x + W, p.y + H / 2];
        });
        var dpath = "M" + pts[0][0] + " " + pts[0][1];
        for (var k = 1; k < pts.length; k++) {
          var a0 = pts[k - 1], b0 = pts[k], mx = (a0[0] + b0[0]) / 2;
          dpath += " C" + mx + " " + a0[1] + " " + mx + " " + b0[1] + " " + b0[0] + " " + b0[1];
        }
        var e = ch.edge, cls = "modchart-edge " + (e.ok ? "ok" : "bad") + (e.alt ? " alt" : "") +
          (e.core ? " core" : "");
        s.appendChild(svg("path", {d: dpath, "class": cls,
                                   "marker-end": "url(#mc-arrow-" + (e.ok ? "ok" : "bad") + ")"}));
      });

      nodes.forEach(function (id) {
        var p = pos[id], m = mods[id];
        var cls = st.sel[id] ? (st.status[id] === "ok" ? "ok" : "bad") : "missing";
        var g = svg("g", {"class": "modchart-node " + cls + (id === coreId ? " core" : ""),
                          tabindex: 0, role: "button"});
        g.appendChild(svg("rect", {x: p.x, y: p.y, width: W, height: H, rx: 7, ry: 7}));
        var lab = (cls === "missing" ? "+ " : "") + m.label;
        var t = svg("text", {x: p.x + W / 2, y: p.y + H / 2 + 4.5, "text-anchor": "middle"});
        t.textContent = lab.length > 27 ? lab.slice(0, 26) + "…" : lab;
        var tt = svg("title");
        tt.textContent = m.label + (m.summary ? ": " + m.summary : "") +
          (cls === "missing" ? " (click to add)" : id === coreId ? "" : " (click to remove)");
        g.appendChild(tt);
        g.appendChild(t);
        if (id !== coreId) {
          g.addEventListener("click", function () { toggle(id, !st.sel[id]); });
          g.addEventListener("keydown", function (ev) {
            if (ev.key === "Enter" || ev.key === " ") { ev.preventDefault(); toggle(id, !st.sel[id]); }
          });
        }
        s.appendChild(g);
      });
      graph.appendChild(s);
      var legend = el("div", {"class": "modchart-legend"});
      [["ok", "selected, complete"], ["bad", "selected, needs something"],
       ["missing", "needed, not selected (click to add)"]].forEach(function (k) {
        legend.appendChild(el("span", {"class": "modchart-key " + k[0]}));
        legend.appendChild(el("span", {text: " " + k[1] + "   "}));
      });
      legend.appendChild(el("span", {text: "Arrows point to what a module needs; dashed: one of several options."}));
      graph.appendChild(legend);
    }

    function fmt(x, unit) {
      if (x >= 1e6) return (x / 1e6).toPrecision(3) + " M" + (unit || "");
      if (x >= 1e3) return (x / 1e3).toPrecision(3) + " k" + (unit || "");
      return (x < 10 ? x.toPrecision(2) : Math.round(x)) + (unit ? " " + unit : "");
    }
    function years(y) {
      if (y >= 1e6) return +(y / 1e6).toPrecision(3) + " Myr";
      if (y >= 1e3) return +(y / 1e3).toPrecision(3) + " kyr";
      return y + " yr";
    }

    // ---- examples closest to the selection (same kind of model first)
    function writeExamples(st) {
      var ex = data.examples || [];
      var want = Object.keys(st.sel).filter(function (id) { return !automatic(id) || picked[id]; });
      var chosen = Object.keys(picked);
      if (!chosen.length || !ex.length) return;
      var best = {};
      ex.forEach(function (r) {
        if ((r.group === "global") !== (domain === "global")) return;
        var have = {};
        r.modules.forEach(function (i) { have[i] = true; });
        var hit = want.filter(function (i) { return have[i]; }).length;
        var score = hit / Math.max(1, want.length) - 0.01 * Math.max(0, r.modules.length - hit);
        if (!best[r.example] || score > best[r.example].score) best[r.example] = {r: r, score: score, hit: hit};
      });
      var list = Object.keys(best).map(function (k) { return best[k]; })
        .filter(function (b) { return b.hit > 1; })
        .sort(function (a, b) { return b.score - a.score; }).slice(0, 5);
      var box = el("div", {"class": "modchart-panel"});
      box.appendChild(el("h4", {text: "Closest examples (" + (domain === "global" ? "global" : "regional") + ")"}));
      if (!list.length) {
        box.appendChild(el("p", {"class": "modchart-hint",
          text: "No " + domain + " example uses these modules yet."}));
        extras.appendChild(box);
        return;
      }
      var ul = el("ul", {"class": "modchart-examples"});
      list.forEach(function (b) {
        var r = b.r, li = el("li");
        var a = el("a", {href: (window.GOSPL_DOCROOT || "../") + r.page + ".html"});
        a.appendChild(el("strong", {text: r.example.replace(/_/g, " ")}));
        li.appendChild(a);
        li.appendChild(el("code", {"class": "modchart-block", text: r.input}));
        var mesh = r.nodes ? fmt(r.nodes) + " nodes, spacing " +
          (r.spacing_m >= 1000 ? (r.spacing_m / 1000).toPrecision(2) + " km" : Math.round(r.spacing_m) + " m") +
          (r.extent_km ? ", " + r.extent_km[0] + " x " + r.extent_km[1] + " km" : ", whole Earth") : "";
        li.appendChild(el("div", {"class": "modchart-why", text: [mesh, years(r.duration_yr) + " in steps of " +
          years(r.dt_yr) + " (" + r.steps + " steps)", r.strat_yr ? "a layer every " + years(r.strat_yr) : ""]
          .filter(Boolean).join("; ")}));
        var missing = want.filter(function (i) { return r.modules.indexOf(i) < 0; })
          .map(function (i) { return mods[i].label; });
        li.appendChild(el("div", {"class": "modchart-why", text: "Uses " + b.hit + " of your " + want.length +
          " modules" + (missing.length ? "; not in it: " + missing.join(", ") : "")}));
        var use = el("button", {type: "button", "class": "modchart-add", text: "use its mesh and time"});
        use.addEventListener("click", function () { resources.fromExample(r); });
        li.appendChild(use);
        ul.appendChild(li);
      });
      box.appendChild(ul);
      extras.appendChild(box);
    }

    // ---- main parameters of the selected modules, with their usual values
    function writeParams(st) {
      var ids = Object.keys(st.sel).filter(function (id) { return mods[id].params; });
      if (!Object.keys(picked).length || !ids.length) return;
      var box = el("div", {"class": "modchart-panel"});
      box.appendChild(el("h4", {text: "Main parameters and usual values"}));
      ids.forEach(function (id) {
        var d = el("details", {"class": "modchart-params"});
        d.appendChild(el("summary", {text: mods[id].label + " (" + mods[id].params.length + ")"}));
        var t = el("table");
        t.appendChild(el("tr", {}, [el("th", {text: "key"}), el("th", {text: "unit"}),
          el("th", {text: "default"}), el("th", {text: "usual values"})]));
        mods[id].params.forEach(function (p) {
          t.appendChild(el("tr", {}, [el("td", {}, [el("code", {text: p.key})]), el("td", {text: p.unit}),
            el("td", {text: p.default}), el("td", {text: p.range})]));
        });
        d.appendChild(t);
        box.appendChild(d);
      });
      box.appendChild(el("p", {"class": "modchart-hint", html: ""}));
      var a = el("a", {href: (window.GOSPL_DOCROOT || "../") + "user_guide/reference/inputs.html"});
      a.textContent = "every key, with how it depends on the resolution";
      var p = box.lastChild;
      p.appendChild(document.createTextNode("All units and ranges: "));
      p.appendChild(a);
      extras.appendChild(box);
    }

    // ---- resources: nodes, steps, time and memory from the measured costs
    function makeResources() {
      // base step: 1.67 s at 1 rank for 47,672 nodes (cost page), of which
      // the marine deposition is ~64 %
      var NODE_S = 1.6664 / 47672, MARINE = 0.64, LIB_GB = 0.56, KB_NODE = 6.5, ROOT_KB = 1.8;
      var node = el("div", {"class": "modchart-panel modchart-res"});
      node.appendChild(el("h4", {text: "Your run: size and resources"}));
      var form = el("div", {"class": "modchart-form"}), out = el("div", {"class": "modchart-out"});
      var fields = {};
      function field(key, label, val, unit, cls) {
        var inp = el("input", {type: "number", min: "0", step: "any", value: val});
        inp.addEventListener("input", function () { update(); });
        fields[key] = inp;
        var w = el("label", {"class": "modchart-field " + (cls || "")}, [el("span", {text: label}), inp,
          el("em", {text: unit})]);
        form.appendChild(w);
        return w;
      }
      field("lx", "Domain width", 200, "km", "reg");
      field("ly", "Domain length", 100, "km", "reg");
      field("dx", "Mesh spacing", 500, "m", "reg");
      field("dxg", "Mesh spacing", 25, "km", "glo");
      field("dur", "Duration", 1, "Myr");
      field("dt", "Time step", 5, "kyr");
      field("strat", "Layer interval", 5, "kyr", "strat");
      field("ranks", "Processors", 0, "(0: suggested)");
      node.appendChild(form);
      node.appendChild(out);
      var last = null;

      function setDomain(d) {
        if (d === "global") { fields.dur.value = 10; fields.dt.value = 50; fields.strat.value = 100; }
        else { fields.dur.value = 1; fields.dt.value = 5; fields.strat.value = 5; }
        update();
      }
      function fromExample(r) {
        if (r.group === "global") fields.dxg.value = +(r.spacing_m / 1000).toPrecision(3);
        else if (r.extent_km) {
          fields.lx.value = r.extent_km[0]; fields.ly.value = r.extent_km[1];
          fields.dx.value = Math.round(r.spacing_m);
        }
        fields.dur.value = +(r.duration_yr / 1e6).toPrecision(4);
        fields.dt.value = +(r.dt_yr / 1e3).toPrecision(4);
        if (r.strat_yr) fields.strat.value = +(r.strat_yr / 1e3).toPrecision(4);
        update();
        node.scrollIntoView({behavior: "smooth", block: "start"});
      }
      function update(st) {
        if (st) last = st;
        st = last;
        form.classList.toggle("is-global", domain === "global");
        var sel = st ? st.sel : {};
        form.classList.toggle("has-strat", !!sel.strata);
        out.innerHTML = "";
        if (!st || !Object.keys(picked).length) {
          out.appendChild(el("p", {"class": "modchart-hint",
            text: "Tick at least one process to estimate the size of the run."}));
          return;
        }
        var v = function (k) { return Math.max(0, parseFloat(fields[k].value) || 0); };
        var nodes = domain === "global"
          ? 4 * Math.PI * Math.pow(6371, 2) / (0.866 * Math.pow(v("dxg"), 2))
          : v("lx") * v("ly") * 1e6 / (0.866 * Math.pow(Math.max(v("dx"), 1), 2));
        var steps = v("dt") > 0 ? Math.ceil(v("dur") * 1e3 / v("dt")) : 0;
        var factor = (sel.marine_deposition ? 1 : 1 - MARINE);
        Object.keys(sel).forEach(function (id) { factor += (mods[id].cost || 0) / 100; });
        var cpuStep = nodes * NODE_S * factor;                      // core-seconds per step
        var coreH = cpuStep * steps / 3600;
        var ranks = v("ranks") >= 1 ? Math.round(v("ranks")) : Math.max(1, Math.round(nodes / 50000));
        var eff = ranks === 1 ? 1 : ranks <= 128 ? 0.85 : 0.7;
        var wall = cpuStep * steps / (ranks * eff);
        var layers = sel.strata && v("strat") > 0 ? Math.ceil(v("dur") * 1e3 / v("strat")) : 0;
        var perLayer = 32 + (sel.dual ? 16 : 0) + (sel.carbonate ? 48 : 0) + (sel.lithologies ? 16 : 0) +
          (sel.provenance ? 8 * (3 + (sel.carbonate || sel.pelagic ? 1 : 0)) : 0) +
          (sel.duricrust ? 24 : 0);
        var strataGB = layers * nodes * perLayer / 1e9;
        var memGB = ranks * LIB_GB + nodes * (KB_NODE + ROOT_KB) / 1e6 + strataGB;
        var perRank = LIB_GB + (nodes * KB_NODE / 1e6 + strataGB) / ranks;
        var rootGB = perRank + nodes * ROOT_KB / 1e6;
        var level = coreH < 1 && memGB < 8 ? ["light", "a laptop, minutes"] :
          coreH < 50 && memGB < 64 ? ["moderate", "a workstation, hours"] :
          coreH < 2000 && memGB < 1000 ? ["heavy", "a cluster, tens to hundreds of processors"] :
          ["very heavy", "a national HPC facility (hundreds of processors)"];
        function wallTxt(s) {
          return s < 120 ? Math.round(s) + " s" : s < 7200 ? Math.round(s / 60) + " min" :
            s < 172800 ? (s / 3600).toPrecision(2) + " h" : (s / 86400).toPrecision(2) + " days";
        }
        var rows = [
          ["Mesh", fmt(nodes) + " nodes" + (domain === "global" ? " on the sphere" : "")],
          ["Steps", steps + (layers ? ", " + layers + " stratigraphic layers" : "")],
          ["Time per step", (cpuStep).toPrecision(2) + " s on 1 processor (" + factor.toFixed(2) +
            " x the base step per node)"],
          ["Total", (coreH < 10 ? coreH.toPrecision(2) : Math.round(coreH)) + " core-hours; about " +
            wallTxt(wall) + " on " + ranks + " processor" + (ranks > 1 ? "s" : "")],
          ["Memory", (memGB < 10 ? memGB.toPrecision(2) : Math.round(memGB)) + " GB in total, " +
            perRank.toPrecision(2) + " GB per processor (" + rootGB.toPrecision(2) + " GB on the first)" +
            (strataGB > 0.05 ? "; stratigraphy " + strataGB.toPrecision(2) + " GB" : "")]
        ];
        var t = el("table");
        rows.forEach(function (r) { t.appendChild(el("tr", {}, [el("th", {text: r[0]}), el("td", {text: r[1]})])); });
        out.appendChild(t);
        var badge = el("p", {"class": "modchart-level " + level[0].replace(" ", "-")});
        badge.appendChild(el("strong", {text: level[0].charAt(0).toUpperCase() + level[0].slice(1)}));
        badge.appendChild(document.createTextNode(": " + level[1] + "."));
        out.appendChild(badge);
        var notes = [];
        if (nodes > 2e5 && ranks === 1) notes.push("More than ~200,000 nodes: run with MPI (mpirun -np N).");
        if (layers * nodes > 5e8) notes.push("Many layers on a large mesh: a longer layer interval (time: strat) saves memory.");
        if (domain === "global" && v("dxg") < 10) notes.push("Below 10 km, global runs are HPC jobs: see the HPC scaling page.");
        notes.push("Order of magnitude only: scaled from the costs measured on a 47,672-node margin; the real cost depends on the landscape and the solvers.");
        notes.forEach(function (n) { out.appendChild(el("p", {"class": "modchart-hint", text: n})); });
      }
      return {node: node, update: update, setDomain: setDomain, fromExample: fromExample};
    }

    render();
  }

  // ---- filter of the reference pages (inputs and outputs): hides the
  // entries that do not match, the sections left empty and the contents box,
  // counts the matches and scrolls to the first one; a key name that matches
  // exactly is listed first by scrolling to it
  function initFilter(input) {
    var entries = Array.prototype.slice.call(document.querySelectorAll(".gospl-ref-entry"));
    var article = document.querySelector(".bd-article") || document.body;
    var sections = Array.prototype.slice.call(article.querySelectorAll("section"));
    var toc = article.querySelectorAll(".contents, nav.contents, aside.topic");
    var count = document.createElement("span");
    count.className = "gospl-ref-count";
    input.parentNode.insertBefore(count, input.nextSibling);
    var timer = null;
    function apply(scroll) {
      var q = input.value.trim().toLowerCase(), first = null, exact = null, n = 0;
      entries.forEach(function (e) {
        var head = e.querySelector(".gospl-ref-head");
        var key = head ? (head.querySelector("code") || head).textContent.trim().toLowerCase() : "";
        var hit = !q || e.textContent.toLowerCase().indexOf(q) >= 0;
        e.style.display = hit ? "" : "none";
        if (hit && q) {
          n++;
          if (!first) first = e;
          if (!exact && key === q) exact = e;
        }
      });
      sections.forEach(function (sec) {
        if (!q) { sec.style.display = ""; return; }
        var own = sec.querySelectorAll(".gospl-ref-entry");
        if (!own.length) return;
        var any = Array.prototype.some.call(own, function (e) { return e.style.display !== "none"; });
        sec.style.display = any ? "" : "none";
      });
      Array.prototype.forEach.call(toc, function (t) { t.style.display = q ? "none" : ""; });
      count.textContent = q ? (n ? n + (n > 1 ? " matches" : " match") : "no match") : "";
      var target = exact || first;
      if (scroll && target) target.scrollIntoView({behavior: "smooth", block: "center"});
    }
    input.addEventListener("input", function () {
      apply(false);
      clearTimeout(timer);
      timer = setTimeout(function () { apply(true); }, 500);
    });
    input.addEventListener("keydown", function (ev) {
      if (ev.key === "Enter") { ev.preventDefault(); apply(true); }
      if (ev.key === "Escape") { input.value = ""; apply(false); }
    });
  }

  document.addEventListener("DOMContentLoaded", function () {
    document.querySelectorAll(".gospl-modchart").forEach(init);
    document.querySelectorAll(".gospl-ref-filter").forEach(initFilter);
  });
})();
