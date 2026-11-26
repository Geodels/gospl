/* goSPL module chooser: pick the processes you want, see what they need.

   Data: window.GOSPL_MODULES = {modules: [{id, label, block, page, summary,
   group, always, requires: [{id, why, enforced}], requires_any: [[{id, why}],
   ...], excludes: [{id, why, enforced}], recommends: [{id, why}], mesh}]}
   written into the page by docs/user_guide/build_reference.py.

   The core is selected as soon as anything is. A selected module is green
   when everything it needs is selected too, red while something is missing;
   a missing requirement is drawn as a dashed "+ add" block (click it, or
   its checkbox, to add it). The graph is laid out in rows (what is
   needed above), ordered to limit crossings, and the links are curves
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
    // the always-on modules (the core and the rain it needs) are always
    // selected, even before anything is ticked
    function selected() {
      var s = {};
      Object.keys(picked).forEach(function (id) { s[id] = true; });
      if (coreId) {
        s[coreId] = true;
        if (rainId && !(oroId && s[oroId])) s[rainId] = true;
      }
      return s;
    }
    function automatic(id) {
      return id === coreId || (id === rainId && !picked[rainId]);
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
    var skeleton = el("div", {"class": "modchart-panel modchart-skel"});
    root.appendChild(skeleton);
    var lastSt = null;
    resources.onChange(function () { if (lastSt) writeSkeleton(lastSt); });

    function toggle(id, on) {
      if (id === coreId) return;
      if (id === rainId && !on && !picked[rainId]) {
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
      lastSt = st;
      writeSkeleton(st);
    }

    // ---- report under the graph
    function writeReport(st) {
      report.innerHTML = "";
      var ids = Object.keys(st.sel);
      if (!Object.keys(picked).length) {
        report.appendChild(el("p", {"class": "modchart-hint",
          text: "The graph shows what every run has (the core and its rain). Tick the " +
                "processes you want to use: a selected module turns green when everything " +
                "it needs is selected, red while something is missing; click a dashed " +
                "block to add it."}));
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

      // the core underlies every module: link it to those with no other
      // requirement shown (layout and the "what it makes possible" arrows)
      if (coreId && shown[coreId]) {
        nodes.forEach(function (id) {
          if (id === coreId || id === rainId || !st.sel[id]) return;
          var has = edges.some(function (e) { return e.from === id; });
          if (!has) edges.push({from: id, to: coreId, ok: true, alt: false, core: true});
        });
      }

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

      // neighbours of every point (one row up and one row down)
      var nbr = {}, up = {}, down = {};
      function push(m, k, v) { (m[k] = m[k] || []).push(v); }
      chains.forEach(function (ch) {
        // ids[k] is one row below ids[k + 1]
        for (var k = 0; k + 1 < ch.ids.length; k++) {
          push(nbr, ch.ids[k], ch.ids[k + 1]); push(nbr, ch.ids[k + 1], ch.ids[k]);
          push(up, ch.ids[k], ch.ids[k + 1]); push(down, ch.ids[k + 1], ch.ids[k]);
        }
      });

      // order of every row: barycentre sweeps down and up, then swaps of
      // neighbours, keeping the order with the fewest link crossings
      function index() {
        var idx = {};
        cols.forEach(function (row) { row.forEach(function (n, i) { idx[n.id] = i; }); });
        return idx;
      }
      function crossings() {
        var idx = index(), total = 0;
        for (var k = 1; k < cols.length; k++) {
          var segs = [];
          cols[k].forEach(function (n) {
            (up[n.id] || []).forEach(function (u) { segs.push([idx[u], idx[n.id]]); });
          });
          for (var i = 0; i < segs.length; i++) {
            for (var j = i + 1; j < segs.length; j++) {
              if ((segs[i][0] - segs[j][0]) * (segs[i][1] - segs[j][1]) < 0) total++;
            }
          }
        }
        return total;
      }
      function snapshot() { return cols.map(function (row) { return row.slice(); }); }
      function sweep(k, nb) {
        var idx = index(), row = cols[k];
        row.forEach(function (n, i) {
          var ns = nb[n.id] || [];
          n.key = ns.length ? ns.reduce(function (t, x) { return t + idx[x]; }, 0) / ns.length : i;
        });
        row.sort(function (p, q) { return p.key - q.key; });
      }
      var best = snapshot(), bestC = crossings();
      for (var it = 0; it < 24 && bestC > 0; it++) {
        if (it % 2 === 0) { for (var k1 = 1; k1 < cols.length; k1++) sweep(k1, up); }
        else { for (var k2 = cols.length - 2; k2 >= 0; k2--) sweep(k2, down); }
        var cc = crossings();
        if (cc < bestC) { bestC = cc; best = snapshot(); }
      }
      cols = best;
      for (var improved = true, guard = 0; improved && bestC > 0 && guard < 50; guard++) {
        improved = false;
        cols.forEach(function (row) {
          for (var i = 0; i + 1 < row.length; i++) {
            var t = row[i]; row[i] = row[i + 1]; row[i + 1] = t;
            var c2 = crossings();
            if (c2 < bestC) { bestC = c2; improved = true; }
            else { row[i + 1] = row[i]; row[i] = t; }
          }
        });
      }

      // geometry: rows from top to bottom, what a module needs above it, the
      // modules that need it below; a dummy point takes a narrow slot
      var W = 230, H = 46, GX = 28, ROWH = 110, DW = 22;
      function slot(n) { return n.dummy ? DW : W; }
      // x of every point: start packed in the barycentre order, then move each
      // point towards the mean of its neighbours (link points pulled harder, so
      // long links run straight) and push apart the ones that overlap, keeping
      // the order of every row
      var cx = {};
      cols.forEach(function (row) {
        var x = 0;
        row.forEach(function (n) { cx[n.id] = x + slot(n) / 2; x += slot(n) + GX; });
      });
      function separate(row) {
        for (var i = 1; i < row.length; i++) {          // left to right
          var a = row[i - 1], b = row[i], m = (slot(a) + slot(b)) / 2 + GX;
          if (cx[b.id] < cx[a.id] + m) cx[b.id] = cx[a.id] + m;
        }
        var shift = 0;                                   // keep the row centred on its wish
        row.forEach(function (n) { shift += n.want - cx[n.id]; });
        shift /= row.length;
        row.forEach(function (n) { cx[n.id] += shift; });
        for (var j = row.length - 2; j >= 0; j--) {      // right to left
          var c = row[j], e = row[j + 1], m2 = (slot(c) + slot(e)) / 2 + GX;
          if (cx[c.id] > cx[e.id] - m2) cx[c.id] = cx[e.id] - m2;
        }
      }
      for (var pass = 0; pass < 30; pass++) {
        cols.forEach(function (row) {
          row.forEach(function (n) {
            var ns = nbr[n.id] || [];
            if (!ns.length) { n.want = cx[n.id]; return; }
            var sw = 0, sx = 0;
            ns.forEach(function (x) { var w = x.charAt(0) === "_" ? 2 : 1; sw += w; sx += w * cx[x]; });
            n.want = sx / sw;
            cx[n.id] = n.dummy ? n.want : 0.5 * (cx[n.id] + n.want);
          });
          separate(row);
        });
      }
      // close the empty vertical bands left between unconnected groups
      var spans = [];
      cols.forEach(function (row) {
        row.forEach(function (n) { spans.push([cx[n.id] - slot(n) / 2, cx[n.id] + slot(n) / 2]); });
      });
      spans.sort(function (p, q) { return p[0] - q[0]; });
      var cuts = [], reach = -Infinity;
      spans.forEach(function (sp) {
        if (reach > -Infinity && sp[0] - reach > 2 * GX) cuts.push([reach, sp[0] - reach - 2 * GX]);
        reach = Math.max(reach, sp[1]);
      });
      cols.forEach(function (row) {
        row.forEach(function (n) {
          var x = cx[n.id], dx = 0;
          cuts.forEach(function (c) { if (x > c[0]) dx += c[1]; });
          cx[n.id] = x - dx;
        });
      });
      var xmin = Infinity, xmax = -Infinity;
      cols.forEach(function (row) {
        row.forEach(function (n) {
          xmin = Math.min(xmin, cx[n.id] - slot(n) / 2);
          xmax = Math.max(xmax, cx[n.id] + slot(n) / 2);
        });
      });
      var width = Math.max(xmax - xmin, W) + 24;
      var height = ncol * ROWH - (ROWH - H) + 24;
      var pos = {};
      cols.forEach(function (row, k) {
        var y = 12 + k * ROWH;
        row.forEach(function (n) {
          pos[n.id] = {x: 12 + cx[n.id] - slot(n) / 2 - xmin, y: y, w: slot(n), dummy: !!n.dummy};
        });
      });

      // drawn at its natural size (never enlarged, never shrunk): the box
      // scrolls sideways when the graph is wider than the page
      var s = svg("svg", {viewBox: "0 0 " + width + " " + height, width: width, height: height,
                          "class": "modchart-svg", role: "img",
                          "aria-label": "Selected modules and what they need"});
      var defs = svg("defs");
      ["ok", "bad"].forEach(function (k) {
        var mk = svg("marker", {id: "mc-arrow-" + k, viewBox: "0 0 10 10", refX: 10, refY: 5,
                                markerWidth: 11, markerHeight: 11, markerUnits: "userSpaceOnUse",
                                orient: "auto"});
        mk.appendChild(svg("path", {d: "M0 0 L10 5 L0 10 z", fill: k === "ok" ? "#6b8f71" : "#c0392b"}));
        defs.appendChild(mk);
      });
      s.appendChild(defs);

      // links: from the bottom of the needed module to the top of the one
      // that needs it (the arrow shows what a module makes possible), straight
      // down through the slots of the rows they cross. The links leave and
      // reach a box at ports spread along its edge, in the order of the other
      // ends, so they do not cross at the box, and end with a short vertical
      // run so the arrow follows the line.
      function ports(id, others) {
        var o = others.slice().sort(function (p, q) { return p.x - q.x; });
        var p = pos[id], m = o.length, res = {};
        o.forEach(function (it, i) {
          res[it.key] = p.dummy ? p.x + p.w / 2 : p.x + p.w * (0.5 + (m > 1 ? (i / (m - 1) - 0.5) * 0.6 : 0));
        });
        return res;
      }
      var topPort = {}, botPort = {};
      chains.forEach(function (ch, c) {
        var a = ch.ids[0], b = ch.ids[ch.ids.length - 1];
        var na = pos[ch.ids[1]], nb = pos[ch.ids[ch.ids.length - 2]];
        push(topPort, a, {key: c, x: nb.x + nb.w / 2});
        push(botPort, b, {key: c, x: na.x + na.w / 2});
      });
      var tp = {}, bp = {};
      Object.keys(topPort).forEach(function (id) { tp[id] = ports(id, topPort[id]); });
      Object.keys(botPort).forEach(function (id) { bp[id] = ports(id, botPort[id]); });
      var STUB = 10;
      chains.forEach(function (ch, c) {
        // from the needed module (top) to the one that needs it (bottom)
        var ids = ch.ids.slice().reverse(), last = ids.length - 1;
        var pts = [];
        ids.forEach(function (id, k) {
          var p = pos[id];
          if (k === 0) { var x0 = bp[id][c]; pts.push([x0, p.y + H], [x0, p.y + H + STUB]); }
          else if (k === last) { var x1 = tp[id][c]; pts.push([x1, p.y - STUB], [x1, p.y]); }
          else { var xd = p.x + p.w / 2; pts.push([xd, p.y], [xd, p.y + H]); }
        });
        var dpath = "M" + pts[0][0] + " " + pts[0][1];
        for (var k = 1; k < pts.length; k++) {
          var a0 = pts[k - 1], b0 = pts[k];
          if (k % 2 === 1) { dpath += " L" + b0[0] + " " + b0[1]; continue; }   // vertical run
          var my = (a0[1] + b0[1]) / 2;
          dpath += " C" + a0[0] + " " + my + " " + b0[0] + " " + my + " " + b0[0] + " " + b0[1];
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
        g.appendChild(svg("rect", {x: p.x, y: p.y, width: p.w, height: H, rx: 7, ry: 7}));
        // short label in the box (no bracketed part); the full one in the tooltip
        var lab = (cls === "missing" ? "+ " : "") + (id === coreId ? "Core: mesh, flow, rivers"
          : m.label.replace(/\s*\([^)]*\)\s*/g, " ").trim());
        var t = svg("text", {x: p.x + p.w / 2, y: p.y + H / 2 + 5.5, "text-anchor": "middle"});
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
      legend.appendChild(el("span", {text: "Arrows go down from a module to the modules that need it; dashed: one of several options."}));
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

    // ---- skeleton input file of the selection, and where to read more
    var BLOCK_ORDER = ["", "domain", "time", "output", "spl", "diffusion", "sea", "climate",
      "orography", "sedfactor", "soil", "compaction", "strata", "provenance", "gravity",
      "carbonate", "sst", "waves", "coastal_erosion", "pelagic", "carbonate_supply", "karst",
      "groundwater", "ice", "tectonics", "flexure", "temap"];
    function num(x) {
      var s = String(+x.toPrecision(6));
      if (/e/.test(s)) return s.replace(/e\+?/, ".e").replace(".e-", ".e-").replace("..", ".");
      return s.indexOf(".") < 0 ? s + "." : s;
    }
    function buildYaml(st) {
      var r = resources.values(), end = r.dur, dt = r.dt || 1;
      var tout = Math.max(dt, Math.round(end / 10 / dt) * dt);
      var fill = function (v) {
        return v.replace(/\{end\}/g, num(end)).replace(/\{dt\}/g, num(dt))
          .replace(/\{tout\}/g, num(tout)).replace(/\{strat\}/g, num(r.strat || dt));
      };
      var ids = Object.keys(st.sel);
      ids.sort(function (a, b) { return a === coreId ? -1 : b === coreId ? 1 : 0; });
      // tree: block path -> ordered keys
      var tree = {}, order = [];
      function node(path) {
        if (!tree[path]) { tree[path] = {keys: [], vals: {}}; order.push(path); }
        return tree[path];
      }
      ids.forEach(function (id) {
        (mods[id].yaml || []).forEach(function (e) {
          if (e.key === "bc" && domain === "global") return;
          // a top-level list (climate, tectonics, pelagic, ...) is a block of its own
          var path = e.block === "" && e.key !== "name" ? "#" + e.key : e.block;
          var n = node(path), v = fill(e.value);
          if (path === "#climate" && st.sel.evaporation)
            v += "\n    evap_uniform: 0.2   # evaporation (m/yr)";
          if (!(e.key in n.vals)) n.keys.push(e.key);
          n.vals[e.key] = {v: v, c: e.comment, id: id};
        });
      });
      // on a sphere the waves need a fetch (no domain edge brings them in)
      if (tree.waves && domain === "global" && !("fetch" in tree.waves.vals)) {
        tree.waves.keys.unshift("fetch");
        tree.waves.vals.fetch = {v: "5.e5", c: "m; growth length to full energy over open water (required on a global mesh)"};
      }
      // goSPL needs both creep coefficients whenever the diffusion block is
      // there (marine keys or channelQ alone made it stop at start-up)
      if (tree.diffusion) {
        [["hillslopeKm", "Linear diffusivity of the seafloor creep"],
         ["hillslopeKa", "Linear hillslope diffusivity on land"]].forEach(function (k) {
          if (!(k[0] in tree.diffusion.vals)) {
            tree.diffusion.keys.unshift(k[0]);
            tree.diffusion.vals[k[0]] = {v: "0.", c: "m2/yr; " + k[1] +
              " (required in this block; 0 = none)"};
          }
        });
      }
      // the block of a sub-block exists even when the parent has no key
      order.slice().forEach(function (p) {
        var parts = p.split(".");
        for (var k = 1; k < parts.length; k++) node(parts.slice(0, k).join("."));
      });
      var tops = order.filter(function (p) { return p.indexOf(".") < 0; });
      tops.sort(function (a, b) {
        var ia = BLOCK_ORDER.indexOf(a.replace("#", "")), ib = BLOCK_ORDER.indexOf(b.replace("#", ""));
        return (ia < 0 ? 99 : ia) - (ib < 0 ? 99 : ib);
      });
      var L = [];
      L.push("# goSPL input file: skeleton for a " + (domain === "global" ? "global (spherical)" : "regional (flat)") + " model");
      L.push("# made by 'Which modules do I need?' in the goSPL documentation, with the minimal keys of:");
      L.push("#   " + ids.map(function (id) { return mods[id].label; }).join(", "));
      L.push("# Values are starting points: check each one (units and ranges: the input reference page).");
      L.push("# The mesh (domain: npdata) and any map ([file, key]) are npz files you build first.");
      if (domain === "global") L.push("# Global model: the mesh is a sphere; there are no edge conditions (domain: bc).");
      L.push("");
      function emit(path, indent) {
        var n = tree[path], pad = new Array(indent + 1).join(" ");
        n.keys.forEach(function (k) {
          var e = n.vals[k];
          if (e.v.charAt(0) === "\n") {
            L.push(pad + k + ":" + (e.c ? "   # " + e.c : ""));
            e.v.slice(1).split("\n").forEach(function (line) { L.push(pad + line); });
          } else {
            L.push(pad + k + ": " + e.v + (e.c ? "   # " + e.c : ""));
          }
        });
        order.filter(function (p) {
          return p.indexOf(path + ".") === 0 && p.split(".").length === path.split(".").length + 1;
        }).forEach(function (sub) {
          L.push(pad + sub.split(".").pop() + ":");
          emit(sub, indent + 4);
        });
      }
      tops.forEach(function (t) {
        if (t === "" || t.charAt(0) === "#") { emit(t, 0); L.push(""); return; }
        L.push(t + ":");
        emit(t, 4);
        L.push("");
      });
      return L.join("\n");
    }
    // where to get the mesh (domain: npdata) for the chosen kind of model
    function meshHelp() {
      var ml = data.meshes || {}, r = window.GOSPL_DOCROOT || "../";
      var box = el("div", {"class": "modchart-mesh"});
      box.appendChild(el("strong", {text: "The mesh (domain: npdata). "}));
      function a(href, txt) { var x = el("a", {href: href}); x.textContent = txt; return x; }
      function t(x) { return document.createTextNode(x); }
      if (domain === "global") {
        [t("Download a uniform global mesh (5, 8, 10 or 20 km spacing) from "),
         a(ml.global_download, "the goSPL meshes folder"),
         t(", or build your own (any resolution, refined where you need it) with "),
         a(r + ml.global_build + ".html", "the continental_flux build_inputs notebook"),
         t(". To put your NetCDF data (elevation, rain, uplift, ...) on the mesh, see "),
         a(r + ml.global_data + ".html", "Global data on a goSPL mesh"),
         t(" (longitudes -180 to 180 or 0 to 360, latitudes -90 to 90, both ends included).")]
          .forEach(function (n) { box.appendChild(n); });
      } else {
        [t("Build a regional mesh (from a DEM, refined or uniform) with "),
         a(r + ml.regional_build + ".html", "the Building 2D meshes example"),
         t(".")].forEach(function (n) { box.appendChild(n); });
      }
      return box;
    }

    function writeSkeleton(st) {
      skeleton.innerHTML = "";
      if (!Object.keys(picked).length) return;
      skeleton.appendChild(el("h4", {text: "Your input file (skeleton)"}));
      var yaml = buildYaml(st);
      var bar = el("div", {"class": "modchart-skelbar"});
      var dl = el("button", {type: "button", "class": "modchart-add", text: "Download input.yml"});
      dl.addEventListener("click", function () {
        var blob = new Blob([yaml + "\n"], {type: "text/yaml"});
        var a = document.createElement("a");
        a.href = URL.createObjectURL(blob);
        a.download = "input.yml";
        document.body.appendChild(a);
        a.click();
        setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 0);
      });
      var cp = el("button", {type: "button", "class": "modchart-add", text: "Copy"});
      cp.addEventListener("click", function () {
        if (navigator.clipboard) navigator.clipboard.writeText(yaml + "\n").then(function () {
          cp.textContent = "Copied"; setTimeout(function () { cp.textContent = "Copy"; }, 1500);
        });
      });
      bar.appendChild(dl);
      bar.appendChild(cp);
      bar.appendChild(el("span", {"class": "modchart-hint",
        text: "  End time, step, output and layer intervals come from the size form above."}));
      skeleton.appendChild(bar);
      skeleton.appendChild(meshHelp());
      skeleton.appendChild(el("pre", {"class": "modchart-yaml", text: yaml}));
      if (Object.keys(st.missing).length)
        skeleton.appendChild(el("p", {"class": "modchart-hint modchart-bad",
          text: "Some modules still need others (in red above): the file is not complete until they are added."}));

      // where to read more about each selected module
      var root_ = window.GOSPL_DOCROOT || "../";
      skeleton.appendChild(el("h4", {text: "Read more about your modules"}));
      var t = el("table", {"class": "modchart-links"});
      t.appendChild(el("tr", {}, [el("th", {text: "module"}), el("th", {text: "theory"}),
        el("th", {text: "how to set it"}), el("th", {text: "its keys"})]));
      Object.keys(st.sel).forEach(function (id) {
        var m = mods[id], anchors = [];
        (m.yaml || []).forEach(function (e) { if (e.anchor && anchors.indexOf(e.anchor) < 0) anchors.push(e.anchor); });
        var cell = function (href, txt) {
          if (!href) return el("td", {text: "-"});
          var a = el("a", {href: href});
          a.textContent = txt;
          return el("td", {}, [a]);
        };
        var keys = el("td");
        anchors.forEach(function (an, i) {
          var a = el("a", {href: root_ + "user_guide/reference/inputs.html#" + an});
          a.textContent = an.replace(/^in-/, "").replace(/-/g, ".");
          if (i) keys.appendChild(document.createTextNode(", "));
          keys.appendChild(a);
        });
        if (!anchors.length) keys.textContent = "-";
        var row = el("tr", {}, [el("td", {text: m.label}),
          cell(m.tech ? root_ + m.tech + ".html" : "", m.tech ? m.tech.split("/").pop().replace("ero", "erosion").replace("dep", "depressions").replace("hill", "hillslope").replace("strat", "stratigraphy").replace("tecto", "tectonics").replace("carb", "carbonates") : ""),
          cell(m.page ? root_ + m.page + ".html" : "", "user guide"), keys]);
        t.appendChild(row);
      });
      skeleton.appendChild(t);
    }

    // ---- resources: nodes, steps, time and memory from the measured costs
    function makeResources() {
      // base step: 1.67 s at 1 rank for 47,672 nodes (cost page), of which
      // the marine deposition is ~64 %
      var NODE_S = 1.6664 / 47672, MARINE = 0.64, LIB_GB = 0.56, KB_NODE = 6.5, ROOT_KB = 1.8;
      var onChange = null;
      var node = el("div", {"class": "modchart-panel modchart-res"});
      node.appendChild(el("h4", {text: "Your run: size and resources"}));
      var form = el("div", {"class": "modchart-form"}), out = el("div", {"class": "modchart-out"});
      var fields = {};
      function field(key, label, val, unit, cls) {
        var inp = el("input", {type: "number", min: "0", step: "any", value: val});
        inp.addEventListener("input", function () { update(); if (onChange) onChange(); });
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
        if (onChange) onChange();
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
        if (onChange) onChange();
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
      function values() {
        var v = function (k) { return Math.max(0, parseFloat(fields[k].value) || 0); };
        return {dur: v("dur") * 1e6, dt: v("dt") * 1e3, strat: v("strat") * 1e3};
      }
      return {node: node, update: update, setDomain: setDomain, fromExample: fromExample,
              values: values, onChange: function (f) { onChange = f; }};
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
