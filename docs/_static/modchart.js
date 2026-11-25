/* goSPL module chooser: pick the processes you want, see what they need.

   Data: window.GOSPL_MODULES = {modules: [{id, label, block, page, summary,
   group, requires: [{id, why, enforced}], requires_any: [[{id, why}], ...],
   excludes: [{id, why, enforced}], recommends: [{id, why}], mesh}]}
   written into the page by docs/user_guide/build_reference.py.

   Selecting a module adds everything it requires (transitively); the graph
   shows the selected modules, what they pull in, and conflicts. */
(function () {
  "use strict";

  function el(tag, attrs, children) {
    var e = document.createElement(tag);
    if (attrs) {
      Object.keys(attrs).forEach(function (k) {
        if (k === "text") e.textContent = attrs[k];
        else if (k === "html") e.innerHTML = attrs[k];
        else e.setAttribute(k, attrs[k]);
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
    var mods = {};
    data.modules.forEach(function (m) { mods[m.id] = m; });
    var picked = {};

    // closure of the requirements of the picked modules: id -> list of reasons
    function closure() {
      var need = {}, stack = Object.keys(picked);
      stack.forEach(function (id) { need[id] = need[id] || []; });
      while (stack.length) {
        var id = stack.pop(), m = mods[id];
        if (!m) continue;
        (m.requires || []).forEach(function (r) {
          if (!mods[r.id]) return;
          var first = !(r.id in need);
          need[r.id] = need[r.id] || [];
          need[r.id].push({by: id, why: r.why, enforced: r.enforced});
          if (first) stack.push(r.id);
        });
      }
      return need;
    }

    // ---- controls: one checkbox per module, by group
    var groups = {};
    data.modules.forEach(function (m) {
      if (m.always) return;
      var g = m.group || "Other";
      (groups[g] = groups[g] || []).push(m);
    });
    var controls = el("div", {"class": "modchart-controls"});
    Object.keys(groups).forEach(function (g) {
      var box = el("fieldset", {"class": "modchart-group"}, [el("legend", {text: g})]);
      groups[g].forEach(function (m) {
        var cb = el("input", {type: "checkbox", id: "mc-" + m.id});
        cb.addEventListener("change", function () {
          if (cb.checked) picked[m.id] = true; else delete picked[m.id];
          render();
        });
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
    root.appendChild(controls);
    root.appendChild(graph);
    root.appendChild(report);

    function link(m) {
      if (!m.page) return el("strong", {text: m.label});
      var a = el("a", {href: (window.GOSPL_DOCROOT || "../") + m.page + ".html"});
      a.appendChild(el("strong", {text: m.label}));
      return a;
    }

    function render() {
      var need = closure();
      var ids = Object.keys(need);
      // mark the controls of modules pulled in by others
      data.modules.forEach(function (m) {
        var item = root.querySelector("#mci-" + m.id);
        if (!item) return;
        item.classList.toggle("is-required", m.id in need && !picked[m.id]);
      });

      // ---- report
      report.innerHTML = "";
      if (!ids.length) {
        report.appendChild(el("p", {"class": "modchart-hint",
          text: "Tick the processes you want to use: the chart adds what they need " +
                "and lists what does not go with them."}));
        drawGraph([], need);
        return;
      }
      var conflicts = [], recs = [], anys = [], meshes = [];
      ids.forEach(function (id) {
        var m = mods[id];
        (m.excludes || []).forEach(function (x) {
          if (x.id in need) conflicts.push({a: m, b: mods[x.id], why: x.why});
        });
        (m.recommends || []).forEach(function (r) {
          if (!(r.id in need) && mods[r.id]) recs.push({a: m, b: mods[r.id], why: r.why});
        });
        (m.requires_any || []).forEach(function (alts) {
          var ok = alts.some(function (r) { return r.id in need; });
          if (!ok) anys.push({a: m, alts: alts});
        });
        if (m.mesh && m.mesh !== "both") meshes.push({a: m, mesh: m.mesh});
      });

      var ul = el("ul");
      ids.sort(function (a, b) { return (picked[b] ? 1 : 0) - (picked[a] ? 1 : 0); });
      ids.forEach(function (id) {
        var m = mods[id], li = el("li", {}, [link(m)]);
        if (m.block) li.appendChild(el("code", {"class": "modchart-block", text: m.block}));
        if (picked[id]) li.appendChild(el("span", {"class": "modchart-tag picked", text: "selected"}));
        else {
          li.appendChild(el("span", {"class": "modchart-tag required", text: "required"}));
          var why = el("ul", {"class": "modchart-why"});
          need[id].forEach(function (r) {
            why.appendChild(el("li", {text: "by " + mods[r.by].label + ": " + r.why +
              (r.enforced ? " (" + r.enforced + ")" : "")}));
          });
          li.appendChild(why);
        }
        ul.appendChild(li);
      });
      report.appendChild(el("h4", {text: "Modules to set in the input file"}));
      report.appendChild(ul);

      function section(title, cls, items, fmt) {
        if (!items.length) return;
        report.appendChild(el("h4", {"class": cls, text: title}));
        var u = el("ul");
        items.forEach(function (it) { u.appendChild(el("li", {text: fmt(it)})); });
        report.appendChild(u);
      }
      section("Not compatible: remove one of them", "modchart-bad", conflicts, function (c) {
        return c.a.label + " and " + c.b.label + ": " + c.why;
      });
      section("Choose at least one of", "modchart-warn", anys, function (c) {
        return c.a.label + ": " + c.alts.map(function (r) {
          return (mods[r.id] ? mods[r.id].label : r.id) + (r.why ? " (" + r.why + ")" : "");
        }).join(" or ");
      });
      section("Mesh", "modchart-info", meshes, function (c) {
        return c.a.label + ": " + c.mesh;
      });
      section("Often used with", "modchart-info", recs, function (c) {
        return c.b.label + " (for " + c.a.label + "): " + c.why;
      });
      drawGraph(ids, need);
    }

    // ---- graph: modules in columns by depth of requirement, arrows from a
    // module to what it requires
    function drawGraph(ids, need) {
      graph.innerHTML = "";
      if (!ids.length) return;
      var depth = {};
      function d(id, seen) {
        if (id in depth) return depth[id];
        seen = seen || {};
        if (seen[id]) return 0;
        seen[id] = true;
        var m = mods[id], k = 0;
        (m.requires || []).forEach(function (r) {
          if (r.id in need) k = Math.max(k, d(r.id, seen) + 1);
        });
        depth[id] = k;
        return k;
      }
      ids.forEach(function (id) { d(id); });
      var cols = {};
      ids.forEach(function (id) { (cols[depth[id]] = cols[depth[id]] || []).push(id); });
      var ncol = Math.max.apply(null, Object.keys(cols).map(Number)) + 1;
      var W = 190, H = 34, GX = 60, GY = 14;
      var nrow = Math.max.apply(null, Object.keys(cols).map(function (k) { return cols[k].length; }));
      var width = ncol * W + (ncol - 1) * GX + 20, height = nrow * (H + GY) + 20;
      var s = svg("svg", {viewBox: "0 0 " + width + " " + height, width: "100%",
                          "class": "modchart-svg", role: "img",
                          "aria-label": "Selected modules and their requirements"});
      var defs = svg("defs");
      var mk = svg("marker", {id: "mc-arrow", viewBox: "0 0 10 10", refX: 10, refY: 5,
                              markerWidth: 8, markerHeight: 8, markerUnits: "userSpaceOnUse", orient: "auto"});
      mk.appendChild(svg("path", {d: "M0 0 L10 5 L0 10 z", fill: "#898781"}));
      defs.appendChild(mk);
      s.appendChild(defs);
      var pos = {};
      Object.keys(cols).forEach(function (k) {
        // deepest requirements on the left, selected modules on the right
        var c = ncol - 1 - Number(k);
        cols[k].forEach(function (id, i) {
          pos[id] = {x: 10 + c * (W + GX), y: 10 + i * (H + GY) + (nrow - cols[k].length) * (H + GY) / 2};
        });
      });
      ids.forEach(function (id) {
        (mods[id].requires || []).forEach(function (r) {
          if (!(r.id in pos)) return;
          var a = pos[id], b = pos[r.id];
          var x1 = a.x, y1 = a.y + H / 2, x2 = b.x + W, y2 = b.y + H / 2, mx = (x1 + x2) / 2;
          s.appendChild(svg("path", {d: "M" + x1 + " " + y1 + " C" + mx + " " + y1 + " " + mx + " " +
                                       y2 + " " + x2 + " " + y2, "class": "modchart-edge",
                                     "marker-end": "url(#mc-arrow)"}));
        });
      });
      ids.forEach(function (id) {
        var p = pos[id], g = svg("g", {"class": "modchart-node " + (picked[id] ? "picked" : "required")});
        g.appendChild(svg("rect", {x: p.x, y: p.y, width: W, height: H, rx: 6, ry: 6}));
        var t = svg("text", {x: p.x + W / 2, y: p.y + H / 2 + 4, "text-anchor": "middle"});
        var lab = mods[id].label;
        t.textContent = lab.length > 26 ? lab.slice(0, 25) + "…" : lab;
        var tt = svg("title");
        tt.textContent = mods[id].label + (mods[id].summary ? ": " + mods[id].summary : "");
        g.appendChild(tt);
        g.appendChild(t);
        s.appendChild(g);
      });
      graph.appendChild(s);
      graph.appendChild(el("p", {"class": "modchart-legend",
        text: "Arrows point from a module to what it needs. Dark: selected; light: required."}));
    }

    render();
  }

  // ---- search filter for the reference pages (inputs and outputs)
  function initFilter(input) {
    var target = document.querySelectorAll(".gospl-ref-entry");
    input.addEventListener("input", function () {
      var q = input.value.trim().toLowerCase();
      target.forEach(function (e) {
        e.style.display = !q || e.textContent.toLowerCase().indexOf(q) >= 0 ? "" : "none";
      });
    });
  }

  document.addEventListener("DOMContentLoaded", function () {
    document.querySelectorAll(".gospl-modchart").forEach(init);
    document.querySelectorAll(".gospl-ref-filter").forEach(initFilter);
  });
})();
