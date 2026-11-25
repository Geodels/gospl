/* Section navigation (left sidebar): one branch open at a time. Opening a
   chevron folds the other open branches at the same level, so a section
   folds back when you move on to another one. */
document.addEventListener("DOMContentLoaded", function () {
  var nav = document.querySelector(".bd-sidenav");
  if (!nav) return;
  nav.querySelectorAll("details").forEach(function (d) {
    d.addEventListener("toggle", function () {
      if (!d.open) return;
      var li = d.parentElement, list = li && li.parentElement;
      if (!list) return;
      Array.prototype.forEach.call(list.children, function (other) {
        if (other === li) return;
        other.querySelectorAll(":scope > details[open]").forEach(function (o) {
          o.open = false;
        });
      });
    });
  });
});

/* "Contents" button below the top bar (computers): opens the same section
   navigation panel as the header menu icon, which the CSS hides there. */
document.addEventListener("DOMContentLoaded", function () {
  var toggle = document.querySelector(".primary-toggle");
  // only on pages with a section list (not the main page, search, index)
  if (!toggle || !document.querySelector("#pst-primary-sidebar .bd-sidenav")) return;
  document.documentElement.classList.add("gospl-has-contents");
  var btn = document.createElement("button");
  btn.type = "button";
  btn.className = "gospl-contents-btn";
  btn.setAttribute("aria-label", "Section navigation");
  btn.innerHTML = '<i class="fa-solid fa-bars" aria-hidden="true"></i><span>Contents</span>';
  btn.addEventListener("click", function () { toggle.click(); });
  document.body.appendChild(btn);
  // just below the top bar, whose height grows when its links wrap
  var header = document.querySelector(".bd-header");
  function place() {
    if (header) btn.style.top = (header.getBoundingClientRect().height + 16) + "px";
  }
  place();
  window.addEventListener("resize", place);
});
