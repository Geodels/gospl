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
