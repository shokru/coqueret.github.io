/* Small enhancements; every page still works without JavaScript. */
(function () {
  var root = document.documentElement;

  // 1. Underline the current section in the top menu.
  //    Sub-pages (course or research material) point to their parent
  //    with the small link above their title (class "crumb").
  var page = location.pathname.split("/").pop() || "index.html";
  var crumb = document.querySelector(".crumb a");
  var current = crumb ? crumb.getAttribute("href") : page;
  document.querySelectorAll(".site-nav a").forEach(function (a) {
    if (a.getAttribute("href") === current) a.setAttribute("aria-current", "page");
  });

  // 2. Light / dark switch. The choice is remembered in the visitor's browser;
  //    until they use it, the site follows their system setting.
  var btn = document.querySelector(".theme-toggle");
  if (btn) {
    var media = window.matchMedia ? window.matchMedia("(prefers-color-scheme: dark)") : null;
    var effective = function () {
      return root.getAttribute("data-theme") || (media && media.matches ? "dark" : "light");
    };
    var label = function () {
      var text = effective() === "dark" ? "Switch to light mode" : "Switch to dark mode";
      btn.setAttribute("aria-label", text);
      btn.title = text;
    };
    btn.addEventListener("click", function () {
      var next = effective() === "dark" ? "light" : "dark";
      root.setAttribute("data-theme", next);
      try { localStorage.setItem("theme", next); } catch (e) {}
      label();
    });
    if (media && media.addEventListener) media.addEventListener("change", label);
    label();
    btn.hidden = false;
  }

  // 3. Links to other sites, notebooks, slides and papers open in a new tab
  document.querySelectorAll("a[href]").forEach(function (a) {
    var href = a.getAttribute("href");
    var external = /^https?:/.test(href) && a.hostname !== location.hostname;
    if (external || /^(files|papers)\//.test(href)) {
      a.target = "_blank";
      a.rel = "noopener";
    }
  });

  // 4. Topic filter for published articles (research page).
  //    Buttons are built from the "topic" values in data/articles.yml,
  //    so a new topic gets its own button automatically.
  var list = document.querySelector("#listing-articles .pubs");
  if (!list) return;
  var items = Array.prototype.slice.call(list.querySelectorAll(".pub"));
  var topics = [];
  items.forEach(function (it) {
    var t = it.getAttribute("data-topic");
    if (t && topics.indexOf(t) < 0) topics.push(t);
  });
  topics.sort(function (a, b) { return a.localeCompare(b); });

  var hasSelected = items.some(function (it) { return it.hasAttribute("data-selected"); });
  var filters = (hasSelected ? ["Selected"] : []).concat(["All"], topics);

  function matches(it, f) {
    if (f === "All") return true;
    if (f === "Selected") return it.hasAttribute("data-selected");
    return it.getAttribute("data-topic") === f;
  }

  var bar = document.createElement("div");
  bar.className = "filters";
  bar.setAttribute("role", "group");
  bar.setAttribute("aria-label", "Filter articles by topic");
  var chips = filters.map(function (f) {
    var b = document.createElement("button");
    b.type = "button";
    b.className = "chip";
    b.textContent = f;
    var n = document.createElement("span");
    n.className = "n";
    n.textContent = items.filter(function (it) { return matches(it, f); }).length;
    b.appendChild(n);
    b.addEventListener("click", function () { apply(f); });
    bar.appendChild(b);
    return b;
  });
  list.parentNode.insertBefore(bar, list);

  function apply(f) {
    items.forEach(function (it) { it.hidden = !matches(it, f); });
    chips.forEach(function (c, i) { c.setAttribute("aria-pressed", String(filters[i] === f)); });
  }
  apply(filters[0]);
})();
