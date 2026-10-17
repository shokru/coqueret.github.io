/* Small enhancements; every page still works without JavaScript. */
(function () {
  // 1. Underline the current page in the top menu
  var page = location.pathname.split("/").pop() || "index.html";
  document.querySelectorAll(".site-nav a").forEach(function (a) {
    if (a.getAttribute("href") === page) a.setAttribute("aria-current", "page");
  });

  // 2. Links to other sites open in a new tab
  document.querySelectorAll('a[href^="http"]').forEach(function (a) {
    if (a.hostname !== location.hostname) {
      a.target = "_blank";
      a.rel = "noopener";
    }
  });

  // 3. Topic filter for published articles (research page).
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
