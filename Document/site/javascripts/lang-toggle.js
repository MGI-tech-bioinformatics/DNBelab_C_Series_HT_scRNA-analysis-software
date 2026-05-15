document.addEventListener("DOMContentLoaded", () => {
  // Always sync Material palette with current system preference.
  // This prevents stale localStorage palette values from locking the theme.
  try {
    const isDark = window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
    const scheme = isDark ? "slate" : "default";
    const primary = "indigo";
    const accent = "blue";
    const palette = { color: { scheme, primary, accent } };
    localStorage.setItem("__palette", JSON.stringify(palette));
    document.body.setAttribute("data-md-color-scheme", scheme);
    document.body.setAttribute("data-md-color-primary", primary);
    document.body.setAttribute("data-md-color-accent", accent);
  } catch (e) {
    // no-op
  }

  const langInner = document.querySelector(".md-header__option .md-select__inner");
  if (!langInner) return;
  const option = langInner.closest(".md-header__option");
  if (!option) return;

  const toTargetPath = (pathname) => {
    if (pathname.includes("/en/")) return pathname.replace("/en/", "/");
    if (pathname.includes("/site/")) return pathname.replace("/site/", "/site/en/");
    if (pathname.startsWith("/en")) return pathname.replace(/^\/en/, "/");
    return `/en${pathname.startsWith("/") ? "" : "/"}${pathname}`;
  };

  option.addEventListener(
    "click",
    (event) => {
      event.preventDefault();
      event.stopPropagation();
      const { pathname, search, hash } = window.location;
      const next = toTargetPath(pathname);
      window.location.assign(`${next}${search}${hash}`);
    },
    true
  );
});
